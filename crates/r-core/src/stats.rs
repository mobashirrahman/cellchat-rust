//! Bit-exact port of the summary statistics `computeCommunProb` uses.
//!
//! Parity contract: **Exact** (`docs/SEMANTICS.md` R2, R3, R4, R9, R10, R11).
//!
//! Upstream sources:
//! * `triMean`          — `R/modeling.R:878`, via `collapse::fquantile`
//! * `truncatedMean`    — `R/modeling.R:70`, via `mean(x, trim = trim, na.rm = TRUE)`
//! * `median`           — `R/modeling.R:72`, via `stats::median`
//! * `geometricMean`    — `R/modeling.R:863`, via `exp(mean(log(x), na.rm = TRUE))`
//! * `thresholdedMean`  — `R/modeling.R:889`
//!
//! ## The type-7 quantile discrepancy
//!
//! `collapse::fquantile` and R's `stats::quantile(type = 7)` compute the *same*
//! mathematical quantity but not the *same* floating-point expression:
//!
//! | | `h` | interpolation |
//! |---|---|---|
//! | R `ap.c` | `(n - 1)*p + 1` (1-based) | `lo + g*(hi - lo)`, `g = h - floor(h)` |
//! | `collapse` | `(l - 1)*Q` (0-based) | `a + (h - ih)*(b - a)`, `ih = (int)h` |
//!
//! The `+ 1` in R's form is a rounding opportunity that collapse avoids. Measured over
//! 50 000 random vectors, **7.4 % of cases differ, by up to 4.44e-16 (2 ulp)**. So
//! porting R's formula would be wrong: this module implements *collapse's*.
//!
//! Two more details of `fquantileC` that are load-bearing:
//!
//! * the boundary guard is `ih == n - 1` where `n` is the **original** `length(x)`,
//!   not the post-`na.rm` count `l`. Benign (see [`fquantile_type7`]), but replicated.
//! * `fquantile` switches to a radix-order fast path when `length(x) > 1e5 &&
//!   length(probs) > log(length(x))` (`R/my_RcppExports.R:107`). The result is identical
//!   either way, so this port always takes the direct path and only asserts the
//!   threshold is understood.

/// `collapse`'s default `eps = 10 * DBL_EPSILON` (from `fnth_fmedian_fquantile.c:16`).
/// Only used by the weighted paths, which CellChat never reaches.
pub const COLLAPSE_EPS: f64 = 10.0 * f64::EPSILON;

/// Which `fquantile` algorithm family `collapse` dispatches to.
///
/// Reproduced from `R/my_RcppExports.R:107` so the threshold is pinned by a test
/// rather than being an assumption.
#[inline]
pub fn collapse_uses_radix_order(len_x: usize, len_probs: usize) -> bool {
    (len_x > 100_000) && (len_probs as f64 > (len_x as f64).ln())
}

/// `collapse::fquantile(x, probs, na.rm = TRUE, type = 7)`, unweighted.
///
/// Port of **`dquickselect`** in `src/fnth_fmedian_fquantile.c:414` — *not* of the
/// `FQUANTILE_ORDVEC` macro that also lives in that file. The distinction is not
/// cosmetic:
///
/// * `fquantile` reaches `dquickselect` whenever `o = NULL`, which is always the case
///   for `triMean`: its default `o` (`radixorder(x)`) requires
///   `length(probs) > log(length(x))`, and with `length(probs) == 4` that needs
///   `length(x) < e^4`, contradicting the `length(x) > 1e5` the same condition demands.
///   [`collapse_uses_radix_order`] pins this.
/// * `dquickselect` subtracts the integer part from `h` **before** the zero test:
///   ```c
///   RETQSWITCH(n);            // type 7: h = (n - 1) * Q
///   elem = h; h -= elem;      // h is now the interpolation weight in [0, 1)
///   QUICKSELECT(dswap);       // places the order statistic at x[elem]
///   if (... || elem == n-1 || h <= 0.0) return a;
///   b = x[elem+1]; /* ... */;
///   return a + h*(b-a);
///   ```
///   The `FQUANTILE_ORDVEC` macro instead tests the *un*-reduced `h <= 0.0`. That
///   difference is observable whenever `b` is infinite: with a zero weight the macro
///   would evaluate `0.0 * Inf = NaN`, whereas `dquickselect` returns `a` before ever
///   touching `b`. Concretely, for `c(1, Inf, 2)` at `Q = 0.5` collapse returns `2.0`,
///   which is only reachable via the `dquickselect` guard.
///
/// * `n` here is the count of **non-missing** values, and `QUICKSELECT` leaves the
///   selected order statistic at `x[elem]` with `x[elem+1]` the next one up, so passing
///   a fully sorted slice is equivalent — and much cheaper.
pub fn fquantile_type7(x_sorted: &[f64], probs: &[f64], _n: usize) -> Vec<f64> {
    let n = x_sorted.len();
    let mut out = Vec::with_capacity(probs.len());
    let mut prev = f64::NEG_INFINITY;
    for &q in probs {
        assert!(
            (0.0..=1.0).contains(&q),
            "probabilities need to be in range [0, 1]"
        );
        assert!(
            q >= prev,
            "probabilities need to be passed in ascending order"
        );
        prev = q;
        if n == 0 {
            // `dquickselect`: `if (n == 0) return NA_REAL;`
            out.push(f64::NAN);
            continue;
        }
        // RETQSWITCH, type 7. `n` is a count, promoted to double for the multiply.
        let h_full = (n as f64 - 1.0) * q;
        // `elem = h;` — C truncation toward zero; h >= 0 so this is floor.
        let elem = h_full as usize;
        let h = h_full - elem as f64; // the interpolation weight
        let a = x_sorted[elem.min(n - 1)];
        // `elem == n-1` guards the last position, which quickselect can only select for
        // Q == 1; the `min` above mirrors the fact that `elem <= n-1` always.
        if elem + 1 >= n || h <= 0.0 {
            out.push(a);
        } else {
            let b = x_sorted[elem + 1];
            // NB the exact association collapse uses. Do NOT "simplify" this to
            // (1-h)*a + h*b: that changes the last bit.
            out.push(a + h * (b - a));
        }
    }
    out
}

/// `collapse::fquantile(x, probs, na.rm, type = 7)` for an unsorted input.
///
/// Sorts a copy, drops missing values, then calls [`fquantile_type7`]. Returns `None`
/// per probability if every value was missing, matching collapse's `NA` result.
pub fn fquantile_type7_unsorted(x: &[f64], probs: &[f64]) -> Vec<f64> {
    let n = x.len();
    let mut sorted: Vec<f64> = x.iter().copied().filter(|v| !v.is_nan()).collect();
    // R's `order`/`sort` places NA last and removes them under na.rm; after filtering
    // NaNs the remainder is a plain ascending sort of the finite values.
    sorted.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    fquantile_type7(&sorted, probs, n)
}

/// `triMean(x, na.rm = TRUE)`:
/// `mean(collapse::fquantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = TRUE))`.
///
/// The duplicate `0.50` is deliberate: the trimean is the mean of
/// `{Q1, Q2, Q2, Q3}` weighted by the `c(.25,.5,.5,.75)` probability vector, which
/// collapse's `probs` argument only supplies as a *selection* mechanism — the mean of
/// the four returned values is what forms the trimean.
///
/// Note that the division by 4 uses the accumulation order of R's `mean`, which
/// accumulates in `LONG_DOUBLE`. See [`r_mean`] for what that costs us.
pub fn tri_mean(x: &[f64]) -> f64 {
    const PROBS: [f64; 4] = [0.25, 0.50, 0.50, 0.75];
    let q = fquantile_type7_unsorted(x, &PROBS);
    // NB upstream's inner `mean()` takes R's **default** `na.rm = FALSE`; the `na.rm`
    // argument belongs to `fquantile`, not to `mean`. So a `NaN` among the four
    // quantiles propagates instead of being dropped. This is observable: for
    // `c(1, -Inf, 2)` the quartiles are `(NaN, 1, 1, 1.5)`, and upstream's trimean is
    // `NaN` while a `na.rm = TRUE` mean would return 1.1666...
    r_mean_no_rm(&q)
}

/// R's `mean(x, na.rm = TRUE)`.
///
/// Faithful port of `real_mean` in R 4.3.3's `src/main/summary.c:476`, which is a
/// **two-pass corrected mean** evaluated in `LDOUBLE` (x87 80-bit, 64-bit
/// significand):
///
/// ```c
/// LDOUBLE s = 0.0;
/// for (k = 0; k < n; k++) s += dx[k];        // pass 1
/// if (R_FINITE((double) s)) { s /= n; }
/// else { s = 0.; for (k) s += dx[k]/n; }
/// if (finite_s && R_FINITE((double) s)) {
///     LDOUBLE t = 0.0;
///     for (k = 0; k < n; k++) t += (dx[k] - s);  // pass 2: residuals
///     s += t/n;
/// } else if (R_FINITE((double) s)) {
///     LDOUBLE t = 0.0;
///     for (k = 0; k < n; k++) t += (dx[k] - s)/n;
///     s += t;
/// }
/// return ScalarReal((double) s);
/// ```
///
/// Three details that are easy to get wrong and each of which changes the last bit:
///
/// 1. `R_FINITE` is applied to **`(double) s`**, i.e. after narrowing, not to the
///    80-bit value. A finite `s` that overflows `f64` therefore takes the recovery
///    branch. [`F80`] has no infinities precisely so this ordering is expressible.
/// 2. `na.rm` is applied by R's `mean.default` *before* `real_mean` sees the vector,
///    so the C code never encounters a missing value. Filtering here matches that.
/// 3. The recovery branch's first form is `dx[k] / n` in **`f64`**, then widened,
///    whereas the second is `(dx[k] - s) / n` in 80-bit. They are not interchangeable.
pub fn r_mean(x: &[f64]) -> f64 {
    let vals: Vec<f64> = x.iter().copied().filter(|v| !v.is_nan()).collect();
    if vals.is_empty() {
        // R warns "no non-missing arguments to mean" and returns NaN.
        return f64::NAN;
    }
    r_mean_impl(&vals)
}

/// `real_mean` proper, on a vector that is already free of missing values.
fn r_mean_impl(vals: &[f64]) -> f64 {
    let n = vals.len();
    let nf = n as f64;

    // A non-finite *term* makes pass 1 infinite (or NaN) in 80-bit, which takes R's
    // recovery branches and then returns that value unchanged. `F80` deliberately has
    // no infinity representation (see its module docs), so mirror the same arithmetic
    // in `f64`. This is exact rather than an approximation: `Inf` and `NaN` propagate
    // identically at any precision, and the recovery branches below only ever combine
    // a non-finite accumulator with a finite `n`.
    //
    // Note this is *not* the same as the finite-input overflow case below, where the
    // 80-bit sum can exceed `f64` range and R still takes the recovery branch.
    if vals.iter().any(|v| !v.is_finite()) {
        return r_mean_non_finite(vals);
    }

    let mut s = crate::longdouble::F80::ZERO;
    for &v in vals {
        s = s.add(crate::longdouble::F80::from_f64(v));
    }
    let finite_s = s.to_f64().is_finite();

    if finite_s {
        s = s.div_int(n as u64);
    } else {
        // Overflow recovery: rescale each term in f64 before widening.
        s = crate::longdouble::F80::ZERO;
        for &v in vals {
            s = s.add(crate::longdouble::F80::from_f64(v / nf));
        }
    }

    if finite_s && s.to_f64().is_finite() {
        let mut t = crate::longdouble::F80::ZERO;
        for &v in vals {
            t = t.add(crate::longdouble::F80::from_f64(v).sub(s));
        }
        s = s.add(t.div_int(n as u64));
    } else if s.to_f64().is_finite() {
        let mut t = crate::longdouble::F80::ZERO;
        for &v in vals {
            t = t.add(crate::longdouble::F80::from_f64(v).sub(s).div_int(n as u64));
        }
        s = s.add(t);
    }

    s.to_f64()
}

/// R's `mean(x)` with `na.rm = FALSE`: `NaN` propagates, as in R.
///
/// Distinct from [`r_mean`] and both are needed — upstream `triMean` calls
/// `mean()` without `na.rm`, while `geometricMean` and `thresholdedMean` pass it
/// explicitly.
pub fn r_mean_no_rm(x: &[f64]) -> f64 {
    if x.is_empty() {
        return f64::NAN;
    }
    r_mean_impl(x)
}

/// R's `real_mean` recovery path for vectors containing `Inf`/`NaN`, in `f64`.
///
/// Split out so [`r_mean`] reads as a straight transcription of the C source.
fn r_mean_non_finite(vals: &[f64]) -> f64 {
    let nf = vals.len() as f64;
    let mut s = 0.0f64;
    for &v in vals {
        s += v;
    }
    let finite_s = s.is_finite();
    if finite_s {
        s /= nf;
    } else {
        s = 0.0;
        for &v in vals {
            s += v / nf;
        }
    }
    if finite_s && s.is_finite() {
        let mut t = 0.0f64;
        for &v in vals {
            t += v - s;
        }
        s += t / nf;
    } else if s.is_finite() {
        let mut t = 0.0f64;
        for &v in vals {
            t += (v - s) / nf;
        }
        s += t;
    }
    s
}

/// `median(x, na.rm = TRUE)`, a port of R's `median.default`
/// (`src/library/stats/R/median.R`):
///
/// ```r
/// if (na.rm) x <- x[!is.na(x)] else if (any(is.na(x))) return(NA)
/// n <- length(x)
/// if (n == 0L) return(NA)
/// half <- (n + 1L) %/% 2L
/// if (n %% 2L == 1L) sort(x, partial = half)[half]
/// else mean(sort(x, partial = half + 0L:1L)[half + 0L:1L])
/// ```
///
/// The even-`n` branch averages the two central values with [`r_mean`], i.e. R's
/// two-pass 80-bit mean — **not** a plain `(a + b) / 2`, which differs in the last bit.
pub fn r_median(x: &[f64], na_rm: bool) -> f64 {
    let mut v: Vec<f64> = x.to_vec();
    if na_rm {
        v.retain(|z| !z.is_nan());
    } else if v.iter().any(|z| z.is_nan()) {
        return f64::NAN;
    }
    let n = v.len();
    if n == 0 {
        return f64::NAN;
    }
    v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
    // R: half <- (n + 1L) %/% 2L, then 1-based index `half`.
    let half = n.div_ceil(2);
    if n % 2 == 1 {
        v[half - 1]
    } else {
        r_mean(&[v[half - 1], v[half]])
    }
}

/// `truncatedMean(x, trim, na.rm = TRUE)`: R's `mean(x, trim = trim, na.rm = TRUE)`.
///
/// Port of `mean.default` (`src/library/base/R/mean.R`):
///
/// ```r
/// if (isTRUE(na.rm)) x <- x[!is.na(x)]
/// n <- length(x)
/// if (trim > 0 && n) {
///     if (trim >= 0.5) return(stats::median(x, na.rm = FALSE))
///     lo <- floor(n*trim)+1
///     hi <- n+1-lo
///     x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
/// }
/// .Internal(mean(x))
/// ```
///
/// Note the `trim >= 0.5` short-circuit routes to `median` of the **already
/// `na.rm`-filtered** vector, and that `lo`/`hi` are 1-based inclusive.
pub fn r_trimmed_mean(x: &[f64], trim: f64, na_rm: bool) -> f64 {
    let mut v: Vec<f64> = x.to_vec();
    if na_rm {
        v.retain(|z| !z.is_nan());
    }
    let n = v.len();
    if trim > 0.0 && n > 0 {
        if trim >= 0.5 {
            // Upstream passes `na.rm = FALSE`, and NAs are already gone.
            return r_median(&v, false);
        }
        let lo = (n as f64 * trim).floor() as usize + 1; // 1-based inclusive
        let hi = n + 1 - lo;
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        v = v[lo - 1..hi].to_vec();
    }
    r_mean(&v)
}

/// R's `prod(x, na.rm)`: an `LDOUBLE` accumulation of `s = 1.0; s *= x[k]`.
///
/// Port of `rprod` in `src/main/summary.c:367`. Needed because
/// `computeExpr_coreceptor`, `computeExpr_agonist` and `computeExpr_antagonist` all end
/// in `apply(..., 2, prod)` whenever a cofactor has more than one subunit, and a plain
/// `f64` left-to-right product is **1 ulp off on the first random trial**.
///
/// R's explicit `DBL_MAX` clamp to infinity is subsumed by `to_f64()`, which also
/// yields an infinity for any finite 80-bit value beyond `f64` range.
pub fn r_prod(x: &[f64], na_rm: bool) -> f64 {
    if x.is_empty() {
        return 1.0;
    }
    let mut acc = crate::longdouble::F80::from_f64(1.0);
    for &v in x {
        if na_rm && v.is_nan() {
            continue;
        }
        acc = acc.mul(crate::longdouble::F80::from_f64(v));
    }
    acc.to_f64()
}

/// `geometricMean(x, na.rm = TRUE)`: `exp(mean(log(x), na.rm = TRUE))`.
///
/// For a **vector** input (the column case) this is
/// `exp((log x_1 + ... + log x_k) / k)` with missing values dropped.
///
/// A single zero makes the whole thing `0`, because `log 0 = -Inf` (R-ism 4 in
/// `docs/SEMANTICS.md`). The Rust `ln()` returns `-Inf` for `0.0` too, so this
/// behaviour comes for free — but it is worth an explicit test, because a naive
/// implementation that filters zeros would silently differ from upstream.
pub fn geometric_mean(x: &[f64]) -> f64 {
    if x.iter().all(|v| v.is_nan()) {
        return f64::NAN;
    }
    // `log` is applied *before* the na.rm filtering in R too, but NaN.ln() is NaN, so
    // filtering first and last are equivalent here -- and filtering first avoids
    // computing a log we would throw away.
    let logs: Vec<f64> = x
        .iter()
        .copied()
        .filter(|v| !v.is_nan())
        .map(|v| v.ln())
        .collect();
    r_mean(&logs).exp()
}

/// `thresholdedMean(x, trim, na.rm = TRUE)`:
///
/// ```r
/// percent <- Matrix::nnzero(x)/length(x)
/// if (percent < trim) return(0) else return(mean(x, na.rm = na.rm))
/// ```
///
/// `Matrix::nnzero(x)`, as R actually computes it.
///
/// It does **not** "count the non-zero entries and skip the missing ones". Measured against
/// `Matrix::nnzero`:
///
/// ```text
/// all_nan   nnzero=NA     all_zero  nnzero=0     one_nan  nnzero=NA
/// one_na    nnzero=NA     all_pos   nnzero=15
/// ```
///
/// A single missing entry makes the whole answer `NA`, not a smaller count. Upstream's
/// `thresholdedMean` then evaluates `if (percent < trim)` with `percent == NA`, and R raises
/// `missing value where TRUE/FALSE needed` *inside the mean*, before any `Prob` exists.
pub fn r_nnzero(x: &[f64]) -> f64 {
    let mut n = 0.0f64;
    for &v in x {
        if v.is_nan() {
            return f64::NAN;
        }
        if v != 0.0 {
            n += 1.0;
        }
    }
    n
}

/// `thresholdedMean(x, trim, na.rm = TRUE)`:
///
/// ```r
/// percent <- Matrix::nnzero(x)/length(x)
/// if (percent < trim) return(0) else return(mean(x, na.rm = na.rm))
/// ```
///
/// Note the *asymmetry* with the other means: the threshold is evaluated on the raw
/// vector including missing slots, but the `trim` comparison is against the fraction of
/// non-zero entries — so a vector of all zeros has `nnzero == 0` and returns `0`.
///
/// When `percent` is `NA` this returns `NaN` rather than choosing a branch. R cannot
/// evaluate `NA < trim` either, and it raises there. Letting the `NaN` travel into `P1` puts it
/// in `sum(P1 * P.spatial)`, where the kernel's `if (NA)` check raises upstream's
/// `missing value where TRUE/FALSE needed` with the identical message. Choosing `0` here instead
/// is what let the port return an all-zero `Prob` where upstream refuses to return anything, and it
/// is invisible to every fixture without a missing value.
pub fn thresholded_mean(x: &[f64], trim: f64) -> f64 {
    let percent = r_nnzero(x) / x.len() as f64;
    if percent.is_nan() {
        return f64::NAN;
    }
    if percent < trim {
        return 0.0;
    }
    r_mean(x)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Measured against `Matrix::nnzero` on the same vectors; the port returned a count for
    /// every one of these and upstream returns `NA` for three of the five.
    #[test]
    fn nnzero_is_na_when_any_entry_is_missing() {
        let all_nan = vec![f64::NAN; 15];
        let all_zero = vec![0.0; 15];
        let mut one_nan = vec![0.0; 15];
        one_nan[3] = f64::NAN;
        let mut one_na = vec![0.0; 15];
        one_na[3] = f64::NAN;
        assert!(r_nnzero(&all_nan).is_nan());
        assert_eq!(r_nnzero(&all_zero), 0.0);
        assert!(r_nnzero(&one_nan).is_nan());
        assert!(r_nnzero(&one_na).is_nan());
        assert_eq!(
            r_nnzero(&(1..=15).map(|i| i as f64).collect::<Vec<_>>()),
            15.0
        );
    }

    /// `thresholdedMean` must not return `0` for a vector with a missing entry: upstream's
    /// `if (percent < trim)` raises, and a silent `0` turns a hard error into an all-zero `Prob`.
    #[test]
    fn thresholded_mean_propagates_nnzero_na() {
        let all_nan = vec![f64::NAN; 15];
        assert!(thresholded_mean(&all_nan, 0.1).is_nan());
        let mut one_nan = vec![0.0; 15];
        one_nan[3] = f64::NAN;
        assert!(thresholded_mean(&one_nan, 0.1).is_nan());
        // A vector with no missing value keeps the threshold behaviour.
        assert_eq!(thresholded_mean(&[0.0; 15], 0.1), 0.0);
        assert_eq!(
            thresholded_mean(&(1..=15).map(|i| i as f64).collect::<Vec<_>>(), 0.1),
            8.0
        );
    }

    #[test]
    fn type7_boundary_probs_take_endpoints() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let q = fquantile_type7(&x, &[0.0, 0.25, 0.5, 0.75, 1.0], 4);
        assert_eq!(q, vec![1.0, 1.75, 2.5, 3.25, 4.0]);
    }

    #[test]
    fn type7_single_observation() {
        let x = [42.0];
        let q = fquantile_type7(&x, &[0.25, 0.5, 0.75], 1);
        assert_eq!(q, vec![42.0; 3]);
    }

    #[test]
    fn type7_empty_input_is_nan() {
        let q = fquantile_type7(&[], &[0.25, 0.5], 0);
        assert!(q.iter().all(|v| v.is_nan()));
    }

    #[test]
    fn type7_drops_nan() {
        let x = [1.0, f64::NAN, 3.0, f64::NAN, 5.0];
        let q = fquantile_type7_unsorted(&x, &[0.0, 0.5, 1.0]);
        assert_eq!(q, vec![1.0, 3.0, 5.0]);
    }

    #[test]
    fn tri_mean_of_uniform_equals_the_mean() {
        // For a symmetric distribution the trimean coincides with the mean.
        let x: Vec<f64> = (1..=101).map(|i| i as f64).collect();
        assert!((tri_mean(&x) - 51.0).abs() < 1e-12);
    }

    #[test]
    fn tri_mean_is_zero_when_three_quarters_are_zero() {
        // The L1 zero short-circuit in the aggregation relies on this being exactly 0.0.
        let x = [0.0; 8];
        let mut v = x.to_vec();
        v.push(1.0);
        v.push(2.0);
        assert_eq!(tri_mean(&x), 0.0);
        assert_eq!(tri_mean(&v), 0.0);
    }

    #[test]
    fn geometric_mean_of_a_zero_is_zero() {
        assert_eq!(geometric_mean(&[1.0, 0.0, 4.0]), 0.0);
    }

    #[test]
    fn geometric_mean_is_the_log_space_arithmetic_mean() {
        assert!((geometric_mean(&[1.0, 2.0, 3.0, 4.0]) - 24f64.powf(0.25)).abs() < 1e-12);
    }

    #[test]
    fn geometric_mean_drops_nan_but_not_zero() {
        let with_nan = geometric_mean(&[1.0, f64::NAN, 4.0]);
        let without = geometric_mean(&[1.0, 4.0]);
        assert_eq!(with_nan, without);
    }

    #[test]
    fn thresholded_mean_returns_zero_below_threshold() {
        // 1 non-zero out of 100 = 1% < trim = 10%.
        let mut x = vec![0.0; 100];
        x[0] = 5.0;
        assert_eq!(thresholded_mean(&x, 0.10), 0.0);
        // At exactly 10% it must take the mean branch.
        let mut y = vec![0.0; 10];
        y[0] = 5.0;
        assert!((thresholded_mean(&y, 0.10) - 0.5).abs() < 1e-12);
    }

    /// This test used to assert that `c(NaN, 1, 3, 5)` averages to `3.0`, i.e. that a missing
    /// entry is skipped and the mean branch still runs. It does not: upstream **raises**.
    ///
    /// ```text
    /// c(NaN,1,3,5) trim=0.5 : ERROR: missing value where TRUE/FALSE needed
    /// c(0,0,0,0)   trim=0.1 : 0
    /// c(1,3,5,7)   trim=0.5 : 4
    /// ```
    ///
    /// `Matrix::nnzero` is `NA` for any vector containing a missing value, so `percent` is `NA`
    /// and R's `if` cannot choose a branch. The old assertion was a plausible reading of the
    /// function -- it is the `na.rm = TRUE` half of the signature taken in isolation -- and it is
    /// precisely why the port returned an all-zero `Prob` on six matrix configurations.
    #[test]
    fn thresholded_mean_raises_rather_than_skipping_a_missing_value() {
        let x = [f64::NAN, 1.0, 3.0, 5.0];
        assert!(thresholded_mean(&x, 0.5).is_nan());
        assert!(thresholded_mean(&x, 0.9).is_nan());
        // No missing value: the threshold branch is decided by the non-zero count alone.
        assert_eq!(thresholded_mean(&[1.0, 3.0, 5.0, 7.0], 0.5), 4.0);
        assert_eq!(thresholded_mean(&[0.0; 4], 0.1), 0.0);
    }

    /// The radix-order fast path is **unreachable for `triMean`**, which is the only
    /// caller. `collapse` switches when
    /// `length(x) > 1e5 && length(probs) > log(length(x))`; with `length(probs) == 4`
    /// the second clause needs `length(x) < e^4 ~= 54.6`, which contradicts the first.
    ///
    /// Worth pinning as a test: if a future `collapse` release changes the condition,
    /// the assumption that the direct path is always taken would silently rot.
    #[test]
    fn radix_order_threshold_matches_collapse() {
        assert!(!collapse_uses_radix_order(100_000, 4));
        assert!(!collapse_uses_radix_order(1_000_000, 4));
        for len in [1usize, 2, 10, 1000, 100_000, 1_000_000] {
            assert!(
                !collapse_uses_radix_order(len, 4),
                "triMean's 4 probs must never take the radix path (len={len})"
            );
        }
        // The path is reachable with many probs, so the predicate is not vacuous:
        // log(200000) = 12.2, so 13 probs qualifies.
        assert!(collapse_uses_radix_order(200_000, 13));
        assert!(!collapse_uses_radix_order(200_000, 12));
    }
}
