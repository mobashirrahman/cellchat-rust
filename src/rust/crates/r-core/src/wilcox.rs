//! `rank()`, `wilcox.test`'s two-sample branch, `p.adjust`, and the pieces of
//! `identifyOverExpressedGenes` that surround them.
//!
//! Parity contract: **Exact** (`docs/SEMANTICS.md`).
//!
//! Upstream: `R/utilities.R:376` (`identifyOverExpressedGenes`, `do.fast = FALSE` path),
//! `stats::wilcox.test.default`, `stats::p.adjust`, `stats::rank`.
//!
//! ## Why Wilcoxon is in scope at all
//!
//! `identifyOverExpressedGenes` is on the critical path for `createCellChat` and is the
//! second-largest term in the Amdahl decomposition the objective asks for. With
//! `do.fast = FALSE` — the setting the objective's correctness claims rest on, because
//! `do.fast = TRUE` routes to `presto::wilcoxauc`, a *different* test with different logFC
//! values — it calls `stats::wilcox.test` once per gene per cell group.
//!
//! ## The two branches
//!
//! `wilcox.test`'s default is `exact = NULL`, which resolves to
//! `(n.x < 50) && (n.y < 50)`, and the exact branch additionally requires **no ties**:
//!
//! ```r
//! if (exact && !TIES) { ... exact ... } else { ... normal approximation ... }
//! ```
//!
//! Real single-cell data has hundreds to thousands of cells per group and integer counts
//! with many zeros, so ties are essentially always present and the normal approximation
//! is what actually runs. The exact branch is nevertheless implemented, because a
//! bit-identical port cannot be one that silently takes a different branch on small
//! input.
//!
//! ## Continuity correction
//!
//! `correct = TRUE` by default, and for `alternative = "two.sided"` the correction is
//! `sign(z) * 0.5` where `z = W - n.x * n.y / 2`. For `z == 0` exactly, `sign(0) == 0`, so
//! **no** correction is applied — a discontinuity that a fixture with a perfectly balanced
//! split will hit.

use crate::mathfn::{choose, pnorm};
use std::collections::HashMap;

/// `rank(x)` with `ties.method = "average"` (the default) and `na.last = "keep"`.
///
/// `digits.rank = Inf`, so no `signif` is applied.
///
/// C99/R's `sign(x)`, which returns `+/-0.0` for zero.
///
/// **Not** Rust's `f64::signum`, which is `1.0f64.copysign(self)` and therefore `+1.0` for
/// `+0.0`. Using it makes `sign(z) * 0.5 == 0.5` at `z == 0`, so a perfectly balanced split
/// gets a *half* correction it should not have, and the p-value moves in the middle of the
/// range. This is the discontinuity the module notes call out, and it caused a real corpus
/// failure: an all-ties gene has `SIGMA == 0`, so `z = 0/0` must be `NaN` (R returns `NaN`),
/// but with the spurious `-0.5` numerator it became `-inf` and the p-value `0` -- which
/// would have reported every all-ties gene as maximally significant.
#[inline]
fn r_sign(x: f64) -> f64 {
    if x > 0.0 {
        1.0
    } else if x < 0.0 {
        -1.0
    } else {
        // Preserves -0.0, like C's signbit-aware sign().
        x
    }
}

/// ## `na.last = "keep"` means NA is *ranked*, not that the result is NA
///
/// This is the part that surprises, and it is R's actual behaviour:
///
/// ```r
/// rank(c(1, NA, 2))   # 1 3 2
/// rank(c(NA, NA, 1))  # 2 3 1
/// ```
///
/// NA sorts *last* and takes the highest ranks, in input order among themselves. Returning
/// `NA` for those positions -- which is the obvious reading of the name -- gives
/// `c(1, NA, 3)` and is wrong for every vector with a missing value.
///
/// NB this does not affect `wilcox.test`, which drops NAs before ranking.
pub fn rank_average(x: &[f64]) -> Vec<f64> {
    let n = x.len();
    // `order(..., na.last = TRUE)` -- non-NA ascending, then the NAs in input order.
    let mut idx: Vec<usize> = (0..n).filter(|&i| !x[i].is_nan()).collect();
    idx.sort_by(|&a, &b| x[a].partial_cmp(&x[b]).unwrap_or(std::cmp::Ordering::Equal));
    let mut out = vec![0.0f64; n];
    let mut i = 0;
    while i < idx.len() {
        let mut j = i + 1;
        while j < idx.len() && x[idx[j]] == x[idx[i]] {
            j += 1;
        }
        // 1-based ranks i+1 ..= j averaged.
        let avg = ((i + 1 + j) as f64) / 2.0;
        for &k in &idx[i..j] {
            out[k] = avg;
        }
        i = j;
    }
    // The NAs take the remaining ranks, in input order.
    let mut r = idx.len() as f64 + 1.0;
    for (k, &v) in x.iter().enumerate() {
        if v.is_nan() {
            out[k] = r;
            r += 1.0;
        }
    }
    out
}

/// The two-sample Wilcoxon rank-sum statistic and the p-value, for R's defaults
/// (`alternative = "two.sided"`, `mu = 0`, `exact = NULL`, `correct = TRUE`).
///
/// Returns `(W, p, exact_used)`.
pub fn wilcox_test_two_sided(x: &[f64], y: &[f64]) -> (f64, f64, bool) {
    wilcox_test_two_sided_correct(x, y, true)
}

/// [`wilcox_test_two_sided`] with R's `correct` argument.
///
/// CellChat always uses the default `correct = TRUE`, but `stats::wilcox.test` exposes the
/// argument, so the port has to reproduce it to keep the two results distinguishable: for
/// case 2 of the corpus, `W = 8` gives `p = 0.3976...` with the correction and
/// `p = 0.3412...` without. Reading a `correct = FALSE` golden record and comparing it
/// against a `correct = TRUE` implementation is a silent 1-in-13 chance of passing
/// (n.x * n.y even) and otherwise an obvious mismatch.
pub fn wilcox_test_two_sided_correct(x: &[f64], y: &[f64], correct: bool) -> (f64, f64, bool) {
    // R: `y <- y[!is.na(y)]` happens inside the `!is.null(y)` branch, and
    // `x <- x[!is.na(x)]` after. Both are dropped, not imputed.
    let xv: Vec<f64> = x.iter().copied().filter(|v| !v.is_nan()).collect();
    let yv: Vec<f64> = y.iter().copied().filter(|v| !v.is_nan()).collect();
    // R stops with "not enough 'y' observations" if y is empty after NA removal, and
    // "not enough (non-missing) 'x' observations" if x is empty. CellChat guarantees both
    // sides are non-empty (each group and its complement), so the ordering of the two
    // checks is not observable here; they are documented rather than modelled.
    let n_x = xv.len() as f64;
    let n_y = yv.len() as f64;
    assert!(
        !xv.is_empty() && !yv.is_empty(),
        "wilcox_test: empty sample"
    );

    // `r <- rank(c(x - mu, y))` with mu = 0.
    let mut all = Vec::with_capacity(xv.len() + yv.len());
    all.extend_from_slice(&xv);
    all.extend_from_slice(&yv);
    let r = rank_average(&all);
    // `STATISTIC <- c(W = sum(r[seq_along(x)]) - n.x * (n.x + 1) / 2)`
    //
    // NB `r[seq_along(x)]` is x's *ranks within the pooled sample*, and the subtrahend is
    // n.x*(n.x+1)/2 — not n.x*n.y/2. The n.x*n.y/2 appears later, in the null mean.
    let w: f64 = r.iter().take(xv.len()).sum::<f64>() - n_x * (n_x + 1.0) / 2.0;

    // `TIES <- (length(r) != length(unique(r)))`
    let mut seen: HashMap<u64, ()> = HashMap::new();
    let mut n_distinct = 0;
    for v in &r {
        if seen.insert(v.to_bits(), ()).is_none() {
            n_distinct += 1;
        }
    }
    let ties = n_distinct != r.len();

    let exact = (n_x < 50.0) && (n_y < 50.0);
    if exact && !ties {
        // p <- if (W > n.x*n.y/2) pwilcox(W - 1, lower=FALSE) else pwilcox(W)
        let p = if w > n_x * n_y / 2.0 {
            pwilcox(w - 1.0, n_x, n_y, false)
        } else {
            pwilcox(w, n_x, n_y, true)
        };
        return (w, (2.0 * p).min(1.0), true);
    }

    // Normal approximation with the tie correction and the continuity correction.
    //
    //   NTIES <- table(r)
    //   z <- STATISTIC - n.x * n.y / 2
    //   SIGMA <- sqrt((n.x * n.y / 12) * ((n.x + n.y + 1) -
    //                   sum(NTIES^3 - NTIES) / ((n.x + n.y) * (n.x + n.y - 1))))
    //   if (correct) CORRECTION <- sign(z) * 0.5
    //   z <- (z - CORRECTION) / SIGMA
    //   PVAL <- 2 * min(pnorm(z), pnorm(z, lower.tail = FALSE))
    //
    // `sum(NTIES^3 - NTIES)` is an *integer* sum in R (`table()` returns integers, and `^`
    // promotes to double). Computed as f64 from integer cubes, which is exact up to 2^53.
    let mut tie_sum = 0.0f64;
    {
        let mut counts: HashMap<u64, f64> = HashMap::new();
        for v in &r {
            *counts.entry(v.to_bits()).or_insert(0.0) += 1.0;
        }
        for c in counts.values() {
            let n = *c;
            tie_sum += n * n * n - n;
        }
    }
    let nn = n_x + n_y;
    let mut z = w - n_x * n_y / 2.0;
    let sigma = ((n_x * n_y / 12.0) * (nn + 1.0 - tie_sum / (nn * (nn - 1.0)))).sqrt();
    // `alternative = "two.sided"` -> `p <- 2 * min(pnorm(z), pnorm(z, lower.tail = FALSE))`,
    // and with `correct = TRUE` R first does `z <- (z - sign(z) * 0.5)/SIGMA`. R's
    // `sign(0)` is 0, so a perfectly balanced split gets *no* correction -- see [`r_sign`].
    if correct {
        z = (z - r_sign(z) * 0.5) / sigma;
    } else {
        z /= sigma;
    }
    let p = 2.0 * pnorm(z, true).min(pnorm(z, false));
    (w, p, false)
}

/// `cwilcox(k, m, n)` from `src/nmath/wilcox.c`: the number of arrangements of the
/// Mann-Whitney statistic equal to `k`.
///
/// Transcribed including both of R's reductions, which are load-bearing for the
/// *arithmetic* and therefore for the last bits:
///
/// ```c
/// u = m*n; if (k < 0 || k > u) return 0;
/// c = u/2; if (k > c) k = u - k;          /* symmetry */
/// if (m < n) {i = m; j = n;} else {i = n; j = m;}   /* i <= j */
/// if (j == 0) return (k == 0);
/// if (j > 0 && k < j) return cwilcox(k, i, k);      /* small-k reduction */
/// return cwilcox(k - j, i - 1, j) + cwilcox(k, i, j - 1);
/// ```
///
/// Memoised on `(i, j, k)`, as R is. The counts are `double` in R (so `choose` in
/// `pwilcox` divides a double by a double); the values here are integers representable
/// exactly in `f64` for every `m, n < 50`.
fn cwilcox(k: i64, m: i64, n: i64, memo: &mut HashMap<(i64, i64, i64), f64>) -> f64 {
    let u = m * n;
    if k < 0 || k > u {
        return 0.0;
    }
    let c = u / 2;
    let mut k = k;
    if k > c {
        k = u - k;
    }
    let (i, j) = if m < n { (m, n) } else { (n, m) };
    if j == 0 {
        return if k == 0 { 1.0 } else { 0.0 };
    }
    if k < j {
        return cwilcox(k, i, k, memo);
    }
    if let Some(&v) = memo.get(&(i, j, k)) {
        return v;
    }
    let v = cwilcox(k - j, i - 1, j, memo) + cwilcox(k, i, j - 1, memo);
    memo.insert((i, j, k), v);
    v
}

/// `pwilcox(q, m, n, lower_tail)` from `src/nmath/wilcox.c`.
///
/// ```c
/// q = floor(q + 1e-7);
/// if (q < 0) return 0; if (q >= m*n) return 1;
/// c = choose(m + n, n); p = 0;
/// if (q <= m*n/2) { for (i = 0; i <= q; i++) p += cwilcox(i, mm, nn) / c; }
/// else { q = m*n - q; for (i = 0; i < q; i++) p += cwilcox(i, mm, nn) / c;
///        lower_tail = !lower_tail; }
/// return lower_tail ? p : 1 - p;
/// ```
///
/// The summation order and the `1 - p` flip are R's. Note that the upper tail is
/// `1 - p`, computed here exactly as R does — unlike `pnorm`, where R computes the tail
/// directly.
pub fn pwilcox(q: f64, m: f64, n: f64, lower_tail: bool) -> f64 {
    if q.is_nan() || m.is_nan() || n.is_nan() {
        return q + m + n;
    }
    if !m.is_finite() || !n.is_finite() {
        return f64::NAN;
    }
    let (m, n) = (m as i64 as f64, n as i64 as f64);
    if m <= 0.0 || n <= 0.0 {
        return f64::NAN;
    }
    let q = (q + 1e-7).floor();
    // R returns `R_DT_0` / `R_DT_1`, which are `lower_tail ? 0 : 1` and its complement.
    // Returning the value unconditionally is a tail-direction bug that only shows up
    // outside the support, i.e. at q = -1 and q >= m*n -- both in the fixture.
    if q < 0.0 {
        return if lower_tail { 0.0 } else { 1.0 };
    }
    if q >= m * n {
        return if lower_tail { 1.0 } else { 0.0 };
    }
    let (mm, nn) = (m as i64, n as i64);
    let mut memo: HashMap<(i64, i64, i64), f64> = HashMap::new();
    let c = choose(m + n, n);
    let mut p = 0.0f64;
    let mut lower = lower_tail;
    let qi = q as i64;
    if q <= m * n / 2.0 {
        for i in 0..=qi {
            p += cwilcox(i, mm, nn, &mut memo) / c;
        }
    } else {
        let q2 = (m * n) as i64 - qi;
        for i in 0..q2 {
            p += cwilcox(i, mm, nn, &mut memo) / c;
        }
        lower = !lower;
    }
    if lower {
        p
    } else {
        1.0 - p
    }
}

/// `stats::p.adjust(p, method = "bonferroni")`.
///
/// ```r
/// p.adjust.bonferroni <- function(p, n = length(p)) pmin(1.0, n * p)
/// ```
///
/// Two details: the default `n` is `length(p)` — the number of *genes in the matrix*, not
/// the number of surviving features — and `n * p` is evaluated in that order, so a
/// commutative rearrangement is not the same expression.
pub fn p_adjust_bonferroni(p: &[f64], n: f64) -> Vec<f64> {
    p.iter().map(|&v| (n * v).min(1.0)).collect()
}

/// `mean.fxn <- function(x) log(mean(expm1(x)) + 1)` from `identifyOverExpressedGenes`.
///
/// `mean` is R's two-pass corrected 80-bit mean ([`crate::stats::r_mean_no_rm`]),
/// `expm1` and `log` are libm, which R also calls. `apply(X, 1, mean.fxn)` runs over the
/// genes of `data.use[features, cells]`, i.e. over cells within a group.
pub fn mean_fxn(x: &[f64]) -> f64 {
    let expm1: Vec<f64> = x.iter().map(|v| v.exp_m1()).collect();
    (crate::stats::r_mean_no_rm(&expm1) + 1.0).ln()
}

/// One row of `identifyOverExpressedGenes`'s output.
#[derive(Clone, Debug, PartialEq)]
pub struct Marker {
    /// `level.use[i]`, resolved by the caller from `group_index`; empty in
    /// [`identify_over_expressed_one_group`], which does not know the level names.
    pub cluster: String,
    /// Index into the group levels, set by [`identify_over_expressed`].
    pub group_index: usize,
    pub feature: String,
    pub pvalue: f64,
    /// `FC[feature]` — the difference of `mean.fxn` between the group and the rest.
    pub log_fc: f64,
    /// `pct.1` and `pct.2`, each `round(..., 3)`.
    pub pct_1: f64,
    pub pct_2: f64,
    /// `pvalues.adj` — Bonferroni over the whole matrix, `n = nrow(X)`.
    pub pvalue_adj: f64,
}

/// `identifyOverExpressedGenes(..., do.fast = FALSE)`'s single-group inner body.
///
/// `data` is genes x cells column-major; `features` is the full gene universe in
/// `rownames` order; `in_group` / `in_rest` are the cell indices for `level.use[i]` and
/// everything else.
///
/// Reproduces the `do.fast = FALSE` branch of the pinned commit exactly, and specifically
/// **not** the `do.fast = TRUE` (presto) branch. The two are easy to conflate because the
/// presto branch's filters read as if they applied to both:
///
/// ```r
/// markers.all <- dplyr::filter(genes.de, pvalues < thresh.p,
///                              logFC_abs >= thresh.fc, pct.max > thresh.pc*100)
/// ```
///
/// That `pct.max > thresh.pc*100` compares a *percentage* (0-100) against `thresh.pc*100`
/// while the selection step compares a *fraction* (0-1) against `thresh.pc` -- two
/// different scales. The Wilcoxon branch has no such filter at all; it selects with
/// `alpha.min > thresh.pc` and `FC > thresh.fc`, and then filters only on
/// `pvalues < thresh.p` and, when `only.pos`, `logFC > 0`. Porting the presto filter
/// wholesale into the Wilcoxon path silently drops every feature at the default
/// `thresh.pc = 0`... no: it does not, at 0 -- but it does for any `thresh.pc > 0`, and
/// `identifyOverExpressedLigandsReceptors` is routinely called with `thresh.pc = 0.01`.
pub fn identify_over_expressed_one_group(
    data: &[f64],
    features: &[String],
    in_group: &[usize],
    in_rest: &[usize],
    thresh_pc: f64,
    thresh_fc: f64,
    thresh_p: f64,
    only_pos: bool,
) -> Vec<Marker> {
    let n_adjust = features.len();
    identify_over_expressed_one_group_n(
        data, features, in_group, in_rest, thresh_pc, thresh_fc, thresh_p, only_pos, n_adjust,
    )
    .rows
}

/// The per-group result plus the count the R shim needs to reproduce upstream's *schema*.
///
/// Upstream's aggregation is
/// ```r
/// markers.all <- data.frame()
/// for (i in 1:numCluster) {
///   gde <- genes.de[[i]]
///   if (!is.null(gde)) {
///     gde <- gde[order(gde$pvalues, -gde$logFC), ]
///     gde <- subset(gde, subset = pvalues < thresh.p)
///     if (nrow(gde) > 0) markers.all <- rbind(markers.all, gde)
///   }
/// }
/// ```
/// so when **no** group has a row clearing `thresh.p`, `markers.all` is never rbind-ed and
/// stays the zero-column `data.frame()`. The following
/// `markers.all$features <- as.character(markers.all$features)` then *adds* a single
/// zero-length column, and the "marker table" is `0x1` with just `features` rather than
/// `0x7`. `identical()` on the S4 object sees the difference, so the port has to know which
/// of the two shapes it is: this counter is the only place the information still exists.
/// A "just return an empty Vec<Vec>" API cannot express it.
pub struct GroupMarkers {
    pub rows: Vec<Marker>,
    /// Rows that passed `pvalues < thresh.p`, *before* the `only.pos` filter and before
    /// the cross-group `rbind`. `> 0` implies the full seven-column schema.
    pub n_before_only_pos: usize,
}

/// [`identify_over_expressed_one_group`] with the Bonferroni multiplier decoupled from the
/// feature count.
///
/// `p.adjust(pvalues, "bonferroni", n = nrow(X))` uses `nrow(X)` -- the row count of
/// `object@data.signaling` -- **not** `nrow(data.use)`. `data.use <- X[features.use, ]` can
/// be a strict subset whenever the caller passes `features`, and the multiplier then
/// diverges. With the default `nrow(X)` this is invisible; with a `features` argument the
/// adjusted p-values are smaller by exactly the subset ratio, which changes which features
/// clear `thresh.p`.
#[allow(clippy::too_many_arguments)]
pub fn identify_over_expressed_one_group_n(
    data: &[f64],
    features: &[String],
    in_group: &[usize],
    in_rest: &[usize],
    thresh_pc: f64,
    thresh_fc: f64,
    thresh_p: f64,
    only_pos: bool,
    n_adjust: usize,
) -> GroupMarkers {
    let n_genes = features.len();
    let n_cells = data.len() / n_genes;
    assert_eq!(data.len(), n_genes * n_cells, "bad data buffer");
    assert!(!in_group.is_empty(), "the group has no cells");
    // `cell.use2 <- base::setdiff(1:length(labels), cell.use1)` -- every cell not in the
    // group, which is never empty when there is more than one group.
    assert!(!in_rest.is_empty(), "the complement has no cells");

    let col = |g: usize, c: usize| data[c * n_genes + g];

    // Feature selection by percentage:
    //   pct.1 <- round(rowSums(data.use[features, cell.use1] > 0) / length(cell.use1), 3)
    //   alpha.min <- apply(cbind(pct.1, pct.2), 1, max)
    //   features <- names(which(alpha.min > thresh.pc))
    let pct_of = |cells: &[usize]| -> Vec<f64> {
        (0..n_genes)
            .map(|g| {
                let n_pos = cells.iter().filter(|&&c| col(g, c) > 0.0).count();
                round3(n_pos as f64 / cells.len() as f64)
            })
            .collect()
    };
    let pct_1 = pct_of(in_group);
    let pct_2 = pct_of(in_rest);

    // Feature selection by average difference:
    //   data.1 <- apply(data.use[features, cell.use1], 1, mean.fxn)
    //   FC <- data.1 - data.2
    //   features.diff <- names(which(FC > thresh.fc))          [only.pos]
    //   features <- intersect(features, features.diff)
    let row_mean_fxn = |cells: &[usize], g: usize| -> f64 {
        let v: Vec<f64> = cells.iter().map(|&c| col(g, c)).collect();
        mean_fxn(&v)
    };
    let mut out = Vec::new();
    // Upstream computes the two vectors over the *pct-selected* subset and then
    // `intersect`s. R's `intersect(x, y)` is `unique(y[match(x, y, 0L)])`: the surviving
    // order follows `x`, so it is the pct order, and both arguments are already ascending
    // in feature order. `pct_pass`/`diff_pass` reproduce that.
    let pct_pass: Vec<usize> = (0..n_genes)
        .filter(|&g| pct_1[g].max(pct_2[g]) > thresh_pc)
        .collect();
    if pct_pass.is_empty() {
        return GroupMarkers {
            rows: out,
            n_before_only_pos: 0,
        };
    }
    // A single surviving feature makes the group contribute **nothing**, and that is upstream's
    // behaviour rather than an accident of this port. `apply(X, 1, FUN)` on a one-row matrix
    // returns an unnamed *scalar* instead of a named vector of length one:
    //
    // ```r
    // a <- apply(matrix(1:4, nrow = 1), 1, function(x) x[1] + x[2])   # 3, no names
    // FC <- a - a2
    // FC["gene1"]                                                     # NA
    // names(which(FC > 0))                                            # character(0)
    // ```
    //
    // so `features.diff <- names(which(FC > thresh.fc))` is empty, `intersect(features,
    // features.diff)` is empty, and the loop's `next` skips the group. The port, which keeps the
    // feature names attached to the values, would keep the feature and report a marker -- so the
    // two disagree on exactly the inputs where a group has a single marker candidate.
    if pct_pass.len() == 1 {
        return GroupMarkers {
            rows: out,
            n_before_only_pos: 0,
        };
    }
    let m1: Vec<f64> = pct_pass
        .iter()
        .map(|&g| row_mean_fxn(in_group, g))
        .collect();
    let m2: Vec<f64> = pct_pass.iter().map(|&g| row_mean_fxn(in_rest, g)).collect();
    let fc: Vec<f64> = m1.iter().zip(&m2).map(|(a, b)| a - b).collect();
    let diff_pass: Vec<usize> = (0..pct_pass.len())
        .filter(|&i| {
            if only_pos {
                fc[i] > thresh_fc
            } else {
                fc[i].abs() > thresh_fc
            }
        })
        .collect();
    if diff_pass.is_empty() {
        return GroupMarkers {
            rows: out,
            n_before_only_pos: 0,
        };
    }

    // `pvalues <- unlist(my.sapply(1:nrow(data1), function(x) wilcox.test(data1[x,],
    // data2[x,])$p.value))`. `wilcox.test`'s defaults are `alternative = "two.sided"` and
    // `correct = TRUE`; the commented-out `alternative = "greater"` line above it is dead.
    let pvalues: Vec<f64> = diff_pass
        .iter()
        .map(|&i| {
            let g = pct_pass[i];
            let a: Vec<f64> = in_group.iter().map(|&c| col(g, c)).collect();
            let b: Vec<f64> = in_rest.iter().map(|&c| col(g, c)).collect();
            wilcox_test_two_sided(&a, &b).1
        })
        .collect();
    // `pval.adj = stats::p.adjust(pvalues, method = "bonferroni", n = nrow(X))` -- nrow of
    // the *whole* matrix, not of the selected features, so the multiplier is the full
    // number of genes even when only a handful survive selection.
    let adj = p_adjust_bonferroni(&pvalues, n_adjust as f64);

    // genes.de[[i]] <- data.frame(clusters, features, pvalues, logFC, pct.1, pct.2,
    //                             pvalues.adj)
    let rows: Vec<Marker> = diff_pass
        .iter()
        .enumerate()
        .map(|(k, &i)| {
            let g = pct_pass[i];
            Marker {
                cluster: String::new(),
                group_index: 0,
                feature: features[g].clone(),
                pvalue: pvalues[k],
                log_fc: fc[i],
                pct_1: pct_1[g],
                pct_2: pct_2[g],
                pvalue_adj: adj[k],
            }
        })
        .collect();

    // The tail of the function, which is the part that decides the *row order* of the
    // marker table and therefore of `object@var.features[["features"]]` and the
    // `de_foldChange` columns:
    //
    //   markers.all <- data.frame()
    //   for (i in 1:numCluster) {
    //     gde <- genes.de[[i]]
    //     if (!is.null(gde)) {
    //       gde <- gde[order(gde$pvalues, -gde$logFC), ]
    //       gde <- subset(gde, subset = pvalues < thresh.p)
    //       if (nrow(gde) > 0) markers.all <- rbind(markers.all, gde)
    //     }
    //   }
    //   if (only.pos & nrow(markers.all) > 0) markers.all <- subset(markers.all, logFC > 0)
    //
    // So the order is **per group, then by p-value ascending, then by logFC descending** --
    // not feature order. Emitting `diff_pass` in feature order agrees on the *set* of
    // markers and disagrees on the *order*, which changes `features.sig` and hence
    // `identifyOverExpressedLigandsReceptors`'s LR table.
    let mut sorted: Vec<Marker> = rows
        .into_iter()
        .filter(|m| {
            // `subset(gde, pvalues < thresh.p)` -- `NA`/`NaN` comparisons are dropped.
            m.pvalue < thresh_p
        })
        .collect();
    // `order(pvalues, -logFC)`. R's `order` with `method = "auto"` is a stable radix sort
    // for numeric input, so features with equal (pvalue, -logFC) keep feature order.
    sorted.sort_by(|a, b| {
        a.pvalue
            .total_cmp(&b.pvalue)
            .then_with(|| (-a.log_fc).total_cmp(&(-b.log_fc)))
    });
    let n_before_only_pos = sorted.len();
    if only_pos {
        sorted.retain(|m| m.log_fc > 0.0);
    }
    out.extend(sorted);
    GroupMarkers {
        rows: out,
        n_before_only_pos,
    }
}

/// `identifyOverExpressedGenes(..., do.fast = FALSE)`: the per-group loop plus `rbind`.
///
/// `labels` is the per-cell group index in *factor level* order, so the group iteration
/// order is the factor's level order, not the order of first appearance.
#[allow(clippy::too_many_arguments)]
pub fn identify_over_expressed(
    data: &[f64],
    features: &[String],
    labels: &[usize],
    n_groups: usize,
    thresh_pc: f64,
    thresh_fc: f64,
    thresh_p: f64,
    only_pos: bool,
    n_adjust: usize,
) -> IdentifyResult {
    let n_cells = data.len() / features.len();
    assert_eq!(labels.len(), n_cells, "one label per cell");
    let mut select = |level: usize| -> Option<(Vec<usize>, Vec<usize>)> {
        let in_group: Vec<usize> = (0..n_cells).filter(|&c| labels[c] == level).collect();
        // `level.use <- levels(labels)[levels(labels) %in% unique(labels)]` drops a
        // level with no cells before the loop, so `genes.de[[i]]` stays NULL and the
        // aggregation skips it.
        if in_group.is_empty() {
            return None;
        }
        let in_rest: Vec<usize> = (0..n_cells).filter(|&c| labels[c] != level).collect();
        if in_rest.is_empty() {
            return None;
        }
        Some((in_group, in_rest))
    };
    identify_over_expressed_selected(
        data,
        features,
        n_groups,
        thresh_pc,
        thresh_fc,
        thresh_p,
        only_pos,
        n_adjust,
        &mut select,
    )
}

/// [`identify_over_expressed`], with the per-group cell selection supplied by the caller.
///
/// Upstream picks `cell.use1` / `cell.use2` four different ways, and the choice is the *only*
/// thing `group.dataset` changes about the Wilcoxon branch -- the percentage filter, the
/// `mean.fxn` fold change, the `wilcox.test` and the Bonferroni multiplier are identical
/// whichever way the two cell sets are chosen:
///
/// ```r
/// if (is.null(group.dataset)) {
///   cell.use1 <- which(labels == level.use[i])
///   cell.use2 <- base::setdiff(1:length(labels), cell.use1)
/// } else if (!group.DE.combined) {
///   cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
///   cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
/// } else {
///   cell.use1 <- which(labels.dataset == pos.dataset)
///   cell.use2 <- which(labels.dataset != pos.dataset)
/// }
/// ```
///
/// Taking the selection as a parameter keeps the numerics in one place instead of growing a
/// `group.dataset` flag through the whole Wilcoxon path, and it is the same shape as
/// [`identify_over_expressed_one_group_n`], which already takes the two index vectors.
///
/// Returning `None` skips the level, which is how both "no cells in the group" and "no cells
/// outside it" are handled. Upstream's `pct.1` divides by `length(cell.use1)` and `pct.2` by
/// `length(cell.use2)`, so an empty set would be a division by zero followed by a `wilcox.test`
/// with no data; the `group.dataset` branches can produce one, and skipping matches what the
/// `is.null(gde)` check downstream then does.
pub fn identify_over_expressed_selected(
    data: &[f64],
    features: &[String],
    n_groups: usize,
    thresh_pc: f64,
    thresh_fc: f64,
    thresh_p: f64,
    only_pos: bool,
    n_adjust: usize,
    select: &mut dyn FnMut(usize) -> Option<(Vec<usize>, Vec<usize>)>,
) -> IdentifyResult {
    let mut out = Vec::new();
    let mut n_before_only_pos = 0usize;
    for level in 0..n_groups {
        let Some((in_group, in_rest)) = select(level) else {
            continue;
        };
        let mut g = identify_over_expressed_one_group_n(
            data, features, &in_group, &in_rest, thresh_pc, thresh_fc, thresh_p, only_pos, n_adjust,
        );
        for m in &mut g.rows {
            m.cluster = String::new();
            m.group_index = level;
        }
        n_before_only_pos += g.n_before_only_pos;
        out.extend(g.rows);
    }
    IdentifyResult {
        rows: out,
        n_before_only_pos,
    }
}

/// [`identify_over_expressed`]'s result, including the schema counter.
pub struct IdentifyResult {
    pub rows: Vec<Marker>,
    pub n_before_only_pos: usize,
}

/// `identifyOverExpressedGenes(..., do.DE = FALSE)`: the non-Wilcoxon branch.
///
/// ```r
/// markers.all <- data.frame(features = rownames(data.use),
///                           nCells = rowSums(data.use > 0))
/// markers.all <- dplyr::filter(markers.all, nCells >= min.cells)
/// ```
///
/// Note `min.cells` appears **only** here (and in the presto branch). The
/// `do.DE = TRUE` path never applies it -- a fact that is easy to miss because the
/// argument is in the same signature.
///
/// The third element is the feature's 1-based row index in `data.use`, because
/// `dplyr::filter` rewrites the surviving row names to those indices *as character*. The
/// marker table therefore has row names like `"17"`, not `"G17"` and not `"1"`, and a
/// filter that returns a subset of 29 of 30 features does not renumber from 1.
pub fn expressed_in_at_least_n_cells(
    data: &[f64],
    features: &[String],
    min_cells: usize,
) -> Vec<(String, f64, usize)> {
    let n_genes = features.len();
    let n_cells = data.len() / n_genes;
    assert_eq!(data.len(), n_genes * n_cells, "bad data buffer");
    (0..n_genes)
        .map(|g| {
            let n = (0..n_cells)
                .filter(|&c| data[c * n_genes + g] > 0.0)
                .count();
            (features[g].clone(), n as f64, g)
        })
        .filter(|&(_, n, _)| n >= min_cells as f64)
        .collect()
}

/// R's `round(x, 3)`. Exposed for the corpus test, which recomputes `pct.1`/`pct.2`
/// directly rather than trusting the marker table.
pub fn round3_for_tests(x: f64) -> f64 {
    round3(x)
}

/// R's `round(x, 3)`.
fn round3(x: f64) -> f64 {
    // R uses IEC 60559 `rint` for halfway cases, i.e. round-half-to-even, which Rust's
    // `f64::round_tie_even` matches (stable since 1.77).
    (x * 1000.0).round_ties_even() / 1000.0
}

/// R 4.3.3 `dwilcox` reference values. The counts are *not* normalised to 1 over the
/// full support: `cwilcox` applies `k > m*n/2 -> k = m*n - k`, so the upper half
/// mirrors the lower half. `pwilcox` never sums there, which is why its own
/// probabilities are right, but a test that sums the density over the full support
/// and expects 1 fails -- R gives 0.9 for `m=2, n=3`, and so does this port.
#[cfg(test)]
const DWILCOX_REF_: [(i64, i64, i64, f64); 231] = [
    (2, 3, 0, 0.10000000000000001),
    (2, 3, 1, 0.10000000000000001),
    (2, 3, 2, 0.20000000000000001),
    (2, 3, 3, 0.20000000000000001),
    (2, 3, 4, 0.20000000000000001),
    (2, 3, 5, 0.10000000000000001),
    (2, 3, 6, 0.10000000000000001),
    (3, 4, 0, 0.028571428571428571),
    (3, 4, 1, 0.028571428571428571),
    (3, 4, 2, 0.057142857142857141),
    (3, 4, 3, 0.085714285714285715),
    (3, 4, 4, 0.11428571428571428),
    (3, 4, 5, 0.11428571428571428),
    (3, 4, 6, 0.14285714285714285),
    (3, 4, 7, 0.11428571428571428),
    (3, 4, 8, 0.11428571428571428),
    (3, 4, 9, 0.085714285714285715),
    (3, 4, 10, 0.057142857142857141),
    (3, 4, 11, 0.028571428571428571),
    (3, 4, 12, 0.028571428571428571),
    (5, 5, 0, 0.003968253968253968),
    (5, 5, 1, 0.003968253968253968),
    (5, 5, 2, 0.0079365079365079361),
    (5, 5, 3, 0.011904761904761904),
    (5, 5, 4, 0.01984126984126984),
    (5, 5, 5, 0.027777777777777776),
    (5, 5, 6, 0.035714285714285712),
    (5, 5, 7, 0.043650793650793648),
    (5, 5, 8, 0.055555555555555552),
    (5, 5, 9, 0.063492063492063489),
    (5, 5, 10, 0.071428571428571425),
    (5, 5, 11, 0.075396825396825393),
    (5, 5, 12, 0.079365079365079361),
    (5, 5, 13, 0.079365079365079361),
    (5, 5, 14, 0.075396825396825393),
    (5, 5, 15, 0.071428571428571425),
    (5, 5, 16, 0.063492063492063489),
    (5, 5, 17, 0.055555555555555552),
    (5, 5, 18, 0.043650793650793648),
    (5, 5, 19, 0.035714285714285712),
    (5, 5, 20, 0.027777777777777776),
    (5, 5, 21, 0.01984126984126984),
    (5, 5, 22, 0.011904761904761904),
    (5, 5, 23, 0.0079365079365079361),
    (5, 5, 24, 0.003968253968253968),
    (5, 5, 25, 0.003968253968253968),
    (1, 3, 0, 0.25),
    (1, 3, 1, 0.25),
    (1, 3, 2, 0.25),
    (1, 3, 3, 0.25),
    (7, 2, 0, 0.027777777777777776),
    (7, 2, 1, 0.027777777777777776),
    (7, 2, 2, 0.055555555555555552),
    (7, 2, 3, 0.055555555555555552),
    (7, 2, 4, 0.083333333333333329),
    (7, 2, 5, 0.083333333333333329),
    (7, 2, 6, 0.1111111111111111),
    (7, 2, 7, 0.1111111111111111),
    (7, 2, 8, 0.1111111111111111),
    (7, 2, 9, 0.083333333333333329),
    (7, 2, 10, 0.083333333333333329),
    (7, 2, 11, 0.055555555555555552),
    (7, 2, 12, 0.055555555555555552),
    (7, 2, 13, 0.027777777777777776),
    (7, 2, 14, 0.027777777777777776),
    (4, 4, 0, 0.014285714285714285),
    (4, 4, 1, 0.014285714285714285),
    (4, 4, 2, 0.028571428571428571),
    (4, 4, 3, 0.042857142857142858),
    (4, 4, 4, 0.071428571428571425),
    (4, 4, 5, 0.071428571428571425),
    (4, 4, 6, 0.10000000000000001),
    (4, 4, 7, 0.10000000000000001),
    (4, 4, 8, 0.11428571428571428),
    (4, 4, 9, 0.10000000000000001),
    (4, 4, 10, 0.10000000000000001),
    (4, 4, 11, 0.071428571428571425),
    (4, 4, 12, 0.071428571428571425),
    (4, 4, 13, 0.042857142857142858),
    (4, 4, 14, 0.028571428571428571),
    (4, 4, 15, 0.014285714285714285),
    (4, 4, 16, 0.014285714285714285),
    (9, 9, 0, 2.0567667626491154e-05),
    (9, 9, 1, 2.0567667626491154e-05),
    (9, 9, 2, 4.1135335252982309e-05),
    (9, 9, 3, 6.1703002879473463e-05),
    (9, 9, 4, 0.00010283833813245579),
    (9, 9, 5, 0.00014397367338543808),
    (9, 9, 6, 0.00022624434389140272),
    (9, 9, 7, 0.00030851501439736734),
    (9, 9, 8, 0.00045248868778280545),
    (9, 9, 9, 0.00061703002879473468),
    (9, 9, 10, 0.00082270670505964628),
    (9, 9, 11, 0.0010695187165775401),
    (9, 9, 12, 0.0014191690662278898),
    (9, 9, 13, 0.0017893870835047306),
    (9, 9, 14, 0.0022830111065405183),
    (9, 9, 15, 0.0028383381324557796),
    (9, 9, 16, 0.0035170711641299875),
    (9, 9, 17, 0.0042575071986836691),
    (9, 9, 18, 0.0051624845742492802),
    (9, 9, 19, 0.0061085972850678733),
    (9, 9, 20, 0.0072398190045248872),
    (9, 9, 21, 0.0084533113944878658),
    (9, 9, 22, 0.0097902097902097911),
    (9, 9, 23, 0.01120937885643768),
    (9, 9, 24, 0.012793089263677499),
    (9, 9, 25, 0.014376799670917317),
    (9, 9, 26, 0.016083916083916083),
    (9, 9, 27, 0.017832167832167831),
    (9, 9, 28, 0.019621554915672561),
    (9, 9, 29, 0.021390374331550801),
    (9, 9, 30, 0.023200329082682023),
    (9, 9, 31, 0.024886877828054297),
    (9, 9, 32, 0.026573426573426574),
    (9, 9, 33, 0.0281365693130399),
    (9, 9, 34, 0.029555738379267792),
    (9, 9, 35, 0.030830933772110242),
    (9, 9, 36, 0.031982723159193746),
    (9, 9, 37, 0.032867132867132866),
    (9, 9, 38, 0.033566433566433566),
    (9, 9, 39, 0.034060057589469353),
    (9, 9, 40, 0.034286301933360755),
    (9, 9, 41, 0.034286301933360755),
    (9, 9, 42, 0.034060057589469353),
    (9, 9, 43, 0.033566433566433566),
    (9, 9, 44, 0.032867132867132866),
    (9, 9, 45, 0.031982723159193746),
    (9, 9, 46, 0.030830933772110242),
    (9, 9, 47, 0.029555738379267792),
    (9, 9, 48, 0.0281365693130399),
    (9, 9, 49, 0.026573426573426574),
    (9, 9, 50, 0.024886877828054297),
    (9, 9, 51, 0.023200329082682023),
    (9, 9, 52, 0.021390374331550801),
    (9, 9, 53, 0.019621554915672561),
    (9, 9, 54, 0.017832167832167831),
    (9, 9, 55, 0.016083916083916083),
    (9, 9, 56, 0.014376799670917317),
    (9, 9, 57, 0.012793089263677499),
    (9, 9, 58, 0.01120937885643768),
    (9, 9, 59, 0.0097902097902097911),
    (9, 9, 60, 0.0084533113944878658),
    (9, 9, 61, 0.0072398190045248872),
    (9, 9, 62, 0.0061085972850678733),
    (9, 9, 63, 0.0051624845742492802),
    (9, 9, 64, 0.0042575071986836691),
    (9, 9, 65, 0.0035170711641299875),
    (9, 9, 66, 0.0028383381324557796),
    (9, 9, 67, 0.0022830111065405183),
    (9, 9, 68, 0.0017893870835047306),
    (9, 9, 69, 0.0014191690662278898),
    (9, 9, 70, 0.0010695187165775401),
    (9, 9, 71, 0.00082270670505964628),
    (9, 9, 72, 0.00061703002879473468),
    (9, 9, 73, 0.00045248868778280545),
    (9, 9, 74, 0.00030851501439736734),
    (9, 9, 75, 0.00022624434389140272),
    (9, 9, 76, 0.00014397367338543808),
    (9, 9, 77, 0.00010283833813245579),
    (9, 9, 78, 6.1703002879473463e-05),
    (9, 9, 79, 4.1135335252982309e-05),
    (9, 9, 80, 2.0567667626491154e-05),
    (9, 9, 81, 2.0567667626491154e-05),
    (2, 2, 0, 0.16666666666666666),
    (2, 2, 1, 0.16666666666666666),
    (2, 2, 2, 0.33333333333333331),
    (2, 2, 3, 0.16666666666666666),
    (2, 2, 4, 0.16666666666666666),
    (6, 5, 0, 0.0021645021645021645),
    (6, 5, 1, 0.0021645021645021645),
    (6, 5, 2, 0.004329004329004329),
    (6, 5, 3, 0.0064935064935064939),
    (6, 5, 4, 0.010822510822510822),
    (6, 5, 5, 0.015151515151515152),
    (6, 5, 6, 0.021645021645021644),
    (6, 5, 7, 0.025974025974025976),
    (6, 5, 8, 0.034632034632034632),
    (6, 5, 9, 0.041125541125541128),
    (6, 5, 10, 0.049783549783549784),
    (6, 5, 11, 0.054112554112554112),
    (6, 5, 12, 0.062770562770562768),
    (6, 5, 13, 0.064935064935064929),
    (6, 5, 14, 0.069264069264069264),
    (6, 5, 15, 0.069264069264069264),
    (6, 5, 16, 0.069264069264069264),
    (6, 5, 17, 0.064935064935064929),
    (6, 5, 18, 0.062770562770562768),
    (6, 5, 19, 0.054112554112554112),
    (6, 5, 20, 0.049783549783549784),
    (6, 5, 21, 0.041125541125541128),
    (6, 5, 22, 0.034632034632034632),
    (6, 5, 23, 0.025974025974025976),
    (6, 5, 24, 0.021645021645021644),
    (6, 5, 25, 0.015151515151515152),
    (6, 5, 26, 0.010822510822510822),
    (6, 5, 27, 0.0064935064935064939),
    (6, 5, 28, 0.004329004329004329),
    (6, 5, 29, 0.0021645021645021645),
    (6, 5, 30, 0.0021645021645021645),
    (3, 10, 0, 0.0034965034965034965),
    (3, 10, 1, 0.0034965034965034965),
    (3, 10, 2, 0.006993006993006993),
    (3, 10, 3, 0.01048951048951049),
    (3, 10, 4, 0.013986013986013986),
    (3, 10, 5, 0.017482517482517484),
    (3, 10, 6, 0.024475524475524476),
    (3, 10, 7, 0.027972027972027972),
    (3, 10, 8, 0.034965034965034968),
    (3, 10, 9, 0.04195804195804196),
    (3, 10, 10, 0.048951048951048952),
    (3, 10, 11, 0.052447552447552448),
    (3, 10, 12, 0.05944055944055944),
    (3, 10, 13, 0.05944055944055944),
    (3, 10, 14, 0.062937062937062943),
    (3, 10, 15, 0.062937062937062943),
    (3, 10, 16, 0.062937062937062943),
    (3, 10, 17, 0.05944055944055944),
    (3, 10, 18, 0.05944055944055944),
    (3, 10, 19, 0.052447552447552448),
    (3, 10, 20, 0.048951048951048952),
    (3, 10, 21, 0.04195804195804196),
    (3, 10, 22, 0.034965034965034968),
    (3, 10, 23, 0.027972027972027972),
    (3, 10, 24, 0.024475524475524476),
    (3, 10, 25, 0.017482517482517484),
    (3, 10, 26, 0.013986013986013986),
    (3, 10, 27, 0.01048951048951049),
    (3, 10, 28, 0.006993006993006993),
    (3, 10, 29, 0.0034965034965034965),
    (3, 10, 30, 0.0034965034965034965),
];

/// R 4.3.3 `pwilcox` reference values over 10 `(m, n)` shapes and the boundary
/// quanta `-1`, `0`, `1`, `floor(m*n/2)`, `m*n-1`, `m*n`, `m*n+3`, both tails.
#[cfg(test)]
const PWILCOX_REF_: [(i64, i64, f64, bool, f64); 140] = [
    (2, 3, -1.0, true, 0.0),
    (2, 3, -1.0, false, 1.0),
    (2, 3, 0.0, true, 0.10000000000000001),
    (2, 3, 0.0, false, 0.90000000000000002),
    (2, 3, 1.0, true, 0.20000000000000001),
    (2, 3, 1.0, false, 0.80000000000000004),
    (2, 3, 3.0, true, 0.60000000000000009),
    (2, 3, 3.0, false, 0.39999999999999991),
    (2, 3, 5.0, true, 0.90000000000000002),
    (2, 3, 5.0, false, 0.10000000000000001),
    (2, 3, 6.0, true, 1.0),
    (2, 3, 6.0, false, 0.0),
    (2, 3, 9.0, true, 1.0),
    (2, 3, 9.0, false, 0.0),
    (3, 4, -1.0, true, 0.0),
    (3, 4, -1.0, false, 1.0),
    (3, 4, 0.0, true, 0.028571428571428571),
    (3, 4, 0.0, false, 0.97142857142857142),
    (3, 4, 1.0, true, 0.057142857142857141),
    (3, 4, 1.0, false, 0.94285714285714284),
    (3, 4, 6.0, true, 0.5714285714285714),
    (3, 4, 6.0, false, 0.4285714285714286),
    (3, 4, 11.0, true, 0.97142857142857142),
    (3, 4, 11.0, false, 0.028571428571428571),
    (3, 4, 12.0, true, 1.0),
    (3, 4, 12.0, false, 0.0),
    (3, 4, 15.0, true, 1.0),
    (3, 4, 15.0, false, 0.0),
    (5, 5, -1.0, true, 0.0),
    (5, 5, -1.0, false, 1.0),
    (5, 5, 0.0, true, 0.003968253968253968),
    (5, 5, 0.0, false, 0.99603174603174605),
    (5, 5, 1.0, true, 0.0079365079365079361),
    (5, 5, 1.0, false, 0.99206349206349209),
    (5, 5, 12.0, true, 0.5),
    (5, 5, 12.0, false, 0.5),
    (5, 5, 24.0, true, 0.99603174603174605),
    (5, 5, 24.0, false, 0.003968253968253968),
    (5, 5, 25.0, true, 1.0),
    (5, 5, 25.0, false, 0.0),
    (5, 5, 28.0, true, 1.0),
    (5, 5, 28.0, false, 0.0),
    (1, 3, -1.0, true, 0.0),
    (1, 3, -1.0, false, 1.0),
    (1, 3, 0.0, true, 0.25),
    (1, 3, 0.0, false, 0.75),
    (1, 3, 1.0, true, 0.5),
    (1, 3, 1.0, false, 0.5),
    (1, 3, 1.0, true, 0.5),
    (1, 3, 1.0, false, 0.5),
    (1, 3, 2.0, true, 0.75),
    (1, 3, 2.0, false, 0.25),
    (1, 3, 3.0, true, 1.0),
    (1, 3, 3.0, false, 0.0),
    (1, 3, 6.0, true, 1.0),
    (1, 3, 6.0, false, 0.0),
    (7, 2, -1.0, true, 0.0),
    (7, 2, -1.0, false, 1.0),
    (7, 2, 0.0, true, 0.027777777777777776),
    (7, 2, 0.0, false, 0.97222222222222221),
    (7, 2, 1.0, true, 0.055555555555555552),
    (7, 2, 1.0, false, 0.94444444444444442),
    (7, 2, 7.0, true, 0.55555555555555558),
    (7, 2, 7.0, false, 0.44444444444444442),
    (7, 2, 13.0, true, 0.97222222222222221),
    (7, 2, 13.0, false, 0.027777777777777776),
    (7, 2, 14.0, true, 1.0),
    (7, 2, 14.0, false, 0.0),
    (7, 2, 17.0, true, 1.0),
    (7, 2, 17.0, false, 0.0),
    (2, 2, -1.0, true, 0.0),
    (2, 2, -1.0, false, 1.0),
    (2, 2, 0.0, true, 0.16666666666666666),
    (2, 2, 0.0, false, 0.83333333333333337),
    (2, 2, 1.0, true, 0.33333333333333331),
    (2, 2, 1.0, false, 0.66666666666666674),
    (2, 2, 2.0, true, 0.66666666666666663),
    (2, 2, 2.0, false, 0.33333333333333337),
    (2, 2, 3.0, true, 0.83333333333333337),
    (2, 2, 3.0, false, 0.16666666666666666),
    (2, 2, 4.0, true, 1.0),
    (2, 2, 4.0, false, 0.0),
    (2, 2, 7.0, true, 1.0),
    (2, 2, 7.0, false, 0.0),
    (9, 9, -1.0, true, 0.0),
    (9, 9, -1.0, false, 1.0),
    (9, 9, 0.0, true, 2.0567667626491154e-05),
    (9, 9, 0.0, false, 0.99997943233237352),
    (9, 9, 1.0, true, 4.1135335252982309e-05),
    (9, 9, 1.0, false, 0.99995886466474704),
    (9, 9, 40.0, true, 0.49999999999999994),
    (9, 9, 40.0, false, 0.5),
    (9, 9, 80.0, true, 0.99997943233237352),
    (9, 9, 80.0, false, 2.0567667626491154e-05),
    (9, 9, 81.0, true, 1.0),
    (9, 9, 81.0, false, 0.0),
    (9, 9, 84.0, true, 1.0),
    (9, 9, 84.0, false, 0.0),
    (12, 20, -1.0, true, 0.0),
    (12, 20, -1.0, false, 1.0),
    (12, 20, 0.0, true, 4.4288383989501176e-09),
    (12, 20, 0.0, false, 0.99999999557116159),
    (12, 20, 1.0, true, 8.8576767979002353e-09),
    (12, 20, 1.0, false, 0.99999999114232319),
    (12, 20, 120.0, true, 0.50764372776390965),
    (12, 20, 120.0, false, 0.49235627223609035),
    (12, 20, 239.0, true, 0.99999999557116159),
    (12, 20, 239.0, false, 4.4288383989501176e-09),
    (12, 20, 240.0, true, 1.0),
    (12, 20, 240.0, false, 0.0),
    (12, 20, 243.0, true, 1.0),
    (12, 20, 243.0, false, 0.0),
    (4, 4, -1.0, true, 0.0),
    (4, 4, -1.0, false, 1.0),
    (4, 4, 0.0, true, 0.014285714285714285),
    (4, 4, 0.0, false, 0.98571428571428577),
    (4, 4, 1.0, true, 0.028571428571428571),
    (4, 4, 1.0, false, 0.97142857142857142),
    (4, 4, 8.0, true, 0.55714285714285716),
    (4, 4, 8.0, false, 0.44285714285714284),
    (4, 4, 15.0, true, 0.98571428571428577),
    (4, 4, 15.0, false, 0.014285714285714285),
    (4, 4, 16.0, true, 1.0),
    (4, 4, 16.0, false, 0.0),
    (4, 4, 19.0, true, 1.0),
    (4, 4, 19.0, false, 0.0),
    (6, 3, -1.0, true, 0.0),
    (6, 3, -1.0, false, 1.0),
    (6, 3, 0.0, true, 0.011904761904761904),
    (6, 3, 0.0, false, 0.98809523809523814),
    (6, 3, 1.0, true, 0.023809523809523808),
    (6, 3, 1.0, false, 0.97619047619047616),
    (6, 3, 9.0, true, 0.54761904761904756),
    (6, 3, 9.0, false, 0.45238095238095244),
    (6, 3, 17.0, true, 0.98809523809523814),
    (6, 3, 17.0, false, 0.011904761904761904),
    (6, 3, 18.0, true, 1.0),
    (6, 3, 18.0, false, 0.0),
    (6, 3, 21.0, true, 1.0),
    (6, 3, 21.0, false, 0.0),
];

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn rank_averages_ties() {
        assert_eq!(rank_average(&[1.0, 2.0, 3.0]), vec![1.0, 2.0, 3.0]);
        // R: rank(c(1,1,2)) = c(1.5, 1.5, 3)
        assert_eq!(rank_average(&[1.0, 1.0, 2.0]), vec![1.5, 1.5, 3.0]);
        // R: rank(c(1,1,1,1)) = c(2.5, 2.5, 2.5, 2.5)
        assert_eq!(rank_average(&[5.0; 4]), vec![2.5; 4]);
        // R: rank(c(3,1,2)) = c(3, 1, 2)
        assert_eq!(rank_average(&[3.0, 1.0, 2.0]), vec![3.0, 1.0, 2.0]);
    }

    #[test]
    fn rank_keeps_na_by_ranking_it_last() {
        // R 4.3.3: rank(c(1, NA, 2)) = c(1, 3, 2); rank(c(NA, NA, 1)) = c(2, 3, 1);
        // rank(c(1, NA, NA, 2)) = c(1, 3, 4, 2).
        assert_eq!(rank_average(&[1.0, f64::NAN, 2.0]), vec![1.0, 3.0, 2.0]);
        assert_eq!(
            rank_average(&[f64::NAN, f64::NAN, 1.0]),
            vec![2.0, 3.0, 1.0]
        );
        assert_eq!(
            rank_average(&[1.0, f64::NAN, f64::NAN, 2.0]),
            vec![1.0, 3.0, 4.0, 2.0]
        );
        assert_eq!(
            rank_average(&[f64::NAN, 1.0, f64::NAN, 2.0]),
            vec![3.0, 1.0, 4.0, 2.0]
        );
    }

    #[test]
    fn p_adjust_bonferroni_is_pmin_one_n_times_p() {
        assert_eq!(
            p_adjust_bonferroni(&[0.01, 0.5, 1.0], 100.0),
            vec![1.0, 1.0, 1.0]
        );
        assert_eq!(p_adjust_bonferroni(&[0.001, 0.002], 10.0), vec![0.01, 0.02]);
    }

    #[test]
    fn pwilcox_matches_r_bit_for_bit() {
        for &(m, n, q, lower_tail, want) in PWILCOX_REF_.iter() {
            let got = pwilcox(q, m as f64, n as f64, lower_tail);
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "pwilcox({q}, {m}, {n}, lower.tail={lower_tail}) = {got:.17e}, R gives {want:.17e}"
            );
        }
    }

    #[test]
    fn dwilcox_matches_r_bit_for_bit() {
        for &(m, n, k, want) in DWILCOX_REF_.iter() {
            let c = choose(m as f64 + n as f64, n as f64);
            let mut memo = HashMap::new();
            let got = cwilcox(k, m, n, &mut memo) / c;
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "dwilcox({k}, {m}, {n}) = {got:.17e}, R gives {want:.17e}"
            );
        }
    }

    #[test]
    fn r_sign_returns_zero_for_zero_which_rust_signum_does_not() {
        assert_eq!(r_sign(0.0), 0.0);
        assert_eq!(r_sign(-0.0), -0.0);
        assert_eq!(r_sign(3.0), 1.0);
        assert_eq!(r_sign(-3.0), -1.0);
        // The trap: Rust's signum returns +1.0 for +0.0.
        assert_eq!(0.0f64.signum(), 1.0);
        assert_ne!(r_sign(0.0), 0.0f64.signum());
    }

    #[test]
    fn an_all_ties_gene_gives_a_nan_p_value_not_zero() {
        // Every value equal => one tie group of size 80 => SIGMA == 0 => z = 0/0 = NaN.
        // R returns NaN; an implementation that gets the continuity correction wrong here
        // returns -inf and a p-value of 0, which would call every such gene significant.
        let (w, p, exact) = wilcox_test_two_sided(&vec![0.0; 30], &vec![0.0; 50]);
        assert_eq!(w, 750.0);
        assert!(!exact);
        assert!(p.is_nan(), "p = {p}, R gives NaN");
    }

    #[test]
    fn round3_is_rounding_to_three_decimal_places_half_to_even() {
        // Not significant digits: `round(0.5, 3)` is 0.5. And the halfway cases go to even,
        // because R uses `rint` under IEC 60559.
        assert_eq!(round3(0.5), 0.5);
        assert_eq!(round3(0.1235), 0.124);
        assert_eq!(round3(0.1234), 0.123);
        assert_eq!(round3(0.0005), 0.0);
        assert_eq!(round3(0.0015), 0.002);
        // Reference from R: round(1/3, 3) = 0.33300000000000002, which is the *double*
        // nearest 0.333, not a re-rounded string.
        assert_eq!(
            round3(1.0 / 3.0).to_bits(),
            0.33300000000000002f64.to_bits()
        );
        assert_eq!(
            round3(2.0 / 3.0).to_bits(),
            0.66700000000000004f64.to_bits()
        );
        assert_eq!(
            round3(1.0 / 6.0).to_bits(),
            0.16700000000000001f64.to_bits()
        );
    }

    #[test]
    fn mean_fxn_is_log_one_plus_mean_expm1() {
        assert_eq!(mean_fxn(&[0.0, 0.0]), 0.0);
        let x = [0.1f64, 0.2, 0.3];
        let want = (((x[0].exp_m1() + x[1].exp_m1()) + x[2].exp_m1()) / 3.0 + 1.0).ln();
        assert!((mean_fxn(&x) - want).abs() < 1e-15);
    }
}

#[cfg(test)]
mod single_candidate_tests {
    use super::{identify_over_expressed_one_group_n, GroupMarkers};

    fn data(n_genes: usize, n_cells: usize, f: impl Fn(usize, usize) -> f64) -> Vec<f64> {
        (0..n_cells)
            .flat_map(|c| (0..n_genes).map(move |g| (c, g)))
            .map(|(c, g)| f(g, c))
            .collect()
    }

    /// Upstream's `apply(X, 1, FUN)` returns an unnamed scalar for a one-row matrix, so
    /// `FC[features]` is `NA` and the group yields no markers at all. The port has to agree, or
    /// it invents a marker for every group with exactly one percentage-passing feature.
    #[test]
    fn a_single_pct_passing_feature_yields_no_markers() {
        // Gene 0 is expressed in every cell of `in_group` and no cell of `in_rest`, so it is the
        // only feature clearing `thresh.pc = 0`; gene 1 is all zero in both, so `pct.2` and
        // `pct.1` are 0 for it and it is excluded.
        let d = data(2, 4, |g, c| if g == 0 && c < 2 { 1.0 } else { 0.0 });
        let features = vec!["gene0".to_string(), "gene1".to_string()];
        let out: GroupMarkers = identify_over_expressed_one_group_n(
            &d,
            &features,
            &[0, 1],
            &[2, 3],
            0.0,
            0.0,
            0.05,
            true,
            2,
        );
        assert!(
            out.rows.is_empty(),
            "one candidate must not become a marker"
        );
        assert_eq!(out.n_before_only_pos, 0);
    }

    /// Two candidates is the smallest input that produces a marker, and it must still work --
    /// the previous test's rule is about `length(features) == 1`, not about being small.
    ///
    /// The split is 5 against 5 rather than 2 against 2 on purpose: with 2 and 2 the smallest
    /// attainable two-sided Wilcoxon p-value is 1/6, and `p.adjust(bonferroni, n = 3)` turns
    /// that into 0.5, so *no* feature can clear `thresh.p` and the test would pass for the
    /// wrong reason.
    #[test]
    fn two_pct_passing_features_still_produce_markers() {
        let d = data(3, 10, |g, c| if g < 2 && c < 5 { 1.0 } else { 0.0 });
        let features = ["gene0", "gene1", "gene2"]
            .iter()
            .map(|s| s.to_string())
            .collect::<Vec<_>>();
        let out: GroupMarkers = identify_over_expressed_one_group_n(
            &d,
            &features,
            &[0, 1, 2, 3, 4],
            &[5, 6, 7, 8, 9],
            0.0,
            0.0,
            0.05,
            true,
            3,
        );
        assert_eq!(
            out.rows.len(),
            2,
            "both candidates clear pct, fc and thresh.p"
        );
        assert_eq!(out.n_before_only_pos, 2);
    }
}
