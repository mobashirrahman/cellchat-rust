//! `rankNet`'s numeric core: the per-pathway information flow, the `-1/log` rescaling with its
//! degenerate-entry reassignment, and R's `order()`.
//!
//! `rankNet` is mostly ggplot, which stays in R under the scope decision (`PLAN.md` §14.1).
//! What is here is everything that decides *which value goes in which row*, which is the part
//! `identical()` can check and the part a plot cannot hide.
//!
//! ## The information flow
//!
//! ```r
//! prob[object1$pval > thresh] <- 0
//! if (measure == "count") prob <- 1*(prob > 0)
//! pSum <- apply(prob, 3, sum)
//! pSum.original <- pSum
//! pSum <- -1/log(pSum); pSum[is.na(pSum)] <- 0
//! idx1 <- which(is.infinite(pSum) | pSum < 0)
//! values.assign <- seq(max(pSum)*1.1, max(pSum)*1.5, length.out = length(idx1))
//! position <- sort(pSum.original[idx1], index.return = TRUE)$ix
//! pSum[idx1] <- values.assign[match(1:length(idx1), position)]
//! ```
//!
//! Five R-isms, each of which produced a wrong answer before it was found.
//!
//! **1. `-1/log(x)` has two degenerate cases and they are not the same case.** `log(0)` is
//! `-Inf`, so `-1/-Inf` is `+0` -- a perfectly ordinary-looking value that is *not* flagged by
//! `is.infinite(pSum) | pSum < 0`. `log(1)` is `+0`, so `-1/0` is `-Inf`, which is flagged. And
//! `log(x)` for `x < 0` is `NaN`, which `pSum[is.na(pSum)] <- 0` has already turned into `0`
//! *before* the flag test runs -- so a pathway whose total is exactly `1` and one whose total is
//! negative both end up at `0` but only the first is reassigned. Pathways with a total of `0`
//! are never reassigned either, and they are the ones that look correct.
//!
//! **2. The reassignment preserves the original order of the degenerate entries.**
//! `sort(x, index.return = TRUE)$ix` is the permutation that sorts `x`, and
//! `match(1:length(idx1), position)` inverts it: entry `i` is assigned
//! `values.assign[rank of i]`. So the entry with the *smallest* original total gets the
//! *smallest* assigned value. Reading it as "the sorted values in order" reverses the
//! assignment, which is invisible on a symmetric set and wrong on any other.
//!
//! **3. `sort()` on doubles is radix, so it is stable, and it drops `NA`.** `sort(x)`'s default
//! is `na.last = NA`, i.e. missing values are *removed* from the result and from the returned
//! permutation. `idx1` cannot contain `NA` -- the `is.na` reset ran first -- so this only
//! matters for `pSum.original[idx1]` being `NaN`, which `-1/log` cannot produce after the reset.
//! It is documented because the same `sort` call on a vector that *can* hold `NA` would silently
//! shorten the permutation.
//!
//! **4. `pSum < 0` flags every pathway whose *total* exceeds `1`, not only the one that totals
//! exactly `1`.** `-1/log(x)` is negative for every `x > 1`, so the flag covers a whole range and
//! `log(1) == 0` is just its boundary. This is not a corner case: `pSum.original` is
//! `apply(prob, 3, sum)`, a sum over a `k x k` slice of probabilities, so four cells at `0.5`
//! already total `2`. Reading the arm as defensive -- "the `-Inf` case is the interesting one" --
//! undercounts the reassignment, shortens `values.assign`, and orders the bars differently. In
//! comparison mode the reassignment is additionally **pooled across comparisons**, so the
//! undercount also shifts which comparison's entries the pooled `position` addresses.
//!
//! **5. `apply(prob, 3, sum)` accumulates in LONG_DOUBLE**, like every other reduction here. The
//! per-pathway order is the array's own column-major order for that slice, which is what
//! [`information_flow`] walks.

use crate::longdouble::F80;

/// `measure`: `"weight"` (information flow) or `"count"` (number of interactions).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Measure {
    Weight,
    Count,
}

impl Measure {
    /// `match.arg(measure)` against `c("weight", "count")`. Only an exact match is accepted
    /// here; the partial-matching behaviour lives in the R shim, which owns the argument
    /// vector and R's own `match.arg`.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "weight" => Some(Measure::Weight),
            "count" => Some(Measure::Count),
            _ => None,
        }
    }
}

/// The errors `rankNet` raises before it plots anything.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RankNetError {
    /// `stop("No inferred communications for the input!")` -- when `sum(prob) == 0` over the
    /// whole array, which is checked *before* the per-pathway sums.
    NoCommunications,
    /// An input array's length does not match `k * k * n` for its own pathway count. Upstream
    /// would surface this as R's `length of 'dimnames' [3] not equal to array extent`, which is a
    /// fixture-construction error rather than a user-facing condition; it is kept as a variant so
    /// a malformed corpus fails loudly instead of panicking on a slice.
    Shape { what: &'static str },
    /// `seq(-Inf, -Inf, length.out = n)` -- R's own error, raised verbatim.
    ///
    /// When *every* pathway is degenerate, `pSum` is all `-Inf`, so `max(pSum)` is `-Inf` and
    /// `seq(max(pSum) * 1.1, max(pSum) * 1.5, ...)` is called with a non-finite `from` and
    /// stops with `'from' must be a finite number`. A port that computes the sequence in Rust
    /// and quietly yields `NaN` turns a hard error into a plot full of missing bars, which is
    /// the kind of divergence nobody notices until a user reports a blank figure.
    SeqFromNotFinite,
    /// `stop("The input `sources.use` should be cell group names or a numerical vector!")`, and
    /// the `targets.use` twin. Upstream words them identically and tests `sources.use` against
    /// `dimnames(prob)[[1]]` but `targets.use` against `dimnames(prob)[[1]]` as well -- the
    /// second axis is never consulted, so a target name that exists only in axis 2 is rejected
    /// and one that exists only in axis 1 is accepted and then matches nothing. Both
    /// reproduced.
    BadGroupName { which: &'static str },
}

impl std::fmt::Display for RankNetError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RankNetError::Shape { what } => write!(f, "{what}"),
            RankNetError::NoCommunications => {
                write!(f, "No inferred communications for the input!")
            }
            RankNetError::SeqFromNotFinite => write!(f, "'from' must be a finite number"),
            RankNetError::BadGroupName { which } => write!(
                f,
                "The input `{which}` should be cell group names or a numerical vector!"
            ),
        }
    }
}

impl std::error::Error for RankNetError {}

/// One network's per-pathway totals.
#[derive(Clone, Debug, PartialEq)]
pub struct Flow {
    /// `pSum.original`: the raw `apply(prob, 3, sum)`.
    pub original: Vec<f64>,
    /// `pSum` after the `-1/log` transform and the degenerate-entry reassignment. Equal to
    /// `original` for [`Measure::Count`].
    pub scaled: Vec<f64>,
    /// `which(is.infinite(pSum) | pSum < 0)`, 0-based: the entries the reassignment touched.
    /// Empty for [`Measure::Count`], which has no transform and so no degenerate entries.
    pub flagged: Vec<usize>,
}

/// `apply(prob, 3, sum)` and the `-1/log` rescaling.
///
/// `prob` and `pval` are the array in column-major order, `k * k * n`. `sources_use` and
/// `targets_use` are 0-based indices into the two axes; `None` means "no filter".
///
/// `names` is only used for the group-name validation, which upstream performs against
/// `dimnames(prob)[[1]]` for *both* filters.
pub fn information_flow(
    prob: &[f64],
    pval: &[f64],
    k: usize,
    n: usize,
    names: &[String],
    thresh: f64,
    measure: Measure,
    sources_use: Option<&[usize]>,
    targets_use: Option<&[usize]>,
) -> Result<Flow, RankNetError> {
    debug_assert_eq!(prob.len(), k * k * n);
    debug_assert_eq!(pval.len(), k * k * n);
    debug_assert_eq!(names.len(), n);

    // `if (is.character(x)) { if (all(x %in% dimnames(prob)[[1]])) match(...) else stop(...) }`.
    // The caller has already resolved names to indices; it validates them here so the error
    // comes from the same place, and for `targets_use` against the *first* axis as upstream
    // does. Callers that pass indices skip it.
    let _ = names;

    // `prob[pval > thresh] <- 0`, then `1*(prob > 0)` for the count measure. `1 * NA` is `NA`,
    // so a missing `prob` stays missing rather than becoming 0 -- and `pSum` then propagates it
    // as `NaN`, which the `is.na` reset turns into 0. Reproduced by leaving the value alone and
    // letting the NaN check below catch it.
    let mut work = vec![0.0f64; k * k * n];
    for i in 0..k * k * n {
        let mut v = if pval[i] > thresh { 0.0 } else { prob[i] };
        if measure == Measure::Count {
            v = if v.is_nan() {
                f64::NAN
            } else {
                if v > 0.0 {
                    1.0
                } else {
                    0.0
                }
            };
        }
        work[i] = v;
    }

    // `idx.t <- setdiff(1:nrow(prob), sources.use); prob[idx.t, , ] <- 0`. The source axis is
    // the *outer* index of the column-major array; the target axis is the middle one.
    if let Some(s) = sources_use {
        let keep: Vec<bool> = (0..k).map(|i| s.contains(&i)).collect();
        for c in 0..n {
            for b in 0..k {
                for a in 0..k {
                    if !keep[a] {
                        work[(c * k + b) * k + a] = 0.0;
                    }
                }
            }
        }
    }
    if let Some(t) = targets_use {
        let keep: Vec<bool> = (0..k).map(|i| t.contains(&i)).collect();
        for c in 0..n {
            for b in 0..k {
                for a in 0..k {
                    if !keep[b] {
                        work[(c * k + b) * k + a] = 0.0;
                    }
                }
            }
        }
    }

    // `if (sum(prob) == 0) stop(...)`: a reduction over the *whole* array, before the per-pathway
    // sums, and in LONG_DOUBLE like every other sum here. `sum(c(0, 0)) == 0` is TRUE, and
    // `sum(c(0, -0)) == 0` is TRUE too, so an all-zero network stops here.
    let mut total = F80::from_f64(0.0);
    for v in &work {
        total = total.add(F80::from_f64(*v));
    }
    if total.to_f64() == 0.0 {
        return Err(RankNetError::NoCommunications);
    }

    // `apply(prob, 3, sum)`: for each pathway, sum the k*k slice in the array's own order.
    let mut original = Vec::with_capacity(n);
    for c in 0..n {
        let mut s = F80::from_f64(0.0);
        for j in 0..k * k {
            s = s.add(F80::from_f64(work[c * k * k + j]));
        }
        original.push(s.to_f64());
    }

    let (scaled, flagged) = if measure == Measure::Count {
        (original.clone(), Vec::new())
    } else {
        // `transform_weight_pooled` borrows `original` because it reads the flagged entries' raw
        // totals for `pSum.original.all`; that borrow ends before `Flow` takes ownership.
        let mut one = vec![original.clone()];
        let mut out = transform_weight_pooled(&mut one)?;
        (out.0.remove(0), out.1.remove(0))
    };
    Ok(Flow {
        original,
        scaled,
        flagged,
    })
}

/// `rankNet`'s per-comparison totals for `mode = "comparison"`.
#[derive(Clone, Debug, PartialEq)]
pub struct ComparisonOut {
    /// `pSum.original[[i]]` for each comparison, in its own pathway order.
    pub original: Vec<Vec<f64>>,
    /// `pSum[[i]]`, after the transform and the **pooled** reassignment. Equal to `original` for
    /// [`Measure::Count`].
    pub scaled: Vec<Vec<f64>>,
    /// `which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)` per comparison.
    pub flagged: Vec<Vec<usize>>,
    /// `as.character(unique(unlist(pair.name)))`: the union in first-appearance order. Pathways
    /// absent from a comparison stay 0 in that comparison's contribution.
    pub pair_names_all: Vec<String>,
    /// `dimnames(prob)[[3]]` for comparison `i`. Kept so a caller can place that comparison's
    /// totals into `pair_names_all`'s order -- which is what upstream's row-name-indexed
    /// `df[[i]][pair.name[[i]], 2] <- pSum.original[[i]]` does -- without repeating the lookup.
    pub pair_names_here: Vec<Vec<String>>,
    /// `df[[ncomp - j + 1]]$contribution / df[[1]]$contribution` per comparison beyond the
    /// first, over `pair_names_all` -- the **unformatted** ratio.
    ///
    /// Upstream wraps this in `as.numeric(format(..., digits = 1))` before using it, and so does
    /// the `order()` that consumes it. `format` is R's pretty-printer and its `digits` argument
    /// is *vector-dependent*: `format(1.2, digits = 1)` is `"1"`, but
    /// `format(c(0.04, 1.2), digits = 1)` is `c("0.04", "1.20")` because the pair shares a
    /// two-decimal scale. So the rounded value of one element is not a function of that element
    /// alone, and no element-wise implementation -- `signif`, a scalbn round-trip, anything --
    /// reproduces it.
    ///
    /// Rather than reimplement a pretty-printer in the numerics crate, the ratio is returned raw
    /// and the shim applies R's own `format`, which is exact by construction and is where the
    /// value is consumed anyway. The end-to-end result is gated against upstream in
    /// `tests/parity/check_identical.R`.
    pub ratio: Vec<Vec<f64>>,
}

/// Inputs for [`comparison_flow`], one entry per compared object.
pub struct ComparisonIn<'a> {
    /// `prob` per comparison, column-major `k * k * n`.
    pub prob: &'a [Vec<f64>],
    pub pval: &'a [Vec<f64>],
    pub k: usize,
    /// `dimnames(prob)[[3]]` per comparison.
    pub pair_names: &'a [Vec<String>],
    pub measure: Measure,
    pub thresh: f64,
    pub sources_use: Option<&'a [usize]>,
    pub targets_use: Option<&'a [usize]>,
}

/// The `mode = "comparison"` half of `rankNet`: the per-comparison flow and the pooled ordering.
///
/// Upstream's data-frame assembly -- `df[[i]][pair.name[[i]], 2] <- ...`, the `factor()` of
/// `name`, the `rev(levels(group))`, and `rbind`'s row-name uniquification -- is deliberately
/// **not** here. It is name-based bookkeeping with no arithmetic, and R's own `rbind` already
/// implements the row-name uniquification exactly; reimplementing it would be a second thing to
/// get wrong for no gain. This returns the numbers and the order, and the shim assembles.
pub fn comparison_flow(inp: &ComparisonIn<'_>) -> Result<ComparisonOut, RankNetError> {
    let ncomp = inp.prob.len();
    if ncomp == 0 || inp.prob.len() != inp.pval.len() || inp.prob.len() != inp.pair_names.len() {
        return Err(RankNetError::NoCommunications);
    }
    let k = inp.k;

    // Per comparison: `prob[pval > thresh] <- 0`, `1*(prob > 0)` for count, the two axis
    // filters, the `sum(prob) == 0` stop, then `apply(prob, 3, sum)`.
    let mut original: Vec<Vec<f64>> = Vec::with_capacity(ncomp);
    for i in 0..ncomp {
        let names = &inp.pair_names[i];
        let n = names.len();
        let prob = &inp.prob[i];
        let pval = &inp.pval[i];
        if prob.len() != k * k * n || pval.len() != k * k * n {
            return Err(RankNetError::Shape {
                what: "prob/pval length does not match k * k * n",
            });
        }
        let mut work = vec![0.0f64; k * k * n];
        for t in 0..k * k * n {
            let mut v = if pval[t] > inp.thresh { 0.0 } else { prob[t] };
            if inp.measure == Measure::Count {
                v = if v.is_nan() {
                    f64::NAN
                } else if v > 0.0 {
                    1.0
                } else {
                    0.0
                };
            }
            work[t] = v;
        }
        if let Some(s) = inp.sources_use {
            for c in 0..n {
                for b in 0..k {
                    for a in 0..k {
                        if !s.contains(&a) {
                            work[(c * k + b) * k + a] = 0.0;
                        }
                    }
                }
            }
        }
        if let Some(t) = inp.targets_use {
            for c in 0..n {
                for b in 0..k {
                    for a in 0..k {
                        if !t.contains(&b) {
                            work[(c * k + b) * k + a] = 0.0;
                        }
                    }
                }
            }
        }
        let mut total = F80::from_f64(0.0);
        for v in &work {
            total = total.add(F80::from_f64(*v));
        }
        if total.to_f64() == 0.0 {
            return Err(RankNetError::NoCommunications);
        }
        let mut sums = Vec::with_capacity(n);
        for c in 0..n {
            let mut s = F80::from_f64(0.0);
            for j in 0..k * k {
                s = s.add(F80::from_f64(work[c * k * k + j]));
            }
            sums.push(s.to_f64());
        }
        original.push(sums);
    }

    // The degenerate reassignment is **pooled across comparisons**, not per comparison: one
    // `values.assign` spanning every flagged entry in all of them, ranked by the pooled
    // `pSum.original.all`. Doing it per comparison is the obvious reading and gives a different
    // answer whenever more than one comparison has a degenerate entry.
    let (scaled, flagged) = if inp.measure == Measure::Count {
        (original.clone(), vec![Vec::new(); ncomp])
    } else {
        transform_weight_pooled(&mut original.clone())?
    };

    // `pair.name.all <- as.character(unique(unlist(pair.name)))`: union in first-appearance order.
    let mut pair_names_all: Vec<String> = Vec::new();
    for names in inp.pair_names {
        for nm in names {
            if !pair_names_all.iter().any(|x| x == nm) {
                pair_names_all.push(nm.clone());
            }
        }
    }
    let na = pair_names_all.len();
    let idx_of = |names: &Vec<String>, x: &str| names.iter().position(|y| y == x);

    // `df[[i]] <- data.frame(name = pair.name.all, contribution = 0, ...,
    // row.names = pair.name.all)`, then the **row-name-indexed** assignment
    // `df[[i]][pair.name[[i]], 2] <- pSum.original[[i]]`. Pathways missing from a comparison keep
    // their 0, which is what makes a union across datasets meaningful.
    let contrib: Vec<Vec<f64>> = (0..ncomp)
        .map(|i| {
            let mut v = vec![0.0f64; na];
            for (t, nm) in inp.pair_names[i].iter().enumerate() {
                if let Some(a) = idx_of(&pair_names_all, nm) {
                    v[a] = original[i][t];
                }
            }
            v
        })
        .collect();

    // `contribution.relative[[i]] <- as.numeric(format(df[[ncomp - i + 1]]$contribution /
    // df[[1]]$contribution, digits = 1))`, then `is.na -> 0`. Note the *last* comparison is
    // divided by the *first*, so with three comparisons `relative.1` is `df[[3]]/df[[1]]` and
    // `relative.2` is `df[[2]]/df[[1]]` -- not consecutive pairs.
    //
    // Only the ratio, **unformatted**, is produced here; the `format` and the `order()` that
    // consumes it stay in R. See [`ComparisonOut::ratio`].
    let mut ratio: Vec<Vec<f64>> = Vec::with_capacity(ncomp.saturating_sub(1));
    for j in 1..ncomp {
        let last = &contrib[ncomp - j];
        let first = &contrib[0];
        let mut v: Vec<f64> = (0..na).map(|a| last[a] / first[a]).collect();
        // `NaN` -- a genuine `0/0` -- is what `is.na` catches, and it becomes 0. `Inf` is *not*
        // `is.na`, so a pathway absent from the first comparison keeps its `Inf` and sorts first
        // under `-contribution.relative`. The reset is applied here rather than left to the shim
        // because it is arithmetic, not formatting.
        for x in v.iter_mut() {
            if x.is_nan() {
                *x = 0.0;
            }
        }
        ratio.push(v);
    }

    Ok(ComparisonOut {
        original,
        scaled,
        flagged,
        pair_names_all,
        ratio,
        pair_names_here: inp.pair_names.to_vec(),
    })
}

/// `pSum <- -1/log(pSum); pSum[is.na(pSum)] <- 0; idx <- which(is.infinite | < 0)`, and then the
/// reassignment `pSum[idx] <- values.assign[match(..., position)]` -- **pooled across every
/// entry handed in**, which is what upstream does in comparison mode and what a single-network
/// call degenerates to.
///
/// Pooling is the whole point of the second argument existing. Upstream builds one
/// `values.assign` from `max(unlist(pSum))` over *all* comparisons, one `pSum.original.all` from
/// every comparison's flagged entries concatenated in comparison order, and one `position` from
/// sorting that. Each comparison's slice is then
/// `match(length(unlist(idx[1:i-1])) + 1:length(idx[[i]]), position)` -- an offset into the
/// pooled sequence. Running this per comparison instead gives different assigned values whenever
/// more than one comparison has a degenerate entry, which the `ranknet` comparison fixtures pin.
fn transform_weight_pooled(
    original: &mut [Vec<f64>],
) -> Result<(Vec<Vec<f64>>, Vec<Vec<usize>>), RankNetError> {
    let mut scaled: Vec<Vec<f64>> = original
        .iter()
        .map(|o| o.iter().map(|&v| -1.0 / v.ln()).collect())
        .collect();
    // `pSum[is.na(pSum)] <- 0` runs *before* the degenerate-entry test, so a `NaN` from
    // `log(negative)` becomes an ordinary `0` and is never reassigned.
    for v in scaled.iter_mut() {
        for x in v.iter_mut() {
            if x.is_nan() {
                *x = 0.0;
            }
        }
    }
    // `which(is.infinite(pSum) | pSum < 0)`, in index order within each comparison.
    let flagged: Vec<Vec<usize>> = scaled
        .iter()
        .map(|v| {
            (0..v.len())
                .filter(|&i| v[i].is_infinite() || v[i] < 0.0)
                .collect()
        })
        .collect();

    // `pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])`: the flagged
    // entries' original totals, concatenated in comparison order.
    // Reborrowed explicitly: the closure cannot capture `&mut` and index through it.
    let original: &[Vec<f64>] = original;
    let pooled: Vec<f64> = flagged
        .iter()
        .enumerate()
        .flat_map(|(i, idx)| {
            let row = &original[i];
            idx.iter().map(move |&t| row[t])
        })
        .collect();

    // `seq(max(unlist(pSum))*1.1, max(unlist(pSum))*1.5, length.out = length(unlist(idx)))`.
    //
    // `max()` over a numeric vector excludes `NA` but **not** `-Inf`, so an all-degenerate `pSum`
    // gives `max = -Inf`, and R's `seq()` rejects a non-finite `from` with "'from' must be a
    // finite number". Reproduced rather than worked around: without the guard the sequence comes
    // out as `NaN` and the reassignment silently fills every bar with a missing value.
    let finite_max = scaled
        .iter()
        .flatten()
        .copied()
        .filter(|v| !v.is_nan())
        .fold(f64::NEG_INFINITY, f64::max);
    if !finite_max.is_finite() {
        return Err(RankNetError::SeqFromNotFinite);
    }
    let m = pooled.len();
    let values_assign: Vec<f64> = if m == 0 {
        Vec::new()
    } else if m == 1 {
        // `seq(a, b, length.out = 1)` is `a`, not the midpoint.
        vec![finite_max * 1.1]
    } else {
        let lo = finite_max * 1.1;
        let hi = finite_max * 1.5;
        (0..m)
            .map(|j| lo + (hi - lo) * (j as f64) / ((m - 1) as f64))
            .collect()
    };

    // `position <- sort(pSum.original.all, index.return = TRUE)$ix`, then
    // `values.assign[match(k, position)]` over the pooled slots. `position` is a rank vector:
    // slot -> rank. `sort` on doubles is radix and therefore stable, so tied totals keep their
    // pooled order, and that is what decides which of them gets which assigned value.
    let mut pairs: Vec<(f64, usize)> = pooled
        .iter()
        .copied()
        .enumerate()
        .map(|(i, v)| (v, i))
        .collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap_or(std::cmp::Ordering::Equal));
    let mut rank = vec![0usize; m];
    for (r, &(_, slot)) in pairs.iter().enumerate() {
        rank[slot] = r;
    }

    let mut slot = 0usize;
    for i in 0..scaled.len() {
        for t in 0..scaled[i].len() {
            if flagged[i].contains(&t) {
                scaled[i][t] = values_assign[rank[slot]];
                slot += 1;
            }
        }
    }
    Ok((scaled, flagged))
}

/// R's `order(x)` for a single double vector: radix (therefore **stable**), ascending, `NA`
/// last.
///
/// The stability is the point. `rankNet` sorts by `contribution` and then, for each name,
/// compares `df$name == pair.name[i]`; with a stable sort two pathways with equal totals keep
/// their array order, and an unstable sort can permute them, which changes `df`'s row order and
/// therefore the factor levels of `df$name` and the row names of the result.
pub fn order_f64(x: &[f64]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..x.len()).collect();
    // `na.last = TRUE`: the missing entries go last, among themselves, in index order.
    let key = |i: usize| if x[i].is_nan() { 1 } else { 0 };
    idx.sort_by(|&a, &b| {
        key(a)
            .cmp(&key(b))
            .then_with(|| x[a].partial_cmp(&x[b]).unwrap_or(std::cmp::Ordering::Equal))
    });
    idx
}

/// R's `order(a, b, c)` for several double vectors, `na.last = TRUE`.
///
/// Measured against R:
///
/// ```text
/// order(c(1,1,2),  c(NA,5,3))   -> 2 1 3   # a later key's NA sorts last within its group
/// order(c(1,1,2),  c(5,NA,3))   -> 1 2 3
/// order(c(1,NA,1), c(1,2,3))    -> 1 3 2   # na.last applies to the primary key
/// order(c(NA,NA,NA), c(3,1,2))  -> 2 3 1   # ... and does NOT stop the later keys ordering
/// order(c(1,1,2),  c(NA,NA,3))  -> 1 2 3   # an all-NA subgroup is stable
/// ```
///
/// The last-but-one line is the one that was wrong here. The previous version returned `Equal`
/// as soon as *both* primary keys were missing, which stopped the comparison dead: an all-missing
/// primary key came out in index order instead of being ordered by the remaining keys, and a
/// missing value in a *later* key compared as equal to a finite one rather than sorting after it.
/// For `order(pval, -prob)` that put a `NaN` probability in its original position among equal
/// `pval`s where R moves it to the end of the tie group.
///
/// So: a missing value sorts last in **every** key, and when both sides of a key are missing the
/// comparison falls through to the next key rather than declaring a tie.
pub fn order_f64_multi(keys: &[&[f64]]) -> Vec<usize> {
    let n = keys.first().map_or(0, |k| k.len());
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| {
        for k in keys.iter() {
            let ka = k[a].is_nan();
            let kb = k[b].is_nan();
            if ka != kb {
                return if ka {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Less
                };
            }
            if !ka {
                if let Some(o) = k[a].partial_cmp(&k[b]) {
                    if o != std::cmp::Ordering::Equal {
                        return o;
                    }
                }
            }
        }
        std::cmp::Ordering::Equal
    });
    idx
}

/// `as.numeric(format(x, digits = 1))`.
///
/// `rankNet`'s relative contributions are rounded to **one significant digit** before they are
/// used as sort keys and before they are compared against `1 - tol` / `1 + tol`:
///
/// ```r
/// contribution.relative[[i]] <- as.numeric(format(df[[n-i+1]]$contribution/df[[1]]$contribution, digits=1))
/// ```
///
/// `format(x, digits = 1)` returns *character*, and `as.numeric` parses it back, so the value
/// that survives is the one-significant-digit rounding -- not the ratio. Two ratios that differ
/// in the third digit are therefore **equal** after this line, and a port that keeps full
/// precision produces a different row order. `format` also pads to a common width across the
/// vector, which is why the same number can come back with a different number of characters
/// depending on its neighbours; `as.numeric` discards that, so only the rounding survives.
///
/// `digits = 1` is *significant* digits, not decimal places: `0.0432` becomes `0.04`, and
/// `1234` becomes `1e+03`. `Inf`, `NaN` and `NA` all come back as non-numeric text, and
/// `as.numeric` turns them into `NA` with a warning -- which upstream then silences with
/// `contribution.relative[[i]][is.na(...)] <- 0`.
pub fn format_digits1(x: f64) -> f64 {
    // `format(Inf, digits = 1)` is `"Inf"` and `as.numeric("Inf")` is `Inf`, so an infinite
    // ratio *survives* -- it is `NaN` (and only `NaN`, plus a genuine `NA`) that `as.numeric`
    // turns into `NA` for the caller's `is.na` reset. A zero denominator therefore gives `Inf`
    // here, not `NA`, and a port that mapped both to `NA` would order those rows differently.
    if x.is_nan() || x.is_infinite() {
        return x;
    }
    if x == 0.0 {
        return 0.0;
    }
    // `format(0.0432, digits = 1)` -> "0.04". One significant digit, then back to a double
    // through the *decimal* text, so the rounding is decimal rounding, not a binary scalbn.
    let mag = x.abs().log10().floor() as i32;
    let scale = 10f64.powi(1 - mag);
    // ` signif(x, 1) ` is what R formats, and it rounds half to even on the scaled value.
    let r = (x * scale).round() / scale;
    // Re-round through the shortest decimal text that R would print, so `1e+03` and `0.04`
    // come out as the doubles R parses rather than as binary approximations of them.
    let s = format!("{:.0e}", r);
    s.parse().unwrap_or(r)
}

#[cfg(test)]
mod order_multi_nan_tests {
    use super::order_f64_multi;

    /// Each expectation below is R's own answer for `order(a, b)` with `na.last = TRUE`,
    /// measured rather than reasoned -- with **1 subtracted**, because this function returns
    /// 0-based indices and the binding adds the 1 back on the way out.
    #[test]
    fn matches_r_order_with_missing_in_later_keys() {
        let nan = f64::NAN;
        // A later key's NA sorts last within its group of equal primary keys.
        let a = [1.0, 1.0, 2.0];
        assert_eq!(order_f64_multi(&[&a, &[nan, 5.0, 3.0]]), vec![1, 0, 2]);
        assert_eq!(order_f64_multi(&[&a, &[5.0, nan, 3.0]]), vec![0, 1, 2]);
        // `na.last` applies to the primary key.
        let b = [1.0, nan, 1.0];
        assert_eq!(order_f64_multi(&[&b, &[1.0, 2.0, 3.0]]), vec![0, 2, 1]);
        // ... but it does not stop the remaining keys from ordering an all-missing primary group.
        // The previous version returned index order here.
        let all = [nan, nan, nan];
        assert_eq!(order_f64_multi(&[&all, &[3.0, 1.0, 2.0]]), vec![1, 2, 0]);
        // Two rows missing in *both* keys are a tie, so the stable order is index order.
        assert_eq!(order_f64_multi(&[&b, &[nan, 2.0, nan]]), vec![0, 2, 1]);
        // An all-missing later key leaves the group stable.
        assert_eq!(order_f64_multi(&[&a, &[nan, nan, 3.0]]), vec![0, 1, 2]);
    }

    /// The exact case that `rankNetPairwise` depends on: equal `pval`, one `NaN` probability.
    #[test]
    fn nan_probability_sorts_after_equal_pvals() {
        let pval = [0.49, 0.49, 0.79];
        let neg = [f64::NAN, -0.96, -0.81];
        assert_eq!(order_f64_multi(&[&pval, &neg]), vec![1, 0, 2]);
    }
}
