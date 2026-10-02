// --------------------------------------------------------------------------------- filter.rs
//
// `filterCommunication`.
//
// Two independent halves, and the second one is nearly all of the work:
//
// 1. **`min.cells`** (always active): zero `net$prob` on the rows and columns of cell groups
//    with at most `min.cells` cells, and report the percentage of positive entries removed.
//
// 2. **`min.samples >= 2`** (cross-sample consistency): for each sample, recompute the group
//    means of the ligand and receptor genes, score every L-R pair as the outer product of
//    its ligand and receptor profiles, binarise, count in how many samples each
//    (source, target) is present, and keep only those present in at least `min.samples`.
//
// # R-isms
//
// * **`apply(score.LR[, , jj, ], c(1, 2), sum)` collapses the *sample* axis.** `d.call` is
//   the non-`MARGIN` dimension, so the array is reshaped to `(k*k) x n_samples` with the row
//   index `source + k*target` and the columns in sample order, and `sum` runs over the
//   columns. Summing over samples in a different order (say all of sample 1's cells, then
//   sample 2's) gives a different answer, because the values are 0/1 counts and the order
//   changes only the *rounding* -- but it is a different expression, and it is the axis that
//   matters.
// * **`droplevels(group.use)` inside the per-sample loop** changes the *group order* per
//   sample: a sample containing only `g2` and `g5` produces a 2-column aggregate in *its*
//   level order, which then has to be scattered back into a `k`-column matrix via
//   `group.exist <- which(levels(group) %in% unique(group.use))`. Doing the zero-fill once,
//   outside the loop, or using the level order of the *full* factor, gives a different
//   matrix -- and a plausible one, because the values are all non-negative.
// * **`table(object@idents[cell.use])` keeps every level of the full factor.** Subsetting a
//   factor does not drop its levels, so `table(...)` has `k` entries and a group *absent*
//   from a sample counts as 0 cells and is therefore "excluded" for that sample. This is
//   what makes `cell.excludes.sample` include groups that are not in the sample at all, and
//   it is the reason `rare.keep` has anything to preserve.
// * **`cell.excludes.sample <- setdiff(cell.excludes.sample, cell.excludes)` at the end**, so
//   a group already dropped for having too few cells *overall* is not also treated as rare.
// * **`if ("data.smooth" %in% methods::slotNames(object) == FALSE)`** -- `%in%` binds tighter
//   than `==`, so it is `(x %in% y) == FALSE`, which is what was meant. It reads as a
//   precedence bug and happens to be harmless.
// * **`for (jj in 1:length(LR.nonzero))`** with `LR.nonzero` empty is `1:0`, i.e. `c(1, 0)`,
//   and `score.LR[, , 0, i] <- ...` raises **"subscript out of bounds"**. So an all-zero
//   `net$prob` combined with `min.samples >= 2` crashes upstream, on the very first
//   iteration. Reproduced as [`FilterError::EmptyLrNonzero`]; the R shim turns it back into
//   upstream's message.
// * **`LR.nonzero <- LR[which(apply(net$prob, 3, sum) != 0)]`** is computed *after* the
//   `min.cells` zeroing, so the two halves interact: a pair that `min.cells` removed is also
//   removed from the per-sample scoring.
// * **`geneLR <- extractGeneSubset(geneLR, complex_input, geneIfo)`** reduces to
//   `unique(c(intersect(geneSet, geneIfo$Symbol), <subunits of the names absent from
//!   Symbol>))`, because `which(geneSet %in% geneIfo$Symbol == "FALSE")` coerces the logical
//   to character and so really asks "is the name *absent*?". Benign for a well-formed
//!   database, since the names absent from `Symbol` are exactly the complex names -- but it
//!   means `geneInfo$Symbol` must not list the complexes.

use crate::aggregate::{aggregate_1, all_cols, GroupMean};
use crate::db::Database;
use crate::expr::{compute_expr_lr, GroupMeans};
use crate::longdouble::F80;

/// `object@options$parameter$type.mean`, which selects the per-group mean.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum MeanKind {
    TriMean,
    TruncatedMean,
    ThresholdedMean,
    Median,
}

impl MeanKind {
    fn to_group_mean(self, trim: f64) -> GroupMean {
        match self {
            MeanKind::TriMean => GroupMean::TriMean,
            // `mean(x, trim = trim, na.rm = TRUE)`. NB `trim = NULL` is *never* substituted
            // upstream, so a `NULL` trim means R's `trim = 0` -- no trimming -- and the
            // kernel's own 0.1 default is a different function. The shim passes 0 in that
            // case; see `computeAveExpr` for the same trap.
            MeanKind::TruncatedMean => GroupMean::TrimmedMean { trim },
            MeanKind::ThresholdedMean => GroupMean::ThresholdedMean { trim },
            MeanKind::Median => GroupMean::Median,
        }
    }
}

/// Flat index into `net$prob`, an R `k x k x nLR` array.
///
/// **Last dimension fastest**: `prob[source, target, lr]` is at
/// `source + k * target + k * k * lr`. The crate has two conventions in circulation and
/// picking the wrong one is invisible here -- every index stays in range, every value stays
/// plausible, and the corpus fails on a handful of entries. Hence a named function rather
/// than arithmetic at each of the dozen use sites.
///
/// `n_lr` is unused (the last axis is addressed by `lr` alone) but kept in the signature so a
/// caller cannot pass a `k` from a different array and get a silently in-range index.
#[inline]
fn pi(source: usize, target: usize, lr: usize, k: usize, _n_lr: usize) -> usize {
    source + k * target + k * k * lr
}

/// Flat index into `score.LR`, an R `k x k x n_nz x n_samples` array. Same convention.
#[inline]
fn si(source: usize, target: usize, jj: usize, sample: usize, k: usize, n_nz: usize) -> usize {
    source + k * target + k * k * jj + k * k * n_nz * sample
}

/// The upstream failures `filterCommunication` can raise, reproduced rather than repaired.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum FilterError {
    /// `stop(paste0("There are only ", n, " samples in the data. Please change the value of
    /// \`min.samples\`! "))`. **Not raised by [`filter_communication`]** -- it reports
    /// `FilterResult::min_samples_ok = false` instead, because upstream raises it after the
    /// `min.cells` messages have already been printed. This variant exists so a standalone
    /// Rust caller can produce the same text.
    MinSamplesTooLarge {
        n_samples: usize,
        min_samples: usize,
    },
    /// `score.LR[, , 0, i] <- ...` from `for (jj in 1:length(LR.nonzero))` with
    /// `LR.nonzero` empty. Upstream's message is "subscript out of bounds".
    EmptyLrNonzero,
    /// `!is.list(object@var.features)`-style preconditions the shim checks in R.
    MissingDataSmooth,
    /// A ligand or receptor could not be resolved against the expression matrix and the
    /// complex table. Upstream raises `subscript out of bounds` from inside
    /// `computeExpr_LR`; that is the same message, from a different cause.
    UnresolvedEntity { name: String },
}

impl std::fmt::Display for FilterError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            // Upstream's exact text, trailing space included: `paste0("There are only ",
            // length(sample.id), " samples in the data. Please change the value of
            // `min.samples`! ")`. It does *not* mention the value the caller passed, which is
            // the first thing anyone debugging this wants.
            FilterError::MinSamplesTooLarge {
                n_samples,
                min_samples: _,
            } => write!(
                f,
                "There are only {n_samples} samples in the data. Please change the value of \
                 `min.samples`! "
            ),
            // Upstream's own text, byte for byte: it comes from R's subscript handler, not
            // from CellChat, so a port that improves the wording is *less* faithful.
            FilterError::EmptyLrNonzero => write!(f, "subscript out of bounds"),
            FilterError::UnresolvedEntity { .. } => write!(f, "subscript out of bounds"),
            FilterError::MissingDataSmooth => write!(
                f,
                "`object@data.smooth` is missing. Please update the CellChat object via \
                 `updateCellChat`! "
            ),
        }
    }
}

/// What `filterCommunication` computes. The `cat()` text is the shim's job -- it needs
/// `scales::percent()` -- so only the numbers and the index sets come from here.
/// R's `NA_real_`, exactly: quiet NaN, sign 0, payload `0x7a2`.
/// `identical(NA_real_, NaN)` is FALSE in R, so the bit pattern -- not `is_nan()` -- is what
/// distinguishes them, and every comparison below goes through this predicate rather than
/// `f64::is_nan`, which conflates the two.
pub const NA_REAL_BITS: u64 = 0x7ff0_0000_0000_07a2;

#[inline]
pub fn is_na(v: f64) -> bool {
    v.to_bits() == NA_REAL_BITS
}

/// R's `sum(x > 0)` with the default `na.rm = FALSE`, as a count.
///
/// Elementwise `x > 0` yields `NA` for both `NA_real_` and `NaN` operands (R's comparison
/// operators propagate missingness rather than returning FALSE the way IEEE does), and `sum()`
/// over a vector containing any `NA` is `NA_real_` -- not a number. So the result is `None`
/// (R's `NA`) iff any element is `NA_real_` *or* `NaN`, and `Some(count)` otherwise. `Inf`
/// counts (`Inf > 0` is TRUE) and `-Inf` does not, exactly as in R.
fn count_positive(xs: &[f64]) -> Option<usize> {
    let mut n = 0usize;
    for &v in xs {
        if v.is_nan() {
            return None;
        }
        if v > 0.0 {
            n += 1;
        }
    }
    Some(n)
}

/// R's `sum()` over one L-R slice with `na.rm = FALSE`, for the `which()` in
/// `LR.nonzero <- LR[which(apply(net$prob, 3, sum) != 0)]`.
///
/// Returns `None` iff the slice contains an `NA_real_` element (R's sum is `NA`, `NA != 0` is
/// `NA`, `which()` drops it). Any other non-finite element goes through the existing `F80`
/// accumulation: R's sum of a `NaN`-containing vector is `NaN`, and `NaN != 0` is TRUE, so the
/// slice is *included* -- which is what `acc.to_f64() != 0.0` already decides, since a NaN
/// accumulator compares unequal to everything including itself.
fn slice_sum(slice: impl Iterator<Item = f64>) -> Option<f64> {
    let mut acc = F80::ZERO;
    for v in slice {
        if is_na(v) {
            return None;
        }
        acc = acc.add(F80::from_f64(v));
    }
    Some(acc.to_f64())
}

pub struct FilterResult {
    /// The filtered `net$prob`, `k * k * nLR` with the last dimension fastest.
    pub prob: Vec<f64>,
    /// `sum(net$prob > 0)` before any filtering -- `None` (R's `NA`) iff the buffer holds any
    /// `NA_real_` or `NaN`, per `count_positive`.
    pub num_interaction0: Option<usize>,
    /// ... after the `min.cells` half.
    pub num_interaction1: Option<usize>,
    /// ... after the `min.samples` half. Equal to `num_interaction1` when that half is
    /// skipped, which is what upstream leaves it as.
    pub num_interaction2: Option<usize>,
    /// `which(as.numeric(table(object@idents)) <= min.cells)`, 0-based, ascending.
    pub cell_excludes: Vec<usize>,
    /// Per sample, `which(as.numeric(table(object@idents[cell.use])) <= min.cells)`, 0-based.
    /// Includes groups absent from that sample, whose count is 0.
    pub sample_excluded: Vec<Vec<usize>>,
    /// `cell.excludes.sample` after `unique()` and the final `setdiff(..., cell.excludes)`.
    pub cell_excludes_sample: Vec<usize>,
    /// Indices (0-based) of the L-R pairs with a non-zero total, ascending. Empty when
    /// [`FilterError::EmptyLrNonzero`] is returned.
    pub lr_nonzero: Vec<usize>,
    /// `score.LR`, `k * k * len(lr_nonzero) * n_samples`, last dimension fastest. Exposed
    /// because the binarised scores are the only place the per-sample structure is visible,
    /// and a corpus that only checked `net$prob` could not tell a wrong `score.LR` that
    /// happened to cancel.
    pub score: Vec<f64>,
    /// Whether the `min.samples` half ran at all.
    pub ran_sample_filter: bool,
    /// `min.samples <= n_samples`. Upstream raises this *after* the `min.cells` messages, so
    /// the kernel reports it and the shim raises -- otherwise the messages would be missing
    /// from the error path, which `capture.output` shows and `identical()` on the object does
    /// not, so nothing else would notice.
    pub min_samples_ok: bool,
}

pub struct FilterInput<'a> {
    /// `k * k * nLR`, last dimension fastest.
    pub prob: &'a [f64],
    pub n_lr: usize,
    /// `levels(object@idents)`.
    pub n_groups: usize,
    /// Per cell, the level index of `object@idents`.
    pub group_index: &'a [usize],
    pub n_samples: usize,
    /// Per cell, the level index of `object@meta$samples`.
    pub sample_index: &'a [usize],
    /// `object@data.signaling` as a genes x cells column-major buffer, already
    /// `data/max(data)`-scaled by the caller so the value is inspectable.
    pub data: &'a [f64],
    pub gene_names: &'a [String],
    pub min_cells: usize,
    pub min_samples: Option<usize>,
    pub rare_keep: bool,
    pub mean: MeanKind,
    pub trim: f64,
    /// The ligand and receptor of each entry of `prob`'s L-R axis, parallel to it.
    pub ligand: &'a [String],
    pub receptor: &'a [String],
}

/// `filterCommunication`'s numerics. See the module docs for the R-isms.
pub fn filter_communication(
    inp: &FilterInput<'_>,
    db: &Database,
) -> Result<FilterResult, FilterError> {
    let k = inp.n_groups;
    let n_lr = inp.n_lr;
    let kk = k * k;
    assert_eq!(
        inp.prob.len(),
        kk * n_lr,
        "filter: prob buffer must be k*k*nLR"
    );
    assert_eq!(inp.ligand.len(), n_lr, "filter: one ligand per L-R");
    assert_eq!(inp.receptor.len(), n_lr, "filter: one receptor per L-R");
    assert_eq!(
        inp.group_index.len(),
        inp.data.len() / inp.gene_names.len().max(1)
    );

    let mut prob = inp.prob.to_vec();
    let num_interaction0 = count_positive(&prob);

    // ---- half 1: min.cells -------------------------------------------------------------
    let mut counts = vec![0usize; k];
    for &g in inp.group_index {
        counts[g] += 1;
    }
    let cell_excludes: Vec<usize> = (0..k).filter(|&g| counts[g] <= inp.min_cells).collect();
    if !cell_excludes.is_empty() {
        // `net$prob[cell.excludes, , ] <- 0` and `net$prob[, cell.excludes, ] <- 0`.
        //
        // Four axes, four loops. The first version of this had three, with the *source*
        // index in the L-R slot -- so it zeroed L-R pairs 0..k-1 of every (g, t) and left
        // the rest, which is a plausible-looking net that happens to be right for a fixture
        // with `n_lr <= k` and wrong for every real one. The corpus caught it at
        // `net$prob[27]` of 54; `n1` was 27 instead of 0.
        for &g in &cell_excludes {
            for t in 0..k {
                for s in 0..k {
                    for l in 0..n_lr {
                        prob[pi(g, t, l, k, n_lr)] = 0.0;
                        prob[pi(s, g, l, k, n_lr)] = 0.0;
                    }
                }
            }
        }
    }
    let num_interaction1 = count_positive(&prob);
    let mut out = FilterResult {
        prob: prob.clone(),
        num_interaction0,
        num_interaction1,
        num_interaction2: num_interaction1,
        cell_excludes,
        sample_excluded: Vec::new(),
        cell_excludes_sample: Vec::new(),
        lr_nonzero: Vec::new(),
        score: Vec::new(),
        ran_sample_filter: false,
        min_samples_ok: true,
    };

    // ---- half 2: min.samples -----------------------------------------------------------
    let min_samples = inp.min_samples.unwrap_or(1);
    if inp.n_samples == 0 || min_samples > inp.n_samples {
        // Reported, not raised: see `min_samples_ok`. The caller turns it into upstream's
        // `stop()` at the right point in the message sequence.
        out.min_samples_ok = false;
        return Ok(out);
    }
    if !(inp.n_samples >= 2 && min_samples >= 2) {
        return Ok(out);
    }
    out.ran_sample_filter = true;

    // `LR.nonzero <- LR[which(apply(net$prob, 3, sum) != 0)]`, on the *already filtered*
    // prob. LONG_DOUBLE accumulation, as R's `sum` does -- except that a slice containing an
    // `NA_real_` element sums to `NA` and is dropped by `which()`, which `slice_sum` returns as
    // `None` for. (The pre-NA code accumulated every slice through `F80`, so an `NA`-containing
    // slice produced an `F80::NAN` that compared `!= 0.0` and was *included* -- the opposite of
    // R. No fixture covered it because no fixture's `prob` held an `NA`; the tutorial's real
    // object was the first.)
    let mut lr_nonzero = Vec::new();
    for l in 0..n_lr {
        // Collected first so the borrow of `prob` ends before the push; a closure over the
        // buffer would capture it across the loop body instead.
        let mut vals = Vec::with_capacity(kk);
        for t in 0..k {
            for s in 0..k {
                vals.push(prob[pi(s, t, l, k, n_lr)]);
            }
        }
        let sum = slice_sum(vals.into_iter());
        let Some(sum) = sum else { continue };
        if sum != 0.0 {
            lr_nonzero.push(l);
        }
    }
    // `for (jj in 1:length(LR.nonzero))` with an empty `LR.nonzero` is `1:0`, and
    // `score.LR[, , 0, i] <- ...` is a subscript error. Reproduced.
    if lr_nonzero.is_empty() {
        return Err(FilterError::EmptyLrNonzero);
    }
    let n_nz = lr_nonzero.len();
    let n_samples = inp.n_samples;

    // `geneLR <- c(unique(geneL), unique(geneR))`, then `extractGeneSubset`.
    let mut gene_lr: Vec<String> = Vec::new();
    for l in &lr_nonzero {
        if !gene_lr.contains(&inp.ligand[*l]) {
            gene_lr.push(inp.ligand[*l].clone());
        }
    }
    for l in &lr_nonzero {
        if !gene_lr.contains(&inp.receptor[*l]) {
            gene_lr.push(inp.receptor[*l].clone());
        }
    }
    let gene_lr = db.extract_gene_subset(&gene_lr);
    // `data.use <- data.use[rownames(data.use) %in% geneLR, ]` -- `%in%` keeps the matrix's
    // own row order.
    let keep_genes: Vec<usize> = (0..inp.gene_names.len())
        .filter(|&g| gene_lr.iter().any(|x| x == &inp.gene_names[g]))
        .collect();
    let gene_names: Vec<String> = keep_genes
        .iter()
        .map(|&g| inp.gene_names[g].clone())
        .collect();
    // genes x cells, the subset only.
    let n_cells = inp.data.len() / inp.gene_names.len().max(1);
    let mut data_use = vec![0.0f64; keep_genes.len() * n_cells];
    for (new_g, &old_g) in keep_genes.iter().enumerate() {
        for c in 0..n_cells {
            data_use[c * keep_genes.len() + new_g] = inp.data[c * inp.gene_names.len() + old_g];
        }
    }

    let gene_l: Vec<String> = lr_nonzero.iter().map(|&l| inp.ligand[l].clone()).collect();
    let gene_r: Vec<String> = lr_nonzero
        .iter()
        .map(|&l| inp.receptor[l].clone())
        .collect();

    let mut score = vec![0.0f64; kk * n_nz * n_samples];
    let mut sample_excluded: Vec<Vec<usize>> = Vec::with_capacity(n_samples);
    for i in 0..n_samples {
        let cell_use: Vec<usize> = (0..n_cells).filter(|&c| inp.sample_index[c] == i).collect();
        // `group.use <- droplevels(group[cell.use])`: the levels present in this sample, in
        // the full factor's level order. Everything downstream keys off this list.
        let present: Vec<usize> = (0..k)
            .filter(|&g| cell_use.iter().any(|&c| inp.group_index[c] == g))
            .collect();
        let n_pres = present.len();
        // `droplevels()` *renumbers*: the aggregated factor has levels `present` in order and
        // its codes are 1..=n_pres. Passing the full factor's codes with `n_pres` levels
        // indexes past the end the moment a sample is missing a group, which is the common
        // case for a rare population.
        let group_use: Vec<usize> = cell_use
            .iter()
            .map(|&c| {
                present
                    .iter()
                    .position(|&g| g == inp.group_index[c])
                    .expect("every cell's group is in `present` by construction")
            })
            .collect();

        // `cell.excludes.sample.i <- which(as.numeric(table(idents[cell.use])) <= min.cells)`.
        // `table()` on a subsetted factor keeps every level, so an absent group counts 0 and
        // is included. That is what gives `rare.keep` something to preserve.
        let mut per_sample = vec![0usize; k];
        for &c in &cell_use {
            per_sample[inp.group_index[c]] += 1;
        }
        let excl_i: Vec<usize> = (0..k).filter(|&g| per_sample[g] <= inp.min_cells).collect();
        sample_excluded.push(excl_i.clone());

        // `data.use.avg <- t(aggregate(t(data.use[, cell.use]), list(group.use), FUN)[, -1])`
        let n_genes = data_use.len() / n_cells.max(1);
        let mut cells_x_genes = vec![0.0f64; cell_use.len() * n_genes];
        for (new_c, &old_c) in cell_use.iter().enumerate() {
            for g in 0..n_genes {
                cells_x_genes[g * cell_use.len() + new_c] = data_use[old_c * n_genes + g];
            }
        }
        let gmean = inp.mean.to_group_mean(inp.trim);
        // `aggregate_1` returns n_pres x n_genes with flat index `group * n_cols + gene`.
        let agg = aggregate_1(
            &cells_x_genes,
            cell_use.len(),
            &all_cols(n_genes),
            &group_use,
            n_pres,
            gmean,
        );
        // `t(agg[, -1])` -> genes x groups-present, flat `gene * n_pres + group`.
        let mut avg = vec![0.0f64; n_genes * n_pres.max(1)];
        for g in 0..n_genes {
            for j in 0..n_pres {
                avg[g * n_pres + j] = agg[j * n_genes + g];
            }
        }
        // The zero-fill back to `levels(group)`. `group.exist` is in level order, so
        // `present` (already in level order) is exactly it.
        //
        // Two transposes meet here and both are silent. `avg` is genes x groups-present
        // (gene-major, from `t(agg[, -1])`), while `GroupMeans` takes a **group-major**
        // buffer -- `values[group * n_genes + gene]` -- so the scatter has to flip as well as
        // place. Feeding `GroupMeans` the gene-major buffer it already had in the first
        // version produced plausible non-zero averages for the wrong groups, and the only
        // symptom was a handful of `net$prob` entries that the consistency mask zeroed
        // because one sample looked absent.
        let mut avg_k = vec![0.0f64; n_genes * k];
        for g in 0..n_genes {
            for (j, &lev) in present.iter().enumerate() {
                avg_k[lev * n_genes + g] = avg[g * n_pres + j];
            }
        }
        let means = GroupMeans::new(&avg_k, &gene_names, k);
        let data_l =
            compute_expr_lr(&gene_l, &means, db).map_err(|e| FilterError::UnresolvedEntity {
                name: e.to_string(),
            })?;
        let data_r =
            compute_expr_lr(&gene_r, &means, db).map_err(|e| FilterError::UnresolvedEntity {
                name: e.to_string(),
            })?;

        for (jj, _) in lr_nonzero.iter().enumerate() {
            // `score.LR[,,jj,i] <- Matrix::crossprod(matrix(a, nrow = 1), matrix(b, nrow = 1))`
            // which is the outer product a[s] * b[t].
            //
            // `compute_expr_lr` returns an `nLR x n_groups` matrix flattened **group-major**
            // (`out[group * n + entity]`), so the profile of pair `jj` for group `g` is at
            // `g * n_nz + jj`. Reading it as `jj * k + g` gives the *transpose* of every
            // score, which is not a small error: the tensor stays in range, every value stays
            // positive, and a fixture where the ligand and receptor profiles happen to be
            // symmetric in two of three groups still agrees. The `two_samples_ms2` case caught
            // it because the excluded group is a *source* in one sample and a *target* in the
            // other, so the transpose moved the zeroed row onto the wrong axis.
            for t in 0..k {
                for s in 0..k {
                    score[si(s, t, jj, i, k, n_nz)] = data_l[s * n_nz + jj] * data_r[t * n_nz + jj];
                }
            }
        }
        for &g in &excl_i {
            for t in 0..k {
                for s in 0..k {
                    for jj in 0..n_nz {
                        score[si(g, t, jj, i, k, n_nz)] = 0.0;
                        score[si(s, g, jj, i, k, n_nz)] = 0.0;
                    }
                }
            }
        }
    }
    out.sample_excluded = sample_excluded;

    // `cell.excludes.sample <- unique(...)` then `setdiff(..., cell.excludes)`.
    let mut ex: Vec<usize> = out.sample_excluded.iter().flatten().copied().collect();
    ex.sort_unstable();
    ex.dedup();
    ex.retain(|g| !out.cell_excludes.contains(g));
    out.cell_excludes_sample = ex;

    // `score.LR[score.LR > 0] <- 1`
    for v in score.iter_mut() {
        if *v > 0.0 {
            *v = 1.0;
        }
    }

    // Per L-R pair, count the samples in which each (source, target) is present.
    for (jj, &l) in lr_nonzero.iter().enumerate() {
        let mut sum = vec![0.0f64; kk];
        for t in 0..k {
            for s in 0..k {
                // `apply(score.LR[,,jj,], c(1,2), sum)` collapses the *sample* axis, so the
                // accumulation is over `i` in ascending sample order. The values are 0/1
                // counts, so this is exact in f64 -- but the axis is the point, and summing
                // the wrong one is a different expression, not a rounding difference.
                let mut acc = F80::ZERO;
                for i in 0..n_samples {
                    acc = acc.add(F80::from_f64(score[si(s, t, jj, i, k, n_nz)]));
                }
                sum[t * k + s] = acc.to_f64();
            }
        }
        // `if (sum((score.LR.sum > 0) * (score.LR.sum < min.samples)) > 0)`
        let n_partial = sum
            .iter()
            .filter(|v| **v > 0.0 && **v < min_samples as f64)
            .count();
        if n_partial == 0 {
            continue;
        }
        // `score.LR.consitent <- (score.LR.sum >= min.samples) * 1`
        let consitent: Vec<f64> = sum
            .iter()
            .map(|v| if *v >= min_samples as f64 { 1.0 } else { 0.0 })
            .collect();
        let mut consitent = consitent;
        if inp.rare_keep && !out.cell_excludes_sample.is_empty() {
            for &g in &out.cell_excludes_sample {
                for t in 0..k {
                    consitent[t * k + g] = 1.0;
                    consitent[g * k + t] = 1.0;
                }
            }
        }
        for t in 0..k {
            for s in 0..k {
                prob[pi(s, t, l, k, n_lr)] *= consitent[t * k + s];
            }
        }
    }
    out.num_interaction2 = count_positive(&prob);
    out.prob = prob;
    out.lr_nonzero = lr_nonzero;
    out.score = score;
    Ok(out)
}

#[cfg(test)]
mod na_tests {
    use super::{count_positive, is_na, slice_sum, NA_REAL_BITS};

    const NA: f64 = f64::from_bits(NA_REAL_BITS);

    #[test]
    fn na_real_is_detected_and_plain_nan_is_not_na() {
        assert!(is_na(NA));
        assert!(is_na(f64::from_bits(NA_REAL_BITS)));
        assert!(!is_na(f64::NAN));
        assert!(!is_na(f64::NEG_INFINITY));
        assert!(!is_na(0.0));
        // `identical(NA_real_, NaN)` is FALSE in R; the predicate must agree.
        assert_ne!(NA.to_bits(), f64::NAN.to_bits());
    }

    #[test]
    fn counts_are_na_when_any_element_is_na_or_nan() {
        // R: `sum(c(TRUE, NA))` is NA; `sum(c(TRUE, NaN > 0))` is NA because `NaN > 0` is NA.
        assert_eq!(count_positive(&[1.0, 2.0]), Some(2));
        assert_eq!(count_positive(&[1.0, 0.0, -1.0]), Some(1));
        assert_eq!(count_positive(&[1.0, NA]), None);
        assert_eq!(count_positive(&[1.0, f64::NAN]), None);
        assert_eq!(count_positive(&[0.0, NA]), None);
        // Infinities are ordinary values to `>`: `Inf > 0` is TRUE, `-Inf > 0` is FALSE.
        assert_eq!(count_positive(&[f64::INFINITY]), Some(1));
        assert_eq!(count_positive(&[f64::NEG_INFINITY]), Some(0));
        assert_eq!(count_positive(&[]), Some(0));
    }

    #[test]
    fn slice_sums_drop_na_slices_and_keep_nan_ones() {
        // `which(apply(prob, 3, sum) != 0)`: an NA-summing slice is NA, and `which()` drops NAs.
        assert_eq!(slice_sum([1.0, 2.0].into_iter()), Some(3.0));
        assert!(matches!(slice_sum([0.0, 0.0].into_iter()), Some(x) if x == 0.0));
        assert_eq!(slice_sum([1.0, NA].into_iter()), None);
        // R's `sum(c(1, NaN))` is NaN (not NA), and `NaN != 0` is TRUE: included.
        let nan_sum = slice_sum([1.0, f64::NAN].into_iter());
        assert!(matches!(nan_sum, Some(x) if x.is_nan() && !is_na(x)));
    }
}
