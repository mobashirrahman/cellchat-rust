//! `computeCommunProb` / `computeCommunProbPathway`: the bootstrap permutation test that
//! produces `net$prob` and `net$pval`.
//!
//! Parity contract: **Exact** (`docs/SEMANTICS.md`).
//!
//! Upstream: `R/modeling.R:114` (`computeCommunProb`) and `R/modeling.R:259`
//! (`computeCommunProbPathway`, which differs only in which objects it returns).
//!
//! This is 92–95 % of CellChat's runtime on the authors' datasets, so the design goal is
//! both bit-exactness and doing the per-`(interaction, bootstrap)` work without ever
//! touching a name. Every string resolution happens once, in [`LrPlan::resolve`]; the
//! hot loops then work entirely in `usize` gene indices.
//!
//! ## The loop-carried `P.spatial` mutation
//!
//! `R/modeling.R:232` does `if (i > nLR1) P.spatial <- P.spatial * adj.contact` **inside**
//! the interaction loop, and `P.spatial` is what `Pnull` is built from. Reading it
//! literally, `adj.contact` is applied once per iteration past `nLR1`.
//!
//! It is idempotent only because `adj.contact` is 0/1 (and all-ones for RNA), so
//! `x * adj * adj == x * adj`. [`Kernel::prob_arrays`] therefore hoists it to
//! `P_eff(i) = P0 * adj` for `i > nLR1`, which is *not* the same as the value used for
//! `i <= nLR1` and is observable in `Prob`. What is not observable: `P1_Pspatial` in the
//! `sum(...) == 0` early-exit branch is computed *before* the mutation of the same
//! iteration, so for `i > nLR1` it already sees the mutated matrix. Both branches
//! therefore use the same `P_eff(i)`, and the hoisted form is exact.
//!
//! ## The `sum(P1_Pspatial) == 0` early exit
//!
//! When it fires, `Prob[,,i] <- P1_Pspatial`, `Pval[,,i] <- 1` -- note `1`, not
//! `1/nboot`, and *not* the R-ism 7 bootstrap path. So a pair whose observed network is
//! identically zero is never resampled.

use std::collections::HashMap;

use crate::aggregate::{aggregate_1, all_cols, GroupMean};
use crate::db::Database;
#[cfg(test)]
use crate::expr::{compute_expr_coreceptor, compute_expr_lr, CoreceptorKind};
use crate::expr::{EntityRows, ExprError, GroupMeans};
use crate::longdouble::F80;
use crate::rng::MersenneTwister;
use crate::stats::r_prod;

/// The `type.mean` values `match.arg` accepts, in the order upstream declares them, so
/// partial matching behaves the same. R-ism 9.
pub const MEAN_TYPES: [&str; 4] = ["triMean", "truncatedMean", "thresholdedMean", "median"];

/// `match.arg(type)`: exact match, else a unique partial match, else R's error text.
pub fn match_mean_type(type_: &str) -> Result<GroupMean, KernelError> {
    for t in MEAN_TYPES {
        if t == type_ {
            return Ok(mean_for(t, 0.1));
        }
    }
    let hits: Vec<&&str> = MEAN_TYPES.iter().filter(|t| t.starts_with(type_)).collect();
    match hits.len() {
        1 => Ok(mean_for(hits[0], 0.1)),
        0 => Err(KernelError::MatchArg(type_.to_string())),
        _ => Err(KernelError::MatchArg(type_.to_string())),
    }
}

fn mean_for(t: &str, trim: f64) -> GroupMean {
    match t {
        "triMean" => GroupMean::TriMean,
        "truncatedMean" => GroupMean::TrimmedMean { trim },
        "thresholdedMean" => GroupMean::ThresholdedMean { trim },
        _ => GroupMean::Median,
    }
}

/// Upstream's `aggregate` applies `FunMean` to **every gene**, not just the ones the
/// interaction list mentions, and for `thresholdedMean` a single missing value makes the function
/// raise rather than return a number:
///
/// ```text
/// c(NaN,1,3,5) trim=0.5 : ERROR: missing value where TRUE/FALSE needed
/// c(0,0,0,0)   trim=0.1 : 0
/// ```
///
/// The port aggregates all genes too, so it *has* the offending `NaN` -- and then discards it,
/// because a gene no interaction references is never read again. The result was an all-zero `Prob`
/// where upstream refuses to produce one. `thresholded_mean` returns `NaN` in exactly the case
/// where upstream raised (that is its only route to `NaN`), so a `NaN` in the observed aggregate
/// is a faithful raise signal for this one `FUN`.
///
/// The check is deliberately restricted to `thresholdedMean`. For `triMean` and `median` a `NaN`
/// group mean is an ordinary *value*, not an error: upstream's `if (sum(P1_Pspatial) == 0)` is what
/// turns it into `missing value where TRUE/FALSE needed`, and that is handled where it happens.
pub fn check_observed_means(fun: GroupMean, avg: &[f64]) -> Result<(), KernelError> {
    if !matches!(fun, GroupMean::ThresholdedMean { .. }) {
        return Ok(());
    }
    if avg.iter().any(|v| v.is_nan()) {
        return Err(KernelError::IfNa);
    }
    Ok(())
}

#[derive(Debug, Clone, PartialEq)]
pub enum KernelError {
    /// R's `match.arg` failure, with upstream's message.
    MatchArg(String),
    /// `numCluster != length(unique(group))` -- upstream's `stop()`.
    UnusedLevels,
    /// A ligand/receptor that cannot be resolved: upstream aborts with
    /// `subscript out of bounds` (R-ism 11).
    Expr(ExprError),
    /// `max(data)` is non-finite, so `data/max(data)` is all `NaN`. Upstream warns
    /// (`-Inf`) and continues with an all-zero network.
    MaxNotFinite(f64),
    /// R's own error from `if (sum(P1_Pspatial) == 0)` when the sum is `NaN`.
    ///
    /// `NaN == 0` is `NA` in R, and `if (NA)` stops with "missing value where TRUE/FALSE needed".
    /// The port has no `if`, so without this the comparison is simply false and the kernel
    /// proceeds with an all-`NaN` network and *returns* it. That is worse than an error: the caller
    /// gets a silent `Prob = NaN` where upstream refuses, and no later stage can tell the
    /// difference between "no communication" and "the arithmetic went missing".
    IfNa,
}

impl std::fmt::Display for KernelError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            // Curly quotes: R's `match.arg` formats with `sQuote()` (U+201C/U+201D), and
            // the message is compared byte for byte. The trailing newline is R's too.
            KernelError::MatchArg(_) => write!(
                f,
                "'arg' should be one of \u{201c}triMean\u{201d}, \u{201c}truncatedMean\u{201d}, \
                 \u{201c}thresholdedMean\u{201d}, \u{201c}median\u{201d}"
            ),
            KernelError::UnusedLevels => write!(
                f,
                "Please check `unique(object@idents)` and ensure that the factor levels are correct!\n\
                 You may need to drop unused levels using 'droplevels' function. e.g.,\n\
                 `meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`"
            ),
            KernelError::Expr(e) => write!(f, "{e}"),
            KernelError::MaxNotFinite(m) => {
                write!(f, "max(data.use) is {m}; data/max(data) is all NaN")
            }
            KernelError::IfNa => write!(f, "missing value where TRUE/FALSE needed"),
        }
    }
}

impl std::error::Error for KernelError {}

/// R's `max()` over a numeric vector, which is **not** Rust's `f64::max`.
///
/// R's `max` propagates `NaN` -- and therefore `NA_real_`, which is a `NaN` with a payload -- so
/// `max(c(1, NA))` is `NA`. Rust's `f64::max` is a `maxNum`: it *skips* `NaN` and returns the other
/// operand, so `1.0_f64.max(f64::NAN)` is `1.0`. Using it here would make a matrix containing a
/// missing value normalise against a finite maximum instead of producing `NaN` throughout, and the
/// two diverge exactly where it matters: every `data.use` entry becomes `NaN`, every group mean
/// becomes `NaN`, and `if (sum(P1 * ...) == 0)` raises "missing value where TRUE/FALSE needed".
///
/// `-Inf` is a *value*, not a missing one: `max(c(-Inf, 1))` is `1` and `max(-Inf)` is `-Inf`. So
/// only `NaN` short-circuits.
pub fn r_max(v: &[f64]) -> f64 {
    let mut m = f64::NEG_INFINITY;
    for &x in v {
        if x.is_nan() {
            return f64::NAN;
        }
        if x > m {
            m = x;
        }
    }
    m
}

/// `data.use <- data/max(data)` -- upstream's very first numerical step.
///
/// It is easy to read this as a normalisation and skip it, and that is exactly what happened: the
/// port omitted it, and **every fixture in the suite had `max(data) == 1`**, so the omission was a
/// no-op in all of them and the differential gate stayed green. The property it buys is
/// *scale-invariance*: upstream's `Prob` depends only on `data/max(data)`, so scaling the input
/// matrix by any constant leaves `Prob` bit-for-bit unchanged. Without the division the port's
/// `Prob` moves, because the Hill terms `n^Kh * P^Kh / (Kh + P^Kh)` are not scale-invariant.
///
/// The degenerate maxima are not special-cased, and that is deliberate. Each one falls out of the
/// arithmetic and the port must reach the same place by the same route:
///
/// * `max` finite and positive -- the normal case;
/// * `max == 0` (an all-zero matrix) -- `0/0` is `NaN` everywhere, every group mean is `NaN`, and
///   upstream stops with "missing value where TRUE/FALSE needed";
/// * `max` is `NaN` or `NA` -- `x/NaN` is `NaN` everywhere, same stop;
/// * `max == Inf` -- finite entries become `0` and the `Inf` entry becomes `NaN`, so most pathways
///   have a zero total and upstream returns an all-zero `Prob` with no error at all.
pub fn scale_by_max(data: &[f64]) -> Vec<f64> {
    let m = r_max(data);
    data.iter().map(|v| v / m).collect()
}

impl From<ExprError> for KernelError {
    fn from(e: ExprError) -> Self {
        KernelError::Expr(e)
    }
}

/// One L-R pair, with every name already resolved to gene indices.
///
/// Built once per interaction. The bootstrap loop then costs no string hashing, no
/// `data.frame` subsetting and no `intersect` -- which is where upstream spends most of
/// the time it does not spend in `triMean`.
#[derive(Clone, Debug)]
pub struct LrPlan {
    pub ligand: EntityRows,
    pub receptor: EntityRows,
    /// Present subunits of `co_A_receptor`, `None` when the field is empty/unknown.
    pub co_a: Option<Vec<usize>>,
    pub co_i: Option<Vec<usize>>,
    /// Present subunits of `agonist` / `antagonist`, and whether the pair has one at all.
    /// Upstream gates on `which(!is.na(x) & x != "")`, so a *present but unresolvable*
    /// cofactor still enters the branch and yields all ones.
    pub agonist: Option<Vec<usize>>,
    pub antagonist: Option<Vec<usize>>,
    /// `as.numeric(table(group))/nC`, repeated over interactions. Empty unless
    /// `population.size`.
    pub pop: Vec<f64>,
    /// `rownames(pairLRsig)[i]`, for the `dimnames` of `Prob`/`Pval`.
    pub label: String,
}

/// `P2`/`P3` for one (pair, matrix) given already-resolved rows.
#[inline]
fn hill_product(
    rows: &[usize],
    m: &GroupMeans<'_>,
    j: usize,
    kh: f64,
    n: f64,
    agonist: bool,
) -> f64 {
    let khn = kh.powf(n);
    if rows.len() == 1 {
        let x = m.at(rows[0], j);
        let xn = x.powf(n);
        return if agonist {
            1.0 + xn / (khn + xn)
        } else {
            khn / (khn + xn)
        };
    }
    let mut buf = vec![0.0; rows.len()];
    for (t, &r) in rows.iter().enumerate() {
        let xn = m.at(r, j).powf(n);
        buf[t] = if agonist {
            1.0 + xn / (khn + xn)
        } else {
            khn / (khn + xn)
        };
    }
    r_prod(&buf, false)
}

/// `1 + x` for one subunit, or the `LDOUBLE` product of `1 + x` for several.
#[inline]
fn coreceptor_value(rows: &[usize], m: &GroupMeans<'_>, j: usize) -> f64 {
    if rows.len() == 1 {
        return 1.0 + m.at(rows[0], j);
    }
    let mut buf = vec![0.0; rows.len()];
    for (t, &r) in rows.iter().enumerate() {
        buf[t] = 1.0 + m.at(r, j);
    }
    r_prod(&buf, false)
}

/// The per-group ligand and receptor vectors for one (pair, matrix), including the
/// co-agonist/co-inhibitor scaling.
fn expr_pair(plan: &LrPlan, m: &GroupMeans<'_>, kh: f64, n: f64) -> (Vec<f64>, Vec<f64>) {
    let k = m.n_groups();
    let mut l = vec![0.0; k];
    for j in 0..k {
        l[j] = match &plan.ligand {
            EntityRows::Single(r) => m.at(*r, j),
            EntityRows::Complex(rows) => geom_mean_over(rows, m, j),
        };
    }
    let mut r = vec![0.0; k];
    for j in 0..k {
        let mut v = match &plan.receptor {
            EntityRows::Single(rr) => m.at(*rr, j),
            EntityRows::Complex(rows) => geom_mean_over(rows, m, j),
        };
        // `dataRavg * co.A / co.I`, in that association, as upstream writes it.
        if let Some(rows) = &plan.co_a {
            v *= coreceptor_value(rows, m, j);
        }
        if let Some(rows) = &plan.co_i {
            v /= coreceptor_value(rows, m, j);
        }
        r[j] = v;
    }
    // The agonist/antagonist factors are separate outer products upstream (P2, P3), so
    // they are returned alongside rather than folded in here.
    let _ = (kh, n);
    (l, r)
}

#[inline]
fn geom_mean_over(rows: &[usize], m: &GroupMeans<'_>, j: usize) -> f64 {
    let mut logs = Vec::with_capacity(rows.len());
    for &r in rows {
        logs.push(m.at(r, j).ln());
    }
    crate::stats::r_mean_no_rm(&logs).exp()
}

/// Everything `computeCommunProb` needs, already normalised and indexed.
pub struct Kernel {
    /// `data.use.avg` rownames: the gene universe, in `data.use` rownames order.
    pub genes: Vec<String>,
    /// `data.use.avg`, `n_genes x n_groups` column-major, for the observed data.
    pub avg: Vec<f64>,
    /// `data.use.avg.boot`, `nboot * n_genes * n_groups`, replicate-major.
    pub boot: Vec<f64>,
    pub plans: Vec<LrPlan>,
    /// `nboot`.
    pub nboot: usize,
    /// Hill parameters.
    pub kh: f64,
    pub n: f64,
    /// `nlevels(group)`.
    pub n_groups: usize,
    /// `levels(group)`.
    pub groups: Vec<String>,
    /// `P.spatial` as constructed by upstream, `n_groups^2` column-major.
    pub p_spatial: Vec<f64>,
    /// `adj.contact`, same layout. Unused for RNA.
    pub adj_contact: Vec<f64>,
    /// 0-based `nLR1`: interactions with index `> n_lr1` use `P.spatial * adj.contact`.
    pub n_lr1: usize,
    /// `p_spatial * adj.contact`, precomputed. The upstream mutation is idempotent only
    /// because `adj.contact` is 0/1, so this is exactly what the loop would converge to
    /// after the first mutation past `n_lr1`.
    p_scaled: Vec<f64>,
    /// `as.numeric(table(group[permutation[, nE]]))/nC` for each replicate,
    /// `nboot * n_groups`. Only populated when `population.size`.
    pop_boot: Vec<f64>,
}

impl Kernel {
    /// Resolve every interaction against one expression matrix, producing a [`LrPlan`] each.
    ///
    /// This is the only string-touching pass in the whole computation. Upstream performs
    /// the equivalent work `nLR * (nboot + 1)` times, inside data.frame subsetting.
    ///
    /// # Panics
    /// Never, for well-formed input: an unresolvable ligand or receptor is an
    /// `Err(KernelError::Expr)`, reproducing upstream's `subscript out of bounds`.
    pub fn resolve(
        genes: Vec<String>,
        lr: &[LrPair],
        db: &Database,
        group: &[usize],
        n_groups: usize,
        n_cells: usize,
        population_size: bool,
    ) -> Result<Vec<LrPlan>, KernelError> {
        // Name-only: the plans are resolved once and reused for all nboot + 1 matrices,
        // which share rownames but not values.
        let probe = GroupMeans::names_only(&genes, n_groups);
        let pop: Vec<f64> = {
            let mut counts = vec![0.0; n_groups];
            for &g in group {
                counts[g] += 1.0;
            }
            counts.iter().map(|c| c / n_cells as f64).collect()
        };
        let mut out = Vec::with_capacity(lr.len());
        for p in lr {
            let resolve = |name: &str| -> Result<Option<Vec<usize>>, KernelError> {
                if name.is_empty() {
                    return Ok(None);
                }
                if let Some(r) = probe.row(name) {
                    return Ok(Some(vec![r]));
                }
                // `intersect(subunits, rownames(data.use))` -- the cofactor path never
                // errors, it degrades to "no modulation".
                let cells = db.cofactor_cells(name);
                let rows: Vec<usize> = match cells {
                    Some(c) => c
                        .iter()
                        .filter(|s| !s.is_empty())
                        .filter_map(|s| probe.row(s))
                        .collect(),
                    None => Vec::new(),
                };
                Ok(if rows.is_empty() { None } else { Some(rows) })
            };
            out.push(LrPlan {
                ligand: crate::expr::resolve_entity(&p.ligand, &probe, db)?,
                receptor: crate::expr::resolve_entity(&p.receptor, &probe, db)?,
                co_a: resolve(&p.co_a)?,
                co_i: resolve(&p.co_i)?,
                agonist: resolve(&p.agonist)?,
                antagonist: resolve(&p.antagonist)?,
                pop: if population_size {
                    pop.clone()
                } else {
                    Vec::new()
                },
                label: p.label.clone(),
            });
        }
        Ok(out)
    }

    /// `set.seed(seed); permutation <- replicate(nboot, sample.int(nC, nC))`, then
    /// `data.use.avg.boot`.
    ///
    /// `data` is `t(data.use)`, i.e. **cells x genes** column-major, with `group[c]` the
    /// level index of cell `c`.
    ///
    /// Note the transposition. `object@data.signaling` is genes x cells upstream; the
    /// aggregate is `aggregate(t(data.use), list(group), FUN)`, so what arrives here has
    /// already been flipped. `aggregate_columns` then indexes `m[cell, gene]` at
    /// `gene * n_cells + cell`, i.e. the gene index strides.
    pub fn build_boot(
        data: &[f64],
        n_genes: usize,
        n_cells: usize,
        group: &[usize],
        n_groups: usize,
        nboot: usize,
        seed: i32,
        fun: GroupMean,
    ) -> (Vec<f64>, Vec<Vec<usize>>) {
        let cols = all_cols(n_genes);

        // Two phases, and the split is forced by the RNG rather than chosen for convenience.
        //
        // **The permutations must be drawn sequentially.** R's `computeCommunProb` does
        // `set.seed(seed.use); permutation <- replicate(nboot, sample.int(nC, nC))`, so replicate
        // `nE` consumes the `nE`-th slice of one MT19937 stream. Drawing them in parallel would
        // need one stream per replicate, which is a *different* permutation set -- not a
        // reordering, a different answer. So phase one is a sequential loop over the stream.
        //
        // **The aggregation is embarrassingly parallel and bit-safe.** Each replicate's
        // `aggregate_1` reads a disjoint permutation of the same input and writes a disjoint
        // `n_genes * n_groups` slot; nothing is shared and nothing is accumulated across
        // replicates. `collect()` on an indexed parallel iterator preserves input order, so
        // replicate `nE` still lands at offset `nE * n_genes * n_groups` and the flattened buffer
        // is byte-identical to the sequential one -- including the `F80` accumulation order
        // *within* a replicate, which is what the LONG_DOUBLE parity tests pin.
        //
        // Measured on the human-skin fixture at `nboot = 100` (see docs/BENCHMARKS.md): the
        // sequential version took 10.66 s at 1 thread and 10.63 s at 8 -- flat, because this loop
        // *was* the loop and nothing else in the kernel was parallel. The bootstrap aggregation is
        // 92-95% of runtime, so this is the whole opportunity.
        let mut rng = MersenneTwister::new(seed);
        let mut perms: Vec<Vec<usize>> = Vec::with_capacity(nboot);
        for _ in 0..nboot {
            // `sample.int` returns **1-based** cell indices, matching R. `permutation` is
            // then indexed as `group[permutation[, nE]]`, so the **0-based cell** indices are what
            // goes back to the caller: `pop_boot` does `counts[groups[c]] += 1` and needs cells,
            // not levels. Returning the group assignments here would silently count the wrong
            // thing whenever `population.size = TRUE`.
            perms.push(
                rng.sample_int_permutation(n_cells)
                    .into_iter()
                    .map(|c| c as usize - 1)
                    .collect(),
            );
        }
        // `group[permutation[, nE]]`: the group label each cell takes in this replicate, which
        // is what `aggregate_1` is given.
        let gboots: Vec<Vec<usize>> = perms
            .iter()
            .map(|p| p.iter().map(|&c| group[c]).collect())
            .collect();

        use rayon::prelude::*;
        let per_replicate: Vec<Vec<f64>> = gboots
            .par_iter()
            .map(|gboot| {
                // The per-replicate aggregation is `aggregate(t(data.use), list(groupboot), FUN)`
                // transposed -- the same `aggregate_columns` the observed pass uses, so the
                // accumulation order is identical by construction.
                // `aggregate_1` returns the groups x genes intermediate *before* R's `t()`, and
                // that is already R's flat layout for a genes x groups matrix (index
                // `group * n_genes + gene`, which is `GroupMeans::at(gene, group)` verbatim).
                // Going through `aggregate_columns` would transpose it into a different, equally
                // plausible order -- see the note on that function.
                let avg = aggregate_1(data, n_cells, &cols, gboot, n_groups, fun);
                debug_assert_eq!(avg.len(), n_genes * n_groups);
                avg
            })
            .collect();
        let mut out = Vec::with_capacity(nboot * n_genes * n_groups);
        for avg in &per_replicate {
            out.extend_from_slice(avg);
        }
        (out, perms)
    }

    /// Assemble a kernel, precomputing the hoisted `P.spatial * adj.contact`.
    #[allow(clippy::too_many_arguments)]
    pub fn new(
        genes: Vec<String>,
        avg: Vec<f64>,
        boot: Vec<f64>,
        plans: Vec<LrPlan>,
        nboot: usize,
        kh: f64,
        n: f64,
        n_groups: usize,
        groups: Vec<String>,
        p_spatial: Vec<f64>,
        adj_contact: Vec<f64>,
        n_lr1: usize,
        pop_boot: Vec<f64>,
    ) -> Result<Self, KernelError> {
        let kk = n_groups * n_groups;
        if p_spatial.len() != kk || adj_contact.len() != kk {
            return Err(KernelError::UnusedLevels);
        }
        if avg.len() != genes.len() * n_groups {
            return Err(KernelError::UnusedLevels);
        }
        let p_scaled: Vec<f64> = p_spatial
            .iter()
            .zip(&adj_contact)
            .map(|(a, b)| a * b)
            .collect();
        Ok(Self {
            genes,
            avg,
            boot,
            plans,
            nboot,
            kh,
            n,
            n_groups,
            groups,
            p_spatial,
            adj_contact,
            n_lr1,
            p_scaled,
            pop_boot,
        })
    }

    /// `P.spatial` for the observed pass: `P0`, or `P0 * adj.contact` past `n_lr1`.
    #[inline]
    fn p_eff(&self, i: usize) -> &[f64] {
        if i > self.n_lr1 {
            &self.p_scaled
        } else {
            &self.p_spatial
        }
    }

    /// `Prob` and `Pval`, both `n_groups x n_groups x nLR` column-major (R's `array(0,
    /// dim = c(k, k, nLR))` fills last-dimension-fastest, so slice `i` is a
    /// `n_groups x n_groups` column-major matrix).
    ///
    /// `Prob[Prob == 0] -> Pval == 1` is already applied.
    pub fn prob_arrays(&self) -> Result<(Vec<f64>, Vec<f64>), KernelError> {
        let k = self.n_groups;
        let kk = k * k;
        let n_lr = self.plans.len();
        let mut prob = vec![0.0f64; kk * n_lr];
        let mut pval = vec![0.0f64; kk * n_lr];
        let genes = self.genes.clone();
        let observed = GroupMeans::new(&self.avg, &genes, k);
        let khn = self.kh.powf(self.n);

        for (i, plan) in self.plans.iter().enumerate() {
            let psp = self.p_eff(i);
            let (lv, rv) = expr_pair(plan, &observed, self.kh, self.n);
            // P2 / P3 / P4 for the observed matrix.
            let (mut p2, mut p3) = (vec![1.0; k], vec![1.0; k]);
            if let Some(rows) = &plan.agonist {
                for j in 0..k {
                    p2[j] = hill_product(rows, &observed, j, self.kh, self.n, true);
                }
            }
            if let Some(rows) = &plan.antagonist {
                for j in 0..k {
                    p3[j] = hill_product(rows, &observed, j, self.kh, self.n, false);
                }
            }

            // P1 * P.spatial, and the early-exit test.
            let mut p1psp = vec![0.0f64; kk];
            for b in 0..k {
                for a in 0..k {
                    let lr = lv[a] * rv[b];
                    let lrn = lr.powf(self.n);
                    let p1 = lrn / (khn + lrn);
                    p1psp[b * k + a] = p1 * psp[b * k + a];
                }
            }
            // `sum()` in R accumulates in LDOUBLE over the column-major vector. With P1 and
            // P.spatial both non-negative the test is really "all zero", but the 80-bit
            // accumulation is what makes the comparison exact rather than approximately
            // exact, and it costs k^2 adds once per interaction.
            let mut acc = F80::from_f64(0.0);
            for v in &p1psp {
                acc = acc.add(F80::from_f64(*v));
            }
            // `NaN == 0` is `NA` in R and `if (NA)` stops. Checked *before* the zero test, because
            // in Rust the zero test is simply false for `NaN` and the loop would continue.
            let acc = acc.to_f64();
            if acc.is_nan() {
                return Err(KernelError::IfNa);
            }
            if acc == 0.0 {
                // R-ism: `Pval` is set to 1, not 1/nboot, and the pair is never resampled.
                prob[i * kk..(i + 1) * kk].copy_from_slice(&p1psp);
                pval[i * kk..(i + 1) * kk].fill(1.0);
                continue;
            }

            // Pnull = P1 * P2 * P3 * P4 * P.spatial
            //
            // Two details, both observable:
            // * `P2` and `P3` are each `Matrix::crossprod` of a 1 x k row vector, so
            //   `P2[a,b] = p2[a] * p2[b]` and `P3[a,b] = p3[a] * p3[b]`. Writing
            //   `p2[a] * p3[b]` gives a different matrix that happens to agree whenever
            //   p2 == p3 or one of them is 1 -- i.e. for every pair without both an
            //   agonist and an antagonist, which is most of them.
            // * R evaluates `P1*P2*P3*P4*P.spatial` strictly left to right, and each of
            //   the intermediate roundings is a rounding. So the multiplication order is
            //   reproduced literally rather than regrouped.
            let mut pnull = vec![0.0f64; kk];
            for b in 0..k {
                for a in 0..k {
                    let lr = lv[a] * rv[b];
                    let lrn = lr.powf(self.n);
                    let p1 = lrn / (khn + lrn);
                    let p4 = if plan.pop.is_empty() {
                        1.0
                    } else {
                        plan.pop[a] * plan.pop[b]
                    };
                    let mut t = p1 * (p2[a] * p2[b]);
                    t *= p3[a] * p3[b];
                    t *= p4;
                    t *= psp[b * k + a];
                    pnull[b * k + a] = t;
                }
            }
            prob[i * kk..(i + 1) * kk].copy_from_slice(&pnull);

            // Bootstrap: rowSums(Pboot - Pnull > 0) -- strict, so ties are not rejections
            // (R-ism 7).
            let mut n_reject = vec![0i64; kk];
            for e in 0..self.nboot {
                let off = e * self.avg.len();
                let m = GroupMeans::new(&self.boot[off..off + self.avg.len()], &genes, k);
                let (lb, rb) = expr_pair(plan, &m, self.kh, self.n);
                // P2 and P3 for this replicate: full per-group vectors, so the outer
                // products below can index both ends. (Indexing only one end and squaring
                // is right on the diagonal and wrong everywhere else.)
                let mut p2v = vec![1.0f64; k];
                if let Some(rows) = &plan.agonist {
                    for j in 0..k {
                        p2v[j] = hill_product(rows, &m, j, self.kh, self.n, true);
                    }
                }
                let mut p3v = vec![1.0f64; k];
                if let Some(rows) = &plan.antagonist {
                    for j in 0..k {
                        p3v[j] = hill_product(rows, &m, j, self.kh, self.n, false);
                    }
                }
                let p4v: &[f64] = if plan.pop.is_empty() {
                    &[]
                } else {
                    &self.pop_boot[e * k..(e + 1) * k]
                };
                for b in 0..k {
                    for a in 0..k {
                        let lr = lb[a] * rb[b];
                        let lrn = lr.powf(self.n);
                        let p1 = lrn / (khn + lrn);
                        let p4 = if p4v.is_empty() { 1.0 } else { p4v[a] * p4v[b] };
                        // Same left-to-right association as
                        // `P1.boot*P2.boot*P3.boot*P4.boot*P.spatial`.
                        let mut pb = p1 * (p2v[a] * p2v[b]);
                        pb *= p3v[a] * p3v[b];
                        pb *= p4;
                        pb *= psp[b * k + a];
                        if pb - pnull[b * k + a] > 0.0 {
                            n_reject[b * k + a] += 1;
                        }
                    }
                }
            }
            let nboot = self.nboot as f64;
            for t in 0..kk {
                pval[i * kk + t] = n_reject[t] as f64 / nboot;
            }
        }
        // R-ism 6: after the whole loop, `Pval[Prob == 0] <- 1`.
        for t in 0..prob.len() {
            if prob[t] == 0.0 {
                pval[t] = 1.0;
            }
        }
        Ok((prob, pval))
    }

    /// `dimnames(Prob)`: `list(levels(group), levels(group), rownames(pairLRsig))`.
    pub fn dimnames(&self) -> (Vec<String>, Vec<String>, Vec<String>) {
        (
            self.groups.clone(),
            self.groups.clone(),
            self.plans.iter().map(|p| p.label.clone()).collect(),
        )
    }
}

/// One row of `object@LR$LRsig`, with the fields the kernel reads.
#[derive(Clone, Debug, Default)]
pub struct LrPair {
    pub ligand: String,
    pub receptor: String,
    pub agonist: String,
    pub antagonist: String,
    pub co_a: String,
    pub co_i: String,
    /// `rownames(pairLRsig)[i]`, i.e. `"ligand^receptor"`.
    pub label: String,
}

/// Names -> row index, for callers that hold a gene universe as a map.
pub fn index_map(genes: &[String]) -> HashMap<&str, usize> {
    genes
        .iter()
        .enumerate()
        .map(|(i, g)| (g.as_str(), i))
        .collect()
}

/// Unused-import guard for [`compute_expr_lr`] / [`compute_expr_coreceptor`]: the kernel
/// deliberately does *not* call the public wrappers, because they re-resolve names. This
/// function exists so the two are kept in step by a test rather than by eye.
#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn the_public_expr_wrappers_agree_with_the_inlined_hot_path() {
        // `expr_pair` inlines the same arithmetic as `compute_expr_lr` +
        // `compute_expr_coreceptor`. If they ever drift, the hot loop is the one that
        // ships, so pin them together.
        let genes: Vec<String> = (0..8).map(|i| format!("G{i}")).collect();
        let vals: Vec<f64> = (0..8 * 3).map(|i| 0.1 + (i as f64) * 0.07).collect();
        let m = GroupMeans::new(&vals, &genes, 3);
        let mut db = Database::default();
        db.complexes
            .insert("AB".into(), vec!["G0".into(), "G1".into()]);
        db.n_subunit_cols = 2;
        db.cofactors.insert("cA".into(), vec!["G2".into()]);
        db.n_cofactor_cols = 1;

        let names = vec!["G0".to_string(), "AB".to_string()];
        let lr = compute_expr_lr(&names, &m, &db).unwrap();
        for (i, ent) in [EntityRows::Single(0), EntityRows::Complex(vec![0, 1])]
            .iter()
            .enumerate()
        {
            for j in 0..3 {
                let got = match ent {
                    EntityRows::Single(r) => m.at(*r, j),
                    EntityRows::Complex(rows) => geom_mean_over(rows, &m, j),
                };
                assert_eq!(
                    got.to_bits(),
                    lr[j * 2 + i].to_bits(),
                    "entity {i} group {j}"
                );
            }
        }
        let co = compute_expr_coreceptor(&["cA".to_string()], &m, &db, CoreceptorKind::Activator);
        for j in 0..3 {
            // `j * 1 + 0` spells the layout out; `co[j]` would be the same index and would stop
            // testing that it is. Bound to a local so the `* 1`/`+ 0` survive as a value the
            // compiler cannot fold, and the stride stays visible to a reader.
            let cell_j = j;
            assert_eq!(
                coreceptor_value(&[2], &m, j).to_bits(),
                co[cell_j].to_bits()
            );
        }
    }

    #[test]
    fn match_arg_rejects_unknown_types() {
        assert_eq!(match_mean_type("triMean"), Ok(GroupMean::TriMean));
        assert_eq!(match_mean_type("median"), Ok(GroupMean::Median));
        assert!(matches!(
            match_mean_type("nope"),
            Err(KernelError::MatchArg(_))
        ));
    }
}

#[cfg(test)]
mod observed_means_tests {
    use super::*;
    use crate::aggregate::GroupMean;

    /// The port returned an all-zero `Prob` for six matrix configurations because the raise
    /// upstream performs inside `thresholdedMean` happened on a gene no interaction referenced.
    #[test]
    fn thresholded_mean_nan_in_an_unused_gene_is_still_a_raise() {
        let fun = GroupMean::ThresholdedMean { trim: 0.1 };
        // One missing group mean, as `aggregate` over all genes leaves behind.
        let avg = vec![0.0, 1.0, f64::NAN, 2.0];
        assert_eq!(check_observed_means(fun, &avg), Err(KernelError::IfNa));
    }

    /// For the other means a `NaN` group mean is a value, not an error: upstream only turns it
    /// into `missing value where TRUE/FALSE needed` later, at `if (sum(P1_Pspatial) == 0)`.
    #[test]
    fn nan_is_only_a_raise_for_thresholded_mean() {
        let avg = vec![0.0, 1.0, f64::NAN, 2.0];
        for fun in [
            GroupMean::Mean,
            GroupMean::TriMean,
            GroupMean::TrimmedMean { trim: 0.1 },
            GroupMean::Median,
        ] {
            assert_eq!(check_observed_means(fun, &avg), Ok(()), "{fun:?}");
        }
    }

    #[test]
    fn finite_observed_means_never_raise() {
        assert_eq!(
            check_observed_means(GroupMean::ThresholdedMean { trim: 0.1 }, &[0.0, 1.5, 2.5]),
            Ok(())
        );
    }
}
