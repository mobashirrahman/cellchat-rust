//! The kernel pipeline, with no R in it.
//!
//! This is the same sequence the extendr binding runs, and the order is load-bearing rather than
//! incidental:
//!
//! 1. `data.use <- data / max(data)` — once, before *both* consumers.
//! 2. the observed per-group aggregate, then `check_observed_means`, **before** the bootstrap is
//!    built, because upstream aggregates the observed data at `modelling.R:115` and only builds the
//!    replicates at `:208`. Reversing the two changes which error is reported first.
//! 3. the bootstrap replicates, whose permutations must be drawn **sequentially** from one MT19937
//!    stream, then aggregated in parallel. `sample.int(nC, nC)` in replicate `nE` consumes the
//!    `nE`-th slice of one stream, so drawing them in parallel would need one stream per replicate
//!    and produce a *different* permutation set — not a reordering, a different answer.
//! 4. `population.size` recomputes `table(group[permutation[, nE]])/nC` per replicate.
//!
//! The step numbers refer to upstream's `R/modeling.R`.

use crate::input::Input;
use r_core::aggregate::GroupMean;
use r_core::aggregate::{aggregate_1, all_cols};
use r_core::prob::Kernel;
use r_core::stats::{geometric_mean, r_trimmed_mean, thresholded_mean, tri_mean};

/// `FunMean` resolved from the requested name and the caller's `trim`.
///
/// `r_core::prob::match_mean_type` cannot be used here: it hard-codes `trim = 0.1` when it builds
/// the `GroupMean`, because in upstream's flow `match.arg` runs on `type` alone and `trim` is
/// threaded in separately afterwards. So the partial-matching logic is reproduced here against the
/// same constant, and the `trim` that the caller actually asked for is what ends up in the enum.
/// Getting this wrong would be silent -- `trim` has no effect on `triMean` or `median`, so only two
/// of the four `type.mean` values would notice.
const MEAN_TYPES: [&str; 4] = ["triMean", "truncatedMean", "thresholdedMean", "median"];

fn fun_mean(name: &str, trim: f64) -> Result<GroupMean, String> {
    // R's `match.arg` partial matching: exact, else a unique prefix, else upstream's error. The
    // resolved name is what the caller records, so `describe` and the parity gate both report the
    // matched spelling rather than the typed one.
    let resolved = if let Some(&t) = MEAN_TYPES.iter().find(|&&t| t == name) {
        t
    } else {
        let hits: Vec<&&str> = MEAN_TYPES.iter().filter(|t| t.starts_with(name)).collect();
        match hits.len() {
            1 => hits[0],
            _ => {
                return Err(format!(
                    "'arg' should be one of {}",
                    MEAN_TYPES
                        .iter()
                        .map(|s| format!("\"{s}\""))
                        .collect::<Vec<_>>()
                        .join(", ")
                ))
            }
        }
    };
    Ok(match resolved {
        "triMean" => GroupMean::TriMean,
        "truncatedMean" => GroupMean::TrimmedMean { trim },
        "thresholdedMean" => GroupMean::ThresholdedMean { trim },
        _ => GroupMean::Median,
    })
}

pub struct Run {
    pub prob: Vec<f64>,
    pub pval: Vec<f64>,
    pub n_groups: usize,
    pub n_lr: usize,
    /// The `groupMean` name as `options$parameter$type.mean` records it, i.e. after `match.arg`.
    pub type_mean_resolved: String,
    /// Per-gene observed means, so the caller can report or compare them.
    pub avg: Vec<f64>,
}

/// Run `computeCommunProb` for one dataset.
pub fn compute_commun_prob(inp: &Input, db: &r_core::db::Database) -> Result<Run, String> {
    let cfg = &inp.cfg;
    let fun = fun_mean(&cfg.type_mean, cfg.trim)?;
    let n_genes = inp.genes.len();
    let n_cells = inp.cells.len();
    let n_groups = inp.groups.len();
    let n_lr = inp.lr.len();

    // `raw.use` chooses *which* matrix, and nothing else. Upstream is two lines:
    //
    // ```r
    // if (raw.use) data <- as.matrix(object@data.signaling)
    // else          data <- as.matrix(object@data.smooth)
    // ...
    // data.use <- data/max(data)          # modelling.R:111, unconditional
    // ```
    //
    // There is no library-size normalisation anywhere in `computeCommunProb`. A first version of this
    // file divided by the row mean when `raw.use` was `FALSE`, on the reasonable-sounding guess that
    // "raw versus projected" implied a normalisation; the gate caught it -- upstream produced a full
    // `Prob` array and the CLI raised "missing value where TRUE/FALSE needed", because
    // `data.smooth` is a row-normalised matrix whose row mean is 1, so `data / 1` fed the Hill
    // function values several orders of magnitude above `Kh` and `P1` saturated to 1 everywhere. The
    // R shim never had this: it copies upstream's two lines verbatim.
    let base = if cfg.raw_use {
        inp.data.clone()
    } else {
        inp.data_smooth.clone().ok_or_else(|| {
            "raw_use is FALSE but the input carries no data_smooth matrix".to_string()
        })?
    };

    // 1. `data.use <- data / max(data)`.
    let data_scaled = r_core::prob::scale_by_max(&base);

    // 2. The observed aggregate, then the `NaN` check, before any resampling.
    let avg = aggregate_1(
        &data_scaled,
        n_cells,
        &all_cols(n_genes),
        &inp.cell_group,
        n_groups,
        fun,
    );
    r_core::prob::check_observed_means(fun, &avg).map_err(|e| e.to_string())?;

    // 3. Resolve the interactions, then build the replicates.
    let plans = Kernel::resolve(
        inp.genes.clone(),
        &inp.lr,
        db,
        &inp.cell_group,
        n_groups,
        n_cells,
        cfg.population_size,
    )
    .map_err(|e| e.to_string())?;

    let (boot, perms) = Kernel::build_boot(
        &data_scaled,
        n_genes,
        n_cells,
        &inp.cell_group,
        n_groups,
        cfg.nboot,
        cfg.seed,
        fun,
    );

    // 4. `population.size` per replicate.
    let pop_boot: Vec<f64> = if cfg.population_size {
        let mut v = Vec::with_capacity(cfg.nboot * n_groups);
        for p in &perms {
            let mut counts = vec![0.0f64; n_groups];
            for &c in p {
                counts[inp.cell_group[c]] += 1.0;
            }
            v.extend(counts.iter().map(|c| c / n_cells as f64));
        }
        v
    } else {
        Vec::new()
    };

    // Upstream's RNA branch: `P.spatial <- matrix(1, k, k)`, `adj.contact <- matrix(1, k, k)`,
    // `nLR1 <- nLR`. The spatial terms are not part of the standalone path, which is
    // `computeCellChat`-shaped rather than `computeRegionDistance`-shaped.
    let kk = n_groups * n_groups;
    let kernel = Kernel::new(
        inp.genes.clone(),
        avg.clone(),
        boot,
        plans,
        cfg.nboot,
        cfg.kh,
        cfg.n,
        n_groups,
        inp.groups.clone(),
        vec![1.0; kk],
        vec![1.0; kk],
        n_lr,
        pop_boot,
    )
    .map_err(|e| e.to_string())?;
    let (prob, pval) = kernel.prob_arrays().map_err(|e| e.to_string())?;

    Ok(Run {
        prob,
        pval,
        n_groups,
        n_lr,
        type_mean_resolved: match fun {
            GroupMean::TriMean => "triMean".into(),
            GroupMean::TrimmedMean { .. } => "truncatedMean".into(),
            GroupMean::ThresholdedMean { .. } => "thresholdedMean".into(),
            GroupMean::Median => "median".into(),
            GroupMean::Mean => "triMean".into(),
        },
        avg,
    })
}

/// `aggregateNet`'s unfiltered branch: the `K x K` count and weight matrices.
///
/// Ported because the CLI needs them to produce the same interaction table a user's next call would,
/// and because they are exact integer sums and additions rather than a solver's output.
pub fn aggregate_net(
    prob: &[f64],
    pval: &[f64],
    n_groups: usize,
    n_lr: usize,
    thresh: f64,
) -> (Vec<Vec<f64>>, Vec<Vec<f64>>) {
    let mut count = vec![vec![0.0f64; n_groups]; n_groups];
    let mut weight = vec![vec![0.0f64; n_groups]; n_groups];
    for i in 0..n_lr {
        for b in 0..n_groups {
            for a in 0..n_groups {
                let p = prob[a + b * n_groups + i * n_groups * n_groups];
                let pv = pval[a + b * n_groups + i * n_groups * n_groups];
                // `pval[prob == 0] <- 1; prob[pval >= thresh] <- 0` then `apply(prob > 0, c(1,2),
                // sum)`. The NaN-safe comparison matters: upstream's is `prob[pval >= thresh]`, and
                // a `pval` that is `NA` makes the subscript `NA` assignments do nothing, leaving
                // `prob` counted. The port reproduces that by treating an unorderable p-value as
                // surviving, which is what the `!(pv >= thresh)` form does.
                let zeroed = if pv >= thresh { 0.0 } else { p };
                if zeroed > 0.0 {
                    count[a][b] += 1.0;
                    weight[a][b] += zeroed;
                }
            }
        }
    }
    (count, weight)
}

/// `computeAveExpr` for one group, exposed so the CLI can report a per-group summary.
///
/// Present because it is part of the required numeric surface and a CLI that cannot show the
/// averaged expression is a CLI that cannot be checked against `computeAveExpr`'s output.
pub fn compute_ave_expr(inp: &Input) -> Result<Vec<f64>, String> {
    let fun = fun_mean(&inp.cfg.type_mean, inp.cfg.trim)?;
    let n_genes = inp.genes.len();
    // **Always** `data.signaling`, even when `cfg.raw_use` is `FALSE`.
    //
    // `computeAveExpr` has no `raw.use` parameter -- its signature is
    // `computeAveExpr(object, features, group.by, type, trim, slot.name = c("data.signaling", "data"),
    // data.use = NULL)` -- and it defaults to `slot.name = "data.signaling"`. `raw.use` belongs to
    // `computeCommunProb` alone. So the two functions read *different matrices* for the same object
    // when `raw.use = FALSE`, and that is upstream, not an oversight to be tidied up.
    //
    // An intermediate version of this file "fixed" the asymmetry by honouring `raw_use` here, and the
    // gate rejected it: `computeAveExpr means MISMATCH 96 values` on the `raw.use = FALSE`
    // configuration. The gate was right and the change was wrong. The CLI's `compute_commun_prob` does
    // honour `raw_use`; this does not, and the two disagreeing is the correct behaviour.
    let scaled = r_core::prob::scale_by_max(&inp.data);
    Ok(aggregate_1(
        &scaled,
        inp.cells.len(),
        &all_cols(n_genes),
        &inp.cell_group,
        inp.groups.len(),
        fun,
    ))
}

/// The four mean functions, named, for `cellchatrs mean`. Pure `r-core` re-exports so the CLI has
/// no reason to reimplement any of them.
pub fn means(x: &[f64], kind: &str, trim: f64) -> Result<Vec<f64>, String> {
    Ok(match kind {
        "triMean" => vec![tri_mean(x)],
        "truncatedMean" => vec![r_trimmed_mean(x, trim, true)],
        "thresholdedMean" => vec![thresholded_mean(x, trim)],
        "geometricMean" => vec![geometric_mean(x)],
        other => return Err(format!("unknown mean function {other:?}")),
    })
}
