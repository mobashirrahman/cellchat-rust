// extendr bindings for `r_core::prob` and `r_core::net` -- argument marshalling only.
//
// No numerics here. Every value that crosses the boundary is a plain `f64` or a `String`;
// the computation itself is `r_core::prob::Kernel`, so the binding layer cannot drift from
// the tested core.
//
// ## Why everything is flattened
//
// The R side hands over a dense genes x cells matrix as a column-major `NumericVector` plus
// an explicit `dim`, and the L-R table as six parallel string vectors. It does **not** hand
// over R matrix / data.frame objects: indexing those from Rust means calling back into R for
// every element, and the entire point of this port is that the hot path never crosses the
// FFI boundary.
//
// ## Why this file is `include!`d into `lib.rs` rather than a module
//
// `#[extendr]` emits a **private** `meta__<name>` beside each annotated function, and
// `extendr_module!` expands to `crate::meta__<name>` paths. Those paths only resolve if the
// annotated functions live in the crate root, so this file is textually included there.
// That is a real extendr 0.9 constraint, not a style choice.

use r_core::aggregate::{aggregate_1, all_cols, GroupMean};
use r_core::db::Database;
use r_core::prob::{Kernel, LrPair};

/// Unwrap a fallible step, raising a **real R error** instead of unwinding.
///
/// # Why not just return `extendr_api::Result`
///
/// With no `result_*` feature enabled, extendr's `From<Result<T, E>> for Robj` is literally
/// `Err(err) => panic!("{}", err)` (`robj/into_robj.rs:73`). So every `Err` a binding returns
/// reaches R as a Rust panic that extendr catches at the FFI boundary and converts back into an
/// R error. R cannot tell the difference -- the message and the condition class come out right,
/// which is why the differential tests pass -- but two things leak. extendr's own docs say
/// `panic!` there is "discouraged due to memory leakage", and each one prints
/// `thread '<unnamed>' panicked at .../into_robj.rs:73` to stderr, so an ordinary user-facing
/// validation failure looks like a crash.
///
/// `throw_r_error` raises the error through R's own condition system and diverges (`-> !`), so
/// it is a drop-in for both `?` and an early `return Err(..)`.
///
/// # Why not enable `result_condition`
///
/// Because it breaks bit-identity, which outranks tidiness. That impl builds
/// `list(message = "extendr_err", value = x)` with class `c("extendr_error", "error", "condition")`:
/// `conditionMessage()` then returns the literal string `"extendr_err"` and the real message
/// moves into `$value`. `tests/test-contract-inputs.R` asserts that this port's error text
/// matches upstream's exactly, and a `computeCommunProb` validation failure is one of the cases
/// it compares. Trading a stderr line for a wrong error message is a bad trade.
///
/// # Status
///
/// Applied to `compute_commun_prob` as the pilot, since it is the function the error-path tests
/// actually exercise. The other 30 bindings returning `extendr_api::Result` still take the
/// `panic!` route; converting them is the same mechanical change.
fn unwrap_r<T>(r: extendr_api::Result<T>) -> T {
    match r {
        Ok(v) => v,
        Err(e) => extendr_api::throw_r_error(e.to_string()),
    }
}

/// Load a `CellChatDB` export from disk and report what was found.
///
/// Returns a named list rather than the [`Database`] struct: `Database` has no `Robj`
/// representation and does not need one, since the kernel takes the directory path and
/// loads it itself. The shim calls this once to surface a corrupt or stale export as an R
/// error at the point the user can act on it, rather than deep inside the kernel.
///
/// # Arguments
/// * `dir`: a directory written by `tests/parity/export_db.R`, manifest-pinned by MD5.
/// * `species`: `"human"` or `"mouse"`, used only in error messages.
#[extendr]
pub fn db_load(dir: &str, species: &str) -> List {
    let db = unwrap_r(
        Database::load(std::path::Path::new(dir))
            .map_err(|e| extendr_api::Error::Other(format!("{species} CellChatDB at {dir:?}: {e}"))),
    );
    list! {
        species=db.species.clone(),
        n_interactions=db.interactions.len() as i32,
        n_complexes=db.complexes.len() as i32,
        n_cofactors=db.cofactors.len() as i32,
        n_symbols=db.symbols.len() as i32,
        n_lr1=db.nlr1 as i32,
    }
}

/// `data.use.avg`: `aggregate(t(data.use), list(group), FUN)`, in R's genes x groups order.
///
/// `data` is `t(data.use)`, i.e. **cells x genes** column-major, because that is what the
/// upstream aggregate consumes. `dim` is `c(n_genes, n_cells)` so the R shim can pass
/// `as.matrix(object@data.signaling)` through unchanged.
#[extendr]
pub fn average_expression(
    data: Vec<f64>,
    dim: Vec<i32>,
    group: Vec<i32>,
    group_levels: Vec<String>,
    mean_type: &str,
    trim: f64,
) -> Vec<f64> {
    let (n_genes, n_cells) = unwrap_r(dims(&dim));
    let groups = unwrap_r(group_index(&group, &group_levels));
    let fun = unwrap_r(mean_fun(mean_type, trim));
    if data.len() != n_genes * n_cells {
        extendr_api::throw_r_error(format!(
            "data has {} values but dim implies {}",
            data.len(),
            n_genes * n_cells
        ));
    }
    aggregate_1(
        &data,
        n_cells,
        &all_cols(n_genes),
        &groups,
        group_levels.len(),
        fun,
    )
}

/// `computeCommunProb` for RNA data. Returns `(Prob, Pval)` and the array dimensions.
///
/// # Arguments
/// * `data`: `t(data.use)`, cells x genes column-major.
/// * `dim`: `c(n_genes, n_cells)`.
/// * `genes`: `rownames(object@data.signaling)`, the gene universe and its order.
/// * `group`: 0-based level index per cell, indexing `group_levels`.
/// * `group_levels`: `levels(object@idents)`; also the `dimnames` of `Prob`.
/// * `lr_*`: one string per row of `object@LR$LRsig`, all the same length and order.
/// * `complex_*`, `cofactor_*`: the object's own subunit tables, column-major.
/// * `nboot`, `seed`, `kh`, `n`, `mean_type`, `trim`, `population_size`: the upstream
///   parameters, with upstream's names.
/// * `p_spatial` / `adj_contact` / `n_lr1`: spatial-only. Pass an empty `p_spatial` (the
///   default for RNA) to take upstream's RNA branch, where `P.spatial` and `adj.contact`
///   are both all-ones and `nLR1 = nLR`.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn compute_commun_prob(
    data: Vec<f64>,
    dim: Vec<i32>,
    genes: Vec<String>,
    group: Vec<i32>,
    group_levels: Vec<String>,
    lr_ligand: Vec<String>,
    lr_receptor: Vec<String>,
    lr_agonist: Vec<String>,
    lr_antagonist: Vec<String>,
    lr_co_a: Vec<String>,
    lr_co_i: Vec<String>,
    lr_label: Vec<String>,
    complex_names: Vec<String>,
    complex_subunits: Vec<String>,
    complex_n_cols: usize,
    cofactor_names: Vec<String>,
    cofactor_subunits: Vec<String>,
    cofactor_n_cols: usize,
    nboot: usize,
    seed: i32,
    kh: f64,
    n: f64,
    mean_type: &str,
    trim: f64,
    population_size: bool,
    p_spatial: Vec<f64>,
    adj_contact: Vec<f64>,
    n_lr1: i32,
    advance_rng: Function,
) -> List {
    let (n_genes, n_cells) = unwrap_r(dims(&dim));
    if data.len() != n_genes * n_cells {
        extendr_api::throw_r_error(format!(
            "data has {} values but dim implies {}",
            data.len(),
            n_genes * n_cells
        ));
    }
    if genes.len() != n_genes {
        extendr_api::throw_r_error(format!(
            "got {} gene names but dim says {n_genes} genes",
            genes.len()
        ));
    }
    let n_lr = lr_ligand.len();
    for (name, v) in [
        ("lr_receptor", &lr_receptor),
        ("lr_agonist", &lr_agonist),
        ("lr_antagonist", &lr_antagonist),
        ("lr_co_a", &lr_co_a),
        ("lr_co_i", &lr_co_i),
        ("lr_label", &lr_label),
    ] {
        if v.len() != n_lr {
            extendr_api::throw_r_error(format!(
                "{name} has {} entries but lr_ligand has {n_lr}"
            , v.len()));
        }
    }
    let n_groups = group_levels.len();
    if n_groups == 0 {
        extendr_api::throw_r_error("no cell groups");
    }
    let groups = unwrap_r(group_index(&group, &group_levels));
    if groups.len() != n_cells {
        extendr_api::throw_r_error(format!(
            "group has {} entries but data has {n_cells} cells",
            groups.len()
        ));
    }
    // Upstream errors when `nlevels(group) != length(unique(group))`; an unused level
    // therefore must not reach here silently.
    for g in 0..n_groups {
        if !groups.contains(&g) {
            extendr_api::throw_r_error(format!(
                "unused factor level {:?}: Please check `unique(object@idents)` and ensure \\
                 that the factor levels are correct!",
                group_levels[g]
            ));
        }
    }
    let fun = unwrap_r(mean_fun(mean_type, trim));

    let complexes = unwrap_r(subunit_table(&complex_names, &complex_subunits, complex_n_cols));
    let cofactors = unwrap_r(subunit_table(&cofactor_names, &cofactor_subunits, cofactor_n_cols));
    let db = Database::from_parts(
        "object@DB", Vec::new(), complexes, cofactors,
        complex_names, cofactor_names, complex_n_cols, cofactor_n_cols,
        std::collections::HashMap::new(),
    );

    let lr: Vec<LrPair> = (0..n_lr)
        .map(|i| LrPair {
            ligand: lr_ligand[i].clone(),
            receptor: lr_receptor[i].clone(),
            agonist: lr_agonist[i].clone(),
            antagonist: lr_antagonist[i].clone(),
            co_a: lr_co_a[i].clone(),
            co_i: lr_co_i[i].clone(),
            label: lr_label[i].clone(),
        })
        .collect();
    // `data.use <- data/max(data)`, upstream's first numerical step, done once before *both*
    // consumers -- the observed aggregate and the bootstrap replicates.
    //
    // It was missing from the port for the whole life of the kernel and nothing noticed, because
    // **every fixture in the suite had `max(data) == 1`**: `runif(0.01, 1)` reaches its upper bound
    // often enough that both the hand fixture and the configuration matrix normalised to the
    // identity. What the step actually buys is scale-*invariance* -- upstream's `Prob` depends only
    // on `data/max(data)`, so a constant factor on the matrix leaves it bit-for-bit unchanged --
    // and without the division the port's `Prob` moved by more than an order of magnitude when the
    // same data was divided by two. The differential gate could not see it because it never handed
    // the port a matrix whose maximum was not 1, which is why the configuration matrix now carries
    // an explicit `scale` axis.
    let data_scaled = r_core::prob::scale_by_max(&data);

    // Upstream aggregates the observed data (modelling.R:115) *before* it builds the bootstrap
    // replicates (:208), so a raise inside `FunMean` is reported before any resampling happens.
    // The port built the replicates first, which gave it a different error precedence.
    let avg = aggregate_1(&data_scaled, n_cells, &all_cols(n_genes), &groups, n_groups, fun);
    // Upstream applies `FunMean` to every gene; see `check_observed_means` for why a `NaN` here is
    // a raise rather than a value, and why the port used to swallow it.
    unwrap_r(r_core::prob::check_observed_means(fun, &avg).map_err(|e| kernel_err(&e)));

    // The only string-touching pass. Everything after this is `usize` indices.
    let plans = unwrap_r(
        Kernel::resolve(
            genes.clone(),
            &lr,
            &db,
            &groups,
            n_groups,
            n_cells,
            population_size,
        )
        .map_err(|e| kernel_err(&e)),
    );

    // Invoke R only on the main thread, after observed-data validation and before bootstrap.
    // Early errors leave the caller's RNG untouched; errors in the LR loop consume the draws.
    unwrap_r(advance_rng.call(pairlist!()));

    let (boot, perms) =
        Kernel::build_boot(&data_scaled, n_genes, n_cells, &groups, n_groups, nboot, seed, fun);
    // `population.size` recomputes `table(group[permutation[, nE]])/nC` per replicate.
    let pop_boot: Vec<f64> = if population_size {
        let mut v = Vec::with_capacity(nboot * n_groups);
        for p in &perms {
            let mut counts = vec![0.0f64; n_groups];
            for &c in p {
                counts[groups[c]] += 1.0;
            }
            v.extend(counts.iter().map(|c| c / n_cells as f64));
        }
        v
    } else {
        Vec::new()
    };
    let kk = n_groups * n_groups;
    // Upstream's RNA branch: `P.spatial <- matrix(1, numCluster, numCluster)`,
    // `adj.contact <- matrix(1, ...)`, `nLR1 <- nLR`.
    let (psp, adj, n_lr1) = if p_spatial.is_empty() {
        (vec![1.0; kk], vec![1.0; kk], n_lr)
    } else {
        if p_spatial.len() != kk {
            extendr_api::throw_r_error(format!(
                "p_spatial must have {kk} values, got {}",
                p_spatial.len()
            ));
        }
        if adj_contact.len() != kk {
            extendr_api::throw_r_error(format!(
                "adj_contact must have {kk} values, got {}",
                adj_contact.len()
            ));
        }
        if n_lr1 < 0 || n_lr1 as usize > n_lr {
            extendr_api::throw_r_error(format!(
                "n_lr1 must be in [0, {n_lr}], got {n_lr1}"
            ));
        }
        (p_spatial, adj_contact, n_lr1 as usize)
    };

    let kernel = unwrap_r(
        Kernel::new(
            genes,
            avg,
            boot,
            plans,
            nboot,
            kh,
            n,
            n_groups,
            group_levels,
            psp,
            adj,
            n_lr1,
            pop_boot,
        )
        .map_err(|e| kernel_err(&e)),
    );
    // `prob_arrays` is fallible now: R's `if (sum(P1_Pspatial) == 0)` raises
    // "missing value where TRUE/FALSE needed" when the sum is `NaN`, and a kernel that returns an
    // all-`NaN` network instead is strictly worse than one that stops.
    let (prob, pval) = unwrap_r(kernel.prob_arrays().map_err(|e| kernel_err(&e)));
    // `list(prob=, pval=, dim=)`. The two vectors are column-major flat buffers, which is
    // exactly R's array storage order, so the shim only has to `dim<-` them.
    list! {
        prob=r!(prob),
        pval=r!(pval),
        dim=r!(vec![n_groups as i32, n_groups as i32, n_lr as i32]),
    }
}

/// Reassemble a column-major R subunit table without guessing a species or export.
fn subunit_table(
    names: &[String],
    cells: &[String],
    n_cols: usize,
) -> extendr_api::Result<std::collections::HashMap<String, Vec<String>>> {
    let n_rows = names.len();
    if cells.len() != n_rows * n_cols {
        return Err(extendr_api::Error::Other("invalid subunit table dimensions".into()));
    }
    let mut rows = std::collections::HashMap::new();
    for (r, name) in names.iter().enumerate() {
        // R matches the first row when row names are repeated.
        rows.entry(name.clone()).or_insert_with(|| {
            (0..n_cols).map(|c| cells[c * n_rows + r].clone()).collect()
        });
    }
    Ok(rows)
}

/// Resolve a cofactor's present subunits against a gene universe.
///
/// Exposed so the R shim can build the equivalent of a [`r_core::prob::LrPlan`] for
/// diagnostics, and so `computeExpr_*` has a binding of its own. Returns the row indices
/// in `cofactor*` order; an unresolvable cofactor yields an empty vector, which upstream
/// treats as "no modulation" rather than an error.
#[extendr]
pub fn resolve_cofactor_rows(name: &str, genes: Vec<String>, db_dir: &str) -> Vec<i64> {
    if name.is_empty() {
        return Vec::new();
    }
    let db = unwrap_r(
        Database::load(std::path::Path::new(db_dir))
            .map_err(|e| extendr_api::Error::Other(format!("CellChatDB at {db_dir:?}: {e}"))),
    );
    let index: std::collections::HashMap<&str, usize> =
        genes.iter().enumerate().map(|(i, g)| (g.as_str(), i)).collect();
    let Some(cells) = db.cofactor_cells(name) else {
        return Vec::new();
    };
    cells
        .iter()
        .filter(|s| !s.is_empty())
        .filter_map(|s| index.get(s.as_str()).map(|&i| i as i64))
        .collect()
}

/// `match.arg(type)` over `type.mean`, returning the canonical spelling R would store in
/// `options$parameter$type.mean`.
#[extendr]
pub fn match_mean_type(type_: &str, trim: f64) -> String {
    match unwrap_r(mean_fun(type_, trim)) {
        // `GroupMean::Mean` is not reachable from `match.arg`: upstream's
        // `switch(aggregate.fun, mean = ...)` has no `mean` choice, and
        // `computeCommunProb`'s `type` argument does not offer it either.
        GroupMean::TriMean => "triMean",
        GroupMean::TrimmedMean { .. } => "truncatedMean",
        GroupMean::ThresholdedMean { .. } => "thresholdedMean",
        GroupMean::Median => "median",
        GroupMean::Mean => "triMean",
    }
    .to_string()
}

fn mean_fun(type_: &str, trim: f64) -> extendr_api::Result<GroupMean> {
    const TYPES: [&str; 4] = ["triMean", "truncatedMean", "thresholdedMean", "median"];
    if let Some(&t) = TYPES.iter().find(|&&t| t == type_) {
        return Ok(match t {
            "triMean" => GroupMean::TriMean,
            "truncatedMean" => GroupMean::TrimmedMean { trim },
            "thresholdedMean" => GroupMean::ThresholdedMean { trim },
            _ => GroupMean::Median,
        });
    }
    // `match.arg`'s partial matching: a unique prefix is accepted, and the *matched* name
    // is what gets recorded.
    let hits: Vec<&&str> = TYPES.iter().filter(|t| t.starts_with(type_)).collect();
    match hits.len() {
        1 => mean_fun(hits[0], trim),
        _ => Err(extendr_api::Error::Other(format!(
            "'arg' should be one of \"triMean\", \"truncatedMean\", \"thresholdedMean\", \\
             \"median\" (got \"{type_}\")"
        ))),
    }
}

fn dims(dim: &[i32]) -> extendr_api::Result<(usize, usize)> {
    if dim.len() != 2 {
        return Err(extendr_api::Error::Other("dim must have length 2".into()));
    }
    let (a, b) = (i64::from(dim[0]), i64::from(dim[1]));
    if a < 0 || b < 0 {
        return Err(extendr_api::Error::Other("dim must be non-negative".into()));
    }
    Ok((a as usize, b as usize))
}

// ----------------------------------------------------------------------------------- de.rs

/// `computeAveExpr`: average expression per cell group.
///
/// `data` is `t(data.use)`, cells x genes column-major. Returns the flat buffer in R's
/// genes x groups order together with the feature and group names.
#[extendr]
pub fn compute_ave_expr(
    data: Vec<f64>,
    genes: Vec<String>,
    group: Vec<i32>,
    group_levels: Vec<String>,
    features: Vec<String>,
    type_: &str,
    trim: f64,
) -> List {
    if genes.is_empty() {
        extendr_api::throw_r_error("computeAveExpr: no genes in the expression matrix");
    }
    if data.len() % genes.len() != 0 {
        extendr_api::throw_r_error(format!(
            "computeAveExpr: data has {} values, not a multiple of {} genes",
            data.len(),
            genes.len()
        ));
    }
    let groups = unwrap_r(group_index(&group, &group_levels));
    if groups.len() != data.len() / genes.len() {
        extendr_api::throw_r_error(format!(
            "computeAveExpr: group has {} entries but the matrix has {} cells",
            groups.len(),
            data.len() / genes.len()
        ));
    }
    let feats: Option<&[String]> = if features.is_empty() { None } else { Some(&features) };
    let (values, names) = r_core::de::select_features(&data, &genes, feats);
    let avg = unwrap_r(
        r_core::de::compute_ave_expr(
            &data,
            &genes,
            &groups,
            group_levels.len(),
            feats,
            type_,
            trim,
        )
        .map_err(|e| extendr_api::Error::Other(e.to_string())),
    );
    debug_assert_eq!(values.len(), avg.len());
    list! {
        values=r!(avg),
        dim=r!(vec![names.len() as i32, group_levels.len() as i32]),
        features=r!(names),
        groups=r!(group_levels),
    }
}

/// `subsetData`'s gene list: `intersect(gene.use, rownames(data))`, in R's `intersect`
/// order. `features` empty means "use the database's own gene list".
#[extendr]
pub fn subset_data_gene_use(
    gene_use_input: Vec<String>,
    data_rownames: Vec<String>,
    features: Vec<String>,
) -> Vec<String> {
    let feats: Option<&[String]> = if features.is_empty() { None } else { Some(&features) };
    r_core::de::subset_data_gene_use(&gene_use_input, &data_rownames, feats)
}

/// `subsetDB`'s row selection on the `annotation` key, and its default `search`.
///
/// Returns the kept interaction indices, in database order, so the caller can subset its
/// own copy without a second ordering.
#[extendr]
pub fn subset_db_by_annotation(
    annotations: Vec<String>,
    search: Vec<String>,
    key_is_annotation: bool,
    non_protein: bool,
) -> List {
    if !key_is_annotation {
        extendr_api::throw_r_error(
            "Each element of the `key` should be one of the column names of the \
             interaction_input from CellChatDB",
        );
    }
    // An empty `search` selects **nothing**, which is what upstream does for an explicit
    // `search = c()`. The `NULL -> default` substitution is R's `is.null(search)` test and
    // happens in the shim, because only R can distinguish `NULL` from `character(0)`.
    let search: Vec<String> = search;
    // Naming it in `search` flips `non_protein` to TRUE upstream, before the exclusion step.
    let effective_np = non_protein || r_core::de::search_implies_non_protein(&search);
    let keep: Vec<i32> = annotations
        .iter()
        .enumerate()
        .filter(|(_, a)| {
            (effective_np || a.as_str() != "Non-protein Signaling")
                && search.iter().any(|s| s == *a)
        })
        .map(|(i, _)| i as i32 + 1)
        .collect();
    // `keep` is **1-based**. R data.frames are 1-based, and `d[c(0, 1, ..., n)]` silently
    // *drops* the 0 and shifts every index, so returning the 0-based positions the core
    // works in loses a row and returns a table that is one short of correct. This went
    // unnoticed for a long time because the differential test's reference for `subsetDB`
    // was this shim itself: `cellchatrs_upstream_env()` did not source `database.R`, so
    // `get("subsetDB", envir = env)` fell through to the cellchatrs namespace and
    // "upstream" and "Rust" were the same 2238 rows.
    list! {
        keep=r!(keep),
        search=r!(search),
        non_protein=r!(effective_np),
    }
}

// --------------------------------------------------------------------------------- net.rs

/// `aggregateNet` with no filters: returns `net$count` and `net$weight`, both `k x k`
/// column-major with `dimnames` in `group_levels` order.
#[extendr]
pub fn aggregate_net(
    prob: Vec<f64>,
    pval: Vec<f64>,
    group_levels: Vec<String>,
    interaction_names: Vec<String>,
    thresh: f64,
) -> List {
    let net = r_core::net::Net::new(prob, pval, group_levels.clone(), interaction_names);
    let (count, weight) = r_core::net::aggregate_net_default(&net, thresh);
    let k = group_levels.len() as i32;
    let (r, c, _) = net.dimnames();
    list! {
        count=count,
        weight=weight,
        dim=r!(vec![k, k]),
        // Unnamed: R expects list(rownames, colnames). A named list gives
        // `dimnames(count)` the wrong names, which `identical()` catches but a value
        // comparison does not.
        dimnames=r!(list! { r, c }),
    }
}

/// `subsetCommunication` for `slot.name = "net"`, `mode = "single"`.
///
/// Returns a list of parallel column vectors plus `colnames`, `nrow` and the source/target
/// factor levels; the R shim assembles the data.frame from it. The L-R table arrives as
/// parallel string vectors, matching [`compute_commun_prob`]'s convention, and `""` stands
/// for `NA` (upstream's `LRsig` has `evidence = NA` by default).
///
/// # Arguments
/// * `lr_*`: one string per row of `object@LR$LRsig`, all the same length.
/// * `sources_use` / `targets_use`: empty vectors for "no filter", matching the `NULL`
///   default. Upstream resolves numeric indices to names via `cells.level[...]`; the caller
///   does that, since it is the only side that has the level order.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn subset_communication(
    prob: Vec<f64>,
    pval: Vec<f64>,
    group_levels: Vec<String>,
    interaction_names: Vec<String>,
    lr_columns_present: Vec<String>,
    lr_interaction_name: Vec<String>,
    lr_interaction_name_2: Vec<String>,
    lr_pathway_name: Vec<String>,
    lr_ligand: Vec<String>,
    lr_receptor: Vec<String>,
    lr_annotation: Vec<String>,
    lr_evidence: Vec<String>,
    thresh: f64,
    sources_use: Vec<String>,
    targets_use: Vec<String>,
) -> List {
    let net = r_core::net::Net::new(prob, pval, group_levels.clone(), interaction_names);
    let n_lr = lr_interaction_name.len();
    // Every optional column must be n_lr long whether or not the table has it, so index
    // unconditionally and let `columns` record presence.
    for (name, v) in [
        ("lr_interaction_name_2", &lr_interaction_name_2),
        ("lr_pathway_name", &lr_pathway_name),
        ("lr_ligand", &lr_ligand),
        ("lr_receptor", &lr_receptor),
        ("lr_annotation", &lr_annotation),
        ("lr_evidence", &lr_evidence),
    ] {
        if v.len() != n_lr {
            extendr_api::throw_r_error(format!(
                "{name} has {} entries but lr_interaction_name has {n_lr}",
                v.len()
            ));
        }
    }
    let cols: Vec<String> = lr_columns_present.clone();
    let mut lr: Vec<r_core::net::LrMeta> = Vec::with_capacity(n_lr);
    for i in 0..n_lr {
        let na = |v: &str| if v.is_empty() { None } else { Some(v.to_string()) };
        lr.push(r_core::net::LrMeta {
            interaction_name: lr_interaction_name[i].clone(),
            columns: cols.clone(),
            interaction_name_2: lr_interaction_name_2[i].clone(),
            pathway_name: lr_pathway_name[i].clone(),
            ligand: lr_ligand[i].clone(),
            receptor: lr_receptor[i].clone(),
            annotation: lr_annotation[i].clone(),
            evidence: na(&lr_evidence[i]),
        });
    }
    let src: Option<&[String]> = if sources_use.is_empty() { None } else { Some(&sources_use) };
    let tgt: Option<&[String]> = if targets_use.is_empty() { None } else { Some(&targets_use) };
    let t = r_core::net::subset_communication(&net, &lr, thresh, src, tgt, false, true);
    let n = t.len();
    let levels = t.group_levels.clone();
    // The interaction names, in the order the `Prob` array's third dimension had them. Upstream's
    // melt factors this column over that vector, so its factor levels are the *whole* L-R set and
    // not the subset that survived the threshold -- an interaction filtered out of the table still
    // contributes a level. `NetTable` does not carry it, so it comes off the net.
    let interaction_levels = net.interaction_names.clone();

    // `list!` needs a fixed shape, so enumerate every column `subsetCommunication` can
    // return and carry `NULL` for the ones the table does not have. R's `list` drops NULL
    // entries, which is exactly right: the shim sees only the real columns.
    macro_rules! col {
        ($name:expr) => {
            if t.columns.iter().any(|c| c == $name) {
                if $name == "prob" || $name == "pval" {
                    r!((0..n).map(|r| t.num(r, $name).unwrap_or(f64::NAN)).collect::<Vec<f64>>())
                } else {
                    r!((0..n)
                        .map(|r| t.get(r, $name).unwrap_or("NA").to_string())
                        .collect::<Vec<String>>())
                }
            } else {
                r!(extendr_api::NULL)
            }
        };
    }
    list! {
        colnames=r!(t.columns.clone()),
        nrow=r!(n as i32),
        source_levels=r!(levels.clone()),
        target_levels=r!(levels.clone()),
        interaction_levels=r!(interaction_levels),
        source=col!("source"),
        target=col!("target"),
        ligand=col!("ligand"),
        receptor=col!("receptor"),
        prob=col!("prob"),
        pval=col!("pval"),
        interaction_name=col!("interaction_name"),
        interaction_name_2=col!("interaction_name_2"),
        pathway_name=col!("pathway_name"),
        annotation=col!("annotation"),
        evidence=col!("evidence"),
    }
}

/// R's `NA_real_` bit pattern. `NA` and `NaN` are both NaN payloads, and the difference
/// between them is load-bearing throughout `subsetCommunication`: a group whose `pval` is
/// `NA` reports `NA`, a group whose `pval` is `NaN` reports `NaN`, and a threshold comparison
/// blanks the row in both cases while `mean`/`sum` do not. Reading it as "is it a NaN" would
/// conflate them.
///
/// The list a kernel error comes back as: `__error` plus a zero `__nrow`.
///
/// One helper because the shape has to be spelled out at every call site -- with
/// `vec![Some(msg)].into()` inline, `List::from_pairs`' element type is inferred from whichever
/// entry the compiler reaches first, and it reaches the `String` one and settles on
/// `Vec<Option<String>>`, which is not a `KeyValue`. Four copies of that is four chances to get
/// it subtly differently wrong.
/// A kernel error, with an opt-in detail line on stderr.
///
/// The message is a contract: `subscript out of bounds` is upstream's text, byte for byte, and
/// extendr turns the returned `Err` into that R condition. The *name* of the gene that failed to
/// resolve is genuinely useful when a run on real data raises it, and it cannot go in the message
/// without breaking parity -- so it is `Debug` output behind an environment variable, never the
/// error text.
fn kernel_err(e: &r_core::prob::KernelError) -> extendr_api::Error {
    if std::env::var_os("CELLCHATRS_DEBUG_EXPR").is_some() {
        eprintln!("cellchatrs: kernel error detail: {e:?}");
    }
    extendr_api::Error::Other(e.to_string())
}

fn kernel_error(msg: String) -> extendr_api::List {
    use extendr_api::{List, Robj};
    let pairs: Vec<(String, Robj)> = vec![
        ("__error".to_string(), vec![Some(msg)].into()),
        ("__nrow".to_string(), r!(0i32)),
    ];
    List::from_pairs(pairs)
}

/// R's `NA_real_` bit pattern, for the *output* direction: a `NA` result has to come back as
/// `NA_real_` and not as a plain `NaN`, because `subsetCommunication` distinguishes the two
/// and the corpus does too.
///
/// `extendr`'s `REALSXP` -> `Vec<f64>` conversion is a `from_bits`, not a canonicalisation, so
/// the payload survives the crossing and can be recognised by bit pattern on either side.
pub const R_NA_REAL_BITS: u64 = 0x7ff0_0000_0000_07a2;

/// `subsetCommunication` from the melt onwards: the eight DEG thresholds, the `netP`
/// aggregation and the final column selection.
///
/// Two input forms, matching upstream's `if (!is.data.frame(net))`:
///
/// * `columns` non-empty -- an already-melted `net` data frame, which is the only way the DEG
///   columns can exist at all (the melt attaches only `LR`'s six metadata columns, so with an
///   array `net` every threshold raises).
/// * `columns` empty -- melt `prob`/`pval` and join `LR` here, *unselected*, because
///   upstream's `rowSums(is.na(net)) != ncol(net)` counts over the full column set.
///
/// The result is returned as a named `list` of parallel vectors plus `colnames` and `nrow`;
/// the R side assembles the data frame. Numeric columns are returned as `REALSXP` with `NA`
/// and `NaN` preserved as their own payloads, and character columns as `STRSXP` with `NA`.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn subset_communication_deg(
    colnames: Vec<String>,
    cell_text: Vec<String>,
    n_rows: i32,
    prob: Vec<f64>,
    pval: Vec<f64>,
    group_levels: Vec<String>,
    interaction_names: Vec<String>,
    lr_columns_present: Vec<String>,
    lr_interaction_name: Vec<String>,
    lr_interaction_name_2: Vec<String>,
    lr_pathway_name: Vec<String>,
    lr_ligand: Vec<String>,
    lr_receptor: Vec<String>,
    lr_annotation: Vec<String>,
    lr_evidence: Vec<String>,
    thresh: f64,
    datasets: Vec<String>,
    ligand_pvalues: f64,
    ligand_logfc: f64,
    ligand_pct1: f64,
    ligand_pct2: f64,
    receptor_pvalues: f64,
    receptor_logfc: f64,
    receptor_pct1: f64,
    receptor_pct2: f64,
    sources_use: Vec<String>,
    targets_use: Vec<String>,
    slot_name: String,
) -> List {
    use extendr_api::List;

    let mut table = if colnames.is_empty() {
        let net = r_core::net::Net::new(prob, pval, group_levels.clone(), interaction_names);
        let n_lr = lr_interaction_name.len();
        let cols: Vec<String> = lr_columns_present.clone();
        let mut lr: Vec<r_core::net::LrMeta> = Vec::with_capacity(n_lr);
        for i in 0..n_lr {
            let na = |v: &str| if v.is_empty() { None } else { Some(v.to_string()) };
            lr.push(r_core::net::LrMeta {
                interaction_name: lr_interaction_name[i].clone(),
                columns: cols.clone(),
                interaction_name_2: lr_interaction_name_2[i].clone(),
                pathway_name: lr_pathway_name[i].clone(),
                ligand: lr_ligand[i].clone(),
                receptor: lr_receptor[i].clone(),
                annotation: lr_annotation[i].clone(),
                evidence: na(&lr_evidence[i]),
            });
        }
        r_core::net::subset_communication(&net, &lr, thresh, None, None, false, false)
    } else {
        unwrap_r(build_table(&colnames, &cell_text, n_rows))
    };

    // Upstream resolves numeric `sources.use` / `targets.use` against `cells.level` before
    // comparing, and the caller holds the level order, so the resolution stays in R.
    let src: Option<Vec<String>> = (!sources_use.is_empty()).then_some(sources_use);
    let tgt: Option<Vec<String>> = (!targets_use.is_empty()).then_some(targets_use);
    let th = r_core::subset::DegThresholds {
        datasets: (!datasets.is_empty()).then_some(datasets),
        ligand_pvalues: opt(ligand_pvalues),
        ligand_logfc: opt(ligand_logfc),
        ligand_pct1: opt(ligand_pct1),
        ligand_pct2: opt(ligand_pct2),
        receptor_pvalues: opt(receptor_pvalues),
        receptor_logfc: opt(receptor_logfc),
        receptor_pct1: opt(receptor_pct1),
        receptor_pct2: opt(receptor_pct2),
        sources_use: src,
        targets_use: tgt,
        slot: r_core::subset::NetSlot::from_name(&slot_name),
    };
    // The kernel's error travels back as a **value**, not as an `Err`.
    //
    // `extendr` 0.9 turns a `Result::Err` from an `#[extendr]` function into
    // `Error::MustNotBeNA`, whose `Display` is the fixed string "Must not be NA." with the
    // payload discarded -- so the error text is lost, and the text is the whole contract here.
    // Returning it in the list also lets the R side raise it with `stop()`, which is closer to
    // upstream than an extendr condition anyway: upstream's own `stop()` is a simpleError with
    // no call, and `conditionMessage()` is what the differential gate compares.
    table = match r_core::subset::subset_communication_deg(&table, &th) {
        Ok(t) => t,
        Err(e) => {
            let msg = e.to_string();
            return kernel_error(msg);
        }
    };
    // A dynamically named `list`, because the column set is data-dependent. `list!` cannot
    // express it, and the alternative -- a fixed 20-column list full of `NULL` -- loses the
    // distinction between "absent" and "present but all NA", which is exactly the distinction
    // `intersect()` makes.
    let n = table.len();
    let numeric: Vec<&str> = vec!["prob", "pval", "ligand.pvalues", "ligand.logFC",
        "ligand.pct.1", "ligand.pct.2", "receptor.pvalues", "receptor.logFC",
        "receptor.pct.1", "receptor.pct.2"];
    let mut pairs: Vec<(String, extendr_api::Robj)> =
        Vec::with_capacity(table.columns.len() + 3);
    for name in &table.columns {
        let col: Vec<Option<String>> = (0..n).map(|r| table.get(r, name).map(str::to_string)).collect();
        let obj = if numeric.contains(&name.as_str()) {
            // `None` is `NA` and `Some("NaN")` is `NaN`, and they are different values to R:
            // `mean(pval)` over a group containing an `NA` reports `NA`, over a group
            // containing a `NaN` reports `NaN`, and `identical()` tells them apart. So `NA`
            // goes back as `NA_real_`'s exact bit pattern rather than as a canonical NaN.
            let v: Vec<f64> = (0..n)
                .map(|r| match &col[r] {
                    None => f64::from_bits(if slot_name == "netP" && matches!(name.as_str(), "prob" | "pval") {
                        R_NA_REAL_BITS | (1_u64 << 51)
                    } else { R_NA_REAL_BITS }),
                    Some(s) if s == "NaN" => f64::NAN,
                    Some(s) => s.parse().unwrap_or(f64::NAN),
                })
                .collect();
            v.into()
        } else {
            col.into()
        };
        pairs.push((name.clone(), obj));
    }
    pairs.push(("__nrow".to_string(), r!(n as i32)));
    pairs.push(("__levels".to_string(), r!(table.group_levels.clone())));
    List::from_pairs(pairs)
}

/// `NaN` -- **not** `NA_real_` -- is the marker for "argument not supplied".
///
/// `extendr` 0.9 rejects `NA_real_` for an `f64` parameter: the generated wrapper raises
/// `Error::MustNotBeNA`, whose message is the fixed string "Must not be NA." and which is
/// raised during *argument conversion*, before the function body runs. So a sentinel of
/// `NA_real_` does not reach the kernel at all -- it aborts the call with a message that names
/// neither the argument nor the function. `NaN` converts fine, and the thresholds upstream
/// accepts are all finite, so nothing legitimate is lost.
fn opt(v: f64) -> Option<f64> {
    if v.is_nan() {
        None
    } else {
        Some(v)
    }
}

/// The data.frame form: a `nrow x ncol` character matrix, flattened column-major, with
/// `NA` as the missing marker.
///
/// A character matrix rather than parallel typed vectors because a typed round trip is a
/// second thing that can be wrong: `Doubles` is an `f64` newtype and `Strings` an `Rstr` one,
/// and the `NA`/`NaN` distinction has to survive both. Seventeen significant digits round-trip
/// every `f64` exactly, so `f64 -> 17-digit text -> f64` is lossless, and `NA` and `NaN` are
/// distinct tokens. One flat `Vec<String>` is also the one shape `extendr` converts without
/// ceremony.
///
/// The matrix is read column-major, which is R's own order, so the R side can hand over
/// `as.matrix()` output unchanged.
fn build_table(
    colnames: &[String],
    flat: &[String],
    n_rows: i32,
) -> extendr_api::Result<r_core::net::NetTable> {
    let n = n_rows.max(0) as usize;
    if n * colnames.len() != flat.len() {
        return Err(extendr_api::Error::Other(format!(
            "the `net` matrix has {} cells but {n} x {} were declared",
            flat.len(),
            colnames.len()
        )));
    }
    let mut rows: Vec<Vec<Option<String>>> = Vec::with_capacity(n);
    for r in 0..n {
        let mut row: Vec<Option<String>> = Vec::with_capacity(colnames.len());
        for c in 0..colnames.len() {
            row.push(match &flat[c * n + r] {
                // `NA_character_` arrives as the literal string "NA", which is also what R
                // prints for a missing character -- and is *not* how this function spells the
                // numeric `NaN`, so the two cannot be confused here.
                s if s.as_str() == "NA" => None,
                s => Some(s.clone()),
            });
        }
        rows.push(row);
    }
    Ok(r_core::net::NetTable {
        columns: colnames.to_vec(),
        rows,
        group_levels: Vec::new(),
        // Internal shim for the R-side `filterCommunication` path: neither the group nor the
        // interaction factor levels are rebuilt here, the R code does that from its own arguments.
        interaction_levels: Vec::new(),
    })
}

fn group_index(group: &[i32], levels: &[String]) -> extendr_api::Result<Vec<usize>> {
    let index: std::collections::HashMap<&str, usize> =
        levels.iter().enumerate().map(|(i, l)| (l.as_str(), i)).collect();
    group
        .iter()
        .map(|&g| {
            let name = levels
                .get(g as usize)
                .ok_or_else(|| extendr_api::Error::Other(format!("group index {g} out of range")))?;
            index
                .get(name.as_str()).copied()
                .ok_or_else(|| extendr_api::Error::Other(format!("unknown group {name}")))
        })
        .collect()
}

/// `aggregateNet`'s **filtered** branch: `count`, `weight` and the two (possibly different)
/// level vectors.
///
/// Returns a `0 x 0` or all-zero matrix, because that is what upstream produces -- see the
/// long note on [`r_core::net::aggregate_net_filtered`]. The shape depends on
/// `remove_isolate`, which is the observable consequence of upstream's
/// `str_split(key, "|")` regex bug.
#[extendr]
pub fn aggregate_net_filtered(
    table_source: Vec<String>,
    table_target: Vec<String>,
    table_prob: Vec<f64>,
    cells_level: Vec<String>,
    remove_isolate: bool,
) -> List {
    let n = table_source.len();
    if table_target.len() != n || table_prob.len() != n {
        extendr_api::throw_r_error(format!(
            "table columns disagree: source {}, target {}, prob {}",
            n,
            table_target.len(),
            table_prob.len()
        ));
    }
    let t = r_core::net::NetTable {
        columns: vec!["source".into(), "target".into(), "prob".into()],
        interaction_levels: Vec::new(),
        rows: (0..n)
            .map(|i| {
                vec![
                    Some(table_source[i].clone()),
                    Some(table_target[i].clone()),
                    Some(format!("{:?}", table_prob[i])),
                ]
            })
            .collect(),
        group_levels: cells_level.clone(),
    };
    let (count, weight, src_levels, tgt_levels) =
        r_core::net::aggregate_net_filtered(&t, &cells_level, remove_isolate);
    let (r, c) = (src_levels.len() as i32, tgt_levels.len() as i32);
    list! {
        count=count,
        weight=weight,
        dim=r!(vec![r, c]),
        // Unnamed, matching `aggregate_net`: R expects `list(rownames, colnames)`.
        dimnames=r!(list! { src_levels, tgt_levels }),
    }
}

/// `geneInfo$Symbol` as a membership set, with `NA` entries dropped.
///
/// `extendr_api::CanBeNA::is_na` on the `&str` items is exact: it distinguishes the `NA_string_`
/// sentinel from a literal `"NA"`, so no unsafe pointer comparison is needed and a gene actually
/// named `NA` (none in the pinned DBs, but the code must not assume that) would survive.
fn symbols_as_strings(
    v: &extendr_api::Robj,
) -> extendr_api::Result<std::collections::HashMap<String, ()>> {
    use extendr_api::CanBeNA;
    let Some(iter) = v.as_str_iter() else {
        return Err(extendr_api::Error::Other(
            "symbols must be a character vector".into(),
        ));
    };
    Ok(iter
        .filter(|s| !s.is_na())
        .map(|s| (s.to_string(), ()))
        .collect())
}

fn db_from_parts_for_filter(
    complexes: std::collections::HashMap<String, Vec<String>>,
    complex_names: Vec<String>,
    n_sub: usize,
    symbols: std::collections::HashMap<String, ()>,
) -> r_core::db::Database {
    r_core::db::Database::from_parts(
        "object@DB",
        vec![],
        complexes,
        std::collections::HashMap::new(),
        complex_names,
        Vec::new(),
        n_sub,
        0,
        symbols,
    )
}

/// `filterCommunication`'s numerics.
///
/// The `cat()` messages are the R shim's job -- they need `scales::percent()` -- so this
/// returns the counts and the index sets, and the shim formats and prints. The two upstream
/// failures are returned as errors carrying upstream's own text, so the shim can raise them
/// unchanged rather than inventing wording.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn filter_communication(
    prob: extendr_api::Robj,
    group_levels: Vec<String>,
    interaction_names: Vec<String>,
    group_index: Vec<i32>,
    sample_levels: Vec<String>,
    sample_index: Vec<i32>,
    data: Vec<f64>,
    gene_names: Vec<String>,
    ligand: Vec<String>,
    receptor: Vec<String>,
    complex_names: Vec<String>,
    complex_subunits: Vec<String>,
    complex_n_subunit_cols: i32,
    symbols: extendr_api::Robj,
    min_cells: i32,
    min_samples: Option<i32>,
    rare_keep: bool,
    mean_type: &str,
    trim: f64,
) -> List {
    let n_cells = group_index.len();
    // `as.numeric(net$prob)` bit-for-bit, *including* `NA_real_`.
    //
    // The generated `Vec<f64>` conversion rejects `NA_real_` during argument conversion with
    // `Error::MustNotBeNA` -- before the body runs and naming neither the argument nor the
    // function. Upstream's `filterCommunication` tolerates `NA` in `prob` (a real tutorial
    // object has it): `which(apply(net$prob, 3, sum) != 0)` drops `NA`-summing slices, and the
    // three `sum(net$prob > 0)` counts come back `NA`. So the conversion here reinterprets the
    // `REALSXP` payload without a missingness check, and `r_core::filter` implements R's
    // `sum()`/`which()` NA semantics explicitly (see `count_positive`/`slice_sum`).
    //
    // `data` keeps the generated conversion -- and therefore the loud `MustNotBeNA` -- on purpose:
    // `NA` in the expression matrix is not modeled anywhere downstream (R's `max()` would return
    // `NA` where the code below folds with `f64::max`, which *ignores* NaN payloads), so a loud
    // error beats a silent divergence. The asymmetry is the contract: `prob` NA is specified,
    // `data` NA is refused.
    let prob: Vec<f64> = unwrap_r(prob.as_real_slice().map(<[f64]>::to_vec).ok_or_else(|| {
        extendr_api::Error::Other("prob must be a numeric vector".into())
    }));
    let n_genes = gene_names.len();
    if n_genes == 0 || data.len() != n_genes * n_cells {
        extendr_api::throw_r_error(format!(
            "data has {} values, expected {n_genes} genes x {n_cells} cells",
            data.len()
        ));
    }
    fn to_usize(v: &[i32], what: &str) -> extendr_api::Result<Vec<usize>> {
        v.iter()
            .map(|&x| {
                usize::try_from(x)
                    .map_err(|_| extendr_api::Error::Other(format!("{what}: negative index {x}")))
            })
            .collect()
    }
    let group_index = unwrap_r(to_usize(&group_index, "group_index"));
    let sample_index = unwrap_r(to_usize(&sample_index, "sample_index"));

    // `data <- data/max(data)` happens *before* anything else, so the caller hands over the
    // raw matrix and the scaling is done here, where its effect is visible.
    let max = data.iter().copied().fold(f64::NEG_INFINITY, f64::max);
    if !max.is_finite() || max == 0.0 {
        extendr_api::throw_r_error(
            "data/max(data) is undefined for an all-zero or non-finite matrix",
        );
    }
    let scaled: Vec<f64> = data.iter().map(|v| v / max).collect();

    let mut complexes: std::collections::HashMap<String, Vec<String>> =
        std::collections::HashMap::new();
    let n_sub = usize::try_from(complex_n_subunit_cols).unwrap_or(0);
    for (r, name) in complex_names.iter().enumerate() {
        let row: Vec<String> = (0..n_sub)
            .map(|c| {
                complex_subunits
                    .get(r * n_sub + c)
                    .cloned()
                    .unwrap_or_default()
            })
            .collect();
        complexes.insert(name.clone(), row);
    }
    // `geneInfo$Symbol` contains exactly one `NA` in the pinned human DB (audited on the tutorial
    // object: 1 of 26,827). Upstream never notices it: `extractGeneSubset` only ever asks `%in%`
    // and `intersect()` with non-`NA` needles (gene names), and R matches an `NA` table entry
    // against nothing but an `NA` needle. So `NA` entries are dropped rather than converted, which
    // is exactly equivalent for every lookup this map will ever see (see
    // `Database::is_symbol`). Inserting a placeholder instead -- `""`, `"NA"` -- would be wrong:
    // it could match a real lookup that R would not match.
    let db = db_from_parts_for_filter(
        complexes,
        complex_names.clone(),
        n_sub,
        unwrap_r(symbols_as_strings(&symbols)),
    );

    let mean = match mean_type {
        "triMean" => r_core::filter::MeanKind::TriMean,
        "truncatedMean" => r_core::filter::MeanKind::TruncatedMean,
        "thresholdedMean" => r_core::filter::MeanKind::ThresholdedMean,
        "median" => r_core::filter::MeanKind::Median,
        other => {
            extendr_api::throw_r_error(format!(
                "`type.mean` must be one of triMean, truncatedMean, thresholdedMean, median; \
                 got {other:?}"
            ))
        }
    };

    let inp = r_core::filter::FilterInput {
        prob: &prob,
        n_lr: interaction_names.len(),
        n_groups: group_levels.len(),
        group_index: &group_index,
        n_samples: sample_levels.len(),
        sample_index: &sample_index,
        data: &scaled,
        gene_names: &gene_names,
        min_cells: usize::try_from(min_cells).unwrap_or(0),
        min_samples: min_samples.map(|m| usize::try_from(m).unwrap_or(0)),
        rare_keep,
        mean,
        trim,
        ligand: &ligand,
        receptor: &receptor,
    };
    let res = unwrap_r(r_core::filter::filter_communication(&inp, &db).map_err(|e| {
        // `Error::Other` renders as a plain R error whose message is `e.toString()`, which is
        // upstream's own text, so the shim's `tryCatch` sees the same condition.
        extendr_api::Error::Other(e.to_string())
    }));
    // One entry per sample, **including samples the filter did not score**. When
    // `min.samples < 2` the whole per-sample half is skipped and `sample_excluded` is empty;
    // padding here means the R side can index `len[i]` for every `i` without an
    // out-of-bounds `NA` turning into "missing value where TRUE/FALSE needed".
    let n_samples = sample_levels.len();
    let mut sample_excluded_len: Vec<i32> =
        res.sample_excluded.iter().map(|v| v.len() as i32).collect();
    sample_excluded_len.resize(n_samples, 0);
    let mut sample_excluded: Vec<i32> = res
        .sample_excluded
        .iter()
        .flat_map(|v| v.iter().map(|&g| g as i32))
        .collect();
    sample_excluded.shrink_to_fit();
    list! {
        prob=r!(res.prob),
        // `Option<i32>`: `None` is R's `NA_integer_`, which is what upstream's
        // `sum(net$prob > 0)` returns when the buffer holds any `NA` or `NaN`. The R side does
        // arithmetic on these (`(n1 - n2)/n1`), and `NA` propagates through it exactly as in R.
        n_interaction0=r!(res.num_interaction0.map(|n| n as i32)),
        n_interaction1=r!(res.num_interaction1.map(|n| n as i32)),
        n_interaction2=r!(res.num_interaction2.map(|n| n as i32)),
        cell_excludes=r!(res.cell_excludes.iter().map(|&g| g as i32).collect::<Vec<i32>>()),
        cell_excludes_sample=r!(res.cell_excludes_sample.iter().map(|&g| g as i32).collect::<Vec<i32>>()),
        sample_excluded=r!(sample_excluded),
        sample_excluded_len=r!(sample_excluded_len),
        lr_nonzero=r!(res.lr_nonzero.iter().map(|&l| l as i32).collect::<Vec<i32>>()),
        ran_sample_filter=r!(res.ran_sample_filter),
        min_samples_ok=r!(res.min_samples_ok),
        n_samples=r!(n_samples as i32),
        nrow=r!(group_levels.len() as i32),
    }
}

/// `rankNet`'s numeric core: the per-pathway information flow and its `-1/log` rescaling.
///
/// The plot stays in R (`PLAN.md` §14.1); what crosses the boundary is `prob`, `pval` and the
/// two group filters, and what comes back is the pair of vectors upstream calls `pSum.original`
/// and `pSum`.
///
/// # Arguments
/// * `prob` / `pval`: the `slot.name` arrays in column-major order, `k * k * n`.
/// * `group_levels` / `names`: the two `dimnames` axes. `names` is the pathway (or L-R) axis.
/// * `sources_use` / `targets_use`: **0-based** indices, already resolved from names by the
///   caller -- upstream resolves numeric `sources.use` against `cells.level` and validates
///   character `sources.use` against `dimnames(prob)[[1]]`, and the R side owns both the level
///   order and the error message.
/// * `measure`: `"weight"` or `"count"`, already through `match.arg`.
///
/// Returns `pSum_original`, `pSum`, `flagged` (the 0-based `which(is.infinite | < 0)` indices) and
/// `order`, the 0-based permutation `order(pSum.original)`. The two error paths come back as
/// `__error` rather than an `Err`, for the reason documented on `subset_communication_deg`.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn ranknet_information_flow(
    prob: Vec<f64>,
    pval: Vec<f64>,
    group_levels: Vec<String>,
    names: Vec<String>,
    thresh: f64,
    measure: String,
    sources_use: Vec<i32>,
    targets_use: Vec<i32>,
) -> List {
    use extendr_api::List;
    let n = names.len();
    let k = (prob.len() / n.max(1)).isqrt();
    if k * k * n != prob.len() || pval.len() != prob.len() {
        // R's own shape error, from `dim(prob) <- c(k, k, n)`. A length that is not `k*k*n`
        // is not a `rankNet` condition at all -- it cannot happen through the exported function,
        // because `prob` always comes from a slot -- so the port reports what R would.
        return kernel_error(
            "length of 'dimnames' [3] not equal to array extent".to_string(),
        );
    }
    let src: Option<Vec<usize>> = (!sources_use.is_empty())
        .then(|| sources_use.iter().map(|&i| i.max(0) as usize).collect());
    let tgt: Option<Vec<usize>> = (!targets_use.is_empty())
        .then(|| targets_use.iter().map(|&i| i.max(0) as usize).collect());
    let m = r_core::ranknet::Measure::from_name(&measure);
    let Some(m) = m else {
        return kernel_error(
            "'arg' should be one of \"weight\", \"count\"".to_string(),
        );
    };
    match r_core::ranknet::information_flow(
        &prob, &pval, k, n, &group_levels, thresh, m, src.as_deref(), tgt.as_deref(),
    ) {
        Ok(flow) => {
            let ord = r_core::ranknet::order_f64(&flow.original);
            List::from_pairs(vec![
                (
                    "pSum_original".to_string(),
                    flow.original.to_vec().into(),
                ),
                (
                    "pSum".to_string(),
                    flow.scaled.to_vec().into(),
                ),
                // 0-based here, as everywhere in `r-core`; the R side adds 1 when it needs
                // `which()`-style indices.
                (
                    "flagged".to_string(),
                    flow.flagged.iter().map(|&i| i as f64).collect::<Vec<f64>>().into(),
                ),
                ("order".to_string(), ord_to_r(&ord)),
            ])
        }
        Err(e) => kernel_error(e.to_string()),
    }
}

/// A 0-based permutation as a 1-based double vector, which is what R's `order()` returns and
/// therefore what the corpus records. Double rather than integer because the caller indexes
/// with it directly, and an integer vector would need a second conversion on the R side.
fn ord_to_r(ord: &[usize]) -> Robj {
    let v: Vec<f64> = ord.iter().map(|&i| (i + 1) as f64).collect();
    v.into()
}

/// `trimmed_mean`, exposed so the R shim can reuse the one implementation.
///
/// `mean(x, trim =, na.rm = )` has a specific shape -- `f <- floor(trim * n + 1/2)` elements off
/// each end, `sort`ed, and `n` recomputed *before* the mean -- and `fquantile`'s type-7 rules are
/// elsewhere in this crate. Sharing it avoids a second, subtly different copy in R.
#[extendr]
pub fn spatial_trimmed_mean(x: Vec<f64>, trim: f64, na_rm: bool) -> f64 {
    r_core::spatial::trimmed_mean(&x, trim, na_rm)
}

/// `fdist`, exposed for the `computeCellDistance` shim.
///
/// `collapse::fdist` returns a **full** `d x d` matrix, whereas upstream's `computeCellDistance`
/// returns a `dist` object -- so the R side, not here, is what rewraps it. Doing the rewrap in Rust
/// would mean inventing the `Lower`/`Upper`/`Labels` attributes in the binding layer, which is
/// where they do not belong.
#[extendr]
pub fn spatial_fdist(x: Vec<f64>, n: i32, d: i32) -> Vec<f64> {
    let (n, d) = (n as usize, d as usize);
    if n.checked_mul(d) != Some(x.len()) {
        extendr_api::throw_r_error("`coordinates` must have `nrow * ncol` entries");
    }
    r_core::spatial::fdist(&x, n, d)
}

/// `computeRegionDistance`, with the exact k-d tree in place of `AnnoyParam`.
///
/// Marshalling only; every decision lives in `r_core::spatial::region_distance`.
///
/// The R side owns the level order (`levels(factor)`) and the sample order, because
/// `nlevels(group)` and `level.use` disagree in a way the caller has to be able to see, and it
/// passes 0-based indices into them. `ratio` and `tol` come through as plain vectors and are
/// recycled by the core exactly as R's `ratio[k]` recycles them.
///
/// `interaction_range` and `contact_range` are `Option<f64>`: upstream's `dist - NULL < tol`
/// selecting nothing is load-bearing, and `f64::NAN` cannot express "absent" because `NaN` is
/// also a legitimate coordinate value in a fixture.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn spatial_region_distance(
    coords: Vec<f64>,
    n: i32,
    d: i32,
    group_levels: Vec<String>,
    group_index: Vec<i32>,
    samples_levels: Vec<String>,
    sample_index: Vec<i32>,
    interaction_range: Option<f64>,
    ratio: Vec<f64>,
    tol: Vec<f64>,
    k_min: i32,
    contact_dependent: bool,
    contact_range: Option<f64>,
    contact_knn_k: Option<i32>,
    do_symmetric: bool,
) -> List {
    use extendr_api::List;
    let (n, d) = (n as usize, d as usize);
    if n.checked_mul(d) != Some(coords.len()) {
        extendr_api::throw_r_error("`coordinates` must have `nrow * ncol` entries");
    }
    if group_index.len() != n || sample_index.len() != n {
        extendr_api::throw_r_error("`meta` must have one row per coordinate row");
    }
    let k = group_levels.len();
    if k == 0 {
        extendr_api::throw_r_error("`meta$group` must have at least one level");
    }
    let inp = r_core::spatial::RegionInput {
        coords: &coords,
        n,
        d,
        group_levels: &group_levels,
        group_index: &group_index.iter().map(|&i| i.max(0) as usize).collect::<Vec<_>>(),
        samples_levels: &samples_levels,
        sample_index: &sample_index.iter().map(|&i| i.max(0) as usize).collect::<Vec<_>>(),
        interaction_range,
        ratio: &ratio,
        tol: &tol,
        k_min: k_min.max(0) as usize,
        contact_dependent,
        contact_range,
        contact_knn_k: contact_knn_k.map(|v| v.max(0) as usize),
        do_symmetric,
    };
    match r_core::spatial::region_distance(&inp) {
        // The two error paths come back as `__error`, like every other kernel here: `extendr` 0.9
        // drops the text of an `Err(String)`, so a real `Err` would reach R as an empty message.
        Err(e) => kernel_error(e.to_string()),
        Ok(res) => List::from_pairs(vec![
            ("d_spatial".to_string(), res.d_spatial.into()),
            ("adj_contact".to_string(), res.adj_contact.into()),
            ("adj_contact_knn".to_string(), res.adj_contact_knn.into()),
            ("group_levels".to_string(), res.group_levels.into()),
        ]),
    }
}

/// `rankNet(mode = "comparison")`'s numeric core: the per-comparison flow, the **pooled**
/// degenerate reassignment, and the union of pathway vocabularies.
///
/// Marshalling only. The data-frame assembly stays in R because it is name-based bookkeeping with
/// no arithmetic -- `df[[i]][pair.name[[i]], 2] <- ...`, `factor()` of `name`,
/// `rev(levels(group))`, and `rbind`'s row-name uniquification -- and R's own `rbind` already
/// implements the uniquification exactly. Reimplementing it here would be a second thing to get
/// wrong for no gain.
///
/// `prob` / `pval` are one flat column-major `k * k * n` array **per comparison**, and
/// `pair_names` the `dimnames(prob)[[3]]` per comparison. Every per-comparison record comes back
/// as a separate list element with an `n` suffix, because one flattened concatenation would leave
/// the R side guessing at boundaries.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn ranknet_comparison_flow(
    prob: Vec<f64>,
    prob_sizes: Vec<i32>,
    pval: Vec<f64>,
    pval_sizes: Vec<i32>,
    k: i32,
    pair_names: Vec<String>,
    pair_name_sizes: Vec<i32>,
    measure: String,
    thresh: f64,
    sources_use: Vec<i32>,
    targets_use: Vec<i32>,
) -> List {
    use extendr_api::List;
    let m = r_core::ranknet::Measure::from_name(&measure);
    let Some(m) = m else {
        return kernel_error(
            "'arg' should be one of \"weight\", \"count\"".to_string(),
        );
    };
    // The list-of-vectors inputs are flattened on the way in and re-split here, because
    // `extendr` has no `TryFrom<Robj>` for `Vec<Vec<_>>`. Sizes rather than a delimiter: a
    // pathway name can contain anything except the separator, so a length vector is the only
    // encoding that cannot be broken by a fixture -- and `subset_communication_deg` sets the
    // precedent for passing text plus a row count.
    let split = |flat: Vec<f64>, sizes: Vec<i32>, what: &str| -> Option<Vec<Vec<f64>>> {
        let mut out = Vec::with_capacity(sizes.len());
        let mut at = 0usize;
        for n in sizes {
            let n = n.max(0) as usize;
            if at + n > flat.len() {
                return None;
            }
            out.push(flat[at..at + n].to_vec());
            at += n;
        }
        if at != flat.len() {
            return None;
        }
        let _ = what;
        Some(out)
    };
    let split_names = |flat: Vec<String>, sizes: Vec<i32>| -> Option<Vec<Vec<String>>> {
        let mut out = Vec::with_capacity(sizes.len());
        let mut at = 0usize;
        for n in sizes {
            let n = n.max(0) as usize;
            if at + n > flat.len() {
                return None;
            }
            out.push(flat[at..at + n].to_vec());
            at += n;
        }
        if at != flat.len() {
            return None;
        }
        Some(out)
    };
    let (Some(prob), Some(pval), Some(pair_names)) = (
        split(prob, prob_sizes, "prob"),
        split(pval, pval_sizes, "pval"),
        split_names(pair_names, pair_name_sizes),
    ) else {
        return kernel_error(
            "`prob`, `pval` and `pair_names` must be given with matching size vectors".to_string(),
        );
    };
    let ncomp = prob.len();
    if ncomp == 0 || pval.len() != ncomp || pair_names.len() != ncomp {
        return kernel_error("the compared networks must agree in number".to_string());
    }
    let k = k.max(0) as usize;
    let src: Option<Vec<usize>> =
        (!sources_use.is_empty()).then(|| sources_use.iter().map(|&i| i.max(0) as usize).collect());
    let tgt: Option<Vec<usize>> =
        (!targets_use.is_empty()).then(|| targets_use.iter().map(|&i| i.max(0) as usize).collect());
    let inp = r_core::ranknet::ComparisonIn {
        prob: &prob,
        pval: &pval,
        k,
        pair_names: &pair_names,
        measure: m,
        thresh,
        sources_use: src.as_deref(),
        targets_use: tgt.as_deref(),
    };
    match r_core::ranknet::comparison_flow(&inp) {
        Err(e) => kernel_error(e.to_string()),
        Ok(res) => {
            // `List::from_pairs` rather than `set`: `extendr`'s `List` has no `set`, so a
            // dynamically-named element set has to be built as a pair list up front. `c{i}_...`
            // rather than `..._i` so the names sort in comparison order, which is what makes the
            // R side's `names(out)` readable.
            let mut pairs: Vec<(String, extendr_api::Robj)> = vec![
                ("n".to_string(), (ncomp as i32).into()),
                ("pair_names_all".to_string(), res.pair_names_all.clone().into()),
                // The *unformatted* ratio. `format(x, digits = 1)` is R's pretty-printer and is
                // vector-dependent -- `format(1.2, digits = 1)` is "1" but
                // `format(c(0.04, 1.2), digits = 1)` is c("0.04", "1.20") because the pair shares
                // a scale -- so it, and the `order()` that consumes its result, are applied in R
                // where they are exact by construction. The kernel returns the arithmetic.
                ("n_ratio".to_string(), (res.ratio.len() as i32).into()),
            ];
            for i in 0..ncomp {
                let base = format!("c{i}");
                pairs.push((format!("{base}_original"), res.original[i].clone().into()));
                pairs.push((format!("{base}_scaled"), res.scaled[i].clone().into()));
                // `which()` is 1-based upstream and the flag is absent (not empty) for
                // `measure = "count"`, so an empty vector means "no flagged entries".
                pairs.push((format!("{base}_flagged"),
                            res.flagged[i].iter().map(|x| (*x as f64) + 1.0).collect::<Vec<f64>>().into()));
                pairs.push((format!("{base}_names"), res.pair_names_here[i].clone().into()));
            }
            for (j, r) in res.ratio.iter().enumerate() {
                pairs.push((format!("ratio{j}"), r.clone().into()));
            }
            List::from_pairs(pairs)
        }
    }
}

/// `rankNetPairwise`'s ordering, for every `(i, j)` group pair at once.
///
/// Upstream's `rankNetPairwise` is
/// ```r
/// data[with(data, order(pval, -prob)), ]
/// ```
/// inside an `n x n` loop, and the whole body is that one `order` plus data-frame construction.
/// `order` with two numeric keys and the default `method = "auto"` is a **stable** radix sort, so
/// the ordering is a total function of `(pval, -prob)` and the input index -- and
/// [`r_core::ranknet::order_f64_multi`] is already a faithful transcription of it. Reproducing it
/// here rather than calling R's `order` is what lets the shim own the function instead of passing
/// through, and it is the only arithmetic in the body.
///
/// Everything else in upstream's function -- `data.frame()`, its column order, `row.names`, the
/// nested `list()` and the `names(temp) <- colnames(prob)` -- is R-ism and stays in R, for the same
/// reason `format` and `rankNet`'s row assembly do.
///
/// The result is one 1-based permutation per group pair, in column-major `(i, j)` order to match
/// R's `array` storage, so the shim can index it with a single counter.
#[extendr]
pub fn ranknet_pairwise_orders(
    prob: Vec<f64>,
    prob_sizes: Vec<i32>,
    pval: Vec<f64>,
    pval_sizes: Vec<i32>,
) -> List {
    use extendr_api::List;
    // `sizes` is R's `dim()`, not a list of slice lengths: a `k x k x n` array flattens to
    // `k * k` consecutive slices of `n` values each, and reading `dim()` as slice lengths sums to
    // `k + k + n` rather than `k * k * n`. That mistake returns the `__error` list, which the shim
    // treats as "kernel declined" and answers from upstream -- so the differential test passes
    // while nothing is ported. Hence the explicit shape check.
    // `sizes` is R's `dim()`. A `k x k x n` array flattens **column-major**, so the first index
    // varies fastest and a fixed `(i, j)` slice is a *strided* gather with stride `k * k`, not a
    // contiguous run:
    //
    // ```text
    // a <- array(1:8, c(2,2,2)); as.vector(a)
    // # 1 2 3 4 5 6 7 8 == a[1,1,1] a[2,1,1] a[1,2,1] a[2,2,1] a[1,1,2] ...
    // ```
    //
    // So `flat.chunks_exact(n)` -- the obvious reading -- returns the run that varies `i`, and every
    // slice after the first is transposed. That version passed a 1-slice test 200 times and failed
    // 539 of 540 slices at `c(3,3,6)`, which is why the shape is checked rather than assumed.
    let split3 = |flat: &[f64], sizes: &[i32]| -> Option<Vec<Vec<f64>>> {
        if sizes.len() != 3 {
            return None;
        }
        let d: Vec<usize> = sizes.iter().map(|&d| d.max(0) as usize).collect();
        let (k, n) = (d[0], d[2]);
        if d[0].checked_mul(d[1])?.checked_mul(n)? != flat.len() {
            return None;
        }
        let mut out = Vec::with_capacity(k * d[1]);
        for j in 0..d[1] {
            for i in 0..k {
                // Column-major: i + j*k + l*k*k.
                out.push((0..n).map(|l| flat[i + j * k + l * k * d[1]]).collect());
            }
        }
        Some(out)
    };
    let (Some(prob), Some(pval)) = (split3(&prob, &prob_sizes), split3(&pval, &pval_sizes)) else {
        return kernel_error(
            "prob and pval must be 3-dimensional arrays whose lengths match their dimensions"
                .to_string(),
        );
    };
    if prob.len() != pval.len() {
        return kernel_error(
            "prob and pval must agree in their first two dimensions".to_string(),
        );
    }
    let mut orders = Vec::with_capacity(prob.len());
    for (p, v) in prob.iter().zip(pval.iter()) {
        // `-prob` is the second key: descending probability. Negating is not a bijection for the
        // two zeros, but `order` only ever *compares*, and `-0.0` and `0.0` compare equal either
        // way, so the permutation is unchanged.
        let neg: Vec<f64> = p.iter().map(|x| -x).collect();
        let ord = r_core::ranknet::order_f64_multi(&[v, &neg]);
        orders.push(ord.into_iter().map(|i| (i + 1) as i32).collect::<Vec<_>>());
    }
    use extendr_api::Robj;
    let pairs: Vec<(String, Robj)> = orders
        .into_iter()
        .enumerate()
        .map(|(i, v)| (i.to_string(), Robj::from(v)))
        .collect();
    List::from_pairs(pairs)
}

/// The deterministic half of `computeCentralityLocal`: unweighted degrees, strengths and
/// betweenness, from one K×K pathway slice.
///
/// `net` is the slice in **row-major** order -- the caller passes `as.numeric(t(slice))`, because
/// `as.numeric()` on a matrix is column-major and the core indexes `m[row * k + col]`. Upstream's
/// whole computation is orientation-agnostic (it builds an igraph object, not an index
/// expression), so the transpose is unobservable downstream; getting it wrong here would move every
/// cell to a different place and the k-d tree... no, the centrality -- would answer confidently
/// about the wrong graph. Stated rather than assumed, as elsewhere.
///
/// Errors travel back as a **value** (`__error`), not as `Err`, for the reason documented on
/// `subset_communication_deg`: extendr 0.9 renders `Err` as "Must not be NA", discarding the text,
/// and the text -- igraph's exact message down to the `Source:` suffix -- is the contract. The
/// tiny-weights condition travels back as a flag for the same reason a boolean travels better than
/// a side effect: the R side raises igraph's exact warning text, and a warning raised from Rust
/// would carry the wrong call.
#[extendr]
pub fn centrality_deterministic(
    net: Vec<f64>,
    k: i32,
) -> List {
    let k = k as usize;
    if net.len() != k * k {
        return kernel_error(format!(
            "centrality: net has {} values, expected {k}x{k}",
            net.len()
        ));
    }
    let det = match r_core::centrality::degrees_strength(&net, k) {
        Ok(d) => d,
        Err(e) => return kernel_error(e.to_string()),
    };
    let (btw, tiny) = match r_core::centrality::betweenness(&net, k) {
        Ok(v) => v,
        Err(e) => return kernel_error(e.to_string()),
    };
    list! {
        outdeg_unweighted=r!(det.outdeg_unweighted),
        indeg_unweighted=r!(det.indeg_unweighted),
        outdeg=r!(det.outdeg),
        indeg=r!(det.indeg),
        betweenness=r!(btw),
        tiny_weights=r!(tiny),
    }
}
