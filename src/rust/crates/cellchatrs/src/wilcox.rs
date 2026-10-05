// --------------------------------------------------------------------------------- wilcox.rs
//
// extendr bindings for `identifyOverExpressedGenes(..., do.fast = FALSE)`.
//
// Argument marshalling only; the numerics are `r_core::wilcox`. The shim in
// `R/modeling.R` keeps the S4 side (`object@var.features[[features.name]]` and the
// `.info` data.frame), the `group.by`/`idents.use`/`invert`/`features` argument
// resolution, and the "replicate upstream" fallback.

/// `identifyOverExpressedGenes(..., do.fast = FALSE)` for the `do.DE = TRUE` branch.
///
/// # Arguments
/// * `data` — `features.use x cells`, column-major doubles, already subset by the shim
///   (`data.use <- X[features.use, ]`).
/// * `features` — `rownames(data)`, i.e. `features.use` in its own order. Length must be
///   `length(data) / length(group_index)`.
/// * `group_index` — per-cell level index into `group_levels`, **0-based**, matching
///   `compute_ave_expr`'s `group`. (0-based, not 1-based: the shared `group_index()` helper
///   indexes `group_levels` directly, so a 1-based vector raises "group index N out of
///   range" on the last level -- which for a two-group fixture is the *common* case, so it
///   is not a subtle off-by-one to discover later.)
/// * `group_levels` — the factor levels, in level order, restricted by `idents.use`.
/// * `n_adjust` — `nrow(X)`, the Bonferroni multiplier. **Not** `nrow(data.use)`: upstream
///   writes `p.adjust(p, "bonferroni", n = nrow(X))` where `X` is the full
///   `object@data.signaling`, so the two differ whenever the caller passed `features`.
#[extendr]
#[allow(clippy::too_many_arguments)]
pub fn identify_over_expressed_genes(
    data: Vec<f64>,
    features: Vec<String>,
    group_index: Vec<i32>,
    group_levels: Vec<String>,
    thresh_pc: f64,
    thresh_fc: f64,
    thresh_p: f64,
    only_pos: bool,
    n_adjust: i32,
) -> List {
    let n_genes = features.len();
    let n_cells = group_index.len();
    if n_genes == 0 {
        extendr_api::throw_r_error(
            "`features` is empty; upstream would have stopped at \
             `nrow(X) < 3` or produced a zero-row `features` vector",
        );
    }
    if data.len() != n_genes * n_cells {
        extendr_api::throw_r_error(format!(
            "data has {} values, expected features({}) x cells({}) = {}",
            data.len(),
            n_genes,
            n_cells,
            n_genes * n_cells
        ));
    }
    let n_groups = group_levels.len();
    if n_groups < 2 {
        extendr_api::throw_r_error(format!(
            "need at least 2 groups with cells, got {n_groups}"
        ));
    }
    let labels = unwrap_r(self::group_index(&group_index, &group_levels));
    let res = r_core::wilcox::identify_over_expressed(
        &data,
        &features,
        &labels,
        n_groups,
        thresh_pc,
        thresh_fc,
        thresh_p,
        only_pos,
        n_adjust as usize,
    );
    // `clusters` is `level.use[i]`, not an index, and the levels may have been filtered by
    // `idents.use`.
    let markers = res.rows;
    let clusters: Vec<String> = markers
        .iter()
        .map(|m| group_levels[m.group_index].clone())
        .collect();
    let n = markers.len();
    list! {
        clusters=r!(clusters),
        features=r!(markers.iter().map(|m| m.feature.clone()).collect::<Vec<String>>()),
        pvalues=r!(markers.iter().map(|m| m.pvalue).collect::<Vec<f64>>()),
        logFC=r!(markers.iter().map(|m| m.log_fc).collect::<Vec<f64>>()),
        // `list!` needs Rust identifiers, so the dotted upstream column names
        // (`pct.1`, `pct.2`, `pvalues.adj`) travel as `pct1`/`pct2`/`pvalues_adj` and the
        // shim names the data.frame columns. Renaming on the R side rather than
        // working around the macro keeps the upstream spelling in one place.
        pct1=r!(markers.iter().map(|m| m.pct_1).collect::<Vec<f64>>()),
        pct2=r!(markers.iter().map(|m| m.pct_2).collect::<Vec<f64>>()),
        pvalues_adj=r!(markers.iter().map(|m| m.pvalue_adj).collect::<Vec<f64>>()),
        nrow=r!(n as i32),
        // > 0 means the full seven-column schema; 0 means upstream's collapsed 0x1
        // `features`-only frame. See `r_core::wilcox::GroupMarkers`.
        n_before_only_pos=r!(res.n_before_only_pos as i32),
    }
}

/// [`identify_over_expressed_genes`], with upstream's `group.dataset` cell selection.
///
/// ```r
/// if (is.null(group.dataset)) {
///   cell.use1 <- which(labels == level.use[i])
///   cell.use2 <- base::setdiff(1:length(labels), cell.use1)
/// } else if ((!is.null(group.dataset)) & (group.DE.combined == FALSE)) {
///   cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
///   cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
/// } else if ((!is.null(group.dataset)) & (group.DE.combined == TRUE)) {
///   cell.use1 <- which(labels.dataset == pos.dataset)
///   cell.use2 <- which(labels.dataset != pos.dataset)
/// }
/// ```
///
/// `dataset_index` is 0-based, matching `group_index`, and identifies the **positive** dataset;
/// everything else is "not the positive dataset". Upstream does not pass a second label through
/// at all, because it has already collapsed every non-positive dataset into a single level:
///
/// ```r
/// labels.dataset[labels.dataset != pos.dataset] <- toString(setdiff(unique(labels.dataset), pos.dataset))
/// ```
///
/// so "the rest" is one group whatever the original number of datasets was. Index `0` is
/// therefore all the caller needs, and `group_de_combined` selects between pooling cells across
/// groups and restricting to the current one.
///
/// This is a **separate binding** rather than extra arguments on `identify_over_expressed_genes`:
/// extendr's generated wrappers require every non-`Option` argument, so adding parameters to a
/// function the shim already calls by name would break that call. The two share the numerics
/// through [`r_core::wilcox::identify_over_expressed_selected`], which is where the selection
/// actually lives.
#[extendr]
pub fn identify_over_expressed_genes_dataset(
    data: Vec<f64>,
    features: Vec<String>,
    group_index: Vec<i32>,
    group_levels: Vec<String>,
    dataset_index: Vec<i32>,
    thresh_pc: f64,
    thresh_fc: f64,
    thresh_p: f64,
    only_pos: bool,
    n_adjust: i32,
    group_de_combined: bool,
) -> List {
    let n_genes = features.len();
    let n_cells = group_index.len();
    if n_genes == 0 {
        extendr_api::throw_r_error(
            "`features` is empty; upstream would have stopped at \
             `nrow(X) < 3` or produced a zero-row `features` vector",
        );
    }
    if data.len() != n_genes * n_cells {
        extendr_api::throw_r_error(format!(
            "data has {} values, expected features({}) x cells({}) = {}",
            data.len(),
            n_genes,
            n_cells,
            n_genes * n_cells
        ));
    }
    if dataset_index.len() != n_cells {
        extendr_api::throw_r_error(format!(
            "dataset_index has {} entries, one per cell is {}",
            dataset_index.len(),
            n_cells
        ));
    }
    let n_groups = group_levels.len();
    if n_groups < 2 {
        extendr_api::throw_r_error(format!(
            "need at least 2 groups with cells, got {n_groups}"
        ));
    }
    let labels = unwrap_r(self::group_index(&group_index, &group_levels));
    let is_pos = |c: usize| dataset_index[c] == 0;
    // `group.DE.combined = TRUE` drops the `labels == level.use[i]` term from **both** sets, so
    // every group is tested against the same pooled cells and the `cluster` column is the only
    // thing distinguishing the rows. Keeping the term on the `in_rest` side, or on either side
    // for the combined case, makes the selection differ from upstream's and loses every marker
    // in any dataset whose labels are not aligned with the groups.
    let mut select = |level: usize| -> Option<(Vec<usize>, Vec<usize>)> {
        let in_group: Vec<usize> = (0..n_cells)
            .filter(|&c| is_pos(c) && (group_de_combined || labels[c] == level))
            .collect();
        let in_rest: Vec<usize> = (0..n_cells)
            .filter(|&c| !is_pos(c) && (group_de_combined || labels[c] == level))
            .collect();
        if in_group.is_empty() || in_rest.is_empty() {
            return None;
        }
        Some((in_group, in_rest))
    };
    let res = r_core::wilcox::identify_over_expressed_selected(
        &data,
        &features,
        n_groups,
        thresh_pc,
        thresh_fc,
        thresh_p,
        only_pos,
        n_adjust as usize,
        &mut select,
    );
    let markers = res.rows;
    let clusters: Vec<String> = markers
        .iter()
        .map(|m| group_levels[m.group_index].clone())
        .collect();
    let n = markers.len();
    list! {
        clusters=r!(clusters),
        features=r!(markers.iter().map(|m| m.feature.clone()).collect::<Vec<String>>()),
        pvalues=r!(markers.iter().map(|m| m.pvalue).collect::<Vec<f64>>()),
        logFC=r!(markers.iter().map(|m| m.log_fc).collect::<Vec<f64>>()),
        pct1=r!(markers.iter().map(|m| m.pct_1).collect::<Vec<f64>>()),
        pct2=r!(markers.iter().map(|m| m.pct_2).collect::<Vec<f64>>()),
        pvalues_adj=r!(markers.iter().map(|m| m.pvalue_adj).collect::<Vec<f64>>()),
        nrow=r!(n as i32),
        n_before_only_pos=r!(res.n_before_only_pos as i32),
    }
}

/// `identifyOverExpressedGenes(..., do.DE = FALSE)`.
///
/// ```r
/// markers.all <- data.frame(features = as.character(rownames(data.use)),
///                           nCells = rowSums(data.use > 0))
/// markers.all <- dplyr::filter(markers.all, nCells >= min.cells)
/// ```
///
/// This is the *only* branch where `min.cells` is applied; the Wilcoxon branch ignores it
/// despite the argument being in the same signature.
#[extendr]
pub fn expressed_in_cells(
    data: Vec<f64>,
    features: Vec<String>,
    min_cells: i32,
) -> List {
    let n_genes = features.len();
    if n_genes == 0 {
        extendr_api::throw_r_error("`features` is empty");
    }
    if data.len() % n_genes != 0 {
        extendr_api::throw_r_error(format!(
            "data length {} is not a multiple of {} features",
            data.len(),
            n_genes
        ));
    }
    let rows = r_core::wilcox::expressed_in_at_least_n_cells(&data, &features, min_cells as usize);
    let n = rows.len();
    list! {
        features=r!(rows.iter().map(|(f, _, _)| f.clone()).collect::<Vec<String>>()),
        /// **Integer**, not double: `rowSums` on a `Matrix` returns `integer`, so upstream's
        /// `nCells` column is an integer vector. `all.equal` calls `integer` and `double`
        /// columns equal, so only `identical()` -- the actual acceptance gate -- catches it.
        n_cells=r!(rows.iter().map(|(_, c, _)| *c as i32).collect::<Vec<i32>>()),
        /// 1-based row index in `data.use`, for the `rowSums`-names-as-row-names rule.
        row_index=r!(rows.iter().map(|(_, _, i)| (*i + 1) as i32).collect::<Vec<i32>>()),
        nrow=r!(n as i32),
    }
}

/// `computeCommunProbPathway` for the `object = NULL` branch.
///
/// # Arguments
/// * `prob` / `pval` — `k x k x nLR` column-major, last dimension fastest.
/// * `pathway_name` — one string per L-R; the pathway names are taken in first-appearance
///   order, which is R's `unique(pairLR.use$pathway_name)`.
#[extendr]
pub fn compute_commun_prob_pathway(
    prob: Vec<f64>,
    pval: Vec<f64>,
    group_levels: Vec<String>,
    interaction_names: Vec<String>,
    pathway_name: Vec<String>,
    thresh: f64,
) -> List {
    let r = r_core::pathway::compute_commun_prob_pathway(
        &prob,
        &pval,
        group_levels,
        interaction_names,
        &pathway_name,
        thresh,
    );
    let n = r.prob.len();
    list! {
        pathways=r!(r.pathways),
        lrs=r!(r.lr_sig),
        prob=r!(r.prob),
        dim=r!(vec![r.dims.0 as i32, r.dims.1 as i32, r.dims.2 as i32]),
        source_levels=r!(r.group_levels),
        nrow=r!(n as i32),
    }
}
