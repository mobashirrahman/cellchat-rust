//! `computeAveExpr`, `subsetDB` and `subsetData` — the entry points that turn a raw
//! count matrix and a `CellChatDB` into the `data.signaling` the kernel consumes.
//!
//! Parity contract: **Exact** (`docs/SEMANTICS.md`).
//!
//! Upstream: `R/modeling.R:473` (`computeAveExpr`), `R/database.R:137` (`subsetDB`),
//! `R/utilities.R:311` (`subsetData`).
//!
//! ## `computeAveExpr` is `aggregate` again
//!
//! It is the same `aggregate(t(data.use), list(labels), FUN = FunMean)` the kernel
//! performs, with two differences that matter: it takes `truncatedMean` where the kernel
//! spells that `mean(x, trim = trim, na.rm = TRUE)` and `trimmedMean` is *not* offered, and
//! it assigns `rownames(data.use.avg) <- features.use` explicitly, so the row order is the
//! *requested* feature order rather than whatever `features.use` happened to be.
//!
//! ## `features.use <- intersect(features, rownames(data.use))` is load-bearing
//!
//! `intersect` returns its **first** argument's values, in first-appearance order and
//! deduplicated, keeping only those present in the second. So `computeAveExpr(features = c("B",
//! "A", "B"))` orders the output `B, A` — the caller's order, not the matrix's. That is
//! R-ism territory and it is why [`intersect_order`] takes the requested list first.

use std::collections::HashMap;

use crate::aggregate::{aggregate_1, all_cols, GroupMean};

/// `match.arg(type)` over `computeAveExpr`'s three choices.
///
/// NB the kernel offers four (`thresholdedMean` as well) and `computeAveExpr` does not, so
/// the two dispatchers are deliberately different and share only two names.
pub const AVE_EXPR_TYPES: [&str; 3] = ["triMean", "truncatedMean", "median"];

/// `match.arg` with R's partial matching, and the `trim` default.
///
/// Upstream's `trim = NULL` default is never substituted: for `truncatedMean`,
/// `mean(x, trim = NULL, na.rm = TRUE)` is R's `trim = 0`, i.e. no trimming. The default is
/// reproduced as 0.0 rather than 0.1, which is the kernel's default, because they are
/// different.
pub fn match_ave_expr_type(type_: &str, trim: f64) -> Result<GroupMean, KernelErrorLite> {
    for t in AVE_EXPR_TYPES {
        if t == type_ {
            return Ok(ave_fun(t, trim));
        }
    }
    let hits: Vec<&&str> = AVE_EXPR_TYPES
        .iter()
        .filter(|t| t.starts_with(type_))
        .collect();
    match hits.len() {
        1 => Ok(ave_fun(hits[0], trim)),
        _ => Err(KernelErrorLite::MatchArg(type_.to_string())),
    }
}

fn ave_fun(t: &str, trim: f64) -> GroupMean {
    match t {
        "triMean" => GroupMean::TriMean,
        "truncatedMean" => GroupMean::TrimmedMean { trim },
        _ => GroupMean::Median,
    }
}

/// The `trim` upstream actually uses when the caller passes the default `NULL`.
pub const AVE_TRIM_DEFAULT: f64 = 0.0;

#[derive(Debug, Clone, PartialEq)]
pub enum KernelErrorLite {
    MatchArg(String),
    UnknownKey(String),
}

impl std::fmt::Display for KernelErrorLite {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            // NB the curly quotes: R formats `match.arg`'s deparse output with
            // `sQuote()`, which uses U+201C/U+201D, and the message is compared byte for
            // byte, so straight quotes are a mismatch.
            KernelErrorLite::MatchArg(_) => write!(
                f,
                "'arg' should be one of \u{201c}triMean\u{201d}, \u{201c}truncatedMean\u{201d}, \
                 \u{201c}median\u{201d}"
            ),
            KernelErrorLite::UnknownKey(_) => write!(
                f,
                "Each element of the `key` should be one of the column names of the \
                 interaction_input from CellChatDB"
            ),
        }
    }
}

impl std::error::Error for KernelErrorLite {}

/// `intersect(a, b)`: `a`'s values, in first-appearance order, deduplicated, keeping only
/// those in `b`.
///
/// R's `intersect` is *not* symmetric and *not* sorted. `match()` on the second argument
/// gives membership; the order and the deduplication come from the first. This matters
/// because both `computeAveExpr(features = ...)` and `computeExpr_LR`'s subunit filtering
/// route through it, and getting it backwards permutes an output matrix.
pub fn intersect_order(a: &[String], b: &[String]) -> Vec<String> {
    let present: HashMap<&str, ()> = b.iter().map(|s| (s.as_str(), ())).collect();
    let mut seen: HashMap<&str, ()> = HashMap::new();
    let mut out = Vec::new();
    for s in a {
        if present.contains_key(s.as_str()) && !seen.contains_key(s.as_str()) {
            seen.insert(s.as_str(), ());
            out.push(s.clone());
        }
    }
    out
}

/// `computeAveExpr`: average expression per cell group, as a `k x n_features`
/// column-major buffer in R's flat order (group strides).
///
/// `data` is `t(data.use)`, cells x genes column-major, exactly as [`aggregate_1`] wants.
/// `features` is the requested feature list; `None` means "all, in matrix order".
///
/// # Panics
/// If `data.len() != n_cells * genes.len()` or `group.len() != n_cells`.
pub fn compute_ave_expr(
    data: &[f64],
    genes: &[String],
    group: &[usize],
    n_groups: usize,
    features: Option<&[String]>,
    type_: &str,
    trim: f64,
) -> Result<Vec<f64>, KernelErrorLite> {
    let fun = match_ave_expr_type(type_, trim)?;
    if genes.is_empty() || data.len() % genes.len() != 0 {
        panic!("compute_ave_expr: data buffer is not a whole number of cells");
    }
    let n_cells = data.len() / genes.len();
    if group.len() != n_cells {
        panic!(
            "compute_ave_expr: group has {} entries but data has {n_cells} cells",
            group.len()
        );
    }
    // NB the *selected* buffer and the *selected* count. Using the unselected `genes.len()`
    // here is a silent no-op for `features = ...`, which returns a buffer of the wrong width
    // and an average over the wrong columns -- plausible numbers, wrong answer.
    let (values, names) = select_features(data, genes, features);
    let n_feat = names.len();
    let avg = aggregate_1(&values, n_cells, &all_cols(n_feat), group, n_groups, fun);
    // `aggregate_1` returns groups x features, whose column-major flat index is
    // `group * n_features + feature` -- which is already a features x groups matrix in R's
    // order, so no transpose is needed. See the note on `aggregate_columns`.
    Ok(avg)
}

/// Row-select a cells x genes buffer by feature name, in `intersect` order.
///
/// Returns `(values, names)`, both in the selected order.
pub fn select_features(
    data: &[f64],
    genes: &[String],
    features: Option<&[String]>,
) -> (Vec<f64>, Vec<String>) {
    let n_cells = if genes.is_empty() {
        0
    } else {
        data.len() / genes.len()
    };
    match features {
        None => (data.to_vec(), genes.to_vec()),
        Some(f) => {
            let keep = intersect_order(f, genes);
            let index: HashMap<&str, usize> = genes
                .iter()
                .enumerate()
                .map(|(i, g)| (g.as_str(), i))
                .collect();
            let cols: Vec<usize> = keep
                .iter()
                .filter_map(|g| index.get(g.as_str()).copied())
                .collect();
            let mut out = vec![0.0; n_cells * cols.len()];
            for (j, &src) in cols.iter().enumerate() {
                for c in 0..n_cells {
                    out[j * n_cells + c] = data[src * n_cells + c];
                }
            }
            (out, keep)
        }
    }
}

/// `subsetDB`'s `search` default, which depends on `non_protein`.
///
/// ```r
/// if (is.null(search) & non_protein == FALSE & any(key == "annotation"))
///   search <- c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact")
/// ```
///
/// Note `Non-protein Signaling` is **absent** by default: excluding it is the point, and
/// `interaction_input <- subset(interaction_input, annotation != "Non-protein Signaling")`
/// would be redundant otherwise.
pub fn subset_db_default_search(non_protein: bool) -> &'static [&'static str] {
    const BASE: [&str; 3] = ["Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"];
    const ALL: [&str; 4] = [
        "Secreted Signaling",
        "ECM-Receptor",
        "Cell-Cell Contact",
        "Non-protein Signaling",
    ];
    if non_protein {
        &ALL
    } else {
        &BASE
    }
}

/// Was `Non-protein Signaling` requested? Upstream sets `non_protein = TRUE` and emits a
/// message if `"Non-protein Signaling" %in% unlist(search)`, *before* the exclusion step.
pub fn search_implies_non_protein(search: &[String]) -> bool {
    search.iter().any(|s| s == "Non-protein Signaling")
}

/// `subsetData`'s annotation reordering.
///
/// ```r
/// interaction_input$annotation <- factor(annotation,
///     levels = c("Secreted Signaling","ECM-Receptor","Non-protein Signaling","Cell-Cell Contact"))
/// interaction_input <- interaction_input[order(annotation), , drop = FALSE]
/// interaction_input$annotation <- as.character(annotation)
/// ```
///
/// Two R-isms, both already known but load-bearing here: `order()` on a **factor** is level
/// order, while on a character vector it is alphabetical (R-ism 9 in the DB layer), and
/// `order()` is stable so equal keys keep their relative order.
pub const ANNOTATION_LEVELS: [&str; 4] = [
    "Secreted Signaling",
    "ECM-Receptor",
    "Non-protein Signaling",
    "Cell-Cell Contact",
];

/// Rank of an annotation for `order(factor(...))`. Unknown labels sort **last** (R's
/// `factor()` gives them `NA`, and `order()` drops `NA` by default -- so upstream *removes*
/// rows with an unrecognised annotation, which is worth knowing and is why the
/// `level_one_based` helper returns `None` for them).
pub fn annotation_rank(a: &str) -> Option<usize> {
    ANNOTATION_LEVELS.iter().position(|l| *l == a)
}

/// `order(annotation)` indices: rows with a known annotation, in level order, stable
/// within a level.
pub fn order_by_annotation(annotations: &[String]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..annotations.len())
        .filter(|&i| annotation_rank(&annotations[i]).is_some())
        .collect();
    idx.sort_by_key(|&i| annotation_rank(&annotations[i]).unwrap_or(usize::MAX));
    idx
}

/// `subsetData`: the gene list for `object@data.signaling`.
///
/// `gene.use <- intersect(gene.use_input, rownames(object@data))` where `gene.use_input`
/// is `extractGene(DB)`. `features = NULL` takes the DB's list first; otherwise the
/// caller's list, again first. Then the matrix keeps rows whose name is in `gene.use`.
pub fn subset_data_gene_use(
    gene_use_input: &[String],
    data_rownames: &[String],
    features: Option<&[String]>,
) -> Vec<String> {
    match features {
        None => intersect_order(gene_use_input, data_rownames),
        Some(f) => intersect_order(f, data_rownames),
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn intersect_keeps_the_first_arguments_order_and_deduplicates() {
        let a: Vec<String> = ["B", "A", "B", "C"].iter().map(|s| s.to_string()).collect();
        let b: Vec<String> = ["A", "B", "C", "D"].iter().map(|s| s.to_string()).collect();
        assert_eq!(intersect_order(&a, &b), vec!["B", "A", "C"]);
        // Not sorted, and not b's order.
        assert_ne!(intersect_order(&a, &b), vec!["A", "B", "C"]);
    }

    #[test]
    fn intersect_drops_names_absent_from_the_second_argument() {
        let a: Vec<String> = ["A", "Z"].iter().map(|s| s.to_string()).collect();
        let b: Vec<String> = ["A"].iter().map(|s| s.to_string()).collect();
        assert_eq!(intersect_order(&a, &b), vec!["A"]);
    }

    #[test]
    fn match_arg_covers_ave_exprs_three_choices() {
        assert_eq!(match_ave_expr_type("triMean", 0.1), Ok(GroupMean::TriMean));
        assert_eq!(
            match_ave_expr_type("truncatedMean", 0.2),
            Ok(GroupMean::TrimmedMean { trim: 0.2 })
        );
        assert_eq!(match_ave_expr_type("median", 0.1), Ok(GroupMean::Median));
        // `thresholdedMean` is a kernel choice, not a computeAveExpr one.
        assert!(match_ave_expr_type("thresholdedMean", 0.1).is_err());
        // "t" is ambiguous across triMean/truncatedMean.
        assert!(match_ave_expr_type("t", 0.1).is_err());
        // "med" is a unique prefix.
        assert_eq!(match_ave_expr_type("med", 0.1), Ok(GroupMean::Median));
    }

    #[test]
    fn the_default_search_excludes_non_protein_signalling() {
        assert!(!subset_db_default_search(false).contains(&"Non-protein Signaling"));
        assert!(subset_db_default_search(true).contains(&"Non-protein Signaling"));
    }

    #[test]
    fn annotation_order_is_factor_level_order_not_alphabetical() {
        let a: Vec<String> = [
            "Cell-Cell Contact",
            "Non-protein Signaling",
            "Secreted Signaling",
            "ECM-Receptor",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        let order = order_by_annotation(&a);
        assert_eq!(
            order.iter().map(|&i| a[i].as_str()).collect::<Vec<_>>(),
            vec![
                "Secreted Signaling",
                "ECM-Receptor",
                "Non-protein Signaling",
                "Cell-Cell Contact"
            ]
        );
    }

    #[test]
    fn annotation_order_is_stable_within_a_level() {
        let a: Vec<String> = [
            "ECM-Receptor",
            "Secreted Signaling",
            "ECM-Receptor",
            "Secreted Signaling",
        ]
        .iter()
        .map(|s| s.to_string())
        .collect();
        assert_eq!(order_by_annotation(&a), vec![1, 3, 0, 2]);
    }

    #[test]
    fn an_unrecognised_annotation_is_dropped_like_order_does_with_na() {
        // `factor(x, levels = ...)` gives NA for an unknown label and `order()` drops NA by
        // default, so such a row disappears from the reordered table.
        let a: Vec<String> = ["Secreted Signaling", "Bogus"]
            .iter()
            .map(|s| s.to_string())
            .collect();
        assert_eq!(order_by_annotation(&a), vec![0]);
    }

    #[test]
    fn ave_expr_returns_a_features_by_groups_matrix_in_r_order() {
        // 2 genes x 4 cells, 2 groups; t(data.use) is cells x genes.
        let genes = vec!["A".to_string(), "B".to_string()];
        let data = [1.0, 2.0, 3.0, 4.0, 10.0, 20.0, 30.0, 40.0];
        let group = [0usize, 0, 1, 1];
        assert!(
            compute_ave_expr(&data, &genes, &group, 2, None, "thresholdedMean", 0.0).is_err(),
            "thresholdedMean is a kernel choice, not a computeAveExpr one"
        );
        let avg = compute_ave_expr(&data, &genes, &group, 2, None, "triMean", 0.0).unwrap();
        // A = [1,2,3,4] so group 0 = 1.5 and group 1 = 3.5; B = [10,20,30,40] so group 0 =
        // 15 and group 1 = 35.
        //
        // The flat index is `group * n_features + feature`, which is R's column-major order
        // for the features x groups result -- so [1.5, 15, 3.5, 35], *not* the row-major
        // [1.5, 3.5, 15, 35]. Both hold the same four numbers; only one of them is what
        // `t(data.use.avg[,-1])` is, and a matrix handed to R with the other order is
        // silently transposed.
        assert_eq!(avg, vec![1.5, 15.0, 3.5, 35.0]);
        assert_eq!(avg[0], 1.5, "group 0, feature A");
        assert_eq!(avg[1], 15.0, "group 0, feature B");
        assert_eq!(avg[2], 3.5, "group 1, feature A");
        assert_eq!(avg[3], 35.0, "group 1, feature B");
    }

    #[test]
    fn feature_selection_uses_the_callers_order_not_the_matrixs() {
        let genes = vec!["A".to_string(), "B".to_string()];
        let data = [1.0, 2.0, 10.0, 20.0];
        let want = vec!["B".to_string(), "A".to_string()];
        let (vals, names) = select_features(&data, &genes, Some(&want));
        assert_eq!(names, want);
        assert_eq!(vals, vec![10.0, 20.0, 1.0, 2.0]);
    }
}
