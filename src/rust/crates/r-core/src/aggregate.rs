//! `stats::aggregate` for a single grouping factor, plus the mean functions CellChat
//! dispatches to.
//!
//! This is the primitive underneath `computeAveExpr` and the exported-but-unused
//! `computeExprGroup_agonist` / `computeExprGroup_antagonist`.
//!
//! Upstream is `aggregate(matrix, list(group), FUN = FunMean)`, where `group` is a
//! `factor`. Three upstream details this reproduces:
//!
//! 1. **Group order is factor-level order, not order of appearance.** `split()` indexes
//!    by `levels(f)`, so a cell of group "b" appearing before any cell of group "a" still
//!    lands in the "a" row. This is R-ism 5 and it is why every port that indexes groups
//!    by first appearance is quietly wrong for any input whose levels are not already
//!    sorted.
//! 2. **Empty levels are dropped**, not back-filled. `aggregate.data.frame` computes
//!    `y <- y[match(sort(unique(grp)), grp, 0L), ]`, so the result has one row per group
//!    *actually present* and a level with no cells simply has no row -- there is no `NA`
//!    placeholder to find. [`aggregate_1`] instead returns a dense `n_groups`-row buffer
//!    with `NaN` in the empty rows. That is a deliberate, documented divergence: it is
//!    unreachable from CellChat, where `ident` is derived from the cell labels so every
//!    level has at least one cell, and a dense buffer keeps the caller's indexing trivial.
//!    `empty_groups_are_dense_nan_and_that_is_a_documented_divergence` pins it.
//! 3. **Values arrive in ascending cell order within a group**, which is what makes
//!    `fquantile`'s radix-vs-comparison sort (see [`crate::stats::collapse_uses_radix_order`])
//!    and `mean`'s two-pass accumulation order reproduce.
//!
//! Column-major throughout, matching R's storage order — which matters for
//! [`crate::longdouble::F80`] accumulations, where summation order is observable.

use crate::stats::{r_mean_no_rm, r_median, r_trimmed_mean, thresholded_mean, tri_mean};

/// The `FUN` CellChat dispatches to via `switch(aggregate.fun, ...)`, plus
/// `thresholdedMean` for completeness.
#[derive(Clone, Copy, Debug, PartialEq)]
pub enum GroupMean {
    /// `aggregate.fun = "mean"` -> `mean(x, na.rm = TRUE)`.
    Mean,
    /// `aggregate.fun = "triMean"` -> `triMean(x)`.
    TriMean,
    /// `aggregate.fun = "trimmedMean"` -> `trimmedMean(x, trim, na.rm = TRUE)`.
    TrimmedMean { trim: f64 },
    /// `thresholdedMean(x, trim)`.
    ThresholdedMean { trim: f64 },
    /// `aggregate.fun = "median"` -> `median(x, na.rm = TRUE)`.
    Median,
}

impl GroupMean {
    /// Apply the function to one group's values, in the order given.
    #[inline]
    pub fn apply(self, x: &[f64]) -> f64 {
        match self {
            GroupMean::Mean => r_mean_no_rm(x),
            GroupMean::TriMean => tri_mean(x),
            GroupMean::TrimmedMean { trim } => r_trimmed_mean(x, trim, true),
            GroupMean::ThresholdedMean { trim } => thresholded_mean(x, trim),
            GroupMean::Median => r_median(x, true),
        }
    }
}

/// `aggregate(matrix, list(group), FUN)` for one factor.
///
/// `values` is a `n_rows x src_cols` column-major buffer (so its length is
/// `n_rows * values.len() / n_rows`), `cols` names the source columns to aggregate -- all
/// of them for a full aggregation -- and `group[c]` is the level index of row `c`.
///
/// Returns a `n_groups x cols.len()` column-major buffer.
///
/// Passing a column subset matters for [`aggregate_columns`] and for matching R, where
/// `aggregate` is called on a matrix built by `t()` or `matrix(..., ncol = 1)` and so may
/// see fewer columns than the expression matrix has.
///
/// # Panics
/// If `group.len() != n_rows`, if `values.len() % n_rows != 0`, or if a selected column is
/// out of range. These are internal wiring errors; the buffers come from
/// [`crate::db::extract_gene_subset`] and the bindings.
pub fn aggregate_1(
    values: &[f64],
    n_rows: usize,
    cols: &[usize],
    group: &[usize],
    n_groups: usize,
    fun: GroupMean,
) -> Vec<f64> {
    let n_cols = cols.len();
    assert_eq!(
        values.len() % n_rows,
        0,
        "aggregate_1: value buffer is not a whole number of rows"
    );
    let n_src = values.len() / n_rows;
    assert_eq!(
        group.len(),
        n_rows,
        "aggregate_1: group length must match n_rows"
    );
    assert!(
        cols.iter().all(|&c| c < n_src),
        "aggregate_1: column out of range"
    );
    debug_assert!(
        group.iter().all(|&g| g < n_groups),
        "aggregate_1: group index out of range"
    );

    // Bucket the row indices per level once, in ascending order, so the inner loop is a
    // contiguous gather. This is also what makes the accumulation order match `split()`'s.
    let mut buckets: Vec<Vec<usize>> = vec![Vec::new(); n_groups];
    for (c, &g) in group.iter().enumerate() {
        buckets[g].push(c);
    }

    // R's result is `n_groups` rows x `n_cols` columns stored column-major, so the index
    // is `g * n_cols + j` -- the **group** strides. Writing `j * n_groups + g` here instead
    // produces the same multiset of numbers in a different order, which is why it survived
    // review and then scrambled every multi-subunit group mean.
    let mut out = vec![f64::NAN; n_groups * n_cols];
    let mut buf: Vec<f64> = Vec::new();
    for (g, rows) in buckets.iter().enumerate() {
        for (j, &src) in cols.iter().enumerate() {
            if rows.is_empty() {
                continue; // see the empty-level note on `aggregate_1`
            }
            buf.clear();
            buf.extend(rows.iter().map(|&c| values[src * n_rows + c]));
            out[g * n_cols + j] = fun.apply(&buf);
        }
    }
    out
}

/// Every column of a `n_rows x n_cols` buffer, in order.
pub fn all_cols(n_cols: usize) -> Vec<usize> {
    (0..n_cols).collect()
}

/// `aggregate_1` restricted to a set of columns, then transposed to
/// **`n_selected x n_groups`**, indexed `selected * n_groups + group`.
///
/// This is the `data.avg <- t(data.avg[, -1])` step of `computeExprGroup_agonist`: R builds
/// a `n_groups x n_subunits` intermediate and transposes it to subunits x groups, and the
/// callers index that with `[subunit, group]` -- exactly the order this returns.
///
/// ## Do not use this to feed [`crate::expr::GroupMeans`]
///
/// The two layouts differ, and the difference is silent. R stores an `nrow x ncol` matrix
/// column-major, so a genes x groups matrix has flat index `group * n_genes + gene`: the
/// **group** strides. This function's output is a genes x groups matrix whose flat index is
/// `gene * n_groups + group`: the **gene** strides. Both hold the same numbers.
///
/// [`aggregate_1`] already returns R's layout for a genes x groups result (its
/// `n_groups x n_cols` output has flat index `group * n_cols + gene`, which is
/// `GroupMeans::at(gene, group)` verbatim), so anything that wants a genes x groups matrix
/// in R order should call [`aggregate_1`] and not this. `pinned_by_a_unit_test`.
pub fn aggregate_columns(
    values: &[f64],
    n_rows: usize,
    cols: &[usize],
    group: &[usize],
    n_groups: usize,
    fun: GroupMean,
) -> Vec<f64> {
    let n_cols = cols.len();
    let sub = aggregate_1(values, n_rows, cols, group, n_groups, fun);
    // aggregate_1 gives n_groups x n_cols; the caller wants n_cols x n_groups.
    let mut out = vec![f64::NAN; n_cols * n_groups];
    for c in 0..n_cols {
        for g in 0..n_groups {
            out[c * n_groups + g] = sub[g * n_cols + c];
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn groups_are_level_order_not_appearance_order() {
        // rows: cell 0 -> group 1, cell 1 -> group 0. `split` must still put group 0 first.
        let values = [10.0, 20.0, 30.0]; // 3 rows x 1 col
        let group = [1usize, 0, 1];
        let out = aggregate_1(&values, 3, &[0], &group, 2, GroupMean::Mean);
        assert_eq!(out[0], 20.0, "group 0 is cell 1");
        assert_eq!(out[1], 20.0, "group 1 is the mean of cells 0 and 2");
    }

    #[test]
    fn empty_groups_are_dense_nan_and_that_is_a_documented_divergence() {
        // R drops the level entirely (see the module note); we return a dense buffer with
        // NaN. Unreachable from CellChat because every `ident` level has cells.
        let values = [1.0, 3.0];
        let out = aggregate_1(&values, 2, &[0], &[0, 0], 3, GroupMean::Mean);
        assert_eq!(out[0], 2.0);
        assert!(
            out[1].is_nan() && out[2].is_nan(),
            "empty levels stay NaN, never 0"
        );
    }

    #[test]
    fn single_cell_group_gets_that_cell() {
        let values = [5.0, 6.0, 7.0, 8.0];
        let out = aggregate_1(&values, 4, &[0], &[0, 1, 2, 3], 4, GroupMean::TriMean);
        assert_eq!(out, [5.0, 6.0, 7.0, 8.0]);
    }

    #[test]
    fn mean_is_column_wise_over_a_matrix() {
        // 4 rows x 2 cols, column-major: col 0 = [1,2,3,4], col 1 = [10,20,30,40]
        let values = [1.0, 2.0, 3.0, 4.0, 10.0, 20.0, 30.0, 40.0];
        let group = [0, 0, 1, 1];
        let out = aggregate_1(&values, 4, &all_cols(2), &group, 2, GroupMean::Mean);
        // 2 groups x 2 cols, column-major, group-major: index = g * n_cols + j
        assert_eq!(out[0], 1.5); // group 0, col 0 = (1+2)/2
        assert_eq!(out[1], 15.0); // group 0, col 1 = (10+20)/2
        assert_eq!(out[2], 3.5); // group 1, col 0 = (3+4)/2
        assert_eq!(out[3], 35.0); // group 1, col 1 = (30+40)/2
    }

    #[test]
    fn multi_subunit_group_means_keep_the_subunit_axis() {
        // The layout `computeExprGroup_*` depends on: 2 subunits x 3 groups, and the
        // transposed view must be subunit-major.
        let values = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]; // 6 cells x 1 col
        let group = [0, 0, 1, 1, 2, 2];
        // Pretend columns 0 and 5 hold two subunits' cell values.
        let mut m = vec![0.0; 12]; // 6 rows x 2 cols
        for c in 0..6 {
            m[c] = values[c];
            m[6 + c] = values[c] * 10.0;
        }
        let sub = aggregate_1(&m, 6, &all_cols(2), &group, 3, GroupMean::Mean);
        assert_eq!(sub.len(), 6);
        let t = aggregate_columns(&m, 6, &all_cols(2), &group, 3, GroupMean::Mean);
        // R's aggregate gives 3x2 (groups x subunits); t() gives 2x3 (subunits x groups).
        assert_eq!(t[0], 1.5); // subunit 0, group 0 = (1+2)/2
        assert_eq!(t[1], 3.5); // subunit 0, group 1 = (3+4)/2  <- a wrong stride scrambles this
        assert_eq!(t[2], 5.5); // subunit 0, group 2 = (5+6)/2
        assert_eq!(t[3], 15.0); // subunit 1, group 0
        assert_eq!(t[4], 35.0);
        assert_eq!(t[5], 55.0);
    }

    #[test]
    fn aggregate_columns_transposes() {
        // 4 rows x 2 cols; select column 1 only -> 2 groups x 1 -> returned as 1 x 2.
        let values = [1.0, 2.0, 3.0, 4.0, 10.0, 20.0, 30.0, 40.0];
        let out = aggregate_columns(&values, 4, &[1], &[0, 0, 1, 1], 2, GroupMean::Mean);
        assert_eq!(out, [15.0, 35.0]);
    }

    #[test]
    fn aggregate_1_gives_group_strides_and_aggregate_columns_gene_strides() {
        // The two layouts hold the same numbers in a different order, which is exactly why
        // the mix-up is invisible until a Prob slice comes out permuted. 2 genes x 2 groups
        // over 4 cells.
        let values = [1.0, 2.0, 3.0, 4.0, 10.0, 20.0, 30.0, 40.0]; // 4 cells x 2 genes
        let group = [0, 0, 1, 1];
        // 4 cells x 2 genes, gene 0 = [1,2,3,4] and gene 1 = [10,20,30,40]; groups {0,1,2,3}
        // pair up as (1+2)/2 = 1.5, (3+4)/2 = 3.5, (10+20)/2 = 15, (30+40)/2 = 35.
        //
        // `aggregate_1` returns groups x genes, R column-major: index = group * n_genes +
        // gene -- which is *also* a genes x groups matrix in R's order, so it feeds
        // `GroupMeans` unchanged.
        let r_order = aggregate_1(&values, 4, &all_cols(2), &group, 2, GroupMean::Mean);
        assert_eq!(
            r_order,
            [1.5, 15.0, 3.5, 35.0],
            "group 0, gene 0 / gene 1 / group 1, ..."
        );
        // `aggregate_columns` transposes to genes x groups with index = gene * n_groups +
        // group -- a *different* flat order for the same numbers.
        let t_order = aggregate_columns(&values, 4, &all_cols(2), &group, 2, GroupMean::Mean);
        assert_eq!(
            t_order,
            [1.5, 3.5, 15.0, 35.0],
            "gene 0 group 0 / group 1 / gene 1, ..."
        );
        assert_ne!(r_order, t_order, "the two orders must actually differ here");
        // Feeding the transposed order to GroupMeans silently transposes the matrix.
        let names = ["a".to_string(), "b".to_string()];
        let ok = crate::expr::GroupMeans::new(&r_order, &names, 2);
        assert_eq!(
            (ok.at(0, 0), ok.at(1, 0), ok.at(0, 1), ok.at(1, 1)),
            (1.5, 15.0, 3.5, 35.0)
        );
        let bad = crate::expr::GroupMeans::new(&t_order, &names, 2);
        assert_eq!(
            (bad.at(0, 0), bad.at(1, 0), bad.at(0, 1), bad.at(1, 1)),
            (1.5, 3.5, 15.0, 35.0)
        );
        assert_ne!(
            ok.at(1, 0),
            bad.at(1, 0),
            "and the mistake is visible at (gene 1, group 0)"
        );
    }
}
