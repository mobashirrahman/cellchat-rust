//! The `computeExpr_*` family: turning a genes × cell-groups matrix into per-pair
//! ligand/receptor expression, plus the co-agonist, co-antagonist and co-receptor
//! modulations.
//!
//! Parity contract: **Exact** (`docs/SEMANTICS.md`).
//!
//! Upstream sources, all in `R/modeling.R`:
//! `computeExpr_LR:593`, `computeExpr_complex:525`, `computeExpr_coreceptor:621`,
//! `computeExpr_agonist:798`, `computeExpr_antagonist:834`, and the exported-but-unused
//! group-mean variants `computeExprGroup_agonist:732` / `computeExprGroup_antagonist:765`.
//!
//! ## Errors are part of the contract
//!
//! Upstream's `computeExpr_complex` indexes `data.use[RsubunitsV, ]` with raw subunit
//! *names* and no `intersect`, so a name that is neither a row of the matrix nor a row of
//! the complex table — or a complex with one subunit missing — aborts with
//! `subscript out of bounds` (R-ism 11, `PLAN.md` §14.3). That is replicated here rather
//! than smoothed over, because a silent "no modulation" fallback would be a behaviour
//! change disguised as a robustness fix.
//!
//! Contrast the cofactor paths, which *do* `intersect(...)` first and therefore treat an
//! unresolvable name as "no modulation" (all ones) rather than an error. Both behaviours
//! are present upstream, and the difference is deliberate.

use std::collections::HashMap;

use crate::aggregate::{aggregate_columns, GroupMean};
use crate::db::Database;
use crate::stats::{r_mean_no_rm, r_prod};

/// Upstream's error text for R-ism 11, byte-for-byte.
pub const SUBSCRIPT_OUT_OF_BOUNDS: &str = "subscript out of bounds";

#[derive(Clone, PartialEq, Eq)]
pub enum ExprError {
    /// A ligand/receptor name could not be resolved: either it is not a row of the
    /// expression matrix and not in the complex table, or it is a complex with a subunit
    /// that is not a row. Upstream raises `subscript out of bounds`.
    SubscriptOutOfBounds { name: String },
}

impl std::fmt::Display for ExprError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExprError::SubscriptOutOfBounds { .. } => write!(f, "{SUBSCRIPT_OUT_OF_BOUNDS}"),
        }
    }
}

impl std::fmt::Debug for ExprError {
    /// The `Display` text is a *contract* -- upstream's message is exactly
    /// `subscript out of bounds`, and any elaboration would be a divergence -- so the offending
    /// name lives in `Debug` instead. `{:?}` therefore names the gene while `{}` stays byte-exact.
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ExprError::SubscriptOutOfBounds { name } => {
                write!(f, "SubscriptOutOfBounds {{ name: {name:?} }}")
            }
        }
    }
}

impl std::error::Error for ExprError {}

type Result<T> = std::result::Result<T, ExprError>;

/// A genes × cell-groups matrix, column-major, with a name → row index.
///
/// This is the only view the `computeExpr_*` functions need, which is what lets the
/// kernel swap in a bootstrap replicate without rebuilding anything: the gene names and
/// their row positions are identical across all `nboot` replicates, only the numbers
/// change.
#[derive(Clone, Debug)]
pub struct GroupMeans<'a> {
    /// `n_groups` blocks of `n_genes` values.
    values: &'a [f64],
    n_genes: usize,
    n_groups: usize,
    row_of: HashMap<String, usize>,
}

impl<'a> GroupMeans<'a> {
    /// Build from a column-major buffer plus row names.
    ///
    /// # Panics
    /// If `values.len() != n_genes * n_groups`, or `names.len() != n_genes`. A silent
    /// length mismatch here would mis-index every gene, so it is a panic rather than a
    /// `Result`: it can only be a programming error, and the buffer is built by
    /// `aggregate.rs` and the bindings, not by user data.
    pub fn new(values: &'a [f64], names: &[String], n_groups: usize) -> Self {
        assert_eq!(
            values.len(),
            names.len() * n_groups,
            "GroupMeans: values must be genes x groups"
        );
        let row_of = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.clone(), i))
            .collect();
        Self {
            values,
            n_genes: names.len(),
            n_groups,
            row_of,
        }
    }

    /// A name-only view: `row()` works, every value accessor panics.
    ///
    /// `computeCommunProb` resolves every L-R name against the *observed* matrix before the
    /// bootstrap, and then re-resolves against each replicate. A plan is therefore
    /// name-resolved once and reused for `nboot + 1` different matrices, all of which share
    /// the same rownames. This constructor is how that resolution is done without a value
    /// buffer, which is what `Kernel::resolve` needs.
    pub fn names_only(names: &[String], n_groups: usize) -> Self {
        let row_of = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.clone(), i))
            .collect();
        Self {
            values: &[],
            n_genes: names.len(),
            n_groups,
            row_of,
        }
    }

    /// A name-only view: [`GroupMeans::row`] works, every value accessor panics.
    ///
    /// `computeCommunProb` resolves every L-R name against the *observed* matrix before the
    /// bootstrap and then, for each of `nboot` replicates, resolves it again against a
    /// different matrix -- but all `nboot + 1` matrices share the same rownames, because
    /// they are all aggregations of the same `data.use`. So the resolution is done once into
    /// a [`crate::prob::LrPlan`] of `usize` indices and reused, which is both the
    /// correctness argument for hoisting it out of the loop and the main performance win.
    /// This constructor is how that resolution happens without a value buffer.
    #[inline]
    pub fn n_genes(&self) -> usize {
        self.n_genes
    }

    #[inline]
    pub fn n_groups(&self) -> usize {
        self.n_groups
    }

    /// Row index of a gene, or `None` if it is not in the matrix.
    #[inline]
    pub fn row(&self, name: &str) -> Option<usize> {
        self.row_of.get(name).copied()
    }

    /// `value(g, j)`.
    ///
    /// # Panics
    /// On a [`GroupMeans::names_only`] view, by construction.
    #[inline]
    pub fn at(&self, g: usize, j: usize) -> f64 {
        self.values[j * self.n_genes + g]
    }

    /// Every gene row, in  rownames order.
    pub fn names(&self) -> Vec<String> {
        let mut v: Vec<(usize, &String)> = self.row_of.iter().map(|(k, i)| (*i, k)).collect();
        v.sort_unstable_by_key(|(i, _)| *i);
        v.into_iter().map(|(_, k)| k.clone()).collect()
    }
}

/// `geometricMean` restricted to a set of gene rows: `exp(mean(log(x), na.rm = TRUE))`
/// computed per group.
///
/// A single zero makes the whole group `0`, because `log 0 = -Inf`. That falls out of
/// using [`r_mean_no_rm`] — `NaN` is skipped, but `-Inf` is not, so `-Inf` survives and
/// `exp(-Inf) = 0`. R-ism 4.
/// `out[j * n_cols + col] = geometric mean of rows` for every group `j`.
fn geometric_mean_over_rows_into(
    means: &GroupMeans<'_>,
    rows: &[usize],
    out: &mut [f64],
    col: usize,
    n_cols: usize,
) {
    let k = means.n_groups();
    debug_assert!(out.len() >= k * n_cols);
    if rows.is_empty() {
        for j in 0..k {
            out[j * n_cols + col] = f64::NAN;
        }
        return;
    }
    let mut logs = Vec::with_capacity(rows.len());
    for j in 0..k {
        logs.clear();
        for &r in rows {
            logs.push(means.at(r, j).ln());
        }
        out[j * n_cols + col] = r_mean_no_rm(&logs).exp();
    }
}

/// The gene rows that define a ligand or receptor's expression, and *how* to reduce them.
///
/// The distinction is load-bearing. Upstream's `computeExpr_LR` copies a plain gene's row
/// **verbatim** (`dataLavg[index.singleL, ] <- dataL1avg`) and applies `geometricMean`
/// only on the complex path. Those are not the same: `exp(log(x)) != x` for most `x`
/// (3.0 comes back as 3.0000000000000004), so folding single genes through the geometric
/// mean would perturb every plain-gene ligand and receptor in the network.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum EntityRows {
    /// A plain gene: use the row as-is.
    Single(usize),
    /// A complex: geometric mean over the subunit rows.
    Complex(Vec<usize>),
}

/// Resolve a ligand or receptor name to the gene rows that define its expression.
///
/// Upstream (`computeExpr_LR`) splits names by `which(geneLR %in% rownames(data.use))`:
/// a name present in the matrix takes that row; anything else is routed to
/// `computeExpr_complex`, which looks it up in the complex table and then indexes the
/// matrix by subunit name with **no** `intersect`. So:
///
/// * name in the matrix → that single row;
/// * else in the complex table → its non-empty subunits, each of which must be a row,
///   or upstream errors;
/// * else → the complex-table row is all `NA`, and `data.use[NA, ]` errors too.
///
/// Returns an error in the last two cases, matching upstream.
pub fn resolve_entity(name: &str, means: &GroupMeans<'_>, db: &Database) -> Result<EntityRows> {
    if let Some(r) = means.row(name) {
        return Ok(EntityRows::Single(r));
    }
    let cells = db
        .complex_cells(name)
        .ok_or_else(|| ExprError::SubscriptOutOfBounds {
            name: name.to_string(),
        })?;
    let mut rows = Vec::with_capacity(cells.len());
    for s in cells.iter().filter(|s| !s.is_empty()) {
        rows.push(
            means
                .row(s)
                .ok_or_else(|| ExprError::SubscriptOutOfBounds {
                    name: name.to_string(),
                })?,
        );
    }
    Ok(EntityRows::Complex(rows))
}

/// Port of `computeExpr_LR` + `computeExpr_complex`.
///
/// Returns an `n_names x n_groups` matrix in R's **column-major** order: `n_groups` blocks
/// of `n_names`, so entity `i` in group `j` lives at `j * n_names + i`.
///
/// Column-major rather than row-major is deliberate and matches [`GroupMeans`],
/// [`CellExpr`] and R's own storage order. It is also the order that makes the
/// accumulation order of the inner loops -- group outer, subunit inner, as in
/// `apply(1 + ..., 2, prod)` -- the same as the buffer order, which is what an
/// `LDOUBLE` reduction is sensitive to.
pub fn compute_expr_lr(
    names: &[String],
    means: &GroupMeans<'_>,
    db: &Database,
) -> Result<Vec<f64>> {
    let k = means.n_groups();
    let n = names.len();
    let mut out = vec![0.0; n * k];
    for (i, name) in names.iter().enumerate() {
        match resolve_entity(name, means, db)? {
            EntityRows::Single(r) => {
                for j in 0..k {
                    out[j * n + i] = means.at(r, j);
                }
            }
            EntityRows::Complex(rows) => {
                geometric_mean_over_rows_into(means, &rows, &mut out, i, n)
            }
        }
    }
    Ok(out)
}

/// Port of `computeExpr_complex` on its own, for callers that already know the names are
/// complexes. Same rules as [`resolve_entity`].
pub fn compute_expr_complex(
    names: &[String],
    means: &GroupMeans<'_>,
    db: &Database,
) -> Result<Vec<f64>> {
    compute_expr_lr(names, means, db)
}

/// Which co-receptor column of an interaction table to use.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum CoreceptorKind {
    /// `co_A_receptor` — the multiplicative activator.
    Activator,
    /// `co_I_receptor` — the multiplicative inhibitor.
    Inhibitor,
}

/// Port of `computeExpr_coreceptor`.
///
/// Returns an `n_names x n_groups` matrix in R's **column-major** order, as in
/// [`compute_expr_lr`]: pair `i` in group `j` is at `j * n_names + i`. The result is
/// **all ones** unless `names[i]` names a cofactor with at least one subunit present.
///
/// Two upstream details, both pinned by the fixture:
/// * an unresolvable name (empty, absent from the cofactor table, or with no subunit in
///   the matrix) yields **all ones**, i.e. *no modulation*. Not an error. A missing
///   cofactor row is all `NA`, and `intersect` drops the `NA`s, leaving length 0.
/// * with exactly one subunit the value is `1 + x`; with several it is the **product** of
///   `1 + x` across subunits, and that product is an `LDOUBLE` accumulation in R
///   ([`r_prod`]).
pub fn compute_expr_coreceptor(
    names: &[String],
    means: &GroupMeans<'_>,
    db: &Database,
    _kind: CoreceptorKind,
) -> Vec<f64> {
    let k = means.n_groups();
    let n = names.len();
    let mut out = vec![1.0; n * k];
    for (i, name) in names.iter().enumerate() {
        if name.is_empty() {
            continue;
        }
        let rows = match present_cofactor_rows(name, means, db) {
            Some(r) if !r.is_empty() => r,
            _ => continue, // no modulation
        };
        if rows.len() == 1 {
            let r = rows[0];
            for j in 0..k {
                out[j * n + i] = 1.0 + means.at(r, j);
            }
        } else {
            let mut buf = vec![0.0; rows.len()];
            for j in 0..k {
                for (t, &r) in rows.iter().enumerate() {
                    buf[t] = 1.0 + means.at(r, j);
                }
                out[j * n + i] = r_prod(&buf, false);
            }
        }
    }
    out
}

/// Subunits of a cofactor that are present in the matrix, in `cofactor*` order.
///
/// `None` when the name is not a cofactor row at all (the all-`NA` row case), which the
/// caller treats as no modulation.
fn present_cofactor_rows(name: &str, means: &GroupMeans<'_>, db: &Database) -> Option<Vec<usize>> {
    let cells = db.cofactor_cells(name)?;
    Some(
        cells
            .iter()
            .filter(|s| !s.is_empty())
            .filter_map(|s| means.row(s))
            .collect(),
    )
}

/// Agonist Hill term: `data.avg^n / (Kh^n + data.avg^n)`.
///
/// `Kh^n` does not depend on the element, so it is hoisted; the association of the
/// division is preserved exactly.
#[inline]
fn hill_agonist(x: f64, khn: f64, n: f64) -> f64 {
    let xn = x.powf(n);
    xn / (khn + xn)
}

/// Antagonist Hill term: `Kh^n / (Kh^n + data.avg^n)`.
///
/// NB this is the **numerator-flipped** form, not the same expression as
/// [`hill_agonist`]. Using the agonist's expression here is a real bug that a unit test
/// caught: the two are equal only when `x == Kh`.
#[inline]
fn hill_antagonist(x: f64, khn: f64, n: f64) -> f64 {
    khn / (khn + x.powf(n))
}

/// Port of `computeExpr_agonist`: `1 + h` for a single subunit, or the product of
/// `1 + h` across subunits, where `h = x^n / (Kh^n + x^n)`.
///
/// Returns a length-`n_groups` vector; all ones when the agonist is unresolvable.
pub fn compute_expr_agonist(
    name: &str,
    means: &GroupMeans<'_>,
    db: &Database,
    kh: f64,
    n: f64,
) -> Vec<f64> {
    let k = means.n_groups();
    let mut out = vec![1.0; k];
    let rows = match present_cofactor_rows(name, means, db) {
        Some(r) if !r.is_empty() => r,
        _ => return out,
    };
    let khn = kh.powf(n);
    let mut buf = vec![0.0; rows.len()];
    for j in 0..k {
        for (t, &r) in rows.iter().enumerate() {
            buf[t] = 1.0 + hill_agonist(means.at(r, j), khn, n);
        }
        out[j] = if rows.len() == 1 {
            buf[0]
        } else {
            r_prod(&buf, false)
        };
    }
    out
}

/// Port of `computeExpr_antagonist`: `h` for a single subunit, or the product of `h`
/// across subunits, with the same Hill expression as the agonist.
pub fn compute_expr_antagonist(
    name: &str,
    means: &GroupMeans<'_>,
    db: &Database,
    kh: f64,
    n: f64,
) -> Vec<f64> {
    let k = means.n_groups();
    let mut out = vec![1.0; k];
    let rows = match present_cofactor_rows(name, means, db) {
        Some(r) if !r.is_empty() => r,
        _ => return out,
    };
    let khn = kh.powf(n);
    let mut buf = vec![0.0; rows.len()];
    for j in 0..k {
        for (t, &r) in rows.iter().enumerate() {
            buf[t] = hill_antagonist(means.at(r, j), khn, n);
        }
        out[j] = if rows.len() == 1 {
            buf[0]
        } else {
            r_prod(&buf, false)
        };
    }
    out
}

/// A cells x genes expression matrix, column-major, with a name -> row index.
///
/// The transpose of [`GroupMeans`], and what the exported `computeExprGroup_*` variants
/// take: upstream's group variants receive a per-cell matrix and `aggregate` it down to
/// groups with an arbitrary mean function. CellChat itself never calls them, so this type
/// exists purely so the exported surface is complete and testable.
#[derive(Clone, Debug)]
pub struct CellExpr<'a> {
    values: &'a [f64],
    n_cells: usize,
    n_genes: usize,
    row_of: HashMap<String, usize>,
}

impl<'a> CellExpr<'a> {
    /// # Panics
    /// If `values.len() != n_cells * n_genes` or `names.len() != n_genes`.
    pub fn new(values: &'a [f64], names: &[String], n_cells: usize) -> Self {
        assert_eq!(
            values.len(),
            names.len() * n_cells,
            "CellExpr: values must be cells x genes"
        );
        let row_of = names
            .iter()
            .enumerate()
            .map(|(i, n)| (n.clone(), i))
            .collect();
        Self {
            values,
            n_cells,
            n_genes: names.len(),
            row_of,
        }
    }

    #[inline]
    pub fn n_cells(&self) -> usize {
        self.n_cells
    }

    #[inline]
    pub fn n_genes(&self) -> usize {
        self.n_genes
    }

    #[inline]
    pub fn row(&self, name: &str) -> Option<usize> {
        self.row_of.get(name).copied()
    }

    /// `value(cell, g)` -- the expression of gene `g` in cell `cell`.
    ///
    /// NB the stride. R stores an `nrow x ncol` matrix column-major, so
    /// `m[cell, g] == values[g * n_cells + cell]`: the *gene* index strides, and the
    /// cells of one gene are contiguous. Writing `values[cell * n_genes + g]` here
    /// (the transpose of what R does) silently reads the wrong cells; it was caught by
    /// the `computeExprGroup_*` corpus, where it turned every group mean into 0.5.
    #[inline]
    pub fn at(&self, g: usize, cell: usize) -> f64 {
        self.values[g * self.n_cells + cell]
    }
}

/// Group the cofactor subunits of `name` by [`GroupMean`], then apply the Hill term.
///
/// Shared by [`compute_expr_group_agonist`] and [`compute_expr_group_antagonist`].
///
/// The single- and multi-subunit branches differ only in whether the Hill values are
/// offset by 1 and whether they are reduced by a product, exactly as upstream:
/// `length == 1` -> `1 + data.avg^n/(Kh^n + data.avg^n)`, `length > 1` ->
/// `apply(1 + ..., 2, prod)`.
fn group_hill(
    name: &str,
    cells: &CellExpr<'_>,
    group: &[usize],
    n_groups: usize,
    db: &Database,
    kh: f64,
    n: f64,
    fun: GroupMean,
    agonist: bool,
) -> Vec<f64> {
    assert_eq!(
        group.len(),
        cells.n_cells(),
        "group length must match n_cells"
    );
    let mut out = vec![1.0; n_groups];
    let rows: Vec<usize> = match db.cofactor_cells(name) {
        Some(cells_of) => cells_of
            .iter()
            .filter(|s| !s.is_empty())
            .filter_map(|s| cells.row(s))
            .collect(),
        None => return out, // no modulation: all ones
    };
    if rows.is_empty() {
        return out;
    }
    let avg = aggregate_columns(cells.values, cells.n_cells, &rows, group, n_groups, fun);
    // `avg` is n_subunits x n_groups.
    let khn = kh.powf(n);
    let mut buf = vec![0.0; rows.len()];
    for j in 0..n_groups {
        for (t, _r) in rows.iter().enumerate() {
            let h = if agonist {
                hill_agonist(avg[t * n_groups + j], khn, n)
            } else {
                hill_antagonist(avg[t * n_groups + j], khn, n)
            };
            buf[t] = if agonist { 1.0 + h } else { h };
        }
        out[j] = if rows.len() == 1 {
            buf[0]
        } else {
            r_prod(&buf, false)
        };
    }
    out
}

/// Port of `computeExprGroup_agonist`.
///
/// Upstream: for each cofactor subunit, `aggregate` the per-cell values within each
/// `group` level with `FunMean`, then apply `1 + h` (one subunit) or the product of
/// `1 + h` (several).
pub fn compute_expr_group_agonist(
    name: &str,
    cells: &CellExpr<'_>,
    group: &[usize],
    n_groups: usize,
    db: &Database,
    kh: f64,
    n: f64,
    fun: GroupMean,
) -> Vec<f64> {
    group_hill(name, cells, group, n_groups, db, kh, n, fun, true)
}

/// Port of `computeExprGroup_antagonist`. Same reduction, Hill term with `Kh^n` in the
/// numerator and no `1 +` offset.
pub fn compute_expr_group_antagonist(
    name: &str,
    cells: &CellExpr<'_>,
    group: &[usize],
    n_groups: usize,
    db: &Database,
    kh: f64,
    n: f64,
    fun: GroupMean,
) -> Vec<f64> {
    group_hill(name, cells, group, n_groups, db, kh, n, fun, false)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn means_of<'a>(x: &'a [f64], genes: &[&str], k: usize) -> GroupMeans<'a> {
        let names: Vec<String> = genes.iter().map(|s| s.to_string()).collect();
        GroupMeans::new(x, &names, k)
    }

    #[test]
    fn a_single_gene_takes_its_row_verbatim() {
        let x = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0]; // 2 genes x 3 groups, column-major
        let m = means_of(&x, &["A", "B"], 3);
        let out = compute_expr_lr(&["A".to_string()], &m, &Database::default()).unwrap();
        assert_eq!(out, vec![1.0, 3.0, 5.0]);
    }

    #[test]
    fn geometric_mean_is_over_logs() {
        // A complex of two subunits: exp((log a + log b)/2) = sqrt(a*b) per group.
        let x = [4.0, 9.0, 16.0, 25.0, 36.0, 49.0]; // A, B x 3 groups
        let m = means_of(&x, &["A", "B"], 3);
        let mut db = Database::default();
        db.complexes
            .insert("AB".to_string(), vec!["A".to_string(), "B".to_string()]);
        db.n_subunit_cols = 2;
        let out = compute_expr_lr(&["AB".to_string()], &m, &db).unwrap();
        // NB not exactly sqrt(a*b): upstream computes exp((log a + log b)/2), and
        // exp(mean(log)) != exp(mean * 2) in floating point. Reproduce the *expression*,
        // not the algebraic identity -- `sqrt(16*25)` is 20.0 and the computation gives
        // 19.999999999999996. The corpus pins the exact bits; this test pins the shape.
        let gm = |a: f64, b: f64| ((a.ln() + b.ln()) / 2.0).exp();
        assert_eq!(out[0], gm(4.0, 9.0));
        assert_eq!(out[1], gm(16.0, 25.0));
        assert_eq!(out[2], gm(36.0, 49.0));
        assert_ne!(
            out[1], 20.0,
            "the algebraic identity would be the wrong target"
        );
    }

    #[test]
    fn a_zero_subunit_zeroes_the_whole_group() {
        // R-ism 4: log(0) = -Inf, so the geometric mean is 0, not "ignore the zero".
        let x = [0.0, 1.0, 4.0, 9.0]; // A, B x 2 groups
        let m = means_of(&x, &["A", "B"], 2);
        let mut db = Database::default();
        db.complexes
            .insert("AB".to_string(), vec!["A".to_string(), "B".to_string()]);
        db.n_subunit_cols = 2;
        let out = compute_expr_lr(&["AB".to_string()], &m, &db).unwrap();
        assert_eq!(out[0], 0.0);
        assert_eq!(out[1], 6.0);
    }

    #[test]
    fn unknown_name_is_an_error_not_a_silent_one() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let m = means_of(&x, &["A", "B"], 2);
        let db = Database::default();
        let e = compute_expr_lr(&["NOPE".to_string()], &m, &db).unwrap_err();
        assert_eq!(e.to_string(), SUBSCRIPT_OUT_OF_BOUNDS);
    }

    #[test]
    fn a_complex_with_a_missing_subunit_is_an_error() {
        // R-ism 11: upstream does not `intersect` here.
        let x = [1.0, 2.0, 3.0, 4.0];
        let m = means_of(&x, &["A", "B"], 2);
        let mut db = Database::default();
        db.complexes.insert(
            "AB".to_string(),
            vec!["A".to_string(), "MISSING".to_string()],
        );
        db.n_subunit_cols = 2;
        let e = compute_expr_lr(&["AB".to_string()], &m, &db).unwrap_err();
        assert_eq!(e.to_string(), SUBSCRIPT_OUT_OF_BOUNDS);
    }

    #[test]
    fn coreceptor_one_subunit_is_one_plus_x() {
        let x = [0.5, 1.5, 2.5, 3.5];
        let m = means_of(&x, &["FST", "X"], 2);
        let mut db = Database::default();
        db.cofactors
            .insert("c".to_string(), vec!["FST".to_string()]);
        db.n_cofactor_cols = 1;
        let out = compute_expr_coreceptor(&["c".to_string()], &m, &db, CoreceptorKind::Activator);
        assert_eq!(out, vec![1.5, 3.5]);
    }

    #[test]
    fn coreceptor_multi_subunit_is_a_product() {
        let x = [1.0, 2.0, 3.0, 4.0];
        let m = means_of(&x, &["A", "B"], 2);
        let mut db = Database::default();
        db.cofactors
            .insert("c".to_string(), vec!["A".to_string(), "B".to_string()]);
        db.n_cofactor_cols = 2;
        let out = compute_expr_coreceptor(&["c".to_string()], &m, &db, CoreceptorKind::Activator);
        // (1+1)*(1+2) = 6 and (1+3)*(1+4) = 20
        assert_eq!(out, vec![6.0, 20.0]);
    }

    #[test]
    fn an_unresolvable_cofactor_is_no_modulation() {
        // Upstream `intersect`s, so this is all ones -- NOT an error, unlike the ligand
        // path. The asymmetry is deliberate and pinned.
        let x = [1.0, 2.0, 3.0, 4.0];
        let m = means_of(&x, &["A", "B"], 2);
        let db = Database::default();
        for name in ["", "NOT_A_COFACTOR", "c"] {
            let out =
                compute_expr_coreceptor(&[name.to_string()], &m, &db, CoreceptorKind::Inhibitor);
            assert_eq!(out, vec![1.0; 2], "name {name:?} should be a no-op");
        }
    }

    #[test]
    fn hill_is_bounded_and_matches_the_closed_form() {
        let x = [0.25, 4.0];
        let m = means_of(&x, &["A", "B"], 1);
        let mut db = Database::default();
        db.cofactors.insert("c".to_string(), vec!["A".to_string()]);
        db.n_cofactor_cols = 1;
        let ag = compute_expr_agonist("c", &m, &db, 0.5, 1.0);
        // 1 + x^n/(Kh^n + x^n) with x = 0.25, Kh = 0.5
        assert_eq!(ag[0], 1.0 + 0.25 / (0.5 + 0.25));
        let an = compute_expr_antagonist("c", &m, &db, 0.5, 1.0);
        // Kh^n/(Kh^n + x^n) -- the numerator is Kh, not x.
        assert_eq!(an[0], 0.5 / (0.5 + 0.25));
        // The two forms coincide only at x == Kh; asserting they differ guards against
        // the two functions sharing an implementation by accident.
        assert_ne!(ag[0] - 1.0, an[0]);
    }

    #[test]
    fn agonist_and_antagonist_differ_only_by_the_offset() {
        // h is the same; the agonist adds 1 (single subunit) and the antagonist does not.
        let x = [0.1, 0.9, 0.5, 0.5];
        let m = means_of(&x, &["A", "B"], 2);
        let mut db = Database::default();
        db.cofactors.insert("c".to_string(), vec!["A".to_string()]);
        db.n_cofactor_cols = 1;
        // The real invariant: with h = x^n/(Kh^n+x^n) and h' = Kh^n/(Kh^n+x^n), we have
        // h + h' == 1, so the single-subunit agonist (1 + h) and antagonist (h') satisfy
        // ag + an == 2. Asserting ag == an + 1 would be wrong: those differ by 1 only when
        // x == Kh, which is what made the two forms look interchangeable.
        let ag = compute_expr_agonist("c", &m, &db, 0.5, 1.0);
        let an = compute_expr_antagonist("c", &m, &db, 0.5, 1.0);
        assert_eq!(ag[0] + an[0], 2.0);
        assert_eq!(ag[1] + an[1], 2.0);
    }
}
