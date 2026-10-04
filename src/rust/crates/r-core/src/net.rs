//! `aggregateNet`, `subsetCommunication` and the array-to-table melt they share.
//!
//! Parity contract: **Exact** (`docs/SEMANTICS.md`).
//!
//! Upstream: `R/modeling.R:449` (`aggregateNet`), `R/analysis.R:1944`
//! (`subsetCommunication`) and `R/analysis.R:1900` (`subsetCommunication_internal`).
//!
//! ## Why this module is fussier than it looks
//!
//! Both functions are a handful of arithmetic operations wrapped in data-frame plumbing,
//! and all of the difficulty is in the plumbing:
//!
//! * **The melt order is not the array's own order.** `reshape2::melt` on a
//!   `k x k x nLR` array builds its label frame with `expand.grid(dimnames(x))`, and
//!   `expand.grid` varies its **first** argument fastest. So the melted rows run
//!   interaction-major, then target, then source -- the opposite of the column-major
//!   traversal R uses to store the array. The result is a table whose row order is
//!   `(interaction, target, source)`, and that order is the return value.
//! * **`sum()` is not associative, so the accumulation order is part of the result.**
//!   `aggregateNet`'s second branch groups by the *character* key `source|target`
//!   (dplyr sorts it), and `tapply(..., sum)` then adds in that grouped order. Getting the
//!   sort or the within-cell order wrong changes the last ulp of `net$weight`.
//! * **`net$weight` and `net$count` are NA-filled to 0, but only after `tapply`**, and
//!   `tapply` leaves `NA` wherever a cell has no rows -- which is different from summing
//!   zeroes, and different again from leaving the cell absent.
//! * **`prob[pval >= thresh] <- 0` is a `>=` against `thresh`,** so a pair with
//!   `pval == thresh` is dropped. `filterCommunication` and `aggregateNet` both rely on
//!   it, and the boundary is reachable with `thresh = 1/nboot`.
//!
//! ## `Prob` can exceed 1 (R-ism 21)
//!
//! Neither function clamps. The co-agonist product makes `Prob > 1` legitimate, and
//! `net$weight` inherits it. Any "sanitising" here would change published CellChat results.

use std::collections::HashMap;

use crate::longdouble::F80;

/// One `k x k x nLR` network: `Prob` and `Pval`, column-major like R's `array(0, dim =)`.
#[derive(Clone, Debug, PartialEq)]
pub struct Net {
    pub prob: Vec<f64>,
    pub pval: Vec<f64>,
    /// `k`.
    pub n_groups: usize,
    /// `levels(object@idents)`.
    pub group_levels: Vec<String>,
    /// `dimnames(prob)[[3]]`, i.e. `rownames(object@LR$LRsig)`.
    pub interaction_names: Vec<String>,
}

impl Net {
    /// # Panics
    /// If the two buffers are not `n_groups^2 * n_lr`, or the level/name vectors disagree
    /// in length. Both are internal wiring errors.
    pub fn new(
        prob: Vec<f64>,
        pval: Vec<f64>,
        group_levels: Vec<String>,
        interaction_names: Vec<String>,
    ) -> Self {
        let k = group_levels.len();
        let n = k * k * interaction_names.len();
        assert_eq!(prob.len(), n, "Net: prob buffer must be k*k*nLR");
        assert_eq!(pval.len(), n, "Net: pval buffer must be k*k*nLR");
        assert_eq!(
            group_levels.len(),
            n_groups_of(&prob, &pval, &interaction_names),
            "Net: group_levels is self-consistent"
        );
        Self {
            prob,
            pval,
            n_groups: k,
            group_levels,
            interaction_names,
        }
    }

    #[inline]
    fn n_lr(&self) -> usize {
        self.interaction_names.len()
    }

    /// `prob[a, b, c]`, R's 1-based indices.
    #[inline]
    pub fn prob_at(&self, a: usize, b: usize, c: usize) -> f64 {
        self.prob[self.idx(a, b, c)]
    }

    /// `pval[a, b, c]`.
    #[inline]
    pub fn pval_at(&self, a: usize, b: usize, c: usize) -> f64 {
        self.pval[self.idx(a, b, c)]
    }

    #[inline]
    fn idx(&self, a: usize, b: usize, c: usize) -> usize {
        c * self.n_groups * self.n_groups + b * self.n_groups + a
    }

    /// `dimnames(prob)`.
    pub fn dimnames(&self) -> (Vec<String>, Vec<String>, Vec<String>) {
        (
            self.group_levels.clone(),
            self.group_levels.clone(),
            self.interaction_names.clone(),
        )
    }
}

fn n_groups_of(prob: &[f64], _pval: &[f64], names: &[String]) -> usize {
    let n = prob.len() / names.len().max(1);
    // Round-trip of k from k*k = n; the assert above catches a genuine mismatch, and this
    // keeps the constructor total for the empty case.
    let mut k = (n as f64).sqrt().round() as usize;
    while k * k < n {
        k += 1;
    }
    k
}

/// The columns `subsetCommunication` returns for `slot.name = "net"`, in upstream's order.
///
/// `intersect(c(...), colnames(net))` over the DEG-augmented column list when
/// `identifyOverExpressedGenes` has been run, and over the shorter list when it has not.
/// Both orders are reproduced; `has_deg` selects.
pub const COLUMNS_PLAIN: [&str; 11] = [
    "source",
    "target",
    "ligand",
    "receptor",
    "prob",
    "pval",
    "interaction_name",
    "interaction_name_2",
    "pathway_name",
    "annotation",
    "evidence",
];

/// Every column the melt and the L-R join can produce, in build order, for
/// `final_select = false`. `source, target, interaction_name, prob, pval` come from the melt
/// and are in *melt* order, not in the order any `col.use` lists them; the rest are attached
/// in `LR`'s column order.
pub const ALL_COLUMNS: [&str; 12] = [
    "source",
    "target",
    "interaction_name",
    "prob",
    "pval",
    "interaction_name_2",
    "pathway_name",
    "ligand",
    "receptor",
    "annotation",
    "evidence",
    "datasets",
];

/// The `intersect()` of the longer list with the columns that exist: identical to
/// [`COLUMNS_PLAIN`] but with the DEG block in the middle, so the *order* differs.
pub const COLUMNS_DEG: [&str; 20] = [
    "source",
    "target",
    "ligand",
    "receptor",
    "prob",
    "pval",
    "interaction_name",
    "interaction_name_2",
    "pathway_name",
    "annotation",
    "evidence",
    "datasets",
    "ligand.logFC",
    "ligand.pct.1",
    "ligand.pct.2",
    "ligand.pvalues",
    "receptor.logFC",
    "receptor.pct.1",
    "receptor.pct.2",
    "receptor.pvalues",
];

/// The optional columns `LR` can supply, in upstream's `col.use` order.
pub const LR_OPTIONAL: [&str; 6] = [
    "interaction_name_2",
    "pathway_name",
    "ligand",
    "receptor",
    "annotation",
    "evidence",
];

/// One row of the L-R table, the columns `subsetCommunication` attaches.
///
/// ## `columns` is the point
///
/// `subsetCommunication` picks output columns with
/// `intersect(c("source", "target", "ligand", ...), colnames(net))`, and `colnames(net)`
/// lists what the **table** has, not what its cells are non-`NA`. So `evidence` appears in
/// the output with `NA` in every row, and it is *absent* only if the user's `LRsig` has no
/// such column at all. Testing the values instead of the column set drops `evidence` from
/// the result, which silently changes the output's column count.
///
/// The same distinction is what makes upstream's
/// `net[rowSums(is.na(net)) != ncol(net), ]` necessary: a row with `evidence = NA` must
/// survive, which it would not if a single `NA` dropped it.
#[derive(Clone, Debug, PartialEq)]
pub struct LrMeta {
    /// `rownames(LR)`, matched against `dimnames(prob)[[3]]`.
    pub interaction_name: String,
    /// Which of [`LR_OPTIONAL`] this table actually has. Defaults to all of them, which is
    /// what a `CellChatDB`-derived `LRsig` has.
    pub columns: Vec<String>,
    pub interaction_name_2: String,
    pub pathway_name: String,
    pub ligand: String,
    pub receptor: String,
    pub annotation: String,
    /// `NA` unless the user filled it in; the column still exists.
    pub evidence: Option<String>,
}

impl Default for LrMeta {
    fn default() -> Self {
        Self {
            interaction_name: String::new(),
            columns: LR_OPTIONAL.iter().map(|s| (*s).to_string()).collect(),
            interaction_name_2: String::new(),
            pathway_name: String::new(),
            ligand: String::new(),
            receptor: String::new(),
            annotation: String::new(),
            evidence: None,
        }
    }
}

impl LrMeta {
    /// The upstream default: `evidence` is `NA` but the column exists.
    pub fn minimal(
        interaction_name: &str,
        interaction_name_2: &str,
        ligand: &str,
        receptor: &str,
        annotation: &str,
    ) -> Self {
        Self {
            interaction_name: interaction_name.into(),
            interaction_name_2: interaction_name_2.into(),
            ligand: ligand.into(),
            receptor: receptor.into(),
            annotation: annotation.into(),
            ..Default::default()
        }
    }
}

/// The `data.frame` `subsetCommunication` returns.
#[derive(Clone, Debug, Default, PartialEq)]
pub struct NetTable {
    /// Column names, in the order `subsetCommunication` selects them.
    pub columns: Vec<String>,
    /// One entry per cell; `None` is `NA`. Indexed as `rows[r][columns.index_of(name)]`.
    pub rows: Vec<Vec<Option<String>>>,
    /// The source/target factor levels, which the melted order follows.
    pub group_levels: Vec<String>,
    /// The `interaction_name` factor's levels: the `Prob` array's whole third dimnames, in order.
    ///
    /// Not the subset of interaction names that survived the threshold. Upstream melts the array and
    /// `var.convert`s the result, which factors every column that came from a dimname, so
    /// `interaction_name` is a factor over the *full* L-R set and an interaction filtered out of the
    /// table still contributes an unused level. R's `subsetCommunication` needs these levels to
    /// build the factor; carrying them here means the R-isms stay in one place and the level
    /// semantics are testable from Rust.
    pub interaction_levels: Vec<String>,
}

impl NetTable {
    #[inline]
    fn col(&self, name: &str) -> usize {
        self.columns
            .iter()
            .position(|c| c == name)
            .unwrap_or_else(|| panic!("NetTable has no column {name:?}"))
    }

    /// `nrow()`.
    pub fn len(&self) -> usize {
        self.rows.len()
    }

    pub fn is_empty(&self) -> bool {
        self.rows.is_empty()
    }

    /// Cell `(r, name)`, or `None` for `NA`.
    pub fn get(&self, r: usize, name: &str) -> Option<&str> {
        self.rows[r][self.col(name)].as_deref()
    }

    /// `f64` cell `(r, name)`.
    pub fn num(&self, r: usize, name: &str) -> Option<f64> {
        self.get(r, name).and_then(|v| v.parse().ok())
    }
}

/// `subsetCommunication(object)` for `slot.name = "net"`, `mode = "single"`.
///
/// `thresh` is the `pval` cut: `prob[pval >= thresh] <- 0`, then rows with `prob == 0` are
/// dropped. `sources_use` / `targets_use` are applied last, in that order, after the
/// threshold, exactly as upstream does.
pub fn subset_communication(
    net: &Net,
    lr: &[LrMeta],
    thresh: f64,
    sources_use: Option<&[String]>,
    targets_use: Option<&[String]>,
    has_deg: bool,
    final_select: bool,
) -> NetTable {
    let k = net.n_groups;
    let n_lr = net.n_lr();
    let by_name: HashMap<&str, usize> = lr
        .iter()
        .enumerate()
        .map(|(i, m)| (m.interaction_name.as_str(), i))
        .collect();

    // `reshape2::melt` order: interaction-major, then target, then source. See the module
    // note -- this is NOT the array's own order and it is the return value's row order.
    // `flat[r]` is the array index of surviving row `r`. It cannot be recovered from `r`:
    // the melt emits a *contiguous* enumeration of the array and then `subset(prob > 0)`
    // drops rows, so the row ordinal and the array index diverge as soon as anything is
    // filtered. Carrying it explicitly is the only way to keep the source/target/interaction
    // labels attached to the right cell.
    let mut flat: Vec<usize> = Vec::new();
    let mut prob_col: Vec<f64> = Vec::new();
    let mut pval_col: Vec<f64> = Vec::new();
    let mut meta_col: Vec<Option<usize>> = Vec::new();
    for c in 0..n_lr {
        for b in 0..k {
            for a in 0..k {
                let p = net.prob_at(a, b, c);
                let pv = net.pval_at(a, b, c);
                // `prob[pval >= thresh] <- 0` happens *before* the melt, so it applies to
                // every cell, and then `subset(net, prob > 0)` drops what became zero.
                if !(pv < thresh) {
                    continue;
                }
                if !(p > 0.0) {
                    continue;
                }
                flat.push(net.idx(a, b, c));
                prob_col.push(p);
                pval_col.push(pv);
                // `idx <- match(net$interaction_name, rownames(pairLR))`; a name absent
                // from `LR` gives `NA` for every attached column, which is what upstream's
                // `rowSums(is.na(net)) != ncol(net)` guard exists to tolerate.
                meta_col.push(by_name.get(net.interaction_names[c].as_str()).copied());
            }
        }
    }
    let n_rows = flat.len();

    // `rowSums(is.na(net)) != ncol(net)`: drop rows that are NA in *every* column. With the
    // plain column set that is impossible (source/target/prob/pval are never NA), but a
    // row whose `interaction_name` is absent from `LR` has NA ligand/receptor/annotation,
    // which is what this guard is for upstream. Reproduced: a row is dropped when every
    // metadata column is NA.
    let keep: Vec<usize> = (0..n_rows).filter(|&r| meta_col[r].is_some()).collect();

    // sources/targets filters, applied after the threshold and in that order. Numeric
    // `sources.use` has already been resolved to names by the caller, exactly as upstream
    // does with `cells.level[sources.use]`.
    let wanted = |u: Option<&[String]>, level: usize| match u {
        None => true,
        Some(names) => names
            .iter()
            .any(|n| net.group_levels.iter().position(|l| l == n) == Some(level)),
    };
    let keep: Vec<usize> = keep
        .into_iter()
        .filter(|&r| {
            wanted(sources_use, src_idx(net, flat[r])) && wanted(targets_use, tgt_idx(net, flat[r]))
        })
        .collect();

    // `final_select = false` returns the melted-and-joined table *unfiltered by column*, so
    // that [`crate::subset::subset_communication_deg`] can run the DEG thresholds and the
    // final `col.use` itself. It has to be the unselected table, not the `slot = "net"` one:
    // upstream's `rowSums(is.na(net)) != ncol(net)` is evaluated on the full column set, so
    // selecting first would change which rows count as all-`NA`.
    let (want, final_select): (&[&str], bool) = match (has_deg, final_select) {
        (true, true) => (&COLUMNS_DEG, true),
        (false, true) => (&COLUMNS_PLAIN, true),
        // Every column the table actually has, in the order the melt and the join built them.
        (_, false) => (&ALL_COLUMNS[..], false),
    };
    // `intersect(want, colnames(net))`: `net` carries the melt columns
    // (source, target, interaction_name, prob, pval) plus whatever `LR` supplied among
    // {interaction_name_2, pathway_name, ligand, receptor, annotation, evidence}, and the
    // DEG block only when it exists.
    let mut columns: Vec<String> = Vec::new();
    for c in want {
        if melt_has(c)
            || meta_has(c, lr)
            || (!final_select && deg_has(c))
            || (has_deg && deg_has(c))
        {
            columns.push((*c).to_string());
        }
    }
    let mut table = NetTable {
        columns,
        rows: Vec::with_capacity(keep.len()),
        group_levels: net.group_levels.clone(),
        interaction_levels: net.interaction_names.clone(),
    };
    for &r in &keep {
        let mut row: Vec<Option<String>> = Vec::with_capacity(table.columns.len());
        for name in table.columns.clone() {
            row.push(match name.as_str() {
                "source" => Some(net.group_levels[src_idx(net, flat[r])].clone()),
                "target" => Some(net.group_levels[tgt_idx(net, flat[r])].clone()),
                "prob" => Some(fmt_f64(prob_col[r])),
                "pval" => Some(fmt_f64(pval_col[r])),
                "interaction_name" => Some(net.interaction_names[c_idx(net, flat[r])].clone()),
                _ => match meta_col[r].and_then(|i| lr.get(i)) {
                    Some(m) => lr_field(m, &name),
                    None => None,
                },
                // `evidence` is NA by default, which is why the DEG guard exists upstream.
            });
        }
        table.rows.push(row);
    }
    table
}

fn melt_has(c: &str) -> bool {
    matches!(
        c,
        "source" | "target" | "prob" | "pval" | "interaction_name"
    )
}

fn deg_has(c: &str) -> bool {
    matches!(
        c,
        "datasets"
            | "ligand.logFC"
            | "ligand.pct.1"
            | "ligand.pct.2"
            | "ligand.pvalues"
            | "receptor.logFC"
            | "receptor.pct.1"
            | "receptor.pct.2"
            | "receptor.pvalues"
    )
}

/// Does the attached L-R table have column `c`? Presence, not non-`NA`.
///
/// Deliberately independent of which rows survived the threshold. `colnames(net)` is fixed
/// by the *structure* -- the melt columns plus `intersect(col.use, colnames(LR))` -- so a
/// query that returns zero rows still has all eleven columns. Deriving the columns from
/// the surviving rows instead collapses them to the five melt columns, which is a
/// different-shaped zero-row answer.
fn meta_has(c: &str, lr: &[LrMeta]) -> bool {
    lr.iter().any(|m| m.columns.iter().any(|x| x == c))
}

fn lr_field(m: &LrMeta, name: &str) -> Option<String> {
    Some(match name {
        "ligand" => m.ligand.clone(),
        "receptor" => m.receptor.clone(),
        "annotation" => m.annotation.clone(),
        "interaction_name_2" => m.interaction_name_2.clone(),
        "pathway_name" => m.pathway_name.clone(),
        "evidence" => m.evidence.clone()?,
        _ if deg_has(name) => return None,
        _ => return None,
    })
}

fn fmt_f64(v: f64) -> String {
    // 17 significant digits round-trips every finite double, which is what R's
    // `as.character` does for the columns that come from a numeric.
    format!("{v:.17e}")
}

// The melt row index encodes (interaction, target, source) as
// `c * k^2 + b * k + a`, which is also the array's own column-major index. That is not a
// coincidence -- `as.vector(data)` is column-major, and the labels are permuted to match --
// so the two agree and the helpers below are just a restatement.
fn src_idx(net: &Net, r: usize) -> usize {
    r % net.n_groups
}
fn tgt_idx(net: &Net, r: usize) -> usize {
    (r / net.n_groups) % net.n_groups
}
fn c_idx(net: &Net, r: usize) -> usize {
    r / (net.n_groups * net.n_groups)
}
/// `aggregateNet`'s default branch, with no filters.
///
/// ```r
/// prob <- net$prob; pval <- net$pval
/// pval[prob == 0] <- 1
/// prob[pval >= thresh] <- 0
/// net$count  <- apply(prob > 0, c(1,2), sum)
/// net$weight <- apply(prob, c(1,2), sum)
/// ```
///
/// Two details that are easy to miss:
///
/// * `pval[prob == 0] <- 1` runs **first**, so a zero-probability cell is rescued to
///   `pval = 1` and then dropped by `pval >= thresh` — the same outcome as leaving it, but
///   only because `thresh <= 1`. The ordering is reproduced anyway.
/// * `apply(prob, c(1,2), sum)` is an R `sum()`, i.e. an **LDOUBLE** accumulation in
///   column-major order over the `k x k` slice. With `k` in the tens this is exactly the
///   case where an `f64` accumulator drifts, so [`F80`] is used.
///
/// # Panics
/// If the buffers are not `k^2 * n_lr`.
pub fn aggregate_net_default(net: &Net, thresh: f64) -> (Vec<f64>, Vec<f64>) {
    let k = net.n_groups;
    let n_lr = net.n_lr();
    let kk = k * k;
    let mut count = vec![0.0f64; kk];
    // The accumulators stay in F80 for the whole loop and are narrowed once at the end.
    // Narrowing after each add -- which a `weight[pos] = (weight[pos] + p)` shape does
    // implicitly -- makes it an f64 sum with extra steps, and that is 1 ulp off on
    // `apply(prob, c(1,2), sum)` at 8 terms, which is the *only* place the difference
    // shows up: at 3 or 4 terms both give the same answer.
    let mut acc: Vec<F80> = vec![F80::from_f64(0.0); kk];
    for c in 0..n_lr {
        for b in 0..k {
            for a in 0..k {
                let t = c * kk + b * k + a;
                let mut p = net.prob[t];
                let mut pv = net.pval[t];
                if p == 0.0 {
                    pv = 1.0;
                }
                if !(pv < thresh) {
                    p = 0.0;
                }
                let pos = b * k + a;
                if p > 0.0 {
                    count[pos] += 1.0;
                }
                if p != 0.0 {
                    acc[pos] = acc[pos].add(F80::from_f64(p));
                }
            }
        }
    }
    let mut weight: Vec<f64> = acc.iter().map(|a| a.to_f64()).collect();
    // `net$weight[is.na(...)] <- 0` and likewise for count. Neither can be NaN from an
    // f64 sum of a dense array, but the NA path exists upstream because `apply` on an
    // all-`NA` slice yields NA, and `prob` can be NA if `P.spatial` was. Preserved as a
    // no-op rather than dropped, with a comment, so the two are not confused later.
    for v in weight.iter_mut().chain(count.iter_mut()) {
        if v.is_nan() {
            *v = 0.0;
        }
    }
    (count, weight)
}

/// `aggregateNet`'s **filtered** branch, given the `df.net` that `subsetCommunication`
/// already produced.
///
/// Upstream is
/// ```r
/// df.net$source_target <- paste(df.net$source, df.net$target, sep = "|")
/// df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = "drop")
/// df.net3 <- df.net %>% group_by(source_target) %>% summarize(prob = sum(prob), .groups = "drop")
/// df.net2$prob <- df.net3$prob
/// a <- stringr::str_split(df.net2$source_target, "|", simplify = TRUE)
/// df.net2$source <- as.character(a[, 1])
/// df.net2$target <- as.character(a[, 2])
/// count <- tapply(df.net2[["count"]], list(df.net2$source], df.net2$target]), sum)
/// prob  <- tapply(df.net2[["prob"]],  list(df.net2$source], df.net2$target]), sum)
/// net$weight[is.na(net$weight)] <- 0
/// net$count[is.na(net$count)] <- 0
/// ```
///
/// # The filtered branch is broken upstream, and this reproduces that
///
/// `stringr::str_split(x, "|", simplify = TRUE)` takes a **regular expression**, and `|` is
/// the alternation metacharacter -- so it matches the empty string at every position rather
/// than the literal pipe:
///
/// ```r
/// stringr::str_split(c("g1|g2", "g10|g3"), "|", simplify = TRUE)
/// #      "" "g" "1" "|" "g" "2" "" ""
/// #      "" "g" "1" "0" "|" "g" "3" "" ""
/// ```
///
/// The result has `nchar(key) + 1` columns -- `""`, one per **character**, `""` -- so
/// `a[, 1]` is the empty string for every row and `a[, 2]` is a *single character*, not a
/// group name. Upstream then does
/// `df.net2$source <- factor(a[, 1], levels = cells.level[cells.level %in% unique(...)])`,
/// and no cell group is named `""`, so:
///
/// * `remove.isolate = TRUE` restricts the levels to `character(0)` and `tapply` returns a
///   **0 x 0** matrix;
/// * `remove.isolate = FALSE` keeps `cells.level`, every entry is `NA`, and the two
///   `is.na` fills turn it into a **k x k matrix of zeros**.
///
/// Either way `net$count` and `net$weight` carry no information: the aggregation upstream
/// intended never reaches the result. A bit-identical port has to return the same empty or
/// all-zero matrices, so this does. The `df.net2` grouping, which *is* correct, is still
/// computed and still ordered by the byte-wise sort of the joined key -- the interesting and
/// non-obvious part of the branch, and the part that would matter if upstream ever fixed the
/// regex.
///
/// Pinned by `src/rust/crates/r-core/tests/netfiltered_parity.rs` and re-checked through the
/// installed shim by `tests/parity/check_identical.R`.
pub fn aggregate_net_filtered(
    t: &NetTable,
    cells_level: &[String],
    remove_isolate: bool,
) -> (Vec<f64>, Vec<f64>, Vec<String>, Vec<String>) {
    assert!(!t.is_empty(), "subsetCommunication returned no rows");
    // `dplyr::group_by` orders by the key. `vctrs`'s character proxy collation defaults to
    // "C", i.e. byte order, which is what `as_bytes().cmp()` gives -- *not* the locale-aware
    // order that `sort()` on a character vector uses in a UTF-8 locale. The corpus pins
    // this: the key order is `g10|g1, g10|g10, g10|g2, ..., g1|g1, ...`, i.e. `"g10|g1"`
    // before `"g1|g1"` because `'0' < '|'`.
    let mut keys: Vec<(String, usize)> = Vec::with_capacity(t.len());
    for r in 0..t.len() {
        let s = t.get(r, "source").unwrap_or("NA");
        let tg = t.get(r, "target").unwrap_or("NA");
        keys.push((format!("{s}|{tg}"), r));
    }
    keys.sort_by(|a, b| a.0.as_bytes().cmp(b.0.as_bytes()));

    // `str_split(key, "|", simplify = TRUE)`: columns are "", one per character, "".
    // `a[, 1]` is therefore "" and `a[, 2]` is the first character; `a[, 2]` is what upstream
    // reads for *both* axes, so `target` is the first character too (upstream's `a[, 2]`,
    // not the third character -- the two reads are both wrong, and both are reproduced
    // below by taking column 1 for `source` and column 2 for `target`, which is literally
    // what the code says).
    let col = |key: &str, idx: usize| -> String {
        // Column 0 is the leading "". Column i (1-based, as `a[, i]`) is character i-1.
        if idx == 0 {
            return String::new();
        }
        match key.chars().nth(idx - 1) {
            Some(c) => c.to_string(),
            None => String::new(),
        }
    };
    let pairs: Vec<(String, String)> = keys.iter().map(|(k, _)| (col(k, 1), col(k, 2))).collect();

    let mut present_src: Vec<String> = Vec::new();
    let mut present_tgt: Vec<String> = Vec::new();
    for (s, tg) in &pairs {
        if !present_src.contains(s) {
            present_src.push(s.clone());
        }
        if !present_tgt.contains(tg) {
            present_tgt.push(tg.clone());
        }
    }
    let (src_levels, tgt_levels) = if remove_isolate {
        (
            cells_level
                .iter()
                .filter(|l| present_src.iter().any(|s| s.as_str() == l.as_str()))
                .cloned()
                .collect::<Vec<_>>(),
            cells_level
                .iter()
                .filter(|l| present_tgt.iter().any(|s| s.as_str() == l.as_str()))
                .cloned()
                .collect::<Vec<_>>(),
        )
    } else {
        (cells_level.to_vec(), cells_level.to_vec())
    };
    let ns = src_levels.len();
    let nt = tgt_levels.len();

    // `tapply(x, list(source, target), sum)`, index (source, target) with source fastest.
    // No cell group is named `""` or a single character in practice, so every entry stays
    // `NA`; the scatter is written out anyway so a level that *is* a single character (or an
    // empty name) behaves as upstream does rather than being special-cased away.
    let mut count = vec![f64::NAN; ns * nt];
    let mut weight = vec![f64::NAN; ns * nt];
    for ((key, _), (s, tg)) in keys.iter().zip(pairs.iter()) {
        let i = match src_levels.iter().position(|l| l == s) {
            Some(i) => i,
            None => continue,
        };
        let j = match tgt_levels.iter().position(|l| l == tg) {
            Some(j) => j,
            None => continue,
        };
        let mut n = 0.0f64;
        let mut p = 0.0f64;
        for r in 0..t.len() {
            let rs = t.get(r, "source").unwrap_or("NA");
            let rt = t.get(r, "target").unwrap_or("NA");
            if format!("{rs}|{rt}") == *key {
                n += 1.0;
                p += t.num(r, "prob").unwrap_or(f64::NAN);
            }
        }
        let c = i + ns * j;
        count[c] = n;
        weight[c] = p;
    }
    // `net$weight[is.na(net$weight)] <- 0` then `net$count[is.na(net$count)] <- 0`.
    for v in weight.iter_mut() {
        if v.is_nan() {
            *v = 0.0;
        }
    }
    for v in count.iter_mut() {
        if v.is_nan() {
            *v = 0.0;
        }
    }
    (count, weight, src_levels, tgt_levels)
}

#[cfg(test)]
mod tests {
    use super::*;

    fn net_2x2x2() -> Net {
        let prob = vec![
            0.0, 0.5, 0.25, 0.0, // interaction 1: (1,1)=0 (1,2)=.5 (2,1)=.25 (2,2)=0
            1.0, 0.0, 0.0, 2.0, // interaction 2
        ];
        let pval = vec![1.0, 0.01, 0.2, 0.5, 0.5, 0.01, 0.5, 0.5];
        Net::new(
            prob,
            pval,
            vec!["a".into(), "b".into()],
            vec!["i1".into(), "i2".into()],
        )
    }

    #[test]
    fn aggregate_counts_nonzero_cells_and_sums_weights() {
        let net = net_2x2x2();
        let (count, weight) = aggregate_net_default(&net, 0.05);
        // Interaction 1: (1,1) prob 0 -> pval rescued to 1 -> dropped. (1,2) 0.5 kept.
        //                      (2,1) 0.25 pval 0.2 >= 0.05 -> dropped. (2,2) 0 -> dropped.
        // Interaction 2: (1,1) prob 1 pval 0.5 -> dropped. (1,2) 0 dropped. (2,1) 0 dropped.
        //                      (2,2) prob 2 pval 0.5 -> dropped.
        assert_eq!(count, vec![0.0, 1.0, 0.0, 0.0]);
        assert_eq!(weight, vec![0.0, 0.5, 0.0, 0.0]);
    }

    #[test]
    fn the_weight_sum_is_accumulated_in_eighty_bit() {
        // 1 ulp of 1230 is 2.27e-13, which is exactly the difference an f64 accumulator
        // makes here. Guards against reintroducing a per-step narrowing.
        let terms = [0.28325667885710665, 1230.0234071066059, 0.39647518008571769];
        let net = Net::new(
            vec![
                terms[0], 0.0, 0.0, 0.0, terms[1], 0.0, 0.0, 0.0, terms[2], 0.0, 0.0, 0.0,
            ],
            vec![0.01; 12],
            vec!["a".into(), "b".into()],
            vec!["i1".into(), "i2".into(), "i3".into()],
        );
        let (_, weight) = aggregate_net_default(&net, 0.05);
        let mut want = 0.0f64;
        for t in terms {
            want += t;
        }
        let mut e80 = F80::from_f64(0.0);
        for t in terms {
            e80 = e80.add(F80::from_f64(t));
        }
        assert_eq!(weight[0].to_bits(), e80.to_f64().to_bits(), "80-bit sum");
        if want != e80.to_f64() {
            // Not an assertion that it differs -- just that the guard above is meaningful.
            eprintln!(
                "note: the f64 sum is {} and the 80-bit sum is {}",
                want,
                e80.to_f64()
            );
        }
    }

    #[test]
    fn a_pval_exactly_at_thresh_is_dropped() {
        // `prob[pval >= thresh] <- 0` is a `>=`, so the boundary is *not* kept.
        let net = Net::new(
            vec![0.4, 0.0, 0.0, 0.0],
            vec![0.05, 1.0, 1.0, 1.0],
            vec!["a".into(), "b".into()],
            vec!["i".into()],
        );
        let (count, _) = aggregate_net_default(&net, 0.05);
        assert_eq!(count[0], 0.0, "pval == thresh must be dropped");
        let (count, _) = aggregate_net_default(&net, 0.05000000000000001);
        assert_eq!(count[0], 1.0, "pval just below thresh is kept");
    }

    #[test]
    fn a_zero_probability_cell_is_rescued_to_pval_one_then_dropped() {
        let net = Net::new(
            vec![0.0, 0.0, 0.0, 0.0],
            vec![0.0, 0.0, 0.0, 0.0],
            vec!["a".into(), "b".into()],
            vec!["i".into()],
        );
        // Without `pval[prob == 0] <- 1` these would all pass `pval >= 0.05`... no, they
        // would *fail* it and be kept. The rescue is what drops them.
        let (count, weight) = aggregate_net_default(&net, 0.05);
        assert_eq!(count, vec![0.0; 4]);
        assert_eq!(weight, vec![0.0; 4]);
    }

    #[test]
    fn subset_rows_run_interaction_major_then_target_then_source() {
        let net = Net::new(
            // All four cells of interaction 1 are positive and below thresh; interaction 2
            // has only (2,2).
            vec![0.1, 0.2, 0.3, 0.4, 0.0, 0.0, 0.0, 0.9],
            vec![0.01; 8],
            vec!["a".into(), "b".into()],
            vec!["i1".into(), "i2".into()],
        );
        let lr = vec![
            LrMeta::minimal("i1", "L1^R1", "L1", "R1", "Secreted Signaling"),
            LrMeta::minimal("i2", "L2^R2", "L2", "R2", "ECM-Receptor"),
        ];
        let t = subset_communication(&net, &lr, 0.05, None, None, false, true);
        let got: Vec<(String, String, String, f64)> = (0..t.len())
            .map(|r| {
                (
                    t.get(r, "source").unwrap().to_string(),
                    t.get(r, "target").unwrap().to_string(),
                    t.get(r, "interaction_name").unwrap().to_string(),
                    t.num(r, "prob").unwrap(),
                )
            })
            .collect();
        assert_eq!(
            got,
            vec![
                // interaction 1, target a, sources a then b
                ("a".into(), "a".into(), "i1".into(), 0.1),
                ("b".into(), "a".into(), "i1".into(), 0.2),
                // interaction 1, target b
                ("a".into(), "b".into(), "i1".into(), 0.3),
                ("b".into(), "b".into(), "i1".into(), 0.4),
                // interaction 2: the only positive cell is (source b, target b)
                ("b".into(), "b".into(), "i2".into(), 0.9),
            ]
        );
    }

    #[test]
    fn subset_columns_follow_the_intersect_order() {
        let net = net_2x2x2();
        let lr = vec![LrMeta::minimal(
            "i1",
            "L1^R1",
            "L1",
            "R1",
            "Secreted Signaling",
        )];
        let t = subset_communication(&net, &lr, 0.05, None, None, false, true);
        // `evidence` is NA in every row and `pathway_name` is "", yet both stay: `intersect`
        // looks at `colnames`, not at the values. Dropping a column because its cells are NA
        // is the bug this pins.
        assert_eq!(t.columns, COLUMNS_PLAIN.to_vec());
        assert!(t.get(0, "evidence").is_none(), "evidence is present but NA");
        assert!(t.columns.contains(&"evidence".to_string()));
    }
}
