//! `subsetCommunication`'s **post-melt** half: the eight DEG thresholds, the all-`NA` row
//! drop, the `netP` aggregation and the final column selection.
//!
//! `subset_communication` in [`crate::net`] covers everything up to and including the
//! `reshape2::melt`, the `thresh` cut and the L-R join, and works from the 3-d `prob`/`pval`
//! arrays. This module starts where that one stops, on an already-melted `net` data frame,
//! which is also the form a caller passes in as `subsetCommunication(net = <data.frame>)`.
//! Keeping the two halves separate is not an aesthetic choice: the melt path needs the
//! arrays, and the DEG path cannot use them, because the DEG columns only exist once
//! `netMappingDEG` has annotated the table.
//!
//! ## The `NA` row rule, which is the whole reason this is not eight `filter()` calls
//!
//! Upstream filters with `net[net$ligand.pvalues <= ligand.pvalues, , drop = FALSE]`. For a
//! row whose value is `NA` -- or `NaN`, which compares `NA` in R just as `NA` does -- the
//! index entry is `NA`, and indexing a data frame by a logical vector containing `NA`
//! **inserts a row of `NA`s** rather than dropping the row. So a threshold does not
//! silently discard rows with a missing value: it blanks the whole row, and the row is
//! removed later by
//!
//! ```r
//! net <- net[rowSums(is.na(net)) != ncol(net), , drop = FALSE]
//! ```
//!
//! which is reached unconditionally. A later threshold applied to an all-`NA` row yields
//! `NA` again, so the row stays blank through the rest of the chain and is dropped exactly
//! once, at the end. The net effect is "rows with a missing value in any *applied* threshold
//! column are dropped" -- but the all-`NA` drop is a **separate** step, because a row that
//! was already entirely `NA` on input is dropped even when no threshold is supplied. Both
//! are implemented below; they are not the same test.
//!
//! `datasets` and the two `sources.use` / `targets.use` filters are the exceptions: `%in%`
//! is `match()`, which returns `FALSE` for a missing value and never `NA`, so those drop
//! the row outright and insert nothing.

use crate::longdouble::F80;
use crate::net::NetTable;
use std::collections::HashMap;

/// Which slot's column list `subsetCommunication` finishes with.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum NetSlot {
    /// `slot.name = "net"`: the L-R-level table.
    Net,
    /// `slot.name = "netP"`: aggregated to `(source, target, pathway)`.
    NetP,
}

impl NetSlot {
    /// Upstream compares with `==` against the argument, so an unrecognised value falls
    /// through both branches of the final `if/else` and keeps every column. `None` models
    /// that.
    pub fn from_name(name: &str) -> Option<Self> {
        match name {
            "net" => Some(NetSlot::Net),
            "netP" => Some(NetSlot::NetP),
            _ => None,
        }
    }
}

/// The eight DEG thresholds plus the two cell-group filters, all optional as upstream has
/// them: each is an independent `if (!is.null(...))`, so supplying several applies all.
#[derive(Clone, Debug, Default)]
pub struct DegThresholds {
    pub datasets: Option<Vec<String>>,
    pub ligand_pvalues: Option<f64>,
    pub ligand_logfc: Option<f64>,
    pub ligand_pct1: Option<f64>,
    pub ligand_pct2: Option<f64>,
    pub receptor_pvalues: Option<f64>,
    pub receptor_logfc: Option<f64>,
    pub receptor_pct1: Option<f64>,
    pub receptor_pct2: Option<f64>,
    pub sources_use: Option<Vec<String>>,
    pub targets_use: Option<Vec<String>>,
    pub slot: Option<NetSlot>,
}

/// The errors `subsetCommunication_internal` raises, each with upstream's exact text.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum SubsetCommError {
    /// One of the eight shared `stop()`s.
    NeedsDegThreshold { column: String },
    /// `stop("Please run `identifyOverExpressedGenes` and `netMappingDEG` before selecting
    /// 'datasets'")` -- worded differently from the eight, and the test asserts it.
    NeedsDegDatasets,
    /// `stop("No significant signaling interactions are inferred based on the input!")`.
    EmptyAfterThreshold,
}

impl std::fmt::Display for SubsetCommError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            SubsetCommError::NeedsDegThreshold { column } => write!(
                f,
                "Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the \
                 threshold '{column}'"
            ),
            SubsetCommError::NeedsDegDatasets => write!(
                f,
                "Please run `identifyOverExpressedGenes` and `netMappingDEG` before selecting \
                 'datasets'"
            ),
            SubsetCommError::EmptyAfterThreshold => write!(
                f,
                "No significant signaling interactions are inferred based on the input!"
            ),
        }
    }
}

impl std::error::Error for SubsetCommError {}

/// A cell as a number, treating `NaN` the way R's comparisons do.
///
/// `NetTable` stores `NA` as `None` and lets a literal `NaN` through as `Some(NaN)`, but
/// `NaN >= x` is `NA` in R exactly as `NA >= x` is, so both must take the blank-the-row
/// branch. Getting this wrong keeps every `NaN` row instead of dropping it.
#[inline]
fn num(of: &NetTable, r: usize, name: &str) -> Option<f64> {
    match of.num(r, name) {
        Some(v) if !v.is_nan() => Some(v),
        _ => None,
    }
}

/// `rowSums(is.na(x)) != ncol(x)`: drop rows whose every cell is `NA`.
fn drop_all_na(net: &NetTable) -> NetTable {
    let ncol = net.columns.len();
    let rows = net
        .rows
        .iter()
        .filter(|r| r.iter().filter(|c| c.is_none()).count() != ncol)
        .cloned()
        .collect();
    NetTable {
        rows,
        ..net.clone()
    }
}

/// `net[net[[name]] <cmp> t, ]`, where `cmp` can yield `NA` -- and an `NA` index entry
/// becomes an all-`NA` row, not a dropped one.
fn filter_cmp<F>(net: &NetTable, name: &str, cmp: F) -> Result<NetTable, SubsetCommError>
where
    F: Fn(f64) -> bool,
{
    if !net.columns.iter().any(|c| c == name) {
        return Err(SubsetCommError::NeedsDegThreshold {
            column: name.to_string(),
        });
    }
    let ncol = net.columns.len();
    let mut rows: Vec<Vec<Option<String>>> = Vec::new();
    for (r, row) in net.rows.iter().enumerate() {
        match num(net, r, name) {
            Some(v) if cmp(v) => rows.push(row.clone()),
            Some(_) => {}
            None => rows.push(vec![None; ncol]),
        }
    }
    Ok(NetTable {
        rows,
        ..net.clone()
    })
}

/// `dplyr::group_by` on character columns: sorted **byte-wise**, not by locale collation.
///
/// The same ordering as `aggregateNet`'s `group_by`, and for the same reason: with values
/// `g1`, `g10`, `g2` R orders `"g10"` before `"g1"`, because the pinned environment collates
/// byte-wise. A fixture whose group keys distinguish the two orders is what keeps this
/// honest -- a port that sorted by (pathway, source, target) would pass on three groups.
fn group_key_order(keys: &[(String, String)]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..keys.len()).collect();
    idx.sort_by(|&a, &b| {
        keys[a]
            .0
            .as_bytes()
            .cmp(keys[b].0.as_bytes())
            .then_with(|| keys[a].1.as_bytes().cmp(keys[b].1.as_bytes()))
    });
    idx
}

/// `stringr::str_split(x, "sourceTotarget", simplify = TRUE)`, then `a[, 1]` / `a[, 2]`.
fn split_source_target<'a>(v: &'a str, sep: &str) -> (&'a str, &'a str) {
    match v.split_once(sep) {
        Some((a, b)) => (a, b),
        // `str_split` on a string with no separator yields a length-1 vector, so `a[, 2]`
        // is `NA`. `paste()` always inserts the separator, so this is unreachable for
        // upstream's own output; it exists so a hand-built fixture cannot diverge silently.
        None => (v, ""),
    }
}

/// `subsetCommunication_internal` from the already-melted `net` onwards.
pub fn subset_communication_deg(
    net: &NetTable,
    th: &DegThresholds,
) -> Result<NetTable, SubsetCommError> {
    let mut cur = net.clone();

    // Upstream's order: datasets, then the four ligand thresholds, then the four receptor
    // ones. Each is an independent `if`, so the order is observable: ask for two
    // thresholds whose columns are both absent and the *first* one's message is the one
    // raised. The fixtures pin that.
    if let Some(d) = &th.datasets {
        cur = filter_in_set(&cur, "datasets", d)?;
    }
    if let Some(t) = th.ligand_pvalues {
        cur = filter_cmp(&cur, "ligand.pvalues", |v| v <= t)?;
    }
    if let Some(t) = th.ligand_logfc {
        // `if (ligand.logFC >= 0) >= else <=`: a *negative* threshold keeps the
        // downregulated genes. The sign of the **argument** chooses, not the data -- so
        // `ligand.logFC = 0` takes the `>=` branch and keeps only `logFC >= 0`.
        cur = if t >= 0.0 {
            filter_cmp(&cur, "ligand.logFC", |v| v >= t)?
        } else {
            filter_cmp(&cur, "ligand.logFC", |v| v <= t)?
        };
    }
    if let Some(t) = th.ligand_pct1 {
        cur = filter_cmp(&cur, "ligand.pct.1", |v| v >= t)?;
    }
    if let Some(t) = th.ligand_pct2 {
        cur = filter_cmp(&cur, "ligand.pct.2", |v| v >= t)?;
    }
    if let Some(t) = th.receptor_pvalues {
        cur = filter_cmp(&cur, "receptor.pvalues", |v| v <= t)?;
    }
    if let Some(t) = th.receptor_logfc {
        cur = if t >= 0.0 {
            filter_cmp(&cur, "receptor.logFC", |v| v >= t)?
        } else {
            filter_cmp(&cur, "receptor.logFC", |v| v <= t)?
        };
    }
    if let Some(t) = th.receptor_pct1 {
        cur = filter_cmp(&cur, "receptor.pct.1", |v| v >= t)?;
    }
    if let Some(t) = th.receptor_pct2 {
        cur = filter_cmp(&cur, "receptor.pct.2", |v| v >= t)?;
    }

    cur = drop_all_na(&cur);
    if cur.rows.is_empty() {
        return Err(SubsetCommError::EmptyAfterThreshold);
    }

    if th.slot == Some(NetSlot::NetP) {
        cur = aggregate_netp(&cur)?;
    }

    if let Some(s) = &th.sources_use {
        cur = filter_membership(&cur, "source", s);
    }
    if let Some(s) = &th.targets_use {
        cur = filter_membership(&cur, "target", s);
    }

    let Some(slot) = th.slot else {
        // Neither branch of upstream's final `if/else` runs, so every column survives.
        return Ok(cur);
    };
    let keep: Vec<&str> = match slot {
        NetSlot::Net => {
            const BASE: [&str; 11] = [
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
            const DEG: [&str; 8] = [
                "ligand.logFC",
                "ligand.pct.1",
                "ligand.pct.2",
                "ligand.pvalues",
                "receptor.logFC",
                "receptor.pct.1",
                "receptor.pct.2",
                "receptor.pvalues",
            ];
            // The `datasets` branch needs `ligand.logFC` *as well*, so a table carrying
            // `datasets` but no `ligand.logFC` falls to the plain `else` and drops
            // `datasets`. A fixture pins that, because it looks like a bug and is not.
            let has = |n: &str| cur.columns.iter().any(|c| c == n);
            // `datasets` sits *between* `evidence` and the DEG list, not after it. The
            // `dplyr::select` renames to the requested order, so the column order is this
            // vector's order, and the corpus compares column names as well as cells.
            if has("ligand.logFC") && has("datasets") {
                BASE.iter()
                    .copied()
                    .chain(["datasets"])
                    .chain(DEG.iter().copied())
                    .collect()
            } else if has("ligand.logFC") {
                BASE.iter().chain(DEG.iter()).copied().collect()
            } else {
                BASE.to_vec()
            }
        }
        NetSlot::NetP => vec!["source", "target", "pathway_name", "prob", "pval"],
    };
    // `intersect(requested, colnames(net))`: the *requested* order, minus the names the
    // table does not have.
    let idx: Vec<usize> = keep
        .iter()
        .filter_map(|n| cur.columns.iter().position(|c| c == n))
        .collect();
    let columns: Vec<String> = idx.iter().map(|&i| cur.columns[i].clone()).collect();
    let rows = cur
        .rows
        .iter()
        .map(|r| idx.iter().map(|&i| r[i].clone()).collect())
        .collect();
    Ok(NetTable {
        columns,
        rows,
        ..cur
    })
}

/// `subset(net, source %in% sources.use)`. `%in%` is `match()`: a missing value is simply
/// not a member, so the row is dropped and nothing is inserted.
fn filter_membership(net: &NetTable, name: &str, keep: &[String]) -> NetTable {
    let rows = net
        .rows
        .iter()
        .enumerate()
        .filter(|(r, _)| match net.get(*r, name) {
            Some(v) => keep.iter().any(|k| k == v),
            None => false,
        })
        .map(|(_, r)| r.clone())
        .collect();
    NetTable {
        rows,
        ..net.clone()
    }
}

fn filter_in_set(net: &NetTable, name: &str, keep: &[String]) -> Result<NetTable, SubsetCommError> {
    if !net.columns.iter().any(|c| c == name) {
        return Err(SubsetCommError::NeedsDegDatasets);
    }
    Ok(filter_membership(net, name, keep))
}

/// The `slot.name = "netP"` aggregation.
///
/// The column *order* is reproduced too, and it is not the order anyone would write:
/// `summarize` returns the group keys then the aggregate, `source`/`target` are appended
/// after that, `source_target` is dropped, and `pval` is assigned last. So the result is
/// `pathway_name, prob, source, target, pval` -- `pval` ends up in the **last** position
/// even though it is the first thing grouped.
fn aggregate_netp(net: &NetTable) -> Result<NetTable, SubsetCommError> {
    // `col.use <- intersect(c("source","target","pathway_name","prob","pval","annotation"), ...)`
    const PRE: [&str; 6] = [
        "source",
        "target",
        "pathway_name",
        "prob",
        "pval",
        "annotation",
    ];
    let present: Vec<(usize, usize)> = PRE
        .iter()
        .enumerate()
        .filter_map(|(p, n)| net.columns.iter().position(|c| c == n).map(|i| (p, i)))
        .collect();
    // `summarize` needs the two group keys plus both aggregates; a frame missing any of them
    // fails upstream inside dplyr with a different message, and the fixtures always supply
    // all four, so this arm only guards a malformed hand-built table.
    if present.len() < 5 {
        return Err(SubsetCommError::EmptyAfterThreshold);
    }
    let col = |name: &str| -> usize {
        let p = PRE.iter().position(|c| *c == name).unwrap();
        present.iter().find(|(pp, _)| *pp == p).unwrap().1
    };
    let (i_src, i_tgt, i_pw) = (col("source"), col("target"), col("pathway_name"));
    let (i_prob, i_pval) = (col("prob"), col("pval"));

    // `paste(source, target, sep = "sourceTotarget")` on an `NA` gives the literal string
    // "NA", *not* `NA_character_`: `paste` coerces with `as.character`, and
    // `as.character(NA)` is "NA". So a missing cell becomes a real group key. That is the
    // one place in this function where a missing value survives as data, and it is why the
    // all-`NA` row drop has to run *before* the aggregation -- it does.
    const SEP: &str = "sourceTotarget";
    let cell = |r: &[Option<String>], i: usize| -> String {
        r[i].clone().unwrap_or_else(|| "NA".to_string())
    };
    let keys: Vec<(String, String)> = net
        .rows
        .iter()
        .map(|r| {
            (
                format!("{}{}{}", cell(r, i_src), SEP, cell(r, i_tgt)),
                cell(r, i_pw),
            )
        })
        .collect();
    let order = group_key_order(&keys);
    // One row per *distinct* key, in sorted-key order -- `dplyr::summarize` collapses each group.
    // `group_key_order` returns every row index sorted, so without this the emission loop below
    // writes one row per input row and a table with duplicate `(source_target, pathway)` keys
    // comes back unaggregated. No fixture caught it because every corpus case has unique keys;
    // the tutorial's real `netP` was the first table with two interactions in the same pathway
    // between the same groups (576 distinct keys in 1042 rows, emitted as 1042).
    let mut seen: std::collections::HashSet<&(String, String)> =
        std::collections::HashSet::with_capacity(order.len());
    let order: Vec<usize> = order
        .into_iter()
        .filter(|&r| seen.insert(&keys[r]))
        .collect();

    // `mean(pval)` and `sum(prob)` both accumulate in LONG_DOUBLE and round once at the end.
    // Neither has an `na.rm`, so a missing member poisons the group -- but `NA` and `NaN` are
    // poisoned *differently*, and the corpus distinguishes them: a group containing an `NA`
    // reports `<NA>`, a group containing only `NaN` reports `<NaN>`. Both are NaN payloads
    // to `f64::is_nan`, so a single `is_nan` check cannot tell them apart and the two have to
    // be tracked as flags alongside the accumulator. `NA` wins when both are present, which
    // is what R's arithmetic does: `NA_real_` has a payload that survives addition.
    let mut acc: HashMap<(String, String), Acc> = HashMap::new();
    for (r, key) in keys.iter().enumerate() {
        let e = acc.entry(key.clone()).or_default();
        e.add(&net.rows[r][i_pval], &net.rows[r][i_prob]);
    }

    let mut rows: Vec<Vec<Option<String>>> = Vec::with_capacity(order.len());
    for &r in &order {
        let key = &keys[r];
        let a = &acc[key];
        let (src, tgt) = split_source_target(&key.0, SEP);
        rows.push(vec![
            Some(key.1.clone()),
            a.prob_text(),
            Some(src.to_string()),
            Some(tgt.to_string()),
            a.pval_text(),
        ]);
    }
    Ok(NetTable {
        columns: vec![
            "pathway_name".into(),
            "prob".into(),
            "source".into(),
            "target".into(),
            "pval".into(),
        ],
        rows,
        group_levels: net.group_levels.clone(),
        // This is the pathway table, so it has no `interaction_name` column and nothing needs the
        // levels. The field is carried anyway so the struct has one shape; the pathway factor's
        // levels are the table's own `pathway_name` values, which `pathway_parity` checks directly.
        interaction_levels: Vec::new(),
    })
}

/// One cell as `(value, is_na)`. Both `None` and the literal `"NA"` are `NA`; `NaN` is not.
fn cell_f64(cell: &Option<String>) -> (f64, bool) {
    match cell {
        None => (f64::NAN, true),
        Some(s) if s.trim() == "NA" => (f64::NAN, true),
        Some(s) => (s.trim().parse().unwrap_or(f64::NAN), false),
    }
}

/// One `(source_target, pathway_name)` group's running `sum(prob)` and `mean(pval)`.
struct Acc {
    /// `mean(pval)`: the LONG_DOUBLE accumulator, its divisor, and the two missing-value
    /// flags. Named for the quantity, not for a position, because the two accumulators were
    /// once swapped and the flags travelled with them -- which showed up as a `sum(prob)`
    /// column reporting `NaN` whenever a group's *pval* was `NaN`.
    msum: F80,
    mn: u64,
    m_na: bool,
    m_nan: bool,
    /// `sum(prob)`: no divisor, and the same two flags.
    ssum: F80,
    s_na: bool,
    s_nan: bool,
}

impl Default for Acc {
    fn default() -> Self {
        Self {
            msum: F80::from_f64(0.0),
            mn: 0,
            m_na: false,
            m_nan: false,
            ssum: F80::from_f64(0.0),
            s_na: false,
            s_nan: false,
        }
    }
}

impl Acc {
    fn add(&mut self, pval: &Option<String>, prob: &Option<String>) {
        // `mn` counts *every* member, missing ones included: R divides by the full length even
        // when the result is `NA`, and it matters for the value in the (impossible) case of a
        // group whose members are all present but sum to an exact zero.
        self.mn += 1;
        let (v, is_na) = cell_f64(pval);
        self.m_na |= is_na;
        self.m_nan |= !is_na && v.is_nan();
        self.msum = self.msum.add(F80::from_f64(v));

        let (v, is_na) = cell_f64(prob);
        self.s_na |= is_na;
        self.s_nan |= !is_na && v.is_nan();
        self.ssum = self.ssum.add(F80::from_f64(v));
    }

    /// `None` for `NA`, `Some("NaN")` for `NaN`. **Not** the strings `"<NA>"` / `"<NaN>"`:
    /// those are the *corpus's* escape for a tab-separated file, and a `NetTable` cell is a
    /// typed value. Spelling the escape into the cell meant the binding's `"NA"` / `"NaN"`
    /// match fell through to `parse()`, both became `f64::NAN`, and a `mean` that upstream
    /// reports as `NA` came back as `NaN` -- a difference `identical()` sees and a numeric
    /// comparison does not.
    fn prob_text(&self) -> Option<String> {
        if self.s_na {
            None
        } else if self.s_nan {
            Some("NaN".into())
        } else {
            Some(fmt_f64(self.ssum.to_f64()))
        }
    }

    fn pval_text(&self) -> Option<String> {
        if self.m_na {
            None
        } else if self.m_nan {
            Some("NaN".into())
        } else {
            Some(fmt_f64(self.msum.div_int(self.mn.max(1)).to_f64()))
        }
    }
}

/// `as.character(x)` for a double, at 17 significant digits.
///
/// Only used by the golden corpus and by the R side, both of which re-parse the text, so
/// the requirement is an exact round trip rather than R's 15-digit display width. 17 digits
/// is the shortest width that is guaranteed to round-trip every `f64`.
pub fn fmt_f64(v: f64) -> String {
    if v.is_nan() {
        return "NaN".to_string();
    }
    if v.is_infinite() {
        return if v > 0.0 { "Inf".into() } else { "-Inf".into() };
    }
    format!("{v:.17e}")
}

#[cfg(test)]
mod netp_tests {
    use super::aggregate_netp;
    use crate::net::NetTable;

    fn table(rows: Vec<Vec<Option<String>>>) -> NetTable {
        NetTable {
            columns: vec![
                "source".into(),
                "target".into(),
                "pathway_name".into(),
                "prob".into(),
                "pval".into(),
                "annotation".into(),
            ],
            rows,
            group_levels: vec!["g1".into(), "g2".into()],
            interaction_levels: vec![],
        }
    }

    fn row(src: &str, tgt: &str, pw: &str, prob: &str, pval: &str) -> Vec<Option<String>> {
        vec![
            Some(src.into()),
            Some(tgt.into()),
            Some(pw.into()),
            Some(prob.into()),
            Some(pval.into()),
            Some("Secreted Signaling".into()),
        ]
    }

    #[test]
    fn duplicate_keys_collapse_to_one_row_per_group() {
        // `dplyr::summarize` emits one row per group. The emission loop used to walk every *row*
        // index in sorted order, so a table with two interactions in the same pathway between the
        // same groups came back unaggregated -- 1042 rows for 576 distinct keys on the tutorial's
        // real `netP`, where no corpus case has a duplicate key.
        let t = table(vec![
            row("g1", "g2", "WNT", "0.5", "0.01"),
            row("g1", "g2", "WNT", "0.25", "0.03"),
            row("g2", "g1", "WNT", "1.0", "0.0"),
        ]);
        let out = aggregate_netp(&t).expect("aggregates");
        assert_eq!(
            out.rows.len(),
            2,
            "one row per distinct (source, target, pathway)"
        );
        // Output order is `[pathway_name, prob, source, target, pval]` -- `pval` last, per the
        // doc comment on `aggregate_netp`.
        let keys: Vec<String> = out
            .rows
            .iter()
            .map(|r| {
                format!(
                    "{}|{}|{}",
                    r[2].clone().unwrap(),
                    r[3].clone().unwrap(),
                    r[0].clone().unwrap()
                )
            })
            .collect();
        assert!(keys.contains(&"g1|g2|WNT".to_string()), "{keys:?}");
        assert!(keys.contains(&"g2|g1|WNT".to_string()), "{keys:?}");
        // `summarize(prob = sum(prob))`: 0.5 + 0.25, accumulated in LONG_DOUBLE like R's sum.
        let g1g2 = out
            .rows
            .iter()
            .find(|r| r[2].clone().unwrap() == "g1")
            .unwrap();
        assert_eq!(g1g2[1].clone().unwrap().parse::<f64>().unwrap(), 0.75);
        // `summarize(pval = mean(pval))`: (0.01 + 0.03) / 2.
        assert_eq!(g1g2[4].clone().unwrap().parse::<f64>().unwrap(), 0.02);
    }

    #[test]
    fn groups_come_out_in_sorted_key_order() {
        // `dplyr::group_by` sorts character keys byte-wise; the second group here sorts first.
        let t = table(vec![
            row("g2", "g1", "WNT", "1.0", "0.0"),
            row("g1", "g2", "BMP", "0.5", "0.5"),
        ]);
        let out = aggregate_netp(&t).expect("aggregates");
        assert_eq!(out.rows.len(), 2);
        let pws: Vec<String> = out.rows.iter().map(|r| r[0].clone().unwrap()).collect();
        assert_eq!(pws, vec!["BMP".to_string(), "WNT".to_string()]);
    }
}
