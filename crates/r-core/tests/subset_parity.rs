//! Parity for `subsetCommunication_internal`'s post-melt half: the eight DEG thresholds,
//! the all-`NA` row drop, the `netP` aggregation and the final column selection.
//!
//! The corpus is `tests/fixtures/subset_golden.txt`, written by
//! `tests/parity/gen_subset_golden.R` from the pinned upstream's own body. There is no RNG
//! in this branch -- the input table is hand-built -- so the corpus is reproducible from the
//! generator alone and the failures it pins are all deterministic.
//!
//! The comparison is on the *text* of every cell rather than on parsed doubles, because the
//! distinctions that matter here are exactly the ones a numeric comparison would erase: an
//! `NA` and a `NaN` compare equal to nothing, and the string `"NA"` that `paste` produces
//! from a missing cell group is not the same value as a missing cell.

use r_core::net::NetTable;
use r_core::subset::{DegThresholds, NetSlot, SubsetCommError};
use std::collections::HashMap;
use std::path::Path;

#[derive(Debug, Default, Clone)]
struct Case {
    in_cols: Vec<String>,
    in_rows: Vec<Vec<Option<String>>>,
    out_cols: Option<Vec<String>>,
    out_rows: Vec<Vec<Option<String>>>,
    out_rn: Option<String>,
    error: Option<String>,
    warn: Vec<String>,
    slot: String,
    arg: HashMap<String, String>,
}

/// Warnings raised by *CellChat* itself, as opposed to by the versions of `dplyr` and
/// `tidyselect` installed alongside it.
///
/// The `netP` path goes through `dplyr::group_by`/`summarize`/`select`, and modern
/// tidyselect warns "Using an external vector in selections was deprecated" and "Setting row
/// names on a tibble is deprecated". Those are the *dependency's* messages, pinned to a
/// dplyr release rather than to CellChat, and reproducing them from Rust would mean emitting
/// another package's version-specific deprecation text -- which would be wrong on the next
/// dplyr release. The shim therefore does not emit them, and they are recorded here only so
/// that the corpus is a faithful dump. Everything else must match exactly.
const DEPENDENCY_WARNINGS: [&str; 2] = [
    "Using an external vector in selections was deprecated",
    "Setting row names on a tibble is deprecated.",
];

/// Corpus text -> a cell. `<NA>` is a missing cell, `NA` is the *string* `NA`.
fn cell(s: &str) -> Option<String> {
    match s {
        "<NA>" => None,
        "<NaN>" => Some("NaN".to_string()),
        other => Some(other.replace("<TAB>", "\t")),
    }
}

fn parse_num(s: &str) -> Option<f64> {
    match s {
        "-" => None,
        other => other.parse().ok(),
    }
}

fn parse_set(s: &str) -> Option<Vec<String>> {
    match s {
        "-" => None,
        other => Some(other.split(',').map(|x| x.to_string()).collect()),
    }
}

fn load() -> Vec<(String, Case)> {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("../../tests/fixtures/subset_golden.txt");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    let mut out: Vec<(String, Case)> = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.is_empty() {
            continue;
        }
        match f[0] {
            "case" => out.push((f[1].to_string(), Case::default())),
            "in_cols" => {
                out.last_mut().unwrap().1.in_cols = f[2].split(',').map(str::to_string).collect()
            }
            "in_row" => {
                let c = &mut out.last_mut().unwrap().1;
                assert_eq!(
                    c.in_rows.len() + 1,
                    f[2].parse::<usize>().unwrap(),
                    "row order"
                );
                c.in_rows.push(f[3..].iter().map(|x| cell(x)).collect());
            }
            "slot" => out.last_mut().unwrap().1.slot = f[2].to_string(),
            "arg" => {
                out.last_mut()
                    .unwrap()
                    .1
                    .arg
                    .insert(f[2].to_string(), f[3].to_string());
            }
            "error" => out.last_mut().unwrap().1.error = Some(f[2].replace("<TAB>", "\t")),
            "out_cols" => {
                out.last_mut().unwrap().1.out_cols =
                    Some(f[2].split(',').map(str::to_string).collect())
            }
            "out_nrow" => {}
            "out_row" => {
                let c = &mut out.last_mut().unwrap().1;
                assert_eq!(
                    c.out_rows.len() + 1,
                    f[2].parse::<usize>().unwrap(),
                    "row order"
                );
                c.out_rows.push(f[3..].iter().map(|x| cell(x)).collect());
            }
            // `out_rn` and `warn` are (tag, case, value): three fields, so the value is
            // `f[2]`. A zero-row result has no row names and the field is then empty, which
            // is why `f[2]` is read as an `Option` rather than indexed.
            "out_rn" => {
                out.last_mut().unwrap().1.out_rn = Some(f.get(2).copied().unwrap_or("").to_string())
            }
            "warn" => {
                let c = &mut out.last_mut().unwrap().1;
                if f.get(2).copied().unwrap_or("-") != "-" {
                    c.warn = f[2]
                        .split(" || ")
                        .map(|s| s.replace("<TAB>", "\t").replace("<NL>", "\n"))
                        .collect();
                }
            }
            "n_cases" | "in_nrow" => {}
            other => panic!("unknown corpus record {other:?}"),
        }
    }
    out
}

fn thresholds(c: &Case) -> DegThresholds {
    let g = |k: &str| c.arg.get(k).map(|s| s.as_str()).unwrap_or("-");
    DegThresholds {
        datasets: parse_set(g("datasets")),
        ligand_pvalues: parse_num(g("ligand.pvalues")),
        ligand_logfc: parse_num(g("ligand.logFC")),
        ligand_pct1: parse_num(g("ligand.pct.1")),
        ligand_pct2: parse_num(g("ligand.pct.2")),
        receptor_pvalues: parse_num(g("receptor.pvalues")),
        receptor_logfc: parse_num(g("receptor.logFC")),
        receptor_pct1: parse_num(g("receptor.pct.1")),
        receptor_pct2: parse_num(g("receptor.pct.2")),
        sources_use: parse_set(g("sources.use")),
        targets_use: parse_set(g("targets.use")),
        slot: NetSlot::from_name(&c.slot),
    }
}

fn table(columns: &[String], rows: &[Vec<Option<String>>]) -> NetTable {
    NetTable {
        columns: columns.to_vec(),
        rows: rows.to_vec(),
        group_levels: vec!["g1".into(), "g2".into(), "g10".into()],
        interaction_levels: vec!["LR1".into(), "LR2".into()],
    }
}

fn run_all() -> Vec<(String, Result<NetTable, SubsetCommError>)> {
    load()
        .into_iter()
        .map(|(n, c)| {
            let t = table(&c.in_cols, &c.in_rows);
            let r = r_core::subset::subset_communication_deg(&t, &thresholds(&c));
            (n, r)
        })
        .collect()
}

/// Render a table for comparison, normalising every cell that parses as a double to its
/// **bit pattern**.
///
/// The corpus is written by R's `format(digits = 17)` and the port writes Rust's `{:.17e}`,
/// and those two agree on the value while differing on width and exponent spelling
/// (`4.0000000000000002e-01` against `4.00000000000000022e-1`). Comparing the text would pin
/// the formatting rather than the arithmetic and would fail on a value that is in fact
/// identical. The bit pattern is the thing that must match, and it also makes a 1-ulp
/// difference impossible to overlook, because it changes the pattern.
fn norm(c: &Option<String>) -> String {
    match c {
        None => "<NA>".to_string(),
        Some(s) => match s.trim().parse::<f64>() {
            Ok(v) if v.is_finite() => format!("{:016x}", v.to_bits()),
            Ok(v) if v.is_nan() => {
                // `NA` and `NaN` arrive as distinct markers, but a group that sums to `NaN`
                // for want of a member must not be confused with one that has an `NA` member.
                if s.trim() == "NaN" {
                    "<NaN>".to_string()
                } else {
                    "<NA>".to_string()
                }
            }
            _ => s.clone(),
        },
    }
}

fn dump(t: &NetTable) -> String {
    let mut s = t.columns.join(",");
    s.push('\n');
    for r in &t.rows {
        s.push_str(&r.iter().map(norm).collect::<Vec<_>>().join("|"));
        s.push('\n');
    }
    s
}

/// The corpus must actually contain every branch the implementation has, or a green run
/// means nothing. Cheap, and it is the test that fails when a case is accidentally dropped.
#[test]
fn the_corpus_covers_every_branch() {
    let all = load();
    let names: Vec<&str> = all.iter().map(|(n, _)| n.as_str()).collect();
    for want in [
        "no_thresholds",
        "all_na_row_no_threshold",
        "all_eight_permissive",
        "all_eight_tight",
        "ligand_logfc_neg",
        "ligand_logfc_zero",
        "receptor_logfc_neg",
        "datasets_d1",
        "datasets_d2_nodeg",
        "netp_plain",
        "netp_after_deg",
        "netp_sources_targets",
        "netp_allna_row",
        "unknown_slot",
        "unknown_slot_deg",
    ] {
        assert!(names.contains(&want), "corpus lost the {want:?} case");
    }
    // Named rather than counted: a count silently passes if one error case is replaced by
    // another, and stops passing if an unrelated case starts erroring.
    let mut errors: Vec<&str> = all
        .iter()
        .filter(|(_, c)| c.error.is_some())
        .map(|(n, _)| n.as_str())
        .collect();
    errors.sort_unstable();
    assert_eq!(
        errors,
        [
            "err_datasets_missing",
            "err_empty_result",
            "err_ligand_logfc_missing",
            "err_ligand_pvalues_missing",
            "err_receptor_pct2_missing",
            "err_two_absent_order",
            "ligand_pvalues_none",
        ],
        "the set of error cases changed"
    );
    for want in [
        "single_threshold_blanking",
        "all_eight_permissive",
        "all_eight_tight",
    ] {
        assert!(names.contains(&want), "corpus lost the {want:?} case");
    }
}

/// The headline: for every case, either the same error text or a byte-identical table.
#[test]
fn the_result_matches_upstream_cell_for_cell() {
    let all = load();
    for (name, got) in run_all() {
        let c = &all.iter().find(|(n, _)| *n == name).unwrap().1;
        match (&c.error, &got) {
            (Some(want), Err(e)) => assert_eq!(&e.to_string(), want, "{name}: error text"),
            (Some(want), Ok(_)) => {
                panic!("{name}: upstream errors with {want:?}, port returned a table")
            }
            (None, Err(e)) => panic!("{name}: port errored with {e}, upstream returned a table"),
            (None, Ok(t)) => {
                let want = dump(&table(c.out_cols.as_ref().unwrap(), &c.out_rows));
                assert_eq!(dump(t), want, "{name}: table");
                // `rownames(net) <- 1:nrow(net)` is applied only in the non-empty branch, so
                // a zero-row result keeps the data.frame's default `integer(0)`.
                let rn = c.out_rn.clone().unwrap_or_default();
                assert_eq!(
                    rn.split(',').filter(|x| !x.is_empty()).count(),
                    t.rows.len(),
                    "{name}: rownames are assigned 1:nrow only when nrow > 0"
                );
                assert_cellchat_warnings(&name, c);
            }
        }
    }
}

/// Every warning in the corpus must be either a CellChat message the port also produces, or
/// one of the pinned dependency deprecations. Anything else means the corpus grew a new
/// observable that nothing is checking.
fn assert_cellchat_warnings(name: &str, c: &Case) {
    for w in &c.warn {
        // Prefix, not equality: the tidyselect message embeds the offending call site and a
        // URL, so it changes shape between releases while staying the same deprecation. The
        // *prefix* is the stable part, and matching it exactly would make the suite fail for
        // a reason that has nothing to do with this port.
        let known = DEPENDENCY_WARNINGS.iter().any(|p| w.starts_with(p));
        assert!(
            known || w == CELLCHAT_EMPTY_WARNING,
            "{name}: unrecognised warning {w:?} -- the corpus grew an observable"
        );
    }
    // CellChat's own warning fires on exactly one condition: the table is empty by the time
    // the final `nrow` check runs. Anything else would be a message with no cause.
    let empty = c.out_cols.is_some() && c.out_rows.is_empty();
    assert_eq!(
        c.warn.iter().any(|w| w == CELLCHAT_EMPTY_WARNING),
        empty,
        "{name}: the empty-table warning must fire if and only if the result is empty"
    );
}

/// `warning("No significant signaling interactions are inferred!")` -- no "based on the
/// input", which is the `stop()` two stages earlier. The two are easy to confuse.
const CELLCHAT_EMPTY_WARNING: &str = "No significant signaling interactions are inferred!";

/// `subsetCommunication` **warns** where it does not error, and the two are different paths:
/// the threshold stage stops with `stop()`, while `sources.use`/`targets.use` run after it and
/// an over-restrictive group list leaves an empty table with only a warning.
#[test]
fn an_over_restrictive_group_list_warns_where_a_threshold_errors() {
    let all = load();
    let c = &all.iter().find(|(n, _)| n == "sources_missing").unwrap().1;
    assert_eq!(
        c.warn,
        vec!["No significant signaling interactions are inferred!"]
    );
    let t = table(&c.in_cols, &c.in_rows);
    let got = r_core::subset::subset_communication_deg(&t, &thresholds(c))
        .expect("an empty result after sources.use is a warning, not an error");
    assert!(got.rows.is_empty());
    assert_eq!(
        c.out_rn.as_deref(),
        Some(""),
        "a zero-row result gets no row names"
    );
    // ... while the same emptiness one stage earlier is an error.
    let c = &all.iter().find(|(n, _)| n == "err_empty_result").unwrap().1;
    let t = table(&c.in_cols, &c.in_rows);
    assert!(r_core::subset::subset_communication_deg(&t, &thresholds(c)).is_err());
}

/// `paste` coerces a missing cell with `as.character(NA)`, which is the **string** `"NA"`, so
/// a row with a missing `source` groups under the key `"NAsourceTotarget..."` -- and comes
/// back with `source` as that string, not as a missing cell. It is the one place a missing
/// value becomes data, and it is invisible unless a fixture has one.
#[test]
fn a_missing_source_becomes_the_string_na_as_a_group_key() {
    let all = load();
    let c = &all.iter().find(|(n, _)| n == "netp_plain").unwrap().1;
    let t = table(&c.in_cols, &c.in_rows);
    let got = r_core::subset::subset_communication_deg(&t, &thresholds(c)).unwrap();
    let first = &got.rows[0];
    assert_eq!(
        first[0].as_deref(),
        Some("NA"),
        "`source` for the missing-source group is the string \"NA\", not NA"
    );
    assert_eq!(first[1].as_deref(), Some("g1"));
    // And it sorts *first*: 'N' (0x4e) < 'g' (0x67), byte-wise.
    assert_eq!(got.rows.len(), 12, "no rows merged or lost");
}

/// `group_by` on character keys sorts byte-wise, so `"g10"` precedes `"g1"`. Sorting by
/// (source, target, pathway) as strings, or by group *level* order, gives a different
/// answer, and a fixture with three groups cannot tell the two apart.
#[test]
fn the_group_keys_are_ordered_byte_wise() {
    let all = load();
    let c = &all.iter().find(|(n, _)| n == "netp_plain").unwrap().1;
    let t = table(&c.in_cols, &c.in_rows);
    let got = r_core::subset::subset_communication_deg(&t, &thresholds(c)).unwrap();
    // The sort key is the **pasted** string, not the two halves. That is the whole point:
    // `"g1"` < `"g10"` as separate strings, but `"g10sourceTotargetg1"` < `"g1sourceTotargetg1"`
    // as pasted ones, because the second byte is '0' against 's'. So a port that sorts by
    // `(source, target)` gets a *different* order and is wrong, while looking right.
    let pasted = |r: &[Option<String>]| {
        format!(
            "{}sourceTotarget{}",
            r[0].as_ref().unwrap(),
            r[1].as_ref().unwrap()
        )
    };
    let mut want: Vec<(String, String)> = got
        .rows
        .iter()
        .map(|r| (pasted(r), r[2].clone().unwrap()))
        .collect();
    want.sort_by(|a, b| {
        a.0.as_bytes()
            .cmp(b.0.as_bytes())
            .then_with(|| a.1.as_bytes().cmp(b.1.as_bytes()))
    });
    let have: Vec<(String, String)> = got
        .rows
        .iter()
        .map(|r| (pasted(r), r[2].clone().unwrap()))
        .collect();
    assert_eq!(
        have, want,
        "group keys are not in byte-wise order of the pasted key"
    );

    let i10 = have.iter().position(|(k, _)| k.starts_with("g10")).unwrap();
    let i1 = have
        .iter()
        .position(|(k, _)| k.starts_with("g1source"))
        .unwrap();
    assert!(
        i10 < i1,
        "\"g10sourceTotarget...\" must sort before \"g1sourceTotarget...\""
    );
    // ... and `"g1"` alone still precedes `"g10"`, which is what makes this non-obvious.
    assert!("g1".as_bytes() < "g10".as_bytes());
}

/// The `NA`-blanking rule, isolated: with thresholds that exclude *nothing*, the rows that
/// disappear are exactly the ones with a missing value in an applied threshold column, and
/// a row that was entirely missing on input goes even with no threshold at all.
#[test]
fn a_missing_threshold_value_drops_the_row_via_blanking() {
    let all = load();
    let c = &all
        .iter()
        .find(|(n, _)| n == "all_eight_permissive")
        .unwrap()
        .1;
    let t = table(&c.in_cols, &c.in_rows);
    let got = r_core::subset::subset_communication_deg(&t, &thresholds(c)).unwrap();
    // 12 input rows, 4 of them carrying a missing value in a threshold column that the
    // permissive thresholds apply to: `NA` in ligand.pvalues (row 5), `NaN` in
    // ligand.pvalues (row 9), `NA` in ligand.logFC (row 7), `NA` in receptor.logFC (row 11).
    // None of the comparisons exclude anything, so these four are the whole difference.
    assert_eq!(
        got.rows.len(),
        6,
        "permissive thresholds: only the poisoned rows go"
    );

    // The same rule with nothing else in play: one threshold, set to the column maximum, so
    // no comparison excludes anything and the only rows lost are the two with a missing
    // `ligand.pvalues` -- one `NA` and one `NaN`, which are the same to R's `<=`.
    let c = &all
        .iter()
        .find(|(n, _)| n == "single_threshold_blanking")
        .unwrap()
        .1;
    let t = table(&c.in_cols, &c.in_rows);
    let got = r_core::subset::subset_communication_deg(&t, &thresholds(c)).unwrap();
    assert_eq!(
        got.rows.len(),
        10,
        "exactly the two rows with a missing ligand.pvalues go"
    );

    let c = &all
        .iter()
        .find(|(n, _)| n == "all_na_row_no_threshold")
        .unwrap()
        .1;
    let t = table(&c.in_cols, &c.in_rows);
    let got = r_core::subset::subset_communication_deg(&t, &thresholds(c)).unwrap();
    assert_eq!(
        got.rows.len(),
        11,
        "an all-NA row is dropped by rowSums(is.na(x)) != ncol(x) with no threshold set"
    );
}

/// `ligand.logFC` picks its comparison from the sign of the *argument*, so a negative
/// threshold keeps the downregulated genes and `0` takes the `>=` branch.
#[test]
fn the_sign_of_the_logfc_argument_chooses_the_comparison() {
    let all = load();
    let rows_of = |name: &str| -> Vec<Vec<Option<String>>> {
        let c = &all.iter().find(|(n, _)| n == name).unwrap().1;
        let t = table(&c.in_cols, &c.in_rows);
        r_core::subset::subset_communication_deg(&t, &thresholds(c))
            .unwrap()
            .rows
    };
    // The fixture's `ligand.logFC` is 0.5, -0.5, 1.5, 0, -2, 0.25, NA, 0.75, 0.1, 0, -1, 0.
    let pos = rows_of("ligand_logfc_pos");
    let neg = rows_of("ligand_logfc_neg");
    let zero = rows_of("ligand_logfc_zero");
    assert_eq!(pos.len(), 3, "logFC >= 0.5 keeps 0.5, 1.5, 0.75");
    assert_eq!(neg.len(), 3, "logFC <= -0.5 keeps -0.5, -2, -1");
    // `logFC = 0` takes the `>=` branch, so it keeps the non-negative values *including the
    // three exact zeros* -- and the row with an `NA` logFC goes, blanked.
    assert_eq!(zero.len(), 8, "logFC >= 0");
    assert_ne!(
        pos, neg,
        "the two branches select different rows, not just different counts"
    );
    assert!(pos.len() < zero.len() && neg.len() < zero.len());
}

/// The two `stop()` wordings are not the same string, and the order of the `if`s decides
/// which fires when two columns are absent at once.
#[test]
fn the_two_stop_messages_differ_and_the_first_test_wins() {
    let all = load();
    let err = |name: &str| -> String {
        let c = &all.iter().find(|(n, _)| n == name).unwrap().1;
        let t = table(&c.in_cols, &c.in_rows);
        r_core::subset::subset_communication_deg(&t, &thresholds(c))
            .unwrap_err()
            .to_string()
    };
    assert_eq!(
        err("err_ligand_pvalues_missing"),
        "Please run `identifyOverExpressedGenes` and `netMappingDEG` before using the threshold 'ligand.pvalues'"
    );
    assert_eq!(
        err("err_datasets_missing"),
        "Please run `identifyOverExpressedGenes` and `netMappingDEG` before selecting 'datasets'",
        "the datasets stop() is worded differently from the eight threshold ones"
    );
    assert_eq!(
        err("err_two_absent_order"),
        err("err_datasets_missing"),
        "datasets is tested first, so its message wins over three absent threshold columns"
    );
    assert_eq!(
        err("err_empty_result"),
        "No significant signaling interactions are inferred based on the input!"
    );
}

/// Upstream's final `if/else` reaches the `datasets` branch only when `ligand.logFC` is
/// *also* present, so a table carrying `datasets` and no `ligand.logFC` silently loses the
/// column. That looks like a bug; it is the code.
#[test]
fn the_datasets_column_needs_ligand_logfc_to_survive() {
    let all = load();
    let cols = |name: &str| -> Vec<String> {
        let c = all
            .iter()
            .find(|(n, _)| n == name)
            .unwrap()
            .1
            .out_cols
            .clone()
            .unwrap();
        c
    };
    assert!(
        cols("no_thresholds").contains(&"datasets".to_string()),
        "ligand.logFC + datasets keeps datasets"
    );
    assert!(
        !cols("no_thresholds_nodatasets").contains(&"datasets".to_string()),
        "datasets without ligand.logFC is dropped by the final select"
    );
}

/// An unrecognised `slot.name` makes both branches of upstream's final `if/else` false, so
/// the table comes back with every column it went in with.
#[test]
fn an_unknown_slot_name_keeps_every_column() {
    let all = load();
    let c = &all.iter().find(|(n, _)| n == "unknown_slot").unwrap().1;
    let t = table(&c.in_cols, &c.in_rows);
    let got = r_core::subset::subset_communication_deg(&t, &thresholds(c)).unwrap();
    assert_eq!(got.columns, c.in_cols, "no column selection at all");
}
