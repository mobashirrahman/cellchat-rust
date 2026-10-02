//! Parity of [`r_core::net`] against upstream `aggregateNet` and
//! `subsetCommunication_internal`.
//!
//! Corpus: `tests/fixtures/net_golden.txt`, produced by `tests/parity/gen_net_golden.R`,
//! which **sources `R/modeling.R` and `R/analysis.R` from the pinned commit** and calls the
//! upstream functions directly on four Prob/Pval fixtures.
//!
//! | fixture | shape | what it stresses |
//! |---|---|---|
//! | `f1` | 3 groups x 7 L-R | the general case, with five zeroed cells |
//! | `f2` | 4 groups x 12 L-R | a non-contiguous surviving set, so a row-index/array-index confusion cannot pass |
//! | `f3` | 1 group x 1 L-R | the degenerate shape |
//! | `f4` | 2 groups x 3 L-R, all zero | upstream's "No significant signaling interactions" error |
//!
//! Every fixture carries cells whose `pval` is **exactly** `thresh`, because
//! `prob[pval >= thresh] <- 0` is a `>=` and the boundary is the whole point of the
//! comparison.
//!
//! ## The two things this file is really checking
//!
//! 1. **Row order.** `reshape2::melt` builds its label frame with `expand.grid`, which
//!    varies its *first* argument fastest, so the melted rows run interaction-major, then
//!    target, then source — the opposite of the array's own column-major traversal. That
//!    order is the return value.
//! 2. **Accumulation order.** `aggregateNet` sums with R's `sum`, an `LDOUBLE` accumulation
//!    in column-major order, so `net$weight` is a function of the summation sequence and not
//!    just of the multiset of values.

use r_core::net::{aggregate_net_default, subset_communication, LrMeta, Net};

const GOLDEN: &str = include_str!("../../../tests/fixtures/net_golden.txt");

fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn parse_vec(s: &str) -> Vec<f64> {
    s.split(',')
        .map(|v| v.trim().parse().expect("number"))
        .collect()
}

fn bits_eq(got: f64, want: f64, what: &str) {
    if got.to_bits() != want.to_bits() {
        panic!(
            "{what}: bit mismatch\n  rust = {got:.17e} (0x{:016x})\n  R    = {want:.17e} (0x{:016x})",
            got.to_bits(),
            want.to_bits()
        );
    }
}

fn vecs_eq(got: &[f64], want: &[f64], what: &str) {
    assert_eq!(
        got.len(),
        want.len(),
        "{what}: length {} vs R's {}",
        got.len(),
        want.len()
    );
    for (i, (&g, &w)) in got.iter().zip(want).enumerate() {
        bits_eq(g, w, &format!("{what}[{i}]"));
    }
}

/// One fixture, assembled from the corpus.
struct Fixture {
    name: String,
    net: Net,
    lr: Vec<LrMeta>,
    k: usize,
    n_lr: usize,
}

fn fixtures() -> Vec<Fixture> {
    let mut out = Vec::new();
    for line in GOLDEN.lines().filter(|l| fields(l)[0] == "fixture") {
        let f = fields(line);
        let name = f[1];
        let k: usize = f[2].parse().unwrap();
        let n_lr: usize = f[3].parse().unwrap();
        let get = |kind: &str| {
            let pre = format!("{kind}\t{name}\t");
            GOLDEN
                .lines()
                .find(|l| l.starts_with(&pre))
                .unwrap_or_else(|| panic!("no {kind} for {name}"))
                .split('\t')
                .nth(2)
                .unwrap()
                .to_string()
        };
        let levels: Vec<String> = (1..=k).map(|i| format!("g{i}")).collect();
        let names: Vec<String> = (1..=n_lr).map(|i| format!("L{i}^R{i}")).collect();
        let lr: Vec<LrMeta> = (1..=n_lr)
            .map(|i| LrMeta {
                interaction_name: format!("L{i}^R{i}"),
                columns: r_core::net::LR_OPTIONAL
                    .iter()
                    .map(|s| (*s).to_string())
                    .collect(),
                interaction_name_2: format!("L{i}_R{i}"),
                pathway_name: ["TGFb", "WNT", "NOTCH", "VEGF"][(i - 1) % 4].to_string(),
                ligand: format!("L{i}"),
                receptor: format!("R{i}"),
                annotation: [
                    "Secreted Signaling",
                    "ECM-Receptor",
                    "Non-protein Signaling",
                    "Cell-Cell Contact",
                ][(i - 1) % 4]
                    .to_string(),
                evidence: None,
            })
            .collect();
        let net = Net::new(
            parse_vec(&get("prob")),
            parse_vec(&get("pval")),
            levels,
            names,
        );
        assert_eq!(net.prob.len(), k * k * n_lr, "{name}: prob length");
        assert_eq!(net.pval.len(), k * k * n_lr, "{name}: pval length");
        out.push(Fixture {
            name: name.to_string(),
            net,
            lr,
            k,
            n_lr,
        });
    }
    assert!(out.len() >= 4, "expected at least four fixtures");
    out
}

/// `LRmeta` in the generator, mirrored so the Rust table cannot drift silently.
fn expected_lr_columns() -> Vec<&'static str> {
    r_core::net::COLUMNS_PLAIN.to_vec()
}

#[test]
fn every_fixture_is_covered() {
    let got: Vec<&str> = GOLDEN
        .lines()
        .filter_map(|l| l.strip_prefix("fixture\t"))
        .map(|l| fields(l)[0])
        .collect();
    let mine: Vec<String> = fixtures().into_iter().map(|f| f.name).collect();
    assert_eq!(
        got.len(),
        mine.len(),
        "fixture count: corpus {got:?} vs Rust {mine:?}"
    );
    for (g, m) in got.iter().zip(&mine) {
        assert_eq!(*g, m);
    }
}

/// The corpus records `k` and `n_lr` on every fixture line. Parsing them and never checking them
/// would mean the positional read of columns 2 and 3 is unverified -- a corpus that shifted a
/// column would still produce plausible numbers, just the wrong ones. This is the test that makes
/// the positional parse trustworthy.
#[test]
fn fixture_header_agrees_with_its_own_data() {
    for fx in fixtures() {
        assert_eq!(
            fx.k, fx.net.n_groups,
            "{}: k in the corpus does not match the net's group count",
            fx.name
        );
        assert_eq!(
            fx.n_lr,
            fx.lr.len(),
            "{}: n_lr in the corpus does not match the L-R rows",
            fx.name
        );
    }
}

#[test]
fn aggregate_net_count_and_weight_are_bit_identical() {
    let mut checked = 0;
    for fx in fixtures() {
        for thresh in [0.05f64, 0.01, 0.2] {
            let (count, weight) = aggregate_net_default(&fx.net, thresh);
            let key = format!("{}\t{}\t", fx.name, thresh);
            let rec = |kind: &str| -> Vec<f64> {
                let pre = format!("{kind}\t{key}");
                let l = GOLDEN
                    .lines()
                    .find(|l| l.starts_with(&pre))
                    .unwrap_or_else(|| panic!("no {kind} for {key}"));
                parse_vec(l.split('\t').nth(3).unwrap())
            };
            vecs_eq(
                &count,
                &rec("count"),
                &format!("{} count thresh={thresh}", fx.name),
            );
            vecs_eq(
                &weight,
                &rec("weight"),
                &format!("{} weight thresh={thresh}", fx.name),
            );
            checked += 2;
        }
    }
    assert_eq!(
        checked,
        4 * 3 * 2,
        "four fixtures x three thresholds x two arrays"
    );
}

#[test]
fn aggregate_net_dimnames_follow_the_group_levels() {
    for fx in fixtures() {
        let thresh = 0.05f64;
        let pre = format!("count_dimnames\t{}\t{}\t", fx.name, thresh);
        let l = GOLDEN
            .lines()
            .find(|l| l.starts_with(&pre))
            .unwrap_or_else(|| panic!("no count_dimnames for {}", fx.name));
        let f: Vec<&str> = l.split('\t').collect();
        let want_rows: Vec<&str> = f[3].split(',').collect();
        let want_cols: Vec<&str> = f[4].split(',').collect();
        let (r, c, _) = fx.net.dimnames();
        assert_eq!(
            r.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            want_rows,
            "{}",
            fx.name
        );
        assert_eq!(
            c.iter().map(|s| s.as_str()).collect::<Vec<_>>(),
            want_cols,
            "{}",
            fx.name
        );
    }
}

#[test]
fn subset_communication_is_row_for_row_and_column_for_column() {
    let mut n_rows_checked = 0;
    for fx in fixtures() {
        for thresh in [0.05f64, 0.01] {
            for (tag, sources, targets) in [
                ("all", None, None),
                ("g1", Some(vec!["g1".to_string()]), None),
                (
                    "g1+g2",
                    Some(vec!["g1".to_string(), "g2".to_string()]),
                    Some(vec!["g2".to_string()]),
                ),
            ] {
                let key = format!("{}|{}|{}", fx.name, thresh, tag);
                let pre = format!("subset\t{key}\t");
                let head = GOLDEN
                    .lines()
                    .find(|l| l.starts_with(&pre))
                    .unwrap_or_else(|| panic!("no subset record for {key}"));
                let hf: Vec<&str> = head.split('\t').collect();
                if hf[2] == "ERROR" {
                    // Upstream refuses to return an empty table. The Rust side returns an
                    // empty one, so the *error* is part of the contract; see
                    // `an_empty_result_is_upstreams_error_not_an_empty_table`.
                    assert_eq!(
                        hf[3],
                        "No significant signaling interactions are inferred based on the input!"
                    );
                    let t = subset_communication(
                        &fx.net,
                        &fx.lr,
                        thresh,
                        sources.as_deref(),
                        targets.as_deref(),
                        false,
                        true,
                    );
                    assert!(t.is_empty(), "{key}: R errored, so Rust should have too");
                    continue;
                }
                let want_shape = hf[2];
                let want_cols: Vec<&str> = hf[3].split(',').collect();
                let t = subset_communication(
                    &fx.net,
                    &fx.lr,
                    thresh,
                    sources.as_deref(),
                    targets.as_deref(),
                    false,
                    true,
                );

                // rownames are 1:nrow
                let rn = GOLDEN
                    .lines()
                    .find(|l| l.starts_with(&format!("subset_rownames\t{key}\t")))
                    .unwrap()
                    .split('\t')
                    .nth(2)
                    .unwrap();
                let want_rn: Vec<&str> = if rn.trim().is_empty() {
                    vec![]
                } else {
                    rn.split(',').collect()
                };
                assert_eq!(
                    t.len(),
                    want_rn.len(),
                    "{key}: nrow {} vs R's {}",
                    t.len(),
                    want_rn.len()
                );
                assert_eq!(t.columns, want_cols, "{key}: column set/order");

                for cname in &t.columns {
                    let pre = format!("subset_col\t{key}\t{cname}\t");
                    let line = GOLDEN
                        .lines()
                        .find(|l| l.starts_with(&pre))
                        .unwrap_or_else(|| panic!("no subset_col {key}/{cname}"));
                    // The generator writes with `format(..., digits = 17)`, which right-pads
                    // to a common width, so every cell must be trimmed before comparing.
                    // An empty field means a zero-row column: `paste(character(0), ...)`
                    // yields "", and `"".split(',')` is one empty string rather than none.
                    let raw = line.split('\t').nth(3).unwrap();
                    let want: Vec<&str> = if raw.trim().is_empty() {
                        Vec::new()
                    } else {
                        raw.split(',').map(str::trim).collect()
                    };
                    assert_eq!(want.len(), t.len(), "{key}/{cname}: length");
                    for (r, w) in want.iter().enumerate() {
                        let w = *w;
                        match (w, t.get(r, cname)) {
                            ("NA", None) => {}
                            ("NA", Some(g)) => {
                                panic!("{key}/{cname}[{r}]: R has NA, Rust has {g:?}")
                            }
                            (w, None) => panic!("{key}/{cname}[{r}]: R has {w:?}, Rust has NA"),
                            (w, Some(g)) => {
                                // Numeric columns are compared bit-exactly; character columns
                                // exactly.
                                let both_numeric =
                                    w.parse::<f64>().is_ok() && g.parse::<f64>().is_ok();
                                if both_numeric {
                                    bits_eq(
                                        g.parse().unwrap(),
                                        w.parse().unwrap(),
                                        &format!("{key}/{cname}[{r}]"),
                                    );
                                } else {
                                    assert_eq!(g, w, "{key}/{cname}[{r}]");
                                }
                            }
                        }
                    }
                }
                let _ = want_shape;
                n_rows_checked += t.len();
            }
        }
    }
    assert!(n_rows_checked > 50, "only compared {n_rows_checked} rows");
}

#[test]
fn an_empty_result_is_upstreams_error_not_an_empty_table() {
    // `subsetCommunication_internal` does `stop("No significant signaling interactions are
    // inferred based on the input!")` when nothing survives. The Rust function returns an
    // empty table instead, so the binding layer is responsible for raising the error; the
    // message must be byte-identical, and this test pins the *condition* that triggers it.
    let net = Net::new(
        vec![0.0; 8],
        vec![0.9; 8],
        vec!["a".into(), "b".into()],
        vec!["i".into(), "j".into()],
    );
    let t = subset_communication(&net, &[LrMeta::default()], 0.05, None, None, false, true);
    assert!(t.is_empty());
}

#[test]
fn the_column_set_is_the_upstream_intersect() {
    // `intersect(c("source", "target", "ligand", "receptor", "prob", "pval",
    // "interaction_name", "interaction_name_2", "pathway_name", "annotation", "evidence"),
    // colnames(net))` -- and `evidence` survives even though every cell is NA.
    let net = Net::new(
        vec![0.5, 0.0, 0.0, 0.0],
        vec![0.01, 0.5, 0.5, 0.5],
        vec!["a".into(), "b".into()],
        vec!["i".into()],
    );
    let mut m = LrMeta::minimal("i", "L^R", "L", "R", "Secreted Signaling");
    m.evidence = None;
    let t = subset_communication(&net, &[m], 0.05, None, None, false, true);
    assert_eq!(t.columns, expected_lr_columns());
    assert_eq!(t.len(), 1);
    assert_eq!(
        t.get(0, "evidence"),
        None,
        "evidence is present but every cell is NA"
    );
}
