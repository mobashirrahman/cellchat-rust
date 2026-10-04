//! Parity of [`r_core::de`] against upstream `computeAveExpr`, `subsetDB` and
//! `subsetData`.
//!
//! Corpus: `tests/fixtures/de_golden.txt`, produced by `tests/parity/gen_de_golden.R`,
//! which **sources `modeling.R`, `utilities.R` and `database.R` from the pinned commit**.
//!
//! ## What this pins that a value-only test would miss
//!
//! * **`computeAveExpr`'s row order is the caller's, not the matrix's.**
//!   `features.use <- intersect(features, rownames(data.use))`, and R's `intersect`
//!   returns its *first* argument's values in first-appearance order, deduplicated. The
//!   `featsel` fixture asks for `c("G10","G1","G10","GZZ","G2")` and gets rows
//!   `G10, G1, G2` — a different order from the matrix's `G1, G2, G10` and with the
//!   duplicate collapsed.
//! * **`subsetData` reorders the DB's annotation column**, and the human DB is *not* already
//!   in factor order, so the reordering is a real change rather than a no-op. The fixture
//!   records the annotation sequence before and after.
//! * **`match.arg` and the `key` validation messages**, byte for byte.

use r_core::db::Database;
use r_core::de::{
    annotation_rank, compute_ave_expr, intersect_order, match_ave_expr_type, order_by_annotation,
    search_implies_non_protein, subset_data_gene_use, subset_db_default_search, KernelErrorLite,
};
use std::collections::HashMap;
use std::path::PathBuf;

const GOLDEN: &str = include_str!("../../../../../tests/fixtures/de_golden.txt");

const N_GENES: usize = 24;
const N_CELLS: usize = 60;
const N_GROUPS: usize = 3;

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tests/fixtures")
}

fn human_db() -> Database {
    Database::load(&fixtures().join("db_human")).expect("db_human fixture")
}

fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

/// NA-aware: `median` of an all-missing row is `NA`, and the generator writes it as `NA`
/// rather than as a number, because `sprintf("%.17g", NA_real_)` gives `"NA"`.
fn parse_vec(s: &str) -> Vec<f64> {
    s.split(',')
        .map(|v| {
            let v = v.trim();
            if v == "NA" || v == "NaN" {
                f64::NAN
            } else {
                v.parse().unwrap_or_else(|_| panic!("not a number: {v:?}"))
            }
        })
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

/// `record(kind)` -> the value field, tolerating both 2- and 3-field lines
/// (`sprintf("%s\t%s\n", ...)` vs `sprintf("%s\t%s\t%s\n", ...)`).
fn golden(kind: &str) -> Option<String> {
    let pre = format!("{kind}\t");
    GOLDEN
        .lines()
        .find(|l| l.starts_with(&pre))
        .map(|l| fields(l).last().unwrap().to_string())
}

fn golden_named(kind: &str, name: &str) -> Option<Vec<String>> {
    let pre = format!("{kind}\t{name}\t");
    GOLDEN
        .lines()
        .find(|l| l.starts_with(&pre))
        .map(|l| fields(l).last().unwrap().to_string())
        .map(|v| v.split(',').map(|s| s.to_string()).collect())
}

/// The fixture matrix, regenerated here from the generator's recipe.
///
/// `m[1,] <- 0`, `m[2, 1:30] <- 0`, `m[3,] <- 2`, `m[4, seq(1,60,by=3)] <- 0`,
/// `m[5,] <- NA`, and the rest `runif(0, 5)` from seed 20240303. The generator's
/// `set.seed`/`runif` is *not* re-derived -- see the note in `tests/parity/README.md` --
/// so these constants are only used to check the *shape* and the engineered rows, and the
/// values themselves come from `golden("ave_vals")`.
///
/// Actually the values are needed, so: the matrix is read back from the corpus by
/// inverting nothing -- the corpus stores the *output*, and the input matrix is rebuilt
/// from the same R recipe here. That is the one place a re-derivation is unavoidable,
/// because the input is 1440 numbers and the generator does not dump it. To keep it honest
/// the engineered rows are asserted explicitly.
fn fixture_matrix() -> (Vec<f64>, Vec<String>) {
    let mut mt = r_core::rng::MersenneTwister::new(20240303);
    let genes: Vec<String> = (1..=N_GENES).map(|i| format!("G{i}")).collect();
    // genes x cells, column-major, from runif(0, 5) -> 0.01 + 4.99*u is NOT the same
    // distribution, so scale exactly: runif(a, b) = a + (b - a) * unif_rand().
    let mut m: Vec<f64> = (0..N_GENES * N_CELLS)
        .map(|_| 0.0 + 5.0 * mt.unif_rand())
        .collect();
    for c in 0..N_CELLS {
        m[c * N_GENES] = 0.0;
    }
    for c in 0..30 {
        m[1 + c * N_GENES] = 0.0;
    }
    for c in 0..N_CELLS {
        m[2 + c * N_GENES] = 2.0;
    }
    for c in (0..N_CELLS).step_by(3) {
        m[3 + c * N_GENES] = 0.0;
    }
    for c in 0..N_CELLS {
        m[4 + c * N_GENES] = f64::NAN;
    }
    (m, genes)
}

/// `t(m)`: cells x genes, column-major.
fn cells_by_genes(m: &[f64], genes: &[String]) -> (Vec<f64>, Vec<String>) {
    let mut out = vec![0.0; m.len()];
    for g in 0..N_GENES {
        for c in 0..N_CELLS {
            out[g * N_CELLS + c] = m[g + c * N_GENES];
        }
    }
    (out, genes.to_vec())
}

/// The group assignment. The generator *permutes* it (`sample.int`), so the labels per cell
/// are recovered from the corpus instead: the corpus does not store them, so this
/// reconstructs the same permutation from the same RNG.
fn fixture_group() -> Vec<usize> {
    let base = [0usize; 25];
    let mut labels: Vec<usize> = Vec::with_capacity(N_CELLS);
    labels.extend(base.iter());
    labels.extend(std::iter::repeat_n(1, 20));
    labels.extend(std::iter::repeat_n(2, 15));
    // `c(rep(..), ..)[sample.int(60)]`
    let mut mt = r_core::rng::MersenneTwister::new(20240303);
    // The generator's RNG stream: `runif(24*60, 0, 5)` for the matrix, then the permutation.
    // The matrix draws are already consumed by `fixture_matrix`, so this is a fresh stream
    // offset by exactly that many draws.
    for _ in 0..N_GENES * N_CELLS {
        mt.unif_rand();
    }
    let p = mt.sample_int_permutation(N_CELLS);
    p.into_iter().map(|c| labels[c as usize - 1]).collect()
}

struct AveSpec {
    name: &'static str,
    type_: &'static str,
    trim: f64,
    features: Option<Vec<String>>,
}

fn ave_specs() -> Vec<AveSpec> {
    let f = |v: &[&str]| Some(v.iter().map(|s| s.to_string()).collect());
    vec![
        AveSpec {
            name: "tri",
            type_: "triMean",
            trim: 0.0,
            features: None,
        },
        AveSpec {
            name: "trim01",
            type_: "truncatedMean",
            trim: 0.1,
            features: None,
        },
        AveSpec {
            name: "trim05",
            type_: "truncatedMean",
            trim: 0.5,
            features: None,
        },
        AveSpec {
            name: "median",
            type_: "median",
            trim: 0.0,
            features: None,
        },
        AveSpec {
            name: "featsel",
            type_: "triMean",
            trim: 0.0,
            features: f(&["G10", "G1", "G10", "GZZ", "G2"]),
        },
        AveSpec {
            name: "matcharg",
            type_: "tri",
            trim: 0.0,
            features: None,
        },
        AveSpec {
            name: "median_missing",
            type_: "median",
            trim: 0.0,
            features: f(&["G5", "G1", "G4"]),
        },
    ]
}

#[test]
fn the_fixture_matrix_has_the_engineered_rows() {
    let (m, genes) = fixture_matrix();
    assert_eq!(genes.len(), N_GENES);
    assert_eq!(m.len(), N_GENES * N_CELLS);
    for c in 0..N_CELLS {
        assert_eq!(m[c * N_GENES], 0.0, "gene 1 is all zero");
        assert_eq!(m[2 + c * N_GENES], 2.0, "gene 3 is constant");
        assert!(m[4 + c * N_GENES].is_nan(), "gene 5 is all NA");
    }
    for c in 30..N_CELLS {
        assert_ne!(
            m[1 + c * N_GENES],
            0.0,
            "gene 2 is zero only in the first half"
        );
    }
    for c in 0..30 {
        assert_eq!(m[1 + c * N_GENES], 0.0);
    }
    for c in (0..N_CELLS).step_by(3) {
        assert_eq!(m[3 + c * N_GENES], 0.0, "gene 4 is sparse");
    }
}

#[test]
fn every_ave_spec_is_covered() {
    let got: Vec<String> = GOLDEN
        .lines()
        .filter(|l| l.starts_with("ave\t"))
        .map(|l| fields(l)[1].to_string())
        .collect();
    let mine: Vec<String> = ave_specs()
        .into_iter()
        .map(|s| s.name.to_string())
        .collect();
    assert_eq!(
        got, mine,
        "the Rust spec list has drifted from the generator's"
    );
}

#[test]
fn compute_ave_expr_is_bit_identical_including_row_order() {
    let (m, genes) = fixture_matrix();
    let (data, genes) = cells_by_genes(&m, &genes);
    let group = fixture_group();
    for spec in ave_specs() {
        let want_rows = golden_named("ave_rows", spec.name).expect("ave_rows record");
        let want_cols = golden_named("ave_cols", spec.name).expect("ave_cols record");
        let want_vals: Vec<f64> = parse_vec(
            &golden_named("ave_vals", spec.name)
                .expect("ave_vals record")
                .join(","),
        );

        // The row order is `intersect(features, rownames)`, so check the *names* the Rust
        // side would produce against R's, by comparing the value vectors positionally.
        let got = compute_ave_expr(
            &data,
            &genes,
            &group,
            N_GROUPS,
            spec.features.as_deref(),
            spec.type_,
            spec.trim,
        )
        .unwrap_or_else(|e| panic!("{}: {e}", spec.name));

        assert_eq!(want_cols.len(), N_GROUPS, "{}: ncol", spec.name);
        assert_eq!(got.len(), want_vals.len(), "{}: length", spec.name);
        for (i, (&g, &w)) in got.iter().zip(&want_vals).enumerate() {
            // NaN (gene 5 under `median`) is compared as NaN, not by bits.
            if g.is_nan() || w.is_nan() {
                assert!(
                    g.is_nan() && w.is_nan(),
                    "{}: cell {i}: rust={g} r={w}",
                    spec.name
                );
            } else {
                bits_eq(g, w, &format!("{}: [{i}]", spec.name));
            }
        }
        // And the row order, which the positional comparison above relies on.
        let want: Vec<String> = match &spec.features {
            None => genes.clone(),
            Some(f) => intersect_order(f, &genes),
        };
        assert_eq!(want, want_rows, "{}: row order", spec.name);
    }
}

#[test]
fn intersect_order_matches_r() {
    let a: Vec<String> = ["G10", "G1", "G10", "GZZ", "G2"]
        .iter()
        .map(|s| s.to_string())
        .collect();
    let genes: Vec<String> = (1..=N_GENES).map(|i| format!("G{i}")).collect();
    assert_eq!(intersect_order(&a, &genes), vec!["G10", "G1", "G2"]);
}

#[test]
fn match_arg_matches_r() {
    let want = golden("ave_argerr").expect("ave_argerr record");
    let e = match_ave_expr_type("t", 0.0).unwrap_err();
    assert_eq!(e.to_string(), want, "byte for byte");
    assert!(matches!(e, KernelErrorLite::MatchArg(_)));
    // The partial matches that must succeed.
    assert!(match_ave_expr_type("tri", 0.0).is_ok());
    assert!(match_ave_expr_type("triMean", 0.0).is_ok());
    assert!(match_ave_expr_type("truncated", 0.0).is_ok());
}

#[test]
fn subset_data_reorders_the_annotation_column() {
    let db = human_db();
    let before: Vec<String> = db
        .interactions
        .iter()
        .map(|i| i.annotation.clone())
        .collect();
    let want_after: Vec<String> = golden_named("subset_db_ann", "default").expect("record");
    assert_eq!(want_after.len(), before.len());
    // `r-core`'s DB is *already* in annotation order -- `export_db.R` sorts it so the port
    // never has to reproduce `order()` (PLAN.md, Phase 1d). So `subsetData`'s reordering is
    // a no-op on the Rust side, and the meaningful invariant is the one below: r-core's
    // interaction order equals upstream's order *after* `subsetData`.
    let ranks: Vec<usize> = before
        .iter()
        .map(|a| annotation_rank(a).expect("known"))
        .collect();
    assert!(
        ranks.windows(2).all(|w| w[0] <= w[1]),
        "the TSV export is supposed to be pre-sorted in annotation order"
    );
    assert_eq!(
        before, want_after,
        "r-core's order must equal upstream's post-subsetData order"
    );
    // `order_by_annotation` is still exercised directly, on an unsorted input.
    let mut shuffled = before.clone();
    shuffled.reverse();
    let order = order_by_annotation(&shuffled);
    let got: Vec<String> = order.iter().map(|&i| shuffled[i].clone()).collect();
    assert_eq!(got, want_after, "order_by_annotation on a reversed input");
    assert_eq!(order.len(), before.len(), "no row is dropped");
}

#[test]
fn subset_data_gene_use_is_intersect_of_the_db_list_with_the_matrix() {
    let db = human_db();
    let gene_use_input = db.extract_gene();
    // A matrix rownames set that is a strict subset, in a scrambled order.
    let rownames: Vec<String> = (1..=8).rev().map(|i| format!("G{i}")).collect();
    let got = subset_data_gene_use(&gene_use_input, &rownames, None);
    // R's order: the DB's gene order, filtered.
    let want: Vec<String> = gene_use_input
        .iter()
        .filter(|g| rownames.contains(g))
        .cloned()
        .collect();
    assert_eq!(got, want);
    // An explicit feature list reorders, because `intersect` keeps the first argument.
    let feats: Vec<String> = vec!["G5".into(), "G1".into(), "G5".into(), "G9".into()];
    let got2 = subset_data_gene_use(&gene_use_input, &rownames, Some(&feats));
    assert_eq!(got2, vec!["G5", "G1"]);
}

#[test]
fn subset_db_filters_by_annotation() {
    let db = human_db();
    let total = db.interactions.len();
    let count = |ann: &str| {
        db.interactions
            .iter()
            .filter(|i| i.annotation == ann)
            .count()
    };

    // Default: everything except Non-protein Signaling.
    let def = subset_db_default_search(false);
    let want_default = total - count("Non-protein Signaling");
    assert_eq!(
        db.interactions
            .iter()
            .filter(|i| def.contains(&i.annotation.as_str()))
            .count(),
        want_default
    );
    // With `non_protein = TRUE` the fourth annotation joins the search list.
    let all = subset_db_default_search(true);
    assert!(all.contains(&"Non-protein Signaling"));
    assert_eq!(
        db.interactions
            .iter()
            .filter(|i| all.contains(&i.annotation.as_str()))
            .count(),
        total
    );
    // Naming it in `search` flips `non_protein` even when the flag is FALSE.
    let search: Vec<String> = subset_db_default_search(false)
        .iter()
        .map(|s| s.to_string())
        .collect();
    assert!(!search_implies_non_protein(&search));
    let mut with_np = search.clone();
    with_np.push("Non-protein Signaling".into());
    assert!(search_implies_non_protein(&with_np));
}

#[test]
fn subset_db_matches_the_corpus_row_counts() {
    let db = human_db();
    let total = db.interactions.len();
    let cases: Vec<(&str, Option<&[&str]>, bool)> = vec![
        ("default", None, false),
        ("nonprotein", None, true),
        (
            "explicit_all",
            Some(&[
                "Secreted Signaling",
                "ECM-Receptor",
                "Cell-Cell Contact",
                "Non-protein Signaling",
            ]),
            false,
        ),
        ("contact_only", Some(&["Cell-Cell Contact"]), false),
        ("empty_search", Some(&[]), false),
    ];
    for (name, search, non_protein) in cases {
        let line = GOLDEN
            .lines()
            .find(|l| l.starts_with(&format!("subset_db\t{name}\t")))
            .unwrap_or_else(|| panic!("no subset_db record for {name}"));
        let want_rows: usize = fields(line)[2].parse().expect("row count");
        let search: Vec<String> = match search {
            None => subset_db_default_search(non_protein)
                .iter()
                .map(|s| s.to_string())
                .collect(),
            Some(v) => v.iter().map(|s| s.to_string()).collect(),
        };
        // Naming "Non-protein Signaling" in `search` *flips* `non_protein` to TRUE, upstream:
        //
        //     if ("Non-protein Signaling" %in% unlist(search)) { non_protein = TRUE; ... }
        //
        // which is why `explicit_all` keeps 3233 rows even though its `non_protein` argument
        // is FALSE. Missing this is a 994-row difference.
        let effective_np = non_protein || search_implies_non_protein(&search);
        let got = db
            .interactions
            .iter()
            .filter(|i| effective_np || i.annotation != "Non-protein Signaling")
            .filter(|i| search.iter().any(|s| s == &i.annotation))
            .count();
        assert_eq!(got, want_rows, "{name}: interaction count");
        let _ = total;
    }
}

#[test]
fn the_db_key_error_is_reproduced() {
    let want = golden("subset_db_keyerr").expect("subset_db_keyerr record");
    let e = KernelErrorLite::UnknownKey("nope".into());
    assert_eq!(e.to_string(), want, "byte for byte");
}

#[test]
fn the_corpus_has_records_for_everything_the_test_reads() {
    for kind in [
        "ave_argerr",
        "subset_db_keyerr",
        "db_int_ann0",
        "db_int_rows",
    ] {
        assert!(golden(kind).is_some(), "no `{kind}` record in the corpus");
    }
    for name in ["default", "explicit", "empty"] {
        assert!(
            golden_named("subset_data_rows", name).is_some(),
            "no subset_data_rows for {name}"
        );
    }
    // `db_int_ann0` is the *pre*-`subsetData` sequence straight out of the pinned .rda, and
    // `subset_db_ann default` is the same table after `subsetData` has reordered it. The
    // export matches the second, not the first -- and the two must differ, or the
    // reordering is untested.
    let pre: Vec<String> = golden("db_int_ann0")
        .unwrap()
        .split(',')
        .map(|s| s.to_string())
        .collect();
    let after: Vec<String> = golden_named("subset_db_ann", "default").unwrap();
    assert_eq!(pre.len(), after.len());
    assert_ne!(
        pre, after,
        "the fixture would be vacuous if subsetData changed nothing"
    );
    let db = human_db();
    assert_eq!(after.len(), db.interactions.len());
    for (i, x) in after.iter().enumerate() {
        assert_eq!(x, &db.interactions[i].annotation, "row {i}");
    }
}

#[test]
fn group_levels_are_derived_not_assumed() {
    // The group vector is reconstructed from the same RNG recipe, so assert its shape and
    // that every level is populated -- otherwise the aggregate is vacuous.
    let g = fixture_group();
    assert_eq!(g.len(), N_CELLS);
    for lvl in 0..N_GROUPS {
        assert!(g.contains(&lvl), "level {lvl} is empty");
    }
    let mut counts = HashMap::new();
    for &x in &g {
        *counts.entry(x).or_insert(0usize) += 1;
    }
    assert_eq!(counts.len(), N_GROUPS);
}
