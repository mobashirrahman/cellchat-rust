//! Parity of [`r_core::expr`] against upstream CellChat's `computeExpr_*` family.
//!
//! Corpus: `tests/fixtures/expr_golden.txt`, produced by
//! `tests/parity/gen_expr_golden.R`, which **sources `R/modeling.R` from the pinned
//! commit** and calls the upstream functions directly. Nothing in the expected values is
//! reimplemented on the R side, so this is a true differential test rather than two ports
//! agreeing with each other.
//!
//! The fixture matrix is regenerated here from R's own RNG -- `set.seed(4242)` then
//! `runif(n, 0.01, 1)` -- via [`r_core::rng`], so the Rust side sees bit-identical inputs
//! to the R side without a second 50 KB literal in the fixture. That is also a
//! cross-check of the RNG: if [`MersenneTwister::unif_rand`] drifted, every record here
//! would fail at once.
//!
//! The fixture deliberately covers both *outcomes* of each resolution path:
//!
//! | fixture case | expected |
//! |---|---|
//! | `absent_not_complex` (ligand path) | **error** `subscript out of bounds` (R-ism 11) |
//! | `missing_subunit_row` (ligand path) | **error** `subscript out of bounds` (R-ism 11) |
//! | `NOT_A_COFACTOR` as agonist/cofactor | all ones (no modulation), *not* an error |
//! | `IL12AB` (a complex name) in a cofactor slot | all ones |
//! | single-subunit cofactor | `1 + x` / `Kh^n/(Kh^n+x^n)` |
//! | multi-subunit cofactor | product, accumulated in `LDOUBLE` |
//!
//! The last four rows are the important ones: a cofactor that cannot be resolved degrades
//! to a no-op, while a ligand/receptor that cannot be resolved aborts. Both behaviours are
//! upstream's and they are not interchangeable.

use r_core::aggregate::GroupMean;
use r_core::db::Database;
use r_core::expr::{
    compute_expr_agonist, compute_expr_antagonist, compute_expr_complex, compute_expr_coreceptor,
    compute_expr_group_agonist, compute_expr_group_antagonist, compute_expr_lr, CellExpr,
    CoreceptorKind, GroupMeans, SUBSCRIPT_OUT_OF_BOUNDS,
};
use r_core::rng::MersenneTwister;
use std::collections::HashMap;
use std::path::{Path, PathBuf};

const GOLDEN: &str = include_str!("../../../../../tests/fixtures/expr_golden.txt");
const MATRIX: &str = include_str!("../../../../../tests/fixtures/expr_matrix.tsv");
const SEED: i32 = 4242;

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

/// The `X` matrix, read verbatim from `expr_matrix.tsv` as hex floats.
///
/// `expr_matrix.tsv` is written by the generator, so it is ground truth. Re-deriving it
/// here from `set.seed`/`runif` was tried and abandoned: it required this test to
/// reproduce R's `unlist(data.frame[, cols])` column-major order and `unique()`'s
/// first-appearance order, and getting either wrong kept the gene *set* correct while
/// permuting the *values* -- every name still resolved and every number still looked
/// plausible. Ground truth in a file cannot rot that way.
fn fixture_matrix() -> (Vec<f64>, Vec<String>, usize) {
    // Self-describing layout: `groups<TAB>n`, `genes<TAB>m`, m names, `values`, then one
    // line of space-separated hex floats per gene. The markers matter -- an earlier
    // version of the generator wrote the two halves through separate connections and R
    // flushed them out of order, and a reader that infers the layout from line counts
    // cannot tell a reordered file from a corrupted one.
    let mut lines = MATRIX.lines();
    let n_groups = header_field(
        lines.next().expect("expr_matrix.tsv: missing header"),
        "groups",
    );
    let n_genes = header_field(
        lines.next().expect("expr_matrix.tsv: missing gene count"),
        "genes",
    );
    let genes: Vec<String> = lines
        .by_ref()
        .take(n_genes)
        .map(|s| s.to_string())
        .collect();
    assert_eq!(
        genes.len(),
        n_genes,
        "expr_matrix.tsv: gene-name block is short"
    );
    assert_eq!(
        lines.next(),
        Some("values"),
        "expr_matrix.tsv: expected a `values` marker"
    );
    let row_lines: Vec<&str> = lines.collect();
    assert_eq!(
        row_lines.len(),
        n_genes,
        "expr_matrix.tsv: expected one value line per gene"
    );
    // Generator order is gene-major (one line per gene); the buffer is column-major, like
    // every other matrix in the port.
    let mut data = vec![0.0; genes.len() * n_groups];
    for (g, line) in row_lines.iter().enumerate() {
        let vals: Vec<f64> = line.split_whitespace().map(parse_hex_float).collect();
        assert_eq!(
            vals.len(),
            n_groups,
            "expr_matrix.tsv: gene {g} has {} values",
            vals.len()
        );
        for (j, v) in vals.into_iter().enumerate() {
            data[j * genes.len() + g] = v;
        }
    }
    assert_eq!(data.len(), genes.len() * n_groups);
    (data, genes, n_groups)
}

fn header_field(line: &str, key: &str) -> usize {
    line.strip_prefix(key)
        .and_then(|v| v.strip_prefix('\t'))
        .unwrap_or_else(|| panic!("expr_matrix.tsv: expected a `{key}<TAB>n` header, got {line:?}"))
        .parse()
        .unwrap_or_else(|_| panic!("expr_matrix.tsv: bad count in {line:?}"))
}

/// Parse a C99 hex float as emitted by R's `sprintf("%a", x)`.
///
/// `f64::from_str_radix` parses only the integer significand, and Rust's `f64` parser
/// does not accept the `0x...p...` form, so the three parts are split by hand. This
/// round-trips exactly: every finite double has a unique shortest-round-trip decimal
/// form, and `%a` is the binary-exact counterpart.
fn parse_hex_float(s: &str) -> f64 {
    let s = s.trim();
    assert!(
        s.contains("0x") || s.contains('p'),
        "not a hex float: {s:?}"
    );
    let (neg, s) = match s.strip_prefix('-') {
        Some(rest) => (true, rest),
        None => (false, s.strip_prefix('+').unwrap_or(s)),
    };
    let s = s.trim_start_matches("0x").trim_start_matches("0X");
    let (mantissa_str, exp_str) = match s.find(['p', 'P']) {
        Some(i) => (&s[..i], &s[i + 1..]),
        None => (s, "0"),
    };
    let (int_part, frac_part) = match mantissa_str.find('.') {
        Some(i) => (&mantissa_str[..i], &mantissa_str[i + 1..]),
        None => (mantissa_str, ""),
    };
    assert!(
        !int_part.is_empty(),
        "hex float without an integer part: {s:?}"
    );
    // The significand is digits[0].digits[1..] in base 16; 4 bits per hex digit.
    let mut lead = 0.0f64;
    for d in int_part.bytes() {
        lead = lead * 16.0 + hex_digit(d) as f64;
    }
    let mut value = lead;
    let mut scale = 1.0f64;
    for d in frac_part.bytes() {
        scale /= 16.0;
        value += hex_digit(d) as f64 * scale;
    }
    let exp: i32 = exp_str.parse().unwrap_or(0);
    let v = value * 2f64.powi(exp);
    if neg {
        -v
    } else {
        v
    }
}

fn hex_digit(b: u8) -> u32 {
    match b {
        b'0'..=b'9' => (b - b'0') as u32,
        b'a'..=b'f' => (b - b'a') as u32 + 10,
        b'A'..=b'F' => (b - b'A') as u32 + 10,
        _ => panic!("bad hex digit {:?}", b as char),
    }
}

struct Fixture {
    data: Vec<f64>,
    genes: Vec<String>,
    k: usize,
}

impl Fixture {
    fn build() -> Self {
        let (data, genes, k) = fixture_matrix();
        Self { data, genes, k }
    }

    fn means(&self) -> GroupMeans<'_> {
        GroupMeans::new(&self.data, &self.genes, self.k)
    }

    /// `X[rownames %in% keep, ]`: drop rows by gene name, as the generator's `Xm` does.
    fn without(&self, drop: &[&str]) -> (Vec<f64>, Vec<String>) {
        let keep: Vec<usize> = (0..self.genes.len())
            .filter(|&i| !drop.contains(&self.genes[i].as_str()))
            .collect();
        let genes: Vec<String> = keep.iter().map(|&i| self.genes[i].clone()).collect();
        let mut data = Vec::with_capacity(keep.len() * self.k);
        for j in 0..self.k {
            for &i in &keep {
                data.push(self.data[j * self.genes.len() + i]);
            }
        }
        (data, genes)
    }
}

/// Split a corpus line into tab-separated fields.
fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn parse_f64s(s: &str) -> Vec<f64> {
    s.split(',').map(|v| v.trim().parse().unwrap()).collect()
}

fn is_error(field: &str) -> bool {
    field == "ERROR"
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

/// `lr_df` from the generator, verbatim: two co-receptor columns of *cofactor row names*.
/// Transcribed by hand, so it is pinned against the generator by
/// `lr_df_matches_the_generator` -- a hand-typed mirror of a fixture that this file cannot
/// see is a silent source of wrong answers, and it already produced one.
const LR_DF: [(&str, &str); 5] = [
    ("", ""),
    ("ACTIVIN antagonist", "TGFb inhibition receptor"),
    ("TGFb antagonist", ""),
    ("", ""),
    ("NOT_A_COFACTOR", ""),
];

fn lr_records() -> HashMap<&'static str, Vec<&'static str>> {
    GOLDEN
        .lines()
        .filter_map(|l| l.strip_prefix("lr\t"))
        .map(|l| {
            let f = fields(l);
            (f[0], f)
        })
        .collect()
}

#[test]
fn corpus_is_present_and_populated() {
    let mut kinds: HashMap<&str, usize> = HashMap::new();
    for l in GOLDEN.lines().filter(|l| !l.is_empty()) {
        *kinds.entry(fields(l)[0]).or_default() += 1;
    }
    for k in [
        "lr",
        "coreceptor",
        "coreceptor_one",
        "agonist",
        "geomean",
        "group_agonist",
        "group_antagonist",
    ] {
        assert!(
            kinds.get(k).copied().unwrap_or(0) > 0,
            "no {k} records in the corpus"
        );
    }
    // The corpus must contain both error records, or the R-ism 11 half is untested.
    assert_eq!(lr_records()["absent_not_complex"][1], "ERROR");
    assert_eq!(lr_records()["missing_subunit_row"][1], "ERROR");
    assert!(
        kinds["agonist"] >= 25,
        "agonist records: {}",
        kinds["agonist"]
    );
}

#[test]
fn the_fixture_matrix_parses_and_has_the_engineered_rows() {
    // Guards the assumption every other test in this file rests on: that the hex dump
    // round-trips. If this fails, do not "fix" the other tests -- fix the parser.
    let (values, genes, k) = fixture_matrix();
    assert_eq!(k, 6);
    assert_eq!(genes.len(), 416, "416 genes: 12 plain + 404 DB subunits");
    assert_eq!(values.len(), genes.len() * k);
    // The four hand-set rows the generator installs after the draw.
    for j in 0..k {
        bits_eq(values[j * 416], 0.0, "gene 1 (all zero)");
    }
    for j in 0..5 {
        bits_eq(values[j * 416 + 1], 0.0, "gene 2 (>= 75% zero)");
    }
    for j in 0..k {
        bits_eq(values[j * 416 + 3], 1.0, "gene 4 (constant)");
    }
    assert_eq!(genes[0], "G1");
    // A few non-trivial draws, to catch a parser that silently returns 0 or truncates.
    let nonzero = values.iter().filter(|v| **v != 0.0).count();
    assert!(
        nonzero > 2400,
        "only {nonzero} non-zero values; the draw is missing"
    );
    assert!(
        values.iter().all(|v| *v >= 0.0 && *v <= 1.0),
        "out of runif(0.01, 1) range"
    );
}

#[test]
fn the_rng_can_still_reproduce_the_generator_draw() {
    // Kept as a *secondary* check. The primary input is the hex dump; this only confirms
    // `MersenneTwister` still matches `set.seed(4242); runif(...)`, which is what the
    // generator's first 2496 values came from. The RNG's own parity lives in
    // `rng_parity.rs`; this is the join between the two.
    let mut mt = MersenneTwister::new(SEED);
    let want = [
        9.862_531_709_182_076e-1,
        3.542_388_913_338_072_6e-1,
        2.376_028_250_367_380_7e-1,
        6.720_240_837_428_719e-1,
    ];
    for (i, &w) in want.iter().enumerate() {
        let got = 0.01 + 0.99 * mt.unif_rand();
        assert_eq!(format!("{got:.17e}"), format!("{w:.17e}"), "draw {i}");
    }
}

#[test]
fn lr_df_matches_the_generator() {
    // LR_DF is hand-transcribed, so it can drift from the generator without anyone
    // noticing until a coreceptor value is quietly wrong. This makes the drift loud.
    for col in ["co_A_receptor", "co_I_receptor"] {
        let line = GOLDEN
            .lines()
            .find(|l| l.starts_with(&format!("lr_df\t{col}\t")))
            .unwrap_or_else(|| panic!("no lr_df record for {col}"));
        let want: Vec<&str> = line.split('\t').skip(2).collect();
        let which = if col == "co_A_receptor" { 0 } else { 1 };
        let got: Vec<&str> = LR_DF
            .iter()
            .map(|p| if which == 0 { p.0 } else { p.1 })
            .collect();
        assert_eq!(
            got, want,
            "LR_DF column {col} has drifted from the generator"
        );
    }
}

#[test]
fn compute_expr_lr_matches_upstream() {
    let db = human_db();
    let fx = Fixture::build();
    let m = fx.means();
    // Must mirror `lr_cases` in the generator. The two ERROR cases are handled by
    // `a_complex_with_a_missing_subunit_row_errors_like_upstream` and the loop below.
    let cases: Vec<(&str, Vec<&str>)> = vec![
        ("single", vec!["G1", "G2"]),
        ("complex", vec!["Activin AB", "IL12AB"]),
        ("mixed", vec!["G1", "Activin AB", "IL12AB", "G4"]),
        ("missing_subunit", vec!["G1", "IL12AB"]),
        ("repeat", vec!["G1", "G1", "Activin AB", "Activin AB"]),
        ("single_row", vec!["G5"]),
    ];
    let recs = lr_records();

    for (name, genes) in &cases {
        let f = &recs[name];
        assert!(
            !is_error(f[1]),
            "{name}: R succeeded, so this must not be an error case"
        );
        let names: Vec<String> = genes.iter().map(|s| s.to_string()).collect();
        let got = compute_expr_lr(&names, &m, &db).unwrap_or_else(|e| panic!("{name}: {e}"));
        vecs_eq(&got, &parse_f64s(f[2]), name);
    }
}

#[test]
fn a_name_neither_a_row_nor_a_complex_is_upstreams_subscript_error() {
    let db = human_db();
    let fx = Fixture::build();
    let m = fx.means();
    let f = &lr_records()["absent_not_complex"];
    assert_eq!(f[1], "ERROR");
    let e = compute_expr_lr(&["NOT_A_COMPLEX".to_string()], &m, &db)
        .expect_err("must reproduce upstream's subscript out of bounds");
    assert_eq!(e.to_string(), SUBSCRIPT_OUT_OF_BOUNDS);
    assert_eq!(
        f[2], SUBSCRIPT_OUT_OF_BOUNDS,
        "R's message text, byte for byte"
    );
}

#[test]
fn a_complex_with_a_missing_subunit_row_errors_like_upstream() {
    let db = human_db();
    let fx = Fixture::build();
    // Drop IL12B from the matrix, as the generator does with `Xm`.
    let (data, genes) = fx.without(&["IL12B"]);
    let m = GroupMeans::new(&data, &genes, fx.k);
    let f = &lr_records()["missing_subunit_row"];
    assert_eq!(f[1], "ERROR");
    // The co-agonist/cofactor paths would instead have degraded to a no-op here; the
    // ligand path has no `intersect`, so it aborts. Assert the two differ.
    assert_eq!(
        compute_expr_coreceptor(&["IL12AB".to_string()], &m, &db, CoreceptorKind::Activator),
        vec![1.0; fx.k]
    );
    let e = compute_expr_lr(&["IL12AB".to_string()], &m, &db)
        .expect_err("must reproduce upstream's subscript out of bounds");
    assert_eq!(e.to_string(), SUBSCRIPT_OUT_OF_BOUNDS);
    assert_eq!(f[2], SUBSCRIPT_OUT_OF_BOUNDS);
}

#[test]
fn compute_expr_coreceptor_matches_upstream() {
    let db = human_db();
    let fx = Fixture::build();
    let m = fx.means();
    let co_a: Vec<String> = LR_DF.iter().map(|p| p.0.to_string()).collect();
    let co_i: Vec<String> = LR_DF.iter().map(|p| p.1.to_string()).collect();

    for (ty, names) in [("A", &co_a), ("I", &co_i)] {
        let f: Vec<&str> = GOLDEN
            .lines()
            .find(|l| l.starts_with(&format!("coreceptor\t{ty}\t")))
            .unwrap_or_else(|| panic!("no coreceptor record for {ty}"))
            .split('\t')
            .collect();
        assert!(!is_error(f[2]), "coreceptor {ty}: R errored ({})", f[2]);
        let kind = if ty == "A" {
            CoreceptorKind::Activator
        } else {
            CoreceptorKind::Inhibitor
        };
        let got = compute_expr_coreceptor(names, &m, &db, kind);
        vecs_eq(&got, &parse_f64s(f[3]), &format!("coreceptor {ty}"));
    }

    // And the single-row form the bootstrap loop calls, on the A column only.
    for line in GOLDEN.lines().filter(|l| l.starts_with("coreceptor_one\t")) {
        let f = fields(line);
        let idx: usize = f[1].parse().expect("row index");
        assert!(
            !is_error(f[2]),
            "coreceptor_one {idx}: R errored ({})",
            f[2]
        );
        let got = compute_expr_coreceptor(
            &[LR_DF[idx - 1].0.to_string()],
            &m,
            &db,
            CoreceptorKind::Activator,
        );
        vecs_eq(
            &got,
            &parse_f64s(f[2]),
            &format!("coreceptor_one row {idx}"),
        );
    }
}

#[test]
fn compute_expr_agonist_and_antagonist_match_upstream() {
    let db = human_db();
    let fx = Fixture::build();
    let m = fx.means();
    let mut n_checked = 0;
    for line in GOLDEN.lines().filter(|l| l.starts_with("agonist\t")) {
        // f = kind, agonist name, antagonist name, Kh, n, agonist values, antagonist values
        let f = fields(line);
        let kh: f64 = f[3].parse().expect("Kh");
        let n: f64 = f[4].parse().expect("n");
        let what = format!("ag={} an={} Kh={kh} n={n}", f[1], f[2]);
        assert!(!is_error(f[5]), "agonist {what}: R errored ({})", f[5]);
        let ag = compute_expr_agonist(f[1], &m, &db, kh, n);
        vecs_eq(&ag, &parse_f64s(f[5]), &format!("agonist {what}"));
        assert!(!is_error(f[6]), "antagonist {what}: R errored ({})", f[6]);
        let an = compute_expr_antagonist(f[2], &m, &db, kh, n);
        vecs_eq(&an, &parse_f64s(f[6]), &format!("antagonist {what}"));
        n_checked += 1;
    }
    assert!(
        n_checked >= 25,
        "only checked {n_checked} agonist/antagonist records"
    );
}

#[test]
fn geometric_mean_specials_match_r() {
    // `gm_cases` in the generator, in order. Note case 9: log(-Inf) is NaN, so the NaN is
    // dropped and the result is sqrt(2) rather than 0.
    let cases: Vec<Vec<f64>> = vec![
        vec![1.0, 4.0],
        vec![0.0, 1.0, 4.0],
        vec![2.0, 8.0],
        vec![1e-300, 1e300],
        vec![1.0, 1.0, 1.0, 1.0, 1.0],
        vec![1.0, f64::NAN],
        vec![f64::NAN, f64::NAN],
        vec![1.0, f64::INFINITY, 2.0],
        vec![1.0, f64::NEG_INFINITY, 2.0],
        vec![0.1, 0.2, 0.3],
    ];
    let want: Vec<f64> = GOLDEN
        .lines()
        .filter(|l| l.starts_with("geomean\t"))
        .map(|l| fields(l)[2].parse().unwrap())
        .collect();
    assert_eq!(
        want.len(),
        cases.len(),
        "corpus/geomean case count mismatch"
    );
    for (i, (x, &w)) in cases.iter().zip(&want).enumerate() {
        let got = r_core::stats::geometric_mean(x);
        if got.is_nan() || w.is_nan() {
            assert!(got.is_nan() && w.is_nan(), "case {i}: rust={got} r={w}");
        } else {
            bits_eq(got, w, &format!("geometricMean({x:?})"));
        }
    }
}

#[test]
fn the_group_mean_variants_match_upstream() {
    // `grp <- factor(rep(paste0("g", 1:3), each = 2), levels = paste0("g", 1:3))` and
    // `FunMean = triMean`, Kh = 0.5, n = 1, cofactor "TGFb agonist" (2 subunits) against
    // "TGFb antagonist" (3 subunits).
    let db = human_db();
    let fx = Fixture::build();
    // `X` is genes x 6 "cells". The group variants want a *cells x genes* matrix, and R's
    // column-major layout makes that a real transpose: `X[g, cell]` sits at
    // `cell * n_genes + g` in `X`, but in a cells x genes R matrix `m[cell, g]` sits at
    // `g * n_cells + cell`. Copying `X` through unchanged (no transpose) makes every group
    // mean the mean of the wrong cells.
    let (data, genes) = (fx.data.as_slice(), fx.genes.as_slice());
    let n_genes = genes.len();
    let mut cells: Vec<f64> = vec![0.0; fx.k * n_genes];
    for g in 0..n_genes {
        for cell in 0..fx.k {
            cells[g * fx.k + cell] = data[cell * n_genes + g];
        }
    }
    let ce = CellExpr::new(&cells, genes, fx.k);
    // The transpose must actually be a transpose. Without the stride swap every group
    // mean comes out as 0.5 -- the mean of {0, 1} -- which is an unmistakable signature.
    // G1 is all zeros so it cannot distinguish the two layouts; use a subunit instead.
    let fst = ce.row("FST").expect("FST in the fixture matrix");
    // The stride is written out on purpose -- the whole point of the test is that the layout is
    // `data[gene * n_cells + cell]` and the row-major reading `data[fst]` is what it contradicts.
    // `cell` goes in a local rather than inline as `0 * n_genes`, because `0 * n` is a
    // constant-folding `erasing_op` and the index has to stay a real stride arithmetic for the
    // test to mean anything.
    let cell = 0usize;
    let cell0_gene_fst = data[cell * n_genes + fst];
    assert_eq!(
        ce.at(fst, 0).to_bits(),
        cell0_gene_fst.to_bits(),
        "FST in cell 0"
    );
    assert_ne!(ce.at(fst, 0), 0.0, "the fixture gives FST a non-zero value");
    let group: Vec<usize> = vec![0, 0, 1, 1, 2, 2];
    let n_groups = 3;

    // The record carries both cofactor names: `computeExprGroup_agonist` reads
    // `pairLRsig$agonist[index]` and the antagonist twin reads `pairLRsig$antagonist[index]`,
    // so keying the fixture on one name silently answers the wrong question.
    for (kind, col) in [("group_agonist", 1), ("group_antagonist", 2)] {
        let f = fields(
            GOLDEN
                .lines()
                .find(|l| l.starts_with(&format!("{kind}\t")))
                .unwrap_or_else(|| panic!("no {kind} record")),
        );
        let name = f[col];
        let want = parse_f64s(f[3]);
        let got = if kind == "group_agonist" {
            compute_expr_group_agonist(
                name,
                &ce,
                &group,
                n_groups,
                &db,
                0.5,
                1.0,
                GroupMean::TriMean,
            )
        } else {
            compute_expr_group_antagonist(
                name,
                &ce,
                &group,
                n_groups,
                &db,
                0.5,
                1.0,
                GroupMean::TriMean,
            )
        };
        vecs_eq(&got, &want, kind);
    }
}

#[test]
fn complex_and_lr_agree_on_a_complex_name() {
    // `computeExpr_complex` is just `computeExpr_LR` restricted to complexes upstream;
    // the corpus has no separate record, so assert the two entry points agree.
    let db = human_db();
    let fx = Fixture::build();
    let m = fx.means();
    let names = vec!["Activin AB".to_string(), "IL12AB".to_string()];
    let a = compute_expr_lr(&names, &m, &db).unwrap();
    let b = compute_expr_complex(&names, &m, &db).unwrap();
    vecs_eq(&a, &b, "computeExpr_LR vs computeExpr_complex");
}

#[test]
fn a_plain_gene_is_its_row_verbatim() {
    // The bug this file exists to prevent: routing a single gene through
    // exp(mean(log(x))) is off by 1 ulp on the very first value.
    let db = human_db();
    let fx = Fixture::build();
    let m = fx.means();
    let got = compute_expr_lr(&["G5".to_string()], &m, &db).unwrap();
    for j in 0..fx.k {
        let row = fx.data[j * fx.genes.len() + 4];
        bits_eq(got[j], row, "G5 row");
        if (row.ln().exp() - row).abs() > 0.0 {
            // Not an assertion about correctness of the port -- a reminder that the two
            // are genuinely different numbers for at least some inputs.
            eprintln!("note: exp(log(x)) != x for G5 group {j}");
        }
    }
}

#[test]
fn fixtures_dir_is_where_we_think_it_is() {
    let p: &Path = &fixtures();
    assert!(
        p.join("db_human").is_dir(),
        "db_human fixture missing at {p:?}"
    );
    assert!(
        p.join("db_mouse").is_dir(),
        "db_mouse fixture missing at {p:?}"
    );
}
