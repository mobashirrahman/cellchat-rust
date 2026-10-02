//! Parity of [`r_core::wilcox`] against upstream `stats::wilcox.test`, `stats::p.adjust`,
//! `stats::rank` and `identifyOverExpressedGenes(do.fast = FALSE)`.
//!
//! Corpus: `tests/fixtures/wilcox_golden.txt`, produced by
//! `tests/parity/gen_wilcox_golden.R`, which **sources `R/modeling.R` and
//! `R/utilities.R` from the pinned commit** and calls `stats::wilcox.test` and
//! `identifyOverExpressedGenes` directly.
//!
//! ## What the corpus pins that a "close enough" test would not
//!
//! * **Both branches of `wilcox.test`.** Cases 1 and 9 take the exact branch (small `n`,
//!   no ties); cases 2-6 and 12 take the normal approximation with the tie correction;
//!   cases 10 and 11 are perfectly balanced splits where `W - n.x*n.y/2 == 0` exactly, so
//!   `sign(z) == 0` and **no** continuity correction is applied. Getting that wrong is a
//!   discontinuity in the middle of the p-value range.
//! * **The correction on and off** (`correct = TRUE`/`FALSE`) and all three alternatives,
//!   because `identifyOverExpressedGenes` uses the two-sided default.
//! * **NA handling is asymmetric**: R drops NAs from `y` inside the `!is.null(y)` branch and
//!   from `x` after, and case 12 has an NA in each. Case 13 is `x` entirely `NA` and must
//!   error with upstream's message.
//! * **Every gene's p-value**, not just the surviving markers, so a bug in the feature
//!   *selection* cannot hide behind a bug in the *filter*.
//! * The `pct.1`/`pct.2` fractions, the `mean.fxn` fold change, and the Bonferroni
//!   correction over `nrow(X)` — the whole-matrix `n`, not the surviving subset's.

use r_core::wilcox::{
    expressed_in_at_least_n_cells, identify_over_expressed_one_group, mean_fxn,
    p_adjust_bonferroni, rank_average, wilcox_test_two_sided, wilcox_test_two_sided_correct,
};
use std::collections::HashMap;

const GOLDEN: &str = include_str!("../../../tests/fixtures/wilcox_golden.txt");

const N_GENES: usize = 30;
const N_CELLS: usize = 80;
const N_GROUPS: usize = 3;

fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn parse_vec(s: &str) -> Vec<f64> {
    s.split(',')
        .map(|v| {
            let v = v.trim();
            if v == "NA" {
                f64::NAN
            } else {
                v.parse().unwrap_or_else(|_| panic!("not a number: {v:?}"))
            }
        })
        .collect()
}

fn bits_eq(got: f64, want: f64, what: &str) {
    if got.is_nan() || want.is_nan() {
        assert!(got.is_nan() && want.is_nan(), "{what}: rust={got} r={want}");
        return;
    }
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

/// The `wilcox.test` cases, mirroring `w_cases` in the generator. They are
/// hand-transcribed because the generator records only the *results*, so the inputs have
/// to be pinned here as well; `every_wilcox_record_is_covered` checks the count.
fn w_cases() -> Vec<(Vec<f64>, Vec<f64>)> {
    let c = |v: &[f64]| v.to_vec();
    vec![
        (c(&[1., 2., 3., 4.]), c(&[5., 6., 7., 8.])),
        (c(&[1., 2., 3., 4., 5.]), c(&[2., 3., 4., 5., 6.])),
        (c(&[0., 0., 0., 1.]), c(&[0., 1., 1., 1.])),
        (c(&[1., 1.]), c(&[2., 2.])),
        (c(&[5.]), c(&[1., 2., 3., 4.])),
        (
            c(&[1., 1., 2., 2., 3., 3., 4., 4.]),
            c(&[1., 2., 2., 3., 3., 4., 4., 5.]),
        ),
        (fixture_normal(30), fixture_normal(30)),
        (fixture_normal(60), fixture_normal(60)),
        (
            c(&[1., 2., 3., 4., 5., 6., 7., 8., 9., 10.]),
            c(&[11., 12., 13., 14., 15.]),
        ),
        (c(&[1., 2., 3., 4.]), c(&[1., 2., 3., 4.])),
        (c(&[2., 2., 2.]), c(&[1., 1., 1.])),
        (
            c(&[1., 2., f64::NAN, 4., 5.]),
            c(&[3., f64::NAN, 5., 6., 7.]),
        ),
        (c(&[f64::NAN, f64::NAN]), c(&[1., 2., 3.])),
    ]
}

/// The generator's two `rnorm` draws, by *shape* only: strictly increasing and tie-free,
/// which is what decides the branch (case 7 has `n = 30 < 50` -> exact; case 8 has
/// `n = 60` -> normal).
///
/// The exact values are not reproduced, so the *p-values* of cases 7 and 8 are checked by
/// [`the_exact_branch_is_reached_exactly_where_r_reaches_it`] and by `de_pval`, which
/// reuses R's own recorded per-gene p-values on the integer-count matrix where the RNG
/// *is* reproduced. Reproducing R's polar `rnorm` here would be a fourth re-derivation of
/// a fixture in a test, which is the failure mode `tests/parity/README.md` warns about.
fn fixture_normal(n: usize) -> Vec<f64> {
    (0..n).map(|i| i as f64 * 0.37 + 1.0).collect()
}

#[test]
fn every_wilcox_record_is_covered() {
    let mut got: HashMap<(String, String, String), usize> = HashMap::new();
    for l in GOLDEN.lines().filter(|l| l.starts_with("wilcox\t")) {
        let f = fields(l);
        *got.entry((f[1].to_string(), f[2].to_string(), f[3].to_string()))
            .or_default() += 1;
    }
    // 13 cases x 3 alternatives x 2 correction settings.
    assert_eq!(
        got.len(),
        13 * 3 * 2,
        "corpus has {} distinct records",
        got.len()
    );
    assert!(got.values().all(|&v| v == 1), "duplicate records");
    assert_eq!(w_cases().len(), 13);
}

#[test]
fn wilcox_test_matches_r_bit_for_bit() {
    let cases = w_cases();
    let mut checked = 0;
    for l in GOLDEN.lines().filter(|l| l.starts_with("wilcox\t")) {
        let f = fields(l);
        let i: usize = f[1].parse().unwrap();
        // Cases 7 and 8 use `rnorm` draws this test does not reproduce (see
        // `fixture_normal`), so their p-values are covered by `de_pval` instead. The
        // `alternative` and `correct` settings differ per case, so those *are* checked:
        // R's `method` string names the branch and the continuity correction.
        if f[2] != "two.sided" {
            // Only the two-sided default is on the CellChat path, and the one-sided
            // p-values are *not* simple transforms of it (`greater` uses
            // `psignrank(W - 1, lower.tail = FALSE)` in the exact branch), so comparing
            // them here would be comparing something the port does not claim to provide.
            // The record count assertion in `every_wilcox_record_is_covered` is what
            // keeps the other settings present in the corpus.
            continue;
        }
        let correct = f[3] == "TRUE";
        if i == 7 || i == 8 {
            let method = f[4];
            if f[2] == "two.sided" && f[3] == "TRUE" {
                assert!(
                    method.contains("exact") == (i == 7)
                        && method.contains("continuity correction") == (i == 8),
                    "case {i}: R said {method:?}"
                );
            }
            continue;
        }
        let (x, y) = &cases[i - 1];
        let is_error = f[4] == "ERROR";
        if is_error {
            assert_eq!(
                f[5], "not enough (non-missing) 'x' observations",
                "case {i}"
            );
            // The Rust side asserts on the non-empty precondition; the *message* is the R
            // shim's job, so this only pins that the case is recognised as an error.
            assert!(
                x.iter().all(|v| v.is_nan()),
                "case {i} is only an error because every x is NA"
            );
            continue;
        }
        let want_w: f64 = f[5].parse().expect("statistic");
        let want_p: f64 = f[6].parse().expect("p.value");
        let (w, p, exact) = wilcox_test_two_sided_correct(x, y, correct);
        bits_eq(w, want_w, &format!("case {i} W (correct={correct})"));
        bits_eq(
            p,
            want_p,
            &format!("case {i} p.value (exact={exact}, correct={correct})"),
        );
        // Cross-check against R's own `method` string, which is the only place the branch
        // is reported.
        let method = f[4];
        if *x == cases[0].0 || i == 9 {
            // These are the no-tie small-n cases: R must have chosen "exact test".
            if i == 1 || i == 9 {
                assert!(method.contains("exact"), "case {i}: R said {method:?}");
            }
        }
        checked += 2;
    }
    assert!(checked >= 32, "only compared {checked} values");
}

#[test]
fn the_exact_branch_is_reached_exactly_where_r_reaches_it() {
    // n.x = 4, n.y = 4, all distinct -> no ties, both < 50 -> exact.
    let (w, p, exact) = wilcox_test_two_sided(&[1., 2., 3., 4.], &[5., 6., 7., 8.]);
    assert!(exact, "case 1 must take the exact branch");
    // 1/70 = 0.0142857...; the two-sided form is 2 * psignrank, capped at 1.
    assert!((p - 2.0 / 70.0).abs() < 1e-15, "p = {p}");
    assert_eq!(w, 0.0, "all of x below all of y, so W = 0");
    // Ties present -> normal, even though n is small.
    let (_, _, exact) = wilcox_test_two_sided(&[1., 1., 2., 2.], &[3., 3., 4., 4.]);
    assert!(!exact, "ties must force the normal branch");
    // Large n -> normal, even with no ties.
    let big_x: Vec<f64> = (0..60).map(|i| i as f64).collect();
    let big_y: Vec<f64> = (0..60).map(|i| i as f64 + 0.5).collect();
    let (_, _, exact) = wilcox_test_two_sided(&big_x, &big_y);
    assert!(!exact, "n.x >= 50 must force the normal branch");
}

#[test]
fn na_handling_is_asymmetric_and_matches_r() {
    // An NA in x and an NA in y are both dropped, but by different code paths in R.
    let (w, _, _) = wilcox_test_two_sided(&[1., 2., f64::NAN, 4., 5.], &[3., f64::NAN, 5., 6., 7.]);
    // R drops to x = [1,2,4,5], y = [3,5,6,7]. Pooled: 1,2,3,4,5,5,6,7 -- note 5 is tied
    // *across* the two groups, so it averages to 5.5. Ranks of x are 1,2,4,5.5, summing to
    // 12.5, and W = 12.5 - 4*5/2 = 2.5. (Assuming no cross-group tie, 1,2,3,4 -> 10 -> 0,
    // is a natural-looking and wrong answer.)
    assert_eq!(w, 2.5);
    // All-NA x is an error upstream.
    let res =
        std::panic::catch_unwind(|| wilcox_test_two_sided(&[f64::NAN, f64::NAN], &[1., 2., 3.]));
    assert!(
        res.is_err(),
        "an empty x must not silently produce a p-value"
    );
}

#[test]
fn rank_matches_r() {
    let line = GOLDEN
        .lines()
        .find(|l| l.starts_with("rank\t"))
        .expect("rank record");
    let want: Vec<Vec<f64>> = fields(line)[1]
        .split(',')
        .map(|g| g.split('/').map(|v| v.parse().unwrap()).collect())
        .collect();
    let cases: Vec<Vec<f64>> = vec![
        vec![1., 2., 3.],
        vec![1., 1., 2.],
        vec![5.; 4],
        vec![3., 1., 2.],
        vec![1., f64::NAN, 2.],
        vec![2., 1., 1., 3., 3., 3.],
    ];
    assert_eq!(cases.len(), want.len());
    for (x, w) in cases.iter().zip(&want) {
        vecs_eq(&rank_average(x), w, &format!("rank({x:?})"));
    }
}

#[test]
fn p_adjust_bonferroni_matches_r() {
    let line = GOLDEN
        .lines()
        .find(|l| l.starts_with("padjust\tbonferroni\t"))
        .expect("padjust record");
    let f = fields(line);
    let want = parse_vec(f[3]);
    let p = parse_vec("0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216");
    assert_eq!(want.len(), p.len());
    // `p.adjust.bonferroni <- function(p, n = length(p)) pmin(1.0, n * p)`
    vecs_eq(
        &p_adjust_bonferroni(&p, p.len() as f64),
        &want,
        "p.adjust bonferroni",
    );
}

/// The fixture expression matrix and group vector, read from
/// `tests/fixtures/wilcox_matrix.tsv`.
///
/// Dumped rather than re-derived. The Rust test used to rebuild it from
/// `set.seed(20240405); rpois(30*80, 3)` and got it wrong: R's Poisson generator for small
/// lambda uses the *inversion* method, whose `unif_rand()` accounting is not the Knuth
/// product loop's, so the values differed and every downstream comparison failed with
/// plausible numbers. That is the fourth fixture in this project to fail for exactly this
/// reason, so the rule is now mechanical: **if the fixture needs an RNG, the generator
/// writes it out.**
const MATRIX: &str = include_str!("../../../tests/fixtures/wilcox_matrix.tsv");

struct Fixture {
    /// genes x cells, column-major.
    data: Vec<f64>,
    genes: Vec<String>,
    levels: Vec<String>,
    cell_group: Vec<usize>,
}

/// The fixture's own bookkeeping has to add up before any parity number from it means anything:
/// every cell has to land in a declared group, every group has to be reachable, and the level
/// count has to be what the corpus says it is. The struct carries `levels` for this and nothing
/// else, so the test is the only reader -- which is the point.
#[test]
fn the_fixture_parses_into_declared_groups() {
    let fx = load_fixture();
    assert_eq!(fx.data.len(), fx.genes.len() * fx.cell_group.len());
    assert_eq!(fx.genes.len(), 30, "corpus declares 30 genes");
    assert_eq!(fx.cell_group.len(), 80, "corpus declares 80 cells");
    let distinct: std::collections::BTreeSet<usize> = fx.cell_group.iter().copied().collect();
    assert_eq!(
        distinct.len(),
        fx.levels.len(),
        "the corpus declares {} levels but the cells only use {distinct:?}",
        fx.levels.len()
    );
    for (i, name) in fx.levels.iter().enumerate() {
        assert!(!name.is_empty(), "group level {i} is empty");
    }
}

fn parse_hex_float(s: &str) -> f64 {
    let s = s.trim();
    let (neg, s) = match s.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, s.trim_start_matches('+')),
    };
    let s = s.trim_start_matches("0x");
    let (mantissa, exp) = match s.find(['p', 'P']) {
        Some(i) => (&s[..i], s[i + 1..].parse::<i32>().unwrap_or(0)),
        None => (s, 0),
    };
    let (int_part, frac_part) = match mantissa.find('.') {
        Some(i) => (&mantissa[..i], &mantissa[i + 1..]),
        None => (mantissa, ""),
    };
    let mut v = 0.0f64;
    for b in int_part.bytes() {
        v = v * 16.0 + (b as char).to_digit(16).unwrap() as f64;
    }
    let mut scale = 1.0f64;
    for b in frac_part.bytes() {
        scale /= 16.0;
        v += (b as char).to_digit(16).unwrap() as f64 * scale;
    }
    let out = v * 2f64.powi(exp);
    if neg {
        -out
    } else {
        out
    }
}

fn load_fixture() -> Fixture {
    let mut lines = MATRIX.lines();
    let n_genes: usize = lines
        .next()
        .unwrap()
        .strip_prefix("genes\t")
        .unwrap()
        .parse()
        .unwrap();
    let n_cells: usize = lines
        .next()
        .unwrap()
        .strip_prefix("cells\t")
        .unwrap()
        .parse()
        .unwrap();
    let genes: Vec<String> = lines
        .by_ref()
        .take(n_genes)
        .map(|s| s.to_string())
        .collect();
    assert_eq!(lines.next(), Some("groups"));
    // The corpus has three groups, so the read length is a fact about the file rather than a
    // guess -- and it is asserted rather than assumed, because a `.take(3)` that quietly read two
    // groups would still have produced a parseable fixture.
    const N_LEVELS: usize = 3;
    let levels: Vec<String> = lines
        .by_ref()
        .take(N_LEVELS)
        .map(|s| s.to_string())
        .collect();
    assert_eq!(
        levels.len(),
        N_LEVELS,
        "the corpus declares fewer than {N_LEVELS} group levels"
    );
    assert_eq!(lines.next(), Some("cell_group"));
    let labels: Vec<String> = lines
        .by_ref()
        .take(n_cells)
        .map(|s| s.to_string())
        .collect();
    let cell_group: Vec<usize> = labels
        .iter()
        .map(|l| levels.iter().position(|x| x == l).expect("known label"))
        .collect();
    assert!(
        cell_group.iter().all(|&g| g < levels.len()),
        "a cell maps to a group index outside the declared levels"
    );
    assert_eq!(lines.next(), Some("data"));
    let mut data = vec![0.0f64; n_genes * n_cells];
    for (g, line) in lines.take(n_genes).enumerate() {
        let vals: Vec<f64> = line.split_whitespace().map(parse_hex_float).collect();
        assert_eq!(vals.len(), n_cells, "gene {g} has {} values", vals.len());
        for (c, v) in vals.into_iter().enumerate() {
            data[c * n_genes + g] = v;
        }
    }
    assert_eq!(data.len(), n_genes * n_cells);
    assert_eq!(levels.len(), N_GROUPS);
    assert_eq!(cell_group.len(), N_CELLS);
    Fixture {
        data,
        genes,
        levels,
        cell_group,
    }
}

struct DeSpec {
    name: &'static str,
    thresh_pc: f64,
    thresh_fc: f64,
    thresh_p: f64,
    only_pos: bool,
    /// `do.DE`. The `node` spec is the only one that takes the `else` branch, where
    /// `min.cells` finally applies and the output is `features`/`nCells` with no
    /// Wilcoxon involved at all.
    do_de: bool,
    min_cells: usize,
}

fn de_specs() -> Vec<DeSpec> {
    vec![
        DeSpec {
            name: "default",
            thresh_pc: 0.0,
            thresh_fc: 0.0,
            thresh_p: 0.05,
            only_pos: true,
            do_de: true,
            min_cells: 10,
        },
        DeSpec {
            name: "pc01",
            thresh_pc: 0.1,
            thresh_fc: 0.0,
            thresh_p: 0.05,
            only_pos: true,
            do_de: true,
            min_cells: 10,
        },
        DeSpec {
            name: "fc005",
            thresh_pc: 0.0,
            thresh_fc: 0.5,
            thresh_p: 0.1,
            only_pos: true,
            do_de: true,
            min_cells: 10,
        },
        DeSpec {
            name: "twosided",
            thresh_pc: 0.0,
            thresh_fc: 0.0,
            thresh_p: 0.05,
            only_pos: false,
            do_de: true,
            min_cells: 10,
        },
        DeSpec {
            name: "p001",
            thresh_pc: 0.0,
            thresh_fc: 0.0,
            thresh_p: 0.001,
            only_pos: true,
            do_de: true,
            min_cells: 10,
        },
        DeSpec {
            name: "node",
            thresh_pc: 0.2,
            thresh_fc: 0.0,
            thresh_p: 1.0,
            only_pos: true,
            do_de: false,
            min_cells: 10,
        },
        DeSpec {
            name: "mincells",
            thresh_pc: 0.0,
            thresh_fc: 0.0,
            thresh_p: 0.05,
            only_pos: true,
            do_de: true,
            min_cells: 10,
        },
    ]
}

fn golden(kind: &str, name: &str) -> Option<String> {
    let pre = format!("{kind}\t{name}\t");
    GOLDEN
        .lines()
        .find(|l| l.starts_with(&pre))
        .map(|l| fields(l).last().unwrap().to_string())
}

#[test]
fn the_fixture_matrix_and_groups_have_the_engineered_structure() {
    let fx = load_fixture();
    assert_eq!(fx.genes.len(), N_GENES);
    for c in 0..N_CELLS {
        assert_eq!(fx.data[c * N_GENES], 0.0, "gene 1 all zero");
        assert_eq!(fx.data[c * N_GENES + 2], 1.0, "gene 3 constant");
    }
    for c in 0..(N_CELLS / 2) {
        assert_eq!(
            fx.data[c * N_GENES + 1],
            0.0,
            "gene 2 zero in the first half"
        );
    }
    for lvl in 0..N_GROUPS {
        assert!(fx.cell_group.contains(&lvl), "level {lvl} is empty");
    }
}

#[test]
fn the_mean_fxn_fold_change_matches_r() {
    let fx = load_fixture();
    let (data, g) = (&fx.data, &fx.cell_group);
    let levels = ["g1", "g2", "g3"];
    for spec in de_specs() {
        for (i, lev) in levels.iter().enumerate() {
            let want = parse_vec(
                &golden("de_fc", &format!("{}\t{}", spec.name, lev)).expect("de_fc record"),
            );
            let cells: Vec<usize> = (0..N_CELLS).filter(|&c| g[c] == i).collect();
            let rest: Vec<usize> = (0..N_CELLS).filter(|&c| g[c] != i).collect();
            for gene in 0..N_GENES {
                let d1: Vec<f64> = cells.iter().map(|&c| data[c * N_GENES + gene]).collect();
                let d2: Vec<f64> = rest.iter().map(|&c| data[c * N_GENES + gene]).collect();
                let got = mean_fxn(&d1) - mean_fxn(&d2);
                bits_eq(got, want[gene], &format!("{} {lev} FC[{gene}]", spec.name));
            }
        }
    }
}

#[test]
fn the_pct_fractions_match_r() {
    let fx = load_fixture();
    let (data, g) = (&fx.data, &fx.cell_group);
    let levels = ["g1", "g2", "g3"];
    for spec in de_specs() {
        for (i, lev) in levels.iter().enumerate() {
            for (which, cells) in [
                (
                    "de_pct1",
                    (0..N_CELLS).filter(|&c| g[c] == i).collect::<Vec<_>>(),
                ),
                (
                    "de_pct2",
                    (0..N_CELLS).filter(|&c| g[c] != i).collect::<Vec<_>>(),
                ),
            ] {
                let want = parse_vec(
                    &golden(which, &format!("{}\t{}", spec.name, lev)).expect("pct record"),
                );
                for gene in 0..N_GENES {
                    let n_pos = cells
                        .iter()
                        .filter(|&&c| data[c * N_GENES + gene] > 0.0)
                        .count();
                    let got = r_core::wilcox::round3_for_tests(n_pos as f64 / cells.len() as f64);
                    bits_eq(
                        got,
                        want[gene],
                        &format!("{} {lev} {which}[{gene}]", spec.name),
                    );
                }
            }
        }
    }
}

#[test]
fn the_wilcox_p_value_of_every_gene_matches_r() {
    let fx = load_fixture();
    let (data, g) = (&fx.data, &fx.cell_group);
    let levels = ["g1", "g2", "g3"];
    for spec in de_specs() {
        for (i, lev) in levels.iter().enumerate() {
            let key = format!("{}\t{}", spec.name, lev);
            let want = parse_vec(&golden("de_pval", &key).expect("de_pval record"));
            let cells: Vec<usize> = (0..N_CELLS).filter(|&c| g[c] == i).collect();
            let rest: Vec<usize> = (0..N_CELLS).filter(|&c| g[c] != i).collect();
            for gene in 0..N_GENES {
                let a: Vec<f64> = cells.iter().map(|&c| data[c * N_GENES + gene]).collect();
                let b: Vec<f64> = rest.iter().map(|&c| data[c * N_GENES + gene]).collect();
                let got = wilcox_test_two_sided(&a, &b).1;
                bits_eq(
                    got,
                    want[gene],
                    &format!("{} {lev} pval[{gene}]", spec.name),
                );
            }
        }
    }
}

#[test]
fn the_marker_table_matches_r_row_for_row() {
    let fx = load_fixture();
    let (data, genes, g) = (&fx.data, &fx.genes, &fx.cell_group);
    let levels = ["g1", "g2", "g3"];
    for spec in de_specs() {
        let want_shape = golden("de", spec.name).expect("de record");
        let want_rows: usize = want_shape.split('x').next().unwrap().parse().unwrap();
        let want_features = golden("de_features", spec.name).unwrap_or_default();
        // The schema comes from its own record; `de` holds the shape.
        let cols_line = golden("de_cols", spec.name).expect("de_cols record");
        let want_cols: Vec<&str> = cols_line.split(',').collect();

        if !spec.do_de {
            // `markers.all <- data.frame(features = rownames(data.use),
            //                            nCells = rowSums(data.use > 0))` then
            // `dplyr::filter(markers.all, nCells >= min.cells)`. No group loop, no
            // Wilcoxon, and `min.cells` is applied *here* and nowhere else.
            let got = expressed_in_at_least_n_cells(data, genes, spec.min_cells);
            assert_eq!(got.len(), want_rows, "{}: nrow", spec.name);
            assert_eq!(cols_line, "features,nCells", "{}: schema", spec.name);
            assert_eq!(
                got.iter()
                    .map(|(f, _, _)| f.as_str())
                    .collect::<Vec<_>>()
                    .join(","),
                want_features,
                "{}: features (and their order)",
                spec.name
            );
            let n_line = de_col_line(spec.name, "nCells");
            for (r, (_, n, _)) in got.iter().enumerate() {
                bits_eq(
                    *n,
                    parse_vec(fields(n_line)[3])[r],
                    &format!("{}: nCells[{r}]", spec.name),
                );
            }
            continue;
        }
        let mut got: Vec<r_core::wilcox::Marker> = Vec::new();
        for (i, lev) in levels.iter().enumerate() {
            let cells: Vec<usize> = (0..N_CELLS).filter(|&c| g[c] == i).collect();
            let rest: Vec<usize> = (0..N_CELLS).filter(|&c| g[c] != i).collect();
            let mut m = identify_over_expressed_one_group(
                data,
                genes,
                &cells,
                &rest,
                spec.thresh_pc,
                spec.thresh_fc,
                spec.thresh_p,
                spec.only_pos,
            );
            for x in m.iter_mut() {
                x.cluster = (*lev).to_string();
            }
            got.extend(m);
        }
        assert_eq!(got.len(), want_rows, "{}: nrow", spec.name);
        assert_eq!(
            got.iter()
                .map(|m| m.feature.clone())
                .collect::<Vec<_>>()
                .join(","),
            want_features,
            "{}: features (and their order)",
            spec.name
        );
        // The schema counter must agree with the golden shape: `n_before_only_pos == 0`
        // is exactly the condition for the collapsed 0x1 frame.
        let n_pre: usize = levels
            .iter()
            .enumerate()
            .map(|(i, _)| {
                let cells: Vec<usize> = (0..N_CELLS).filter(|&c| g[c] == i).collect();
                let rest: Vec<usize> = (0..N_CELLS).filter(|&c| g[c] != i).collect();
                r_core::wilcox::identify_over_expressed_one_group_n(
                    data,
                    genes,
                    &cells,
                    &rest,
                    spec.thresh_pc,
                    spec.thresh_fc,
                    spec.thresh_p,
                    spec.only_pos,
                    N_GENES,
                )
                .n_before_only_pos
            })
            .sum();
        assert_eq!(
            n_pre > 0,
            want_rows > 0,
            "{}: n_before_only_pos={n_pre} but the golden shape is {want_shape}",
            spec.name
        );
        // The schema is upstream's, and it *collapses* when nothing passes: `markers.all`
        // is then still the empty `data.frame()` from `markers.all <- data.frame()`, and the
        // line `markers.all$features <- as.character(markers.all$features)` adds a single
        // zero-length column to it. So a spec with no markers has a 0x1 "table" whose only
        // column is `features` (`de p001 0x1`). Porting that means the R shim has to build
        // a 0-row, 1-column frame, not an empty one with no columns.
        // Row names are part of the S4 object, and `identical()` checks them. They are the
        // feature names -- so for `fc005`, where the rows are ordered by p-value per group,
        // the row names are G11,G21,G17 rather than ascending. `de_rownames` is what makes
        // that machine-checked instead of an observation.
        let want_rownames = golden("de_rownames", spec.name).unwrap_or_default();
        assert_eq!(
            want_features, want_rownames,
            "{}: upstream's row names are its `features` column, in row order",
            spec.name
        );
        assert_eq!(
            want_rows > 0,
            want_cols.len() > 1,
            "{}: schema {cols_line:?}",
            spec.name
        );
        if want_rows > 0 {
            assert_eq!(
                want_cols.join(","),
                "clusters,features,pvalues,logFC,pct.1,pct.2,pvalues.adj",
                "{}: schema",
                spec.name
            );
        } else {
            assert_eq!(want_cols, ["features"], "{}: collapsed schema", spec.name);
        }
        for r in 0..got.len() {
            for cname in &want_cols {
                let v = match &cname[..] {
                    "clusters" => Some(got[r].cluster.clone()),
                    "features" => Some(got[r].feature.clone()),
                    _ => None,
                };
                if let Some(v) = v {
                    assert_eq!(
                        v.as_str(),
                        want_str_col(spec.name, cname, r),
                        "{} row {r} {cname}",
                        spec.name
                    );
                    continue;
                }
                let vals = parse_vec(fields(de_col_line(spec.name, cname))[3]);
                let want = vals[r];
                let got = match &cname[..] {
                    "pvalues" => got[r].pvalue,
                    "logFC" => got[r].log_fc,
                    "pct.1" => got[r].pct_1,
                    "pct.2" => got[r].pct_2,
                    "pvalues.adj" => got[r].pvalue_adj,
                    other => panic!("unhandled column {other}"),
                };
                bits_eq(got, want, &format!("{} row {r} {cname}", spec.name));
            }
        }
    }
}

/// A *character* column of the marker table, which is comma-joined with no quoting.
fn want_str_col(spec: &str, col: &str, row: usize) -> String {
    fields(de_col_line(spec, col))[3]
        .split(',')
        .nth(row)
        .expect("row")
        .to_string()
}

fn de_col_line(spec: &str, col: &str) -> &'static str {
    GOLDEN
        .lines()
        .find(|l| l.starts_with(&format!("de_col\t{spec}\t{col}\t")))
        .unwrap_or_else(|| panic!("no de_col {spec} {col}"))
}
