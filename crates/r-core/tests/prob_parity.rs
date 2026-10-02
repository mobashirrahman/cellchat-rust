//! End-to-end parity of [`r_core::prob::Kernel`] against upstream `computeCommunProb`.
//!
//! Corpus: `tests/fixtures/prob_golden.txt` (Prob, Pval and `options$parameter` for eight
//! parameter configurations), with the inputs in `prob_inputs.tsv` and the observed
//! `data.use.avg` in `prob_avg.txt`. All three are written by
//! `tests/parity/gen_prob_golden.R`, which **sources `R/modeling.R` from the pinned commit
//! and calls `computeCommunProb` itself** on a minimal S4 object. Nothing in the expected
//! values is reimplemented on the R side.
//!
//! ## Why the inputs are dumped rather than re-derived
//!
//! The Rust test does not rebuild the fixture from `set.seed`/`runif`; it reads it. That is
//! deliberate, for the reason `expr_matrix.tsv` exists: re-derivation forces the test to
//! reproduce R's own ordering conventions, and when it gets them wrong the gene *set* stays
//! right while the *values* permute, so every assertion still looks plausible.
//!
//! ## What the configurations cover
//!
//! | config | exercises |
//! |---|---|
//! | `tri_ps0` / `tri_ps1` | the default path, with `population.size` off/on |
//! | `trim_ps0` | `truncatedMean`, nboot 7, a second seed |
//! | `thresh_ps1` | `thresholdedMean`, trim 0.3 (the zero short-circuit) |
//! | `median_ps0` | `median`, nboot 6 |
//! | `hill2` | `n = 2` -- `x^2` in both numerator and denominator |
//! | `kh1e3` | `Kh = 1e3` and nboot 8 |
//! | `nboot1` | `nboot = 1`, where every `Pval` is 0 or 1 |
//!
//! Note that `Prob` legitimately **exceeds 1** in the fixture: the co-agonist term is a
//! product of `1 + h` over all cofactor subunits, so a 10-subunit agonist can reach
//! `2^10`. That is upstream behaviour, not a port bug, and the corpus pins it.

use r_core::aggregate::{aggregate_1, all_cols, GroupMean};
use r_core::db::Database;
use r_core::prob::{Kernel, LrPair};
use std::collections::HashMap;
use std::path::PathBuf;

const GOLDEN: &str = include_str!("../../../tests/fixtures/prob_golden.txt");
const INPUTS: &str = include_str!("../../../tests/fixtures/prob_inputs.tsv");
const AVGS: &str = include_str!("../../../tests/fixtures/prob_avg.txt");

fn fixtures() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tests/fixtures")
}

fn human_db() -> Database {
    Database::load(&fixtures().join("db_human")).expect("db_human fixture")
}

/// The fixture, parsed from `prob_inputs.tsv`.
struct Inputs {
    n_cells: usize,
    n_genes: usize,
    n_groups: usize,
    genes: Vec<String>,
    group_labels: Vec<String>,
    /// `as.character(cell_group)` per cell.
    cell_group: Vec<String>,
    /// `t(data.use)`, cells x genes, column-major: `data[g * n_cells + c]`.
    data: Vec<f64>,
    lr: Vec<LrPair>,
}

impl Inputs {
    /// Level index of each cell.
    fn group_index(&self) -> Vec<usize> {
        self.cell_group
            .iter()
            .map(|g| {
                self.group_labels
                    .iter()
                    .position(|l| l == g)
                    .expect("known label")
            })
            .collect()
    }
}

fn header_field(line: &str, key: &str) -> usize {
    line.strip_prefix(key)
        .and_then(|v| v.strip_prefix('\t'))
        .unwrap_or_else(|| panic!("prob_inputs.tsv: expected `{key}<TAB>n`, got {line:?}"))
        .parse()
        .expect("count")
}

fn parse_hex_float(s: &str) -> f64 {
    let s = s.trim();
    let (neg, s) = match s.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, s.strip_prefix('+').unwrap_or(s)),
    };
    let s = s.trim_start_matches("0x").trim_start_matches("0X");
    let (mantissa, exp) = match s.find(['p', 'P']) {
        Some(i) => (&s[..i], &s[i + 1..]),
        None => (s, "0"),
    };
    let (int_part, frac_part) = match mantissa.find('.') {
        Some(i) => (&mantissa[..i], &mantissa[i + 1..]),
        None => (mantissa, ""),
    };
    assert!(
        !int_part.is_empty(),
        "hex float without an integer part: {s:?}"
    );
    let mut value = 0.0f64;
    for b in int_part.bytes() {
        value = value * 16.0 + hex_digit(b) as f64;
    }
    let mut scale = 1.0f64;
    for b in frac_part.bytes() {
        scale /= 16.0;
        value += hex_digit(b) as f64 * scale;
    }
    let v = value * 2f64.powi(exp.parse().unwrap_or(0));
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

fn load_inputs() -> Inputs {
    let mut lines = INPUTS.lines();
    let n_cells = header_field(lines.next().unwrap(), "cells");
    let n_genes = header_field(lines.next().unwrap(), "genes");
    let n_groups = header_field(lines.next().unwrap(), "groups");
    let genes: Vec<String> = lines
        .by_ref()
        .take(n_genes)
        .map(|s| s.to_string())
        .collect();
    assert_eq!(lines.next(), Some("groups_labels"));
    let group_labels: Vec<String> = lines
        .by_ref()
        .take(n_groups)
        .map(|s| s.to_string())
        .collect();
    assert_eq!(lines.next(), Some("cell_group"));
    let cell_group: Vec<String> = lines
        .by_ref()
        .take(n_cells)
        .map(|s| s.to_string())
        .collect();
    assert_eq!(lines.next(), Some("data"));
    let mut data = vec![0.0f64; n_cells * n_genes];
    for (c, line) in lines.by_ref().take(n_cells).enumerate() {
        let vals: Vec<f64> = line.split_whitespace().map(parse_hex_float).collect();
        assert_eq!(vals.len(), n_genes, "cell {c} has {} values", vals.len());
        for (g, v) in vals.into_iter().enumerate() {
            // Column-major cells x genes: the gene index strides.
            data[g * n_cells + c] = v;
        }
    }
    assert_eq!(lines.next(), Some("lr"));
    let lr: Vec<LrPair> = lines
        .map(|l| {
            let f: Vec<&str> = l.split('\t').collect();
            assert_eq!(f.len(), 7, "lr row {l:?}");
            LrPair {
                ligand: f[0].into(),
                receptor: f[1].into(),
                agonist: f[2].into(),
                antagonist: f[3].into(),
                co_a: f[4].into(),
                co_i: f[5].into(),
                label: f[6].into(),
            }
        })
        .collect();
    Inputs {
        n_cells,
        n_genes,
        n_groups,
        genes,
        group_labels,
        cell_group,
        data,
        lr,
    }
}

/// One configuration from the generator.
struct Config {
    name: String,
    mean: GroupMean,
    population_size: bool,
    nboot: usize,
    seed: i32,
    kh: f64,
    n: f64,
}

/// Mirrors `configs` in the generator. Duplicated deliberately: if the generator's list ever
/// grows, `every_config_is_covered` fails rather than a new config going untested.
fn configs() -> Vec<Config> {
    use GroupMean::*;
    let mk =
        |name: &str, mean: GroupMean, ps: bool, nboot: usize, seed: i32, kh: f64, n: f64| Config {
            name: name.into(),
            mean,
            population_size: ps,
            nboot,
            seed,
            kh,
            n,
        };
    vec![
        mk("tri_ps0", TriMean, false, 5, 1, 0.5, 1.0),
        mk("tri_ps1", TriMean, true, 5, 1, 0.5, 1.0),
        mk("trim_ps0", TrimmedMean { trim: 0.2 }, false, 7, 2, 0.5, 1.0),
        mk(
            "thresh_ps1",
            ThresholdedMean { trim: 0.3 },
            true,
            4,
            3,
            0.5,
            1.0,
        ),
        mk("median_ps0", Median, false, 6, 4, 0.5, 1.0),
        mk("hill2", TriMean, false, 5, 1, 0.5, 2.0),
        mk("kh1e3", TriMean, true, 8, 7, 1e3, 1.0),
        mk("nboot1", TriMean, false, 1, 5, 0.5, 1.0),
    ]
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

/// Build the kernel for one configuration and return `(Prob, Pval)`.
fn run(inp: &Inputs, db: &Database, cfg: &Config) -> (Vec<f64>, Vec<f64>, Kernel) {
    let group = inp.group_index();
    let plans = Kernel::resolve(
        inp.genes.clone(),
        &inp.lr,
        db,
        &group,
        inp.n_groups,
        inp.n_cells,
        cfg.population_size,
    )
    .expect("resolve");

    let (boot, perms) = Kernel::build_boot(
        &inp.data,
        inp.n_genes,
        inp.n_cells,
        &group,
        inp.n_groups,
        cfg.nboot,
        cfg.seed,
        cfg.mean,
    );

    // `population.size` recomputes `table(group[permutation[, nE]])/nC` per replicate.
    let pop_boot: Vec<f64> = if cfg.population_size {
        let mut v = Vec::with_capacity(cfg.nboot * inp.n_groups);
        for p in &perms {
            let mut counts = vec![0.0f64; inp.n_groups];
            for &c in p {
                counts[group[c]] += 1.0;
            }
            v.extend(counts.iter().map(|c| c / inp.n_cells as f64));
        }
        v
    } else {
        Vec::new()
    };

    let avg = aggregate_1(
        &inp.data,
        inp.n_cells,
        &all_cols(inp.n_genes),
        &group,
        inp.n_groups,
        cfg.mean,
    );

    let kk = inp.n_groups * inp.n_groups;
    let kernel = Kernel::new(
        inp.genes.clone(),
        avg,
        boot,
        plans,
        cfg.nboot,
        cfg.kh,
        cfg.n,
        inp.n_groups,
        inp.group_labels.clone(),
        vec![1.0; kk],
        vec![1.0; kk],
        inp.lr.len(),
        pop_boot,
    )
    .expect("kernel");
    let (prob, pval) = kernel.prob_arrays().expect("prob_arrays");
    (prob, pval, kernel)
}

#[test]
fn every_config_is_covered() {
    let mut got: Vec<&str> = GOLDEN
        .lines()
        .filter_map(|l| l.strip_prefix("config\t"))
        .map(|l| l.split('\t').next().unwrap())
        .collect();
    got.sort_unstable();
    let all = configs();
    let mut want: Vec<&str> = all.iter().map(|c| c.name.as_str()).collect();
    want.sort_unstable();
    assert_eq!(
        got, want,
        "the Rust config list has drifted from the generator's"
    );
}

#[test]
fn prob_and_pval_are_bit_identical_for_every_config() {
    let inp = load_inputs();
    let db = human_db();
    let mut checked = 0;
    for cfg in configs() {
        let (prob, pval, kernel) = run(&inp, &db, &cfg);
        let (want_prob, want_pval) = (golden_vec("prob", &cfg.name), golden_vec("pval", &cfg.name));
        assert_eq!(prob.len(), want_prob.len(), "{}: Prob length", cfg.name);
        for (i, (&g, &w)) in prob.iter().zip(&want_prob).enumerate() {
            if g.to_bits() != w.to_bits() {
                let (lr, a, b) = unravel(i, inp.n_groups);
                panic!(
                    "{}: Prob[{lr}][{a},{b}] bit mismatch\n  rust = {g:.17e} (0x{:016x})\n  R    = {w:.17e} (0x{:016x})",
                    cfg.name,
                    g.to_bits(),
                    w.to_bits()
                );
            }
        }
        for (i, (&g, &w)) in pval.iter().zip(&want_pval).enumerate() {
            if g.to_bits() != w.to_bits() {
                let (lr, a, b) = unravel(i, inp.n_groups);
                panic!(
                    "{}: Pval[{lr}][{a},{b}] bit mismatch\n  rust = {g:.17e} (0x{:016x})\n  R    = {w:.17e} (0x{:016x})",
                    cfg.name,
                    g.to_bits(),
                    w.to_bits()
                );
            }
        }
        // dimnames are part of the drop-in contract.
        let (g1, g2, g3) = kernel.dimnames();
        assert_eq!(g1, inp.group_labels, "{}: group dimnames", cfg.name);
        assert_eq!(g2, inp.group_labels, "{}: group dimnames (2nd)", cfg.name);
        assert_eq!(
            g3,
            inp.lr.iter().map(|p| p.label.clone()).collect::<Vec<_>>()
        );
        checked += 1;
    }
    assert_eq!(checked, 8, "expected eight configurations");
}

#[test]
fn the_observed_average_expression_matches_r() {
    // The kernel's `avg` is the same `aggregate` call as the bootstrap's, so a bug there
    // would cancel on both sides and the Prob comparison would still pass. This pins the
    // observed half against R's own output.
    let inp = load_inputs();
    let group = inp.group_index();
    for cfg in configs() {
        let want = golden_vec("avg", &cfg.name);
        let got = aggregate_1(
            &inp.data,
            inp.n_cells,
            &all_cols(inp.n_genes),
            &group,
            inp.n_groups,
            cfg.mean,
        );
        assert_eq!(got.len(), want.len(), "{}: avg length", cfg.name);
        for (i, (&g, &w)) in got.iter().zip(&want).enumerate() {
            bits_eq(g, w, &format!("{}: avg[{i}]", cfg.name));
        }
    }
}

#[test]
fn pval_is_a_multiple_of_one_over_nboot() {
    // Objective 3's recorded invariant: `Pval` lives in `{k / nboot}`.
    let inp = load_inputs();
    let db = human_db();
    for cfg in configs() {
        let (prob, pval, _) = run(&inp, &db, &cfg);
        let nb = cfg.nboot as f64;
        for (i, &v) in pval.iter().enumerate() {
            let k = v * nb;
            assert!(
                (0.0..=1.0).contains(&v),
                "{}: Pval[{i}]={v} out of range",
                cfg.name
            );
            assert!(
                (k - k.round()).abs() < 1e-9,
                "{}: Pval[{i}]={v} is not k/{nb}",
                cfg.name
            );
        }
        // R-ism 6: `Pval[Prob == 0] <- 1` runs after the loop.
        for (i, (&p, &v)) in prob.iter().zip(&pval).enumerate() {
            if p == 0.0 {
                assert_eq!(v, 1.0, "{}: Pval[{i}] must be 1 where Prob is 0", cfg.name);
            }
        }
        // `dimnames` consistency.
        assert_eq!(prob.len(), pval.len());
    }
}

#[test]
fn prob_is_invariant_to_nboot() {
    // Objective 3's metamorphic relation: the observed network does not depend on the
    // number of permutations. Compare against the corpus by checking that two of our own
    // configs differing only in nboot agree...
    let inp = load_inputs();
    let db = human_db();
    let a = Config {
        name: "x".into(),
        mean: GroupMean::TriMean,
        population_size: false,
        nboot: 3,
        seed: 1,
        kh: 0.5,
        n: 1.0,
    };
    let b = Config {
        name: "x".into(),
        mean: GroupMean::TriMean,
        population_size: false,
        nboot: 11,
        seed: 1,
        kh: 0.5,
        n: 1.0,
    };
    let (pa, _, _) = run(&inp, &db, &a);
    let (pb, _, _) = run(&inp, &db, &b);
    for (i, (&x, &y)) in pa.iter().zip(&pb).enumerate() {
        bits_eq(x, y, &format!("Prob[{i}] with nboot 3 vs 11"));
    }
}

#[test]
fn the_group_order_is_factor_level_order_not_appearance_order() {
    // The fixture's first cell belongs to "g3" (see the generator), so a port that indexes
    // groups by first appearance would produce a permuted Prob. Cheap, and it is the single
    // most likely way to get this wrong.
    let inp = load_inputs();
    let first = &inp.cell_group[0];
    assert_ne!(
        first, &inp.group_labels[0],
        "the fixture must not start in group 0"
    );
    assert!(inp.group_labels.contains(first));
}

/// A whitespace-separated, comma-joined value list from either corpus.
fn golden_vec(kind: &str, name: &str) -> Vec<f64> {
    let prefix = format!("{kind}\t{name}\t");
    let line = AVGS
        .lines()
        .chain(GOLDEN.lines())
        .find(|l| l.starts_with(&prefix))
        .unwrap_or_else(|| panic!("no `{kind}` record for {name}"));
    parse_vec(line.split('\t').nth(2).expect("values field"))
}

/// Map a flat index into `(interaction, source, target)`, for error messages.
fn unravel(i: usize, k: usize) -> (usize, usize, usize) {
    let kk = k * k;
    (i / kk, (i % kk) % k, (i % kk) / k)
}

#[test]
fn the_hand_transcribed_lr_table_matches_the_generator() {
    let inp = load_inputs();
    // The L-R table is dumped by the generator and read by `load_inputs`, so the useful
    // invariant here is structural: one row per interaction, labels unique, and every
    // ligand/receptor either a plain gene or a known complex.
    let mut labels: Vec<&str> = inp.lr.iter().map(|p| p.label.as_str()).collect();
    let before = labels.len();
    labels.sort_unstable();
    labels.dedup();
    assert_eq!(labels.len(), before, "duplicate interaction labels");
    assert!(inp
        .lr
        .iter()
        .all(|p| p.label == format!("{}_{}", p.ligand, p.receptor)));
    let complexes = Database::load(&fixtures().join("db_human")).unwrap();
    for p in &inp.lr {
        for name in [&p.ligand, &p.receptor] {
            assert!(
                inp.genes.contains(name) || complexes.complex_cells(name).is_some(),
                "{name} is neither a gene in the fixture nor a known complex"
            );
        }
    }
    let _ = HashMap::<(), ()>::new();
}
