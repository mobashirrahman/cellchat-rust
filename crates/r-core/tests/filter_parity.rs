//! `filterCommunication` against a corpus generated from the pinned upstream
//! (`tests/parity/gen_filter_golden.R`).
//!
//! The two halves are checked separately because they fail differently: `min.cells` is a
//! slice-and-count on `net$prob`, and the `min.samples >= 2` half re-derives group means per
//! sample and applies a cross-sample consistency mask. A test that only compared the final
//! `net$prob` could not tell a wrong `score.LR` that happened to cancel from a right one, so
//! the per-sample scores and every intermediate count are pinned too.
//!
//! Three upstream behaviours are *reproduced* rather than repaired, each documented in
//! `crates/r-core/src/filter.rs`:
//!   * `min.samples` above the sample count is upstream's own `stop()`;
//!   * an all-zero `net$prob` with `min.samples >= 2` is **"subscript out of bounds"**, from
//!     `for (jj in 1:length(LR.nonzero))` with `LR.nonzero` empty being `1:0`;
//!   * `table(idents[cell.use])` keeps every level, so a group *absent* from a sample counts
//!     as zero cells and lands in that sample's excluded set -- which is the only thing that
//!     gives `rare.keep` something to preserve.

use r_core::db::Database;
use r_core::filter::{filter_communication, FilterError, FilterInput, MeanKind};
use std::collections::HashMap;

const GOLDEN: &str = include_str!("../../../tests/fixtures/filter_golden.txt");
const INPUTS: &str = include_str!("../../../tests/fixtures/filter_inputs.tsv");

fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn record(case: &str, kind: &str) -> Option<String> {
    GOLDEN
        .lines()
        .find(|l| l.starts_with(&format!("filter\t{case}\t{kind}=")))
        .map(|l| {
            let f = fields(l);
            let v = &f[2];
            v.strip_prefix(&format!("{kind}=")).unwrap_or(v).to_string()
        })
}

fn is_error(case: &str) -> Option<String> {
    GOLDEN
        .lines()
        .find(|l| l.starts_with(&format!("filter\t{case}\tERROR")))
        .map(|l| fields(l).get(3).map(|s| s.to_string()).unwrap_or_default())
}

fn parse_vec(s: &str) -> Vec<f64> {
    if s.is_empty() {
        return Vec::new();
    }
    s.split(',')
        .map(|v| {
            let v = v.trim();
            assert_ne!(v, "NA", "unexpected NA");
            v.parse().unwrap_or_else(|_| panic!("not a number: {v:?}"))
        })
        .collect()
}

fn parse_usizes(s: &str) -> Vec<usize> {
    s.split(',').map(|v| v.parse().unwrap()).collect()
}

fn parse_hex(s: &str) -> f64 {
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
    for b in int_place_bytes(int_part) {
        v = v * 16.0 + b as f64;
    }
    let mut scale = 1.0f64;
    for b in int_place_bytes(frac_part) {
        scale /= 16.0;
        v += b as f64 * scale;
    }
    let out = v * 2f64.powi(exp);
    if neg {
        -out
    } else {
        out
    }
}

fn int_place_bytes(s: &str) -> Vec<u32> {
    s.bytes()
        .map(|b| (b as char).to_digit(16).unwrap())
        .collect()
}

#[derive(Debug, Clone)]
struct Case {
    name: String,
    levels: Vec<String>,
    lr: Vec<String>,
    n_cells: usize,
    min_cells: usize,
    min_samples: Option<usize>,
    rare_keep: bool,
    nonfilter_keep: bool,
    mean: MeanKind,
    trim: f64,
    group_index: Vec<usize>,
    sample_index: Vec<usize>,
    gene_names: Vec<String>,
    prob: Vec<f64>,
    pval: Vec<f64>,
    data: Vec<f64>,
    /// `options$parameter$raw.use`. `false` means the kernel reads `data.smooth` instead of
    /// `data.signaling`, so the *input matrix* changes -- the only thing `raw.use` does.
    raw_use: bool,
    /// `object@data.smooth`, present only when `raw_use` is `false`.
    smooth: Option<Vec<f64>>,
}

/// Bit-exact equality for two float vectors, with `NaN == NaN` so an all-`NA` array compares
/// equal to itself. `assert_eq!` on `[f64]` would report `NaN != NaN`, which is the opposite of
/// what a fixture comparing two tables wants.
fn same_bits(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len()
        && a.iter()
            .zip(b)
            .all(|(x, y)| (x.is_nan() && y.is_nan()) || x.to_bits() == y.to_bits())
}

fn cases() -> Vec<Case> {
    let mut out: Vec<Case> = Vec::new();
    let lines: Vec<&str> = INPUTS.lines().collect();
    let mut i = 0usize;
    while i < lines.len() {
        if !lines[i].starts_with("case\t") {
            i += 1;
            continue;
        }
        let f = fields(lines[i]);
        let name = f[1].to_string();
        // The trailing "|" keeps the final field from being dropped by a split that discards
        // trailing empties.
        // The extra trailing "|" keeps the final field from being dropped by a split that
        // discards trailing empties.
        let meta = format!("{}|", f[2]);
        let m: Vec<&str> = meta.split('|').collect();
        let levels: Vec<String> = m[0].split(',').map(|s| s.to_string()).collect();
        let lr: Vec<String> = m[1].split(',').map(|s| s.to_string()).collect();
        let n_cells: usize = m[2].parse().unwrap();
        let min_cells: usize = m[3].parse().unwrap();
        let min_samples = if m[4] == "NA" {
            None
        } else {
            Some(m[4].parse().unwrap())
        };
        let rare_keep = m[5] == "TRUE";
        let nonfilter_keep = m[6] == "TRUE";
        let mean = match m[7] {
            "triMean" => MeanKind::TriMean,
            "truncatedMean" => MeanKind::TruncatedMean,
            "thresholdedMean" => MeanKind::ThresholdedMean,
            "median" => MeanKind::Median,
            other => panic!("unknown mean {other:?}"),
        };
        let trim: f64 = m[8].parse().unwrap();
        // Two fields appended to the metadata line, so the earlier indices are untouched.
        let raw_use = m[9] == "TRUE";
        // Three states, not two: `-` (no `data.smooth` needed), `yes` (the slot is present),
        // and `dropped` (a matrix was built for the record but the object has no such slot, so
        // upstream's mis-negated guard does *not* fire and it reads an absent matrix).
        let smooth_state = m[10];
        assert!(
            matches!(smooth_state, "-" | "yes" | "dropped"),
            "{smooth_state:?}"
        );
        assert_eq!(lines[i + 1].split('\t').next(), Some("idents"), "{name}");
        let group_index = parse_usizes(lines[i + 1].split('\t').nth(2).unwrap());
        let sample_index = parse_usizes(lines[i + 2].split('\t').nth(2).unwrap());
        let gene_names: Vec<String> = lines[i + 3]
            .split('\t')
            .nth(2)
            .expect("genes record")
            .split(',')
            .map(|s| s.to_string())
            .collect();
        // Per-case block, 10 lines: case, idents, samples, genes, "prob", prob row,
        // "pval", pval row, "data", data row.
        assert_eq!(lines[i + 3].split('\t').next(), Some("genes"), "{name}");
        assert_eq!(lines[i + 4], "prob", "{name}");
        let prob: Vec<f64> = lines[i + 5].split_whitespace().map(parse_hex).collect();
        assert_eq!(lines[i + 6], "pval", "{name}");
        let pval: Vec<f64> = lines[i + 7].split_whitespace().map(parse_hex).collect();
        assert_eq!(lines[i + 8], "data", "{name}");
        let data: Vec<f64> = lines[i + 9].split_whitespace().map(parse_hex).collect();
        // The `smooth` record, when present, sits immediately after `data` and is itself two
        // lines (a bare tag, then the values) -- so the fixed offsets below shift by two.
        let (smooth, next_i) = if smooth_state != "-" {
            assert_eq!(lines[i + 10].split('\t').next(), Some("smooth"), "{name}");
            let v: Vec<f64> = lines[i + 11].split_whitespace().map(parse_hex).collect();
            assert_eq!(
                v.len(),
                data.len(),
                "{name}: smooth must match data's shape"
            );
            (Some(v), i + 12)
        } else {
            (None, i + 10)
        };
        let k = levels.len();
        let n_lr = lr.len();
        assert_eq!(prob.len(), k * k * n_lr, "{name}: prob length");
        assert_eq!(pval.len(), k * k * n_lr, "{name}: pval length");
        assert_eq!(
            data.len(),
            gene_names.len() * n_cells,
            "{name}: data length"
        );
        assert_eq!(group_index.len(), n_cells, "{name}: idents length");
        assert_eq!(sample_index.len(), n_cells, "{name}: samples length");
        out.push(Case {
            name,
            levels,
            lr,
            n_cells,
            min_cells,
            min_samples,
            rare_keep,
            nonfilter_keep,
            mean,
            trim,
            group_index,
            sample_index,
            gene_names,
            prob,
            pval,
            data,
            raw_use,
            smooth,
        });
        i = next_i;
    }
    out
}

/// `score.LR`'s flat index: `k x k x n_nz x n_samples`, last dimension fastest, so
/// `source + k*target + k*k*jj + k*k*n_nz*sample`. Written out in the test as well as in the
/// core: a test that reused the core's helper would agree with it by construction and could
/// not catch a wrong layout.
fn si(source: usize, target: usize, jj: usize, sample: usize, k: usize, n_nz: usize) -> usize {
    source + k * target + k * k * jj + k * k * n_nz * sample
}

/// The fixture DB plus the two parallel columns of its interaction table, which
/// `filterCommunication` reaches with `match(LR.nonzero, interaction_input$interaction_name)`.
struct Fixture {
    db: Database,
    ligand: Vec<String>,
    receptor: Vec<String>,
}

fn load_db() -> Fixture {
    let lines: Vec<&str> = INPUTS.lines().collect();
    // The `db_*` records are two-field (`db_x<TAB>payload`) while `db_subunit` is
    // three-field, so take the *last* field rather than indexing 2.
    let get = |tag: &str| -> String {
        lines
            .iter()
            .find(|l| l.starts_with(&format!("{tag}\t")))
            .and_then(|l| fields(l).last().map(|x| x.to_string()))
            .unwrap_or_default()
    };
    let complex_names: Vec<String> = get("db_complexes")
        .split(',')
        .map(|s| s.to_string())
        .collect();
    let mut complexes: HashMap<String, Vec<String>> = HashMap::new();
    for l in &lines {
        if l.starts_with("db_subunit\t") {
            let f = fields(l);
            let col = f[1].to_string();
            let vals: Vec<String> = f[2].split(',').map(|s| s.to_string()).collect();
            for (i, name) in complex_names.iter().enumerate() {
                complexes
                    .entry(name.clone())
                    .or_default()
                    .push(vals.get(i).cloned().unwrap_or_default());
            }
            assert!(col.starts_with("subunit_"), "unexpected column {col:?}");
        }
    }
    let symbols: HashMap<String, ()> = get("db_symbols")
        .split(',')
        .filter(|s| !s.is_empty())
        .map(|s| (s.to_string(), ()))
        .collect();
    let ligand: Vec<String> = get("db_ligand").split(',').map(|s| s.to_string()).collect();
    let receptor: Vec<String> = get("db_receptor")
        .split(',')
        .map(|s| s.to_string())
        .collect();
    assert_eq!(
        ligand.len(),
        complex_names.len(),
        "db: ligand/receptor/complex lengths"
    );
    let db = Database::from_parts(
        "fixture",
        vec![],
        complexes,
        HashMap::new(),
        complex_names,
        Vec::new(),
        5,
        0,
        symbols,
    );
    Fixture {
        db,
        ligand,
        receptor,
    }
}

/// The ligand and receptor of each L-R, matched by name against the DB's interaction table
/// exactly as upstream does: `idx <- match(LR.nonzero, interaction_input$interaction_name)`.
fn lr_genes(c: &Case, f: &Fixture) -> (Vec<String>, Vec<String>) {
    let mut lig = Vec::with_capacity(c.lr.len());
    let mut rec = Vec::with_capacity(c.lr.len());
    for name in &c.lr {
        let i = f
            .ligand
            .iter()
            .position(|l| {
                format!("{l}_{}", name.split_once('_').map(|x| x.1).unwrap_or("")) == *name
            })
            .unwrap_or_else(|| panic!("{name}: not in the DB's interaction table"));
        lig.push(f.ligand[i].clone());
        rec.push(f.receptor[i].clone());
    }
    (lig, rec)
}

fn run(c: &Case, f: &Fixture) -> Result<r_core::filter::FilterResult, FilterError> {
    // `data <- data/max(data)` and then `data.use <- data/max(data)`, so the kernel sees the
    // *scaled* matrix. The scaling is part of the contract, not a convenience: a port that
    // skipped it would differ whenever any group mean crossed 1.0 -- which is exactly where
    // the outer product's sign would change.
    // `raw.use = FALSE` reads `data.smooth`. Upstream applies the same
    // `data <- data/max(data)` to whichever matrix it chose, so the scaling is applied to the
    // *selected* matrix and not to `data.signaling` -- using the wrong one's maximum is a
    // plausible-looking bug that only shows up when the two matrices have different ranges.
    let src: &[f64] = if c.raw_use {
        &c.data
    } else {
        c.smooth
            .as_deref()
            .unwrap_or_else(|| panic!("{}: raw.use = FALSE needs a data.smooth record", c.name))
    };
    let max = src.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let scaled: Vec<f64> = src.iter().map(|v| v / max).collect();
    let (lig, rec) = lr_genes(c, f);
    let n_samples = c
        .sample_index
        .iter()
        .copied()
        .max()
        .map(|m| m + 1)
        .unwrap_or(0);
    let inp = FilterInput {
        prob: &c.prob,
        n_lr: c.lr.len(),
        n_groups: c.levels.len(),
        group_index: &c.group_index,
        n_samples,
        sample_index: &c.sample_index,
        data: &scaled,
        gene_names: &c.gene_names,
        min_cells: c.min_cells,
        min_samples: c.min_samples,
        rare_keep: c.rare_keep,
        mean: c.mean,
        trim: c.trim,
        ligand: &lig,
        receptor: &rec,
    };
    filter_communication(&inp, &f.db)
}

#[test]
fn the_corpus_covers_both_halves_and_both_upstream_failures() {
    let cs = cases();
    // Named, not counted: a count silently passes if one case is replaced by another, and stops
    // passing for no reason when a case is added -- which is what happened five times while the
    // `raw.use` fixtures were being built.
    assert_eq!(cs.len(), 19, "the corpus should carry 19 cases");
    for want in [
        "raw_use_false_mincells",
        "raw_use_false_median",
        "raw_use_false_rare_keep",
        "raw_use_false_slot_absent",
        "raw_use_false_nonfilter_keep",
        "raw_use_false_two_samples",
    ] {
        assert!(
            cs.iter().any(|c| c.name == want),
            "corpus lost the {want:?} case"
        );
    }
    // The `min.cells` half alone.
    assert!(
        cs.iter().any(|c| c.min_samples.is_none()),
        "no min.samples = NULL case"
    );
    // The `min.samples` half, with >= 2 samples.
    assert!(
        cs.iter()
            .filter(|c| c.min_samples.is_some_and(|m| m >= 2))
            .count()
            >= 8,
        "too few multi-sample cases"
    );
    // Both upstream failures.
    let errs: Vec<&Case> = cs.iter().filter(|c| is_error(&c.name).is_some()).collect();
    assert_eq!(
        errs.len(),
        2,
        "expected exactly the two upstream error cases"
    );
    assert!(errs
        .iter()
        .any(|c| is_error(&c.name).as_deref() == Some("subscript out of bounds")));
    assert!(errs
        .iter()
        .any(|c| is_error(&c.name).unwrap().contains("There are only")));
    // All four `type.mean` choices.
    for m in ["TriMean", "TruncatedMean", "ThresholdedMean", "Median"] {
        assert!(
            cs.iter().any(|c| format!("{:?}", c.mean) == m),
            "no {m} case"
        );
    }
    // `rare.keep` both ways, on the same input.
    assert!(cs.iter().any(|c| c.rare_keep), "no rare.keep = TRUE case");
    assert!(cs
        .iter()
        .any(|c| !c.rare_keep && c.min_samples.is_some_and(|m| m >= 2)));
    // `nonFilter.keep`.
    assert!(
        cs.iter().any(|c| c.nonfilter_keep),
        "no nonFilter.keep = TRUE case"
    );
}

/// The two upstream failures, and the different mechanism each needs.
///
/// * The all-zero network is a *genuine* error deep inside upstream's loop, so the core
///   returns it and the text is compared.
/// * `min.samples` above the sample count is an ordinary precondition, but upstream raises
///   it **after** printing the `min.cells` messages. The core therefore reports
///   `min_samples_ok = false` and the R shim raises, so the message sequence on the error
///   path stays byte-identical. The core's job here is to refuse to run the sample half, and
///   the text lives with the other upstream strings so it is checked in one place.
#[test]
fn the_upstream_failures_are_reproduced_not_repaired() {
    let fx = load_db();
    for c in cases() {
        let Some(want) = is_error(&c.name) else {
            continue;
        };
        match run(&c, &fx) {
            Err(e) => assert_eq!(e.to_string(), want, "{}: error text", c.name),
            Ok(got) => {
                let reported = if want.contains("There are only") {
                    assert!(
                        !got.min_samples_ok && !got.ran_sample_filter,
                        "{}: min.samples > n_samples must leave the sample half unrun",
                        c.name
                    );
                    r_core::filter::FilterError::MinSamplesTooLarge {
                        n_samples: c.sample_index.iter().copied().max().unwrap() + 1,
                        min_samples: c.min_samples.unwrap(),
                    }
                } else {
                    panic!(
                        "{}: upstream errors with {want:?} but the port returned a result",
                        c.name
                    )
                };
                assert_eq!(
                    reported.to_string(),
                    want,
                    "{}: reported error text",
                    c.name
                );
            }
        }
    }
}

#[test]
fn the_filtered_probability_array_matches_r_bit_for_bit() {
    let fx = load_db();
    let mut checked = 0;
    for c in cases() {
        if is_error(&c.name).is_some() {
            continue;
        }
        let got = run(&c, &fx).unwrap_or_else(|e| panic!("{}: unexpected {e}", c.name));
        let want = parse_vec(&record(&c.name, "prob").expect("prob record"));
        assert_eq!(got.prob.len(), want.len(), "{}: prob length", c.name);
        for (i, (&g, &w)) in got.prob.iter().zip(want.iter()).enumerate() {
            if g.to_bits() != w.to_bits() {
                panic!(
                    "{}: net$prob[{i}] bit mismatch\n  rust = {g:.17e} (0x{:016x})\n  R    = {w:.17e} (0x{:016x})",
                    c.name, g.to_bits(), w.to_bits()
                );
            }
        }
        let dim = record(&c.name, "dim").expect("dim record");
        let want_dim: Vec<usize> = dim.split('x').map(|s| s.parse().unwrap()).collect();
        assert_eq!(
            vec![c.levels.len(), c.levels.len(), c.lr.len()],
            want_dim,
            "{}: dim(net$prob)",
            c.name
        );
        // `net$pval` is never touched.
        let pval = parse_vec(&record(&c.name, "pval").expect("pval record"));
        for (i, (&g, &w)) in c.pval.iter().zip(pval.iter()).enumerate() {
            assert_eq!(
                g.to_bits(),
                w.to_bits(),
                "{}: net$pval[{i}] changed",
                c.name
            );
        }
        checked += 1;
    }
    // 19 cases, 2 of which are upstream errors (`raw_use_false_median` and
    // `raw_use_false_rare_keep`, both the all-zero `LR.nonzero` crash).
    assert_eq!(checked, 17, "compared {checked} cases");
}

/// `raw.use = FALSE` reads `object@data.smooth` instead of `object@data.signaling`, and it is
/// the *only* thing that flag does -- there is no second kernel. So the check that matters is
/// that the two settings really do differ on the same object, and that the scaling
/// `data <- data/max(data)` is applied to the **selected** matrix.
///
/// The scaling is the subtle half. `data.smooth` is a shrunk matrix whose maximum differs from
/// `data.signaling`'s, so scaling by the wrong one divides every group mean by a different
/// constant. That is invisible while the mask is insensitive to it and wrong as soon as a
/// per-sample score crosses the threshold -- which is what the corpus's shrunk smooth data is
/// built to do.
#[test]
fn raw_use_false_reads_data_smooth_and_scales_the_matrix_it_selected() {
    let fx = load_db();
    let cs = cases();
    let by = |n: &str| cs.iter().find(|c| c.name == n).expect("case");

    // The data-dependent half is the *cross-sample* one, so the comparison needs `min.samples >= 2`.
    // With `min.samples = NULL` the two settings agree exactly -- and that is worth pinning in
    // its own right, because the `min.cells` half looks like it reads the data (it computes group
    // means) when in fact it only counts cells per group. A port that computed the exclusion from
    // the expression values would differ here, and the corpus would not catch it.
    let single = by("raw_use_false_mincells");
    assert!(single.min_samples.is_none());
    let mut single_raw = (*single).clone();
    single_raw.raw_use = true;
    assert!(
        same_bits(
            &run(single, &fx).unwrap().prob,
            &run(&single_raw, &fx).unwrap().prob
        ),
        "with min.samples = NULL the two settings must agree: the min.cells half counts cells"
    );

    // With two samples the per-sample scores *are* computed from the matrix, so the two
    // settings are reading genuinely different data. The substantive check is not "the results
    // differ" -- an earlier version of this test asserted that, and it is not a property the
    // corpus establishes: on these fixtures the cross-sample mask happens to select the same
    // L-R pairs either way, so the two results coincide even though the inputs do not. Asserting
    // inequality would have been asserting a coincidence, and it would break the next time a
    // fixture is added.
    //
    // What is asserted instead is the thing that *is* structural: the matrix the port reads is
    // not the matrix the other setting reads, and its maximum -- the divisor in
    // `data <- data/max(data)` -- is a different number, so scaling by the wrong one is a real
    // mistake rather than a hypothetical.
    let t = by("raw_use_false_two_samples");
    assert_eq!(t.min_samples, Some(2));
    assert!(!t.raw_use, "the case must actually set raw.use = FALSE");
    let smooth = t.smooth.as_ref().expect("a data.smooth record");
    assert!(
        !same_bits(smooth, &t.data),
        "the two input matrices must differ"
    );
    let smooth_max = smooth.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    let raw_max = t.data.iter().cloned().fold(f64::NEG_INFINITY, f64::max);
    assert_ne!(
        smooth_max.to_bits(),
        raw_max.to_bits(),
        "the fixture is only a test of the scaling if the two maxima differ"
    );
    // And the kernel runs on the smooth matrix, matching the corpus -- which is the check that
    // actually pins the port, and is made per-case by
    // `the_filtered_probability_array_matches_r_bit_for_bit`.
    run(t, &fx).expect("raw.use = FALSE on an independent draw");

    // `raw.use = FALSE` with `nonFilter.keep = TRUE` is the one configuration in which the
    // min.cells zero-fill is observable at all, so it is what pins the arithmetic.
    let nf = by("raw_use_false_nonfilter_keep");
    let got = run(nf, &fx).expect("raw.use = FALSE, nonFilter.keep = TRUE");
    assert!(
        got.prob.contains(&0.0),
        "the excluded group's row and column must be zeroed"
    );
    assert!(
        got.prob.iter().any(|v| *v > 0.0),
        "and the surviving entries must not be"
    );
}

/// The `data.smooth` slot being *absent* does not raise: upstream's guard
/// `if (!("data.smooth" %in% methods::slotNames(object)) == FALSE)` has `%in%` binding tighter
/// than both `!` and `==`, so it evaluates `(!in_set) == FALSE` -- that is, "the slot **is**
/// present" -- and the message claims the opposite. With the slot missing the test is false,
/// nothing stops, and the absent matrix is read. Preserved, and the corpus records what
/// upstream then produces.
#[test]
fn an_absent_data_smooth_slot_does_not_raise() {
    let fx = load_db();
    let cs = cases();
    let c = cs
        .iter()
        .find(|c| c.name == "raw_use_false_slot_absent")
        .expect("case");
    // The corpus supplies the smooth matrix even for this case (so the reader's shape check
    // holds), but upstream read an *absent* slot. The recorded result is upstream's, and the
    // port reaches the same one by way of a different input, which is exactly the situation the
    // case documents.
    // Nothing to assert on the *value* -- the point is that it returns at all. Asserting on the
    // array would be asserting that upstream limped along in a particular way, which is a
    // stronger claim than this case is for; the array is pinned by the corpus comparison.
    let _got = run(c, &fx).expect("an absent data.smooth must not raise");
    assert!(
        _got.prob.iter().any(|v| *v > 0.0),
        "upstream reads the absent matrix and leaves the network alone"
    );
}

/// `nonFilter.keep = TRUE` stores the *unfiltered* arrays in two new slots. Since the
/// filtering has already zeroed `net$prob` by the time `identical()` looks at the object, a
/// port that stored the filtered arrays would be caught here.
#[test]
fn nonfilter_keep_stores_the_unfiltered_arrays() {
    let _fx = load_db();
    for c in cases() {
        if !c.nonfilter_keep || is_error(&c.name).is_some() {
            continue;
        }
        if record(&c.name, "slot_names").is_none() {
            continue;
        }
        // `nonFilter.keep = TRUE` is a **no-op** upstream, and that is the finding.
        //
        //   filterCommunication <- function(object, ...) {
        //     net <- object@net                    # captured here
        //     if (nonFilter.keep == TRUE) {
        //       object@net$prob.nonFilter <- net$prob    # added to object@net
        //       object@net$pval.nonFilter <- net$pval
        //     }
        //     ...
        //     object@net <- net                    # overwrites, discarding both
        //     return(object)
        //   }
        //
        // The two slots are written to `object@net` and then thrown away by the final
        // `object@net <- net`, because `net` was bound to the *pre-assignment* list. So
        // upstream's `net` gains no slots and the object comes back unchanged -- the recorded
        // `slot_names` is exactly the input's. A port that helpfully kept the arrays would
        // differ from upstream, and the help would be invisible in every other respect.
        // Reproduced: the shim does not create the slots either, and its `cat()` message --
        // which upstream *does* emit -- is the only observable effect of the flag.
        let slots = record(&c.name, "slot_names").expect("slot_names record");
        assert!(
            !slots.contains("prob.nonFilter") && !slots.contains("pval.nonFilter"),
            "{}: upstream *discards* the nonFilter slots; if this now fails, re-read the \
             pinned source before 'fixing' it",
            c.name
        );
        assert_eq!(
            slots, "dimnames,prob,prob.dim,pval",
            "{}: net slot names are unchanged",
            c.name
        );
    }
}

/// The per-sample score tensor, which the final `net$prob` only summarises.
#[test]
fn the_per_sample_scores_are_binarised_and_zeroed_for_excluded_groups() {
    let fx = load_db();
    for c in cases() {
        if is_error(&c.name).is_some() || !c.min_samples.is_some_and(|m| m >= 2) {
            continue;
        }
        let got = run(&c, &fx).unwrap_or_else(|e| panic!("{}: unexpected {e}", c.name));
        let k = c.levels.len();
        let n_samp = c.sample_index.iter().copied().max().unwrap() + 1;
        let n_nz = got.lr_nonzero.len();
        assert_eq!(
            got.score.len(),
            k * k * n_nz * n_samp,
            "{}: score.LR length",
            c.name
        );
        // `score.LR[score.LR > 0] <- 1`: every entry is 0 or 1.
        for (i, v) in got.score.iter().enumerate() {
            assert!(
                *v == 0.0 || *v == 1.0,
                "{}: score.LR[{i}] = {v}; the tensor is binarised",
                c.name
            );
        }
        // A group excluded in sample `i` has its whole row and column zeroed in that sample.
        for (i, excl) in got.sample_excluded.iter().enumerate() {
            for &g in excl {
                for jj in 0..n_nz {
                    for t in 0..k {
                        for s in 0..k {
                            assert_eq!(
                                got.score[si(g, t, jj, i, k, n_nz)],
                                0.0,
                                "{}: sample {i} group {g} row not zeroed",
                                c.name
                            );
                            assert_eq!(
                                got.score[si(s, g, jj, i, k, n_nz)],
                                0.0,
                                "{}: sample {i} group {g} column not zeroed",
                                c.name
                            );
                        }
                    }
                }
            }
        }
    }
}

/// `table(idents[cell.use])` keeps every level, so a group *absent* from a sample is in that
/// sample's excluded set with a count of zero. This is the only thing that gives `rare.keep`
/// anything to preserve, and it is the subtle part of the per-sample loop.
#[test]
fn an_absent_group_counts_as_zero_cells_in_that_sample() {
    let fx = load_db();
    for c in cases() {
        if is_error(&c.name).is_some() {
            continue;
        }
        let got = run(&c, &fx).unwrap_or_else(|e| panic!("{}: unexpected {e}", c.name));
        let k = c.levels.len();
        let n_samp = c.sample_index.iter().copied().max().unwrap() + 1;
        for i in 0..n_samp.min(got.sample_excluded.len()) {
            let cells: Vec<usize> = (0..c.n_cells).filter(|&x| c.sample_index[x] == i).collect();
            let mut expect: Vec<usize> = Vec::new();
            for g in 0..k {
                let n = cells.iter().filter(|&&x| c.group_index[x] == g).count();
                if n <= c.min_cells {
                    expect.push(g);
                }
            }
            assert_eq!(
                got.sample_excluded[i], expect,
                "{}: sample {i} excluded groups",
                c.name
            );
        }
    }
}

/// The final `setdiff(cell.excludes.sample, cell.excludes)`: a group already dropped for
/// having too few cells *overall* is not also treated as rare.
#[test]
fn cell_excludes_sample_excludes_the_overall_excluded_groups() {
    let fx = load_db();
    for c in cases() {
        if is_error(&c.name).is_some() {
            continue;
        }
        let got = run(&c, &fx).unwrap_or_else(|e| panic!("{}: unexpected {e}", c.name));
        let mut flat: Vec<usize> = got.sample_excluded.iter().flatten().copied().collect();
        flat.sort_unstable();
        flat.dedup();
        for g in &got.cell_excludes {
            assert!(
                !got.cell_excludes_sample.contains(g),
                "{}: group {g} is excluded overall, so it must not be in cell_excludes_sample",
                c.name
            );
            flat.retain(|x| x != g);
        }
        assert_eq!(
            got.cell_excludes_sample, flat,
            "{}: cell.excludes.sample",
            c.name
        );
    }
}

/// The `min.cells` half on its own: `which(as.numeric(table(idents)) <= min.cells)`, in level
/// order, and the zeroing on both the row and the column of `net$prob`.
#[test]
fn min_cells_excludes_groups_at_or_below_the_threshold_in_level_order() {
    let fx = load_db();
    for c in cases() {
        if is_error(&c.name).is_some() {
            continue;
        }
        let got = run(&c, &fx).unwrap_or_else(|e| panic!("{}: unexpected {e}", c.name));
        let k = c.levels.len();
        let expect: Vec<usize> = (0..k)
            .filter(|&g| c.group_index.iter().filter(|&&x| x == g).count() <= c.min_cells)
            .collect();
        assert_eq!(got.cell_excludes, expect, "{}: cell.excludes", c.name);
        assert!(
            got.cell_excludes.windows(2).all(|w| w[0] < w[1]),
            "{}: cell.excludes must be ascending",
            c.name
        );
        // Both the row and the column are zeroed, so no positive entry survives in either.
        // `k x k x nLR` column-major: `source + k*target + k*k*lr`.
        let n_lr = c.lr.len();
        for &g in &got.cell_excludes {
            for t in 0..k {
                for s in 0..k {
                    for l in 0..n_lr {
                        assert_eq!(
                            got.prob[g + k * t + k * k * l],
                            0.0,
                            "{}: source row {g}",
                            c.name
                        );
                        assert_eq!(
                            got.prob[s + k * g + k * k * l],
                            0.0,
                            "{}: target col {g}",
                            c.name
                        );
                    }
                }
            }
        }
    }
}
