//! Differential tests for the deterministic half of `computeCentralityLocal`, against
//! installed igraph 2.3.4.
//!
//! The oracle is `tests/fixtures/centrality_golden.txt`, written by
//! `tests/parity/gen_centrality_golden.R` from igraph itself: 60 randomized graphs engineered for
//! ties, zeros, loops and cancellation-sensitive magnitudes, 7 reciprocal/edge cases, 2 error
//! cases, and 6 real pathway slices from the tutorial object. Every value is `%.17g`, so the
//! comparisons below are on bit patterns, not decimals.
//!
//! What is pinned, and why each piece is shaped the way it is:
//!
//! * Degrees and strengths: `%.17g` text of igraph's output against `%.17g` text of the Rust
//!   values, compared as *strings*. String comparison is deliberate -- it distinguishes `-0.0`
//!   from `0.0`, which a parsed-float comparison would conflate, and `%.17g` uniquely encodes a
//!   double so equal strings mean equal bits.
//! * Betweenness: same treatment. The values are Brandes accumulations over epsilon-compared
//!   Dijkstra distances; any difference in heap order, tie rule, or accumulation order shows up
//!   as a different last bit, usually in the first few cases.
//! * Error cases: the port must raise igraph's exact message, `Source:` suffix included, because
//!   that suffix is part of `conditionMessage()` and the gates compare it.
//!
//! Deliberately absent: any assertion on `hub_score` / `authority_score` / `eigen_centrality` /
//! `page_rank`. Those are iterative solvers whose output differs run to run on identical input
//! (measured), so a golden pins nothing and a passing test would prove nothing.

use r_core::centrality::{betweenness, degrees_strength};
use std::collections::HashMap;

const GOLDEN: &str = include_str!("../../../../../tests/fixtures/centrality_golden.txt");

/// Compare golden text against computed values by *parsed bits*, not by string equality.
///
/// R's `%.17g` prints the shortest round-trip spelling (`4`, `17`, `1e+16`), while Rust's `{:.17e}`
/// always uses scientific notation (`4.00000000000000000e0`) -- so equal values have different
/// spellings and string comparison fails on every integral result. Parsing both sides and comparing
/// bit patterns keeps the distinction that matters (`-0.0` versus `0.0`) while ignoring the one
/// that does not (spelling). `%.17g` round-trips every double uniquely, so parsing is exact.
///
/// The one deliberate blind spot: `NaN` payloads. R prints every NaN as `NaN`, so a payload the
/// port produced differently is invisible here -- but NaN outputs do not occur in these measures
/// (degrees are counts, strengths of finite inputs are finite-or-`Inf`, betweenness likewise),
/// and any `NaN` on either side fails loudly rather than passing quietly.
struct Case {
    name: String,
    k: usize,
    matrix: Vec<f64>,
    /// Construction failed: both entry points must fail with this text.
    construct_error: Option<String>,
    /// Only betweenness failed (reciprocal checks): degrees/strength records are present and
    /// must match; betweenness must fail with this text.
    betweenness_error: Option<String>,
    records: HashMap<String, String>,
}

fn parse() -> Vec<Case> {
    let mut cases = Vec::new();
    let mut lines = GOLDEN.lines().peekable();
    while let Some(line) = lines.next() {
        if line.trim().is_empty() {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        assert_eq!(f[0], "case", "expected a case record, found {line:?}");
        let name = f[1].to_string();
        let k: usize = f[2].parse().expect("k");
        let marker = lines.next().expect("matrix marker");
        assert_eq!(marker, "matrix", "{name}: expected the matrix marker");
        let mut matrix = Vec::with_capacity(k * k);
        for _ in 0..k {
            let row = lines.next().expect("matrix row");
            for v in row.split_whitespace() {
                // `%.17g` never emits these, but the generator writes R's own spellings for the
                // non-finite entries (`NA`, `NaN`, `Inf`), and Rust's float parser accepts none
                // of the missing-value spellings -- so they are mapped by hand, with `NA_real_`
                // getting its exact bit pattern rather than a canonical NaN.
                let x = match v {
                    "NA" => f64::from_bits(0x7ff0_0000_0000_07a2),
                    "NaN" => f64::NAN,
                    "Inf" | "+Inf" => f64::INFINITY,
                    "-Inf" => f64::NEG_INFINITY,
                    other => other.parse::<f64>().expect("a double"),
                };
                matrix.push(x);
            }
        }
        assert_eq!(matrix.len(), k * k, "{name}: matrix shape");
        // Records run until the next `case` line or EOF. Tags: `nedges\tN` (informational,
        // skipped), bare measure tags with a one-line payload, `error\t<line1>` with optional
        // `error2\t<rest>` (an error message spans lines verbatim, so it cannot fit the
        // tag-plus-one-line shape), and the bare `betweenness_error` marker after which the
        // `error` lines belong to betweenness rather than construction.
        let mut construct_error: Option<String> = None;
        let mut betweenness_error: Option<String> = None;
        let mut pending_betweenness_error = false;
        let mut error_buf: Option<String> = None;
        let mut records = HashMap::new();
        // Local flush: an `error` line without a following `error2` is complete, but that is
        // only known when the next record starts -- so completion is decided inline below.
        while let Some(&next) = lines.peek() {
            if next.trim().is_empty() {
                lines.next();
                continue;
            }
            if next.starts_with("case\t") {
                break;
            }
            if next == "betweenness_error" {
                lines.next();
                pending_betweenness_error = true;
                continue;
            }
            if let Some(rest) = next.strip_prefix("error\t") {
                let rest = rest.to_string();
                lines.next();
                match error_buf.take() {
                    // A second `error` line with no intervening `error2` means the previous one
                    // was complete on its own (single-line message).
                    Some(prev) => {
                        if pending_betweenness_error {
                            betweenness_error = Some(prev);
                        } else {
                            construct_error = Some(prev);
                        }
                        pending_betweenness_error = false;
                        error_buf = Some(rest);
                    }
                    None => error_buf = Some(rest),
                }
                continue;
            }
            if let Some(rest) = next.strip_prefix("error2\t") {
                lines.next();
                let first = error_buf
                    .take()
                    .expect("error2 without a preceding error line");
                let full = format!("{first}\n{rest}");
                if pending_betweenness_error {
                    betweenness_error = Some(full);
                } else {
                    construct_error = Some(full);
                }
                pending_betweenness_error = false;
                continue;
            }
            if next.starts_with("nedges\t") {
                lines.next();
                continue;
            }
            if next == "warning" {
                lines.next();
                records.insert("warning".to_string(), "present".to_string());
                continue;
            }
            // A bare measure tag; payload is the single following line. A pending lone `error`
            // line is complete at this point (no `error2` followed it).
            if let Some(buf) = error_buf.take() {
                if pending_betweenness_error {
                    betweenness_error = Some(buf);
                } else {
                    construct_error = Some(buf);
                }
                pending_betweenness_error = false;
            }
            let tag = lines.next().unwrap().to_string();
            let payload = lines.next().unwrap_or("").to_string();
            records.insert(tag, payload);
        }
        // End of case (or file) completes a trailing lone `error` line.
        if let Some(buf) = error_buf.take() {
            if pending_betweenness_error {
                betweenness_error = Some(buf);
            } else {
                construct_error = Some(buf);
            }
        }
        cases.push(Case {
            name,
            k,
            matrix,
            construct_error,
            betweenness_error,
            records,
        });
    }
    cases
}

fn parse_token(tok: &str) -> u64 {
    match tok {
        "NA" => 0x7ff0_0000_0000_07a2,
        "NaN" => f64::NAN.to_bits(),
        "Inf" | "+Inf" => f64::INFINITY.to_bits(),
        "-Inf" => f64::NEG_INFINITY.to_bits(),
        other => other.parse::<f64>().expect("a double").to_bits(),
    }
}

/// Asserts `got` parses to the golden text's bits, elementwise. A golden `NaN` matches any NaN
/// payload (R's spelling carries none); everything else must agree bit for bit.
fn assert_bits_eq(case: &str, measure: &str, got: &[f64], want_text: &str) {
    let want: Vec<&str> = want_text.split_whitespace().collect();
    assert_eq!(
        got.len(),
        want.len(),
        "{case}: {measure} length {} vs {}",
        got.len(),
        want.len()
    );
    for (i, (g, w)) in got.iter().zip(want.iter()).enumerate() {
        let gb = g.to_bits();
        let wb = parse_token(w);
        let both_nan = g.is_nan() && f64::from_bits(wb).is_nan();
        assert!(
            both_nan || gb == wb,
            "{case}: {measure}[{i}] differs\n  got  {g:.17e} ({gb:#x})\n  want {w} ({wb:#x})"
        );
    }
}

/// The whole corpus, one assertion per measure per case. Split into two tests so a strength
/// failure and a betweenness failure report independently rather than masking each other.
#[test]
fn degrees_and_strength_match_igraph_bit_for_bit() {
    let mut checked = 0;
    for c in parse() {
        if c.construct_error.is_some() {
            continue;
        }
        let d = degrees_strength(&c.matrix, c.k).expect("constructed test matrix");
        for (measure, got) in [
            ("outdeg_unweighted", &d.outdeg_unweighted),
            ("indeg_unweighted", &d.indeg_unweighted),
            ("outdeg", &d.outdeg),
            ("indeg", &d.indeg),
        ] {
            let want = c.records.get(measure).unwrap_or_else(|| {
                panic!("{}: no {} record", c.name, measure);
            });
            assert_bits_eq(&c.name, measure, got, want);
        }
        checked += 1;
    }
    assert!(
        checked >= 60,
        "only {checked} cases carried measure records"
    );
}

#[test]
fn betweenness_matches_igraph_bit_for_bit() {
    let mut checked = 0;
    for c in parse() {
        if c.construct_error.is_some() || c.betweenness_error.is_some() {
            continue;
        }
        let (got, tiny) = betweenness(&c.matrix, c.k).expect("successful test matrix");
        // The warning flag tracks igraph's warning exactly: it fires iff a reciprocal weight is
        // at or below 1e-10, and the corpus records its presence per case. A flag that disagreed
        // would mean the R shim warns where upstream stays silent or vice versa.
        let want_warn = c.records.contains_key("warning");
        assert_eq!(
            tiny,
            want_warn,
            "{}: tiny-weights flag {tiny} but golden {} a warning",
            c.name,
            if want_warn { "records" } else { "records no" }
        );
        let want = c
            .records
            .get("betweenness")
            .unwrap_or_else(|| panic!("{}: no betweenness record", c.name));
        assert_bits_eq(&c.name, "betweenness", &got, want);
        checked += 1;
    }
    assert!(
        checked >= 60,
        "only {checked} cases carried betweenness records"
    );
}

#[test]
fn construction_errors_raise_igraphs_exact_text() {
    // `graph_from_adjacency_matrix` raises before any edge exists, so both entry points must fail
    // with the construction text -- and upstream's `computeCentralityLocal` never reaches degrees,
    // strength, or betweenness on such input either.
    let mut checked = 0;
    for c in parse() {
        let Some(msg) = &c.construct_error else {
            continue;
        };
        let e1 = degrees_strength(&c.matrix, c.k).unwrap_err();
        let e2 = betweenness(&c.matrix, c.k).unwrap_err();
        assert_eq!(e1.to_string(), *msg, "{}: degrees error text", c.name);
        assert_eq!(e2.to_string(), *msg, "{}: betweenness error text", c.name);
        assert_eq!(e1, e2, "{}: error variants agree", c.name);
        checked += 1;
    }
    assert!(checked >= 2, "only {checked} construction-error cases");
}

#[test]
fn reciprocal_errors_leave_degrees_and_strength_intact() {
    // Negative (or NaN-reciprocal) weights die in the reciprocal checks, *after* construction:
    // degrees and strength on the same matrices must still match, and betweenness must fail with
    // the reciprocal text. A single early-return error record could not express this split, which
    // is why the golden records measures and betweenness outcomes separately.
    let mut checked = 0;
    for c in parse() {
        let Some(msg) = &c.betweenness_error else {
            continue;
        };
        let d = degrees_strength(&c.matrix, c.k)
            .unwrap_or_else(|e| panic!("{}: degrees must succeed: {e}", c.name));
        for (measure, got) in [("outdeg", &d.outdeg), ("indeg", &d.indeg)] {
            let want = c.records.get(measure).unwrap_or_else(|| {
                panic!("{}: no {} record", c.name, measure);
            });
            assert_bits_eq(&c.name, measure, got, want);
        }
        let e = betweenness(&c.matrix, c.k).unwrap_err();
        assert_eq!(e.to_string(), *msg, "{}: betweenness error text", c.name);
        checked += 1;
    }
    assert!(checked >= 2, "only {checked} reciprocal-error cases");
}
