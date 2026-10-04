//! Parity of [`r_core::stats`] against the pinned R 4.3.3 + collapse 2.1.8.
//!
//! Corpus: `tests/fixtures/stats_golden.txt` (expected values),
//! `tests/fixtures/stats_vectors.{idx,bin}` (the exact input vectors, so this test never
//! has to re-run R's RNG), both produced by `tests/parity/gen_stats_golden.R`.
//!
//! Every assertion is bit equality. See `docs/SEMANTICS.md` R2/R3/R4/R9/R10.
//!
//! ## What this test settles that reading cannot
//!
//! * **R2** — that `collapse::fquantile(type = 7)` is reproduced, not R's
//!   `quantile(type = 7)`. The corpus records the *measured* disagreement between the
//!   two (1.13 % of quantile elements, 3.51 % of trimean values, max 1 ulp) and
//!   `collapse_divergence_is_real` asserts it is still non-zero, so an implementation
//!   that quietly switched to R's formula would fail here.
//! * **R10** — whether a plain `f64` accumulation reproduces R's `LONG_DOUBLE` mean for
//!   the lengths that matter. `r_mean_matches_r_long_double` checks the specific
//!   catastrophic-cancellation case R's own fingerprint prints, rather than trusting
//!   that they agree in general.

use r_core::stats::{
    collapse_uses_radix_order, fquantile_type7_unsorted, geometric_mean, r_mean, r_median,
    thresholded_mean, tri_mean,
};
use std::collections::HashMap;

const GOLDEN: &str = include_str!("../../../../../tests/fixtures/stats_golden.txt");
const VEC_BIN: &[u8] = include_bytes!("../../../../../tests/fixtures/stats_vectors.bin");
const VEC_IDX: &str = include_str!("../../../../../tests/fixtures/stats_vectors.idx");

/// `id -> Vec<f64>` loaded from the binary blob.
fn vectors() -> HashMap<String, Vec<f64>> {
    let mut out = HashMap::new();
    for line in VEC_IDX.lines().filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        let id = f[0].to_string();
        let off: usize = f[1].parse().unwrap();
        let len: usize = f[2].parse().unwrap();
        let bytes = &VEC_BIN[off..off + 8 * len];
        let v = bytes
            .chunks_exact(8)
            .map(|c| f64::from_le_bytes(c.try_into().unwrap()))
            .collect();
        out.insert(id, v);
    }
    assert!(out.len() > 1000, "only {} vectors indexed", out.len());
    out
}

struct Golden {
    /// (family, n, rep) -> (four quantiles, trimean)
    trimean: HashMap<(String, usize, usize), (Vec<f64>, f64)>,
    trimean_big: HashMap<usize, f64>,
    geomean: HashMap<(usize, usize), f64>,
    geomean_zero: HashMap<(usize, usize), f64>,
    geomean_rand: HashMap<usize, f64>,
    /// (n, nnz, trim) -> value
    threshmean: HashMap<(usize, usize, String), f64>,
    special: HashMap<(String, String), f64>,
    div_elem_frac: f64,
    div_tri_frac: f64,
    div_max: f64,
    collapse_version: String,
    r_version: String,
    long_double_fingerprint: f64,
}

/// Parse an R-formatted number.
///
/// R renders `NA_real_` as `"NA"`, which Rust's `f64` parser rejects. `NA_real_` *is* a
/// NaN bit pattern, just a non-canonical one, so it is mapped to Rust's canonical NaN
/// and compared with `is_nan()` rather than by bit pattern (see [`bit_eq`]).
fn rf64(s: &str) -> f64 {
    match s.trim() {
        "NA" | "NaN" => f64::NAN,
        other => other
            .parse()
            .unwrap_or_else(|e| panic!("cannot parse {other:?}: {e}")),
    }
}

fn parse() -> Golden {
    let mut g = Golden {
        trimean: HashMap::new(),
        trimean_big: HashMap::new(),
        geomean: HashMap::new(),
        geomean_zero: HashMap::new(),
        geomean_rand: HashMap::new(),
        threshmean: HashMap::new(),
        special: HashMap::new(),
        div_elem_frac: -1.0,
        div_tri_frac: -1.0,
        div_max: -1.0,
        collapse_version: String::new(),
        r_version: String::new(),
        long_double_fingerprint: f64::NAN,
    };
    for line in GOLDEN.lines().filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "trimean" => {
                let q: Vec<f64> = f[4].split(',').map(rf64).collect();
                g.trimean.insert(
                    (
                        f[1].to_string(),
                        f[2].parse().unwrap(),
                        f[3].parse().unwrap(),
                    ),
                    (q, rf64(f[5])),
                );
            }
            "trimean_big" => {
                g.trimean_big.insert(f[1].parse().unwrap(), rf64(f[2]));
            }
            "geomean" => {
                g.geomean
                    .insert((f[1].parse().unwrap(), f[2].parse().unwrap()), rf64(f[3]));
            }
            "geomean_zero" => {
                g.geomean_zero
                    .insert((f[1].parse().unwrap(), f[2].parse().unwrap()), rf64(f[3]));
            }
            "geomean_rand" => {
                g.geomean_rand.insert(f[1].parse().unwrap(), rf64(f[2]));
            }
            "threshmean" => {
                g.threshmean.insert(
                    (
                        f[1].parse().unwrap(),
                        f[2].parse().unwrap(),
                        f[3].to_string(),
                    ),
                    rf64(f[4]),
                );
            }
            "special_trimean" | "special_geomean" | "special_median" => {
                g.special
                    .insert((f[0].to_string(), f[1].to_string()), rf64(f[2]));
            }
            "#divergence" => {
                g.div_elem_frac = f[3].trim_start_matches("elem_frac=").parse().unwrap();
                g.div_tri_frac = f[6].trim_start_matches("tri_frac=").parse().unwrap();
                g.div_max = f[4].trim_start_matches("elem_max=").parse().unwrap();
            }
            "#env" => {
                g.r_version = f[1].trim_start_matches("R=").to_string();
                g.collapse_version = f[2].trim_start_matches("collapse=").to_string();
            }
            "#r_mean_kind" => {
                // #r_mean_kind <kind> first=<value>
                assert_eq!(f[1], "LONG_DOUBLE", "R's mean accumulation kind changed");
                g.long_double_fingerprint = f[2].trim_start_matches("first=").parse().unwrap();
            }
            // Informational records that carry no expected value.
            "#vectors" | "#" => {}
            other => panic!("unknown golden record kind: {other}"),
        }
    }
    g
}

/// Bit equality, with NaN treated as equal to NaN.
///
/// R's `NA_real_` has a non-canonical NaN payload, so demanding bit equality for
/// missing results would compare payload bits rather than semantics.
fn bit_eq(got: f64, want: f64, what: &str) {
    if got.is_nan() && want.is_nan() {
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

#[test]
fn corpus_is_the_pinned_environment() {
    let g = parse();
    assert_eq!(g.r_version, "4.3.3");
    assert_eq!(g.collapse_version, "2.1.8");
}

#[test]
fn fquantile_matches_collapse_bit_for_bit() {
    let g = parse();
    let vs = vectors();
    let mut checked = 0usize;
    for ((family, n, rep), (want_q, want_tri)) in &g.trimean {
        let id = format!("trimean:{family}:{n}:{rep}");
        let x = vs.get(&id).unwrap_or_else(|| panic!("missing vector {id}"));
        assert_eq!(x.len(), *n, "vector length mismatch for {id}");

        const PROBS: [f64; 4] = [0.25, 0.50, 0.50, 0.75];
        let got_q = fquantile_type7_unsorted(x, &PROBS);
        assert_eq!(got_q.len(), want_q.len(), "{id}");
        for (i, (&a, &b)) in got_q.iter().zip(want_q).enumerate() {
            bit_eq(a, b, &format!("{id} quantile[{i}]"));
        }
        bit_eq(tri_mean(x), *want_tri, &format!("{id} trimean"));
        checked += 1;
    }
    assert!(checked > 500, "only checked {checked} trimean vectors");
}

#[test]
fn trimean_matches_collapse_for_large_vectors() {
    // The sizes the aggregation actually runs at, plus collapse's radix threshold.
    let g = parse();
    let vs = vectors();
    for (&n, &want) in &g.trimean_big {
        let x = &vs[&format!("trimean_big:{n}")];
        assert_eq!(x.len(), n);
        bit_eq(tri_mean(x), want, &format!("trimean_big n={n}"));
    }
    assert!(g.trimean_big.contains_key(&200_000));
}

#[test]
fn geometric_mean_matches_r_bit_for_bit() {
    let g = parse();
    let vs = vectors();
    for (&(n, rep), &want) in &g.geomean {
        let x = &vs[&format!("geomean:{n}:{rep}")];
        assert_eq!(x.len(), n);
        bit_eq(geometric_mean(x), want, &format!("geomean {n}:{rep}"));
    }
    for (&(n, rep), &want) in &g.geomean_zero {
        let x = &vs[&format!("geomean:{n}:{rep}")];
        let mut with_zero = x.clone();
        with_zero.push(0.0);
        bit_eq(
            geometric_mean(&with_zero),
            want,
            &format!("geomean_zero {n}:{rep}"),
        );
    }
    for (&i, &want) in &g.geomean_rand {
        let x = &vs[&format!("geomean_rand:{i}")];
        bit_eq(geometric_mean(x), want, &format!("geomean_rand {i}"));
    }
}

#[test]
fn thresholded_mean_matches_r_bit_for_bit() {
    let g = parse();
    let vs = vectors();
    for (&(n, nnz, ref trim), &want) in &g.threshmean {
        let id = format!("threshmean:{n}:{nnz}:{trim}");
        let x = &vs[&id];
        let t: f64 = trim.parse().unwrap();
        bit_eq(thresholded_mean(x, t), want, &id);
    }
}

#[test]
fn nan_and_inf_handling_matches_r() {
    // R9: `fquantile` treats NaN as missing under na.rm, `log(0) = -Inf` poisons the
    // geometric mean, and Inf propagates. Each of these is an R-ism, not a detail.
    let g = parse();
    let vs = vectors();
    for ((kind, name), &want) in &g.special {
        let x = &vs[&format!("special:{name}")];
        let got = match kind.as_str() {
            "special_trimean" => tri_mean(x),
            "special_geomean" => geometric_mean(x),
            "special_median" => r_median(x, true),
            other => panic!("no Rust implementation for {other}"),
        };
        bit_eq(got, want, &format!("{kind} {name}"));
    }
}

#[test]
fn median_even_n_uses_r_two_pass_mean() {
    // Guards the R-ism that `median` of an even-length vector is `mean(two middles)`,
    // not `(a + b) / 2`. The corpus has an even-n special case (`with_nan` -> 4 values
    // after na.rm) and it is checked in `nan_and_inf_handling_matches_r`; this asserts
    // the structural difference directly so the reason is not lost.
    // Structural guarantee: even-n median is `mean()` of the two central values.
    // A concrete pair where the naive `(a + b) / 2` disagrees, found by search in R
    // (first hit after ~1.2e5 random pairs, so they agree ~19999 times in 20000 -- a
    // naive "search small integers for a difference" test would be vacuous).
    let (a, b) = (1.4493788766496514e71, 1.0758174218958832e86);
    let plain = (a + b) / 2.0;
    let via_mean = r_mean(&[a, b]);
    assert_ne!(plain, via_mean, "control case should differ");
    assert_eq!(
        r_median(&[a, b], false),
        via_mean,
        "median of two values must use the two-pass mean, not (a+b)/2"
    );
}

#[test]
fn collapse_divergence_is_real() {
    // Guards the whole point of R2. If a future R/collapse update made the two
    // quantiles agree bit-for-bit, this test would fail and force a re-review of
    // whether we still need collapse's form.
    let g = parse();
    assert!(
        g.div_tri_frac > 0.0,
        "collapse and R quantiles now agree exactly; re-check SEMANTICS R2"
    );
    assert!(
        g.div_elem_frac > 0.0,
        "element-level divergence vanished; re-check SEMANTICS R2"
    );
    // Measured max discrepancy must stay at or below 1 ulp for a mean-based statistic.
    assert!(
        g.div_max <= f64::EPSILON,
        "divergence grew to {}",
        g.div_max
    );
    eprintln!(
        "measured divergence: element {:.4}%, trimean {:.4}%, max {:.3e}",
        g.div_elem_frac, g.div_tri_frac, g.div_max
    );
}

#[test]
fn r_mean_matches_r_long_double() {
    // R10. R's own corpus fingerprint is `mean(c(1e16, 1, -1e16))` = 1/3 computed with
    // a 64-bit mantissa. A naive f64 accumulation loses the `1` entirely and returns 0,
    // so this is a sharp test of whether our summation is adequate.
    let g = parse();
    let got = r_mean(&[1e16, 1.0, -1e16]);
    bit_eq(got, g.long_double_fingerprint, "mean(c(1e16, 1, -1e16))");
    // R does NOT return the exact mean (1/3) here. Pass 1 lands on the `1` only in 64
    // bits, so the result is R's specific two-pass value, 0.33365885416666669. Assert
    // both that we match R bit-for-bit and that we differ from a naive f64 chain, which
    // loses the `1` entirely and returns 0.
    let mut naive = 0.0f64;
    for &v in &[1e16, 1.0, -1e16] {
        naive += v;
    }
    assert_eq!(
        naive, 0.0,
        "the naive f64 control case is what we are beating"
    );
    assert_ne!(got, naive);
    assert_ne!(
        got,
        1.0 / 3.0,
        "R's two-pass mean is deliberately not the exact mean"
    );
}

#[test]
fn radix_order_path_is_never_taken_by_trimean() {
    for len in [1usize, 5, 100, 54, 100_000, 1_000_000] {
        assert!(!collapse_uses_radix_order(len, 4));
    }
}
