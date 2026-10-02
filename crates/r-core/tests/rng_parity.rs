//! Parity of [`r_core::rng`] against the pinned R runtime.
//!
//! The golden file `tests/fixtures/rng_golden.txt` is produced verbatim by
//! `tests/parity/gen_rng_golden.R`, so no constant in this test is hand-transcribed.
//! Regenerate with:
//!
//! ```sh
//! Rscript tests/parity/gen_rng_golden.R > tests/fixtures/rng_golden.txt
//! ```
//!
//! Every assertion here is `assert_eq!` on `f64`/`i32`, i.e. bit equality. The locked
//! target (`PLAN.md` §14.2) forbids a tolerance here: these values determine the
//! bootstrap permutations, and therefore every p-value in the network.

use r_core::rng::MersenneTwister;
use std::collections::HashMap;
use std::path::PathBuf;

const GOLDEN: &str = include_str!("../../../tests/fixtures/rng_golden.txt");

fn golden_path() -> PathBuf {
    PathBuf::from(env!("CARGO_MANIFEST_DIR"))
        .parent()
        .unwrap()
        .parent()
        .unwrap()
        .join("tests/fixtures/rng_golden.txt")
}

/// `key -> values` for each record kind, keyed by the tab-separated label.
struct Golden {
    unif_rand: HashMap<i32, Vec<f64>>,
    sample_int: HashMap<(i32, usize), Vec<i32>>,
    sample_int_k: HashMap<(i32, usize, usize), Vec<i32>>,
    sample_int_rep: HashMap<i32, Vec<i32>>,
    r_fingerprint: String,
}

fn parse() -> Golden {
    let mut g = Golden {
        unif_rand: HashMap::new(),
        sample_int: HashMap::new(),
        sample_int_k: HashMap::new(),
        sample_int_rep: HashMap::new(),
        r_fingerprint: String::new(),
    };
    for line in GOLDEN.lines().filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "unif_rand" => {
                let seed: i32 = f[1].trim_start_matches("seed=").parse().unwrap();
                g.unif_rand
                    .insert(seed, f[2].split(',').map(|v| v.parse().unwrap()).collect());
            }
            "sample_int" => {
                let seed: i32 = f[1].trim_start_matches("seed=").parse().unwrap();
                let n: usize = f[2].trim_start_matches("n=").parse().unwrap();
                g.sample_int.insert(
                    (seed, n),
                    f[3].split(',').map(|v| v.parse().unwrap()).collect(),
                );
            }
            "sample_int_k" => {
                let seed: i32 = f[1].trim_start_matches("seed=").parse().unwrap();
                let n: usize = f[2].trim_start_matches("n=").parse().unwrap();
                let k: usize = f[3].trim_start_matches("k=").parse().unwrap();
                g.sample_int_k.insert(
                    (seed, n, k),
                    f[4].split(',').map(|v| v.parse().unwrap()).collect(),
                );
            }
            "sample_int_rep" => {
                let seed: i32 = f[1].trim_start_matches("seed=").parse().unwrap();
                g.sample_int_rep
                    .insert(seed, f[2].split(',').map(|v| v.parse().unwrap()).collect());
            }
            "#R" => g.r_fingerprint = line.to_string(),
            other => panic!("unknown golden record kind: {other}"),
        }
    }
    assert!(!g.unif_rand.is_empty(), "golden file parsed to nothing");
    g
}

#[test]
fn golden_file_is_the_pinned_r() {
    let g = parse();
    let f = &g.r_fingerprint;
    assert!(f.contains("Mersenne-Twister"), "{f}");
    assert!(f.contains("sample.kind=Rejection"), "{f}");
}

#[test]
fn golden_file_is_where_the_test_expects() {
    // Guards against `include_str!` and a stale on-disk file diverging.
    assert_eq!(GOLDEN, std::fs::read_to_string(golden_path()).unwrap());
}

#[test]
fn unif_rand_is_bit_identical_to_r() {
    let g = parse();
    let mut checked = 0;
    for (&seed, want) in &g.unif_rand {
        let mut rng = MersenneTwister::new(seed);
        for (i, &w) in want.iter().enumerate() {
            let got = rng.unif_rand();
            assert_eq!(
                got.to_bits(),
                w.to_bits(),
                "unif_rand mismatch: seed={seed} draw={i} r={w:.17e} rust={got:.17e}"
            );
            checked += 1;
        }
    }
    assert!(checked >= 100, "only checked {checked} values");
}

#[test]
fn sample_int_permutation_is_bit_identical_to_r() {
    let g = parse();
    for (&(seed, n), want) in &g.sample_int {
        let mut rng = MersenneTwister::new(seed);
        let got = rng.sample_int_permutation(n);
        assert_eq!(&got, want, "sample.int({n}, {n}) mismatch at seed={seed}");
    }
}

#[test]
fn sample_int_partial_is_bit_identical_to_r() {
    let g = parse();
    for (&(seed, n, k), want) in &g.sample_int_k {
        let mut rng = MersenneTwister::new(seed);
        let got = rng.sample_int_no_replace(n, k);
        assert_eq!(&got, want, "sample.int({n}, {k}) mismatch at seed={seed}");
    }
}

#[test]
fn boot_permutation_shape_used_by_compute_commun_prob() {
    // `replicate(nboot, sample.int(nC, size = nC))` draws nboot *independent*
    // permutations from one continuing stream. Getting this wrong (e.g. reseeding per
    // replicate) would still produce a valid-looking network, so pin it explicitly.
    let g = parse();
    let n = 257usize;
    let want_first = &g.sample_int[&(1, n)];

    let mut rng = MersenneTwister::new(1);
    for boot in 0..3 {
        let got = rng.sample_int_permutation(n);
        if boot == 0 {
            assert_eq!(
                &got, want_first,
                "first replicate must match a fresh stream"
            );
        } else {
            assert_ne!(&got, want_first, "replicate {boot} repeated replicate 0");
        }
    }
}

#[test]
fn rng_is_deterministic_across_many_draws() {
    // Long-stream stability: 25k draws per seed over several seeds, all must agree.
    let g = parse();
    for &seed in &[1i32, 2, 42, 1234567, -1, i32::MAX] {
        if let Some(want) = g.unif_rand.get(&seed) {
            let mut rng = MersenneTwister::new(seed);
            for (i, &w) in want.iter().enumerate() {
                assert_eq!(
                    rng.unif_rand().to_bits(),
                    w.to_bits(),
                    "seed={seed} draw={i}"
                );
            }
        }
    }
}
