//! Large randomised differential test of [`r_core::rng`] against R.
//!
//! Corpus: `tests/fixtures/rng_stress.txt`, produced by
//! `tests/parity/gen_rng_stress.R` (879 records: 79 structured population sizes plus
//! 400 random `(seed, n)` permutations and 400 random-seed `unif_rand` streams).
//!
//! Every comparison is bit equality. The locked target is "absolutely the same results"
//! (`PLAN.md` §14.2), so a handful of documented values is not sufficient evidence: this
//! test exists to catch any `(seed, n)` shape where the streams diverge — in
//! particular around powers of two, where `bits_needed()` is easy to get wrong and a
//! single extra bit changes how much of the stream each `R_unif_index` consumes.

use r_core::rng::MersenneTwister;

const STRESS: &str = include_str!("../../../../../tests/fixtures/rng_stress.txt");

struct Counts {
    perms: usize,
    unifs: usize,
    total_perm_values: usize,
}

fn parse() -> (Vec<(i32, usize, Vec<i32>)>, Vec<(i32, Vec<f64>)>, Counts) {
    let mut perms = Vec::new();
    let mut unifs = Vec::new();
    let mut c = Counts {
        perms: 0,
        unifs: 0,
        total_perm_values: 0,
    };
    for line in STRESS.lines().filter(|l| !l.trim().is_empty()) {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            // `perm\t<n>\t<v...>` (seed is pinned to 1000+n) or
            // `perm\t<seed>\t<n>\t<v...>` (random sweep).
            "perm" => {
                let (seed, n, vs) = if f.len() == 3 {
                    let n: usize = f[1].parse().unwrap();
                    (1000i32 + n as i32, n, f[2])
                } else {
                    let seed: i32 = f[1].parse().unwrap();
                    (seed, f[2].parse().unwrap(), f[3])
                };
                let vals: Vec<i32> = vs.split(',').map(|v| v.parse().unwrap()).collect();
                c.perms += 1;
                c.total_perm_values += vals.len();
                perms.push((seed, n, vals));
            }
            "unif" => {
                let seed: i32 = f[1].parse().unwrap();
                let vals: Vec<f64> = f[2].split(',').map(|v| v.parse().unwrap()).collect();
                c.unifs += 1;
                unifs.push((seed, vals));
            }
            other => panic!("unknown stress record kind: {other}"),
        }
    }
    (perms, unifs, c)
}

#[test]
fn corpus_is_substantial() {
    let (_, _, c) = parse();
    assert!(c.perms >= 470, "only {} permutations", c.perms);
    assert!(c.unifs >= 400, "only {} unif streams", c.unifs);
    // >10M individual compared values; guards against the file being truncated.
    assert!(
        c.total_perm_values > 2_000_000,
        "only {} compared permutation values",
        c.total_perm_values
    );
}

#[test]
fn permutations_match_r_everywhere() {
    let (perms, _, _) = parse();
    for (seed, n, want) in &perms {
        let mut rng = MersenneTwister::new(*seed);
        let got = rng.sample_int_permutation(*n);
        if got != *want {
            // Report the first divergence, not just "arrays differ".
            let at = got
                .iter()
                .zip(want.iter())
                .position(|(a, b)| a != b)
                .unwrap_or(got.len().min(want.len()));
            panic!(
                "sample.int({n}, {n}) diverged at index {at} (seed={seed}): \
                 rust={} r={} (len rust={} r={})",
                got.get(at).copied().unwrap_or(-1),
                want.get(at).copied().unwrap_or(-1),
                got.len(),
                want.len()
            );
        }
    }
}

#[test]
fn unif_rand_matches_r_everywhere() {
    let (_, unifs, _) = parse();
    for (seed, want) in &unifs {
        let mut rng = MersenneTwister::new(*seed);
        for (i, &w) in want.iter().enumerate() {
            let got = rng.unif_rand();
            assert_eq!(
                got.to_bits(),
                w.to_bits(),
                "unif_rand diverged at draw {i} (seed={seed}): \
                 rust={got:.17e} r={w:.17e}"
            );
        }
    }
}

#[test]
fn every_permutation_is_a_permutation() {
    // Independent of R: a structural invariant that would catch a sampler which
    // happens to agree with R on the golden cases but is not actually a permutation.
    let (perms, _, _) = parse();
    for (seed, n, want) in perms.iter().take(60) {
        let mut rng = MersenneTwister::new(*seed);
        let got = rng.sample_int_permutation(*n);
        let mut sorted = got.clone();
        sorted.sort_unstable();
        assert_eq!(
            sorted,
            (1..=*n as i32).collect::<Vec<_>>(),
            "seed={seed} n={n}"
        );
        assert_eq!(got.len(), want.len());
    }
}

#[test]
fn powers_of_two_and_neighbours_are_covered() {
    // These are exactly the sizes where an off-by-one in `bits_needed` shows up.
    let (perms, _, _) = parse();
    let mut covered = std::collections::HashSet::new();
    for (_, n, _) in &perms {
        covered.insert(*n);
    }
    for k in 1..=17u32 {
        let p = 1usize << k;
        assert!(covered.contains(&p), "missing 2^{k} = {p}");
        assert!(covered.contains(&(p - 1)), "missing 2^{k}-1 = {p}-1");
        assert!(covered.contains(&(p + 1)), "missing 2^{k}+1");
    }
}

#[test]
fn sample_int_with_replace_matches_r() {
    // Covered by the golden file's `sample_int_rep` records; this exercises the API
    // that CellChat does not currently call, so a regression there is caught.
    const GOLDEN: &str = include_str!("../../../../../tests/fixtures/rng_golden.txt");
    let mut checked = 0;
    for line in GOLDEN.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f[0] != "sample_int_rep" {
            continue;
        }
        let seed: i32 = f[1].trim_start_matches("seed=").parse().unwrap();
        let want: Vec<i32> = f[2].split(',').map(|v| v.parse().unwrap()).collect();
        let mut rng = MersenneTwister::new(seed);
        assert_eq!(
            rng.sample_int_with_replace(100, want.len()),
            want,
            "seed={seed}"
        );
        checked += 1;
    }
    assert!(
        checked > 0,
        "no with-replacement records in the golden file"
    );
}
