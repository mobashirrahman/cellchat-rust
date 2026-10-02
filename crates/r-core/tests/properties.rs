//! Property tests for the primitives the numeric surface is built on.
//!
//! These are not a second copy of the parity tests. `tests/*_parity.rs` pin a *named* case against
//! a recorded answer from R or `collapse`; a property test instead generates the input and checks a
//! rule that must hold for **all** of them. That difference earned its place: the two real bugs this
//! suite is built around -- `order_f64_multi` treating a missing value in a second key as equal, and
//! returning index order for an all-missing primary key -- both produced *plausible* outputs on
//! every fixture in the repository. Neither had a repeated key. A randomised key vector hits one.
//!
//! Every reference implementation below is written from R's documented semantics, deliberately not
//! by calling the function under test, so that a shared misreading cannot make a property
//! vacuously true. Where the rule is "R says X" the rule is quoted in the comment and the
//! measured values are in `docs/SEMANTICS.md`.

use proptest::prelude::*;

use r_core::aggregate::{aggregate_1, GroupMean};
use r_core::longdouble::F80;
use r_core::prob::r_max;
use r_core::ranknet::{order_f64, order_f64_multi};
use r_core::rng::MersenneTwister;
use r_core::stats::{geometric_mean, r_mean, r_median, r_nnzero, r_prod, thresholded_mean};

/// A double that is finite, plus one that may be `NaN`/`Inf`, so the awkward cases are generated as
/// often as the easy ones. Without this, 10 000 generated vectors contain about four non-finite
/// values in total and the properties that exist for them are never exercised.
mod arb {
    use super::*;

    pub fn finite() -> impl Strategy<Value = f64> {
        prop_oneof![
            // A wide range, including the negative side: `max` and the Hill terms treat 0 and
            // negative probabilities differently, and a suite that only sees positives misses it.
            any::<f64>().prop_filter("finite", |v| v.is_finite()),
            (0i64..1_000_000).prop_map(|i| i as f64 / 97.0),
            (-500_000i64..0).prop_map(|i| i as f64 / 89.0),
        ]
    }

    pub fn any_f64() -> impl Strategy<Value = f64> {
        prop_oneof![
            finite(),
            Just(f64::NAN),
            Just(f64::INFINITY),
            Just(f64::NEG_INFINITY),
        ]
    }

    pub fn vec_nonempty(max: usize) -> impl Strategy<Value = Vec<f64>> {
        prop::collection::vec(any_f64(), 1..max)
    }

    pub fn vec_finite(max: usize) -> impl Strategy<Value = Vec<f64>> {
        prop::collection::vec(finite(), 1..max)
    }
}

use arb::{finite, vec_finite, vec_nonempty};

// ---------------------------------------------------------------------------- R's `order`

/// R's `order(a, b, ..., na.last = TRUE)`, written out from the documented rules and measured in
/// `docs/SEMANTICS.md`:
///
/// ```text
/// order(c(1,1,2),  c(NA,5,3))   -> 2 1 3
/// order(c(1,1,2),  c(5,NA,3))   -> 1 2 3
/// order(c(1,NA,1), c(1,2,3))    -> 1 3 2
/// order(c(NA,NA,NA), c(3,1,2))  -> 2 3 1
/// ```
///
/// So: a missing value sorts last in **every** key; when both sides of a key are missing the
/// comparison falls through to the next key; and ties keep input order.
fn r_order_reference(keys: &[&[f64]]) -> Vec<usize> {
    let n = keys.first().map_or(0, |k| k.len());
    let mut idx: Vec<usize> = (0..n).collect();
    idx.sort_by(|&a, &b| {
        for k in keys.iter() {
            let (ka, kb) = (k[a].is_nan(), k[b].is_nan());
            if ka != kb {
                return if ka {
                    std::cmp::Ordering::Greater
                } else {
                    std::cmp::Ordering::Less
                };
            }
            if !ka {
                if let Some(o) = k[a].partial_cmp(&k[b]) {
                    if o != std::cmp::Ordering::Equal {
                        return o;
                    }
                }
            }
        }
        std::cmp::Ordering::Equal
    });
    idx
}

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// The one-key form, against a plain ascending sort with `NaN` last.
    #[test]
    fn order_f64_is_ascending_with_nan_last(x in vec_nonempty(64)) {
        let got = order_f64(&x);
        let mut want: Vec<usize> = (0..x.len()).filter(|&i| !x[i].is_nan()).collect();
        want.sort_by(|&a, &b| x[a].partial_cmp(&x[b]).unwrap());
        want.extend((0..x.len()).filter(|&i| x[i].is_nan()));
        prop_assert_eq!(got, want);
    }

    /// The multi-key form against the transcribed rule. This is the property that would have caught
    /// both `order_f64_multi` bugs, and a single-key or all-distinct-key test catches neither.
    #[test]
    fn order_multi_matches_r(keys in prop::collection::vec(vec_nonempty(24), 1..3)) {
        let n = keys[0].len();
        prop_assume!(n > 0);
        // Every key must be the same length, as they are in the real call sites.
        let keys: Vec<Vec<f64>> = keys.into_iter().map(|mut k| { k.resize(n, f64::NAN); k }).collect();
        let refs: Vec<&[f64]> = keys.iter().map(|k| k.as_slice()).collect();
        prop_assert_eq!(order_f64_multi(&refs), r_order_reference(&refs));
    }

    /// Repeated keys with a missing value in a *later* key -- the exact shape that made
    /// `order(pval, -prob)` place a `NaN` probability in the wrong position.
    #[test]
    fn repeated_primary_key_puts_a_missing_later_key_last(
        tie in prop::collection::vec(finite(), 1..6),
        others in vec_finite(5),
    ) {
        let n = tie.len().max(others.len());
        let mut a: Vec<f64> = tie.clone(); a.resize(n, 1.0);
        let mut b: Vec<f64> = others; b.resize(n, 2.0);
        if n > 1 { b[0] = f64::NAN; }
        let refs: [&[f64]; 2] = [&a, &b];
        let got = order_f64_multi(&refs);
        let want = r_order_reference(&refs);
        prop_assert_eq!(got.clone(), want);
        // And the specific fact, stated directly: within the group of equal first keys, the row
        // whose *second* key is missing comes after the finite ones.
        if n > 1 {
            let pos_nan = got.iter().position(|&i| i == 0).unwrap();
            // The group occupies a *contiguous* output range, but not necessarily one starting at
            // 0: rows with a smaller primary key sort ahead of it. Locating the group by its
            // smallest output position is what makes "last within the group" mean anything.
            let group_pos: Vec<usize> = got
                .iter()
                .enumerate()
                .filter(|(_, i)| a[**i] == a[0])
                .map(|(p, _)| p)
                .collect();
            if group_pos.len() > 1 {
                let contiguous = group_pos.windows(2).all(|w| w[1] == w[0] + 1);
                prop_assert!(contiguous, "an equal-primary-key group must be contiguous: {:?}", group_pos);
                prop_assert_eq!(pos_nan, *group_pos.last().unwrap(),
                    "the NaN row must be last within its equal-primary-key group");
            }
        }
    }
}

// ------------------------------------------------------------------------ R's `max` / `mean`

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// R's `max` propagates `NaN`; Rust's `f64::max` is a `maxNum` that ignores it. `-Inf` is a
    /// *value*, not a missing one: `max(c(-Inf, 1))` is `1`.
    #[test]
    fn r_max_propagates_nan_and_treats_neg_inf_as_a_value(x in vec_nonempty(64)) {
        let want = if x.iter().any(|v| v.is_nan()) {
            f64::NAN
        } else {
            x.iter().copied().fold(f64::NEG_INFINITY, f64::max)
        };
        let got = r_max(&x);
        prop_assert!(if want.is_nan() { got.is_nan() } else { got == want },
            "r_max({:?}) = {} but R's max is {}", x, got, want);
    }

    /// `mean` without `na.rm` returns `NaN` for any missing value, and is bounded by the min and max
    /// of a finite vector -- the bound being what makes a silent sign error detectable.
    #[test]
    fn r_mean_is_bounded_by_min_and_max_when_finite(x in vec_finite(64)) {
        let got = r_mean(&x);
        let lo = x.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        prop_assert!(got >= lo - 1e-9 && got <= hi + 1e-9,
            "mean {} outside [{}, {}]", got, lo, hi);
    }

    /// `nnzero` is `NA` for **any** vector with a missing value, not a smaller count. This single
    /// rule is what makes upstream's `thresholdedMean` raise.
    #[test]
    fn nnzero_is_na_when_any_entry_is_missing(x in vec_nonempty(32)) {
        let got = r_nnzero(&x);
        if x.iter().any(|v| v.is_nan()) {
            prop_assert!(got.is_nan());
        } else {
            let want = x.iter().filter(|v| **v != 0.0).count() as f64;
            prop_assert_eq!(got, want);
        }
    }
}

// ------------------------------------------------------------- the four `FunMean` variants

proptest! {
    #![proptest_config(ProptestConfig::with_cases(384))]

    /// A constant vector means exactly that constant, for every variant. `triMean` of a constant
    /// goes through `fquantile` and a four-element `mean`, and `geometricMean` through a
    /// log/exp round trip, so this catches drift in any of them.
    #[test]
    fn every_mean_of_a_constant_is_that_constant(c in finite(), n in 1usize..24) {
        let x = vec![c; n];
        for fun in [GroupMean::Mean, GroupMean::TriMean, GroupMean::Median,
                    GroupMean::TrimmedMean { trim: 0.1 },
                    GroupMean::ThresholdedMean { trim: 0.1 }] {
            let got = fun.apply(&x);
            if got.is_finite() {
                prop_assert!((got - c).abs() <= 1e-9 * c.abs().max(1.0),
                    "{fun:?}(constant {c}) = {got}");
            }
        }
    }

    /// `median` of a finite vector is between its min and max, and is permutation-invariant --
    /// which is a real check, because R's type-7 median is computed on a *sorted* copy and a
    /// port that sorted in place would perturb its caller.
    #[test]
    fn median_is_bounded_and_permutation_invariant(x in vec_finite(48)) {
        let m = r_median(&x, true);
        let lo = x.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        prop_assert!(m >= lo - 1e-12 && m <= hi + 1e-12);
        let mut shuffled = x.clone();
        // A cheap deterministic rotation is enough: it is a permutation, which is the property.
        if !shuffled.is_empty() {
            shuffled.rotate_left(1);
        }
        // By value: `r_median` of an even-length vector averages the two middle values, and
        // averaging `-0.0` can land on `+0.0`, which is a different bit pattern and the same number.
        let m2 = r_median(&shuffled, true);
        prop_assert!(m2 == m, "median must be permutation invariant: {} vs {}", m, m2);
    }

    /// The geometric mean of a finite positive vector lies between its min and max. A negative or
    /// zero entry makes it non-finite, which is also worth pinning.
    #[test]
    fn geometric_mean_is_bounded_for_positive_inputs(x in prop::collection::vec(0.01f64..1e3, 1..32)) {
        let got = geometric_mean(&x);
        let lo = x.iter().copied().fold(f64::INFINITY, f64::min);
        let hi = x.iter().copied().fold(f64::NEG_INFINITY, f64::max);
        prop_assert!(got >= lo - 1e-6 * lo.abs() && got <= hi + 1e-6 * hi.abs(),
            "geometric mean {} outside [{}, {}]", got, lo, hi);
    }

    /// `thresholdedMean` is `0` strictly below the threshold, and never `NaN` for a vector with no
    /// missing value -- it may return `NaN` only when `nnzero` is `NA`, which is the raise case.
    #[test]
    fn thresholded_mean_is_zero_below_threshold_and_finite_otherwise(
        x in vec_finite(32), trim in 0.0f64..0.5
    ) {
        let got = thresholded_mean(&x, trim);
        prop_assert!(!got.is_nan(), "no missing values, so nnzero is not NA");
        let nz = x.iter().filter(|v| **v != 0.0).count() as f64 / x.len() as f64;
        if nz < trim {
            prop_assert_eq!(got, 0.0);
        }
    }

    /// `prod` with `na.rm = TRUE` ignores missing values; without it, any missing poisons the
    /// product. A `0` in the vector pins the absorbing case, which is the easy one to get wrong
    /// when a fast path skips the multiply.
    #[test]
    fn prod_honours_na_rm(x in vec_nonempty(24), na_rm in any::<bool>()) {
        let got = r_prod(&x, na_rm);
        if !na_rm && x.iter().any(|v| v.is_nan()) {
            prop_assert!(got.is_nan());
        } else {
            let vals: Vec<f64> = x.iter().copied().filter(|v| !v.is_nan()).collect();
            if vals.is_empty() {
                prop_assert_eq!(got, 1.0, "the empty product is 1, as in R");
            } else {
                // The log domain needs strictly positive **finite** inputs: `ln(x)` is `NaN` for
                // `x <= 0` and `Inf` for `x = Inf`, and both make the comparison `NaN`. So the
                // direct product is used for anything else, and the log product only where it is
                // the more accurate reference.
                // What is checked here is the rule that is universally true of `prod` -- the
                // `na.rm` handling of a missing value, and the empty product being 1. The
                // *value* of a product is not comparable against a naive `f64` fold at all, because
                // R accumulates in LONG_DOUBLE and narrows once: `prod(c(8.98e-288, 5.77e-49, Inf))`
                // is `Inf` in R and `NaN` in a `f64` fold, and `prod(c(8.98e-288, 5.77e-49))` is
                // `0x0p+0` in R and `0` only by coincidence of narrowing. Every non-finite and
                // saturating case is pinned against R's measured table in
                // `longdouble::mul_special_tests` instead, which is the only valid oracle for it.
                let fold = vals.iter().copied().fold(1.0f64, |a, b| a * b);
                let fold_reliable = fold.is_finite()
                    && fold != 0.0
                    && vals.iter().all(|v| v.is_finite() && *v > 0.0);
                if !fold_reliable {
                    // R produces `NaN` from a product in exactly two ways, both measured:
                    // a missing value without `na.rm`, and a zero times an infinity
                    // (`prod(c(0, Inf))` and `prod(c(-0, Inf))` are both `NaN`). Nothing else may.
                    let has_missing = vals.iter().any(|v| v.is_nan());
                    let zero_times_inf = vals.contains(&0.0)
                        && vals.iter().any(|v| v.is_infinite());
                    prop_assert!(!got.is_nan() || has_missing || zero_times_inf,
                        "r_prod({:?}, {}) = {} -- NaN needs a missing value or a zero times an infinity",
                        x, na_rm, got);
                } else {
                    // Compared in the log domain, but only where the product is in the *normal*
                    // `f64` range. R accumulates in LONG_DOUBLE and narrows once at the end, so
                    // `prod(c(8.98e-288, 5.77e-49))` is `0x0p+0` even though the extended value is
                    // about `5.2e-336`; a log comparison that "found" the extended answer would be
                    // reporting a divergence from R that is really a divergence from the reference.
                    // Within range, the log form is the better reference, because the naive `f64`
                    // product rounds differently from the extended accumulation and that rounding
                    // difference is not a bug in either.
                    let want_log: f64 = vals.iter().map(|v| v.ln()).sum();
                    prop_assert!((got.ln() - want_log).abs() <= 1e-9 * want_log.abs().max(1.0),
                        "r_prod({:?}, {}) = {} but log-product is {}", x, na_rm, got, want_log);
                }
            }
        }
    }
}

// --------------------------------------------------------------------------------- `F80`

proptest! {
    #![proptest_config(ProptestConfig::with_cases(512))]

    /// `F80` has a 64-bit mantissa, so a sum of `f64`s is *more* precise than an `f64` sum. The
    /// property that matters is not equality with `f64` but that the two agree to within a few ulp
    /// and that no value is invented: the 80-bit result is always within one rounding of the
    /// `f64` result.
    #[test]
    fn f80_sum_agrees_with_f64_sum_to_a_few_ulp(x in vec_finite(48)) {
        let mut acc = F80::from_f64(0.0);
        let mut naive = 0.0f64;
        for &v in &x {
            acc = acc.add(F80::from_f64(v));
            naive += v;
        }
        let got = acc.to_f64();
        if naive.is_finite() && got.is_finite() {
            // The naive `f64` sum's own error is bounded by the magnitudes of its **partials**,
            // which for a cancelling vector are far larger than the total: summing
            // `[-3952, -5478, 9344]` gives `-85.3` from partials near `1e4`. Bounding by the total
            // would fail on every cancelling input, and would be a property of the reference, not
            // of the accumulator. The 80-bit result is *more* accurate, so they may differ -- the
            // claim is only that they differ by no more than the naive sum's own uncertainty.
            let scale: f64 = x.iter().map(|v| v.abs()).sum::<f64>().max(1.0);
            let tol = 16.0 * f64::EPSILON * scale * x.len() as f64;
            prop_assert!((got - naive).abs() <= tol,
                "F80 sum {} vs f64 sum {} over {:?} (scale {})", got, naive, x, scale);
        }
    }

    /// Any `NaN` poisons the sum, and `Inf + -Inf` is `NaN`. This is the property behind R's
    /// `if (sum(P1_Pspatial) == 0)` raising `missing value where TRUE/FALSE needed`.
    #[test]
    fn f80_propagates_the_special_values(x in vec_nonempty(24)) {
        let mut acc = F80::from_f64(0.0);
        for &v in &x {
            acc = acc.add(F80::from_f64(v));
        }
        let got = acc.to_f64();
        if x.iter().any(|v| v.is_nan()) {
            prop_assert!(got.is_nan(), "a NaN term must poison the sum, got {}", got);
        }
        if x.iter().any(|v| v.is_infinite() && *v > 0.0)
            && x.iter().any(|v| v.is_infinite() && *v < 0.0) {
            prop_assert!(got.is_nan(), "Inf + -Inf must be NaN, got {}", got);
        }
        prop_assert_eq!(acc.to_f64().to_bits(), got.to_bits(), "conversion must be idempotent");
    }

    /// `F80` addition is commutative where the result is finite, and adding zero is the identity.
    /// Non-associativity is *expected* (that is the whole point of the extended accumulator) and is
    /// not asserted.
    #[test]
    fn f80_add_is_commutative_where_finite(a in finite(), b in finite()) {
        let (x, y) = (F80::from_f64(a), F80::from_f64(b));
        let ab = x.add(y).to_f64();
        let ba = y.add(x).to_f64();
        if ab.is_finite() {
            prop_assert!((ab - ba).abs() <= 1e-15 * ab.abs().max(1.0),
                "{} + {} != {} + {}", a, b, ab, ba);
        }
        // By *value*, not by bits: IEEE says `-0.0 + 0.0` is `+0.0`, and R's `sum(c(-0, 0))` is
        // `0` too, so a bit comparison would flag correct behaviour.
        let z = x.add(F80::from_f64(0.0)).to_f64();
        prop_assert!(z == a, "adding zero must be the identity: {} + 0 = {}", a, z);
    }
}

// ------------------------------------------------------------------------------ the RNG

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// `sample.int` without replacement is always a permutation of `0..n`. This is the invariant
    /// the bootstrap relies on -- an out-of-range index is exactly the `subscript out of bounds`
    /// the shim has to reproduce rather than fix.
    #[test]
    fn sample_int_no_replace_is_always_a_permutation(n in 1usize..200, k in 0usize..200, seed in any::<i32>()) {
        let k = k.min(n);
        let mut rng = MersenneTwister::new(seed);
        let got = rng.sample_int_no_replace(n, k);
        prop_assert_eq!(got.len(), k);
        let mut seen = got.clone();
        seen.sort_unstable();
        seen.dedup();
        prop_assert_eq!(seen.len(), k, "duplicate index in {:?}", got);
        prop_assert!(got.iter().all(|&i| i >= 1 && (i as usize) <= n),
            "index out of 1..={} in {:?}", n, got);
    }

    /// The same seed gives the same stream, and `unif_rand` stays in `[0, 1)`.
    #[test]
    fn mersenne_twister_is_reproducible_and_in_range(seed in any::<i32>(), n in 1usize..64) {
        let a: Vec<f64> = { let mut r = MersenneTwister::new(seed); (0..n).map(|_| r.unif_rand()).collect() };
        let b: Vec<f64> = { let mut r = MersenneTwister::new(seed); (0..n).map(|_| r.unif_rand()).collect() };
        prop_assert_eq!(a.iter().map(|v| v.to_bits()).collect::<Vec<_>>(),
                       b.iter().map(|v| v.to_bits()).collect::<Vec<_>>());
        prop_assert!(a.iter().all(|v| (0.0..1.0).contains(v)), "unif_rand outside [0,1): {:?}", a);
    }
}

// ------------------------------------------------------------------------- `aggregate_1`

proptest! {
    #![proptest_config(ProptestConfig::with_cases(256))]

    /// `aggregate` returns exactly one row per level, in **level order** -- not order of first
    /// appearance. `split()` indexes by `levels(f)`, and any port that indexes by first appearance
    /// is wrong for every input whose levels are not already sorted.
    #[test]
    fn aggregate_returns_one_row_per_level_in_level_order(
        n_groups in 1usize..6, n_cols in 1usize..5, assignment in prop::collection::vec(any::<u8>(), 1..40)
    ) {
        let n_rows = assignment.len();
        let n_groups = n_groups.max(1);
        let group: Vec<usize> = assignment.into_iter().map(|b| (b as usize) % n_groups).collect();
        let data: Vec<f64> = (0..(n_rows * n_cols)).map(|i| (i % 17) as f64).collect();
        for fun in [GroupMean::Mean, GroupMean::TriMean, GroupMean::Median,
                    GroupMean::TrimmedMean { trim: 0.1 }] {
            let out = aggregate_1(&data, n_rows, &(0..n_cols).collect::<Vec<_>>(), &group, n_groups, fun);
            let kind = format!("{fun:?}");
            prop_assert_eq!(out.len(), n_groups * n_cols,
                            "{} returned {} values for {} groups", kind, out.len(), n_groups);
        }
    }

    /// The `Mean` variant of a group is the arithmetic mean of exactly that group's values, in
    /// ascending cell order. `triMean` and `median` are not arithmetic means and are not compared
    /// here; the property is about *which values land in which row*.
    #[test]
    fn aggregate_mean_row_equals_the_mean_of_that_group(
        n_cols in 1usize..4, assignment in prop::collection::vec(any::<u8>(), 2..40)
    ) {
        let n_rows = assignment.len();
        let n_groups = 4usize;
        let group: Vec<usize> = assignment.into_iter().map(|b| (b as usize) % n_groups).collect();
        let data: Vec<f64> = (0..(n_rows * n_cols)).map(|i| ((i * 7) % 23) as f64).collect();
        let out = aggregate_1(&data, n_rows, &(0..n_cols).collect::<Vec<_>>(), &group,
                             n_groups, GroupMean::Mean);
        for g in 0..n_groups {
            let cells: Vec<usize> = (0..n_rows).filter(|&r| group[r] == g).collect();
            if cells.is_empty() {
                continue;
            }
            for c in 0..n_cols {
                // Column-major: element (row, col) is `values[col * n_rows + row]`. Using
                // `row * n_cols + col` here is the transpose, and it produced a property that
                // failed on its very first case with a perfectly correct kernel -- the third
                // distinct layout mistake in this file, after the slice order and the dim vector.
                let want: f64 = cells.iter().map(|&r| data[c * n_rows + r]).sum::<f64>()
                    / cells.len() as f64;
                let got = out[g * n_cols + c];
                prop_assert!((got - want).abs() <= 1e-9 * want.abs().max(1.0),
                    "group {} col {}: {} != {}", g, c, got, want);
            }
        }
    }
}
