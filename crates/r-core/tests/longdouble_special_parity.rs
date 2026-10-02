//! `F80` must carry `NaN` and `Inf` rather than collapsing them to zero.
//!
//! R's `sum()` propagates both: `sum(c(1, NaN))` is `NaN`, `sum(c(1, Inf))` is `Inf`, and
//! `sum(c(Inf, -Inf))` is `NaN`. A long-double accumulator that treats a non-finite term as zero
//! cannot produce those answers, and the difference is observable far beyond `sum()`: R's
//! `if (sum(P1_Pspatial) == 0)` raises "missing value where TRUE/FALSE needed" on the `NaN` case,
//! while the `Inf` case is merely non-zero and proceeds.
//!
//! Every reduction in this crate goes through `F80`, so this is not a corner of `sum()` -- it is
//! the behaviour of every accumulated mean, `trimmed_mean`, `tri_mean` and per-pathway total.

use r_core::longdouble::F80;

fn acc(v: &[f64]) -> f64 {
    v.iter()
        .fold(F80::ZERO, |a, x| a.add(F80::from_f64(*x)))
        .to_f64()
}

#[test]
fn r_max_propagates_nan_but_not_negative_infinity() {
    let mx = |v: &[f64]| r_core::prob::r_max(v);
    assert_eq!(mx(&[1.0, 2.0, 3.0]), 3.0);
    assert!(
        mx(&[1.0, f64::NAN]).is_nan(),
        "R's max propagates NaN: max(c(1, NA)) is NA"
    );
    assert!(
        mx(&[f64::NEG_INFINITY, 1.0]) == 1.0,
        "-Inf is a value, not a missing one"
    );
    assert!(mx(&[f64::NEG_INFINITY]).is_infinite(), "max(-Inf) is -Inf");
    assert!(mx(&[f64::INFINITY, 1.0]).is_infinite());
    // Rust's own `f64::max` is a `maxNum` and would return 1.0 for the first NaN case, which is
    // exactly why this function exists rather than a fold over `f64::max`.
    assert_eq!(
        1.0_f64.max(f64::NAN),
        1.0,
        "the Rust behaviour this must not copy"
    );
}

#[test]
fn sum_propagates_nan_and_infinity_like_r() {
    assert_eq!(acc(&[1.0, 2.0]), 3.0);
    assert!(acc(&[1.0, f64::NAN]).is_nan(), "sum(c(1, NA)) is NA in R");
    assert!(
        acc(&[0.0, f64::NAN]).is_nan(),
        "and sum(c(0, NA)) is NA too -- 0 + NaN is NaN"
    );
    assert!(
        acc(&[1.0, f64::INFINITY]).is_infinite() && acc(&[1.0, f64::INFINITY]).is_sign_positive()
    );
    assert!(
        acc(&[1.0, f64::NEG_INFINITY]).is_infinite()
            && acc(&[1.0, f64::NEG_INFINITY]).is_sign_negative()
    );
    assert!(
        acc(&[f64::INFINITY, f64::NEG_INFINITY]).is_nan(),
        "Inf + -Inf is NaN in R"
    );
    assert_eq!(
        acc(&[0.0, 0.0]),
        0.0,
        "and 0 is still 0, with the sign R gives it"
    );
}

#[test]
fn scaling_by_a_non_finite_max_produces_nan_everywhere() {
    // `data/max(data)` with `max == 0` is `0/0`, which is `NaN` in every entry. This is the step
    // that makes an all-zero matrix raise upstream's "missing value where TRUE/FALSE needed".
    let scaled = r_core::prob::scale_by_max(&[0.0, 0.0, 0.0]);
    assert!(scaled.iter().all(|v| v.is_nan()));
    let scaled = r_core::prob::scale_by_max(&[0.0, 1.0, 2.0]);
    assert_eq!(scaled, vec![0.0, 0.5, 1.0]);
    // `Inf` in the data with a finite maximum: the `Inf` entry becomes `NaN` and the rest 0.
    let scaled = r_core::prob::scale_by_max(&[1.0, 2.0, f64::INFINITY]);
    assert_eq!(scaled[0], 0.0);
    assert!(scaled[2].is_nan());
}

#[test]
fn finite_arithmetic_is_unchanged() {
    // The sentinels must not have disturbed the finite path: exact 80-bit accumulation, rounded
    // once to `f64`.
    assert_eq!(
        acc(&[1e16, 1.0, -1e16]),
        1.0,
        "80-bit significand keeps the middle term"
    );
    // 0.1 + 0.2 + 0.3 in 80 bits, rounded once: the correctly-rounded double is 0.6, not the
    // 0.6000000000000001 that naive f64 accumulation gives.
    assert_eq!(
        acc(&[0.1, 0.2, 0.3]),
        0.6,
        "80-bit accumulation rounds once, at the end"
    );
    // The point of the test above: this is the answer an `f64` accumulator gives, and the 80-bit
    // one must *not* reproduce it. `clippy::assertions_on_constants` reads a constant truth as a
    // mistake; here it is the specification. Bound to a local because a `#[allow]` on the statement
    // does not reach spans the macro's own expansion produces.
    let f64_accumulated = 0.1f64 + 0.2 + 0.3;
    assert!(f64_accumulated > 0.6, "the f64 result this must not copy");
    assert!(F80::ZERO.is_zero());
    assert!(!F80::NAN.is_zero(), "NaN is not zero -- that was the bug");
    assert!(!F80::inf(true).is_zero());
    assert_eq!(F80::inf(false).to_f64(), f64::INFINITY);
    assert_eq!(F80::inf(true).to_f64(), f64::NEG_INFINITY);
    assert!(
        F80::NAN.neg().to_f64().is_nan(),
        "negating NaN is still NaN"
    );
    assert_eq!(F80::inf(false).neg().to_f64(), f64::NEG_INFINITY);
}

#[test]
fn round_trips_survive_the_specials() {
    // `-0.0` is deliberately absent: `F80::ZERO` has no sign, so `from_f64(-0.0)` is `+0.0`. That
    // matches R, where `sum(c(-0.0))` is `0` -- `sum` starts from `0` and `0 + -0` is `+0` under
    // round-to-nearest -- so an accumulator that preserved the sign would be *less* faithful.
    for v in [
        f64::NAN,
        f64::INFINITY,
        f64::NEG_INFINITY,
        0.0,
        1.5,
        -2.25,
        f64::MIN_POSITIVE,
    ] {
        let back = F80::from_f64(v).to_f64();
        if v.is_nan() {
            assert!(back.is_nan(), "NaN must round-trip as NaN");
        } else {
            assert_eq!(back.to_bits(), v.to_bits(), "{v:e} did not round-trip");
        }
    }
}
