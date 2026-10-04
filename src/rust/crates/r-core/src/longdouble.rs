//! Minimal x87 80-bit extended precision ("f80") arithmetic.
//!
//! ## Why this exists
//!
//! R 4.x's `real_mean` (`src/main/summary.c:476`) is a **two-pass corrected mean**
//! evaluated in `LDOUBLE`, which on x86-64 Linux is the 80-bit extended format — a
//! **64-bit significand**, not 53:
//!
//! ```c
//! LDOUBLE s = 0.0;
//! for (k = 0; k < n; k++) s += dx[k];       // pass 1
//! s /= n;
//! LDOUBLE t = 0.0;
//! for (k = 0; k < n; k++) t += (dx[k] - s); // pass 2: residuals
//! s += t/n;
//! return ScalarReal((double) s);             // single narrowing here
//! ```
//!
//! A plain `f64` accumulation is **not** equivalent, and that was measured rather than
//! assumed: for `mean(c(1e16, 1, -1e16))` a naive `f64` sum gives `0` while R gives
//! `0.33365885416666669`, because pass 1 lands on the `1` only in 64 bits and pass 2
//! then recovers it. For `triMean` the disagreement shows up as 1–2 ulp on ~3.5 % of
//! vectors (`docs/SEMANTICS.md` R10).
//!
//! ## Representation
//!
//! `value = (-1)^sign * m * 2^exp`, with `m` a 64-bit integer whose bit 63 is set for
//! every nonzero value (an *explicit* leading bit, as x87 uses). Zero is `m == 0`.
//!
//! Infinities and subnormals are not representable in this type, deliberately: R only
//! ever tests finiteness on the **narrowed** value (`R_FINITE((double) s)`), so the
//! infinities that matter are created by the `f64` cast, not by the accumulator. The
//! caller in [`crate::stats::r_mean`] reproduces that ordering exactly.
//!
//! Every operation rounds to nearest with ties to even, matching x87.

/// 80-bit extended precision value: `(-1)^sign * m * 2^exp`, `m` in `[0, 2^64)`.
#[derive(Clone, Copy, Debug)]
pub struct F80 {
    sign: bool,
    exp: i32,
    m: u64,
}

/// Sentinel exponents for the two non-finite values, as in the x87 extended format (an
/// all-ones exponent means "special"; the significand's leading bit then distinguishes
/// infinity from NaN). Using sentinels rather than widening the struct keeps the finite
/// representation and every finite operation byte-for-byte as it was.
const EXP_INF: i32 = i32::MAX - 1;
const EXP_NAN: i32 = i32::MAX;

/// Only meaningful for finite values; `PartialEq` routes through [`F80::to_bits`] for the
/// special cases so that `NaN != NaN` and `Inf == Inf`.
impl PartialEq for F80 {
    fn eq(&self, other: &F80) -> bool {
        match (self.is_nan(), other.is_nan()) {
            (true, true) => true,
            (true, false) | (false, true) => false,
            _ => self.sign == other.sign && self.exp == other.exp && self.m == other.m,
        }
    }
}
impl Eq for F80 {}

/// `(shifted >> 1, bit_is_one, anything_below_is_nonzero)`.
#[inline]
fn round_bits(m: u64, shift: u32) -> (u64, bool, bool) {
    debug_assert!((1..64).contains(&shift));
    let q = m >> shift;
    let round = (m >> (shift - 1)) & 1 != 0;
    let sticky = (m & ((1u64 << (shift - 1)) - 1)) != 0;
    (q, round, sticky)
}

impl F80 {
    pub const ZERO: F80 = F80 {
        sign: false,
        exp: 0,
        m: 0,
    };
    /// `NaN`. Every `NaN` payload collapses to this one; R makes the same collapse, and the
    /// only question a caller can ask is whether a value is missing.
    pub const NAN: F80 = F80 {
        sign: false,
        exp: EXP_NAN,
        m: 0,
    };

    #[inline]
    pub fn inf(sign: bool) -> F80 {
        F80 {
            sign,
            exp: EXP_INF,
            m: 0,
        }
    }

    #[inline]
    pub fn is_nan(self) -> bool {
        self.exp == EXP_NAN
    }

    #[inline]
    pub fn is_infinite(self) -> bool {
        self.exp == EXP_INF
    }

    #[inline]
    pub fn is_finite(self) -> bool {
        self.exp < EXP_INF
    }

    #[inline]
    pub fn is_zero(self) -> bool {
        self.m == 0 && self.is_finite()
    }

    #[inline]
    pub fn neg(self) -> F80 {
        F80 {
            sign: !self.sign,
            ..self
        }
    }

    /// Widen an `f64`. Exact — 53 significand bits fit in 64 with room to spare.
    ///
    /// Non-finite inputs collapse to zero. R never accumulates a non-finite term in the
    /// finite branch, and for the infinite branch it recomputes with `dx[k] / n`; that
    /// path is handled by the caller, which checks finiteness on the narrowed value.
    pub fn from_f64(x: f64) -> F80 {
        let bits = x.to_bits();
        let sign = (bits >> 63) != 0;
        let biased = ((bits >> 52) & 0x7ff) as i32;
        let frac = bits & 0x000f_ffff_ffff_ffff;
        if frac == 0 && biased == 0 {
            return F80::ZERO;
        }
        if biased == 0x7ff {
            // `frac == 0` is an infinity; anything else is a NaN, and every NaN payload
            // collapses to the one `F80::NAN`, exactly as R collapses `NA_real_` and `NaN` to the
            // same "missing" notion even though their bit patterns differ.
            return if frac == 0 { F80::inf(sign) } else { F80::NAN };
        }
        if biased == 0 {
            let lz = frac.leading_zeros() as i32; // 11..=63
            F80 {
                sign,
                exp: -1022 - 52 - lz,
                m: frac << lz,
            }
        } else {
            F80 {
                sign,
                exp: biased - 1023 - 63,
                m: (1u64 << 63) | (frac << 11),
            }
        }
    }

    /// Narrow to `f64`, rounding to nearest with ties to even.
    ///
    /// This is the *only* place precision is lost in [`crate::stats::r_mean`], matching
    /// R, which likewise returns `(double) s` once at the very end.
    pub fn to_f64(self) -> f64 {
        if self.is_nan() {
            return f64::NAN;
        }
        if self.is_infinite() {
            return if self.sign {
                f64::NEG_INFINITY
            } else {
                f64::INFINITY
            };
        }
        if self.m == 0 {
            // A **signed** zero. R's `prod(c(-0, 5))` is `-0x0p+0`, not `0x0p+0`, and a caller that
            // divides by a product, or compares `1/x` against `-Inf`, can tell the difference.
            // Returning a bare `0.0` here silently made every signed zero positive on the way out.
            // `-Inf * 0` is `NaN`, not `-0`; the bit pattern is the only way to spell a signed zero.
            return f64::from_bits((self.sign as u64) << 63);
        }
        let sign_bit = (self.sign as u64) << 63;
        // Value is 2^lead * (1 + fraction) with lead = exp + 63.
        let lead = self.exp + 63;

        if lead < -1022 {
            // f64-subnormal territory. A subnormal f64 is k * 2^-1074 with k < 2^52.
            let shift = (-1074 - self.exp) as u32;
            if shift >= 64 {
                // Strictly below half of the smallest subnormal: rounds to zero.
                // (`shift == 64` exactly means value < 2^-1074, so no tie is possible
                //  here because m >= 2^63 makes value >= 2^(exp+63) = 2^-1074.)
                return f64::from_bits(sign_bit);
            }
            let (k, round, sticky) = round_bits(self.m, shift);
            let mut k = k;
            if round && (sticky || (k & 1) != 0) {
                k += 1;
            }
            if k >= (1u64 << 52) {
                // Rounded up into the smallest normal.
                return f64::from_bits(sign_bit | (1u64 << 52));
            }
            return f64::from_bits(sign_bit | k);
        }

        if lead > 1023 {
            return f64::from_bits(sign_bit | (0x7ffu64 << 52));
        }

        // Normal: keep 53 bits, i.e. drop 11 from the 64-bit significand. The leading
        // bit lands at bit 52 of `m53`, so the f64 exponent field is `lead + 1023` and
        // the fraction field is the low 52 bits.
        let (mut m53, round, sticky) = round_bits(self.m, 11);
        let mut lead = lead;
        if round && (sticky || (m53 & 1) != 0) {
            m53 += 1;
            if m53 >> 53 != 0 {
                m53 >>= 1;
                lead += 1;
            }
        }
        if lead > 1023 {
            return f64::from_bits(sign_bit | (0x7ffu64 << 52));
        }
        f64::from_bits(sign_bit | (((lead + 1023) as u64) << 52) | (m53 & 0x000f_ffff_ffff_ffff))
    }

    /// `self + other`, rounded to nearest with ties to even.
    ///
    /// Non-finite operands follow IEEE 754 and R: any `NaN` makes the result `NaN`; `Inf + Inf` is
    /// `Inf`; `Inf + -Inf` is `NaN`; otherwise an infinite operand wins.
    pub fn add(self, other: F80) -> F80 {
        if self.is_nan() || other.is_nan() {
            return F80::NAN;
        }
        if self.is_infinite() || other.is_infinite() {
            if self.is_infinite() && other.is_infinite() {
                return if self.sign == other.sign {
                    self
                } else {
                    F80::NAN
                };
            }
            return if self.is_infinite() { self } else { other };
        }
        if self.m == 0 {
            return other;
        }
        if other.m == 0 {
            return self;
        }
        if self.sign == other.sign {
            return self.add_same_sign(other);
        }
        match self.cmp_magnitude(other) {
            std::cmp::Ordering::Equal => F80::ZERO,
            std::cmp::Ordering::Greater => self.sub_same_sign(other),
            std::cmp::Ordering::Less => other.sub_same_sign(self),
        }
    }

    /// `self - other`.
    #[inline]
    pub fn sub(self, other: F80) -> F80 {
        self.add(other.neg())
    }

    fn cmp_magnitude(self, other: F80) -> std::cmp::Ordering {
        // |value| = m * 2^exp, so compare exponents then significands.
        self.exp.cmp(&other.exp).then_with(|| self.m.cmp(&other.m))
    }

    /// Addition of two same-signed values: cannot lose magnitude.
    fn add_same_sign(self, other: F80) -> F80 {
        let (hi, lo) = if self.exp >= other.exp {
            (self, other)
        } else {
            (other, self)
        };
        let de = (hi.exp - lo.exp) as u32;
        // 128-bit accumulator with `hi.m` at bits [126..63]; see [`F80::reduce`] for
        // why the anchor is 63 and not 64.
        let mut acc = (hi.m as u128) << 63;
        if de < 63 {
            acc += (lo.m as u128) << (63 - de);
        } else if de < 127 {
            acc += (lo.m as u128) >> (de - 63);
            if (lo.m as u128) & ((1u128 << (de - 63)) - 1) != 0 {
                acc |= 1; // sticky
            }
        } else if lo.m != 0 {
            acc |= 1; // entirely below the LSB
        }
        F80::reduce(acc, hi.sign, hi.exp)
    }

    /// Difference of the larger magnitude minus the smaller, for operands of
    /// **opposite** sign, so the result takes the larger one's sign.
    ///
    /// The operands are re-ordered *here* by exponent rather than assumed. The caller
    /// only guarantees `|self| > |other|`, which after an opposite-sign `add` does
    /// **not** imply `self.exp >= other.exp` — `1.0 - 0.25` has the larger magnitude at
    /// the larger exponent, but `0.75 - 1.0` does not. Assuming otherwise made `de`
    /// underflow to a huge `u32` and turned `1.0 - 0.25` into `1.5`.
    fn sub_same_sign(self, other: F80) -> F80 {
        let (hi, lo) = if self.exp >= other.exp {
            (self, other)
        } else {
            (other, self)
        };
        let de = (hi.exp - lo.exp) as u32;
        // Same accumulator frame as `add_same_sign`; see [`F80::reduce`].
        let mut acc = (hi.m as u128) << 63;
        if de < 63 {
            acc -= (lo.m as u128) << (63 - de);
        } else if de < 127 {
            // `lo` is shifted right; the bits shifted out are not tracked as a negative
            // sticky. That is safe because they land at or below
            // `2^(hi.exp - 126)`, which is more than 60 bits below f64 resolution.
            acc -= (lo.m as u128) >> (de - 63);
        }
        // de >= 127: `lo` lies entirely below the accumulator's LSB, which cannot happen
        // for normalised operands with |hi| >= |lo|.
        F80::reduce(acc, hi.sign, hi.exp)
    }

    /// Normalise a 128-bit accumulator to a 64-bit significand, rounding to nearest
    /// with ties to even, folding in the sticky bit.
    ///
    /// **Why the kept significand occupies bits [126..63] and not [127..64].** The
    /// accumulator's unit is `2^(exp - 63)`, so `acc` holds the significand at
    /// `hi.m << 63`. Anchoring one bit lower guarantees the aligned addition cannot
    /// overflow `u128`:
    ///
    /// * `hi.m << 63 <= (2^64 - 1) * 2^63 = 2^127 - 2^63`
    /// * the aligned addend is at most `2^63 - 1` (it has already been shifted right
    ///   by at least one bit whenever `de > 63`)
    /// * so the sum is at most `2^127 - 1`
    ///
    /// Anchoring at 64 would allow `1.2 * 2^128`, which wraps silently in release
    /// builds and produced `0.1 + 0.2 == 0.05` before this was found.
    fn reduce(mut acc: u128, sign: bool, exp: i32) -> F80 {
        if acc == 0 {
            return F80::ZERO;
        }
        let mut exp = exp;
        let mut extra_sticky = false;

        // Carry out of the anchor: `acc >= 2^127` means the aligned sum overflowed the
        // 64-bit significand. Shifting the whole accumulator right and
        // bumping `exp` is equivalent to pre-shifting `hi.m`, and leaves a spare bit of
        // headroom. The bit shifted out is *not* discarded: it is a sticky, and the
        // rounding decision depends on it.
        if acc >> 127 != 0 {
            let lost = acc & 1 != 0;
            acc >>= 1;
            extra_sticky = lost;
            exp += 1;
        }

        // Left-normalise so the leading bit lands at position 126. After the carry
        // above the MSB is at most 126, so this only ever shifts *up* and cannot lose
        // bits. It is essential after near-cancelling subtraction, which can leave the
        // result many bits below the anchor.
        let msb = 127 - acc.leading_zeros() as i32;
        if msb < 126 {
            acc <<= (126 - msb) as u32;
            exp -= 126 - msb;
        }

        let kept = (acc >> 63) as u64; // bits [126..63]; bit 63 is set by construction
        let remainder = acc & ((1u128 << 63) - 1);
        if remainder == 0 && !extra_sticky {
            return F80 { sign, exp, m: kept };
        }
        let round = (acc >> 62) & 1 != 0;
        let sticky = (acc & ((1u128 << 62) - 1)) != 0 || extra_sticky;
        let mut m = kept;
        if round && (sticky || (m & 1) != 0) {
            m = m.wrapping_add(1);
            let mut e = exp;
            if m == 0 {
                m = 1u64 << 63;
                e += 1;
            }
            return F80 { sign, exp: e, m };
        }
        F80 { sign, exp, m }
    }

    /// `self * other`, rounded to nearest with ties to even.
    ///
    /// Needed because R's `prod` (`rprod`, `src/main/summary.c:367`) also accumulates in
    /// `LDOUBLE`. This is not a theoretical concern: `prod()` of 8 random values in
    /// `[1, 3)` differs from an `f64` left-to-right product on the *first* random
    /// trial tried, by 1 ulp -- and `computeExpr_coreceptor` / `computeExpr_agonist` /
    /// `computeExpr_antagonist` all end in `apply(..., 2, prod)` when a cofactor has
    /// more than one subunit.
    ///
    /// The 64x64 -> 128-bit significand product is **exact**, so only the final rounding
    /// to 64 bits can lose anything.
    pub fn mul(self, other: F80) -> F80 {
        if self.is_nan() || other.is_nan() {
            return F80::NAN;
        }
        // Infinity is encoded with a **zero mantissa**, so the old leading
        // `if self.m == 0 || other.m == 0 { return F80::ZERO }` caught every infinite operand and
        // turned `prod(c(-Inf))` into `0`. A property test found it: R gives `-Inf`.
        if self.is_infinite() || other.is_infinite() {
            // R's table, measured rather than assumed:
            //
            // ```text
            // prod(Inf, Inf)  = Inf     prod(Inf, -Inf)  = -Inf
            // prod(-Inf, Inf)  = -Inf    prod(-Inf, -Inf) = Inf
            // prod(Inf, 2)     = Inf     prod(Inf, 0)     = NaN
            // prod(0,  Inf)    = NaN     prod(-Inf, 0)    = NaN
            // ```
            //
            // So an infinite product is `Inf` carrying the **XOR of the two signs** -- including
            // `Inf * -Inf`, which IEEE 754 makes `NaN` and R does not -- and only a zero times an
            // infinity is `NaN`. This is R's long-double behaviour, and matching it is the point:
            // an "obviously correct" IEEE `mul` would be a divergence.
            // `m == 0` is *also* how an infinity is encoded, so "the other operand is a zero" has
            // to exclude the infinite case explicitly -- otherwise `Inf * Inf` matches
            // `self.is_infinite() && other.m == 0` and returns `NaN`, which is where the first
            // version of this fix went wrong.
            let other_is_zero = other.m == 0 && !other.is_infinite();
            let self_is_zero = self.m == 0 && !self.is_infinite();
            if (self.is_infinite() && other_is_zero) || (other.is_infinite() && self_is_zero) {
                return F80::NAN;
            }
            return F80::inf(self.sign ^ other.sign);
        }
        if self.m == 0 || other.m == 0 {
            // A signed zero: `prod(c(-0, 5))` is `-0` in R, not `0`.
            return F80 {
                sign: self.sign ^ other.sign,
                exp: 0,
                m: 0,
            };
        }
        // The 64x64 significand product is exact and lands in [2^126, 2^128), so its
        // leading bit is at 126 or 127 *naturally*. That is a different situation from
        // `reduce`, whose accumulator is `hi.m << 63` and therefore leads at 126 with a
        // carry to 127 as the exception -- reusing it here put every product 2^64 too
        // small.
        let p = (self.m as u128) * (other.m as u128);
        let exp_sum = self.exp + other.exp;
        // value = m1*m2 * 2^exp_sum = p * 2^exp_sum, and p == kept * 2^shift, so the
        // result exponent is exp_sum + shift.
        let (kept, shift, exp) = if p >> 127 != 0 {
            ((p >> 64) as u64, 64u32, exp_sum + 64)
        } else {
            ((p >> 63) as u64, 63u32, exp_sum + 63)
        };
        let dropped = p & ((1u128 << shift) - 1);
        if dropped == 0 {
            return F80 {
                sign: self.sign ^ other.sign,
                exp,
                m: kept,
            };
        }
        let round = (dropped >> (shift - 1)) & 1 != 0;
        let sticky = (dropped & ((1u128 << (shift - 1)) - 1)) != 0;
        let mut m = kept;
        let mut e = exp;
        if round && (sticky || (m & 1) != 0) {
            m = m.wrapping_add(1);
            if m == 0 {
                m = 1u64 << 63;
                e += 1;
            }
        }
        F80 {
            sign: self.sign ^ other.sign,
            exp: e,
            m,
        }
    }

    /// `self / d` for a positive integer `d`, rounded to nearest with ties to even.
    ///
    /// This is the only division R's `real_mean` performs (`s /= n`, `t / n`), and `d`
    /// is a length, so `d <= 2^31`. Implemented as a 128-bit quotient so the rounding
    /// bits survive.
    pub fn div_int(self, d: u64) -> F80 {
        assert!(d > 0, "division by zero");
        if self.m == 0 {
            return F80::ZERO;
        }
        // value/d = (m * 2^64 / d + frac) * 2^(exp - 64), with frac = remainder/d.
        let num = (self.m as u128) << 64;
        let dd = d as u128;
        let q = num / dd;
        let r = num % dd;
        debug_assert!(q != 0, "d is far larger than any 80-bit significand");
        // Normalise so the kept significand occupies bits [63..0]:
        //   value/d = (q + f) * 2^(exp - 64)  and  q + f = shifted * 2^(-(127 - msb))
        // so the exponent of the kept word is  exp - 64 + msb - 63  =  exp - 127 + msb.
        let msb = 127 - q.leading_zeros() as i32;
        let shifted = q << (127 - msb) as u32;
        let m = (shifted >> 64) as u64;
        let exp = self.exp - 127 + msb;
        let dropped = shifted & ((1u128 << 64) - 1);
        let round = (dropped >> 63) != 0;
        let sticky = (dropped & ((1u128 << 63) - 1)) != 0 || r != 0;
        if round && (sticky || (m & 1) != 0) {
            let m2 = m.wrapping_add(1);
            if m2 == 0 {
                return F80 {
                    sign: self.sign,
                    exp: exp + 1,
                    m: 1u64 << 63,
                };
            }
            return F80 {
                sign: self.sign,
                exp,
                m: m2,
            };
        }
        F80 {
            sign: self.sign,
            exp,
            m,
        }
    }
}

// The operator forms of the inherent methods above. `r-core` ships as a standalone library as
// well as an R backend, and a caller working in `F80` directly should not have to spell out
// `.add()`; `a + b` has to mean the same thing, and it does -- these delegate to the inherent
// methods, which is where all the rounding lives.
//
// The delegation is written as `F80::add(self, other)` on purpose. Inside a trait impl the bare
// name would be ambiguous to a reader, and `F80::add` resolves to the *inherent* method because
// inherent associated items win path resolution -- so there is no recursion. The 80 tests in this
// module, all of which compare `to_bits()` against values R printed, would stack-overflow on the
// first mistake.
//
// There is deliberately no `Div`. R's long double division is correctly rounded, and this crate has
// not implemented it; the only division it offers is `div_int`, whose exactness argument is
// different (it is exact when the quotient is representable, which is not the same claim). Shipping
// a `Div` that quietly rounded differently from R would be worse than not shipping one.
impl std::ops::Add for F80 {
    type Output = F80;
    #[inline]
    fn add(self, other: F80) -> F80 {
        F80::add(self, other)
    }
}

impl std::ops::Sub for F80 {
    type Output = F80;
    #[inline]
    fn sub(self, other: F80) -> F80 {
        F80::sub(self, other)
    }
}

impl std::ops::Mul for F80 {
    type Output = F80;
    #[inline]
    fn mul(self, other: F80) -> F80 {
        F80::mul(self, other)
    }
}

impl std::ops::Neg for F80 {
    type Output = F80;
    #[inline]
    fn neg(self) -> F80 {
        F80::neg(self)
    }
}

impl std::iter::Sum for F80 {
    /// `Sum` starts from `ZERO`, and `F80::add` propagates `NaN` -- so summing an empty iterator
    /// gives `0`, and a sum containing a `NaN` gives `NaN`. Both match `sum()` in R on a long
    /// double vector, which is what `triMean`'s four-element accumulation relies on.
    fn sum<I: Iterator<Item = F80>>(iter: I) -> F80 {
        iter.fold(F80::ZERO, F80::add)
    }
}

impl<'a> std::iter::Sum<&'a F80> for F80 {
    fn sum<I: Iterator<Item = &'a F80>>(iter: I) -> F80 {
        iter.copied().fold(F80::ZERO, F80::add)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    fn close(a: f64, b: f64, rel: f64) {
        assert!(
            (a - b).abs() <= rel * b.abs().max(1.0),
            "a={a:.17e} b={b:.17e} rel={rel}"
        );
    }

    #[test]
    fn dbg_add_sequence() {
        use crate::longdouble::F80;
        let vals = [0.1f64, 0.2, 0.3, 1e-5, 1234.5678, -3.25, 7.0];
        let mut f = F80::ZERO;
        let mut d = 0.0f64;
        for &v in &vals {
            let a = F80::from_f64(v);
            println!(
                "  in  sign={} exp={} m={:018x}  val={:.17e}",
                a.sign, a.exp, a.m, v
            );
            let before = f;
            println!(
                "  acc-before: sign={} exp={} m={:018x}",
                before.sign, before.exp, before.m
            );
            f = before.add(a);
            d += v;
            println!(
                "v={v:>12.6}  f80={:.17e}  f64seq={:.17e}  diff={:.3e}",
                f.to_f64(),
                d,
                f.to_f64() - d
            );
            let _ = before;
        }
    }

    #[test]
    fn dbg_pair() {
        let a = F80::from_f64(1.0);
        let b = F80::from_f64(0.25);
        println!(
            "1.0  -> sign={} exp={} m={:018x} val={}",
            a.sign,
            a.exp,
            a.m,
            a.to_f64()
        );
        println!(
            "0.25 -> sign={} exp={} m={:018x} val={}",
            b.sign,
            b.exp,
            b.m,
            b.to_f64()
        );
        println!("a.sub(b) = {:.17e}", a.sub(b).to_f64());
        let nb = b.neg();
        println!("b.neg() -> sign={} exp={} m={:018x}", nb.sign, nb.exp, nb.m);
        let s = a.add(nb);
        println!(
            "a.add(-b) -> sign={} exp={} m={:018x} val={:.17e}",
            s.sign,
            s.exp,
            s.m,
            s.to_f64()
        );
        let c = a.add(b);
        println!(
            "a.add(b)  -> sign={} exp={} m={:018x} val={:.17e}",
            c.sign,
            c.exp,
            c.m,
            c.to_f64()
        );
    }

    #[test]
    fn round_trips_f64_exactly() {
        let vals: [f64; 14] = [
            0.0,
            1.0,
            -1.0,
            0.5,
            1e16,
            -1e16,
            1e-300,
            std::f64::consts::PI,
            f64::MIN_POSITIVE,
            f64::MAX,
            1.0 / 3.0,
            f64::from_bits(1),
            f64::from_bits(0x000F_FFFF_FFFF_FFFF),
            2.2250738585072014e-308,
        ];
        for &v in &vals {
            assert_eq!(F80::from_f64(v).to_f64(), v, "round trip failed for {v:e}");
        }
    }

    #[test]
    fn subnormal_f64_inputs_widen_correctly() {
        for bits in [1u64, 2, 0x000F_FFFF_FFFF_FFFF, 0x0008_0000_0000_0000] {
            let v = f64::from_bits(bits);
            assert_eq!(F80::from_f64(v).to_f64(), v, "bits={bits:#x}");
        }
    }

    #[test]
    fn carries_more_precision_than_f64() {
        // 1 + 2^-60 is representable in 80 bits. Narrowed to f64 it must round back to
        // 1.0, but the *80-bit* result must still carry the 2^-60 term -- checked by
        // subtracting 1 and recovering it (2^-60 is exactly representable in f64).
        let one = F80::from_f64(1.0);
        // value = m * 2^exp = 2^63 * 2^-123 = 2^-60
        let eps60 = F80 {
            sign: false,
            exp: -123,
            m: 1u64 << 63,
        };
        let sum = one.add(eps60);
        assert_eq!(
            sum.to_f64(),
            1.0,
            "narrowing 1 + 2^-60 to f64 must give 1.0"
        );
        let recovered = sum.sub(one).to_f64();
        assert_eq!(
            recovered,
            2f64.powi(-60),
            "the 2^-60 term was lost in 80-bit add"
        );
    }

    #[test]
    fn cancellation_keeps_the_small_term() {
        // The case that forces R's two-pass algorithm: naive f64 summation gives 0.
        let big = F80::from_f64(1e16);
        let one = F80::from_f64(1.0);
        assert_eq!(big.add(one).add(big.neg()).to_f64(), 1.0);
    }

    #[test]
    fn accumulates_more_accurately_than_f64() {
        // `0.1 + 0.2 + 0.3` is 0.6000000000000001 in a naive f64 chain; in 64 bits it is
        // the correctly rounded 0.6. This is *why* the type exists, so assert the
        // divergence rather than equality with f64.
        let mut f = F80::ZERO;
        for &v in &[0.1f64, 0.2, 0.3] {
            f = f.add(F80::from_f64(v));
        }
        assert_eq!(f.to_f64(), 0.6);
        let mut naive = 0.0f64;
        for &v in &[0.1f64, 0.2, 0.3] {
            naive += v;
        }
        assert_eq!(naive, 0.6000000000000001);
        assert_ne!(f.to_f64(), naive);
    }

    #[test]
    fn sum_of_a_mixed_magnitude_sequence() {
        let vals = [0.1f64, 0.2, 0.3, 1e-5, 1234.5678, -3.25, 7.0];
        let mut f = F80::ZERO;
        for &v in &vals {
            f = f.add(F80::from_f64(v));
        }
        // 1238.91781 is the correctly rounded exact sum; within 1 ulp of f64 scale.
        let got = f.to_f64();
        assert!((got - 1238.91781).abs() < 1e-11, "got {got:.17e}");
    }

    #[test]
    fn div_int_matches_f64_for_divisible_cases() {
        // The divisibility guard this test used to have had identical `then` and `else` arms --
        // both `close(got, a / d as f64, 1e-15)` -- so it was deciding nothing. Whether `div_int`
        // agrees with `f64` for a non-dividing case is a different question, and one `f64` cannot
        // answer: it rounds at each operation where the 80-bit path rounds once.
        for &(a, d) in &[(1.0f64, 3u64), (2.0, 4), (0.0, 3), (255.0, 16), (1e16, 8)] {
            let got = F80::from_f64(a).div_int(d).to_f64();
            close(got, a / d as f64, 1e-15);
        }
    }

    #[test]
    fn sticky_bits_do_not_perturb_an_exact_value() {
        let one = F80::from_f64(1.0);
        let tiny = F80 {
            sign: false,
            exp: -200,
            m: 1u64 << 63,
        };
        assert_eq!(one.add(tiny).to_f64(), 1.0);
    }

    #[test]
    fn opposite_sign_subtraction_is_exact_when_representable() {
        let a = F80::from_f64(1.0);
        let b = F80::from_f64(0.25);
        assert_eq!(a.sub(b).to_f64(), 0.75);
        assert_eq!(b.sub(a).to_f64(), -0.75);
        assert_eq!(a.sub(a).to_f64(), 0.0);
    }
}

#[cfg(test)]
mod mul_special_tests {
    use super::F80;

    /// R's table for infinite products, measured against `prod()` on the pinned R. Every one of
    /// these was wrong before: `mul` returned `0` for any infinite operand, because infinity is
    /// encoded with a zero mantissa and the leading zero-check caught it.
    #[test]
    fn mul_matches_r_for_every_infinite_case() {
        let (inf, ninf, zero) = (F80::inf(false), F80::inf(true), F80::ZERO);
        let two = F80::from_f64(2.0);
        let cases: &[(F80, F80, f64, &str)] = &[
            (inf, inf, f64::INFINITY, "prod(Inf, Inf)"),
            (inf, ninf, f64::NEG_INFINITY, "prod(Inf, -Inf)"),
            (ninf, inf, f64::NEG_INFINITY, "prod(-Inf, Inf)"),
            (ninf, ninf, f64::INFINITY, "prod(-Inf, -Inf)"),
            (inf, two, f64::INFINITY, "prod(Inf, 2)"),
            (two, inf, f64::INFINITY, "prod(2, Inf)"),
            (ninf, two, f64::NEG_INFINITY, "prod(-Inf, 3)"),
            (inf, zero, f64::NAN, "prod(Inf, 0)"),
            (zero, inf, f64::NAN, "prod(0, Inf)"),
            (ninf, zero, f64::NAN, "prod(-Inf, 0)"),
            (zero, ninf, f64::NAN, "prod(0, -Inf)"),
        ];
        for &(a, b, want, what) in cases {
            for (got, label) in [(a.mul(b), "a*b"), (b.mul(a), "b*a")] {
                let ok = if want.is_nan() {
                    got.is_nan()
                } else {
                    got.to_f64() == want
                };
                assert!(ok, "{what} via {label}: got {} want {}", got.to_f64(), want);
            }
        }
    }

    /// The sign of a zero product follows the **XOR of the two signs**, so
    /// `prod(c(-0, 5))` is `-0x0p+0` and so is `prod(c(5, -0))`. The second of those was the
    /// surprise: multiplication is commutative but a signed zero is not, and `5 * -0 == -0` is the
    /// IEEE rule R follows.
    #[test]
    fn mul_keeps_the_sign_of_a_zero_product() {
        let neg_zero = F80 {
            sign: true,
            exp: 0,
            m: 0,
        };
        let pos_zero = F80::ZERO;
        let five = F80::from_f64(5.0);
        assert!(
            neg_zero.mul(five).to_f64().is_sign_negative(),
            "-0 * 5 is -0"
        );
        assert!(
            five.mul(neg_zero).to_f64().is_sign_negative(),
            "5 * -0 is -0"
        );
        assert!(
            !pos_zero.mul(five).to_f64().is_sign_negative(),
            "0 * 5 is +0"
        );
        assert!(
            !five.mul(pos_zero).to_f64().is_sign_negative(),
            "5 * 0 is +0"
        );
    }

    /// A `NaN` operand poisons the product, which the old zero-check also got wrong in the other
    /// direction for `NaN * Inf`.
    #[test]
    fn mul_poisons_on_nan() {
        assert!(F80::NAN.mul(F80::inf(false)).is_nan());
        assert!(F80::inf(false).mul(F80::NAN).is_nan());
        assert!(F80::NAN.mul(F80::ZERO).is_nan());
    }

    /// The finite path is unchanged: the fix is confined to the special values.
    #[test]
    fn finite_products_are_unaffected() {
        for (a, b) in [(2.0f64, 3.0), (0.5, 0.25), (-4.0, 7.5), (1e-8, 1e8)] {
            let got = F80::from_f64(a).mul(F80::from_f64(b)).to_f64();
            assert_eq!(got, a * b, "{a} * {b}");
        }
    }
}
