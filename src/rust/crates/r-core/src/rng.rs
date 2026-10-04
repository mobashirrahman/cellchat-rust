//! Bit-exact port of R's Mersenne-Twister RNG, as used by
//! `computeCommunProb`'s `set.seed(seed.use); replicate(nboot, sample.int(nC, nC))`.
//!
//! Parity contract: **Exact** (`docs/SEMANTICS.md` R1). This module exists because the
//! bootstrap permutations determine every p-value, so a different permutation stream
//! would yield a statistically equivalent but numerically *different* network — not
//! acceptable under the locked bit-identical target (`PLAN.md` §14.2).
//!
//! Sources (R 4.3 branch, matching the pinned upstream's R runtime):
//!
//! * `src/main/RNG.c` — `RNG_Init`, `FixupSeeds`, `MT_sgenrand`, `MT_genrand`,
//!   `unif_rand`, `fixup`, `rbits`, `R_unif_index`
//! * `src/main/random.c` — `do_sample`, the uniform without-replacement branch
//!
//! ## Three subtleties that a textbook MT19937 gets wrong
//!
//! 1. **The seed is not fed to `MT_sgenrand`.** `RNG_Init` scrambles with 50 steps of
//!    `seed = 69069*seed + 1` and then fills *625* words (not 624) with the same LCG.
//!    `i_seed[0]` holds the index `mti`, and `FixupSeeds(MERSENNE_TWISTER, initial =
//!    TRUE)` then **overwrites it with the constant 624**. The state array is
//!    `i_seed[1..=624]`.
//!
//!    This is the most error-prone part: a port that keeps the LCG's `i_seed[0]` reads
//!    `mt[negative]` and silently produces a degenerate stream that still *looks* like
//!    a permutation.
//! 2. **`Int32` is `unsigned int`, not `int`.** `RNG.c` says so in a comment
//!    (`/* typedef unsigned int Int32; in Random.h */`), so the right shifts are
//!    *logical* and R's stream is the standard MT19937 — not a signed-shift variant.
//!    Getting this backwards is easy to do by accident and produces a perfectly
//!    uniform but completely wrong stream.
//!
//!    Verified empirically against R 4.3.3: all 2e6 draws of `set.seed(1); runif(2e6)`
//!    lie exactly on the `k * 2^-32` grid, and none equals `fixup`'s fallback value —
//!    both of which are only possible with logical shifts and non-negative words.
//! 3. **Two different 2^-32 constants.** `unif_rand` scales by
//!    `2.3283064365386963e-10` (2^-32) whereas `fixup`'s fallback uses
//!    `i2_32m1 = 2.328306437080797e-10` (1/(2^32-1)). Substituting one for the other
//!    shifts every value.

/// `1/(2^32 - 1)`, R's `i2_32m1`. Only used by [`fixup`].
const I2_32M1: f64 = 2.328_306_437_080_797e-10;
/// `2^-32`, R's MT19937 reals scaling.
const TWO_POW_M32: f64 = 2.328_306_436_538_696_3e-10;

const N: usize = 624;
const M: usize = 397;
const MATRIX_A: u32 = 0x9908_b0df;
const UPPER_MASK: u32 = 0x8000_0000;
const LOWER_MASK: u32 = 0x7fff_ffff;
const TEMPERING_MASK_B: u32 = 0x9d2c_5680;
const TEMPERING_MASK_C: u32 = 0xefc6_0000;

/// R's Mersenne-Twister state, layout-compatible with R's `.Random.seed` entry for a
/// 625-word Mersenne-Twister generator.
///
/// `mt` is `i_seed[1..=624]`; `mti` is `i_seed[0]`. Both are `u32` because R's `Int32`
/// is `unsigned int` (module docs, subtlety 2).
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct MersenneTwister {
    mt: [u32; N],
    mti: u32,
}

impl MersenneTwister {
    /// Equivalent of R's `set.seed(seed)` for `RNGkind() == "Mersenne-Twister"`.
    ///
    /// R's `do_setseed` truncates the requested seed to `Int32` and passes it to
    /// `RNG_Init(MERSENNE_TWISTER, seed)`. Negative seeds arrive as their two's
    /// complement bit pattern, so the `i32 -> u32` cast is the faithful conversion.
    pub fn new(seed: i32) -> Self {
        // `RNG_Init`: "Initial scrambling" — 50 LCG steps BEFORE the fill loop.
        let mut s = seed as u32;
        for _ in 0..50 {
            s = s.wrapping_mul(69069).wrapping_add(1);
        }
        // Fill `i_seed[0..=624]`: 625 LCG steps.
        let mut i_seed = [0u32; N + 1];
        for slot in i_seed.iter_mut() {
            s = s.wrapping_mul(69069).wrapping_add(1);
            *slot = s;
        }
        // `FixupSeeds(MERSENNE_TWISTER, initial = TRUE)`: `i_seed[0]` (= `mti`) := 624.
        // The "all zeroes" guard cannot fire: the LCG orbit is nonzero, so 625
        // consecutive zero words are unreachable.
        i_seed[0] = N as u32;
        let mut mt = [0u32; N];
        mt.copy_from_slice(&i_seed[1..]);
        Self { mt, mti: N as u32 }
    }

    /// Port of R's `MT_sgenrand`. Reached by [`MersenneTwister::uninitialized`].
    pub fn from_mt_sgenrand(seed: i32) -> Self {
        let mut mt = [0u32; N];
        let mut s = seed as u32;
        for slot in mt.iter_mut() {
            *slot = s & 0xffff_0000;
            s = s.wrapping_mul(69069).wrapping_add(1);
            *slot |= (s & 0xffff_0000) >> 16;
            s = s.wrapping_mul(69069).wrapping_add(1);
        }
        Self { mt, mti: N as u32 }
    }

    /// The state R has when `unif_rand` is called before any `set.seed`:
    /// `MT_genrand` observes `mti == N+1` and self-initialises from 4357.
    pub fn uninitialized() -> Self {
        Self {
            mt: [0; N],
            mti: N as u32 + 1,
        }
    }

    /// Port of `FixupSeeds` for a state reloaded from a possibly user-corrupted
    /// `.Random.seed`.
    pub fn from_raw(mt: [u32; N], mut mti: u32) -> Self {
        if mti == 0 {
            mti = N as u32;
        }
        if mt.iter().all(|&w| w == 0) {
            return Self::from_mt_sgenrand(4357);
        }
        Self { mt, mti }
    }

    /// The full state vector as R stores it in `.Random.seed` (`i_seed[0..=624]`).
    pub fn to_raw(&self) -> ([u32; N], u32) {
        (self.mt, self.mti)
    }

    /// Port of R's `MT_genrand`, returning the raw scaled double *before* [`fixup`].
    fn genrand(&mut self) -> f64 {
        if self.mti >= N as u32 {
            if self.mti == N as u32 + 1 {
                // `MT_genrand`: "a default initial seed is used" (4357).
                self.mt = Self::from_mt_sgenrand(4357).mt;
            }
            for kk in 0..(N - M) {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] = self.mt[kk + M] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            }
            for kk in (N - M)..(N - 1) {
                let y = (self.mt[kk] & UPPER_MASK) | (self.mt[kk + 1] & LOWER_MASK);
                self.mt[kk] =
                    self.mt[kk + M - N] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            }
            let y = (self.mt[N - 1] & UPPER_MASK) | (self.mt[0] & LOWER_MASK);
            self.mt[N - 1] = self.mt[M - 1] ^ (y >> 1) ^ if y & 1 != 0 { MATRIX_A } else { 0 };
            self.mti = 0;
        }
        let mut y = self.mt[self.mti as usize];
        self.mti += 1;
        y ^= y >> 11;
        y ^= (y << 7) & TEMPERING_MASK_B;
        y ^= (y << 15) & TEMPERING_MASK_C;
        y ^= y >> 18;
        y as f64 * TWO_POW_M32
    }

    /// R's `unif_rand()`: `fixup(MT_genrand())`, in `[0, 1)`.
    pub fn unif_rand(&mut self) -> f64 {
        fixup(self.genrand())
    }

    /// R's `rbits(bits)`: `ceil(bits/16) + 1` draws of 16 bits each, masked to `bits`.
    ///
    /// R's loop is `for (n = 0; n <= bits; n += 16)`, so `bits <= 15` consumes one draw,
    /// `bits == 16` consumes two, and so on. `bits` is at most 31 on every path
    /// CellChat exercises. The number of draws is load-bearing: it decides how much of
    /// the stream a single `R_unif_index` call consumes.
    fn rbits(&mut self, bits: u32) -> f64 {
        let mut v: u64 = 0;
        let mut n = 0u32;
        while n <= bits {
            let v1 = (self.unif_rand() * 65536.0).floor() as i64;
            v = 65_536u64.wrapping_mul(v).wrapping_add(v1 as u64);
            n += 16;
        }
        (v & ((1u64 << bits) - 1)) as f64
    }

    /// R's `R_unif_index(dn)` under `sample.kind = "Rejection"` (the R >= 3.6 default):
    /// rejection sampling from integers below the next larger power of two.
    pub fn r_unif_index(&mut self, dn: f64) -> f64 {
        if dn <= 0.0 {
            return 0.0;
        }
        let bits = bits_needed(dn);
        loop {
            let dv = self.rbits(bits);
            if dn > dv {
                return dv;
            }
        }
    }

    /// Port of R's `sample.int(n, k, replace = FALSE)` for the shape CellChat uses:
    /// `n <= i32::MAX`, `k >= 2`, uniform, without replacement.
    ///
    /// R's `do_sample` uniform branch:
    /// ```c
    /// int *x = R_alloc(n, sizeof(int));
    /// for (int i = 0; i < n; i++) x[i] = i;
    /// for (int i = 0; i < k; i++) {
    ///     int j = (int)(R_unif_index(n));
    ///     iy[i] = x[j] + 1;
    ///     x[j] = x[--n];
    /// }
    /// ```
    /// This is selection sampling (a partial Fisher-Yates shuffle): the pool `x[0..n)`
    /// shrinks by one per step. Returned values are **1-based**, matching R.
    pub fn sample_int_no_replace(&mut self, n: usize, k: usize) -> Vec<i32> {
        assert!(k <= n, "cannot take a sample larger than the population");
        let mut pool: Vec<i32> = (0..n as i32).collect();
        let mut remaining = n as i32;
        let mut out = Vec::with_capacity(k);
        for _ in 0..k {
            // NB: R's loop variable `n` is decremented by `x[j] = x[--n]`, and it is
            // that *shrinking* `n` that is handed to `R_unif_index` on the next
            // iteration. Sampling against the original population size yields
            // duplicates — caught immediately by `sample_int_is_a_permutation`.
            let j = self.r_unif_index(remaining as f64) as i32;
            debug_assert!(j >= 0 && j < remaining);
            out.push(pool[j as usize] + 1);
            remaining -= 1;
            pool[j as usize] = pool[remaining as usize];
        }
        out
    }

    /// `sample.int(n, n)` — the full permutation `computeCommunProb` builds.
    pub fn sample_int_permutation(&mut self, n: usize) -> Vec<i32> {
        self.sample_int_no_replace(n, n)
    }

    /// Port of R's `sample.int(n, k, replace = TRUE)`: `iy[i] = R_unif_index(n) + 1`.
    ///
    /// Not on the `computeCommunProb` path, but implemented so this module is a complete
    /// stand-in for the RNG surface CellChat can reach, and covered by the golden file.
    pub fn sample_int_with_replace(&mut self, n: usize, k: usize) -> Vec<i32> {
        let dn = n as f64;
        (0..k).map(|_| self.r_unif_index(dn) as i32 + 1).collect()
    }
}

/// R's `fixup`: guarantees a usable uniform in `[0, 1)`.
///
/// Dead code for Mersenne-Twister, whose `genrand` output is in `(0, 1)` by
/// construction. Retained because R applies it unconditionally, and because the
/// `MT_sgenrand` fallback path can legitimately produce a zero word.
#[inline]
pub fn fixup(x: f64) -> f64 {
    if x <= 0.0 {
        return 0.5 * I2_32M1;
    }
    if (1.0 - x) <= 0.0 {
        return 1.0 - 0.5 * I2_32M1;
    }
    x
}

/// `ceil(log2(dn))` for `dn > 0`, computed exactly.
///
/// This is the smallest `b` with `2^b >= dn`, which equals `bit_length(ceil(dn) - 1)`.
/// Using `bit_length(dn)` is off by one at every power of two; using `floor(dn)` is off
/// by one for every non-integer. Both mistakes shift the whole downstream stream
/// because they change how many `unif_rand` draws a `rbits` call consumes.
#[inline]
fn bits_needed(dn: f64) -> u32 {
    debug_assert!(dn > 0.0 && dn.is_finite());
    let m = (dn.ceil() as u64).wrapping_sub(1);
    u64::BITS - m.leading_zeros()
}

#[cfg(test)]
mod tests {
    use super::*;

    /// The exhaustive parity check lives in `tests/rng_parity.rs`, which reads the
    /// golden file straight out of R. These two smoke tests assert the documented
    /// R 4.3.3 values for `set.seed(1)` so a broken build fails fast here.
    #[test]
    fn first_unif_rand_matches_r() {
        let mut rng = MersenneTwister::new(1);
        assert_eq!(rng.unif_rand(), 0.265508663142099977);
    }

    #[test]
    fn first_permutation_matches_r() {
        let mut rng = MersenneTwister::new(1);
        assert_eq!(
            rng.sample_int_permutation(10),
            vec![9, 4, 7, 1, 2, 5, 3, 10, 6, 8]
        );
    }

    #[test]
    fn sample_int_is_a_permutation() {
        for n in [1usize, 2, 3, 7, 64, 1000] {
            let mut rng = MersenneTwister::new(42);
            let s = rng.sample_int_permutation(n);
            let mut sorted = s.clone();
            sorted.sort_unstable();
            assert_eq!(sorted, (1..=n as i32).collect::<Vec<_>>(), "n = {n}");
        }
    }

    #[test]
    fn partial_sample_is_without_replacement() {
        let mut rng = MersenneTwister::new(7);
        let s = rng.sample_int_no_replace(100, 40);
        let mut sorted = s.clone();
        sorted.sort_unstable();
        sorted.dedup();
        assert_eq!(sorted.len(), 40);
        assert!(s.iter().all(|&v| (1..=100).contains(&v)));
    }

    #[test]
    fn unif_rand_stays_in_unit_interval() {
        let mut rng = MersenneTwister::new(3);
        for _ in 0..200_000 {
            let u = rng.unif_rand();
            assert!((0.0..1.0).contains(&u), "out of range: {u}");
        }
    }

    #[test]
    fn mt_outputs_are_always_strictly_positive() {
        // Documents why `fixup`'s fallback is dead code for Mersenne-Twister: `Int32`
        // is unsigned, so `y as f64 * 2^-32 >= 0` always. If this fails, the port has
        // slipped into signed-shift semantics and the stream is wrong.
        let mut rng = MersenneTwister::new(1);
        for _ in 0..200_000 {
            assert!(rng.genrand() > 0.0);
        }
    }

    #[test]
    fn fixup_clamps_the_full_double_range() {
        assert_eq!(fixup(-0.5), 0.5 * I2_32M1);
        assert_eq!(fixup(0.0), 0.5 * I2_32M1);
        assert_eq!(fixup(1.0), 1.0 - 0.5 * I2_32M1);
        assert_eq!(fixup(0.25), 0.25);
    }

    #[test]
    fn rbits_bit_width_is_exact() {
        assert_eq!(bits_needed(1.0), 0);
        assert_eq!(bits_needed(1.0 + f64::EPSILON), 1);
        assert_eq!(bits_needed(2.0), 1);
        assert_eq!(bits_needed(3.0), 2);
        assert_eq!(bits_needed(4.0), 2);
        assert_eq!(bits_needed(4.5), 3);
        assert_eq!(bits_needed(5.0), 3);
        assert_eq!(bits_needed(21557.0), 15);
        assert_eq!(bits_needed(65536.0), 16);
        assert_eq!(bits_needed(65537.0), 17);
        assert_eq!(bits_needed(i32::MAX as f64), 31);
    }

    /// `.Random.seed` round-trip: R stores exactly these 625 words, so a save/restore
    /// cycle must not perturb the stream.
    #[test]
    fn raw_state_round_trips() {
        let mut a = MersenneTwister::new(11);
        for _ in 0..37 {
            a.unif_rand();
        }
        let (mt, mti) = a.to_raw();
        let mut b = MersenneTwister::from_raw(mt, mti);
        for _ in 0..1000 {
            assert_eq!(a.unif_rand(), b.unif_rand());
        }
    }

    /// The seed LCG must wrap on 32 bits. A non-wrapping implementation drifts after
    /// the first few steps, so this pins the behaviour.
    #[test]
    fn seed_lcg_wraps_at_32_bits() {
        let first = MersenneTwister::new(i32::MAX).unif_rand();
        assert!(first.is_finite() && (0.0..1.0).contains(&first));
        assert_eq!(MersenneTwister::new(i32::MAX).unif_rand(), first);
        assert_ne!(
            MersenneTwister::new(i32::MAX).unif_rand(),
            MersenneTwister::new(i32::MIN).unif_rand()
        );
    }
}
