//! Faithful ports of the R 4.3.3 C routines that `identifyOverExpressedGenes` reaches
//! through `stats::wilcox.test`.
//!
//! These are **not** "a good implementation of the normal CDF". They are R's, transcribed
//! from `src/nmath/pnorm.c`, `src/nmath/lgamma.c`, `src/nmath/lbeta.c` and
//! `src/nmath/choose.c` at tag `R-4.3.3`, coefficient arrays and all, because
//! `identifyOverExpressedGenes` puts a `wilcox.test` p-value straight into a
//! Bonferroni-corrected threshold and a marker list, and a last-ulp difference in
//! `pnorm` becomes a different set of significant genes.
//!
//! The sources are GPL-2-or-later (R's licence), and this crate is GPL-3, so the derivation
//! is compliant; the upstream commit, the R version and the file names are recorded in each
//! item's doc comment so the provenance is checkable.
//!
//! ## What is *not* here, and why it does not matter
//!
//! `pnorm`, `expm1`, `log`, `log1p`, `sqrt`, `trunc` and `ldexp` come from the platform's
//! libm, which R also calls, so the Rust `f64::` versions are the same code. `gammafn` is
//! used by `lbeta` only when *both* arguments are below 10, and `lbeta` is reached from
//! `choose` only for `k >= 30`; [`lgamma_small`] therefore uses `libm::lgamma` there, which
//! is the same call R makes (`Rmath.h0.in` maps `lgamma` to libm).
//!
//! [`lgammacor`] — R's rational correction to Stirling's series for `x > 10` — *is* ported,
//! because `lbeta`'s `p >= 10` branch uses it and that branch is the one the exact
//! Wilcoxon p-value takes.

/// `sqrt(32)`, R's `M_SQRT_32`. The switch point between the central and tail regions of
/// `pnorm`.
const M_SQRT_32: f64 = 5.656854249492380195206754896838;
/// `1/sqrt(2*pi)`, R's `M_1_SQRT_2PI`.
const M_1_SQRT_2PI: f64 = 0.398942280401432677939946059934;
/// `log(sqrt(2*pi))`, R's `M_LN_SQRT_2PI`.
const M_LN_SQRT_2PI: f64 = 0.918938533204672741780329736406;
/// `log(sqrt(pi/2))`, R's `M_LN_SQRT_PId2`.
#[allow(non_upper_case_globals)]
const M_LN_SQRT_PId2: f64 = 0.225791352644727432363097614947;

/// R 4.3.3 `pnorm` reference values, from `format(x, digits=17)`.
#[cfg(test)]
const PNORM_REF_: [(f64, f64); 30] = [
    (0.0, 0.5),
    (1.0, 0.84134474606854293),
    (-1.0, 0.15865525393145705),
    (2.0, 0.97724986805182079),
    (-2.0, 0.022750131948179212),
    (0.5, 0.69146246127401301),
    (3.0, 0.9986501019683699),
    (-3.0, 0.0013498980316300946),
    (5.0, 0.99999971334842808),
    (-5.0, 2.8665157187919391e-07),
    (0.67448975, 0.74999999993768973),
    (5.6568542494923797, 0.99999999229137104),
    (8.0, 0.99999999999999933),
    (10.0, 1.0),
    (37.0, 1.0),
    (-37.0, 5.7255712225245771e-300),
    (1e-04, 0.50003989422797368),
    (0.67448975010000001, 0.74999999996946753),
    (1.5, 0.93319279873114191),
    (-1.5, 0.066807201268858071),
    (4.5, 0.99999660232687526),
    (-4.5, 3.3976731247300598e-06),
    (6.0, 0.9999999990134123),
    (-6.0, 9.8658764503769809e-10),
    (12.0, 1.0),
    (-12.0, 1.776482112077679e-33),
    (20.0, 1.0),
    (-20.0, 2.7536241186062337e-89),
    (0.29999999999999999, 0.61791142218895267),
    (-0.29999999999999999, 0.38208857781104733),
];

/// R 4.3.3 `pnorm(x, lower.tail = FALSE)` reference values. Compared bit for bit as
/// well: R evaluates the upper tail through the `swap_tail` branch rather than as
/// `1 - lower`, and at `x = 8` the lower tail is one ulp below 1, so subtracting
/// throws away the entire answer.
/// R 4.3.3 `pnorm(x, lower.tail = FALSE)` reference values. Compared bit for bit as
/// well: R evaluates the upper tail through the `swap_tail` branch rather than as
/// `1 - lower`, and at `x = 8` the lower tail is one ulp below 1, so subtracting throws
/// away the entire answer.
#[cfg(test)]
const PNORM_UP_REF_: [(f64, f64); 30] = [
    (0.0, 0.5),
    (1.0, 0.15865525393145705),
    (-1.0, 0.84134474606854293),
    (2.0, 0.022750131948179212),
    (-2.0, 0.97724986805182079),
    (0.5, 0.30853753872598694),
    (3.0, 0.0013498980316300946),
    (-3.0, 0.9986501019683699),
    (5.0, 2.8665157187919391e-07),
    (-5.0, 0.99999971334842808),
    (0.67448975, 0.25000000006231027),
    (5.6568542494923797, 7.7086289501400348e-09),
    (8.0, 6.2209605742717849e-16),
    (10.0, 7.6198530241605269e-24),
    (37.0, 5.7255712225245771e-300),
    (-37.0, 1.0),
    (1e-04, 0.49996010577202632),
    (0.67448975010000001, 0.25000000003053252),
    (1.5, 0.066807201268858071),
    (-1.5, 0.93319279873114191),
    (4.5, 3.3976731247300598e-06),
    (-4.5, 0.99999660232687526),
    (6.0, 9.8658764503769809e-10),
    (-6.0, 0.9999999990134123),
    (12.0, 1.776482112077679e-33),
    (-12.0, 1.0),
    (20.0, 2.7536241186062337e-89),
    (-20.0, 1.0),
    (0.29999999999999999, 0.38208857781104733),
    (-0.29999999999999999, 0.61791142218895267),
];

/// R 4.3.3 `choose` reference values.
#[cfg(test)]
const CHOOSE_REF_: [(f64, f64, f64); 12] = [
    (5.0, 2.0, 10.0),
    (10.0, 3.0, 120.0),
    (49.0, 24.0, 63205303218876.0),
    (98.0, 49.0, 2.5477612258980968e+28),
    (30.0, 15.0, 155117520.0),
    (60.0, 30.0, 118264581564861152.0),
    (2.0, 1.0, 2.0),
    (1.0, 0.0, 1.0),
    (45.0, 22.0, 4116715363800.0),
    (51.0, 25.0, 247959266474052.0),
    (52.0, 26.0, 495918532948104.0),
    (60.0, 31.0, 114449595062769136.0),
];

/// R 4.3.3 `lgamma` reference values, spanning the Chebyshev range, the `xbig` branch
/// point (`2^26.5`) and the `x > 1e17` asymptotic branch.
#[cfg(test)]
const LGAMMA_REF_: [(f64, f64); 12] = [
    (11.0, 15.104412573075519),
    (20.0, 39.339884187199495),
    (100.0, 359.1342053695754),
    (1000.0, 5905.2204232091808),
    (1e+05, 1051287.7089736569),
    (1e+06, 12815504.569147611),
    (1e+08, 1742068066.1038349),
    (1e+17, 3814394658089878016.0),
    (1e+18, 40446531673892823040.0),
    (94906265.624251559, 1648370002.6359525),
    (94906265.700000003, 1648370004.0273299),
    (1e+20, 4505170185988091674624.0),
];

/// R's `R_forceint`: `round()` in C, i.e. half-away-from-zero, unlike Rust's `round()`
/// (half-to-even).
#[inline]
fn r_forceint(x: f64) -> f64 {
    // C99 `round` rounds half away from zero.
    if x >= 0.0 {
        (x + 0.5).floor()
    } else {
        (x - 0.5).ceil()
    }
}

/// R's `R_nonint`: `(fabs(x - R_forceint(x)) > 1e-7 * fmax2(1., fabs(x)))`.
#[inline]
fn r_nonint(x: f64) -> bool {
    (x - r_forceint(x)).abs() > 1e-7 * f64::max(1.0, x.abs())
}

/// `pnorm5(x, mu, sigma, lower_tail, log_p)` from `src/nmath/pnorm.c` (R 4.3.3), with
/// `mu = 0` and `sigma = 1` folded in — which is the only shape `wilcox.test` uses.
///
/// Transcribed from `pnorm_both`, which is W. J. Cody's rational Chebyshev approximation
/// (ACM TOMS 19:22-32, ALGORITHM 715). The coefficient arrays are Cody's and must not be
/// "improved": a more accurate `pnorm` would *stop* matching R.
pub fn pnorm(x: f64, lower_tail: bool) -> f64 {
    if x.is_nan() {
        return x;
    }
    // sigma = 1 > 0, so the degenerate branch is unreachable.
    let (cum, ccum) = pnorm_both(x, true, !lower_tail);
    if lower_tail {
        cum
    } else {
        ccum
    }
}

/// `pnorm_both(x, &cum, &ccum, i_tail, log_p = FALSE)`.
///
/// `want_lower` / `want_upper` select which tails are computed, mirroring `i_tail`. R
/// computes only the one it needs, which matters not for the value but for the code
/// structure, so the branches are kept.
fn pnorm_both(x: f64, want_lower: bool, want_upper: bool) -> (f64, f64) {
    const A: [f64; 5] = [
        2.2352520354606839287,
        161.02823106855587881,
        1067.6894854603709582,
        18154.981253343561249,
        0.065682337918207449113,
    ];
    const B: [f64; 4] = [
        47.20258190468824187,
        976.09855173777669322,
        10260.932208618978205,
        45507.789335026729956,
    ];
    const C: [f64; 9] = [
        0.39894151208813466764,
        8.8831497943883759412,
        93.506656132177855979,
        597.27027639480026226,
        2494.5375852903726711,
        6848.1904505362823326,
        11602.651437647350124,
        9842.7148383839780218,
        1.0765576773720192317e-8,
    ];
    const D: [f64; 8] = [
        22.266688044328115691,
        235.38790178262499861,
        1519.377599407554805,
        6485.558298266760755,
        18615.571640885098091,
        34900.952721145977266,
        38912.003286093271411,
        19685.429676859990727,
    ];
    const P: [f64; 6] = [
        0.21589853405795699,
        0.1274011611602473639,
        0.022235277870649807,
        0.001421619193227893466,
        2.9112874951168792e-5,
        0.02307344176494017303,
    ];
    const Q: [f64; 5] = [
        1.28426009614491121,
        0.468238212480865118,
        0.0659881378689285515,
        0.00378239633202758244,
        7.29751555083966205e-5,
    ];

    if x.is_nan() {
        return (x, x);
    }
    let eps = f64::EPSILON * 0.5;
    let mut cum = 0.0;
    let mut ccum = 0.0;
    let y = x.abs();

    if y <= 0.67448975 {
        // qnorm(3/4)
        let (xnum, xden) = if y > eps {
            let xsq = x * x;
            let mut xnum = A[4] * xsq;
            let mut xden = xsq;
            for i in 0..3 {
                xnum = (xnum + A[i]) * xsq;
                xden = (xden + B[i]) * xsq;
            }
            (xnum, xden)
        } else {
            (0.0, 0.0)
        };
        let temp = x * (xnum + A[3]) / (xden + B[3]);
        if want_lower {
            cum = 0.5 + temp;
        }
        if want_upper {
            ccum = 0.5 - temp;
        }
        (cum, ccum)
    } else if y <= M_SQRT_32 {
        // 0.674.. < |x| <= sqrt(32) ~= 5.657
        let mut xnum = C[8] * y;
        let mut xden = y;
        for i in 0..7 {
            xnum = (xnum + C[i]) * y;
            xden = (xden + D[i]) * y;
        }
        let temp = (xnum + C[7]) / (xden + D[7]);
        let (c, cc) = do_del(y, temp);
        cum = c;
        ccum = cc;
        // swap_tail
        if x > 0.0 {
            let t = cum;
            if want_lower {
                cum = ccum;
            }
            ccum = t;
        }
        (cum, ccum)
    } else if (-37.5193 < x && x < 8.2924) || (-8.2924 < x && x < 37.5193) {
        // |x| > sqrt(32) but not so large that the tail is exactly 0 or 1
        let xsq = 1.0 / (x * x);
        let mut xnum = P[5] * xsq;
        let mut xden = xsq;
        for i in 0..4 {
            xnum = (xnum + P[i]) * xsq;
            xden = (xden + Q[i]) * xsq;
        }
        let mut temp = xsq * (xnum + P[4]) / (xden + Q[4]);
        temp = (M_1_SQRT_2PI - temp) / y;
        let (c, cc) = do_del(x, temp);
        cum = c;
        ccum = cc;
        if x > 0.0 {
            let t = cum;
            if want_lower {
                cum = ccum;
            }
            ccum = t;
        }
        (cum, ccum)
    } else {
        if x > 0.0 {
            return (1.0, 0.0);
        }
        (0.0, 1.0)
    }
}

/// R's `do_del` macro (non-`log_p` branch), with `d_2(t) = t/2`.
///
/// `xsq = trunc(x * 16) / 16` — the classic Cody split into a "nice" part and a correction —
/// and the two exponentials are kept **separate** (`exp(-xsq/2) * exp(-del/2) * temp`),
/// because grouping them into one `exp` changes the last bits.
#[inline]
fn do_del(x: f64, temp: f64) -> (f64, f64) {
    let xsq = (x * 16.0).trunc() / 16.0;
    let del = (x - xsq) * (x + xsq);
    let cum = (-xsq * (xsq / 2.0)).exp() * (-del / 2.0).exp() * temp;
    let ccum = 1.0 - cum;
    (cum, ccum)
}

/// R's `chebyshev_eval(x, a, n)` from `src/nmath/chebyshev.c`.
///
/// The coefficients are read in **reverse** (`a[n - i]`), which is not a detail: reading
/// them forward gives a different polynomial and a `lgammacor` that is wrong in the 15th
/// digit, which is a different `lbeta` and a different `choose`.
fn chebyshev_eval(x: f64, a: &[f64], n: usize) -> f64 {
    let twox = x * 2.0;
    let mut b0 = 0.0;
    let mut b1 = 0.0;
    let mut b2 = 0.0;
    for i in 1..=n {
        b2 = b1;
        b1 = b0;
        b0 = twox * b1 - b2 + a[n - i];
    }
    (b0 - b2) * 0.5
}

/// R's `lgammacor(x)` from `src/nmath/lgammacor.c`: Fullerton's Chebyshev correction to
/// Stirling's series, for `x >= 10`.
///
/// The `#define`s from the C source, which are the branch points:
///
/// | constant | value | role |
/// |---|---|---|
/// | `xbig` | `2^26.5` = 94906265.62425156 | below this, the Chebyshev series |
/// | `xmax` | `DBL_MAX / 48` | above this, the 1/(12x) asymptotic form |
#[allow(non_snake_case)]
const XBIG: f64 = 94906265.62425156;
#[allow(non_snake_case)]
const XMAX: f64 = 3.745194030963158e306;
const NALGM: usize = 5;

pub fn lgammacor(x: f64) -> f64 {
    // Only the first 5 of R's 15 are used; R's comment says so explicitly.
    const ALGMCS: [f64; 15] = [
        0.1666389480451863247205729650822e+0,
        -0.1384948176067563840732986059135e-4,
        0.9810825646924729426157171547487e-8,
        -0.1809129475572494194263306266719e-10,
        0.6221098041892605227126015543416e-13,
        -0.3399615005417721944303330599666e-15,
        0.2683181998482698748957538846666e-17,
        -0.2868042435334643284144622399999e-19,
        0.3962837061046434803679306666666e-21,
        -0.6831888753985766870111999999999e-23,
        0.1429227355942498147573333333333e-24,
        -0.3547598158101070547199999999999e-26,
        0.1025680058010470912000000000000e-27,
        -0.3401102254316748799999999999999e-29,
        0.1276642195630062933333333333333e-30,
    ];
    if x < 10.0 {
        return f64::NAN;
    } else if x >= XMAX {
        // R warns about underflow and then falls through to `1/(x*12)`.
    } else if x < XBIG {
        let tmp = 10.0 / x;
        return chebyshev_eval(tmp * tmp * 2.0 - 1.0, &ALGMCS, NALGM) / x;
    }
    1.0 / (x * 12.0)
}

/// `log|gamma(x)|`.
///
/// Implemented as the g = 607/128, n = 15 Lanczos series (the set glibc uses, and the same
/// one Boost and the JVM use) rather than an FFI call to libm's `lgamma`: `f64::ln_gamma`
/// is still unstable, and a `extern "C"` declaration to `lgamma` would not port to the
/// Windows CRT without a different symbol.
///
/// **Unreachable from `pwilcox`.** It is reached only from `lbeta`'s `p < 10 && q < 10`
/// branch, and [`choose`] reaches that branch only after its symmetry reductions have
/// already rewritten `n - k < 2` as `choose(n, n - k)`. So on the exact Wilcoxon path the
/// `lbeta(p, q)` call always has `p >= 10` and takes the `lgammacor` branch instead. It
/// exists so `lbeta` is total and so the function is testable on its own; a difference from
/// R here would not reach a CellChat result.
pub fn gamma_small(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x <= 0.0 && x == x.trunc() {
        // Negative integer: |gamma| has a pole, so log|gamma| = +Inf. Matches R's
        // `lgammafn`, which returns ML_POSINF here without warning.
        return f64::INFINITY;
    }
    let y = x.abs();
    if x < 0.5 {
        // Reflection, gamma(x) * gamma(1 - x) = pi / sin(pi x). Two cases, and the sign of
        // `1 - x` differs:
        //   * 0 < x < 0.5: `1 - x` is in (0.5, 1), so the series applies to it directly.
        //   * x < 0: `1 - x` is `1 + |x|`, **not** `1 - |x|`. Using `1 - |x|` evaluates
        //     gamma(|x|) instead of gamma(1 + |x|) and returns log(gamma(0.5)) for
        //     x = -0.5, which is a plausible-looking number and wrong.
        let one_minus_x = if x > 0.0 { 1.0 - x } else { 1.0 + y };
        let r = lanczos_series(one_minus_x);
        let s = (std::f64::consts::PI * x).sin().abs();
        return (std::f64::consts::PI / (s * r)).ln();
    }
    lanczos_series(x).ln()
}

/// The raw Lanczos product, `|gamma(y)|` for `y >= 0.5`.
fn lanczos_series(y: f64) -> f64 {
    const G: f64 = 607.0 / 128.0;
    const C: [f64; 15] = [
        0.99999999999999709182,
        57.156235665862923517,
        -59.597960355475491248,
        14.136097974741747174,
        -0.49191381609762019978,
        0.33994649984811888699e-4,
        0.46523628927048575665e-4,
        -0.98374475304879564677e-4,
        0.15808870322491248884e-3,
        -0.21026444172410488319e-3,
        0.21743961811521264320e-3,
        -0.16431810653676389022e-3,
        0.84418223983852743293e-4,
        -0.26190838401581408670e-4,
        0.36899182659531622704e-5,
    ];
    let z = y - 1.0;
    let mut a = C[0];
    // NB the divisor is `z + i`, not `z + 1`: `for i in 1..15 { a += c[i] / (z + i) }`.
    // A constant divisor still gives a plausible-looking number -- it just is not gamma.
    for (i, c) in C.iter().enumerate().skip(1) {
        a += c / (z + i as f64);
    }
    let t = z + G + 0.5;
    (2.0 * std::f64::consts::PI).sqrt() * t.powf(z + 0.5) * (-t).exp() * a
}

/// `lgammafn(x)` for `x > 10`: `M_LN_SQRT_2PI + (x - 0.5) * log(x) - x + lgammacor(x)`.
pub fn lgammafn(x: f64) -> f64 {
    if x.is_nan() {
        return x;
    }
    if x <= 0.0 && x == x.trunc() {
        return f64::INFINITY;
    }
    let y = x.abs();
    if y < 1e-306 {
        return -y.ln();
    }
    if y <= 10.0 {
        return gamma_small(x);
    }
    if x > 0.0 {
        if x > 1e17 {
            return x * (x.ln() - 1.0);
        } else if x > 4_934_720.0 {
            return M_LN_SQRT_2PI + (x - 0.5) * x.ln() - x;
        }
        return M_LN_SQRT_2PI + (x - 0.5) * x.ln() - x + lgammacor(x);
    }
    // x < -10
    let sinpiy = (std::f64::consts::PI * y).sin().abs();
    M_LN_SQRT_PId2 + (x - 0.5) * y.ln() - x - sinpiy.ln() - lgammacor(y)
}

/// `lbeta(a, b)` from `src/nmath/lbeta.c`.
pub fn lbeta(a: f64, b: f64) -> f64 {
    if a.is_nan() || b.is_nan() {
        return a + b;
    }
    let mut p = a;
    let mut q = a;
    if b < p {
        p = b;
    }
    if b > q {
        q = b;
    }
    if p < 0.0 {
        return f64::NAN;
    } else if p == 0.0 {
        return f64::INFINITY;
    } else if !q.is_finite() {
        return f64::NEG_INFINITY;
    }
    if p >= 10.0 {
        let corr = lgammacor(p) + lgammacor(q) - lgammacor(p + q);
        q.ln() * -0.5
            + M_LN_SQRT_2PI
            + corr
            + (p - 0.5) * (p / (p + q)).ln()
            + q * (-p / (p + q)).ln_1p()
    } else if q >= 10.0 {
        let corr = lgammacor(q) - lgammacor(p + q);
        lgammafn(p) + corr + p - p * (p + q).ln() + (q - 0.5) * (-p / (p + q)).ln_1p()
    } else {
        if p < 1e-306 {
            return gamma_small(p) + (gamma_small(q) - gamma_small(p + q));
        }
        ((gamma_small(p) + gamma_small(q)) - gamma_small(p + q))
            .exp()
            .ln()
    }
}

/// `lfastchoose(n, k) = -log(n + 1) - lbeta(n - k + 1, k + 1)`, R's `choose.c`.
fn lfastchoose(n: f64, k: f64) -> f64 {
    -(n + 1.0).ln() - lbeta(n - k + 1.0, k + 1.0)
}

/// `choose(n, k)` from `src/nmath/choose.c` (R 4.3.3), for integer `n >= 0`.
///
/// The `k < 30` branch is a plain product loop; above that R uses
/// `exp(lfastchoose(n, k))`. Both are reproduced. The symmetry reductions
/// (`n - k < k` and `n - k < 2`) are load-bearing for *which* arithmetic is used and
/// therefore for the last bits.
pub fn choose(n: f64, k: f64) -> f64 {
    if n.is_nan() || k.is_nan() {
        return n + k;
    }
    const K_SMALL_MAX: f64 = 30.0;
    let mut k = r_forceint(k);
    let mut n = n;
    if k < K_SMALL_MAX {
        if n - k < k && n >= 0.0 && !r_nonint(n) {
            k = r_forceint(n - k);
        }
        if k < 0.0 {
            return 0.0;
        }
        if k == 0.0 {
            return 1.0;
        }
        let mut r = n;
        let mut j = 2.0;
        while j <= k {
            r *= (n - j + 1.0) / j;
            j += 1.0;
        }
        return if !r_nonint(n) { r_forceint(r) } else { r };
    }
    // k >= 30
    if n < 0.0 {
        let mut r = choose(-n + k - 1.0, k);
        if k != 2.0 * (k / 2.0).floor() {
            r = -r;
        }
        return r;
    }
    if !r_nonint(n) {
        n = r_forceint(n);
        if n < k {
            return 0.0;
        }
        if n - k < K_SMALL_MAX {
            return choose(n, n - k);
        }
        return r_forceint(lfastchoose(n, k).exp());
    }
    choose(n, k)
}

#[cfg(test)]
mod tests {
    use super::*;

    /// Reference values from R 4.3.3 on the host this project targets, produced with
    /// `format(x, digits = 17)`. Regenerate with:
    ///
    /// ```sh
    /// R --vanilla -q -e 'for (x in XS) cat(sprintf("(%s, %s),\n", format(x, digits=17),
    ///     format(pnorm(x), digits=17)))'
    /// ```
    ///
    /// These are compared **bit for bit**: a more accurate `pnorm` would stop matching R,
    /// and a `pnorm` that differs in the last ulp puts a different set of genes through
    /// CellChat's significance threshold.
    const PNORM_REF: &[(f64, f64)] = &PNORM_REF_;
    const PNORM_UP_REF: &[(f64, f64)] = &PNORM_UP_REF_;
    const CHOOSE_REF: &[(f64, f64, f64)] = &CHOOSE_REF_;
    const LGAMMA_REF: &[(f64, f64)] = &LGAMMA_REF_;

    #[test]
    fn pnorm_matches_r_bit_for_bit() {
        for &(x, want) in PNORM_REF {
            let got = pnorm(x, true);
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "pnorm({x}) = {got:.17e}, R gives {want:.17e}"
            );
        }
    }

    #[test]
    fn pnorm_upper_tail_matches_r_bit_for_bit() {
        for &(x, want) in PNORM_UP_REF {
            let got = pnorm(x, false);
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "pnorm({x}, lower=FALSE) = {got:.17e}, R gives {want:.17e}"
            );
        }
        // And the substantive point: the tail is *computed*, not `1 - lower`. At x = 8 the
        // lower tail is one ulp below 1, so subtracting throws away the entire answer.
        let x = 8.0;
        let lo = pnorm(x, true);
        let up = pnorm(x, false);
        assert_ne!(up, 1.0 - lo, "the tail must not be 1 - lower");
        assert!(
            (1.0 - lo) / up > 1.05,
            "and the difference is large: {}",
            (1.0 - lo) / up
        );
    }

    #[test]
    fn choose_matches_r_bit_for_bit() {
        for &(n, k, want) in CHOOSE_REF {
            let got = choose(n, k);
            assert_eq!(
                got.to_bits(),
                want.to_bits(),
                "choose({n}, {k}) = {got:.17e}, R gives {want:.17e}"
            );
        }
    }

    #[test]
    fn lgammafn_matches_r_over_the_whole_che_by_and_asymptotic_range() {
        for &(x, want) in LGAMMA_REF {
            let got = lgammafn(x);
            assert!(
                (got - want).abs() < 1e-13 * want.abs().max(1.0),
                "lgammafn({x}) = {got}, R gives {want}"
            );
        }
        // The branch point: xbig = 2^26.5. Just below uses the Chebyshev series, just
        // above the 1/(12x) asymptotic form, and R's `lgamma` switches at the same place.
        assert!(lgammafn(94906265.624251559).is_finite());
        assert!(lgammafn(94906265.7).is_finite());
    }

    /// R 4.3.3 `lgamma`, from `format(lgamma(x), digits = 17)`.
    const LGAMMA_SMALL_REF: &[(f64, f64)] = &[
        (0.5, 0.57236494292470008),
        (1.0, 0.0),
        (1.5, -0.12078223763524518),
        (2.0, 0.0),
        (3.0, 0.69314718055994529),
        (5.0, 3.1780538303479458),
        (9.5, 11.689333420797269),
        (10.0, 12.801827480081469),
        (-0.5, 1.2655121234846454),
        (-1.5, 0.86004701537648109),
        (-2.5, -0.056243716497674033),
    ];

    #[test]
    fn gamma_small_matches_r() {
        for &(x, want) in LGAMMA_SMALL_REF {
            let got = gamma_small(x);
            assert!(
                (got - want).abs() < 1e-13 * want.abs().max(1.0),
                "gamma_small({x}) = {got}, R gives {want}"
            );
        }
        // Negative integer: |gamma| has a pole, so log|gamma| = +Inf, as R returns.
        assert!(gamma_small(-3.0).is_infinite() && gamma_small(-3.0) > 0.0);
        assert!(gamma_small(0.0).is_infinite() && gamma_small(0.0) > 0.0);
    }

    #[test]
    fn r_forceint_matches_c99_round() {
        // Rust's `f64::round` is *also* half-away-from-zero, so the two agree; the
        // difference from R is that `R_forceint` is `round(x)` and not `trunc(x)` or
        // `x as i64` (which saturates and truncates toward zero). Pinned so the choice of
        // primitive stays deliberate.
        for &x in &[2.5f64, -2.5, 2.4, -2.4, 0.5, -0.5, 3.0, 1e7 + 0.5] {
            assert_eq!(r_forceint(x), x.round(), "R_forceint({x})");
        }
        // `as i64` would differ: it truncates toward zero.
        assert_ne!(r_forceint(-2.5), -2.5f64 as i64 as f64);
    }

    #[test]
    fn chebyshev_reads_its_coefficients_in_reverse() {
        // R's loop is `b0 = twox*b1 - b2 + a[n - i]` for i = 1..n, so a 2-term series over
        // [2, 3] reads a[1] then a[0] and evaluates to 3*x + 1. Reading forward would give
        // 2*x + 3, which is a different polynomial and a `lgammacor` wrong in the 15th digit.
        let a = [2.0, 3.0];
        assert_eq!(chebyshev_eval(0.0, &a, 2), 1.0);
        assert_eq!(chebyshev_eval(1.0, &a, 2), 4.0);
    }
}
