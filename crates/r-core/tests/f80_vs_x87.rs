//! Differential test of `F80` against C's real x87 `long double`.
//!
//! The oracle is produced by `tests/parity/gen_f80_ref.c`, compiled with gcc on this
//! host. Comparing against the real thing beats reasoning about the bit layout, and it
//! is how the `sub_same_sign` exponent-underflow bug was actually found.
use r_core::longdouble::F80;

const REF: &str = include_str!("../../../tests/fixtures/f80_ref.txt");

fn approx_bits(a: f64, b: f64) -> bool {
    if a.is_nan() || b.is_nan() {
        return a.is_nan() && b.is_nan();
    }
    if a == b {
        return true;
    }
    (a - b).abs() <= 4.0 * f64::EPSILON * b.abs().max(1e-300).max(1e-30)
}

#[test]
fn matches_x87_long_double() {
    let mut n_add = 0;
    let mut n_sub = 0;
    let mut n_mul = 0;
    let mut n_prod = 0;
    let mut n_chain = 0;
    let mut lines = 0;
    for line in REF
        .lines()
        .filter(|l| !l.trim().is_empty() && !l.starts_with('#'))
    {
        let f: Vec<&str> = line.split('\t').collect();
        lines += 1;
        match f[0] {
            "add" | "sub" | "mul" => {
                let a: f64 = f[1].parse().unwrap();
                let b: f64 = f[2].parse().unwrap();
                let want: f64 = f[3].parse().unwrap();
                let got = match f[0] {
                    "add" => {
                        n_add += 1;
                        F80::from_f64(a).add(F80::from_f64(b)).to_f64()
                    }
                    "sub" => {
                        n_sub += 1;
                        F80::from_f64(a).sub(F80::from_f64(b)).to_f64()
                    }
                    _ => {
                        n_mul += 1;
                        F80::from_f64(a).mul(F80::from_f64(b)).to_f64()
                    }
                };
                assert!(
                    approx_bits(got, want),
                    "{} {} {}: rust={got:.17e} x87={want:.17e}  (ulp delta {})",
                    f[0],
                    a,
                    b,
                    ((got.to_bits() as i64) - (want.to_bits() as i64)).abs()
                );
            }
            "prod" => {
                // R's `prod` (rprod) is an LDOUBLE accumulation, which is what
                // computeExpr_coreceptor/_agonist/_antagonist end in for a multi-subunit
                // cofactor. An f64 left-to-right product is 1 ulp off on the first
                // random trial, so this branch is not optional.
                let vals: Vec<f64> = f[1].split(',').map(|v| v.parse().unwrap()).collect();
                let want: f64 = f[2].parse().unwrap();
                let mut acc = F80::from_f64(1.0);
                for &v in &vals {
                    acc = acc.mul(F80::from_f64(v));
                }
                n_prod += 1;
                let got = acc.to_f64();
                assert!(
                    approx_bits(got, want),
                    "prod of {} values: rust={got:.17e} x87={want:.17e}",
                    vals.len()
                );
            }
            "chain" => {
                // R's real_mean pass 1: repeated addition of a whole vector.
                let vals: Vec<f64> = f[1].split(',').map(|v| v.parse().unwrap()).collect();
                let want: f64 = f[2].parse().unwrap();
                let mut acc = F80::ZERO;
                for &v in &vals {
                    acc = acc.add(F80::from_f64(v));
                }
                n_chain += 1;
                assert!(
                    approx_bits(acc.to_f64(), want),
                    "chain of {} values: rust={:.17e} x87={:.17e}",
                    vals.len(),
                    acc.to_f64(),
                    want
                );
            }
            // Adversarial chains: same accumulation, but the order is observable in the answer.
            // `no_fma.rs` uses these to pin the order; here they are only checked against x87.
            "chain_cancel" => {
                let vals: Vec<f64> = f[1].split(',').map(|v| v.parse().unwrap()).collect();
                let want: f64 = f[2].parse().unwrap();
                let mut acc = F80::ZERO;
                for &v in &vals {
                    acc = acc.add(F80::from_f64(v));
                }
                n_chain += 1;
                assert!(
                    approx_bits(acc.to_f64(), want),
                    "chain_cancel of {} values: rust={:.17e} x87={:.17e}",
                    vals.len(),
                    acc.to_f64(),
                    want
                );
            }
            other => panic!("unknown ref record: {other}"),
        }
    }
    assert!(lines > 9000, "only {lines} reference cases");
    assert!(n_add > 1000 && n_sub > 1000 && n_mul > 1000 && n_prod > 100 && n_chain > 100);
    eprintln!(
        "f80 vs x87: {n_add} adds, {n_sub} subs, {n_mul} muls, {n_prod} prod chains, \
         {n_chain} sum chains"
    );
}
