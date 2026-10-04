//! `computeCommunProbPathway` against a corpus generated from the pinned upstream
//! (`tests/parity/gen_pathway_golden.R`). Inputs are dumped to
//! `tests/fixtures/pathway_inputs.tsv` as `%a` hex rather than regenerated from
//! `set.seed`/`runif`, per the rule in `tests/parity/README.md`.
//!
//! Three R-isms are load-bearing and each has a spec aimed at it:
//!   * the LONG_DOUBLE `sum` -- `bigmag`, whose upstream totals are `30000000000000004`
//!     where an `f64` loop gives `30000000000000000`;
//!   * the two *different* summation orders of the two `apply` calls, since `LR.sig` and
//!     the pathway ranking come from sums taken in opposite orders;
//!   * `sort(decreasing = TRUE)` on tied totals -- `tied` and `tied3`, which settle whether
//!     a stable descending sort reproduces R's shell sort.

use r_core::pathway::compute_commun_prob_pathway;

const GOLDEN: &str = include_str!("../../../../../tests/fixtures/pathway_golden.txt");
const INPUTS: &str = include_str!("../../../../../tests/fixtures/pathway_inputs.tsv");

fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn parse_vec(s: &str) -> Vec<f64> {
    if s.is_empty() {
        return Vec::new();
    }
    s.split(',')
        .map(|v| {
            let v = v.trim();
            assert_ne!(v, "NA", "unexpected NA in the pathway corpus");
            v.parse().unwrap_or_else(|_| panic!("not a number: {v:?}"))
        })
        .collect()
}

/// `sprintf("%a", x)` for a non-negative-or-signed double, e.g. `0x1.8p+1`, `0x0p+0`.
fn parse_hex(s: &str) -> f64 {
    let s = s.trim();
    let (neg, s) = match s.strip_prefix('-') {
        Some(r) => (true, r),
        None => (false, s.trim_start_matches('+')),
    };
    let s = s.trim_start_matches("0x");
    let (mantissa, exp) = match s.find(['p', 'P']) {
        Some(i) => (&s[..i], s[i + 1..].parse::<i32>().unwrap_or(0)),
        None => (s, 0),
    };
    let (int_part, frac_part) = match mantissa.find('.') {
        Some(i) => (&mantissa[..i], &mantissa[i + 1..]),
        None => (mantissa, ""),
    };
    let mut v = 0.0f64;
    for b in int_part.bytes() {
        v = v * 16.0 + (b as char).to_digit(16).unwrap() as f64;
    }
    let mut scale = 1.0f64;
    for b in frac_part.bytes() {
        scale /= 16.0;
        v += (b as char).to_digit(16).unwrap() as f64 * scale;
    }
    let out = v * 2f64.powi(exp);
    if neg {
        -out
    } else {
        out
    }
}

struct Spec {
    name: String,
    k: usize,
    n_lr: usize,
    levels: Vec<String>,
    lr: Vec<String>,
    pathways: Vec<String>,
    thresh: f64,
    prob: Vec<f64>,
    pval: Vec<f64>,
}

fn specs() -> Vec<Spec> {
    let mut out: Vec<Spec> = Vec::new();
    let mut it = INPUTS.lines();
    while let Some(l) = it.next() {
        if !l.starts_with("spec\t") {
            continue;
        }
        let f = fields(l);
        let name = f[1].to_string();
        let k: usize = f[2].parse().unwrap();
        let n_lr: usize = f[3].parse().unwrap();
        let levels: Vec<String> = f[4].split(',').map(|s| s.to_string()).collect();
        let lr: Vec<String> = f[5].split(',').map(|s| s.to_string()).collect();
        let pathways = it
            .next()
            .and_then(|x| x.strip_prefix(&format!("pathways\t{name}\t")))
            .expect("pathways record")
            .split(',')
            .map(|s| s.to_string())
            .collect();
        let thresh: f64 = it
            .next()
            .and_then(|x| x.strip_prefix(&format!("thresh\t{name}\t")))
            .expect("thresh record")
            .trim()
            .parse()
            .unwrap();
        assert_eq!(it.next(), Some("prob"), "{name}: expected prob");
        let prob: Vec<f64> = it
            .next()
            .expect("prob row")
            .split_whitespace()
            .map(parse_hex)
            .collect();
        assert_eq!(it.next(), Some("pval"), "{name}: expected pval");
        let pval: Vec<f64> = it
            .next()
            .expect("pval row")
            .split_whitespace()
            .map(parse_hex)
            .collect();
        assert_eq!(prob.len(), k * k * n_lr, "{name}: prob length");
        assert_eq!(pval.len(), k * k * n_lr, "{name}: pval length");
        assert_eq!(levels.len(), k);
        assert_eq!(lr.len(), n_lr);
        out.push(Spec {
            name,
            k,
            n_lr,
            levels,
            lr,
            pathways,
            thresh,
            prob,
            pval,
        });
    }
    out
}

impl Spec {
    /// `pathway_name` for each L-R: two consecutive L-Rs per pathway, in the order the
    /// `pathways` record gives.
    fn pathway_of(&self, l: usize) -> String {
        self.pathways[l / 2].clone()
    }
    fn run(&self) -> r_core::pathway::NetPathway {
        let pw: Vec<String> = (0..self.n_lr).map(|l| self.pathway_of(l)).collect();
        compute_commun_prob_pathway(
            &self.prob,
            &self.pval,
            self.levels.clone(),
            self.lr.clone(),
            &pw,
            self.thresh,
        )
    }
}

/// `k1` is the one spec upstream *errors* on, so it has an `ERROR` record instead of the
/// per-quantity records. See `single_pathway_errors_like_upstream`.
fn upstream_errors(spec: &str) -> bool {
    spec == "k1"
}

fn record(spec: &str, kind: &str) -> Option<String> {
    GOLDEN
        .lines()
        .find(|l| l.starts_with(&format!("pathway\t{spec}\t{kind}=")))
        .map(|l| {
            let f = fields(l);
            let v = &f[2];
            v.strip_prefix(&format!("{kind}=")).unwrap_or(v).to_string()
        })
}

fn bits_eq(got: &[f64], want: &[f64], what: &str) {
    assert_eq!(
        got.len(),
        want.len(),
        "{what}: length {} vs R's {}",
        got.len(),
        want.len()
    );
    for (i, (&g, &w)) in got.iter().zip(want).enumerate() {
        if g.is_nan() || w.is_nan() {
            assert!(g.is_nan() && w.is_nan(), "{what}[{i}]: rust={g} r={w}");
            continue;
        }
        if g.to_bits() != w.to_bits() {
            panic!(
                "{what}[{i}]: bit mismatch\n  rust = {g:.17e} (0x{:016x})\n  R    = {w:.17e} (0x{:016x})",
                g.to_bits(),
                w.to_bits()
            );
        }
    }
}

#[test]
fn the_corpus_covers_everything_the_test_reads() {
    let specs = specs();
    assert_eq!(specs.len(), 9, "the corpus should carry 9 specs");
    for sp in &specs {
        if upstream_errors(&sp.name) {
            assert!(
                GOLDEN
                    .lines()
                    .any(|l| l.starts_with(&format!("pathway\t{}\tERROR\t", sp.name))),
                "{}: expected an ERROR record",
                sp.name
            );
            continue;
        }
        for kind in [
            "pathways",
            "lr_sig",
            "dim",
            "prob",
            "dimnames",
            "bare",
            "bare_prob",
            "lr_sums",
            "pw_sums",
            "pwp_aperm",
            "pwp_dim",
        ] {
            assert!(
                record(&sp.name, kind).is_some(),
                "{}: no `{kind}=` record",
                sp.name
            );
        }
    }
    // The specs that make the three R-isms observable must all be present, or the test
    // silently stops testing them.
    for needed in ["bigmag", "tied", "tied3", "allzero", "k1", "rand4"] {
        assert!(
            specs.iter().any(|s| s.name == needed),
            "missing spec {needed}"
        );
    }
}

#[test]
fn the_pathway_ranking_matches_r() {
    for sp in specs() {
        if upstream_errors(&sp.name) {
            continue;
        }
        let got = sp.run();
        let want = record(&sp.name, "pathways").expect("pathways record");
        assert_eq!(got.pathways.join(","), want, "{}: netP$pathways", sp.name);
        let dim = record(&sp.name, "dim").expect("dim record");
        let want_dim: Vec<usize> = dim.split('x').map(|s| s.parse().unwrap()).collect();
        assert_eq!(
            vec![got.dims.0, got.dims.1, got.dims.2],
            want_dim,
            "{}: dim(netP$prob)",
            sp.name
        );
    }
}

#[test]
fn the_lrsig_list_matches_r() {
    for sp in specs() {
        if upstream_errors(&sp.name) {
            continue;
        }
        let got = sp.run();
        let want = record(&sp.name, "lr_sig").expect("lr_sig record");
        assert_eq!(got.lr_sig.join(","), want, "{}: net$LRs", sp.name);
    }
}

#[test]
fn the_pathway_probability_array_matches_r_bit_for_bit() {
    for sp in specs() {
        if upstream_errors(&sp.name) {
            continue;
        }
        let got = sp.run();
        let want = parse_vec(&record(&sp.name, "prob").expect("prob record"));
        bits_eq(&got.prob, &want, &format!("{}: netP$prob", sp.name));
        // The `object = NULL` branch must agree with the `object` branch.
        let bare = parse_vec(&record(&sp.name, "bare_prob").expect("bare_prob record"));
        bits_eq(&got.prob, &bare, &format!("{}: bare list()$prob", sp.name));
        let bare_pw = record(&sp.name, "bare").expect("bare record");
        assert_eq!(
            got.pathways.join(","),
            bare_pw,
            "{}: bare pathways",
            sp.name
        );
    }
}

#[test]
fn the_dimnames_are_the_group_levels_and_the_pathways() {
    for sp in specs() {
        if upstream_errors(&sp.name) {
            continue;
        }
        let got = sp.run();
        let want = record(&sp.name, "dimnames").expect("dimnames record");
        // `pathway|source|target` per dimension, comma-joined across dimensions.
        let want_dims: Vec<Vec<&str>> = want.split(',').map(|d| d.split('|').collect()).collect();
        assert_eq!(want_dims.len(), 3, "{}: dimnames record shape", sp.name);
        // Upstream leaves `netP$prob`'s dimnames as the *permuted* `(k, k, nPathways)`
        // array's: `apply(prob, c(1,2), by, group, sum)` produces an array whose first two
        // dims are the `k x k` flattening of the group matrix, so after
        // `aperm(..., c(2, 3, 1))` the *pathways* land in the **third** dimnames slot and
        // the two group levels in the first two. Not `list(pathways, g, g)`, which is the
        // reading the code's variable names suggest.
        let lev: Vec<&str> = sp.levels.iter().map(|s| s.as_str()).collect();
        assert_eq!(want_dims[0], &lev[..], "{}: dimnames[[1]]", sp.name);
        assert_eq!(want_dims[1], &lev[..], "{}: dimnames[[2]]", sp.name);
        // With no surviving pathway, R's `dimnames` third entry is a zero-length vector,
        // which the corpus serialises as an empty field -- so expect one empty string, not
        // zero elements. Pinned so the empty-network case stays distinguishable from a
        // network whose pathways were all dropped for a different reason.
        if got.pathways.is_empty() {
            assert_eq!(
                want_dims[2],
                &[""],
                "{}: dimnames[[3]] for an empty result",
                sp.name
            );
        } else {
            assert_eq!(
                want_dims[2],
                &got.pathways[..],
                "{}: dimnames[[3]]",
                sp.name
            );
        }
    }
}

/// The two `apply` sums, pinned separately.
///
/// `lr_sums` is `apply(prob, 3, sum)`: for each L-R, the `k*k` values in `(c, r)` order.
/// `pw_sums` is `apply(apply(prob, c(1,2), by, group, sum), 3, sum)`: for each pathway, the
/// L-R columns in increasing index, each column in `(r, c)` order.
///
/// Checking them independently is what makes the two-orders claim testable: if the port
/// used one order for both, one of these two lines would move.

#[test]
fn the_two_apply_sums_match_r_in_their_own_orders() {
    use r_core::longdouble::F80;
    use r_core::pathway::r_sum;
    for sp in specs() {
        if upstream_errors(&sp.name) {
            continue;
        }
        let n_lr = sp.n_lr;
        let k = sp.k;
        let mut prob = sp.prob.clone();
        for i in 0..prob.len() {
            if sp.pval[i] > sp.thresh {
                prob[i] = 0.0;
            }
        }
        // apply(prob, 3, sum): c outer, r inner.
        let mut lr_sums = Vec::with_capacity(n_lr);
        for l in 0..n_lr {
            let mut acc = F80::ZERO;
            for c in 0..k {
                for r in 0..k {
                    acc = acc.add(F80::from_f64(prob[r + k * c + k * k * l]));
                }
            }
            lr_sums.push(acc.to_f64());
        }
        bits_eq(
            &lr_sums,
            &parse_vec(&record(&sp.name, "lr_sums").expect("lr_sums record")),
            &format!("{}: apply(prob, 3, sum)", sp.name),
        );
        // Upstream is
        //   pwp <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))
        //   pw_sums <- apply(pwp, 3, sum)
        // and the `aperm` puts the pathway on the **last** axis, so `pwp` is `k x k x
        // nPathways` -- the same layout as `netP$prob`, only without the total-based
        // reordering. (Before the `aperm` the pathway is *first*.)
        //
        // `apply(pwp, 3, sum)` then sums dims 1 and 2: for fixed `p`, the k*k values at
        // `r + k * c + k * k * p` in index order, i.e. `c` outer and `r` inner.
        //
        // Note this is a *flat* sum of the k*k already-summed values, not a sum of the
        // per-(r, c) sums: nesting the long-double accumulator gives a different (also
        // defensible) result, and R does not do it.
        let n_pw = sp.pathways.len();
        let mut pwp = vec![0.0f64; n_pw * k * k];
        for p in 0..n_pw {
            for c in 0..k {
                for r in 0..k {
                    let mut inner = F80::ZERO;
                    for l in 0..n_lr {
                        if sp.pathway_of(l) == sp.pathways[p] {
                            inner = inner.add(F80::from_f64(prob[r + k * c + k * k * l]));
                        }
                    }
                    pwp[r + k * c + k * k * p] = inner.to_f64();
                }
            }
        }
        let mut pw_sums = vec![0.0f64; n_pw];
        for (p, out) in pw_sums.iter_mut().enumerate() {
            let mut acc = F80::ZERO;
            for c in 0..k {
                for r in 0..k {
                    acc = acc.add(F80::from_f64(pwp[r + k * c + k * k * p]));
                }
            }
            *out = acc.to_f64();
        }
        // The intermediate, pinned directly so a `netP$prob` mismatch can be localised to
        // the summation or to the permutation/reordering.
        bits_eq(
            &pwp,
            &parse_vec(&record(&sp.name, "pwp_aperm").expect("pwp_aperm record")),
            &format!(
                "{}: aperm(apply(prob, c(1,2), by, group, sum), c(2,3,1))",
                sp.name
            ),
        );
        bits_eq(
            &pw_sums,
            &parse_vec(&record(&sp.name, "pw_sums").expect("pw_sums record")),
            &format!(
                "{}: apply(apply(prob, c(1,2), by, group, sum), 3, sum)",
                sp.name
            ),
        );
        // And the long-double sum really is needed: at least one spec must have a total
        // that an f64 loop gets wrong, or this test proves nothing about it.
        if sp.name == "bigmag" {
            let naive: f64 = lr_sums.iter().sum();
            assert_ne!(
                naive, lr_sums[0],
                "bigmag's lr_sums must differ between f64 and LDOUBLE accumulation"
            );
        }
        let _ = r_sum(&lr_sums);
    }
}

/// `k1` has a single L-R, hence a single pathway, and upstream then fails:
/// `apply(prob, c(1, 2), by, group, sum)` returns a **matrix** (not an array) when there
/// is one group, so `aperm(x, c(2, 3, 1))` rejects the permutation:
///
///     Error in aperm(...) : 'perm' is of wrong length 3 (!= 2)
///
/// That is a real, reachable upstream failure -- `computeCommunProbPathway` on a
/// one-pathway `LRsig` -- and a drop-in replacement has to fail the same way, so the R
/// shim routes it to upstream rather than guessing. Pinned here so the behaviour is
/// recorded rather than discovered later.

#[test]
fn single_pathway_errors_like_upstream() {
    let sp = specs()
        .into_iter()
        .find(|s| s.name == "k1")
        .expect("k1 spec");
    assert_eq!(sp.n_lr, 1);
    assert_eq!(sp.pathways.len(), 1);
    let got = GOLDEN
        .lines()
        .find(|l| l.starts_with("pathway\tk1\tERROR\t"))
        .map(|l| fields(l)[3].to_string());
    assert_eq!(
        got.as_deref(),
        Some("'perm' is of wrong length 3 (!= 2)"),
        "upstream's single-pathway failure changed; re-read the pinned source"
    );
    // The core itself does not raise: it computes the value. The shim is what must not
    // claim this case, which `check_identical.R` checks against upstream directly.
    let r = sp.run();
    assert_eq!(r.pathways, vec!["PW1"]);
}
