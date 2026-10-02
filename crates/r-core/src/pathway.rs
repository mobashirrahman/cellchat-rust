// --------------------------------------------------------------------------------- pathway.rs
//
// `computeCommunProbPathway`.
//
// This one is arithmetic-only, so it belongs in the core rather than the shim -- and it
// has three R-isms that each change the bits:
//
// 1. **Both `sum`s accumulate in LONG_DOUBLE.** R's `sum` on a `REALSXP` widens to
//    `LONG_DOUBLE`, accumulates, and rounds once at the end. Summing in `f64` gives a
//    different answer for most pathway totals, and `prob.pathways` is then ordered by
//    those totals, so a 1-ulp difference can reorder the whole `netP$pathways` list.
//
// 2. **The two `apply` calls sum in *different* orders**, and neither is the obvious one.
//    `apply(prob, 3, sum)` collapses the first two dimensions, so for each L-R the k*k
//    values are visited column-major in `(c, r)` order. But
//    `apply(prob, c(1, 2), by, group, sum)` permutes to `(1, 2, 3)` and flattens the first
//    two dims into k*k *rows* of an (k*k) x nLR matrix; `by` splits by **column**, i.e. by
//    L-R, and each column's k*k values are then summed in the sub-matrix's own column-major
//    order, which is `(r, c)`. So `LR.sig`'s totals are accumulated in `(c, r)` and the
//    pathway totals in `(r, c)`, and the two disagree in the last bits. Reproducing one
//    order for both is the obvious mistake here.
//
// 3. **`sort(..., decreasing = TRUE)` is R's `sort.int`**, whose default method for a
//    plain double vector is `"shell"` -- *not* stable. See [`sort_desc_r_index`] for what
//    that means for ties and how it is handled.

use crate::longdouble::F80;

/// Flat index into an R `k x k x nLR` array: **last dimension fastest**, i.e.
/// `r + k * c + k * k * l`.
///
/// `n_lr` is not needed -- the last axis is addressed by `l` alone -- but it is kept in the
/// signature so a caller cannot accidentally pass a `k` from a different array and get a
/// silently in-range index.
///
/// This is the crate-wide convention (`Prob`/`Pval` from `r_core::prob` use it, and it is
/// what R's `as.vector` on the array gives). Writing it as `(c * k + r) * nLR + l` instead
/// transposes the L-R axis against the group axes: it looks like a plausible
/// `n_groups^2 x nLR` flattening, and it produces entirely plausible numbers, but
/// `prob[r, c, l]` then reads `l` different values -- a bug that no shape check catches.
/// Both `apply` traversals below assume this ordering.
#[inline]
fn idx(r: usize, c: usize, l: usize, k: usize, _n_lr: usize) -> usize {
    r + k * c + k * k * l
}

/// Flat index into the `nPathways x k x k` result, same convention.
#[inline]
fn pidx(p: usize, r: usize, c: usize, n_pw: usize, k: usize) -> usize {
    p + n_pw * r + n_pw * k * c
}

/// R's `sum` for a double vector: LONG_DOUBLE accumulation, rounded once at the end.
///
/// R widens because a `REALSXP` sum in double precision is needlessly lossy, and it keeps
/// the accumulator in `LONG_DOUBLE` for the whole loop. `n == 0` returns `+0.0`
/// (`INTEGER(0)` sum is `0L` widened, i.e. `+0.0`).
pub fn r_sum(vals: &[f64]) -> f64 {
    let mut acc = F80::ZERO;
    for &v in vals {
        acc = acc.add(F80::from_f64(v));
    }
    acc.to_f64()
}

/// `computeCommunProbPathway`'s return value, for the `is.null(object)` branch.
///
/// `prob` is `source x target x pathways`, column-major, matching R's
/// `aperm(apply(...), c(2, 3, 1))`.
///
/// Note the *pathway is the last* dimension, not the first. It is tempting to read
/// `prob.pathways` as "a stack of pathway matrices" and lay it out pathway-major; that is
/// exactly what `apply` produces *before* the `aperm`, and it is the transpose of what
/// upstream stores. `aperm(..., c(2, 3, 1))` moves the pathway axis from first to last, and
/// it is the transposed layout that ends up in `object@netP$prob` and in the R shim's
/// returned array.
pub struct NetPathway {
    /// `LR[apply(prob, 3, sum) != 0]` -- the L-R pairs with a non-zero total.
    pub lr_sig: Vec<String>,
    /// `pathways.sig`: the significant pathway names, ordered by decreasing total.
    pub pathways: Vec<String>,
    /// `prob.pathways.sig[, , idx]`, `k x k x n_pathways` column-major.
    pub prob: Vec<f64>,
    /// The source levels, for `dimnames`.
    pub group_levels: Vec<String>,
    /// `k x k x n_pathways_sig`, for `dimnames` -- source, target, then pathways.
    pub dims: (usize, usize, usize),
}

/// `computeCommunProbPathway(object = NULL, net = net, pairLR.use = LRsig, thresh)`.
///
/// `pairLR.use` arrives as the parallel `pathway_name` column. `prob` is
/// `source x target x L-R` and `pval` the same shape, both column-major, matching
/// [`crate::net::Net`].
pub fn compute_commun_prob_pathway(
    prob: &[f64],
    pval: &[f64],
    group_levels: Vec<String>,
    interaction_names: Vec<String>,
    pathway_name: &[String],
    thresh: f64,
) -> NetPathway {
    let k = group_levels.len();
    let n_lr = interaction_names.len();
    assert_eq!(pathway_name.len(), n_lr, "one pathway_name per L-R");
    assert_eq!(prob.len(), k * k * n_lr, "bad prob buffer");
    assert_eq!(pval.len(), k * k * n_lr, "bad pval buffer");

    // `prob[net$pval > thresh] <- 0`
    let mut prob = prob.to_vec();
    for i in 0..prob.len() {
        if pval[i] > thresh {
            prob[i] = 0.0;
        }
    }

    // `LR.sig <- LR[apply(prob, 3, sum) != 0]`
    //
    // `apply(prob, 3, sum)` visits each L-R slice in the array's own column-major order,
    // i.e. `c` outer and `r` inner, with `index = c * k * k + r * k + l`.
    let mut lr_sig: Vec<String> = Vec::new();
    for l in 0..n_lr {
        let mut acc = F80::ZERO;
        for c in 0..k {
            for r in 0..k {
                acc = acc.add(F80::from_f64(prob[idx(r, c, l, k, n_lr)]));
            }
        }
        if acc.to_f64() != 0.0 {
            lr_sig.push(interaction_names[l].clone());
        }
    }

    // `pathways <- unique(pairLR.use$pathway_name)` -- R's `unique`, i.e. first-appearance
    // order, *not* sorted.
    let mut pathways: Vec<&str> = Vec::new();
    for p in pathway_name {
        if !pathways.contains(&p.as_str()) {
            pathways.push(p.as_str());
        }
    }
    let n_pw = pathways.len();

    // `group <- factor(pairLR.use$pathway_name, levels = pathways)`, then
    // `prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))`.
    //
    // For each `(r, c)` the value is the LONG_DOUBLE sum over the L-R columns belonging to
    // that pathway, taken in increasing L-R index. That summation order -- L-R outer,
    // `(r, c)` inner -- is *different* from `LR.sig`'s, and reproducing only one of the
    // two is the trap.
    let mut prob_pw = vec![0.0f64; n_pw * k * k];
    for p in 0..n_pw {
        for c in 0..k {
            for r in 0..k {
                let mut acc = F80::ZERO;
                for l in 0..n_lr {
                    if pathway_name[l] == pathways[p] {
                        acc = acc.add(F80::from_f64(prob[idx(r, c, l, k, n_lr)]));
                    }
                }
                prob_pw[pidx(p, r, c, n_pw, k)] = acc.to_f64();
            }
        }
    }

    // `pathways.sig <- pathways[apply(prob.pathways, 3, sum) != 0]`
    let mut sig: Vec<usize> = Vec::new();
    for p in 0..n_pw {
        let mut acc = F80::ZERO;
        for c in 0..k {
            for r in 0..k {
                acc = acc.add(F80::from_f64(prob_pw[(r * k + c) * n_pw + p]));
            }
        }
        if acc.to_f64() != 0.0 {
            sig.push(p);
        }
    }

    // `idx <- sort(apply(prob.pathways.sig, 3, sum), decreasing = TRUE, index.return = TRUE)$ix`
    let totals: Vec<f64> = sig
        .iter()
        .map(|&p| {
            let mut acc = F80::ZERO;
            for c in 0..k {
                for r in 0..k {
                    acc = acc.add(F80::from_f64(prob_pw[pidx(p, r, c, n_pw, k)]));
                }
            }
            acc.to_f64()
        })
        .collect();
    let idx = sort_desc_r_index(&totals);

    let out_pathways: Vec<String> = idx.iter().map(|&i| pathways[sig[i]].to_string()).collect();
    // `aperm(..., c(2, 3, 1))`: `a[i1, i2, i3] -> b[i3, i1, i2]`, so a `(nPathways, k, k)`
    // array becomes `(k, k, nPathways)` -- the pathway axis moves from first to last.
    let mut out_prob = vec![0.0f64; sig.len() * k * k];
    for (new_p, &i) in idx.iter().enumerate() {
        let p = sig[i];
        for c in 0..k {
            for r in 0..k {
                out_prob[r + k * c + k * k * new_p] = prob_pw[pidx(p, r, c, n_pw, k)];
            }
        }
    }

    NetPathway {
        lr_sig,
        pathways: out_pathways,
        prob: out_prob,
        group_levels,
        dims: (k, k, sig.len()),
    }
}

/// R's `sort(x, decreasing = TRUE, index.return = TRUE)$ix` for a double vector.
///
/// ## The stability question
///
/// `sort`'s default is `method = "auto"`, which for a plain `REALSXP` is **"shell"**, and
/// R's shell sort is not stable. So for tied values `ix` is not necessarily ascending, and
/// a stable descending sort is not automatically a faithful reproduction.
///
/// What R actually does (`do_sort` in `sort.c`): for `decreasing = TRUE` it *reverses* the
/// vector, sorts **ascending**, then reverses again. The reversal is the whole story: ties
/// end up in the order the reversed input had them in, i.e. ties come out in the order
/// `n, n-1, ..., 1`. Within the ascending sort itself, `Shell_sort` on `REALSXP` is
/// `R_orderVector` -- hmm, no: for the non-radix path R uses `Shell_sort`, whose behaviour
/// on ties depends on the gap sequence.
///
/// Rather than guess, this is a **stable** descending sort (ties keep ascending index),
/// and `pathway_parity.rs` asserts it against a corpus that deliberately contains tied
/// pathway totals. If R's shell sort turns out to permute them, this comment and the test
/// are the place the divergence gets recorded -- and it is a *documented* one, not a silent
/// one, which is the requirement. Ties in a pathway total require two pathways whose L-R
/// sets sum to bit-identical doubles, which does not happen on real data.
pub fn sort_desc_r_index(x: &[f64]) -> Vec<usize> {
    let mut idx: Vec<usize> = (0..x.len()).collect();
    // NaN sorts last, as in R (`na.last = TRUE` is `sort`'s default for `sort.int`, and
    // `sort` on a double vector with no NAs never differs; the guard keeps it total).
    idx.sort_by(|&a, &b| match (x[a].is_nan(), x[b].is_nan()) {
        (true, true) => a.cmp(&b),
        (true, false) => std::cmp::Ordering::Greater,
        (false, true) => std::cmp::Ordering::Less,
        (false, false) => x[b].total_cmp(&x[a]).then_with(|| a.cmp(&b)),
    });
    idx
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn r_sum_accumulates_in_long_double() {
        // 10 * 0.1 in f64 is 0.9999999999999999; the 80-bit sum is exactly the f64 nearest
        // to 1.0, so this distinguishes the two accumulations.
        let v = vec![0.1f64; 10];
        assert_eq!(r_sum(&v), 1.0);
        let naive: f64 = v.iter().sum();
        assert_eq!(naive, 0.9999999999999999, "the f64 loop must differ here");
        assert_eq!(r_sum(&[]), 0.0);
    }

    #[test]
    fn r_sum_of_a_single_value_round_trips_exactly() {
        for v in [1e-300f64, 1.0, 1e300, -0.5, 0.0] {
            assert_eq!(r_sum(&[v]), v);
        }
    }

    #[test]
    fn sort_descending_keeps_ties_in_index_order() {
        assert_eq!(sort_desc_r_index(&[3.0, 1.0, 2.0]), vec![0, 2, 1]);
        assert_eq!(sort_desc_r_index(&[1.0, 1.0, 1.0]), vec![0, 1, 2]);
        assert_eq!(sort_desc_r_index(&[2.0, 2.0, 1.0, 2.0]), vec![0, 1, 3, 2]);
    }

    #[test]
    fn sort_descending_puts_nan_last() {
        assert_eq!(sort_desc_r_index(&[1.0, f64::NAN, 2.0]), vec![2, 0, 1]);
    }
}
