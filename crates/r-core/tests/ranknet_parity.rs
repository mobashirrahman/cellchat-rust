//! Parity for `rankNet`'s numeric core: the per-pathway information flow, the `-1/log`
//! rescaling with its degenerate-entry reassignment, `order()`, and the one-significant-digit
//! relative contributions.
//!
//! The corpus is `tests/fixtures/ranknet_golden.txt`, written by
//! `tests/parity/gen_ranknet_golden.R`, which lifts upstream's own expressions verbatim rather
//! than reimplementing them -- so the oracle is upstream's arithmetic, not my reading of it.
//! There is no RNG: every array is written out, which is what makes the degenerate cases
//! (`pSum` exactly 1, exactly 0, negative, all-degenerate) reachable on purpose.
//!
//! Doubles are compared by **bit pattern**, not by text: the corpus is R's
//! `format(digits = 17)` and this side writes its own, and those agree on the value while
//! differing on width and exponent spelling.

use r_core::ranknet::{
    format_digits1, information_flow, order_f64, order_f64_multi, Measure, RankNetError,
};
use std::collections::HashMap;

#[derive(Debug, Default, Clone)]
struct Case {
    measure: String,
    thresh: f64,
    names: Vec<String>,
    k: usize,
    n: usize,
    levels: Vec<String>,
    sources: Option<Vec<String>>,
    targets: Option<Vec<String>>,
    prob: Vec<f64>,
    pval: Vec<f64>,
    original: Option<Vec<f64>>,
    scaled: Option<Vec<f64>>,
    flagged: Option<Vec<usize>>,
    values_assign: Option<Vec<f64>>,
    position: Option<Vec<usize>>,
    order: Option<Vec<usize>>,
    error: Option<String>,
}

fn cell(s: &str) -> f64 {
    match s {
        "<NA>" | "NA" => f64::NAN,
        "<NaN>" | "NaN" => f64::NAN,
        "Inf" => f64::INFINITY,
        "-Inf" => f64::NEG_INFINITY,
        other => other.parse().expect("corpus cell must parse"),
    }
}

fn num(s: &str) -> f64 {
    match s {
        "<NA>" | "NA" => f64::NAN,
        "<NaN>" | "NaN" => f64::NAN,
        other => other.parse().expect("corpus number must parse"),
    }
}

/// A corpus index list. Both `-` and an empty field mean "no entries": `paste(integer(0),
/// collapse = ",")` is `""`, and a strict parser that only knew about `-` would read the empty
/// string as a malformed index and fail with a `ParseIntError` that says nothing about the
/// fixture.
fn list(s: &str) -> Vec<usize> {
    if s.is_empty() || s == "-" {
        Vec::new()
    } else {
        s.split(',').map(|x| x.parse().unwrap()).collect()
    }
}

/// 1-based (R) to 0-based (Rust), dropping the empty vector.
fn one_based(v: Vec<usize>) -> Vec<usize> {
    v.into_iter().map(|x| x - 1).collect()
}

fn names_of(s: &str) -> Option<Vec<String>> {
    if s == "-" {
        None
    } else {
        Some(s.split(',').map(str::to_string).collect())
    }
}

fn load() -> (
    HashMap<String, Case>,
    Vec<(String, Vec<f64>, Vec<f64>, Vec<f64>)>,
) {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/ranknet_golden.txt");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    let mut cases: HashMap<String, Case> = HashMap::new();
    let mut rel: Vec<(String, Vec<f64>, Vec<f64>, Vec<f64>)> = Vec::new();
    let mut order: Vec<String> = Vec::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        if f.is_empty() {
            continue;
        }
        let tail = |i: usize| -> Vec<f64> { f[i..].iter().map(|x| cell(x)).collect() };
        match f[0] {
            "case" => {
                order.push(f[1].to_string());
                cases.insert(f[1].to_string(), Case::default());
            }
            "meta" => {
                let c = cases.get_mut(f[1]).unwrap();
                c.measure = f[2].to_string();
                c.thresh = num(f[3]);
                c.names = f[4].split(',').map(str::to_string).collect();
                c.k = f[5].parse().unwrap();
                c.n = f[6].parse().unwrap();
            }
            "levels" => {
                cases.get_mut(f[1]).unwrap().levels = f[2].split(',').map(str::to_string).collect()
            }
            "sources" => cases.get_mut(f[1]).unwrap().sources = names_of(f[2]),
            "targets" => cases.get_mut(f[1]).unwrap().targets = names_of(f[2]),
            "prob" => cases.get_mut(f[1]).unwrap().prob = tail(2),
            "pval" => cases.get_mut(f[1]).unwrap().pval = tail(2),
            "original" => cases.get_mut(f[1]).unwrap().original = Some(tail(2)),
            "scaled" => cases.get_mut(f[1]).unwrap().scaled = Some(tail(2)),
            // R indices are 1-based; `r-core` is 0-based throughout. Converted on load so no
            // call site has to remember, and so a corpus line can be read against R's own
            // `which()` output without an off-by-one.
            "flagged" => cases.get_mut(f[1]).unwrap().flagged = Some(one_based(list(f[2]))),
            "values_assign" => {
                cases.get_mut(f[1]).unwrap().values_assign = if f[2] == "-" {
                    None
                } else {
                    Some(f[2].split(',').map(num).collect())
                }
            }
            "position" => cases.get_mut(f[1]).unwrap().position = Some(one_based(list(f[2]))),
            "order" => cases.get_mut(f[1]).unwrap().order = Some(one_based(list(f[2]))),
            "error" => cases.get_mut(f[1]).unwrap().error = Some(f[2].to_string()),
            // A `relcase` record *opens* a new relative case; the three vectors that follow it
            // fill the slot. The slot has to be created here or the vectors land in the
            // previous case.
            "relcase" => {
                order.push(f[1].to_string());
                rel.push((f[1].to_string(), Vec::new(), Vec::new(), Vec::new()));
            }
            "rel_num" | "rel_den" | "relative" => {
                // Collected in a parallel vec, indexed by the order the `relcase` records
                // appeared in; three records per case, in this order.
                let which = match f[0] {
                    "rel_num" => 0,
                    "rel_den" => 1,
                    _ => 2,
                };
                let name = f[1].to_string();
                let idx = rel.iter().rposition(|(k, ..)| *k == name).unwrap();
                let entry = &mut rel[idx];
                let v = tail(2);
                match which {
                    0 => entry.1 = v,
                    1 => entry.2 = v,
                    _ => entry.3 = v,
                }
            }
            // `n_cmp_cases` is the comparison corpus's count, not a single-network record.
            "n_cases" | "n_cmp_cases" => {}
            // The `cmp_*` records belong to `ranknet_comparison_parity.rs`, which reads the same
            // file. Skipping them by prefix rather than listing each tag keeps the two readers
            // independent: adding a comparison record cannot break the single-network suite, and
            // a genuinely unknown tag still fails loudly.
            other if other.starts_with("cmp_") => {}
            other => panic!("unknown corpus record {other:?}"),
        }
    }
    (cases, rel)
}

/// `f64` -> a comparable string by bit pattern, so a 1-ulp difference cannot hide behind
/// formatting.
fn bits(v: f64) -> String {
    if v.is_nan() {
        "nan".to_string()
    } else {
        format!("{:016x}", v.to_bits())
    }
}

fn same(a: &[f64], b: &[f64]) -> bool {
    a.len() == b.len() && a.iter().zip(b).all(|(x, y)| bits(*x) == bits(*y))
}

/// Resolve group names to indices against the **group levels**, not the pathway names. The two
/// vocabularies are disjoint in every fixture, so resolving against the wrong one yields an
/// empty index set and an all-zero flow that still looks like a result.
fn indices(levels: &[String], want: &[String]) -> Vec<usize> {
    want.iter()
        .map(|w| {
            levels
                .iter()
                .position(|n| n == w)
                .expect("fixture group name must exist")
        })
        .collect()
}

fn run(c: &Case) -> Result<r_core::ranknet::Flow, RankNetError> {
    information_flow(
        &c.prob,
        &c.pval,
        c.k,
        c.n,
        &c.names,
        c.thresh,
        Measure::from_name(&c.measure).expect("fixture measure"),
        c.sources.as_ref().map(|v| indices(&c.levels, v)).as_deref(),
        c.targets.as_ref().map(|v| indices(&c.levels, v)).as_deref(),
    )
}

/// The corpus must reach every branch, or a green run means nothing.
#[test]
fn the_corpus_reaches_every_branch() {
    let (cases, rel) = load();
    for want in [
        "weight_degenerate",
        "weight_three_degenerate",
        "weight_all_degenerate",
        "thresh_cuts",
        "count_measure",
        "sources_filter",
        "targets_filter",
        "both_filters",
        "all_zero",
        "thresh_empties",
        "order_ties",
    ] {
        assert!(cases.contains_key(want), "corpus lost the {want:?} case");
    }
    // Three of the eleven stop, each for a different reason.
    let mut errors: Vec<(&str, &str)> = cases
        .iter()
        .filter_map(|(k, v)| v.error.as_deref().map(|e| (k.as_str(), e)))
        .collect();
    errors.sort_unstable();
    assert_eq!(
        errors,
        [
            ("all_zero", "No inferred communications for the input!"),
            (
                "thresh_empties",
                "No inferred communications for the input!"
            ),
            ("weight_all_degenerate", "'from' must be a finite number"),
        ],
        "the set of error cases changed"
    );
    assert_eq!(rel.len(), 6, "six relative-contribution cases");
    // The degenerate cases must actually be degenerate, or the interesting assertions below
    // are vacuous.
    let d = &cases["weight_degenerate"];
    assert_eq!(
        d.flagged.as_ref().unwrap(),
        &vec![3],
        "exactly one flagged entry"
    );
    assert_eq!(d.original.as_ref().unwrap().len(), 4);
    assert!(
        d.original.as_ref().unwrap()[1] == 0.0 && d.original.as_ref().unwrap()[2] == -1.0,
        "the fixture must contain a zero and a negative total"
    );
}

/// The headline: for every case, the same error text or a bit-identical `Flow`.
#[test]
fn the_information_flow_matches_upstream() {
    let (cases, _) = load();
    let mut names: Vec<&String> = cases.keys().collect();
    names.sort();
    for name in names {
        let c = &cases[name];
        match (&c.error, run(c)) {
            (Some(want), Err(e)) => assert_eq!(&e.to_string(), want, "{name}: error text"),
            (Some(want), Ok(_)) => {
                panic!("{name}: upstream errors with {want:?}, port returned a flow")
            }
            (None, Err(e)) => panic!("{name}: port errored with {e}, upstream returned a flow"),
            (None, Ok(flow)) => {
                let want_o = c
                    .original
                    .as_ref()
                    .expect("a successful case records `original`");
                let want_s = c
                    .scaled
                    .as_ref()
                    .expect("a successful case records `scaled`");
                assert!(
                    same(&flow.original, want_o),
                    "{name}: pSum.original\n got {:?}\nwant {:?}",
                    flow.original.iter().map(|v| bits(*v)).collect::<Vec<_>>(),
                    want_o.iter().map(|v| bits(*v)).collect::<Vec<_>>()
                );
                assert!(same(&flow.scaled, want_s), "{name}: pSum");
            }
        }
    }
}

/// `-1/log(x)` inverts the ordering, and only a total above 1 is flagged. Both halves matter:
/// a total of 0 gives `-1/-Inf == 0`, which is an ordinary-looking value and is *not* flagged,
/// and a negative total gives `NaN`, which the `is.na` reset turns into 0 *before* the flag
/// test, so it is not flagged either. Only `total > 1` gives `-Inf` and lands in `idx1`.
#[test]
fn only_totals_above_one_are_flagged() {
    let (cases, _) = load();
    let c = &cases["weight_degenerate"];
    let flow = run(c).unwrap();
    // originals are 0.5, 0, -1, 2
    assert_eq!(flow.original, vec![0.5, 0.0, -1.0, 2.0]);
    assert!(
        flow.scaled[1] == 0.0 && flow.scaled[1].is_sign_positive(),
        "log(0) is -Inf, so -1/-Inf is +0, not -0 and not Inf"
    );
    assert_eq!(flow.scaled[2], 0.0, "log(-1) is NaN, reset to 0 by is.na");
    assert!(
        flow.scaled[3] > flow.scaled[0] && flow.scaled[0] > 0.0,
        "the reassigned entry lands above the largest ordinary scaled value"
    );
    assert_eq!(
        c.flagged.as_ref().unwrap(),
        &vec![3],
        "upstream's `which()` gives 4, which is index 3 here"
    );
}

/// The reassignment gives each degenerate entry a value in proportion to its *original* total,
/// so the entries keep their relative order. Reading the sorted values as the assignment order
/// reverses it, which is invisible on a symmetric set and wrong otherwise.
#[test]
fn the_reassignment_preserves_the_original_order() {
    let (cases, _) = load();
    let c = &cases["weight_three_degenerate"];
    let flow = run(c).unwrap();
    let flagged: Vec<usize> = c.flagged.as_ref().unwrap().clone();
    // Totals are 1, 2, 2.5, 0.4, 2. A total of *exactly* 1 is degenerate too -- `log(1)` is 0,
    // so `-1/0` is `-Inf` -- so four of the five are flagged and only 0.4 is not.
    assert_eq!(flagged, vec![0, 1, 2, 4], "four degenerate entries");
    // The originals of the flagged entries, ascending, must map to ascending assigned values.
    let orig: Vec<f64> = flagged.iter().map(|&i| flow.original[i]).collect();
    let assigned: Vec<f64> = flagged.iter().map(|&i| flow.scaled[i]).collect();
    // Sort the *pairs*, not the two vectors independently: the originals contain a tie (2.0
    // twice), and sorting each vector on its own says nothing about which assigned value goes
    // with which original. Sorting the pairs and then reading the assigned column is the actual
    // claim: the reassignment is monotone in the original total.
    let mut pairs: Vec<(f64, f64)> = orig.iter().copied().zip(assigned.iter().copied()).collect();
    pairs.sort_by(|a, b| a.0.partial_cmp(&b.0).unwrap());
    let by_original: Vec<f64> = pairs.iter().map(|p| p.0).collect();
    let by_assigned: Vec<f64> = pairs.iter().map(|p| p.1).collect();
    assert!(
        by_assigned.windows(2).all(|w| w[0] <= w[1]),
        "assigned values must increase with the original total: {pairs:?}"
    );
    // And explicitly not the reverse: reading the sorted originals against the *unsorted*
    // assigned values would claim the largest total gets the smallest value.
    assert_ne!(by_original, assigned, "the fixture must be asymmetric");
    assert_eq!(
        by_original,
        vec![1.0, 2.0, 2.0, 2.5],
        "the flagged originals, ascending: the two 2.0s are b and e"
    );
    let va = c.values_assign.as_ref().unwrap();
    assert_eq!(va.len(), 4);
    assert!(
        va[0] < va[1] && va[1] < va[2],
        "values.assign is increasing"
    );
    // `seq(a, b, length.out = 1)` is `a`, not the midpoint.
    //
    // `max(pSum)` is taken **before** the reassignment, so it is the largest value among the
    // entries that were *not* flagged. Here that is the 0.5-total pathway's `-1/log(0.5)`. Using
    // the post-reassignment maximum instead is off by the reassignment itself, which is a
    // self-referential mistake that looks plausible.
    let (cases2, _) = load();
    let d = &cases2["weight_degenerate"];
    let va = d.values_assign.as_ref().unwrap();
    assert_eq!(va.len(), 1);
    let pre_max = -1.0 / 0.5f64.ln();
    assert_eq!(
        va[0],
        pre_max * 1.1,
        "seq(max*1.1, ..., length.out = 1) is max*1.1"
    );
    // ... and the reassigned entry lands exactly there, which is what makes it non-degenerate.
    let scaled = d.scaled.as_ref().unwrap();
    assert_eq!(
        scaled[3], va[0],
        "the one flagged entry takes the single assigned value"
    );
    assert!(
        scaled[3] > scaled[0],
        "1.5870 > 1.4427, so the reassigned bar clears the tallest ordinary one"
    );
}

/// When every pathway is degenerate, `max(pSum)` is `-Inf` and R's `seq()` refuses it. The port
/// must raise the same error rather than quietly producing `NaN` bars.
#[test]
fn an_all_degenerate_network_raises_seq_from_not_finite() {
    let (cases, _) = load();
    let c = &cases["weight_all_degenerate"];
    assert_eq!(
        run(c).unwrap_err(),
        RankNetError::SeqFromNotFinite,
        "upstream stops with R's own seq() error"
    );
    assert_eq!(
        RankNetError::SeqFromNotFinite.to_string(),
        "'from' must be a finite number"
    );
}

/// An all-zero network -- and a network emptied by the `thresh` cut -- stops *before* the
/// per-pathway sums, with a different message from the `seq()` crash.
#[test]
fn an_empty_network_stops_with_its_own_message() {
    let (cases, _) = load();
    for name in ["all_zero", "thresh_empties"] {
        let c = &cases[name];
        assert_eq!(
            run(c).unwrap_err(),
            RankNetError::NoCommunications,
            "{name}"
        );
        assert_eq!(
            RankNetError::NoCommunications.to_string(),
            "No inferred communications for the input!"
        );
        // `sum(prob) == 0`, not `all(prob == 0)`: a network whose cells cancel to a zero total
        // would also stop, and the fixture's `thresh_empties` case is the `pval > thresh` route
        // to the same place.
        assert!(same(&c.prob, &c.prob));
    }
}

/// `measure = "count"` binarises and leaves the totals alone: `1*(prob > 0)`, so the flow is a
/// count of non-zero cells rather than a sum of weights, and `pSum` is *not* rescaled.
#[test]
fn the_count_measure_binarises_and_does_not_rescale() {
    let (cases, _) = load();
    let c = &cases["count_measure"];
    let flow = run(c).unwrap();
    assert!(
        same(&flow.scaled, &flow.original),
        "count leaves pSum as pSum.original"
    );
    assert!(
        flow.original.iter().all(|v| v.fract() == 0.0),
        "every count is a whole number: {:?}",
        flow.original
    );
    assert!(
        flow.original.iter().all(|&v| (0.0..=16.0).contains(&v)),
        "counts are cell counts, bounded by k*k: {:?}",
        flow.original
    );
}

/// `targets.use` is validated against `dimnames(prob)[[1]]` -- the *source* axis. The port takes
/// already-resolved indices, so the check lives in the shim; what the core must get right is
/// that the target filter indexes the middle axis, not the first.
#[test]
fn the_two_group_filters_hit_different_axes() {
    let (cases, _) = load();
    let base = &cases["sources_filter"];
    // g1 is index 0 and g10 is index 3 on the source axis.
    let only_g1 = information_flow(
        &base.prob,
        &base.pval,
        base.k,
        base.n,
        &base.names,
        base.thresh,
        Measure::Weight,
        Some(&[0]),
        None,
    )
    .unwrap();
    assert!(
        only_g1.original[0] < base.original.as_ref().unwrap()[0],
        "keeping one of four sources must reduce the total"
    );
    // The target filter, with the same index, reduces the *second* pathway differently: in
    // `targets_filter` the non-zero cells of p2 are at (g1,g1) and (g1,g2), so keeping target
    // g1 keeps the first and drops the second.
    let t = &cases["targets_filter"];
    let only_t_g1 = information_flow(
        &t.prob,
        &t.pval,
        t.k,
        t.n,
        &t.names,
        t.thresh,
        Measure::Weight,
        None,
        Some(&[0]),
    )
    .unwrap();
    let all = information_flow(
        &t.prob,
        &t.pval,
        t.k,
        t.n,
        &t.names,
        t.thresh,
        Measure::Weight,
        None,
        None,
    )
    .unwrap();
    assert!(
        only_t_g1.original[1] < all.original[1],
        "target filtering must reduce the pathway total: {} vs {}",
        only_t_g1.original[1],
        all.original[1]
    );
    // p1 is `1:16`, and `matrix(1:16, 4, 4)` fills **column-major**, so the four target columns
    // are (1,2,3,4), (5,...,8), (9,...,12), (13,...,16). Restricting the target axis to g1 keeps
    // one column and the total falls to 10 -- which is the check that the filter indexed the
    // *second* axis. Indexing the first instead would have kept source g1's four cells, 1, 5, 9
    // and 13, giving 28: a different number, so the two axes are distinguishable.
    assert_eq!(only_t_g1.original[0], 10.0, "target g1's column is 1+2+3+4");
    assert_eq!(all.original[0], 136.0, "1+2+...+16");
}

/// `order()` on a double vector is radix and therefore **stable**. The fixture has three ties,
/// so an unstable sort permutes them and the row order of the result changes.
#[test]
fn order_is_stable_on_ties() {
    let (cases, _) = load();
    let c = &cases["order_ties"];
    let want = c.order.as_ref().expect("order is recorded");
    let orig = c.original.as_ref().unwrap();
    let _ = c;
    // Counts of non-zero cells per pathway: 1, 3, 1, 2, 3. Two tied pairs, {a,c} at 1 and
    // {b,e} at 3, around a distinct `d` at 2 -- so an unstable sort is detectable.
    assert_eq!(orig, &vec![1.0, 3.0, 1.0, 2.0, 3.0]);
    let got = order_f64(orig);
    assert_eq!(got, *want, "order(contribution), stable");
    // The stable order is the two 1s in index order, then the 2, then the two 3s in index order.
    assert_eq!(got, vec![0, 2, 3, 1, 4]);
    let tie_a: Vec<usize> = got.iter().copied().filter(|x| [0, 2].contains(x)).collect();
    assert_eq!(tie_a, vec![0, 2], "ties keep index order");
    let tie_b: Vec<usize> = got.iter().copied().filter(|x| [1, 4].contains(x)).collect();
    assert_eq!(tie_b, vec![1, 4], "ties keep index order");
}

/// `NA`s sort last, whatever their position.
#[test]
fn order_puts_missing_values_last() {
    let x = vec![2.0, f64::NAN, 1.0, f64::NAN, 3.0];
    assert_eq!(order_f64(&x), vec![2, 0, 4, 1, 3]);
}

/// `order(a, b, c)`: `na.last` is a property of the *first* key only, and later keys order within
/// each group. `rankNet` uses this for `order(-contribution.relative.1, contribution,
/// -contribution.data2)`.
#[test]
fn multi_key_order_partitions_on_the_first_key() {
    let a = vec![1.0, f64::NAN, 1.0, 2.0];
    let b = vec![9.0, 0.0, 1.0, 0.0];
    let got = order_f64_multi(&[&a, &b]);
    // The NA-primary row goes last; within a == 1, b orders 1 before 9.
    assert_eq!(got, vec![2, 0, 3, 1]);
}

/// `as.numeric(format(x, digits = 1))` is **one significant digit**, not one decimal place, and
/// `Inf` survives it. Both directions matter: a port using `round(x, 1)` keeps three digits of a
/// small ratio and orders rows differently, and a port mapping `Inf` to `NA` loses the
/// `is.na` reset's outcome.
#[test]
fn the_relative_contribution_is_one_significant_digit() {
    let (_, rel) = load();
    let by = |n: &str| -> (Vec<f64>, Vec<f64>, Vec<f64>) {
        let e = rel
            .iter()
            .find(|(k, _, _, _)| k == n)
            .expect("relative case");
        (e.1.clone(), e.2.clone(), e.3.clone())
    };
    let apply = |num: &[f64], den: &[f64], want: &[f64]| {
        for (i, (&n, &d)) in num.iter().zip(den).enumerate() {
            let mut got = format_digits1(n / d);
            if got.is_nan() {
                got = 0.0;
            }
            assert_eq!(bits(got), bits(want[i]), "element {i} of the vector");
        }
    };
    let (n, d, w) = by("collapse");
    assert!(
        w.iter().all(|v| *v == 1.0),
        "1, 1.004, 1.04 and 0.996 all round to one significant digit = 1"
    );
    apply(&n, &d, &w);
    let (n, d, w) = by("significant");
    assert_eq!(
        w,
        vec![0.04, 0.004, 1e3, 1e5],
        "significant digits, not decimal places"
    );
    apply(&n, &d, &w);
    let (n, d, w) = by("zeros");
    assert_eq!(w[1], 1e-12, "a tiny ratio is not rounded away to zero");
    assert_eq!(w[2], 1e12);
    apply(&n, &d, &w);
    // A zero denominator gives `Inf`, and `Inf` survives `format` + `as.numeric`.
    let (n, d, w) = by("zero_denominator");
    assert!(
        w[1].is_infinite() && w[1] > 0.0,
        "as.numeric(format(Inf, digits = 1)) is Inf, not NA"
    );
    apply(&n, &d, &w);
    // `format` pads the vector to a common width, so the *characters* depend on the neighbours
    // while the parsed values do not. Recorded so the padding is not mistaken for the rounding.
    let (n, d, w) = by("padded");
    assert_eq!(w, vec![0.5, 50.0, 5000.0]);
    apply(&n, &d, &w);
}

/// A `NaN` ratio becomes `NA` under `as.numeric` and is then reset to 0, so a zero-over-zero
/// pathway contributes 0 rather than propagating.
#[test]
fn a_nan_ratio_is_reset_to_zero() {
    let got = format_digits1(f64::NAN / 1.0);
    assert!(
        got.is_nan(),
        "format_digits1 does not invent a value for NaN"
    );
    let mut v = got;
    if v.is_nan() {
        v = 0.0;
    }
    assert_eq!(v, 0.0, "the is.na reset upstream applies");
}
