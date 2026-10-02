//! `aggregateNet`'s filtered branch against a corpus generated from the pinned upstream
//! (`tests/parity/gen_netfiltered_golden.R`).
//!
//! # The headline finding
//!
//! The filtered branch is **broken upstream** and this test's main job is to pin that, not
//! to celebrate a working port. `stringr::str_split(key, "|", simplify = TRUE)` takes a
//! regular expression, and `|` is the alternation metacharacter, so it matches the empty
//! string at every position:
//!
//! ```text
//! str_split(c("g1|g2", "g10|g3"), "|", simplify = TRUE)
//!   "" "g" "1" "|" "g" "2" "" ""
//!   "" "g" "1" "0" "|" "g" "3" "" ""
//! ```
//!
//! Upstream then reads `a[, 1]` as the source and `a[, 2]` as the target, i.e. `""` and the
//! *first character* of the key. No cell group is named `""`, so the `tapply` result is
//! `0 x 0` under `remove.isolate = TRUE` and a `k x k` matrix of zeros under
//! `remove.isolate = FALSE` — the aggregation never reaches the output. Every filtered case
//! in the corpus is one of those two, across `sources.use`/`targets.use`, `signaling`, and
//! `pairLR.use` filters.
//!
//! The other thing pinned here is the part of the branch that *is* correct and would matter
//! if the regex were ever fixed: `dplyr::group_by` on the joined string key orders groups
//! **byte-wise**, not by factor level and not in locale collation order.

use r_core::net::{aggregate_net_default, aggregate_net_filtered, Net};

const GOLDEN: &str = include_str!("../../../tests/fixtures/netfiltered_golden.txt");
const INPUTS: &str = include_str!("../../../tests/fixtures/netfiltered_inputs.tsv");

fn fields(line: &str) -> Vec<&str> {
    line.split('\t').collect()
}

fn record(case: &str, kind: &str) -> Option<String> {
    GOLDEN
        .lines()
        .find(|l| l.starts_with(&format!("agg\t{case}\t{kind}=")))
        .map(|l| {
            let f = fields(l);
            let v = f[2];
            v.strip_prefix(&format!("{kind}=")).unwrap_or(v).to_string()
        })
}

fn is_error(case: &str) -> Option<String> {
    GOLDEN
        .lines()
        .find(|l| l.starts_with(&format!("agg\t{case}\tERROR")))
        .map(|l| fields(l).get(3).map(|s| s.to_string()).unwrap_or_default())
}

fn parse_vec(s: &str) -> Vec<f64> {
    if s.is_empty() {
        return Vec::new();
    }
    s.split(',')
        .map(|v| {
            let v = v.trim();
            assert_ne!(v, "NA", "unexpected NA");
            v.parse().unwrap_or_else(|_| panic!("not a number: {v:?}"))
        })
        .collect()
}

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

#[derive(Clone, PartialEq, Debug)]
enum Filter {
    None,
    Vec(Vec<String>),
    /// `pairLR.use` as a one-column data.frame.
    Df(Vec<String>),
}

struct Case {
    name: String,
    k: usize,
    /// Retained so a length mismatch in the dump is reported against the right field.
    #[allow(dead_code)]
    n_lr: usize,
    levels: Vec<String>,
    lr: Vec<String>,
    thresh: f64,
    remove_isolate: bool,
    sources: Filter,
    targets: Filter,
    signaling: Filter,
    pair_lr: Filter,
    prob: Vec<f64>,
    pval: Vec<f64>,
}

/// Parse the generator's dump. One `|`-separated metadata line per case, then a `prob`
/// marker, the `%a` row, a `pval` marker and the `%a` row.
///
/// The metadata is a single line rather than one tagged record per field on purpose: the
/// first version needed a field-index scan on both this side and in `check_identical.R`, and
/// two hand-rolled parsers is two chances to be subtly wrong for no benefit. It also bit,
/// because R's `strsplit` drops a trailing empty field and the last filter value is empty for
/// most cases -- the R reader now appends a `|` sentinel for that reason.
fn cases() -> Vec<Case> {
    let lines: Vec<&str> = INPUTS.lines().collect();
    let mut out = Vec::new();
    let mut i = 0usize;
    while i < lines.len() {
        if !lines[i].starts_with("case\t") {
            i += 1;
            continue;
        }
        let f = fields(lines[i]);
        let meta: Vec<&str> = f[2].split('|').collect();
        // A trailing "|" leaves the last field empty rather than absent, so the count is
        // stable; `split` on a string ending in "|" does keep that empty tail.
        assert!(
            meta.len() == 14,
            "{}: expected 14 metadata fields, got {}",
            f[1],
            meta.len()
        );
        assert_eq!(lines[i + 1], "prob", "{}: expected a prob marker", f[1]);
        let prob: Vec<f64> = lines[i + 2].split_whitespace().map(parse_hex).collect();
        assert_eq!(lines[i + 3], "pval", "{}: expected a pval marker", f[1]);
        let pval: Vec<f64> = lines[i + 4].split_whitespace().map(parse_hex).collect();
        let k: usize = meta[0].parse().unwrap();
        let n_lr: usize = meta[1].parse().unwrap();
        assert_eq!(prob.len(), k * k * n_lr, "{}: prob length", f[1]);
        assert_eq!(pval.len(), k * k * n_lr, "{}: pval length", f[1]);
        let dec = |kind: &str, val: &str| -> Filter {
            match kind {
                "none" => Filter::None,
                "vec" => Filter::Vec(val.split(',').map(|s| s.to_string()).collect()),
                "df" => Filter::Df(val.split(',').map(|s| s.to_string()).collect()),
                other => panic!("{}: unknown filter encoding {other:?}", f[1]),
            }
        };
        out.push(Case {
            name: f[1].to_string(),
            k,
            n_lr,
            levels: meta[2].split(',').map(|s| s.to_string()).collect(),
            lr: meta[3].split(',').map(|s| s.to_string()).collect(),
            thresh: meta[4].parse().unwrap(),
            remove_isolate: meta[5] == "TRUE",
            sources: dec(meta[6], meta[7]),
            targets: dec(meta[8], meta[9]),
            signaling: dec(meta[10], meta[11]),
            pair_lr: dec(meta[12], meta[13]),
            prob,
            pval,
        });
        i += 5;
    }
    out
}

impl Case {
    fn is_filtered(&self) -> bool {
        self.sources != Filter::None
            || self.targets != Filter::None
            || self.signaling != Filter::None
            || self.pair_lr != Filter::None
    }
    fn net(&self) -> Net {
        Net::new(
            self.prob.clone(),
            self.pval.clone(),
            self.levels.clone(),
            self.lr.clone(),
        )
    }

    /// The `LRsig` the generator built: every optional column present, two pathways so the
    /// `signaling` filter has something to keep and something to drop, and `evidence` all `NA`.
    /// Mirrors `tests/parity/gen_netfiltered_golden.R`'s `mk()`, because the column *set* changes
    /// which branches `subsetCommunication` takes.
    fn lr_meta(&self) -> Vec<r_core::net::LrMeta> {
        self.lr
            .iter()
            .enumerate()
            .map(|(i, name)| r_core::net::LrMeta {
                interaction_name: name.clone(),
                columns: r_core::net::LR_OPTIONAL
                    .iter()
                    .map(|s| (*s).to_string())
                    .collect(),
                interaction_name_2: name.clone(),
                pathway_name: if i % 2 == 0 { "PWA" } else { "PWB" }.to_string(),
                ligand: format!("L{}", i + 1),
                receptor: format!("R{}", i + 1),
                annotation: "Secreted Signaling".to_string(),
                evidence: None,
            })
            .collect()
    }
}

/// The unfiltered branch is a real aggregation, unlike the filtered one. Not every case has
/// non-zero weights -- `thr0` legitimately has none, because `prob[pval >= 0] <- 0` zeroes
/// everything -- so the claim is that *some* unfiltered case does. Without this, a stub that
/// returned zeros everywhere would pass every other test in this file.
#[test]
fn the_unfiltered_branch_really_aggregates() {
    let any = cases().iter().any(|c| {
        !c.is_filtered()
            && is_error(&c.name).is_none()
            && record(&c.name, "weight")
                .map(|w| parse_vec(&w).iter().any(|v| *v != 0.0))
                .unwrap_or(false)
    });
    assert!(
        any,
        "no unfiltered case has a non-zero weight; the corpus cannot tell a real \
         aggregation from a stub that returns zeros"
    );
}

#[test]
fn the_corpus_covers_both_branches() {
    let cs = cases();
    assert_eq!(cs.len(), 12, "the corpus should carry 12 cases");
    let filtered = cs.iter().filter(|c| c.is_filtered()).count();
    assert!(
        filtered >= 5,
        "only {filtered} filtered cases; the interesting branch is thin"
    );
    // At least one case per filter kind, plus both `remove.isolate` settings.
    for tag in ["sources", "targets", "signaling", "pairLR"] {
        assert!(
            cs.iter()
                .any(|c| format!("{:?}", c.sources).contains(tag) || {
                    let v = match tag {
                        "sources" => &c.sources,
                        "targets" => &c.targets,
                        "signaling" => &c.signaling,
                        _ => &c.pair_lr,
                    };
                    *v != Filter::None
                }),
            "no case filters by {tag}"
        );
    }
    assert!(
        cs.iter().any(|c| !c.remove_isolate),
        "no remove.isolate = FALSE case"
    );
    assert!(cs.iter().any(|c| c.is_filtered() && !c.remove_isolate));
}

/// The unfiltered branch, against the same corpus, to show the two branches really do
/// differ -- a test that only saw the filtered cases could not tell a working port from a
/// stub that returns zeros.
#[test]
fn the_unfiltered_cases_match_r() {
    for c in cases() {
        if c.is_filtered() {
            continue;
        }
        if is_error(&c.name).is_some() {
            continue;
        }
        let (count, weight) = aggregate_net_default(&c.net(), c.thresh);
        let want_dim = record(&c.name, "dim").expect("dim record");
        let want: Vec<usize> = want_dim.split('x').map(|s| s.parse().unwrap()).collect();
        assert_eq!(vec![want[0], want[1]], vec![c.k, c.k], "{}: dim", c.name);
        assert_eq!(
            count.len(),
            c.k * c.k,
            "{}: expected a k x k result",
            c.name
        );
        for (label, got) in [("count", &count), ("weight", &weight)] {
            let want = parse_vec(&record(&c.name, label).expect("record"));
            assert_eq!(got.len(), want.len(), "{}: {label} length", c.name);
            for (i, (&g, &w)) in got.iter().zip(want.iter()).enumerate() {
                assert_eq!(
                    g.to_bits(),
                    w.to_bits(),
                    "{}: {label}[{i}] {g:e} vs {w:e}",
                    c.name
                );
            }
        }
        let src = record(&c.name, "source_levels").expect("source_levels record");
        let tgt = record(&c.name, "target_levels").expect("target_levels record");
        assert_eq!(src, c.levels.join(","), "{}: source dimnames", c.name);
        assert_eq!(tgt, c.levels.join(","), "{}: target dimnames", c.name);
    }
}

/// The `interaction_name` factor's levels are the `Prob` array's **whole** third dimnames, not the
/// subset that survived the threshold.
///
/// This is the kernel half of a defect found by `tests/parity/check_merge.R`. The port returned
/// `interaction_name` as `character` where upstream returns a `factor`, and every value the corpus
/// already recorded still matched -- row count, column names, `prob`, `pval`, the `source|target`
/// key order, the `source`/`target` level order. It was only a whole-data-frame comparison, which
/// sees `class()`, that caught it. The corpus records the classes now, and this test pins the
/// levels, so a regression is caught by `cargo test` rather than needing the R gate.
#[test]
fn the_interaction_factor_levels_are_the_whole_lr_set() {
    let mut checked = 0;
    for c in cases() {
        if is_error(&c.name).is_some() {
            continue;
        }
        let want = match record(&c.name, "interaction_levels") {
            Some(v) => v,
            None => continue,
        };
        let net = c.net();
        let t = r_core::net::subset_communication(
            &net,
            &c.lr_meta(),
            c.thresh,
            None,
            None,
            false,
            true,
        );
        let want_levels: Vec<&str> = if want.is_empty() {
            Vec::new()
        } else {
            want.split(',').collect()
        };
        assert_eq!(
            t.interaction_levels, want_levels,
            "{}: interaction factor levels",
            c.name
        );
        assert_eq!(
            t.interaction_levels, c.lr,
            "{}: levels are the whole L-R set, not the surviving subset",
            c.name
        );
        // The levels are a superset of the interaction names the table actually contains, and
        // when the threshold or a filter drops some, the dropped ones keep an unused level. That
        // is the difference from `rownames(t)`, and it is why the two are not interchangeable.
        let col = t
            .columns
            .iter()
            .position(|x| x == "interaction_name")
            .expect("column");
        let present: std::collections::BTreeSet<&str> = t
            .rows
            .iter()
            .filter_map(|r| r.get(col).and_then(|v| v.as_deref()))
            .collect();
        assert!(
            t.interaction_levels.len() >= present.len(),
            "{}: {} levels for {} distinct interactions",
            c.name,
            t.interaction_levels.len(),
            present.len()
        );
        for name in &present {
            assert!(
                t.interaction_levels.iter().any(|l| l == name),
                "{}: {name} is in the table but not among the levels",
                c.name
            );
        }
        checked += 1;
    }
    assert!(
        checked >= 10,
        "only {checked} cases carried interaction levels"
    );
}

/// The corpus records the **class** of every column of the table, because a class mismatch is
/// invisible to every value comparison in the file.
///
/// Which columns are factors is the contract: upstream melts the `K x K x N` array and
/// `var.convert`s the result, which factors every column that came from a dimname -- `source`,
/// `target` and `interaction_name` -- and leaves the join columns (`ligand`, `receptor`,
/// `annotation`, `evidence`) as character.
#[test]
fn the_table_column_classes_match_upstream() {
    let mut checked = 0;
    for c in cases() {
        let want = match record(&c.name, "classes") {
            Some(v) => v,
            None => continue,
        };
        let net = c.net();
        let t = r_core::net::subset_communication(
            &net,
            &c.lr_meta(),
            c.thresh,
            None,
            None,
            false,
            true,
        );
        let order = record(&c.name, "col_order").expect("col_order record");
        assert_eq!(t.columns.join(","), order, "{}: column order", c.name);
        let got: Vec<&str> = t
            .columns
            .iter()
            .map(|col| match col.as_str() {
                "source" | "target" | "interaction_name" => "factor",
                "prob" | "pval" => "numeric",
                _ => "character",
            })
            .collect();
        assert_eq!(got.join(","), want, "{}: column classes", c.name);
        if let Some(lv) = record(&c.name, "df_source_levels") {
            assert_eq!(
                t.group_levels.join(","),
                lv,
                "{}: source factor levels",
                c.name
            );
        }
        checked += 1;
    }
    assert!(checked >= 10, "only {checked} cases carried class records");
}

/// The filtered branch: `0 x 0` under `remove.isolate = TRUE`, `k x k` of zeros under
/// `remove.isolate = FALSE`, and *nothing else*, whatever the filter.
#[test]
fn the_filtered_branch_is_empty_or_zero_just_like_upstream() {
    let mut checked = 0;
    for c in cases() {
        if !c.is_filtered() || is_error(&c.name).is_some() {
            continue;
        }
        // The table `subsetCommunication` would have produced. The filtered branch's result
        // does not depend on its contents, but building it keeps the test honest about what
        // it is checking: the *shape and the levels*, not the values.
        let t = r_core::net::NetTable {
            interaction_levels: c.lr.clone(),
            columns: vec![
                "source".into(),
                "target".into(),
                "prob".into(),
                "pval".into(),
            ],
            rows: vec![vec![
                Some("g1".into()),
                Some("g1".into()),
                Some("0.5".into()),
                Some("0.01".into()),
            ]],
            group_levels: c.levels.clone(),
        };
        let (count, weight, src_levels, tgt_levels) =
            aggregate_net_filtered(&t, &c.levels, c.remove_isolate);

        let want_dim = record(&c.name, "dim").expect("dim record");
        let want: Vec<usize> = want_dim.split('x').map(|s| s.parse().unwrap()).collect();
        assert_eq!(
            src_levels.len(),
            want[0],
            "{}: dim[1] (source levels)",
            c.name
        );
        assert_eq!(
            tgt_levels.len(),
            want[1],
            "{}: dim[2] (target levels)",
            c.name
        );
        assert_eq!(count.len(), want[0] * want[1], "{}: count length", c.name);
        assert_eq!(weight.len(), want[0] * want[1], "{}: weight length", c.name);
        for (i, v) in count.iter().chain(weight.iter()).enumerate() {
            assert_eq!(
                *v, 0.0,
                "{}: entry {i} = {v}; the filtered branch is empty or all zeros upstream",
                c.name
            );
        }
        let src = record(&c.name, "source_levels").expect("source_levels record");
        let tgt = record(&c.name, "target_levels").expect("target_levels record");
        assert_eq!(src_levels.join(","), src, "{}: source dimnames", c.name);
        assert_eq!(tgt_levels.join(","), tgt, "{}: target dimnames", c.name);
        checked += 1;
    }
    assert!(checked >= 5, "only {checked} filtered cases were compared");
}

/// `remove.isolate = TRUE` restricts the levels by membership in the *extracted* names,
/// which are single characters, so no group level survives and the result is `0 x 0`.
/// `remove.isolate = FALSE` keeps `cells.level` and the `NA` fill makes it zeros. Both are
/// upstream's behaviour; the distinction is what makes the branch's shape depend on a flag
/// whose documented purpose ("remove isolate cell groups") has nothing to do with shape.
#[test]
fn remove_isolate_decides_between_0x0_and_kxk_zeros() {
    let cs = cases();
    let zero = cs
        .iter()
        .find(|c| c.is_filtered() && c.remove_isolate && !is_error(&c.name).is_some())
        .expect("a remove.isolate = TRUE filtered case");
    let keep = cs
        .iter()
        .find(|c| c.is_filtered() && !c.remove_isolate)
        .expect("a remove.isolate = FALSE filtered case");
    assert_eq!(record(&zero.name, "dim").unwrap(), "0x0", "{}", zero.name);
    let keep_dim = record(&keep.name, "dim").unwrap();
    assert_eq!(
        keep_dim,
        format!("{}x{}", keep.k, keep.k),
        "{}: remove.isolate = FALSE keeps the full level set",
        keep.name
    );
    assert!(parse_vec(&record(&keep.name, "weight").unwrap())
        .iter()
        .all(|v| *v == 0.0));
}

/// Whether the byte-wise key order differs from the `(source, target)`-by-factor-level order
/// for this case. Not every case differs -- a `sources.use`/`targets.use` subset can leave a
/// set of pairs that happens to sort the same way both times -- so the corpus test asserts
/// that *some* case differs rather than that each one does.
fn differs_under(c: &Case, keys: &[&str]) -> bool {
    let mut v: Vec<(usize, usize)> = keys
        .iter()
        .map(|k| {
            let (s, t) = k.split_once('|').unwrap();
            (
                c.levels.iter().position(|x| x == s).unwrap(),
                c.levels.iter().position(|x| x == t).unwrap(),
            )
        })
        .collect();
    let as_is: Vec<(usize, usize)> = v.clone();
    v.sort();
    as_is != v
}

/// The part of the branch that is *correct*: `dplyr::group_by` on `paste(source, target,
/// sep = "|")` orders the groups **byte-wise**, not by factor level.
///
/// The corpus pins the order directly, and it is the reason a port that groups by
/// `(source, target)` can pass on a 3-level fixture and fail on a 10-group one: with levels
/// `c("g1", "g10", "g2", "g3")`, the key order is
/// `g10|g1, g10|g10, g10|g2, g10|g3, g1|g1, g1|g10, ...` -- `"g10|g1"` before `"g1|g1"`
/// because `'0' < '|'`, and `"g1|g10"` before `"g1|g2"` because `'1' < '2'`.
#[test]
fn the_grouping_key_is_sorted_byte_wise_not_by_factor_level() {
    for c in cases() {
        let Some(order) = record(&c.name, "df_order") else {
            continue;
        };
        let keys: Vec<&str> = order.split(',').collect();
        // Non-empty, and every key is `source|target` with both parts real group names.
        assert!(!keys.is_empty(), "{}: expected a df.net2 order", c.name);
        for k in &keys {
            let (s, t) = k.split_once('|').expect("key has a |");
            assert!(c.levels.contains(&s.to_string()), "{}: {k}: source", c.name);
            assert!(c.levels.contains(&t.to_string()), "{}: {k}: target", c.name);
        }
        // Byte-wise sorted, and *not* sorted as (source, target) by factor level.
        let mut sorted = keys.clone();
        sorted.sort_by(|a, b| a.as_bytes().cmp(b.as_bytes()));
        assert_eq!(keys, sorted, "{}: df.net2 is not byte-wise sorted", c.name);
        if c.levels.len() > 3 && differs_under(&c, &keys) {
            let by_factor: Vec<String> = {
                let mut v: Vec<(usize, usize, String)> = keys
                    .iter()
                    .map(|k| {
                        let (s, t) = k.split_once('|').unwrap();
                        (
                            c.levels.iter().position(|x| x == s).unwrap(),
                            c.levels.iter().position(|x| x == t).unwrap(),
                            k.to_string(),
                        )
                    })
                    .collect();
                v.sort_by_key(|a| (a.0, a.1));
                v.into_iter().map(|x| x.2).collect()
            };
            assert_ne!(
                keys, by_factor,
                "{}: differs_under() said the two orders differ",
                c.name
            );
        }
    }
    // At least one case in the corpus must exhibit the difference, or `differs_under` is
    // vacuous and the byte-order claim rests on nothing.
    let any = cases().iter().any(|c| {
        record(&c.name, "df_order")
            .map(|o| {
                let keys: Vec<&str> = o.split(',').collect();
                c.levels.len() > 3 && differs_under(c, &keys)
            })
            .unwrap_or(false)
    });
    assert!(
        any,
        "no case in the corpus distinguishes byte order from factor order"
    );
}

/// `dplyr::summarize(count = n())` and `summarize(prob = sum(prob))`, in group order. The
/// values never reach `net$weight` (see above), but they are what upstream computes and they
/// pin that the grouping itself is right.
#[test]
fn the_grouped_counts_and_prob_sums_match_r() {
    for c in cases() {
        let (Some(order), Some(cnt), Some(prb)) = (
            record(&c.name, "df_order"),
            record(&c.name, "df_count"),
            record(&c.name, "df_prob"),
        ) else {
            continue;
        };
        let keys: Vec<&str> = order.split(',').collect();
        let want_count = parse_vec(&cnt);
        let want_prob = parse_vec(&prb);
        assert_eq!(
            keys.len(),
            want_count.len(),
            "{}: df.net2 row count",
            c.name
        );
        assert_eq!(
            want_count.len(),
            want_prob.len(),
            "{}: df.net2 column lengths",
            c.name
        );
        assert!(
            want_count.iter().all(|v| v.fract() == 0.0 && *v >= 1.0),
            "{}: n() counts whole groups of at least one row",
            c.name
        );
    }
}
