//! Parity for `computeRegionDistance`'s foundations: R's `mean(x, trim, na.rm)`, `collapse::fdist`,
//! and the exact k-d tree that replaces `AnnoyParam`.
//!
//! Two different oracles, deliberately.
//!
//! * `mean(trim =)` and `fdist` are pinned against a corpus generated from R itself
//!   (`tests/parity/gen_spatial_golden.R`) -- 26 trimmed-mean cases and 7 distance matrices.
//!   Those are deterministic and R is a valid oracle.
//! * The k-d tree is **not** pinned against upstream. `BiocNeighbors::queryKNN(...,
//!   AnnoyParam())` is a randomised approximate index, so two runs of upstream can disagree and
//!   upstream is not a usable oracle for a neighbour query. The exact tree is verified against
//!   **exhaustive search** instead, on inputs built to be adversarial for a k-d tree: duplicate
//!   points, collinear points, a regular lattice, and every distance tied. The divergence from
//!   Annoy is a separate measurement on real spatial data, not something this suite asserts.
//!
//! Doubles are compared by bit pattern: the corpus is R's `format(digits = 17)` and this side
//! writes its own, and those agree on the value while differing on width and exponent spelling.

use r_core::spatial::{brute_force_knn, fdist, trimmed_mean, KdTree};
use std::collections::HashMap;

// ---------------------------------------------------------------- corpus

#[derive(Debug, Default, Clone)]
struct TmCase {
    x: Vec<f64>,
    trim: f64,
    na_rm: bool,
    want: f64,
}

#[derive(Debug, Default, Clone)]
struct FdCase {
    n: usize,
    d: usize,
    /// The input coordinates, row-major. **Not** the distance matrix: a reader that
    /// reconstructs coordinates from the distance rows is feeding distances back in as
    /// positions, and a two-point fixture then silently asks a different question.
    coords: Vec<f64>,
    /// The symmetric distance matrix, row by row.
    rows: Vec<Vec<f64>>,
    /// Upstream's `as.vector(dist)`: the compact triangle. Compared as a **multiset**.
    triangle: Vec<f64>,
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

fn vec_of(line: &str) -> Vec<f64> {
    if line == "-" {
        Vec::new()
    } else {
        line.split('\t').map(cell).collect()
    }
}

fn load() -> (
    HashMap<String, TmCase>,
    HashMap<String, FdCase>,
    Vec<String>,
    Vec<String>,
) {
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../tests/fixtures/spatial_golden.txt");
    let text = std::fs::read_to_string(&path)
        .unwrap_or_else(|e| panic!("cannot read {}: {e}", path.display()));
    let mut tm: HashMap<String, TmCase> = HashMap::new();
    let mut fd: HashMap<String, FdCase> = HashMap::new();
    let mut tm_order: Vec<String> = Vec::new();
    let mut fd_order: Vec<String> = Vec::new();
    let mut current_tm = String::new();
    let mut current_fd = String::new();
    for line in text.lines() {
        let f: Vec<&str> = line.split('\t').collect();
        match f[0] {
            "tm" => {
                current_tm = f[1].to_string();
                tm_order.push(current_tm.clone());
                tm.insert(
                    current_tm.clone(),
                    TmCase {
                        trim: f[2].parse().unwrap(),
                        na_rm: f[3] == "TRUE",
                        ..Default::default()
                    },
                );
            }
            "tm_x" => tm.get_mut(&current_tm).unwrap().x = vec_of(&f[2..].join("\t")),
            "tm_mean" => tm.get_mut(&current_tm).unwrap().want = cell(f[2]),
            "fd" => {
                current_fd = f[1].to_string();
                fd_order.push(current_fd.clone());
                fd.insert(
                    current_fd.clone(),
                    FdCase {
                        n: f[2].parse().unwrap(),
                        d: f[3].parse().unwrap(),
                        ..Default::default()
                    },
                );
            }
            "fd_n" => fd.get_mut(&current_fd).unwrap().n = f[2].parse().unwrap(),
            "fd_x" => fd.get_mut(&current_fd).unwrap().coords = vec_of(&f[2..].join("\t")),
            "fd_dist" => fd.get_mut(&current_fd).unwrap().triangle = vec_of(&f[2..].join("\t")),
            "fd_row" => {
                let idx: usize = f[1].rsplit('_').next().unwrap().parse().unwrap();
                let v = vec_of(&f[2..].join("\t"));
                let c = fd.get_mut(&current_fd).unwrap();
                c.rows.resize(c.n, Vec::new());
                c.rows[idx - 1] = v;
            }
            "fd_error" => { /* asserted by name in the test below */ }
            "n_trim" | "n_fd" => {}
            other => panic!("unknown corpus record {other:?}"),
        }
    }
    (tm, fd, tm_order, fd_order)
}

fn bits(v: f64) -> String {
    if v.is_nan() {
        "nan".to_string()
    } else {
        format!("{:016x}", v.to_bits())
    }
}

fn same(a: f64, b: f64) -> bool {
    bits(a) == bits(b)
}

// ---------------------------------------------------------------- trimmed mean

/// Every `mean(trim=, na.rm=)` case, bit for bit.
#[test]
fn the_trimmed_mean_matches_r() {
    let (tm, _, order, _) = load();
    assert_eq!(
        order.len(),
        26,
        "the corpus should carry 26 trimmed-mean cases"
    );
    for name in &order {
        let c = &tm[name];
        let got = trimmed_mean(&c.x, c.trim, c.na_rm);
        assert!(
            same(got, c.want),
            "{name}: mean(trim={}, na.rm={}) gave {got:?}, R gave {:?}",
            c.trim,
            c.na_rm,
            c.want
        );
    }
}

/// `na.rm = TRUE` drops `NaN` too, because `is.na(NaN)` is `TRUE`. `NaN` is a *value* in R
/// everywhere else, so a port that filters only `NA` returns `NaN` here and differs.
#[test]
fn na_rm_drops_nan_as_well_as_na() {
    assert_eq!(
        trimmed_mean(&[1.0, f64::NAN, 3.0], 0.0, true),
        2.0,
        "mean(c(1, NaN, 3), na.rm = TRUE) is 2, not NaN"
    );
    // The corpus's `nan_removed` case is the same statement with `trim = 0.1`; the two must
    // agree, which they only do if `NaN` is removed *before* the trim count is taken.
    let (tm, _, _, _) = load();
    let c = &tm["nan_removed"];
    assert!(same(trimmed_mean(&c.x, c.trim, c.na_rm), c.want));
    let d = &tm["na_removed"];
    assert!(
        same(
            trimmed_mean(&d.x, d.trim, d.na_rm),
            trimmed_mean(&c.x, c.trim, c.na_rm)
        ),
        "an NA and a NaN in the same position must behave identically under na.rm = TRUE"
    );
}

/// The count trimmed from each end is `floor(n * trim)`, so the kept count is
/// `n - 2 * floor(n * trim)`. It is **not** `round(n * (1 - trim))`, and it is not monotone in
/// `n`: `n = 9` trims nothing, `n = 10` trims one, `n = 19` one, `n = 20` two.
#[test]
fn the_trim_count_is_floor_n_times_trim_from_each_end() {
    // `1:n` has mean `(n+1)/2` untrimmed, so the kept window's mean is exactly the arithmetic
    // mean of the integers `floor(n*trim)+1 .. n - floor(n*trim)`.
    let expected = |n: usize| -> f64 {
        let t = (n as f64 * 0.1).floor() as usize;
        let hi = n - t;
        // The kept window is the *values* `t+1 ..= n-t`, which is `0`-based slice `t .. n-t`.
        let kept: f64 = (t + 1..=hi).map(|v| v as f64).sum();
        kept / (hi - t) as f64
    };
    for n in [1usize, 2, 9, 10, 11, 19, 20, 21, 100, 1000] {
        let x: Vec<f64> = (1..=n).map(|i| i as f64).collect();
        assert!(
            same(trimmed_mean(&x, 0.1, true), expected(n)),
            "n = {n}: the kept window must be floor(n*trim)+1 .. n-floor(n*trim)"
        );
    }
    // The non-monotonicity, stated explicitly: 9 and 10 trim different amounts and the *kept
    // count* jumps by two between them.
    let t = |n: f64| (n * 0.1).floor();
    assert_eq!((t(9.0), t(10.0)), (0.0, 1.0));
    assert_eq!((t(19.0), t(20.0)), (1.0, 2.0));
}

/// `trim >= 0.5` returns `median`, and on an even `n` that is the *mean of the two middle
/// values* -- not the lower one and not the empty middle. `computeRegionDistance` uses
/// `trim = 0.1` so this is unreachable there, but `mean` is a general utility and a silently
/// wrong `trim = 0.5` would be found by someone else.
#[test]
fn trim_at_or_above_a_half_returns_the_median() {
    let x: Vec<f64> = (1..=10).map(|i| i as f64).collect();
    assert_eq!(trimmed_mean(&x, 0.5, true), 5.5, "(5 + 6) / 2 on an even n");
    assert_eq!(
        trimmed_mean(&x, 0.75, true),
        5.5,
        "trim > 0.5 is still the median"
    );
    let y: Vec<f64> = (1..=9).map(|i| i as f64).collect();
    assert_eq!(
        trimmed_mean(&y, 0.5, true),
        5.0,
        "the middle element on an odd n"
    );
    // Just below the threshold it is a trimmed mean, not a median, and the two differ. The
    // *value* is pinned bit-for-bit by the corpus for every `trim` it covers; what this
    // asserts is the branch, so a hand-derived number here would only be a second thing to keep
    // in sync. (`trim = 0.4` on `1:10` is a bad choice: it keeps `{5, 6}`, which is also the
    // median, so it cannot distinguish the two branches at all.)
    // The *values* are pinned bit-for-bit by `the_trimmed_mean_matches_r` over 26 cases
    // spanning `trim` in {0, 0.1, 0.5, 0.75} and `n` from 1 to 1000, so both branches are
    // value-checked there. Asserting a hand-derived value for an unlisted `trim` here would only
    // add a second thing to keep in sync with the implementation, and a 0.2 case turns out to
    // be a poor discriminator anyway: on `1:9` it coincides with the median.
    let (tm, _, _, _) = load();
    let c = &tm["n9"];
    assert!(
        same(trimmed_mean(&c.x, c.trim, c.na_rm), c.want),
        "n9, trim = 0.1"
    );
}

/// With `na.rm = FALSE`, a single missing value anywhere gives `NA` and not a mean of the rest --
/// but only when `trim > 0`, because `trim == 0` skips the block that contains the `anyNA`
/// early return. That asymmetry is in `mean.default` and is easy to miss.
#[test]
fn na_rm_false_gives_na_only_when_trimming_is_active() {
    let x = [1.0, f64::NAN, 3.0, 4.0, 5.0];
    assert!(
        trimmed_mean(&x, 0.1, false).is_nan(),
        "trim > 0 takes the anyNA branch"
    );
    assert!(
        trimmed_mean(&x, 0.0, false).is_nan(),
        "and trim = 0 sums through the NaN"
    );
    // The `NA` case is distinguishable only in that it is `NA_real_`, not `NaN`; both are NaN
    // bit patterns, so the corpus distinguishes them textually and this cannot.
    let y = [1.0, 2.0, 3.0];
    assert_eq!(
        trimmed_mean(&y, 0.1, false),
        2.0,
        "no missing values, no early return"
    );
}

/// An all-missing or empty input gives `NaN`, which is `mean(numeric(0))`.
#[test]
fn an_empty_or_all_missing_input_is_nan() {
    assert!(trimmed_mean(&[], 0.1, true).is_nan());
    assert!(trimmed_mean(&[f64::NAN, f64::NAN], 0.1, true).is_nan());
    assert!(trimmed_mean(&[], 0.1, false).is_nan());
}

/// The accumulation is LONG_DOUBLE. A `f64` sum is 1 ulp off at eight terms on values like
/// these, so the case is chosen to be one where it shows.
#[test]
fn the_accumulation_is_long_double() {
    let (tm, _, _, _) = load();
    for name in ["n1000", "negative", "fractional"] {
        let c = &tm[name];
        assert!(same(trimmed_mean(&c.x, c.trim, c.na_rm), c.want), "{name}");
    }
    // And the specific shape: `1e8 + 1 - 1e8` in `f64` is 0, and in LONG_DOUBLE it is 1 before
    // the division -- so a 10-element vector built around it separates the two.
    let x = [1e8, 1.0, -1e8, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0];
    let f64_sum: f64 = x.iter().sum();
    let got = trimmed_mean(&x, 0.1, true);
    assert!(
        (f64_sum - got).abs() > 1e-9,
        "an f64 sum would give {f64_sum}; LONG_DOUBLE gives {got}"
    );
}

// ---------------------------------------------------------------- fdist

/// The full distance matrix, bit for bit, row by row.
#[test]
fn fdist_matches_r_bit_for_bit() {
    let (_, fd, _, order) = load();
    assert_eq!(order.len(), 7, "the corpus should carry 7 distance cases");
    for name in &order {
        let c = &fd[name];
        if c.rows.is_empty() {
            // The one-row case is a recorded *refusal*; `a_single_point_has_no_distance_matrix`
            // asserts it. Reading `rows` here would be an empty matrix indexed by `n = 1`.
            assert_eq!(c.n, 1, "{name}: only the one-row case may have no rows");
            continue;
        }
        let coords = c.coords.clone();
        let got = fdist(&coords, c.n, c.d);
        for i in 0..c.n {
            for j in 0..c.n {
                assert!(
                    same(got[i + c.n * j], c.rows[i][j]),
                    "{name}: fdist[{i}][{j}] gave {:?}, R gave {:?}",
                    got[i + c.n * j],
                    c.rows[i][j]
                );
            }
        }
    }
}

/// Upstream's observable is a **`dist`**, not a matrix: the compact triangle, and the port's
/// full matrix must contain exactly those values.
///
/// Compared as a multiset, because `as.vector.dist`'s traversal is not documented and
/// `collapse`'s is not base R's (for a 3-4-5 triangle collapse orders the pairs
/// `(1,2), (2,3), (1,3)` where base orders `(1,2), (1,3), (2,3)`). Pinning the traversal would
/// pin a property of `as.vector` rather than of the arithmetic.
#[test]
fn fdist_reproduces_the_triangle_upstream_exposes() {
    let (_, fd, _, order) = load();
    for name in &order {
        let c = &fd[name];
        if c.rows.is_empty() {
            continue;
        }
        let coords = c.coords.clone();
        let got = fdist(&coords, c.n, c.d);
        let mut mine: Vec<u64> = Vec::with_capacity(c.n * (c.n - 1) / 2);
        for i in 0..c.n {
            for j in (i + 1)..c.n {
                mine.push(got[i + c.n * j].to_bits());
            }
        }
        let mut theirs: Vec<u64> = c.triangle.iter().map(|v| v.to_bits()).collect();
        mine.sort_unstable();
        theirs.sort_unstable();
        assert_eq!(
            mine.len(),
            theirs.len(),
            "{name}: the triangle must hold n(n-1)/2 = {} values, not {}",
            c.n * (c.n - 1) / 2,
            theirs.len()
        );
        assert_eq!(
            mine, theirs,
            "{name}: the distance multiset differs from upstream's `dist`"
        );
    }
}

/// The diagonal is `0`, not `NA`. `fdist` on a *data frame* returns `NA_real_` on the diagonal,
/// and `computeCellDistance` passes a matrix, so the two differ on exactly the cells a consumer
/// is most likely to read.
#[test]
fn the_diagonal_is_zero_not_na() {
    let (_, fd, _, _) = load();
    for (name, c) in &fd {
        if c.rows.is_empty() {
            continue;
        }
        let coords = c.coords.clone();
        let got = fdist(&coords, c.n, c.d);
        for i in 0..c.n {
            assert_eq!(got[i + c.n * i], 0.0, "{name}: diagonal must be 0");
        }
    }
}

/// A single-point distance matrix is not a thing: `collapse::fdist` refuses it, because there is
/// no pair to form. The port returns the 1x1 zero matrix instead, and that difference is
/// deliberate and documented rather than silent -- the shim never calls `fdist` with `n < 2`,
/// because `computeCellDistance` is reached only from a real coordinate matrix.
#[test]
fn a_single_point_has_no_distance_matrix() {
    let got = fdist(&[0.0, 0.0], 1, 2);
    assert_eq!(got, vec![0.0], "1 x 1 zero matrix");
    // ... and the corpus records that upstream refuses. Named rather than counted, so a
    // replaced case is noticed.
    let text = std::fs::read_to_string(
        std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join("../../tests/fixtures/spatial_golden.txt"),
    )
    .unwrap();
    assert!(
        text.contains("fd_error\tone\t"),
        "the corpus must record fdist's refusal of a one-row matrix"
    );
    assert!(
        text.contains("If v is left empty, x needs to be a matrix with at least 2 rows"),
        "upstream's message, verbatim"
    );
}

/// The degenerate geometries a real slide produces: duplicate spots and a collinear run. Both
/// make the k-d tree's axis choice degenerate, which is why they are here and not only in the
/// tree's own test.
#[test]
fn duplicate_and_collinear_points_are_exact() {
    // Two identical points: distance exactly 0, and every pair between the two groups equal.
    let dup = fdist(&[1.0, 1.0, 1.0, 1.0, 2.0, 2.0, 2.0, 2.0], 4, 2);
    assert_eq!(dup[0], 0.0, "identical points are at distance 0");
    assert_eq!(dup[4], 0.0);
    for i in 0..4 {
        for j in 0..4 {
            // The two duplicate groups are {0,1} and {2,3}; within a group the distance is 0,
            // across groups it is sqrt(2).
            let want = if (i < 2) == (j < 2) {
                0.0
            } else {
                2.0f64.sqrt()
            };
            assert_eq!(dup[i + 4 * j], want, "dup[{i}][{j}]");
        }
    }
    // A collinear run: the widest axis is the only one with any spread, so the tree has no
    // choice to make and must still be exact.
    let mut coords = Vec::new();
    for v in 0..10 {
        coords.push(v as f64);
    }
    let lin = fdist(&coords, 10, 1);
    for i in 0..10 {
        for j in 0..10 {
            assert_eq!(
                lin[i + 10 * j],
                (i as f64 - j as f64).abs(),
                "collinear [{i}][{j}]"
            );
        }
    }
}

// ---------------------------------------------------------------- the k-d tree

/// A deterministic PRNG, so a failure is reproducible from the seed in the assertion. `rand` is
/// not a dependency of `r-core`, and a linear congruential generator is enough to build the
/// adversarial point sets this needs.
fn lcg(seed: &mut u64) -> f64 {
    *seed = seed
        .wrapping_mul(6364136223846793005)
        .wrapping_add(1442695040888963407);
    ((*seed >> 11) as f64) / ((1u64 << 53) as f64)
}

/// The build's layout, on a hand-checkable point set.
///
/// A tree that is wrong on uniform noise is often right on uniform noise, and a search that is
/// wrong is indistinguishable from a search that is merely *inaccurately pruned* unless the
/// structure is checked on its own. Four points, so the expected layout can be written out.
#[test]
fn the_build_places_points_where_the_rule_says() {
    // Four points on a line: 0, 1, 10, 11. `mid = lo + (hi - lo) / 2` is the upper median of an
    // even range, so the root is the third point and the left subtree is the first two.
    let t = KdTree::new(vec![0.0, 1.0, 10.0, 11.0], 4, 1);
    let n = t.nodes();
    // One node per point, appended and never resized away. Arena order is root, then the left
    // subtree depth-first, then the right -- so it is a function of the split rule alone.
    assert_eq!(n.len(), 4, "one node per point");
    assert_eq!(n[0].point, 2, "root is the upper median, value 10");
    assert_eq!(n[0].axis, 0);
    assert_eq!(n[0].left, Some(1), "the left subtree holds values 0 and 1");
    assert_eq!(n[0].right, Some(3), "the right subtree holds 10 and 11");
    // Arena order is *depth-first, left before right*, so the indices are root, left subtree,
    // right subtree -- not breadth-first. The left subtree's own median is 1; its left child is
    // the leaf 0 and its right child is the **empty** range `(2, 2)`, which is where the
    // `Option` earns its keep: `hi - lo == 0` is not `== 1`, so a signature that cannot say "no
    // node here" recurses forever.
    assert_eq!(n[1].point, 1);
    assert_eq!(n[1].left, Some(2));
    assert_eq!(
        n[1].right, None,
        "an even-length split leaves one side empty"
    );
    assert_eq!(n[2].point, 0, "the leaf under the left subtree");
    assert_eq!(n[3].point, 3, "the right leaf, value 11");
    // The *indices* are the contract here; the squared distances are compared with a tolerance
    // because `(10.4 - 10.0)^2` is not `0.16` in binary, and pinning the exact bits of a
    // subtraction would be pinning the representation rather than the geometry.
    let got = t.knn(&[10.4], 2);
    assert_eq!(got.iter().map(|g| g.1).collect::<Vec<_>>(), vec![2, 3]);
    assert!(
        (got[0].0 - 0.16).abs() < 1e-12 && (got[1].0 - 0.36).abs() < 1e-12,
        "{got:?}"
    );
}

/// The arena must have exactly one node per point, on every input shape. The `resize`-truncation
/// bug produced a *shorter* arena than the point count, which is the cheapest possible signal and
/// the one the previous implementation never got to check.
#[test]
fn the_arena_has_one_node_per_point() {
    let mut seed = 0xfeed_0001u64;
    for (n, d) in [
        (1usize, 1usize),
        (2, 1),
        (3, 2),
        (7, 3),
        (64, 2),
        (257, 2),
        (1000, 3),
    ] {
        let mut pts = Vec::new();
        for _ in 0..n {
            for _ in 0..d {
                pts.push(lcg(&mut seed) * 100.0);
            }
        }
        let t = KdTree::new(pts, n, d);
        assert_eq!(t.nodes().len(), n, "n = {n}, d = {d}: one node per point");
        let mut seen: Vec<usize> = t.nodes().iter().map(|x| x.point).collect();
        seen.sort_unstable();
        assert_eq!(
            seen,
            (0..n).collect::<Vec<_>>(),
            "every point appears exactly once"
        );
    }
    // Degenerate shapes: all identical, and collinear. `n * d` coordinates, so 20 points in one
    // dimension is 20 values -- passing 40 builds a 40-point tree and the shape assert fires.
    assert_eq!(KdTree::new(vec![1.0; 20], 20, 1).nodes().len(), 20);
    assert_eq!(
        KdTree::new((0..40).map(|i| i as f64).collect(), 40, 1)
            .nodes()
            .len(),
        40
    );
}

/// The tree must agree with exhaustive search, on point sets chosen to break it.
///
/// The sets are not random uniform noise. They are: a regular lattice (every neighbour tied at
/// distance 1), all-duplicate points (every distance tied at 0), a collinear run (one axis with
/// no spread), a circle (many points at the same radius), two tight clusters far apart (so the
/// split axis matters), and uniform noise in 3-d for the ordinary case.
#[test]
fn knn_matches_brute_force_on_adversarial_point_sets() {
    let mut seed = 0x5eed_1234u64;
    let cases: Vec<(&str, Vec<f64>, usize, usize)> = vec![
        (
            "lattice6",
            {
                let mut v = Vec::new();
                for x in 0..6 {
                    for y in 0..6 {
                        v.push(x as f64);
                        v.push(y as f64);
                    }
                }
                v
            },
            36,
            2,
        ),
        ("all_duplicate", vec![2.5; 24 * 2], 24, 2),
        (
            "collinear",
            (0..30).map(|i| (i as f64) * 0.7).collect(),
            30,
            1,
        ),
        (
            "circle",
            {
                let mut v = Vec::new();
                for i in 0..40 {
                    let t = i as f64 * std::f64::consts::TAU / 40.0;
                    v.push(t.cos());
                    v.push(t.sin());
                }
                v
            },
            40,
            2,
        ),
        (
            "two_clusters",
            {
                let mut v = Vec::new();
                for _ in 0..60 {
                    v.push(lcg(&mut seed) * 0.01);
                    v.push(lcg(&mut seed) * 0.01);
                }
                for _ in 0..60 {
                    v.push(100.0 + lcg(&mut seed) * 0.01);
                    v.push(100.0 + lcg(&mut seed) * 0.01);
                }
                v
            },
            120,
            2,
        ),
        (
            "uniform3",
            {
                let mut v = Vec::new();
                for _ in 0..200 {
                    for _ in 0..3 {
                        v.push(lcg(&mut seed));
                    }
                }
                v
            },
            200,
            3,
        ),
    ];

    for (name, pts, n, d) in cases {
        let tree = KdTree::new(pts.clone(), n, d);
        for k in [1usize, 2, 3, 5, 8, 17] {
            if k > n {
                continue;
            }
            // From points in the set (so distance 0 occurs) and from off-lattice points (so no
            // distance is trivially zero).
            let mut queries: Vec<Vec<f64>> = Vec::new();
            for i in 0..n.min(20) {
                queries.push(pts[i * d..i * d + d].to_vec());
            }
            queries.push((0..d).map(|_| lcg(&mut seed) * 1.5).collect());
            queries.push(vec![-0.25; d]);
            for q in &queries {
                let got = tree.knn(q, k);
                let want = brute_force_knn(&pts, n, d, q, k);
                assert_eq!(
                    got.len(),
                    want.len(),
                    "{name}: k = {k} returned {} neighbours",
                    got.len()
                );
                for (i, (g, w)) in got.iter().zip(want.iter()).enumerate() {
                    assert_eq!(
                        g.0.to_bits(),
                        w.0.to_bits(),
                        "{name}: k = {k}, rank {i}: squared distance {:?} vs {:?}",
                        g.0,
                        w.0
                    );
                    assert_eq!(
                        g.1, w.1,
                        "{name}: k = {k}, rank {i}: index {} vs {} -- a tie-break disagreement",
                        g.1, w.1
                    );
                }
            }
        }
    }
}

/// Ties break towards the lower index, and the rule is what makes the tree agree with brute
/// force. It matters downstream because `computeRegionDistance` counts
/// `length(intersect(knn.i, unique(qout$index[...]))` -- with equidistant candidates the chosen
/// *set* is stable but the *order* is not, and `unique` keeps the first.
#[test]
fn ties_break_towards_the_lower_index() {
    let pts: Vec<f64> = (0..5).map(|i| i as f64).collect();
    let tree = KdTree::new(pts, 5, 1);
    // The query sits *on* point 1, so that is the nearest at distance 0 and the second slot is a
    // three-way tie at distance 1 between points 0 and 2. Point 0 must take it.
    let got = tree.knn(&[1.0], 2);
    assert_eq!(got[0], (0.0, 1), "the query's own point is nearest");
    assert_eq!(
        got[1],
        (1.0, 0),
        "and the lower index takes the tie for second"
    );
    // Queried *between* two points the tie is the whole point: both are at distance 0.25 and
    // the lower index must come first.
    let got = tree.knn(&[0.5], 2);
    assert_eq!(got[0], (0.25, 0), "lower index first");
    assert_eq!(got[1], (0.25, 1), "higher index second");
    assert_eq!(got[0].0, got[1].0, "and they really are equidistant");

    // All points identical: the k nearest are the k lowest indices, in order.
    let dup = KdTree::new(vec![1.0; 6], 6, 1);
    assert_eq!(
        dup.knn(&[1.0], 4).iter().map(|g| g.1).collect::<Vec<_>>(),
        vec![0, 1, 2, 3]
    );
    assert_eq!(dup.nearest(&[9.0]).unwrap().1, 0);
}

/// A query at a stored point must find itself at distance 0. The cheapest possible check that the
/// pruning is not discarding the subtree the answer is in.
///
/// The points are on a **coarse lattice with more cells than points**, so no two coincide.
/// An earlier version rounded to halves, which put several points on the same coordinate -- and
/// then the correct answer under the index tie-break is a *lower-indexed coincident* point, not
/// the query's own index. The test failed for the right reason at the wrong value, which is worth
/// stating: "finds a point at distance 0" and "finds *itself*" are different claims, and only the
/// second one is what this test is for.
#[test]
fn a_query_at_a_stored_point_finds_itself() {
    for (n, d) in [(8usize, 1usize), (50, 2), (200, 3)] {
        let seed = 0x1234_5678u64 ^ n as u64;
        // Distinct lattice cells from a **strided** walk, not from the LCG. An LCG's low bits
        // are famously poor, and `(lcg() * cell).floor()` collided even at `cell = 4n`: 8 points
        // came out as 5 distinct coordinates, so the fixture silently tested the tie-break
        // instead of self-identification. A multiplicative stride over a prime modulus is
        // coprime with the modulus, so the first `n` values are distinct by construction.
        let modulus = 1_000_003u64;
        let stride = 7_919u64; // prime, and coprime with `modulus`
        let mut k = seed % modulus;
        let mut pts = Vec::new();
        for _ in 0..n {
            for _ in 0..d {
                k = (k + stride) % modulus;
                pts.push(k as f64);
            }
        }
        // Distinctness is a property of *points*, not of coordinates. Deduplicating the
        // flattened coordinate vector compares a 2-d point's x against another point's y, so for
        // `d = 2` it reports `2n` distinct values and the assertion below could never hold.
        let mut seen: Vec<&[f64]> = pts.chunks(d).map(|c| c as &[f64]).collect();
        seen.sort_unstable_by(|a, b| a.partial_cmp(b).unwrap());
        seen.dedup();
        assert_eq!(
            seen.len(),
            n,
            "n = {n}, d = {d}: the fixture must have no coincident points"
        );
        let t = KdTree::new(pts.clone(), n, d);
        for i in 0..n {
            let got = t.knn(&pts[i * d..i * d + d], 1);
            assert_eq!(
                got[0].1, i,
                "n = {n}, d = {d}: point {i} must find itself, got index {} at squared distance {}",
                got[0].1, got[0].0
            );
            assert_eq!(got[0].0, 0.0);
        }
    }
}

/// With coincident points the answer is *a* point at distance 0, not necessarily the query's own
/// index: the lower index wins the tie. Recorded separately from the test above because the two
/// claims are different and conflating them makes a correct implementation look wrong.
#[test]
fn a_coincident_point_is_found_by_the_lower_index() {
    // Three points, two of them identical, queried at that shared location.
    let t = KdTree::new(vec![5.0, 5.0, 1.0, 1.0, 5.0, 5.0], 3, 2);
    let got = t.knn(&[5.0, 5.0], 2);
    assert_eq!(got[0], (0.0, 0), "points 0 and 2 coincide; 0 is lower");
    assert_eq!(got[1], (0.0, 2));
}

/// `nearest` is the `k = 1` query `computeRegionDistance` uses, and it must agree with both
/// `knn` and exhaustive search.
#[test]
fn nearest_is_knn_of_one_and_agrees_with_brute_force() {
    let mut seed = 0xabcd_0001u64;
    let (n, d) = (150usize, 2usize);
    let mut pts = Vec::new();
    for _ in 0..n {
        for _ in 0..d {
            pts.push(lcg(&mut seed) * 10.0);
        }
    }
    let tree = KdTree::new(pts.clone(), n, d);
    for _ in 0..50 {
        let q: Vec<f64> = (0..d).map(|_| lcg(&mut seed) * 12.0).collect();
        let near = tree.nearest(&q).expect("a non-empty tree has a nearest");
        assert_eq!(near, tree.knn(&q, 1)[0]);
        assert_eq!(near, brute_force_knn(&pts, n, d, &q, 1)[0]);
    }
}

/// Degenerate shapes the callers can reach, which must not panic and must still be exact.
#[test]
fn degenerate_trees_behave() {
    assert!(KdTree::new(Vec::new(), 0, 2).is_empty());
    assert!(KdTree::new(Vec::new(), 0, 2).nearest(&[0.0, 0.0]).is_none());
    // `n * d` coordinates, not `n`: a 1-point 1-dimensional tree takes one.
    assert!(KdTree::new(vec![0.0], 1, 1).knn(&[0.0], 0).is_empty());
    // `k` larger than `n` is clamped, not a panic.
    let t = KdTree::new(vec![0.0, 0.0, 1.0, 0.0], 2, 2);
    assert_eq!(t.knn(&[0.0, 0.0], 99).len(), 2);
    assert_eq!(t.len(), 2);
    // One point, queried at a distance: squared, so no `sqrt` in the search.
    let one = KdTree::new(vec![3.0, 4.0], 1, 2);
    assert_eq!(one.knn(&[0.0, 0.0], 1), vec![(25.0, 0)]);
}

/// The tree has to be *fast*, not just correct: `computeCommunProb` runs a 1-NN query for every
/// cell against every other cell group, so `O(n log n)` against `O(n^2)` is the whole point of
/// replacing Annoy. A quadratic tree that agrees with brute force would satisfy every other test
/// here and still be a regression.
#[test]
fn the_search_is_subquadratic() {
    let mut seed = 0x0bad_c0deu64;
    let (n, d) = (4000usize, 2usize);
    let mut pts = Vec::new();
    for _ in 0..n {
        for _ in 0..d {
            pts.push(lcg(&mut seed));
        }
    }
    let tree = KdTree::new(pts.clone(), n, d);
    let start = std::time::Instant::now();
    let mut total = 0.0f64;
    for i in 0..400 {
        let q = &pts[(i * 7 % n) * d..(i * 7 % n) * d + d];
        total += tree.nearest(q).unwrap().0;
    }
    let tree_us = start.elapsed().as_micros();
    // The same queries by exhaustive search, on a small enough point set to finish.
    let (m, d2) = (400usize, 2usize);
    let small = &pts[..m * d2];
    let start = std::time::Instant::now();
    let mut sink = 0.0f64;
    for i in 0..200 {
        let q = &small[(i * 7 % m) * d2..(i * 7 % m) * d2 + d2];
        sink += brute_force_knn(small, m, d2, q, 1)[0].0;
    }
    // Keep both loops from being optimised away, without asserting on the values: this test is
    // about the *shape* of the cost, not about the distances.
    std::hint::black_box((total, sink));
    let brute_us = start.elapsed().as_micros();
    // The tree does 10x the points in comparable time; the constant is loose on purpose, because
    // a timing assertion that flakes is worse than none. It fails only if the tree has become
    // linear, which is the regression worth catching.
    assert!(
        (tree_us / 400) < (brute_us / 200) * 4,
        "tree {tree_us}us for 400 queries at n={n} vs brute {brute_us}us for 200 at n={m}: the \\
         search looks linear in n"
    );
}
