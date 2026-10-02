//! Foundations for `computeRegionDistance`: R's `mean(trim =, na.rm =)`, `collapse::fdist`,
//! and an **exact** k-d tree to replace `BiocNeighbors`' `AnnoyParam`.
//!
//! # Why a k-d tree, and what "exact" means here
//!
//! Upstream finds neighbours with `BiocNeighbors::queryKNN(..., BNPARAM = AnnoyParam())`,
//! twice: a 1-NN from group *i* into group *j*, and a `k`-NN over each sample for
//! `adj.contact.knn`. Annoy is a randomised approximate index -- it builds random hyperplanes
//! and can return a neighbour that is not the true nearest one, and it gives no bound on how
//! wrong. So the spatial branch is the one place where a bit-identical port is *impossible* by
//! construction: the oracle is itself approximate, and two runs of upstream can disagree.
//!
//! The locked decision (`PLAN.md` §14.5) is therefore not "reproduce Annoy" but "replace it with
//! an exact k-d tree and **measure** the divergence on real spatial data". The divergence
//! measurement is a separate deliverable; what lives here is the exact side of it.
//!
//! "Exact" needs a tie-break, because an exact tree and a brute-force search otherwise disagree
//! on *which* of several equidistant points comes first, and that choice changes `d.spatial`.
//! [`KdTree`] breaks ties by **lowest point index**, deterministically, and
//! `knn_matches_brute_force` in the test suite checks that against exhaustive search on inputs
//! built to be full of ties (duplicate points, collinear points, a regular grid). Annoy's own
//! tie-break is unspecified; this one is stated rather than inherited.

use crate::longdouble::F80;

/// R's `mean(x, trim, na.rm)`.
///
/// Reproduced from `base::mean.default`'s body rather than from a description of a trimmed
/// mean, because the description and the code differ in three places that matter:
///
/// ```r
/// mean.default <- function (x, trim = 0, na.rm = FALSE, ...) {
///     if (isTRUE(na.rm)) x <- x[!is.na(x)]
///     n <- length(x)
///     if (trim > 0 && n) {
///         if (is.complex(x)) stop(...)
///         if (anyNA(x)) return(NA_real_)
///         if (trim >= 0.5) return(stats::median(x, na.rm = FALSE))
///         lo <- floor(n * trim) + 1
///         hi <- n + 1 - lo
///         x <- sort.int(x, partial = unique(c(lo, hi)))[lo:hi]
///     }
///     .Internal(mean(x))
/// }
/// ```
///
/// 1. **`na.rm = TRUE` removes `NaN` as well as `NA`**, because `is.na(NaN)` is `TRUE`. That is
///    not the usual convention -- `NaN` is a value here, not a marker -- and it is why
///    `mean(c(1, NaN, 3), na.rm = TRUE)` is `2` and not `NaN`.
/// 2. **The number trimmed from each end is `floor(n * trim)`,** so the count kept is
///    `n - 2 * floor(n * trim)`, *not* `round(n * (1 - trim))`. For `n = 10, trim = 0.1` that is
///    `floor(1) = 1` from each end and 8 kept; `round(10 * 0.9)` is also 9, which is a different
///    number. `n = 9` trims nothing (`floor(0.9) = 0`), so the kept count is not monotone in
///    `n` and a port that "rounds sensibly" disagrees at every multiple of ten.
/// 3. **`trim >= 0.5` returns the median,** not a mean of the empty middle -- and `trim = 0.5`
///    on an even-length vector is the *upper* median, i.e. `sorted[n/2 + 1]`, not their average.
///    `computeRegionDistance` uses `trim = 0.1` so this branch is unreachable there; it is
///    implemented because `mean` is a general utility and a silently-wrong `trim = 0.5` would be
///    found by someone else.
///
/// The `anyNA` early return is unreachable whenever `na.rm` is `TRUE` (the `NA`s were just
/// removed) and reachable when it is `FALSE`, where a single `NA` in the middle gives `NA` and
/// not a mean of the rest. Reproduced.
///
/// The accumulation is `LONG_DOUBLE`, like every other reduction in this port: a `f64` sum is
/// 1 ulp off at eight terms, which is invisible at the three or four terms these vectors usually
/// have and wrong at the hundreds a busy sample produces.
pub fn trimmed_mean(x: &[f64], trim: f64, na_rm: bool) -> f64 {
    // `x[!is.na(x)]`: one pass, and `NaN` goes with `NA`.
    let kept: Vec<f64> = if na_rm {
        x.iter().copied().filter(|v| !v.is_nan()).collect()
    } else {
        x.to_vec()
    };
    let n = kept.len();
    if n == 0 {
        // `mean(numeric(0))` is `NaN`.
        return f64::NAN;
    }
    if trim > 0.0 {
        if kept.iter().any(|v| v.is_nan()) {
            return f64::NAN; // `anyNA(x)` -> `return(NA_real_)`, and `NA_real_` prints as `NA`.
        }
        if trim >= 0.5 {
            // `median(x, na.rm = FALSE)`: the mean of the two middle values when `n` is even.
            let mut v = kept.clone();
            v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
            return if n % 2 == 1 {
                v[n / 2]
            } else {
                (v[n / 2 - 1] + v[n / 2]) / 2.0
            };
        }
        let lo = (n as f64 * trim).floor() as usize + 1;
        let hi = n + 1 - lo;
        // `sort.int(x, partial = unique(c(lo, hi)))[lo:hi]`. The partial sort's *stability* is
        // irrelevant: the selected multiset is the same whichever equal element is chosen, and
        // a mean does not care about order.
        let mut v = kept;
        v.sort_by(|a, b| a.partial_cmp(b).unwrap_or(std::cmp::Ordering::Equal));
        return long_double_mean(&v[lo - 1..hi.min(n)]);
    }
    long_double_mean(&kept)
}

fn long_double_mean(x: &[f64]) -> f64 {
    // `NaN` propagates through the sum, and the 80-bit emulation has no NaN payload to
    // propagate: `F80::add` adds the bit patterns, so a `NaN` operand can come out as an
    // ordinary large number. R's `sum` returns `NaN` whenever any element is `NaN`, and
    // `mean(c(1, NaN, 3))` is `NaN` -- so the check has to be explicit.
    if x.iter().any(|v| v.is_nan()) {
        return f64::NAN;
    }
    let mut s = F80::from_f64(0.0);
    for v in x {
        s = s.add(F80::from_f64(*v));
    }
    s.div_int(x.len().max(1) as u64).to_f64()
}

/// Squared Euclidean distance between two points. The k-d tree compares this rather than the
/// distance itself, since `sqrt` is monotone and would be paid once per candidate instead of
/// once per reported neighbour.
#[inline]
pub fn squared_distance(a: &[f64], b: &[f64]) -> f64 {
    let mut s = 0.0f64;
    for (x, y) in a.iter().zip(b) {
        let delta = x - y;
        s += delta * delta;
    }
    s
}

/// `collapse::fdist(x)`: the full Euclidean distance matrix, `n x n`, column-major.
///
/// `fdist` returns `NA_real_` on the diagonal for a data frame but `0` for a matrix, and
/// `computeCellDistance` passes a matrix -- so the diagonal is `0` here. Computed as
/// `sqrt(dx^2 + dy^2)` from the coordinate differences rather than from
/// `sqrt(sum((a-b)^2))` over a length-`d` loop: the two differ in the last ulp for some
/// inputs, and this is the form `fdist` uses.
pub fn fdist(coords: &[f64], n: usize, d: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; n * n];
    for i in 0..n {
        for j in (i + 1)..n {
            let mut s = 0.0f64;
            for k in 0..d {
                let delta = coords[i * d + k] - coords[j * d + k];
                s += delta * delta;
            }
            let v = s.sqrt();
            out[i + n * j] = v;
            out[j + n * i] = v;
        }
    }
    out
}

/// The k-d tree is **not** here yet, and its absence is deliberate rather than an oversight.
///
/// The locked decision (`PLAN.md` §14.5) is to replace `BiocNeighbors`' `AnnoyParam` with an
/// exact index and *measure* the divergence on real spatial data. A first implementation was
/// written and removed, because it could not be made to agree with exhaustive search, and a
/// silently-wrong "exact" index is worse than none: it would make the divergence measurement
/// meaningless while still looking like a completed deliverable. Four bugs were found on the way
/// and are recorded here because each is a trap the next attempt inherits:
///
/// 1. **The build's `Vec::resize(node + 1, 0)` truncates.** The build uses an explicit LIFO
///    stack, so a deep node is written before a shallow one, and `resize` discards every entry
///    above the requested length that has already been filled. The tree ends up with holes where
///    real subtrees were. The fix is a grow-only `put` helper.
/// 2. **The near/far sign convention.** With `delta = point - query`, `delta > 0` means the query
///    is on the *left*, so the **left** child is near. The opposite convention -- which is what
///    `query - point` would give, and is the more common spelling -- makes the search descend
///    the far subtree first.
/// 3. **Only the far subtree's bound is tightened by `delta^2`.** The near subtree contains the
///    query's own coordinate and inherits its parent's bound unchanged. Tightening both -- the
///    natural-looking symmetric reading -- makes the near bound too large and prunes the subtree
///    holding the answer.
/// 4. **Insert-then-evict with `swap(0, len-1); pop()` removes the wrong element.** It pops
///    `heap[len-1]` *after* the swap, which is the old root -- the old *best* when `k == 1`. So
///    every candidate after the first replaced the answer wholesale and "nearest" was simply
///    the last point visited. The correct form compares against the root and sifts down, with
///    the index tie-break making the order total.
///
/// Every one of these returns a *plausible* point rather than failing, which is why the test
/// suite that matters is `knn_matches_brute_force` on adversarial point sets -- a lattice, an
/// all-duplicate set, a collinear run, a circle, and two tight clusters -- and not a spot check
/// on random noise. The design and the tie-break rule (lower index wins) are settled; what is
/// missing is a version that passes.
/// Is `a` strictly better than `b` under the `(distance, index)` order? The index
/// One node of the k-d tree.
///
/// An **arena** of these, addressed by index, rather than the textbook `2v+1 / 2v+2` heap
/// layout. The heap layout is the source of the first implementation's worst bug: filling a
/// `Vec` by node index with `resize(node + 1, 0)` **truncates** whenever a deep node is written
/// before a shallow one, which an explicit stack guarantees happens -- and a truncated arena
/// leaves holes where whole subtrees were, so the search walks into them and returns a plausible
/// point that is not the nearest.
///
/// An arena has no such failure mode: nodes are appended, so every index handed out is written
/// exactly once and `None` is the only "absent" marker.
#[derive(Clone, Copy, Debug)]
pub struct Node {
    /// Index into the point set.
    pub point: usize,
    /// The axis this node splits on -- the widest at build time, and **stored**, because the
    /// search must use the axis the node was built on rather than one derived from its depth.
    pub axis: usize,
    /// Child node indices, absent at a leaf.
    pub left: Option<usize>,
    pub right: Option<usize>,
}

/// Is `a` strictly better than `b` under the `(distance, index)` order?
///
/// The index tie-break makes this a **total** order, which is what lets the search prune
/// soundly: a subtree whose lower bound merely *equals* the current k-th best can still hold a
/// lower-indexed point at the same distance, so the bound test is `>` and never `>=`.
#[inline]
fn strictly_better(a: (f64, usize), b: (f64, usize)) -> bool {
    match a.0.partial_cmp(&b.0) {
        Some(std::cmp::Ordering::Less) => true,
        Some(std::cmp::Ordering::Greater) => false,
        // Equal distances: the lower index wins. This is the rule the whole module rests on --
        // `computeRegionDistance` counts `length(intersect(knn.i, unique(qout$index[...]))`, and
        // `unique` keeps the first occurrence, so the chosen *order* is observable even when the
        // chosen *set* is not.
        _ => a.1 < b.1,
    }
}

/// A max-heap under [`strictly_better`]: the root is the **worst** of the kept candidates, so a
/// full heap prunes in one comparison and a better candidate replaces the root directly.
///
/// The first implementation pushed unconditionally and then did `swap(0, len - 1); pop()`, which
/// removes `heap[len - 1]` *after* the swap -- i.e. the old root, which for `k == 1` is the old
/// *best*. Every candidate after the first therefore replaced the answer wholesale and "nearest"
/// was simply the last point visited. It returned a plausible point, which is why it survived a
/// spot check and failed only where ties are common.
struct BestHeap {
    v: Vec<(f64, usize)>,
    cap: usize,
}

impl BestHeap {
    fn new(cap: usize) -> Self {
        BestHeap {
            v: Vec::with_capacity(cap),
            cap,
        }
    }
    fn full(&self) -> bool {
        self.v.len() >= self.cap
    }
    fn worst(&self) -> (f64, usize) {
        self.v[0]
    }
    fn offer(&mut self, cand: (f64, usize)) {
        if self.v.len() < self.cap {
            self.v.push(cand);
            let mut i = self.v.len() - 1;
            while i > 0 {
                let parent = (i - 1) / 2;
                if strictly_better(self.v[parent], self.v[i]) {
                    self.v.swap(i, parent);
                    i = parent;
                } else {
                    break;
                }
            }
        } else if strictly_better(cand, self.v[0]) {
            self.v[0] = cand;
            self.sift_down(0);
        }
    }
    fn sift_down(&mut self, mut i: usize) {
        loop {
            let (l, r) = (2 * i + 1, 2 * i + 2);
            let mut worst = i;
            for c in [l, r] {
                if c < self.v.len() && strictly_better(self.v[worst], self.v[c]) {
                    worst = c;
                }
            }
            if worst == i {
                return;
            }
            self.v.swap(i, worst);
            i = worst;
        }
    }
    /// Ascending `(distance, index)`, the order `computeRegionDistance` reads the list in.
    fn into_sorted(mut self) -> Vec<(f64, usize)> {
        self.v.sort_by(|a, b| {
            a.0.partial_cmp(&b.0)
                .unwrap_or(std::cmp::Ordering::Equal)
                .then(a.1.cmp(&b.1))
        });
        self.v
    }
}

/// An exact k-d tree over `d`-dimensional points.
///
/// Built by recursive median splits on the widest axis. A median split is always balanced, so the
/// recursion is `O(log n)` deep even for degenerate inputs -- every point identical, or all
/// collinear -- and an explicit stack is not needed. The split comparator falls back to the point
/// index so the partition is total, and a node of identical points splits evenly rather than
/// degenerating into a chain.
pub struct KdTree {
    pts: Vec<f64>,
    n: usize,
    d: usize,
    root: Option<usize>,
    nodes: Vec<Node>,
}

impl KdTree {
    /// Build over `n` points of `d` dimensions, `coords` in **row-major** order
    /// (`coords[i * d + k]` is point `i`'s coordinate `k`).
    pub fn new(coords: Vec<f64>, n: usize, d: usize) -> Self {
        assert_eq!(coords.len(), n * d, "coordinate count must be n * d");
        let mut idx: Vec<usize> = (0..n).collect();
        let mut nodes: Vec<Node> = Vec::with_capacity(n);
        let root = if n == 0 || d == 0 {
            None
        } else {
            build(&coords, d, &mut idx, 0, n, &mut nodes)
        };
        KdTree {
            pts: coords,
            n,
            d,
            root,
            nodes,
        }
    }

    pub fn len(&self) -> usize {
        self.n
    }

    pub fn is_empty(&self) -> bool {
        self.n == 0
    }

    /// The arena, for the tests. Exposed rather than derived, because a test that reconstructs
    /// the layout from the documented rule cannot catch a build that does not follow the rule.
    pub fn nodes(&self) -> &[Node] {
        &self.nodes
    }

    /// The `k` nearest neighbours of `query`, as `(squared distance, index)` ascending, ties
    /// broken by **lower index**.
    pub fn knn(&self, query: &[f64], k: usize) -> Vec<(f64, usize)> {
        if k == 0 || self.n == 0 || self.d == 0 {
            return Vec::new();
        }
        let k = k.min(self.n);
        let mut heap = BestHeap::new(k);
        if let Some(root) = self.root {
            // `NEG_INFINITY` as the root's bound: nothing may be excluded before a first
            // candidate exists, and a bound of `0` would be wrong for a query far from the
            // origin.
            self.search(root, query, f64::NEG_INFINITY, &mut heap);
        }
        heap.into_sorted()
    }

    /// The single nearest point to `query`, as `(squared distance, index)`. This is the `k = 1`
    /// query `computeRegionDistance` makes for every cell of group *i* into group *j*.
    pub fn nearest(&self, query: &[f64]) -> Option<(f64, usize)> {
        self.knn(query, 1).into_iter().next()
    }

    /// Recursive, with the bound update that makes the whole thing work.
    ///
    /// `lower` is a valid lower bound on the squared distance from `query` to *anything* in this
    /// subtree. Two rules, and getting the second one backwards is what broke the first
    /// implementation:
    ///
    /// * the **far** subtree is constrained by the split plane, so it inherits at least
    ///   `delta^2`;
    /// * the **near** subtree contains the query's own coordinate and inherits its parent's
    ///   bound **unchanged**.
    ///
    /// Tightening the near side too is the symmetric-looking reading, and it makes the near bound
    /// too large -- so the search prunes the subtree that holds the answer.
    fn search(&self, node: usize, q: &[f64], lower: f64, heap: &mut BestHeap) {
        if heap.full() && lower > heap.worst().0 {
            return;
        }
        let nd = self.nodes[node];
        heap.offer((
            squared_distance(&self.pts[nd.point * self.d..nd.point * self.d + self.d], q),
            nd.point,
        ));
        let delta = self.pts[nd.point * self.d + nd.axis] - q[nd.axis];
        // `delta` is the *point*'s coordinate minus the *query*'s, so `delta > 0` puts the split
        // plane to the right of the query: the query is on the left and the **left** child is
        // near. The more common `query - point` spelling gives the opposite, and the search then
        // descends the far subtree first -- still returning a neighbour, just the wrong one.
        let (near, far) = if delta > 0.0 {
            (nd.left, nd.right)
        } else {
            (nd.right, nd.left)
        };
        if let Some(n) = near {
            self.search(n, q, lower, heap);
        }
        if let Some(f) = far {
            self.search(f, q, lower.max(delta * delta), heap);
        }
    }
}

/// Recursive median split, appending to the arena. `order[lo..hi]` is permuted in place.
///
/// Returns `None` for an **empty** range, which the median split produces routinely: with
/// `mid = lo + (hi - lo) / 2` the left child is `[lo, mid)` and the right is `(mid, hi)`, so an
/// even-length range of 2 gives a right child of `(1, 2)` and a left child of `[0, 1)` -- and
/// higher up, `[mid + 1, hi)` is empty whenever `mid + 1 == hi`. A signature returning a bare
/// `usize` cannot express "no node here", and the empty case recurses forever: `hi - lo == 0`
/// is not `== 1`, so it takes the internal path, computes `mid == lo`, and calls itself with the
/// same range. It shows up as a stack overflow, not as a wrong answer, which is at least loud.
fn build(
    coords: &[f64],
    d: usize,
    order: &mut [usize],
    lo: usize,
    hi: usize,
    nodes: &mut Vec<Node>,
) -> Option<usize> {
    if lo >= hi {
        return None;
    }
    let me = nodes.len();
    if hi - lo == 1 {
        nodes.push(Node {
            point: order[lo],
            axis: 0,
            left: None,
            right: None,
        });
        return Some(me);
    }
    // The widest axis of this node's points. Scanning the box at each node is `O(n log n)`
    // overall and makes the tree's quality independent of the input order, which a
    // "cycle through the axes" rule does not.
    let mut axis = 0usize;
    let mut best_spread = f64::NEG_INFINITY;
    for k in 0..d {
        let mut lo_v = f64::INFINITY;
        let mut hi_v = f64::NEG_INFINITY;
        for &i in &order[lo..hi] {
            let v = coords[i * d + k];
            if v < lo_v {
                lo_v = v;
            }
            if v > hi_v {
                hi_v = v;
            }
        }
        let spread = hi_v - lo_v;
        // `>` so the lowest axis wins a tie, which makes the layout a function of the point set
        // alone and not of the scan order.
        if spread > best_spread {
            best_spread = spread;
            axis = k;
        }
    }
    order[lo..hi].sort_unstable_by(|&a, &b| {
        coords[a * d + axis]
            .partial_cmp(&coords[b * d + axis])
            .unwrap_or(std::cmp::Ordering::Equal)
            .then_with(|| a.cmp(&b))
    });
    // The **upper** median of an even-length range. Either median is defensible; this one is
    // stated because a reader should not have to guess which is in force.
    let mid = lo + (hi - lo) / 2;
    let point = order[mid];
    nodes.push(Node {
        point,
        axis,
        left: None,
        right: None,
    });
    let left = build(coords, d, order, lo, mid, nodes);
    let right = build(coords, d, order, mid + 1, hi, nodes);
    nodes[me].left = left;
    nodes[me].right = right;
    Some(me)
}

/// Exhaustive search under the same `(distance, index)` order, for the tests.
pub fn brute_force_knn(pts: &[f64], n: usize, d: usize, q: &[f64], k: usize) -> Vec<(f64, usize)> {
    let mut all: Vec<(f64, usize)> = (0..n)
        .map(|i| (squared_distance(&pts[i * d..i * d + d], q), i))
        .collect();
    all.sort_by(|a, b| {
        a.0.partial_cmp(&b.0)
            .unwrap_or(std::cmp::Ordering::Equal)
            .then(a.1.cmp(&b.1))
    });
    all.truncate(k.min(n));
    all
}

// ---------------------------------------------------------------- computeRegionDistance

/// The inputs to [`region_distance`], already resolved.
///
/// `ratio` and `tol` are **per sample**, indexed by the sample's position in
/// `samples_levels`: upstream writes `qout$distance * ratio[k]` and
/// `qout$distance - interaction.range < tol[k]`, so a per-sample vector is not a convenience but
/// the actual contract. A scalar in R recycles against the sample index, and the shim therefore
/// takes vectors and holds the caller to them.
#[derive(Clone, Debug)]
pub struct RegionInput<'a> {
    /// `n * d`, row-major.
    pub coords: &'a [f64],
    pub n: usize,
    pub d: usize,
    /// `nlevels(group)`, and the factor's *level* order.
    pub group_levels: &'a [String],
    /// 0-based group index per cell, into `group_levels`.
    pub group_index: &'a [usize],
    pub samples_levels: &'a [String],
    /// 0-based sample index per cell.
    pub sample_index: &'a [usize],
    pub interaction_range: Option<f64>,
    pub ratio: &'a [f64],
    pub tol: &'a [f64],
    pub k_min: usize,
    pub contact_dependent: bool,
    pub contact_range: Option<f64>,
    pub contact_knn_k: Option<usize>,
    pub do_symmetric: bool,
}

/// What `computeRegionDistance` returns: `d.spatial` and the chosen `adj.contact`.
#[derive(Clone, Debug, PartialEq)]
pub struct RegionOutput {
    /// `k x k`, row-major, in `group_levels` order. `NaN` where the groups are not adjacent.
    pub d_spatial: Vec<f64>,
    /// `k x k`, 0/1.
    pub adj_contact: Vec<f64>,
    /// Kept separate from `adj_contact` because upstream *discards* `adj.contact.knn` when
    /// `contact.knn.k` is supplied -- `if (length(contact.knn.k) > 0) adj.contact = adj.contact.knn`
    /// -- and a caller who wants to compare the two needs the value that was dropped.
    pub adj_contact_knn: Vec<f64>,
    pub group_levels: Vec<String>,
}

/// Extend `v` to `ns` entries the way R's recycling does for `v[k]`.
///
/// A length-1 vector repeats, which is what upstream's own documentation implies when it calls
/// `ratio` and `tol` "the ratio and tolerance for each sample" and then indexes them by sample.
/// Any other short vector is an error rather than a wrap: R would return `NA` for the missing
/// entries and silently produce `NaN` distances, and a `NaN` that means "you passed the wrong
/// length" is indistinguishable from a `NaN` that means "these groups are not adjacent".
fn recycle(v: &[f64], ns: usize) -> Result<Vec<f64>, RegionError> {
    // An empty vector is upstream's `NULL`. `qout$distance * NULL` is `numeric(0)`, so
    // `FunMean(numeric(0))` is `NaN` and both adjacency tests select nothing. Returning empty makes
    // the per-sample loop body never execute, which leaves `d.spatial` at its initialised `NaN` and
    // the adjacencies at 0 -- the same values, reached without special-casing `mean` of nothing.
    // This is the *default* signature (`ratio = NULL, tol = NULL`), so it is not a corner case.
    if v.is_empty() || v.len() == ns {
        return Ok(v.to_vec());
    }
    if v.len() == 1 {
        return Ok(vec![v[0]; ns]);
    }
    Err(RegionError::RatioOrTolTooShort { need: ns })
}

/// `computeRegionDistance`, with the exact k-d tree in place of `AnnoyParam`.
///
/// The arithmetic is upstream's, in upstream's order. The R-isms that matter, all of which the
/// corpus in `tests/parity/gen_region_golden.R` pins:
///
/// * **`numCluster` is `nlevels(group)`, but `level.use` drops the levels absent from the
///   data.** The arrays are sized by the former and indexed by the latter, so a level with no
///   cells leaves a row and column that are never written: `d.spatial` stays `NaN` and the
///   `adj.*` stay `0`. Those survive the per-sample mean -- `mean(x, na.rm = TRUE)` over all-`NaN`
///   is `NaN`, and the `adj` means have no `na.rm` but their entries are `0`/`1` -- and they are
///   what the `> 0` binarisation and the symmetrisation then act on.
/// * **The per-sample `mean` of a 0/1 array is a fraction, and `adj[adj > 0] <- 1` turns it back
///   into 0/1.** So a pair adjacent in *any* sample is adjacent overall, and the fraction is
///   discarded.
/// * **Symmetrisation is `adj * t(adj)`**, so a pair with a zero in either direction becomes zero
///   in both. `d.spatial` is symmetrised as `(d + t(d)) / 2` instead, and that **propagates
///   `NaN`**: a level absent from every sample leaves one direction `NaN`, and the mean is `NaN`
///   in both.
/// * **`adj.spatial[adj.spatial == 0] <- NaN` then `d.spatial <- d.spatial * adj.spatial`,** so
///   `d.spatial` is `NaN` exactly where the groups are not adjacent. The caller turns those into
///   zeros: `P.spatial[is.na(d.spatial)] <- 0`.
/// * **`k.min.contact` is `-1`** when `contact.knn.k` is `NULL`, and `length(...) >= -1` is always
///   true, so `adj.contact.knn` is all ones. It is the documented way upstream disables the test
///   rather than a bug, and reproducing it as "skip the comparison" would give the same answer
///   only by accident.
/// * **`if (contact.dependent) { if (is.null(contact.range) & is.null(contact.knn.k)) stop(...) }`.**
///   The error fires when `contact.dependent` is `TRUE` and *neither* range is given. With
///   `contact.dependent = FALSE`, `contact.range` is overwritten with `10000` -- which passes
///   every query through the range test but does not by itself make any `adj.contact` entry 1. With
///   `contact.range = NULL`, the normal case when `contact.knn.k` is given instead, the test
///   `qout$distance - NULL < tol[k]` evaluates on `numeric(0)` and selects nothing, so
///   `adj.contact` is 0 throughout and only the `adj.contact.knn` swap can produce a 1.
pub fn region_distance(inp: &RegionInput<'_>) -> Result<RegionOutput, RegionError> {
    let k = inp.group_levels.len();
    let ns = inp.samples_levels.len();
    if ns == 0 || k == 0 {
        return Err(RegionError::NoSamplesOrGroups);
    }
    // R *recycles* `ratio[k]` and `tol[k]` over the sample index, so a length-1 vector is the
    // normal way to pass one value for every sample. The corpus has a two-sample fixture with a
    // length-1 `ratio` precisely because upstream accepts it, and indexing `inp.ratio[s]` here
    // panics instead. Recycled here, once, so the hot loop indexes like upstream's.
    let ratio = recycle(inp.ratio, ns)?;
    let tol = recycle(inp.tol, ns)?;
    if inp.contact_dependent && inp.contact_range.is_none() && inp.contact_knn_k.is_none() {
        return Err(RegionError::NeedsContactRangeOrKnn);
    }
    // `level.use <- levels(group); level.use <- level.use[level.use %in% unique(group)]`.
    //
    // This is a **compacted** list, and upstream then walks `for (i in 1:numCluster)` -- the *full*
    // level count -- reading the level name out of `level.use[i]`. The loop counter and the list
    // index are two different counters over two different things, and conflating them is the whole
    // bug:
    //
    //   * the **slot written** is the loop counter (1-based in R, so `i - 1` here), and
    //   * the **level whose cells are read** is the `i`-th element of the compacted list.
    //
    // With `levels = c("C","A","B")` and only A and B populated, `level.use` is `c("A","B")`, so the
    // loop puts **A's row in slot 1**, and `rownames(d.spatial) <- levels(group)` then labels that
    // row "C". The names are off by one level and the last slot is never written at all.
    //
    // Indexing the compacted list by the level's own position instead -- the obvious reading, and
    // what this used to do -- gives the transpose: `NaN` on the *first* level's row instead of the
    // last, with every value correct.
    let level_use: Vec<usize> = (0..k).filter(|g| inp.group_index.contains(g)).collect();

    // `nn.ranked`: per sample, the k nearest neighbours of every cell within that sample.
    // `matrix(1, nrow, 1)` when `contact.knn.k` is NULL -- a literal 1, not a neighbour, and
    // `unique(as.vector(...))` of it is just `1`.
    let kmin_contact: i64 = if inp.contact_knn_k.is_some() {
        inp.k_min as i64
    } else {
        -1
    };
    let nn_ranked: Vec<Vec<usize>> = match inp.contact_knn_k {
        None => vec![vec![0usize; inp.n]; ns],
        Some(kk) => (0..ns)
            .map(|s| {
                let idx: Vec<usize> = (0..inp.n).filter(|&i| inp.sample_index[i] == s).collect();
                let coords: Vec<f64> = idx
                    .iter()
                    .flat_map(|&i| inp.coords[i * inp.d..i * inp.d + inp.d].iter().copied())
                    .collect();
                let tree = KdTree::new(coords, idx.len(), inp.d);
                let mut out = vec![0usize; idx.len() * kk];
                for (r, &i) in idx.iter().enumerate() {
                    let q = &inp.coords[i * inp.d..i * inp.d + inp.d];
                    for (c, (_, j)) in tree.knn(q, kk).into_iter().enumerate() {
                        out[r * kk + c] = j;
                    }
                }
                // `nn.ranked[idx.k, ] <- my.knn$index` stores indices into *that sample's*
                // subset, which is what `intersect(knn.i, qout$index[idx])` compares against --
                // and `qout$index` is likewise an index into `idx.j`. Both are sample-local.
                out
            })
            .collect(),
    };

    // `contact.range <- 10000` when not contact-dependent. The upstream comment says "this
    // produces adj.contact with all elements being 1", and that is true of the *range test*, not
    // of the result: every query then passes the test, but the entry is still
    // `length(unique(qout$index[selected rows])) >= k.min`, a count of **distinct** target cells.
    // Three source cells that all pick the same nearest target give 1, so with `k.min = 2` the
    // entry is 0. `not_contact_dependent` is the fixture that pins this.
    let contact_range = if inp.contact_dependent {
        inp.contact_range
    } else {
        Some(10000.0)
    };

    let mut d_spatial = vec![f64::NAN; k * k * ns];
    let mut adj_spatial = vec![0.0f64; k * k * ns];
    let mut adj_contact = vec![0.0f64; k * k * ns];
    let mut adj_knn = vec![0.0f64; k * k * ns];

    for (s, ratio_s) in ratio.iter().enumerate() {
        let tol_s = tol[s];
        let idx_k: Vec<usize> = (0..inp.n).filter(|&i| inp.sample_index[i] == s).collect();
        if idx_k.is_empty() {
            continue;
        }
        for i_slot in 0..k {
            for j_slot in 0..k {
                // `level.use[i]` past its end is `NA` in R, so `which(group == NA & idx.k)` is empty
                // and upstream `next`s. Nothing is written to that slot.
                let (Some(&gi), Some(&gj)) = (level_use.get(i_slot), level_use.get(j_slot)) else {
                    continue;
                };
                let idx_i: Vec<usize> = idx_k
                    .iter()
                    .copied()
                    .filter(|&i| inp.group_index[i] == gi)
                    .collect();
                let idx_j: Vec<usize> = idx_k
                    .iter()
                    .copied()
                    .filter(|&i| inp.group_index[i] == gj)
                    .collect();
                if idx_i.is_empty() || idx_j.is_empty() {
                    continue;
                }
                // For each cell of group *i*, the 1-NN in group *j*.
                let coords_j: Vec<f64> = idx_j
                    .iter()
                    .flat_map(|&i| inp.coords[i * inp.d..i * inp.d + inp.d].iter().copied())
                    .collect();
                let tree_j = KdTree::new(coords_j, idx_j.len(), inp.d);
                let mut dist: Vec<f64> = Vec::with_capacity(idx_i.len());
                let mut near: Vec<usize> = Vec::with_capacity(idx_i.len());
                for &i in &idx_i {
                    let q = &inp.coords[i * inp.d..i * inp.d + inp.d];
                    // `qout$distance` is the *Euclidean* distance, not squared, and
                    // `qout$index` is a position in `idx_j`.
                    let (sq, j) = tree_j.nearest(q).expect("group j is non-empty");
                    dist.push(sq.sqrt());
                    near.push(j);
                }
                // `qout$distance <- qout$distance * ratio[k]` -- the scaling happens *before*
                // every threshold test, so `interaction.range` and `contact.range` are compared
                // in micrometres.
                for v in dist.iter_mut() {
                    *v *= *ratio_s;
                }
                let irange = inp.interaction_range.unwrap_or(f64::INFINITY);
                let in_long: Vec<usize> = (0..dist.len())
                    .filter(|&r| dist[r] - irange < tol_s)
                    .collect();
                // `idx2 <- qout$distance - contact.range < tol[k]`. With `contact.range = NULL`,
                // which is the *normal* case when `contact.knn.k` is supplied instead, `dist - NULL`
                // is `numeric(0)` and the comparison selects nothing. So `in_contact` is empty and
                // `adj.contact` is 0 for the whole run -- it does not mean "skip this pair".
                // Skipping instead left `adj.spatial` and `d.spatial` unwritten, i.e. NaN, for
                // every pair, which is how every `contact.knn.k` fixture came back all-NaN.
                let in_contact: Vec<usize> = match contact_range {
                    Some(cr) => (0..dist.len()).filter(|&r| dist[r] - cr < tol_s).collect(),
                    None => Vec::new(),
                };
                // Row-major, matching the `k x k` matrix the caller reads. The per-sample arrays
                // are only indexed here and in `merge_samples_mean`, so the two have to agree with
                // each other -- but they also have to agree with the *output*, and the corpus, like
                // every R caller, reads `[i, j]` as row `i`. All three outputs are symmetrised
                // before they leave, so a column-major layout looked right on every fixture that
                // has one; it would not have on a caller reading an unsymmetrised matrix.
                let at = i_slot * k + j_slot + k * k * s;
                adj_spatial[at] = if unique_sorted(&near, &in_long).len() as i64 >= inp.k_min as i64
                {
                    1.0
                } else {
                    0.0
                };
                adj_contact[at] =
                    if unique_sorted(&near, &in_contact).len() as i64 >= inp.k_min as i64 {
                        1.0
                    } else {
                        0.0
                    };
                // `knn.i <- unique(as.vector(nn.ranked[idx.i, ]))`, then
                // `length(intersect(knn.i, unique(qout$index[idx])))`. `nn.ranked`'s columns are
                // sample-local positions, as is `qout$index`, so the two are comparable.
                let kk = inp.contact_knn_k.unwrap_or(1);
                let mut knn_i: Vec<usize> = Vec::new();
                for r in 0..idx_i.len() {
                    let local = idx_k.iter().position(|&x| x == idx_i[r]).unwrap();
                    for c in 0..kk {
                        knn_i.push(nn_ranked[s][local * kk + c]);
                    }
                }
                knn_i.sort_unstable();
                knn_i.dedup();
                let long_set = unique_sorted(&near, &in_long);
                let inter = knn_i.iter().filter(|x| long_set.contains(x)).count() as i64;
                adj_knn[at] = if inter >= kmin_contact { 1.0 } else { 0.0 };
                d_spatial[at] = trimmed_mean(&dist, 0.1, true);
            }
        }
    }

    // `apply(., c(1,2), ...)` over the sample axis.
    let d_spatial = merge_samples_mean(d_spatial, k, ns, true);
    let mut adj_spatial = merge_samples_mean(adj_spatial, k, ns, false);
    let mut adj_contact = merge_samples_mean(adj_contact, k, ns, false);
    let mut adj_knn = merge_samples_mean(adj_knn, k, ns, false);

    // `adj[adj > 0] <- 1`: a pair adjacent in *any* sample is adjacent, and the fraction is
    // discarded. Applied before the symmetrisation, so it is that 0/1 matrix which is
    // multiplied by its transpose.
    for v in adj_spatial
        .iter_mut()
        .chain(adj_contact.iter_mut())
        .chain(adj_knn.iter_mut())
    {
        if *v > 0.0 {
            *v = 1.0;
        }
    }
    let mut d_spatial = d_spatial;
    if inp.do_symmetric {
        adj_spatial = symmetrize_and(adj_spatial, k);
        adj_contact = symmetrize_and(adj_contact, k);
        adj_knn = symmetrize_and(adj_knn, k);
    }
    // `d.spatial <- (d.spatial + t(d.spatial)) / 2` is **outside** the `if (do.symmetric)` block in
    // upstream, so it happens either way. Read as a whole that reads like an oversight, and it is
    // tempting to "fix" it -- but the `not_symmetric` fixture pins it: with `do.symmetric = FALSE`
    // the returned `adj.contact` is genuinely asymmetric while `d.spatial` is symmetric. Moving
    // this inside the guard, as I first had it, produced a `d.spatial` that differed from upstream
    // in three of nine entries on that one fixture.
    //
    // It also **propagates NaN**: a level absent from every sample leaves one direction NaN and the
    // mean is NaN in both.
    d_spatial = symmetrize_mean(d_spatial, k);
    // `adj.spatial[adj.spatial == 0] <- NaN; d.spatial <- d.spatial * adj.spatial`.
    for i in 0..k * k {
        if adj_spatial[i] == 0.0 {
            adj_spatial[i] = f64::NAN;
        }
        d_spatial[i] *= adj_spatial[i];
    }
    // `if (length(contact.knn.k) > 0) adj.contact = adj.contact.knn`. `length(NULL)` is 0, so a
    // NULL `contact.knn.k` leaves `adj.contact` alone; `length(NA)` is 1, which is a different
    // case entirely.
    let chosen = if inp.contact_knn_k.map(|k| k > 0).unwrap_or(false) {
        adj_knn.clone()
    } else {
        adj_contact
    };
    Ok(RegionOutput {
        d_spatial,
        adj_contact: chosen,
        adj_contact_knn: adj_knn,
        group_levels: inp.group_levels.to_vec(),
    })
}

/// `unique(qout$index[idx])`: the distinct neighbour positions among the selected rows, ascending.
///
/// Ascending because `intersect` and `length` do not care, but the *order* is what
/// `unique()` would return and the corpus compares it as text.
fn unique_sorted(all: &[usize], sel: &[usize]) -> Vec<usize> {
    let mut v: Vec<usize> = sel.iter().map(|&r| all[r]).collect();
    v.sort_unstable();
    v.dedup();
    v
}

/// `apply(x, c(1, 2), function(v) mean(v, na.rm = ...))` over a `k x k x ns` array.
fn merge_samples_mean(x: Vec<f64>, k: usize, ns: usize, na_rm: bool) -> Vec<f64> {
    let mut out = vec![0.0f64; k * k];
    for i in 0..k {
        for j in 0..k {
            let col: Vec<f64> = (0..ns).map(|s| x[i * k + j + k * k * s]).collect();
            out[i * k + j] = if na_rm {
                trimmed_mean(&col, 0.0, true)
            } else {
                plain_mean(&col)
            };
        }
    }
    out
}

fn plain_mean(x: &[f64]) -> f64 {
    if x.is_empty() {
        return f64::NAN;
    }
    if x.iter().any(|v| v.is_nan()) {
        return f64::NAN;
    }
    let mut s = crate::longdouble::F80::from_f64(0.0);
    for v in x {
        s = s.add(crate::longdouble::F80::from_f64(*v));
    }
    s.div_int(x.len() as u64).to_f64()
}

/// `m * t(m)`.
fn symmetrize_and(m: Vec<f64>, k: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; k * k];
    for i in 0..k {
        for j in 0..k {
            out[i + k * j] = m[i + k * j] * m[j + k * i];
        }
    }
    out
}

/// `(m + t(m)) / 2`, which propagates `NaN`.
fn symmetrize_mean(m: Vec<f64>, k: usize) -> Vec<f64> {
    let mut out = vec![0.0f64; k * k];
    for i in 0..k {
        for j in 0..k {
            out[i + k * j] = (m[i + k * j] + m[j + k * i]) / 2.0;
        }
    }
    out
}

/// The errors `computeRegionDistance` raises.
#[derive(Clone, Debug, PartialEq, Eq)]
pub enum RegionError {
    /// `stop("Please check the documentation of `computeCommunProb` and provide the value of either `contact.range` or `contact.knn.k`")`
    NeedsContactRangeOrKnn,
    /// `ratio` or `tol` has some length other than 1 or `nlevels(samples)`. A length-1 vector is
    /// *not* this error: R recycles it, so it is the normal way to pass one value for every
    /// sample. Anything else short is refused, where upstream would read `NA` from the recycled
    /// vector and silently produce `NaN` distances.
    RatioOrTolTooShort { need: usize },
    /// No samples or no groups.
    NoSamplesOrGroups,
}

impl std::fmt::Display for RegionError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            RegionError::NeedsContactRangeOrKnn => write!(
                f,
                "Please check the documentation of `computeCommunProb` and provide the value of \
                 either `contact.range` or `contact.knn.k`"
            ),
            RegionError::RatioOrTolTooShort { need } => {
                write!(
                    f,
                    "`ratio` and `tol` must have one entry per sample ({need})"
                )
            }
            RegionError::NoSamplesOrGroups => write!(f, "no samples or no cell groups"),
        }
    }
}

impl std::error::Error for RegionError {}
