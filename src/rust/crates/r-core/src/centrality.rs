//! The deterministic half of upstream's `computeCentralityLocal`.
//!
//! Upstream computes eleven measures per pathway network. Four are iterative solvers whose output
//! is not even stable across runs -- `hub_score` (deprecated in igraph 2.0.3), `authority_score`,
//! `eigen_centrality` (ARPACK) and `page_rank` (PRPACK) all differ run to run on identical input,
//! measured in `check_centrality.R`'s probe section. Bit-parity with a nondeterministic oracle is
//! meaningless, so those four stay in R and call igraph directly. What lives here is everything
//! whose output is a pure function of the input:
//!
//! * `outdeg_unweighted` / `indeg_unweighted` -- integer counts,
//! * `outdeg` / `indeg` -- weighted sums,
//! * `betweenness` -- Dijkstra + Brandes over reciprocal weights.
//!
//! `flowbet` / `infocent` are `sna` calls in a `tryCatch` that returns zeros; `sna` is not
//! installed, so both sides return zeros and the shim mirrors the `tryCatch` structure.
//!
//! Every numeric rule below was read out of the igraph C sources vendored for this exact purpose
//! and then verified empirically, because the installed igraph 2.3.4 differs from main-branch
//! source in observable ways (different error strings, different check order). Where the two
//! disagree, the measurement wins and the comment says so:
//!
//! * Edges are created scanning the matrix **column-major** (`j` outer, `i` inner;
//!   `igraph_i_weighted_adjacency_directed`), skipping `M != 0.0` entries -- so `-0.0` is
//!   dropped and `NA`/`NaN` never get this far (the constructor errors first, exactly as
//!   `graph_from_adjacency_matrix` does).
//! * `strength` iterates **edges in edge-ID order** accumulating into per-vertex totals
//!   (`strength_all`), i.e. plain sequential `f64` `+=` -- *not* R's long-double `sum()`.
//!   Verified 4669/4669 against installed igraph; `rowSums` agrees only 3802/4669.
//! * Dijkstra stores distances **plus one** (`dist[source] = 1.0`, `0.0` = infinity), keys the
//!   heap on negated distances, compares with `igraph_cmp_epsilon` at `eps = 1e-10`, and
//!   accumulates path counts (`nrgeo`) as `f64`.
//! * The 2-way heap breaks ties by exact push/sink/shift code replicated below: `shift_up` swaps
//!   on `>=` (ties rise), `sink` swaps only on strict `<` preferring the left child. Any other
//!   tie rule changes pop order, which changes rounding of the accumulated scores.
//! * `betweenness` normalisation factor is `1.0` for directed + unnormalised -- an exact no-op,
//!   skipped rather than applied.
//! * Weight checks, in this order (probed, because the sources disagree): any `NaN` reciprocal
//!   errors first regardless of position; then non-positive; then `<= 1e-10` warns. The messages
//!   are igraph 2.3.4's verbatim, `Source:` suffix included -- that suffix is part of
//!   `conditionMessage()`, so the gate compares it.
//!
//! Two deliberate divergences from a naive port, both verified:
//!
//! * Loops are **included once** in degrees and strength (the probe graphs keep them) but
//!   **excluded** from the betweenness traversal (`inclist` with `IGRAPH_NO_LOOPS`).
//! * An empty graph (no edges) returns zeros with no warning and no error -- probed, not reasoned.

/// `Cannot create a graph object because the adjacency matrix contains NAs.`
///
/// igraph's message verbatim, including the trailing period and the absence of a `Source:`
/// suffix (unlike the weight messages below, which have one -- measured, not explained).
pub const GRAPH_NA_ERROR: &str =
    "Cannot create a graph object because the adjacency matrix contains NAs.";
/// `Weight vector must not contain NaN values. Invalid value\nSource: centrality/betweenness.c:439`
pub const WEIGHT_NAN_ERROR: &str =
    "Weight vector must not contain NaN values. Invalid value\nSource: centrality/betweenness.c:439";
/// `Weight vector must be positive. Invalid value\nSource: centrality/betweenness.c:437`
pub const WEIGHT_POSITIVE_ERROR: &str =
    "Weight vector must be positive. Invalid value\nSource: centrality/betweenness.c:437";
/// `Some weights are smaller than epsilon, ...\nSource: centrality/betweenness.c:441`
pub const WEIGHT_EPS_WARNING: &str =
    "Some weights are smaller than epsilon, calculations may suffer from numerical precision issues.\nSource: centrality/betweenness.c:441";
/// `IGRAPH_SHORTEST_PATH_EPSILON`. The warning fires at `minweight <= eps` -- probed at 1e-10
/// (warns) versus 1.0000001e-10 (silent).
pub const SHORTEST_PATH_EPSILON: f64 = 1e-10;

#[derive(Debug, Clone, PartialEq)]
pub enum CentralityError {
    /// Any `NA_real_` or `NaN` in the adjacency matrix.
    NonFiniteAdjacency,
    /// A `NaN` reciprocal weight.
    NanWeight,
    /// A non-positive reciprocal weight.
    NonPositiveWeight,
}

impl std::fmt::Display for CentralityError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            CentralityError::NonFiniteAdjacency => f.write_str(GRAPH_NA_ERROR),
            CentralityError::NanWeight => f.write_str(WEIGHT_NAN_ERROR),
            CentralityError::NonPositiveWeight => f.write_str(WEIGHT_POSITIVE_ERROR),
        }
    }
}

/// `igraph_cmp_epsilon`, verbatim from `src/math/utils.c`.
///
/// Returns 0 / negative / positive. `a == b` first (handles infinities); then the zero-or-
/// subnormal branch, the overflow branch, and the relative-error branch, in that order --
/// the order matters because more than one can apply.
fn cmp_epsilon(a: f64, b: f64, eps: f64) -> i32 {
    if a == b {
        return 0;
    }
    let diff = a - b;
    let abs_diff = diff.abs();
    let sum = a.abs() + b.abs();
    if a == 0.0 || b == 0.0 || sum < f64::MIN_POSITIVE {
        if abs_diff < eps * f64::MIN_POSITIVE {
            0
        } else if diff < 0.0 {
            -1
        } else {
            1
        }
    } else if !sum.is_finite() {
        if abs_diff < eps * a.abs() + eps * b.abs() {
            0
        } else if diff < 0.0 {
            -1
        } else {
            1
        }
    } else if abs_diff / sum < eps {
        0
    } else if diff < 0.0 {
        -1
    } else {
        1
    }
}

/// One directed edge, in creation order. Creation scans the matrix column-major (`j` outer,
/// `i` inner) skipping `M != 0.0`, so edge IDs ascend fastest in the target index.
#[derive(Debug, Clone, Copy)]
struct Edge {
    from: usize,
    to: usize,
    weight: f64,
}

/// The graph as igraph sees it: edges in creation order. Per-vertex neighbour lists are built
/// where they are needed (with loops for strength, without for traversal), because the two have
/// different membership and sharing one list would conflate them.
#[derive(Debug)]
struct Graph {
    edges: Vec<Edge>,
}

impl Graph {
    fn build(net: &[f64], n: usize) -> Result<Self, CentralityError> {
        assert_eq!(net.len(), n * n, "centrality: net must be k*k");
        // R's `any(is.na())` catches `NaN` too, and igraph reports both the same way before any
        // edge exists -- so this fires before zero-dropping, loop handling, and everything else,
        // exactly like `graph_from_adjacency_matrix`.
        if net.iter().any(|v| v.is_nan()) {
            return Err(CentralityError::NonFiniteAdjacency);
        }
        let mut edges = Vec::new();
        for j in 0..n {
            for i in 0..n {
                let w = net[i * n + j];
                // `M != 0.0`: `-0.0` is dropped, infinities are kept (probed: `Inf` is an edge).
                if w != 0.0 {
                    edges.push(Edge {
                        from: i,
                        to: j,
                        weight: w,
                    });
                }
            }
        }
        Ok(Graph { edges })
    }

    /// Out-neighbour edge IDs per vertex, ascending -- what `igraph_inclist_init` yields, since
    /// it walks edges `0..E-1`. Loops included or excluded by the caller, because strength counts
    /// them and the traversal does not.
    fn out_lists(&self, n: usize, skip_loops: bool) -> Vec<Vec<usize>> {
        let mut out: Vec<Vec<usize>> = vec![Vec::new(); n];
        for (e, edge) in self.edges.iter().enumerate() {
            if skip_loops && edge.from == edge.to {
                continue;
            }
            out[edge.from].push(e);
        }
        out
    }
}

/// igraph's 2-way heap, replicated operation by operation (`src/core/indheap.c`).
///
/// Max-heap on the stored keys; betweenness stores *negated* distances so the maximum key is the
/// minimum distance. 0-based array with `PARENT(x) = ((x)+1)/2-1`. The tie rules are the whole
/// point and are stated on each operation: get one wrong and pop order changes, which changes the
/// rounding of every accumulated score without changing any parent *set*.
struct TwoHeap {
    /// Heap keys in array order.
    data: Vec<f64>,
    /// Vertex index per heap position.
    index: Vec<usize>,
    /// Heap position per vertex (`pos[v]`), or `None` when absent.
    pos: Vec<Option<usize>>,
}

impl TwoHeap {
    fn new(n: usize) -> Self {
        TwoHeap {
            data: Vec::new(),
            index: Vec::new(),
            pos: vec![None; n],
        }
    }

    fn empty(&self) -> bool {
        self.data.is_empty()
    }

    fn switch(&mut self, e1: usize, e2: usize) {
        self.data.swap(e1, e2);
        self.index.swap(e1, e2);
        let v1 = self.index[e1];
        let v2 = self.index[e2];
        self.pos[v1] = Some(e1);
        self.pos[v2] = Some(e2);
    }

    /// `shift_up`: swaps while `data[elem] >= data[parent]` -- ties RISE (the `else` arm runs on
    /// `==`, since the stop condition is strict `<`). A `<`-only version would leave ties in
    /// place and pop them in a different order.
    fn shift_up(&mut self, mut elem: usize) {
        loop {
            if elem == 0 {
                return;
            }
            // `PARENT(x) = ((x)+1)/2-1`, transcribed rather than simplified: clippy suggests
            // `div_ceil`, but `div_ceil(1, 2) == 1` while `PARENT(1) == 0`, so the suggestion
            // computes a different parent and the heap would pop in a different order.
            #[allow(clippy::manual_div_ceil)]
            let parent = (elem + 1) / 2 - 1;
            if self.data[elem] < self.data[parent] {
                return;
            }
            self.switch(elem, parent);
            elem = parent;
        }
    }

    /// `sink`: descends toward the larger child, preferring LEFT on ties (`>=`), and swaps only
    /// on strict `<`. Ties therefore stay put -- the mirror image of `shift_up`, and both halves
    /// are needed to reproduce pop order.
    fn sink(&mut self, mut head: usize) {
        loop {
            let size = self.data.len();
            let left = (head + 1) * 2 - 1;
            let right = (head + 1) * 2;
            if left >= size {
                return;
            }
            let child = if right == size || self.data[left] >= self.data[right] {
                left
            } else {
                right
            };
            if self.data[head] < self.data[child] {
                self.switch(head, child);
                head = child;
            } else {
                return;
            }
        }
    }

    fn push(&mut self, idx: usize, key: f64) {
        let at = self.data.len();
        self.data.push(key);
        self.index.push(idx);
        self.pos[idx] = Some(at);
        self.shift_up(at);
    }

    fn max_index(&self) -> usize {
        self.index[0]
    }

    /// `delete_max`: swap root with last, pop, mark absent, sink. Returns the old root key.
    fn delete_max(&mut self) -> f64 {
        let tmp = self.data[0];
        let tmpidx = self.index[0];
        let last = self.data.len() - 1;
        self.switch(0, last);
        self.data.pop();
        self.index.pop();
        self.pos[tmpidx] = None;
        if !self.data.is_empty() {
            self.sink(0);
        }
        tmp
    }

    /// `modify`: overwrite in place, then sink AND shift up at the position, in that order.
    /// Both run unconditionally -- `modify` does not compare old against new first.
    fn modify(&mut self, idx: usize, key: f64) {
        let at = self.pos[idx].expect("modify of an absent heap element");
        self.data[at] = key;
        self.sink(at);
        self.shift_up(at);
    }
}

/// The five deterministic measures, in upstream's `centr` order for the ones that live here.
#[derive(Debug, Clone)]
pub struct Deterministic {
    pub outdeg_unweighted: Vec<f64>,
    pub indeg_unweighted: Vec<f64>,
    pub outdeg: Vec<f64>,
    pub indeg: Vec<f64>,
    pub betweenness: Vec<f64>,
}

/// `true` iff any reciprocal weight is at or below the comparison tolerance -- the condition for
/// igraph's warning, checked on the reciprocals (what `E(G)$weight` holds when `betweenness`
/// runs), not on the raw probabilities.
pub fn tiny_weights_present(net: &[f64], n: usize) -> bool {
    // No edges: probed -- igraph returns zeros with no warning and no error, so the weight
    // checks (a minimum over an empty vector) never run.
    let mut any_edge = false;
    let mut tiny = false;
    for j in 0..n {
        for i in 0..n {
            let w = net[i * n + j];
            if w != 0.0 {
                any_edge = true;
                // `1.0 / w`: the exact IEEE division R performs elementwise for
                // `E(G)$weight <- 1 / E(G)$weight`, so the comparison below sees the same bits.
                let r = 1.0 / w;
                // NaN and non-positive reciprocals error before the warning is reachable; this
                // predicate is only about the warning, so anything else is ignored here.
                if r.is_finite() && r <= SHORTEST_PATH_EPSILON {
                    tiny = true;
                }
            }
        }
    }
    any_edge && tiny
}

/// Degrees and strengths. Plain sequential `f64` accumulation in edge-ID order -- verified
/// 4669/4669 against installed igraph, where R's long-double `rowSums` agrees only 3802/4669.
/// Loops count once (probed with nonzero diagonals); `-0.0` entries are not edges.
pub fn degrees_strength(net: &[f64], n: usize) -> Result<Deterministic, CentralityError> {
    let g = Graph::build(net, n)?;
    let mut outdeg_unweighted = vec![0.0f64; n];
    let mut indeg_unweighted = vec![0.0f64; n];
    let mut outdeg = vec![0.0f64; n];
    let mut indeg = vec![0.0f64; n];
    // One pass over edges in ID order, accumulating into per-vertex totals -- `strength_all`.
    // Separate out/in loops, out first, exactly as the source has them. The unweighted counts are
    // `rowSums(net > 0)` / `colSums(net > 0)` -- strictly positive, not merely nonzero: a negative
    // entry *creates an edge* (`M != 0.0` in the builder) but is *not counted* (`-0.1 > 0` is
    // FALSE). Counting `!= 0.0` passed every corpus case and was wrong, because no corpus matrix
    // had a negative entry in a degree position before the neg_weight cases.
    for e in &g.edges {
        if e.weight > 0.0 {
            outdeg_unweighted[e.from] += 1.0;
        }
        outdeg[e.from] += e.weight;
    }
    for e in &g.edges {
        if e.weight > 0.0 {
            indeg_unweighted[e.to] += 1.0;
        }
        indeg[e.to] += e.weight;
    }
    Ok(Deterministic {
        outdeg_unweighted,
        indeg_unweighted,
        outdeg,
        indeg,
        betweenness: vec![0.0; n],
    })
}

/// Weighted directed betweenness, replicating `igraph_betweenness_cutoff` with `cutoff = -1`
/// (exact), `normalized = false`, `directed = true`.
///
/// `net` is the K×K pathway probability slice in row-major order. Reciprocals, checks, Dijkstra
/// with `dist`-plus-one encoding, epsilon comparisons, Brandes accumulation -- each step mirrors
/// the C source line for line, and the comments name the lines. Returns the scores plus whether
/// igraph would have warned about tiny weights (the R side raises the warning with igraph's
/// exact text, since a warning's text is part of the observable behaviour and the corpus gates
/// messages as well as values).
pub fn betweenness(net: &[f64], n: usize) -> Result<(Vec<f64>, bool), CentralityError> {
    let g = Graph::build(net, n)?;
    if g.edges.is_empty() {
        return Ok((vec![0.0; n], false));
    }
    // Reciprocals first, then checks in probed order: NaN anywhere errors before positivity,
    // positivity before the epsilon warning.
    let mut rw: Vec<f64> = Vec::with_capacity(g.edges.len());
    for e in &g.edges {
        rw.push(1.0 / e.weight);
    }
    if rw.iter().any(|w| w.is_nan()) {
        return Err(CentralityError::NanWeight);
    }
    // `igraph_vector_min`: plain linear scan, first minimum wins -- but only the *value* is used
    // here (for the two comparisons), never its position, so scan order is unobservable.
    let minweight = rw.iter().copied().fold(f64::INFINITY, f64::min);
    if minweight <= 0.0 {
        return Err(CentralityError::NonPositiveWeight);
    }
    // N.B. `f64::min` ignores NaN (returns the other operand), so the NaN check above must come
    // first -- which is also the probed order. An `is_nan(minweight)`-style check after `min`
    // would miss a NaN hiding behind a smaller value... except NaN is unordered and `f64::min`
    // returns the non-NaN one, so post-hoc detection is impossible. Check-then-min is the only
    // correct order, and it is what the probes confirm.
    let tiny = minweight <= SHORTEST_PATH_EPSILON;

    // Incidence without loops (`IGRAPH_NO_LOOPS`): loops contribute to strength but the
    // traversal never follows them.
    let out = g.out_lists(n, true);
    // Edge weights by edge ID, for the relaxation step.
    let weights = &rw;

    let mut tmpres = vec![0.0f64; n];
    // Per-source scratch, freshly zeroed each iteration -- allocation does not affect arithmetic,
    // and the source loop re-establishes the documented invariant (empty stack, zero dist/nrgeo/
    // tmpscore, empty parents) by construction rather than by reset.
    for source in 0..n {
        // `dist` holds distances PLUS ONE; `0.0` means infinity. `nrgeo` holds path counts as
        // `f64` (igraph uses doubles to avoid overflow on grids -- see the source comment).
        let mut dist = vec![0.0f64; n];
        let mut nrgeo = vec![0.0f64; n];
        let mut tmpscore = vec![0.0f64; n];
        let mut parents: Vec<Vec<usize>> = vec![Vec::new(); n];
        let mut stack: Vec<usize> = Vec::new();
        let mut heap = TwoHeap::new(n);
        heap.push(source, -1.0);
        dist[source] = 1.0;
        nrgeo[source] = 1.0;

        while !heap.empty() {
            let minnei = heap.max_index();
            // `mindist = -delete_max()`: the negation is exact.
            let mindist = -heap.delete_max();
            // `cutoff >= 0` never holds (`cutoff = -1`), so no cutoff reset.
            stack.push(minnei);
            for &e in &out[minnei] {
                let edge = &g.edges[e];
                let to = edge.to;
                // `IGRAPH_OTHER`: the far endpoint. Loops are already excluded above.
                debug_assert_ne!(to, minnei);
                let altdist = mindist + weights[e];
                let curdist = dist[to];
                // `curdist == 0` means infinity (first discovery), not a comparison.
                let cmp = if curdist == 0.0 {
                    -1
                } else {
                    cmp_epsilon(altdist, curdist, SHORTEST_PATH_EPSILON)
                };
                if curdist == 0.0 {
                    parents[to].resize(1, 0);
                    parents[to][0] = minnei;
                    nrgeo[to] = nrgeo[minnei];
                    dist[to] = altdist;
                    heap.push(to, -altdist);
                } else if cmp < 0 {
                    parents[to].resize(1, 0);
                    parents[to][0] = minnei;
                    nrgeo[to] = nrgeo[minnei];
                    dist[to] = altdist;
                    heap.modify(to, -altdist);
                } else if cmp == 0 {
                    // `cutoff < 0` is always true here, so the cutoff guard is omitted.
                    parents[to].push(minnei);
                    nrgeo[to] += nrgeo[minnei];
                }
            }
        }

        // Brandes accumulation in stack-pop (reverse discovery) order. `coeff` divides by the
        // path count *of this node*, and the operations run in exactly this sequence: the adds
        // into `tmpscore[neighbor]` happen before the node's own score is banked, and the
        // resets after.
        while let Some(actnode) = stack.pop() {
            let coeff = (1.0 + tmpscore[actnode]) / nrgeo[actnode];
            // `neis` borrows immutably while `tmpscore` is written -- disjoint indices in
            // practice (a node's parents are never itself: loops excluded), but the borrow
            // checker needs the split, so collect the parent list first. Order preserved.
            let ps: Vec<usize> = parents[actnode].clone();
            for neighbor in ps {
                tmpscore[neighbor] += nrgeo[neighbor] * coeff;
            }
            if actnode != source {
                tmpres[actnode] += tmpscore[actnode];
            }
            dist[actnode] = 0.0;
            nrgeo[actnode] = 0.0;
            tmpscore[actnode] = 0.0;
            parents[actnode].clear();
        }
    }
    // Normalisation factor is exactly 1.0 for directed + unnormalised -- an exact no-op, skipped
    // rather than applied (`x * 1.0 == x` always, including `-0.0`, so this is not a shortcut).
    Ok((tmpres, tiny))
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn cmp_epsilon_matches_documented_branches() {
        // `a == b` shortcut, infinities included.
        assert_eq!(cmp_epsilon(1.0, 1.0, 1e-10), 0);
        assert_eq!(cmp_epsilon(f64::INFINITY, f64::INFINITY, 1e-10), 0);
        // Zero branch: relative error meaningless near zero, compare against `eps * DBL_MIN`.
        assert_eq!(cmp_epsilon(0.0, 1e-320, 1e-10), 0);
        // ...but 1e-300 exceeds that bound, and 0 is simply less than it.
        assert_eq!(cmp_epsilon(0.0, 1e-300, 1e-10), -1);
        // Overflow branch: `|a| + |b|` is infinite, so the comparison is against
        // `eps*|a| + eps*|b| = 2.1e298`, and `|diff| = 1e307` exceeds it.
        assert_eq!(cmp_epsilon(1e308, 1.1e308, 1e-10), -1);
        // Relative branch.
        assert_eq!(cmp_epsilon(1.0, 1.0 + 5e-11, 1e-10), 0);
        assert_eq!(cmp_epsilon(1.0, 1.0 + 5e-10, 1e-10), -1);
        assert_eq!(cmp_epsilon(1.0 + 5e-10, 1.0, 1e-10), 1);
    }

    #[test]
    fn heap_ties_rise_and_sink_prefers_left() {
        // Push order 0, 1, 2 with equal keys: `shift_up` swaps on `>=`, so the last pushed
        // bubbles to the root -- LIFO among equals. Pops must come out 2, 1, 0.
        let mut h = TwoHeap::new(4);
        h.push(0, -1.0);
        h.push(1, -1.0);
        h.push(2, -1.0);
        assert_eq!(h.max_index(), 2);
        h.delete_max();
        assert_eq!(h.max_index(), 1);
        // Distinct keys pop largest-first regardless of insertion order.
        let mut h = TwoHeap::new(4);
        h.push(0, -3.0);
        h.push(1, -1.0);
        h.push(2, -2.0);
        assert_eq!(h.max_index(), 1);
        h.delete_max();
        assert_eq!(h.max_index(), 2);
        // `modify` to a larger key must surface the element.
        h.modify(0, -0.5);
        assert_eq!(h.max_index(), 0);
    }

    #[test]
    fn adjacency_errors_name_themselves() {
        assert_eq!(
            Graph::build(&[0.5, f64::NAN, 0.0, 0.0], 2).unwrap_err(),
            CentralityError::NonFiniteAdjacency
        );
        assert_eq!(
            Graph::build(&[0.5, f64::from_bits(0x7ff0_0000_0000_07a2), 0.0, 0.0], 2).unwrap_err(),
            CentralityError::NonFiniteAdjacency
        );
    }

    #[test]
    fn empty_graph_is_zeros_without_complaint() {
        // Probed: igraph returns zeros with no warning and no error when there are no edges.
        let (b, tiny) = betweenness(&[0.0; 9], 3).expect("empty graph");
        assert_eq!(b, vec![0.0; 3]);
        assert!(!tiny);
        let d = degrees_strength(&[0.0; 9], 3).expect("empty graph");
        assert_eq!(d.outdeg, vec![0.0; 3]);
    }
}
