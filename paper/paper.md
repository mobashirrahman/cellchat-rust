# cellchatrs: a bit-identical, drop-in accelerator for CellChat inference

*Draft — claims are limited to the evidence cited inline. Every number below names the
artifact that produced it; a claim without a pointer is a claim withdrawn.*

## Abstract

CellChat inference is dominated (92--96% of pipeline time) by its bootstrap aggregation
kernel. We reimplemented the package's numeric surface in Rust (`r-core`, zero R
dependency) behind an R package that exports matching inference entry points, and
verify it with 160 machine-checkable parity quantities: `identical()` on whole S4 objects
across a 200-configuration matrix, an installed-shim differential gate (427 comparisons),
metamorphic and fuzz suites, an independent-RNG stream (120 seeds), a full tutorial
reproduction, and a standalone CLI checked against upstream over 14 configurations (231
comparisons). On the authors' own human-skin (7,563 cells) and mouse-wound (21,557 cells)
datasets at the tutorial default `nboot = 100`, the full pipeline runs **16.7x** and
**9.7x** faster (kernel alone: 42.1x and 26.1x, all with bootstrap CIs), with 1,207,952
`Prob` values at maximum absolute difference exactly zero. We further report two
corrections the evidence forced: upstream's approximate spatial neighbour search carries a
reproducible bias (19.6% of `Prob` entries move, 1.0% of peak), and its `Prob` scores are
unbounded above (max 2025.1), so neither "probability" language nor downsampling for
runtime reasons survives measurement.

## Introduction

CellChat [@jin2021] infers cell-cell communication from single-cell transcriptomics by
permuting cell labels through a bootstrap (`computeCommunProb`, default `nboot = 100`)
and aggregating ligand-receptor communication probabilities. The method is widely used
and computationally heavy: the reference implementation needs ~90 s (7.5k cells) to
~220 s (21.5k cells) per run in R, and its own documentation recommends downsampling
large datasets. That recommendation is a runtime workaround, and runtime workarounds
deserve to be remeasured when the runtime changes.

`cellchatrs` keeps the R API and replaces the numerics: a pure-Rust crate implements the
kernel, an R package exports the matching inference entry points, and everything else --
plotting, the Shiny app, database curation -- stays upstream R code by explicit scope
decision. The contribution is therefore not a new method but a verified claim: *the same
numbers, faster, with the verification machine-readable*. The parity ledger
(`parity.json`, 160 quantities, each naming the test that pins it and the rung achieved)
is part of the deliverable, not an appendix to it.

## Results

### Bit-identity on the full surface, not just the kernel

The acceptance gate throughout is `identical()` on whole S4 objects -- values, dimnames,
row names, options, errors and messages -- never a tolerance. The 200-configuration
fixture matrix (pairwise coverage over mean type, cell count, group count, L-R count,
bootstrap depth, population size, raw use, datatype, L-R structure, Hill parameters,
trim, and pathological inputs including all-zero, single-cell, NA/NaN/Inf and duplicate
rownames) passes 1051/1051 axis pairs with scale invariance. The installed-shim gate
passes 427/427 comparisons over 9 configurations; metamorphic invariances (permutation,
relabelling, scaling, duplication, `nboot`-invariance of `Prob`) pass 40/40; 18
property tests, 4 database fuzz properties, and a 120-seed independent-RNG stream
(Kolmogorov-Smirnov $D = 0$ on 15,360 pooled values) all pass. The vignette pipeline
(`createCellChat` through `netAnalysis_computeCentrality` at `nboot = 100`) reproduces
16/16 stage comparisons `identical()`, and the standalone CLI reproduces upstream on
14 configurations (231 comparisons).

Two classes of guarantee needed more than passing tests. Floating-point contraction
(FMA) and reassociation are disabled by build flag, but a flag is a promise, so the
suite contains a control experiment (a 53-bit accumulator fails 9 tests across 6
binaries, proving the tests can see rounding), a three-way codegen comparison (default /
`+fma` / `target-cpu=native` builds agree bit-for-bit), and 275 adversarial chained
summations against a real x87 oracle pinning accumulation order. Separately, the port
found and fixed real defects the fixtures could not see: a melted-table factor level
(`interaction_name` must be a factor, caught only by whole-frame comparison), `NA`
handling in `filterCommunication` (counts, slice selection and messages), and a pathway
aggregation that emitted one row per input row instead of per group.

### Speedup, stated honestly

| workload | upstream (s) | cellchatrs (s) | speedup | Amdahl limit |
|---|---|---|---|---|
| human skin, kernel, nboot=100 | 101.4 [100.8, 105.8] | 2.41 [2.37, 2.47] | 42.1x [40.8, 44.6] | -- |
| human skin, full pipeline | 100.1 | 6.0 | 16.7x | 24.8x |
| mouse wound, kernel, nboot=100 | 245.3 [229.5, 274.9] | 9.41 [9.31, 9.63] | 26.1x [23.8, 29.5] | -- |
| mouse wound, full pipeline | 202.3 | 20.9 | 9.7x | 13.3x |
| visium spatial, kernel, nboot=100 | 17.1 [16.8, 17.2] | 0.30 [0.28, 0.31] | 57.2x [54.3, 62.0] | -- |
| human skin NL, kernel, nboot=100 | 30.3 [29.5, 31.2] | 0.48 [0.48, 0.51] | 62.8x | -- |
| human skin LS, kernel, nboot=100 | 45.0 [44.6, 47.2] | 0.66 [0.64, 0.71] | 67.9x | -- |
| embryonic E13, kernel, nboot=100 | 78.8 [78.0, 90.3] | 2.17 [2.15, 2.22] | 36.3x | -- |
| embryonic E14, kernel, nboot=100 | 83.9 [79.8, 94.8] | 2.27 [2.19, 2.66] | 36.9x | -- |

A controlled synthetic grid (2k--50k cells, fixed K/nLR/genes/nboot) holds 26--32x with
`identical()` at every point including 50k (4.7 s vs 127.3 s); both sides share log-log
slope 0.64, i.e. the same complexity profile. The 50k point passed an enforced
memory-headroom gate before running.

The remaining author datasets extend the range in both directions. Human skin NL
(2,552 cells) and LS (5,011 cells) run 62.8x and 67.9x with `identical()` nets; mouse
embryonic E13 (12,951 cells) and E14 (12,179 cells) run 36.3x and 36.9x, also
`identical()`. The small-data speedups are *larger* because upstream's per-run fixed
costs dominate at small nC while the kernel stays near-flat -- the same Amdahl shape
seen within each dataset. The NL/LS pair additionally exercises the multi-condition
workflow the data exists for: `mergeCellChat` plus comparison-mode `rankNet` agree
exactly on both sides (merge ~1x, as expected for R-side plumbing; rankNet comparison
data frames identical).

One version-skew note, because it constrains what "the authors' data" means: the
deposited NL/LS objects bundle an older CellChatDB (complex table 339x5) than the pinned
export (338x6), and exactly one pair (`IL2_IL2RA_IL2RB_IL2RG`) resolves under the former
but not the latter. Both sides are correct against their own database, so the benchmark
runs the 492-pair intersection both agree on -- the same treatment as the visium
benchmark's 134 resolvable pairs -- and records the dropped pair rather than silently
skipping it.

Medians over $\geq 5$ timed repeats (one warm-up discarded) with bootstrap 95% CIs,
taskset-pinned to 8 cores of an AMD Zen3 (32 GB); peak RSS 3.3/8.9 GB. The host is shared and
was not quiet during these runs, so each result file records the per-core busy fractions sampled
immediately before timing (`host_state.per_cpu_busy_pct`, 18--35 % on the human-skin run) and the
run was taken with the contention explicitly overridden rather than silently. The consequence is
worth stating because it bounds every absolute time below: upstream's human-skin median moved from
93.0 s to 101.4 s to 107.5 s across three runs of byte-identical code as host load varied, while
the *ratio* stayed within 1 % of 42x (42.45x, 42.14x, 44.24x). Both implementations are timed
back to back in one process on the same machine, so contention moves them together; the speedups
are the robust quantity here, the absolute seconds are the soft ones. A reproduction should expect
to match the ratio more closely than the seconds, and should use a quieter host than this one had.
The visium row
covers the spatial path (`distance.use = TRUE`, exact k-d tree); its `Prob` output is
*not* bit-identical by design (see next section), and the benchmark asserts consistency
with the measured Annoy divergence (max abs diff 0.00177, matching the recorded
0.001775) rather than identity. The pipeline column is quoted first deliberately: the kernel is 92--96% of upstream's pipeline, and
the remaining $\sim$10 s (I/O, preprocessing, presto-backed over-expression analysis at
$\sim$1x on both sides) is a floor no kernel accelerates. The Amdahl column states the
ceiling the measured shares imply. Strong scaling of the kernel is 4.9x at 8 threads
(61% efficiency; 16 threads regress on SMT).

### The spatial correction

Upstream finds spatial neighbours with an approximate index (`BiocNeighbors::queryKNN`
with `AnnoyParam`). The port uses an exact k-d tree (lowest-index tie-break, stated).
On the authors' mouse-cortex visium data (1,073 spots, 8 types, 437 pairs, `nboot=100`):
the exact tree and Annoy differ in 42 of 64 inter-group distances (2.8% max) and 1,682
of 8,576 `Prob` entries (19.6%; 1.0% of peak) -- while Annoy agrees with *itself*
exactly. The approximation is therefore a systematic bias, not noise: a reproducible
wrong answer that no user workflow would surface. The exact tree is also 274x faster on
the query itself (8.5 ms vs 2.34 s), because a 2-D exact index over a thousand points
never needed approximation.

### The downsampling correction

Upstream advises downsampling large datasets for runtime. The measured I/O floor is
seconds (2.3 s skin / 7.1 s wound for load + preprocessing) and the full 21,557-cell
pipeline runs in ~21 s, so subsetting cells to fit a runtime budget is unnecessary. The
memory headroom gate (required before quoting any $n_C \geq 50{,}000$) is implemented
but unexercised here: the largest measured dataset is 21,557 cells, and the paper
claims nothing beyond it.

### What is *not* claimed

`Prob` scores are unbounded above (observed max 2025.1: agonist/antagonist terms are
plain sums, and only the Hill term is bounded by 1), so "probability" language about
them is wrong upstream and remains wrong here -- the port reproduces the scores, not
the name. Four centrality measures (`hub`, `authority`, `eigen`, `page_rank`) call
igraph directly: they are iterative solvers whose output varies run to run, so
reimplementing them would trade identical-by-construction for version-fragile
reimplementation; degrees, strengths and betweenness are ported bit-exact. The fast
over-expression path delegates to presto (a normal approximation, not a faster
Wilcoxon): both branches select the same features on fixtures, with p-values agreeing
to $3.8\times 10^{-3}$ and no threshold crossing -- but they are different statistics,
and the paper does not call them interchangeable.

## Methods

*Port architecture.* Pure-Rust numerics (`src/rust/crates/r-core`, no R dependency, all logic in
pure functions) + extendr marshalling layer (`src/rust/crates/cellchatrs`) + R shim overriding
upstream generics + standalone CLI (`src/rust/crates/cellchatrs-cli`) reading self-describing
hex-float input. R 4.3, Rust stable/nightly per CI.

*Borough of verification.* `PLAN.md` §14 records fourteen locked scope decisions;
`docs/SEMANTICS.md` records every measured R quirk the port reproduces (MT19937 +
`sample.int`/`R_unif_index`, collapse type-7 quantiles differing from R's by up to
2 ulp in 7.4% of cases, 80-bit accumulation, NaN-as-missing, radix fast path, strict
`>` p-value counting, loop-carried `P.spatial`, `match.arg` partial matching, exact
error texts including the missing-subunit crash, FMA off, no reassociation).

## Data and code availability

Pinned upstream: jinworks/CellChat @ `75253cd` (v2.2.0.9001). Datasets: authors'
Figshare human skin, mouse wound, mouse cortex visium. All fixtures generated from
pinned upstream by scripts in `tests/parity/`; all benchmark raw timings in
`bench-runner/results/`; parity ledger in `parity.json`. License: GPL-3.

## References

[@jin2021]: Jin et al., Nat Commun 2021 (CellChat).
