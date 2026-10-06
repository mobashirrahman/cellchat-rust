# cellchat-rs — Plan for a rigorous Rust port of CellChat's inference kernel

**Status:** executed. This document is the plan as it was written and worked through; the
per-phase sections below record the state at the time each was finished, and their counts
(190 tests, 272 tests, 359/359, 124/129, `0.73 MB`) are historical. The section
immediately below is the current state, and where the two disagree, it wins.
**Target:** CellChat v2.2.0.9001 (`jinworks/CellChat`, commit
`75253cd0c9e68410e6e721a6d3a0419a1d7e358f`)
**Scope:** numeric core of the R package, exposed through a drop-in R package, with a
bit-level parity harness and a reproducible benchmark suite.

---

## Current state

### What the repository is

The R package is named **`CellChat`** and installs independently — it does not need the
original package, and does not modify its namespace. A caller loads this package in place
of the original. The plan below describes the earlier shape of the work, in which the
package was called `cellchatrs` and overrode a separately installed CellChat; the rename
is recorded in the repository history (`cea8b18`).

The Rust crates keep their names and live under `src/rust/crates/`: `r-core` (numerics,
zero R dependency), `cellchatrs` (the extendr cdylib R loads), `cellchatrs-cli` (the
standalone binary). `src/rust/r-build/` is a separate Cargo workspace root with its own
lockfile and a sha256-pinned `vendor.tar.xz`; `configure` builds it `--frozen` with
`CARGO_NET_OFFLINE=true`, so an install needs no network.

### Verification, as CI runs it

`.github/workflows/parity.yml`, green on `main`:

| job | establishes |
|---|---|
| `pinned upstream` | the checkout is the literal commit the repository claims |
| `rust stable` / `rust nightly` | 331 Rust tests; fmt; clippy `-D warnings`; the regenerated golden corpora; the codegen-variant check (default / `+fma` / `target-cpu=native` byte-identical) |
| `R 4.3` / `R 4.4` | the installed shim against the pinned upstream: `check_identical` 427/427, the 200-configuration matrix, metamorphic, CLI, centrality, delegation, the full tutorial reproduction |
| `R CMD check` | `--as-cran`, 0 errors and 0 warnings |
| `parity.json` | the machine-checkable ledger, published as an artifact and in the job summary |

`parity.json` currently reports **147 of 148 quantities at a passing rung**, plus 12
alternative-algorithm checks (160 entries). The one untested quantity is
`package.independent_dropin`, which compares two independently installed packages and so
needs the original installed beside this one; its dependency tree requires R ≥ 4.5 while
the R-side rungs require R ≤ 4.4 for bit-identity, so no single job can do both. The gate
itself runs and was verified by hand against an oracle built from the pinned SHA (25 byte
comparisons, 109 public signatures).

`R CMD check --as-cran` reports four NOTEs, all understood and none of them defects:

1. *CRAN incoming feasibility* — `New submission`. Clears on acceptance.
2. *package dependencies* — Suggests not installed in the check environment. CRAN has them.
3. *installed package size* — 7.6 MB, of which `upstream` is 4.2 MB. That tree is what all
   four `tests/test-*.R` diff against; dropping it would make `R CMD check` run no tests.
4. *foreign function calls* — `.Call(gen$address, ...)` in `.onLoad` defeats the check's
   static analysis because `gen` is a local. R sees the correct registered routine; the
   note is about check being unable to evaluate it.

### The distribution budget

The source tarball is **4,958,464 bytes** against CRAN's 5,000,000-byte incoming ceiling
(the threshold is `_R_CHECK_CRAN_INCOMING_TARBALL_THRESHOLD_`, default `5e6`).
`scripts/check_tarball_size.sh` fails the build if it crosses, printing the ten largest
contents; `R CMD check` reports the size only as a NOTE, which no gate would trip on.

What is excluded from the tarball and why, all of it in `.Rbuildignore`:

* `data/` — the five CellChatDB/PPI `.rda` files are byte-identical to the copies inside
  `inst/upstream/`, and shipping both is what put the tarball over the ceiling. The
  bundle's copy ships, so the objects remain in the installed package, and `.onLoad`
  binds the five names to it lazily and exports them, so bare `CellChatDB.human` works as
  it does under upstream's `LazyData`; only the `data(CellChatDB.human)` route is gone. `tests/test-public-surface.R` still asserts the
  datasets are byte-identical to upstream, loading them from the bundle.
* `inst/upstream/` is **not** excluded. It was measured as the alternative saving
  (4,769,468 bytes without it) and rejected: all four test targets need the reference.
* `tests/fixtures/` — cargo-test inputs read by relative path; `R CMD check` reads none of
  them, and `scripts/gen_fixtures.sh` regenerates them.

### Known gaps

* `package.independent_dropin` does not run in CI, as above.
* The four NOTEs above.
* R devel was dropped as a canary: it failed 38 of 427 comparisons on every run, all in
  `identifyOverExpressedGenes` and all with `all.equal()` TRUE, because R-devel changed
  how `data.frame()` names and how `identical()` compares `features.info`. Re-adding it
  is worth doing once that is diagnosed.

---

## 0. Executive summary

CellChat's cost is concentrated in one function, `computeCommunProb()`
(`R/modeling.R:63-327`, inside a 1297-line file). It is a pure numeric kernel with
embarrassingly parallel structure over *(ligand–receptor pair, bootstrap replicate)*.

**Phase 0 has already been executed on the target machine.** Measured with upstream
CellChat code on the authors' own example datasets at the default `nboot = 100`:

| fixture | cells | genes (data.signaling) | L-R pairs | upstream wall time | peak RSS |
|---|---:|---:|---:|---:|---:|
| human skin (Figshare 24470719) | 7 563 | 1094 | 1583 | **89.7 s** | 280 MB |
| mouse wound (Figshare 21896400) | 21 557 | 1101 | 1568 | **206.5 s** | 639 MB |

and, critically, the **cost decomposition**: the bootstrap aggregation
`aggregate(t(data), list(group), triMean)` called `nboot` times is **92–95 % of runtime**;
the `nLR × nboot` loop that *looks* like the parallel target is only ~6 %, because most
L-R pairs take the `sum(P1_Pspatial)==0` early exit on sparse data.

**Projected** Rust speedup on that same input: **~22× (7.5 k cells) to ~40×
(21.5 k cells, 206 s → ~5 s)**, rising with dataset size because R is
`Θ(nGenes·nC·nboot)` while an exact Rust implementation becomes memory-bandwidth bound.
These are hypotheses; §9 designs the experiment that falsifies or confirms them.

The plan is built in five gated phases:

| Phase | Weeks | Gate |
|---|---|---|
| **0** Baseline + fixtures | 0.5 | R reference outputs frozen, benchmark harness reproducible (CV ≤ 2 %) |
| **1** Rust core + parity | 3 | Bit-identical to R on 100 % of the Phase-0 fixture matrix |
| **2** Drop-in R shim | 1 | `computeCommunProb(object, ...)` returns an object identical to upstream |
| **3** Benchmark study | 1 | Speedup measured with CI; scaling model fitted; paper-grade figures |
| **4** Spatial + downstream | 3 | `aggregateNet`, `rankNet`, `subsetCommunication`, `filterCommunication`, `mergeCellChat`, KNN |
| **5** Packaging / CI / docs | 1 | `R CMD check` clean, reproducible CI, docs |

All scope, parity, binding, pinning, distribution and hardware decisions are **locked** —
see §14 and `UPSTREAM.md`.

Everything in §2 is *measured on this machine*, not assumed: the database dimensions,
the toolchain (extendr 0.9 builds and is callable from R), the RNG reproducibility, the
`collapse` vs R quantile discrepancy, the real fixtures, and the R cost decomposition.

---

## 1. Scope: what "fully port" should and should not mean

CellChat is ~14k lines of R:

| File | Lines | Port? | Why |
|---|---:|---|---|
| `R/modeling.R` | 1297 | **Full** | Pure numerics. This is the kernel. |
| `R/utilities.R` | 1328 | **Mostly** | `subsetData`, `extractGeneSubset`, `identifyOverExpressedGenes` (Wilcoxon path; `do.fast = TRUE` is presto and falls back), similarity/centrality — all numeric. |
| `R/analysis.R` | 3148 | **Numeric parts** | `subsetCommunication`, `rankNet`, `mergeCellChat`, comparison analysis. Heavy `data.frame`/`apply` loops → big Rust win. |
| `R/visualization.R` | 4672 | **No** | ggplot2 / circlize / ComplexHeatmap. Reimplementing this is a 10× effort with negative scientific value. |
| `R/app.R` | 1982 | **No** | Shiny UI for the CellChat Explorer. |
| `R/database.R` | 507 | **Partial** | `extractGene`, `extractGeneSubset`, `subsetDB` are small and go to Rust; `updateCellChatDB` stays in R (dplyr-heavy). |

**Decision (proposed, needs your sign-off):** the deliverable is
`cellchat-rs` = a Rust workspace + a thin R package `cellchatrs` that is a
**transparent accelerator for the existing CellChat package**. Upstream CellChat stays
installed and unmodified; `cellchatrs` overrides `computeCommunProb` (and later the
other kernels) in the search path. This preserves the entire R API, all 110 exported
functions, all tutorials, and the S4 object format, while the compute moves to Rust.

A standalone pure-Rust library + CLI (`cellchat-core`) is a **secondary deliverable**
that becomes possible once the numerics are validated; it is deliberately out of
scope for v1.

If instead you want a literal full rewrite including visualization, say so — it is
roughly 10× the effort and I would recommend against it.

---

## 2. Verified facts (measured, not assumed)

### 2.1 The algorithm, exactly as implemented

```
data.use        <- data / max(data)                     # data = as.matrix(object@data.signaling)
avg[g, j]       <- FunMean({data.use[g,c] : group[c]==j})
Lavg[i, ]       <- computeExpr_LR(geneL[i], avg, complex)     # geometric mean over subunits if complex
Ravg[i, ]       <- computeExpr_LR(geneR[i], avg, complex)
Ravg[i, ]       <- Ravg[i,] * coA(i,avg) / coI(i,avg)          # co-activator / co-inhibitor coreceptors

permutation     <- replicate(nboot, sample.int(nC, nC))        # R RNG, seed.use
avgB[g, j, b]   <- FunMean({data.use[g,c] : group[permutation[,b]][c]==j})

for i in 1..nLR:
    P1   <- outer(Lavg[i,], Ravg[i,])^n / (Kh^n + outer(...)^n)      # K x K
    if sum(P1 * P.spatial) == 0: Prob[,,i] <- P1*P.spatial; Pval[,,i] <- 1; next
    if i > nLR1: P.spatial <- P.spatial * adj.contact               # <-- loop-carried mutation
    P2   <- outer(agonist(avg, i), rep(1,K))
    P3   <- outer(antagonist(avg, i), rep(1,K))
    P4   <- outer(popFrac, popFrac) if population.size else ones
    Pnull<- P1 * P2 * P3 * P4 * P.spatial
    Prob[,,i] <- Pnull
    for b in 1..nboot:  Pboot[,,b] <- P1b * P2b * P3b * P4b * P.spatial
    pval <- rowSums(Pboot - Pnull > 0) / nboot
Pval[Prob == 0] <- 1
```

`Pboot` is materialised as a full `K x K x nboot` array per pair and then reduced.
That reduction is exact integer counting, so it parallelises without changing results.

### 2.2 Database dimensions (from `CellChatDB.human.rda`)

| Table | Dim | Notes |
|---|---|---|
| `interaction` | 3233 × 28 | 1280 Secreted, 535 Cell–Cell Contact, 424 ECM-Receptor, 994 Non-protein; 290 pathways |
| `complex` | 338 × 5 | up to 5 subunits; geometric-mean aggregation |
| `cofactor` | 32 × 16 | up to 16 subunits; agonist/antagonist/co-A/co-I |
| `geneInfo` | 26827 × 9 | official-symbol filter |

**`data.signaling` after `subsetData`: 1211 genes** (3-annotation DB) or **1446** (4-annotation).
Usable `LRsig` rows: 2239 / 3233 before over-expression filtering, typically 800–1500 after.

So the kernel's shape is: `nGenes ≈ 1.2k`, `nC = 1e3…1e5`, `K = 2…40`, `nLR ≈ 1e3`,
`nboot = 100`. The dense normalised matrix is 1211 × 21557 × 8 B = **209 MB** — fits in RAM
with room for the 24 MB bootstrap-average tensor (`1211 × 25 × 100 × 8 B = 24 MB`, which
fits in this machine's 32 MB L3 — a deliberate design target, see §6.3).

### 2.3 Environment (this machine)

| Item | Value |
|---|---|
| CPU | AMD Ryzen 7 3700X, 8c/16t, Zen 2, **AVX2 (no AVX-512)**, 32 MB L3 |
| RAM | 32 GB (28 GB free) |
| R | 4.3.3, `RNGkind = Mersenne-Twister`, `sample.kind = Rejection` |
| Rust | 1.98.1, cargo 1.98.1, gcc 13.3 |
| libR | `/usr/lib/R/lib/libR.so`, `R CMD config` functional |

### 2.4 Toolchain validation already performed

* `extendr-api 0.9.0` + `extendr-macros 0.9.0` compile as a `cdylib` and are callable
  from R. **The full path R package → `configure` → cargo → `dyn.load` → R function has
  been built and executed end-to-end** (see `R/`, `configure`, `src/rust/crates/cellchatrs` in
  this repo), including the call `get_num_threads()` returning the real rayon pool size.
* **Five gotchas found and fixed empirically** (each would otherwise cost days):
  1. The umbrella crate `extendr` is **not resolvable from the crates.io sparse index**
     here (HTTP 404). Depend on `extendr-api` + `extendr-macros` directly.
  2. `extendr_module!` syntax changed in 0.9: it is now
     `extendr_module! { mod name; fn foo; }`, **not** `extendr_module!(mod name);`.
  3. extendr emits `R_init_<mod>_extendr`, but `dyn.load()` auto-invokes
     `R_init_<dll-filename>`. A shim with the exact module name is mandatory. Empirically
     R derives the name from the *file name* — `libr_bindings.so` was called as
     `R_init_libr_bindings`, not `R_init_r_bindings`.
  4. **The crate must be named to match the R package** (`cellchatrs`, not `r-bindings`).
     extendr resolves the target DLL by *module* name during
     `R_registerRoutines`; when they disagree, registration silently registers nothing and
     R merely reports "could not find function". This is a very confusing failure mode.
  5. Generated wrappers must be emitted with `use_symbols = FALSE` so they call
     `.Call("wrap__foo", PACKAGE = "cellchatrs")`. With `TRUE` they reference a
     `NativeSymbolInfo` object that is not in scope in a namespace.
* `configure` script + `src/Makevars` (no rules) so that **`R CMD INSTALL .` is the only
  command a user needs** and `R CMD check` exercises the real Rust build.
* **Benchmark gotcha:** `parallel::detectCores(logical = FALSE)` returns **16** on this
  8c/16t Zen 2 host. Thread counts must always be pinned explicitly
  (`CELLCHATRS_THREADS`, or `taskset`), never inferred.
* `rayon::ThreadPoolBuilder::build_global()` succeeds only once per process, so
  `set_num_threads()` after the first kernel call is a silent no-op. Fixed by
  `request_num_threads()` (records the request) + lazy pool construction in
  `ensure_pool()`; the buggy first version was caught by this smoke test.
* **Quantile parity is not free.** `triMean` is
  `mean(collapse::fquantile(x, probs=c(.25,.5,.5,.75), na.rm=TRUE))`, and `collapse`'s
  type-7 implementation is *not* bit-identical to R's `stats::quantile`: over 50 000
  random vectors. Measured on 20 000 of them, **1.13 % of the `fquantile` outputs and
  3.51 % of the final `triMean` values differ from R, max 2.22e-16 (1 ulp)**. We must
  port **collapse's** arithmetic, not R's, to be bit-exact.
  `fquantile` also honours NaN as missing (verified: `fquantile(c(1,2,NaN,3), .5) == 2`),
  and switches to a radix-order fast path when `length(x) > 1e5`. All three behaviours must
  be reproduced.
* R's `sample.int` is reproducible under `set.seed`, as required.
* `install.packages()` into a shared library must be **sequential** — parallel installs
  collide with `failed to lock directory` and leave `00LOCK-*` behind.
* CellChat's `Imports:` list (`NMF`, `BiocNeighbors`, `ComplexHeatmap`, `reticulate`,
  `plotly`, `shiny`, `ggpubr`, `sna`, …) is heavy and several are GitHub-only. **The
  oracle harness does not need it**: sourcing `R/modeling.R` into an environment plus a
  minimal S4 class that declares only the slots `computeCommunProb` touches is enough,
  and that approach is already working (it produced every number in §2.5). This is the
  concrete mitigation for Risk R2.

### 2.5 Measured R baseline (this machine, `taskset -c 0-7`)

Both fixtures are the authors' own example data from Figshare (§2.6), fed through the
real upstream `computeCommunProb()` with `triMean`, `population.size = FALSE`,
`Kh = 0.5`, `n = 1`, `seed.use = 1`, `datatype = "RNA"`.

**Fixture A — human skin (`data_humanSkin_CellChat.rda`, Figshare 24470719)**
17 328 genes × **7 563 cells**, 12 cell types, 4.5 % nnz; `data.signaling` = **1094 genes**,
`LRsig` = **1583** pairs. Log1p-CPM normalised.

| config | wall time (3 repeats) | median |
|---|---:|---:|
| `nLR=1`, `nboot=1` | 2.30 | 2.30 s |
| `nLR=1`, `nboot=10` | 10.79 | 10.79 s |
| `nLR=1`, `nboot=100` | 87.69 / 84.18 / 83.97 | **84.18 s** |
| `nLR=1583`, `nboot=10` | 11.40 | 11.40 s |
| `nLR=1583`, `nboot=100` | 89.10 / 89.66 / 91.76 | **89.66 s** |

Peak RSS 280 MB. Result is non-degenerate: `Prob` 12×12×1583, 1035 non-zeros (0.45 %),
`Pval` min 0, 0.43 % below 0.05.

**Cost decomposition (the single most important result of Phase 0):**

| stage | time | share |
|---|---:|---:|
| setup (`as.matrix`, `max`, first `aggregate`) | ~1.3 s | 1.5 % |
| **bootstrap aggregation** (`aggregate` × `nboot`) | **~83 s** | **92 %** |
| L-R × bootstrap loop (1583 × 100 = 158 300 iterations) | ~5.5 s | 6 % |
| **total** | 89.7 s | 100 % |

Cross-check: an isolated `aggregate(t(X), list(g), triMean)` costs **1.83 s** at
1782 genes × 5000 cells and ~1.27 s at 1094 × 7563, i.e. ~1.27 s per bootstrap →
100 boots ≈ 127 s, consistent once setup is removed.

**Fixture B — mouse wound (`data_wound_CellChat.rda`, Figshare 21896400)**
17 090 genes × **21 557 cells**, 25 cell types, 6.7 % nnz; `CellChatDB.mouse`;
`data.signaling` = **1101 genes**, `LRsig` = **1568** pairs.

| config | wall time |
|---|---:|
| `nboot=1` | 4.41 s |
| `nboot=2` | 7.51 s |
| `nboot=10` | 23.54 s |
| **`nboot=100` (CellChat default)** | **206.54 s** (3.4 min), peak RSS 639 MB |

Non-degenerate: 18 120 non-zero `Prob` cells (1.85 %). Least-squares fit over the four
points gives **2.13 s per bootstrap** for 1101 genes × 21 557 cells, i.e. the aggregation
is again ≈ 95 % of the default-configuration runtime.

Two consequences that reshape the engineering priorities:

1. **The obvious optimisation target is nearly irrelevant.** The `nLR × nboot` loop
   — the part that "looks" parallelisable — is ~6 % of runtime, because on sparse data
   most L-R pairs take the `sum(P1_Pspatial)==0` early exit (only ~7 of 1583 pairs produce
   a non-zero network in fixture A). Effort must go into the **aggregation**.
2. **The speedup headroom is essentially the whole runtime.** If aggregation goes to
   ~4 s and the loop to ~0.5 s, the kernel drops from 206.5 s to ~5 s on fixture B.
   Because R's cost is `Θ(nGenes·nC·nboot)` while the Rust cost becomes
   memory-bandwidth-bound at `Θ(nGenes·nC·nboot / cores)`, **the ratio increases with
   `nC`** — exactly the regime where users report 2–3 h.

Measured variance across repeats is 1.5–2 % once the page cache is warm, but an early
measurement of the same config gave 128 s instead of 84 s when another job shared the
machine — **background-load control is not optional** (see §10.2).

### 2.6 Real fixtures are obtainable

The Figshare project (157272) is reachable; verified downloads:

| Dataset | Size | Use |
|---|---:|---|
| `data_wound_CellChat.rda` | 57 MB | 17090 genes × **21 557 cells**, 25 labels, 6.7 % nnz — primary scRNA fixture |
| `data_humanSkin_CellChat.rda` | 13 MB | 7563 cells, 12 labels, 2 conditions — multi-condition |
| `visium_annotated.RData` | 4.7 MB | spatial fixture (mouse cortex visium) |
| `cellchat_visium_mouse_cortex.rds` | 2.2 MB | **already a CellChat object** — can be fed straight to `computeCommunProb` |
| `Spatial_A1_adult_with_predictions.RDS` | 30 MB | spatial multiomic |
| `cellchat_embryonic_E13/E14.rds` | 95/101 MB | comparison/multi-dataset |

GEO FTP and 10X downloads also reachable for scaling to 50k–100k cells.

---

## 3. The performance model (hypothesis to be tested in Phase 3)

The §2.5 decomposition is the load-bearing result: **the bootstrap aggregation is 92 % of
runtime**. Every design decision below follows from that.

| Stage | Complexity | R cost mechanism | R cost (measured, humanSkin, `nboot=100`) |
|---|---|---|---|
| `aggregate(t(data.use), list(group), triMean)` × `nboot` | Θ(nGenes·nC·nboot) | `t()` copy, matrix→data.frame coercion, `split()`, then `nGenes·K` `collapse::fquantile` calls per bootstrap | **~83 s (92 %)** |
| L-R × bootstrap loop | Θ(nLR_active·nboot·K²) | `data.frame` `$` extraction, `cofactor_input[name, grepl(...)]`, `unlist(df[x,])`, `Matrix::crossprod` on K×K — *all inside the inner loop* | ~5.5 s (6 %) |
| setup | Θ(nGenes·nC) | `as.matrix(dgCMatrix)` (0.11–0.17 s), `max`, one `aggregate` | ~1.3 s (1.5 %) |

Isolated micro-measurements (1782 genes × 5000 cells × 25 groups):

| operation | time |
|---|---:|
| `aggregate(t(X), list(g), triMean)` | **2.149 s** |
| `split` + `vapply(apply(...))` hand-written | **0.679 s** (3.2× faster) |
| `t(X)` alone | 0.309 s |
| `as.matrix(dgCMatrix)` | 0.109 s |

So even a *faithful but naive* Rust aggregation already targets a 2–3× win over R's own
`aggregate`, and the real win comes from eliminating interpreter/`data.frame` overhead and
parallelising.

**Projected Rust cost for the same job**

* Aggregation becomes memory-bandwidth bound. `data.use` is
  `1094 × 7563 × 8 B = 66 MB`; streaming it once per bootstrap = 6.6 GB of reads over
  100 boots → **~1–3 s on 8 cores**, versus 83 s.
  (At `nC = 21 557` this becomes 209 MB × 100 = 21 GB → **~3–6 s**.)
* L-R × bootstrap loop: `nLR_active · nboot · K² ≈ 1.4 × 10⁸` element-ops → **<1 s**
  after the `data.frame` lookups are hoisted (§6.2).

**Projection for fixture B (the realistic default workload): 206.5 s → ~5 s, i.e. ~40×.**

| fixture | `nC` | `nGenes` | R `nboot=100` | Rust (projected) | speedup |
|---|---:|---:|---:|---:|---:|
| humanSkin (measured) | 7 563 | 1094 | 89.7 s | ~4 s | ~22× |
| mouse wound (measured) | 21 557 | 1101 | 206.5 s | ~5 s | ~40× |
| 50 k cells | 50 000 | ~1 200 | ~480 s (extrapolated) | ~12 s | ~40× |
| 100 k cells | 100 000 | ~1 200 | ~960 s (extrapolated) | ~24 s | ~40× |

Rust cost model: streaming `nGenes × nC × 8 B` once per bootstrap
(fixture B: 190 MB × 100 = 19 GB) divided across 8 cores, plus the L-R loop at
`nLR_active · nboot · K² ≈ 1.6 × 10⁸` element-ops → `<1 s`.
`1101 × 7563 × 8 B = 66 MB` for fixture A.

These are **hypotheses, not claims**. The experiment design in §10 is built to falsify
them; the published number is whatever Phase 3 measures. The `nC ≥ 50 000` rows are
extrapolations from a single measured point each and must be re-measured directly
(R needs ~16 min at 100 k cells — affordable exactly once).

Two things the model deliberately does **not** yet include, both of which will be measured
in Phase 3 and reported in the paper's Amdahl section:

* `identifyOverExpressedGenes` and `subsetData` (upstream cost, unaffected by this port);
* the fixed ~1.3 s setup, which does not parallelise.

---

## 4. Semantics contract (the spec we implement against)

Bit-exactness is a **ladder**, decided per quantity, and the achieved rung is recorded in
a machine-readable `parity.json` that CI publishes.

| # | Quantity | Target rung | Notes |
|---|---|---|---|
| R1 | `sample.int(n, n)` permutations | **bit-exact** | MT19937 + `R_unif_index`; trivial to port, high value |
| R2 | `triMean` (collapse type-7) | **bit-exact** | port collapse's exact arithmetic incl. radix fast path |
| R3 | `truncatedMean` (R type-1 cutoffs) | **bit-exact** | |
| R4 | `median`, `geometricMean`, `thresholdedMean` | **bit-exact** | `log/exp` must match libm; test on libm-variance |
| R5 | Hill function, outer products, `P.spatial` | **bit-exact** | pure IEEE ops, same order |
| R6 | `Prob` array | **bit-exact** | follows from R1–R5 |
| R7 | `Pval` array | **bit-exact (exact integers)** | counting, order-independent |
| R8 | `dimnames`, `options$parameter`, slot contents | **bit-exact** | `identical()` on the whole S4 object |
| R9 | `collapse::fquantile` NaN/NA handling | **bit-exact** | NaN treated as missing when `na.rm=TRUE` |
| R10 | `mean()` accumulation over 4 elements | **Exact** | R accumulates in `LONG_DOUBLE`; Rust needs a real 64-bit-mantissa accumulator (`u128` or `long double` FFI) — an `f64` sum is *not* sufficient |
| S1 | `BiocNeighbors::queryKNN` (Annoy, approximate) | **documented divergence, intentional** | replaced by an exact k-d tree (§14.4); divergence measured on real spatial data in Phase 4 |

**R-isms that must be reproduced verbatim** (each gets a dedicated test):

1. `P.spatial <- P.spatial * adj.contact` is a **loop-carried mutation** at
   `R/modeling.R:232`. It fires on the first iteration where `i > nLR1`. Although
   idempotent in practice (`adj.contact ∈ {0,1}`), it must be applied at exactly the
   same iteration because the `sum(P1_Pspatial)==0` early-exit for iteration `nLR1+1`
   tests the *pre*-multiplied matrix.
2. `cofactor_input[missing_name, ]` yields an all-`NA` row; `unlist()` then injects `NA`,
   and `v[NA]` (logical indexing by `NA`) inserts an `NA` element. `intersect()` later
   drops it. Net effect is benign but the code path must be replicated.
3. `intersect(coreceptor.subunits, rownames(data.use))` may yield length 0 → returns
   `matrix(1)`, i.e. **no modulation**, not zero.
4. `geometricMean` uses `exp(mean(log(x), na.rm=TRUE))`; a single zero makes the whole
   group mean `0` (`log 0 = -Inf`).
5. `data.use <- data/max(data)` — `max` over the **entire matrix including zeros**, and
   emits a warning (`-Inf`) when all entries are non-finite.
6. `Pval[Prob == 0] <- 1` is applied *after* the loop, overriding any per-pair value.
7. `Pboot - Pnull > 0` is a **strict** comparison; exact ties count as non-rejections.
8. `aggregate()` returns group levels in factor order; `levels(group)` order defines the
   output `dimnames`. If `nlevels != length(unique(group))`, R `stop()`s — reproduce the error.
9. `type` argument is matched via `match.arg`, so partial names (`"tri"`, `"med"`) work.
10. When `population.size = TRUE`, `table(groupboot)/nC` is recomputed per (i, b) in R;
    we hoist it — but it must give identical values.
11. **`computeExpr_complex` has no bounds check.** `R/modeling.R:537` does
    `data.use[RsubunitsV, ]` where `RsubunitsV` are raw subunit *names*, unlike
    `computeExpr_coreceptor` which does `intersect(..., rownames(data.use))` first. If a
    multi-subunit complex has even one subunit absent from `data.signaling`, upstream
    aborts with `subscript out of bounds`. **Confirmed empirically** while building the
    humanSkin fixture (8 of 1614 L-R pairs were dropped to make it run).
    **Decision (locked):** replicate the error with the identical message, and
    separately file it upstream as a bug — see §14.3.

---

## 5. Architecture

```
┌─────────────────────────── R (unchanged CellChat object model) ──────────────────────────┐
│  object@data.signaling (dgCMatrix)   object@LR$LRsig   object@DB$complex / $cofactor      │
│  object@idents (factor)               object@options$datatype                            │
└───────────────┬──────────────────────────────────────────────────────────────────────────┘
                │  extendr-api, zero-copy borrow of the SEXP (dense + sparse paths)
┌───────────────▼──────────────────────────────────────────────────────────────────────────┐
│  src/rust/crates/r-bindings   argument marshalling only: no numerics, no logic                    │
└───────────────┬──────────────────────────────────────────────────────────────────────────┘
┌───────────────▼──────────────────────────────────────────────────────────────────────────┐
│  src/rust/crates/r-core                                                                             │
│    ├── rng/rng.rs          MT19937 + R_unif_index + R unif_rand                           │
│    ├── stats/quantile.rs   collapse type-7, R type-1, median, trimean, trimmed mean       │
│    ├── db.rs               complex/cofactor subunit tables → flat gene-id arrays          │
│    ├── expr.rs             computeExpr_LR / _complex / _coreceptor / _agonist / _antagonist│
│    ├── aggregate.rs        the bootstrap group-mean tensor (parallel, deterministic)      │
│    └── prob.rs             Hill, outer products, P.spatial, LrPlan, the whole i/b loop      │
│  rayon (work-stealing); determinism by construction: every output cell is written by      │
│  exactly one task, no cross-task reductions.                                              │
└─────────────────────────────────────────────────────────────────────────────────────────┘
```

**Data layout decisions**

* `data.use`: column-major `genes × cells`, `f64`, built once. Sparse SEXP is densified
  in Rust (avoids R's 209 MB `as.matrix` copy — measured 0.109 s, small but free to save).
* `group_idx: Vec<u16>` (cell → original cluster), `K` ≤ 65 535.
* Per bootstrap `b`, precompute `boot_cells: Vec<Vec<u32>>` = cell indices of each
  bootstrap cluster, built as a concatenation of *contiguous slices* of the original
  clusters. This is possible because `group[permutation[,b]]` is a relabelling: bootstrap
  cluster *j* receives cells from each original cluster *k* in one contiguous run of the
  original index list. Consequence: the gather pattern is reused by all `nGenes`, which
  keeps the inner loop in L1/L2.
* `avgB`: `Vec<f64>` of `nGenes · K · nboot`, laid out `[gene][cluster][boot]` so the
  per-pair inner loop reads a 20 KB contiguous slab (L1-resident).

**Parallelisation plan**

| Stage | Parallel over | Why safe |
|---|---|---|
| aggregate | `(gene, boot)` tasks | disjoint writes |
| L-R loop | `i` (ligand–receptor pair) | disjoint writes to `Prob[,,i]`, `Pval[,,i]` |
| p-value reduce | `i` | integer counts, order-independent |

Rayon thread count is configurable (`CELLCHATRS_THREADS`); if the caller runs inside a
`future` multisession with *W* workers, we default to `max(1, cores/W)` to avoid
oversubscription, and detect it via `future::nbrOfWorkers()`.

---

## 6. The algorithm design in Rust

### 6.1 Aggregation — **the hot path (92 % of runtime)**

> Engineering priority, re-derived from §2.5: the aggregation is where the time is.
> The L-R loop (§6.2) is a rounding error on sparse data and a secondary target on dense
> data. Any design that optimises the outer loop and hand-rolls the inner one is backwards.

Per `(gene g, boot b)`: stream the `nC` values of gene `g`, scatter each into one of `K`
per-cluster buckets (using the precomputed `boot_cells[b][j]` lists), then take three
order statistics per bucket.

Optimisation ladder. **Each rung is gated by a criterion benchmark, never assumed:**

* **L0 (baseline, correctness)** — gather into `Vec<f64>` buckets, `select_nth_unstable` ×3,
  type-7 interpolation. Expected ~10–20× over R's `aggregate` on its own (the measured
  hand-written R equivalent is already 3.2× faster).
* **L1 zero short-circuit** — during the gather, count exact zeros per bucket; if
  `n_zero ≥ 0.75·n_bucket`, then Q25, Q50, Q75 are all 0 and `triMean = 0` with **no
  selection at all**. Costs ~1 counter per element.
  This is the highest-leverage rung and it is *fixture-driven, not hand-waved*: the real
  humanSkin fixture is 3.9 % dense, so the median gene-cluster bucket is far past 75 %
  zeros. **Measure the short-circuit hit rate per gene in Phase 1 and publish the
  distribution** — it is a result in its own right.
* **L2 bitmask prefilter** — for each `(gene, original cluster)` precompute a bitset of
  "is zero". A bootstrap bucket's zero count is then `popcount(mask_b[j] & zeros[g][k])`,
  i.e. O(n_bucket/64) instead of O(n_bucket). Converts L1 into a ~64× cheaper test for the
  majority of genes and lets the whole bucket be skipped without touching its values.
* **L3 radix select** — for buckets that survive L1/L2, one 256-bin MSD pass brackets all
  three target ranks at once, then refine only the bracketed bins. O(n) with a small
  constant and a single streaming pass, replacing three quickselects.
* **L4 shared-index exploitation** — the *scatter pattern* is gene-independent. Hoist the
  index lists so the three targets are located in a single fused traversal, and tile
  genes × cells so a tile fits in L2.
* **L5 — REMOVED.** An `f32` fast path for `avgB` (halves memory traffic; plausibly
  doubles the speedup at large `nC`) was considered and **rejected**: under decision 2
  (bit-identical) it forfeits rung R6, so it is deleted rather than shipped behind a flag.
  For the same reason FMA contraction is disabled in the aggregation and `avgB` stays `f64`.

Expected end state: aggregation becomes memory-bandwidth bound
(`Θ(nGenes·nC·nboot / cores)`), which is the theoretical floor for an exact
implementation.

### 6.2 L-R loop — correctness-critical, cost-secondary

For each pair `i` and bootstrap `b`:
`P1b[s,t] = h(Lb[s]·Rb[t])`, `P2b[s]`, `P3b[t]`, `P4b[s]P4b[t]`, `P.spatial[s,t]`.
Compute `W_b = (P1b ∘ Psp) ⊗ (P2b ⊗ P3b ∘ P4b)` and compare to `Pnull`. Since `K ≤ 40`,
this is a 40×40 f64 tile that stays in L1; `nLR·nboot·K² ≈ 1.4 × 10⁸` element-ops total.

Even though it is only ~6 % of runtime on sparse fixtures, it must be **correct and
parallel** because (a) on dense fixtures a much larger fraction of L-R pairs is active
(the `sum(P1_Pspatial)==0` early exit does not fire), and (b) it carries the majority of
the R-isms in §4, so it is where parity bugs will live. Optimisations:
The dominant remaining cost is **not** arithmetic but repeatedly recomputing
`computeExpr_LR`/`_coreceptor` for the same genes across pairs — so:
* precompute a **gene → group-mean vector** lookup keyed on `(gene_id, boot)` and memoise;
* precompute **per-pair LR descriptor** (`geneL id`, `geneR id`, subunit id-lists,
  cofactor id-lists, agonist/antagonist ids) once, outside all loops, as flat `u32` arrays
  with no `Option<String>` lookups.

This removes 100 % of the `data.frame` `$`/`unlist`/`dplyr::select` overhead that makes the
R inner loop expensive.

### 6.3 Why `avgB` fits in cache

`nGenes(1446) × K(30) × nboot(100) × 8 B = 34.7 MB`. At the top of the range it just
exceeds L3, so the L-R loop is streamed rather than resident. Report both regimes;
for the common case (`nGenes ≈ 1211`, `K ≤ 25`) it is L3-resident, which is why the
design is `[gene][cluster][boot]` rather than `[pair][cluster][boot]`.

---

## 7. Repository layout

```
cellchat-rust/
├── Cargo.toml                    # workspace  (r-core + cellchatrs)
├── configure                     # builds the cdylib and stages src/cellchatrs.so
├── DESCRIPTION / NAMESPACE       # the R package: Depends: CellChat, useDynLib(.registration)
├── R/
│   ├── zzz.R                     # .onLoad: register extendr wrappers, size the rayon pool
│   ├── modeling.R                # computeCommunProb() override (Phase 2)
│   └── dataset.R                 # CellChatDB -> Arrow (Phase 5)
├── src/rust/crates/
│   ├── r-core/                   # pure Rust numerics, zero R dependency
│   │   ├── src/{lib,rng,stats,db,expr,aggregate,prob,kernel,spatial}.rs
│   │   ├── benches/aggregate.rs  # criterion micro-benchmarks
│   │   └── tests/                # proptest + golden files
│   └── cellchatrs/               # extendr cdylib; marshalling only, no numerics
├── tests/
│   ├── fixtures/                 # frozen RDS reference outputs + manifest
│   ├── parity/                   # differential tests vs R
│   ├── fuzz/                     # randomised config generator
│   └── bench/                    # memory profilers (heaptrack / valgrind massif)
├── bench-runner/                 # prof_r.R / prof_real.R — the Phase-0 harness
└── docs/
    ├── SEMANTICS.md              # §4 of this plan, expanded into a contract
    └── BENCHMARKS.md             # results, regenerated per release
```

**Status of this scaffold:** the workspace builds (`cargo build --release`),
`R CMD INSTALL .` succeeds, and the R→Rust→R call path is proven
(`cellchatrs_threads()` returns the live rayon pool size and responds to
`set_num_threads()` / `CELLCHATRS_THREADS`). The numerics in `r-core/src/*.rs` are
`prob.rs` now owns what used to be the `kernel.rs` placeholder: there is one kernel, and
it is the one the corpus tests.

---

## 8. Phase-by-phase work plan

### Phase 0 — Baseline & fixtures (3–4 days)

**Goal:** a reference oracle and a trustworthy stopwatch, before any Rust is written.

1. **R baseline install.** Create an isolated library (`.rlib`), install
   `collapse, future, future.apply, pbapply, dplyr, igraph, ggplot2, Rcpp, RcppEigen,
   scales, reshape2, circlize, cowplot, RSpectra, irlba, stringr, magrittr`,
   then `devtools::install_github("jinworks/CellChat")`.
   *Gotcha already hit:* parallel `install.packages` into one library causes
   `failed to lock directory` — install sequentially or use one lib per worker.
   *Known risk:* `NMF (>= 0.23.0)` is in `Imports` and can be painful to build.
   Mitigation: use `pkgload::load_all("CellChat")`, which defers the heavy optional
   imports; the inference path only needs `collapse` + `Matrix`.
2. **Fixture matrix** (see §9.1) → 200–400 configurations.
3. **Oracle generation.** For each config, run upstream `computeCommunProb` and serialise
   `Prob`, `Pval`, `dimnames`, `options$parameter`, `options$run.time`, plus
   `computeCommunProbPathway`, `aggregateNet`, `rankNet`, `net$weight`, `net$count`
   and `subsetCommunication()` frames. Store as `.rds`; index with a manifest
   (config hash → artefact path).
4. **Benchmark harness.** R side uses `bench::mark(min_iterations=…, memory=TRUE)`;
   system metrics via `/usr/bin/time -v` and `getrusage()`. Pin CPUs with `taskset`.
   Warm-up run discarded. ≥5 repetitions, report median + bootstrap 95 % CI.
5. **Cost attribution.** Wrap each stage with timers to publish an R cost breakdown
   (this is the evidence for the paper).

**Gate 0:** benchmark harness CV ≤ 2 % across repeats; oracle manifest complete.

### Phase 1 — Rust core + unit parity (2.5–3 weeks) — **complete**

**Completed and bit-exact** (190 Rust tests, all green across 15 binaries, plus the
installed-shim differential gate at 145/145 `identical()` and `R CMD check`
**Status: OK** on the built 0.73 MB tarball -- no errors, warnings or notes -- installing
and loading with `configure` compiling the Rust cdylib, and `R CMD check` running the R
tests):

| module | what | evidence |
|---|---|---|
| `rng.rs` | MT19937 + `R_unif_index` + `sample.int` | 879-record corpus from R, >2 M compared values |
| `longdouble.rs` | x87 80-bit add/sub/mul/prod for R's `mean`/`prod` | 9 301 cases vs real `long double` in C |
| `stats.rs` | `fquantile` type 7, `triMean`, `geometricMean`, `thresholdedMean`, `median`, `trimmedMean`, `r_prod` | 1 395 vectors, bit equality vs R + collapse |
| `db.rs` | `CellChatDB` load, `extractGene`/`extractGeneSubset` | both species, manifest-pinned by MD5; 1446 genes at 4 annotations, 1211 at 3, order for order |
| `expr.rs` | `computeExpr_LR`/`_complex`/`_coreceptor`/`_agonist`/`_antagonist` and the two exported `computeExprGroup_*` | 416 × 6 fixture, 8 + 3 + 30 + 1 + 1 + 10 records, all bit-exact including the two `subscript out of bounds` aborts |
| `aggregate.rs` | `aggregate(matrix, list(factor), FUN)` + `GroupMean` dispatch | the primitive `computeAveExpr` and `computeExprGroup_*` share |
| `prob.rs` | `computeCommunProb`: the `i` x `b` bootstrap permutation test, end to end | **8 parameter configurations, `Prob` and `Pval` bit-identical to upstream**, including the `sum(P1_Pspatial)==0` early exit, the `nLR1` `P.spatial` mutation, the `match.arg` type dispatch and `population.size` |
| `net.rs` | `aggregateNet` (unfiltered), `subsetCommunication` (`slot.name="net"`, `mode="single"`) | 4 fixtures x 3 thresholds x 3 filters, row order, column set, factor levels, every cell |
| `de.rs` | `computeAveExpr`, `subsetDB`, `subsetData`'s gene list and annotation reordering | 3 `type` values, R's `intersect` order, the `non_protein` flip, both error messages |
| `mathfn.rs` | Cody's `pnint` (`pnorm`, both tails), `choose`, `lgammafn`/`lbeta`/`lgammacor` | `pnorm` bit-exact on 30 lower + 30 upper R values; `choose` on 12; `lgamma` matches to ~1e-15 |
| `net.rs` (filtered branch) | `aggregateNet` with `sources.use`/`targets.use`/`signaling`/`pairLR.use`: the byte-wise `group_by` key order, and upstream's `str_split(key, "|")` regex bug, which makes the result `0 x 0` or `k x k` of zeros for every filter | 12 fixtures, whole-object `identical()`; reproduced rather than "fixed" -- see `docs/SEMANTICS.md` |
| `pathway.rs` | `computeCommunProbPathway`: the LONG_DOUBLE sums, the two different `apply` summation orders, the `aperm` axis layout and the decreasing-total reordering | 9 fixtures, `netP$prob` bit-for-bit; `bigmag` pins the 80-bit accumulation (`30000000000000004` vs `30000000000000000`), `tied`/`tied3` pin `sort`'s tie order, `k1` pins the preserved upstream `aperm` failure |
| `wilcox.rs` | `stats::wilcox.test` (statistic, p-value, `correct`, exact/normal branch), `rank`, `p.adjust("bonferroni")`, `mean.fxn`, `round(x, 3)`, `identifyOverExpressedGenes(do.fast = FALSE)`, and the `do.DE = FALSE` branch | 13 x 3 x 2 = 78 `wilcox.test` records bit-exact; `dwilcox` 231 cases and `pwilcox` 140 cases bit-exact; the marker table row for row including the collapsed `0x1` schema and the feature row names |

`identifyOverExpressedGenes` is the other hot spot the Amdahl model depends on
(risk **R4**), so it is ported rather than deferred. Nine R-isms in it are documented in
`docs/SEMANTICS.md`; three of them would each have made every all-ties gene look
maximally significant, dropped every feature at a non-zero `thresh.pc`, or reordered the
LR table.

The `i` loop is now fully name-free: every ligand, receptor, agonist, antagonist and
cofactor is resolved once into a `LrPlan` of `usize` gene indices, and the per-`(pair,
replicate)` work runs entirely in those indices. That is the whole performance argument —
upstream re-does the equivalent `data.frame` subsetting `nLR * (nboot + 1)` times.
Parallelism over the `i` loop and SIMD over the `n_groups^2` inner product are the next
step, and neither changes any bit: the RNG is consumed entirely before the loop starts, so
the loop is embarrassingly parallel.

### Phase 2a complete — the R shim returns `identical()` objects

`R/modeling.R` overrides `computeCommunProb` and delegates to `r-core`. It marshals and
nothing else; the pinned upstream body is kept verbatim in the same package as
`cellchatrs_upstream_computeCommunProb`, used both as the differential reference and as a
`CELLCHATRS_FALLBACK=1` escape hatch.

`tests/parity/check_identical.R` runs the installed shim and the upstream body on the same
objects and requires `identical()`:

```
tri_ps0    prob=TRUE  pval=TRUE  parameter=TRUE  dimnames=TRUE
tri_ps1    prob=TRUE  pval=TRUE  parameter=TRUE  dimnames=TRUE
trim_ps0   prob=TRUE  pval=TRUE  parameter=TRUE  dimnames=TRUE
thresh_ps1 prob=TRUE  pval=TRUE  parameter=TRUE  dimnames=TRUE
median_ps0 prob=TRUE  pval=TRUE  parameter=TRUE  dimnames=TRUE
hill2      prob=TRUE  pval=TRUE  parameter=TRUE  dimnames=TRUE
kh1e3      prob=TRUE  pval=TRUE  parameter=TRUE  dimnames=TRUE
matcharg   prob=TRUE  pval=TRUE  parameter=TRUE  dimnames=TRUE
ambiguous  error=TRUE   'arg' should be one of "triMean", "truncatedMean",
                        "thresholdedMean", "median"

IDENTICAL: 0 failing comparisons out of 82, over 9 configurations
```

Extended to the data-entry points:

```
ave_tri        dim=TRUE  rownames=TRUE  colnames=TRUE  values=TRUE
ave_trim01     dim=TRUE  rownames=TRUE  colnames=TRUE  values=TRUE
ave_trim05     dim=TRUE  rownames=TRUE  colnames=TRUE  values=TRUE
ave_median     dim=TRUE  rownames=TRUE  colnames=TRUE  values=TRUE
ave_featsel    dim=TRUE  rownames=TRUE  colnames=TRUE  values=TRUE
ave_matcharg   dim=TRUE  rownames=TRUE  colnames=TRUE  values=TRUE
ave_argerr      equal=TRUE   'arg' should be one of “triMean”, “truncatedMean”, “median”
subdb_default  nrow=2238  names=TRUE
subdb_nonprotein nrow=3232  names=TRUE
subdb_contact_only nrow=535   names=TRUE
subdb_empty    nrow=2238  names=TRUE
subdb_keyerr    equal=TRUE
```

Extended to the downstream functions, fed a `net` the *upstream* kernel produced:

```
agg_thresh005  count=TRUE  weight=TRUE  dimnames=TRUE
agg_thresh01   count=TRUE  weight=TRUE  dimnames=TRUE
agg_thresh1    count=TRUE  weight=TRUE  dimnames=TRUE
sub_all        nrow=28  cols=TRUE  levels=TRUE  num=TRUE  chr=TRUE rownames=TRUE
sub_src_g1     nrow=6   cols=TRUE  levels=TRUE  num=TRUE  chr=TRUE rownames=TRUE
sub_src_tgt    nrow=9   cols=TRUE  levels=TRUE  num=TRUE  chr=TRUE rownames=TRUE
```

`matcharg` is `match.arg`'s partial matching (R-ism 9) and `ambiguous` is its error, so
the *message* is part of the contract too, not just the numbers.
`options$run.time` is deliberately excluded: it is wall-clock and cannot be identical.

`aggregateNet`'s unfiltered branch and `subsetCommunication` for `slot.name = "net"` /
`mode = "single"` are also Rust-backed and identical. Cases the shim still routes to
upstream, loudly (`NOT_PORTED` in `R/modeling.R`):
`datatype != "RNA"` (needs `P.spatial` from `computeRegionDistance`, which requires the
exact k-d tree of §14.5), and `raw.use = FALSE` (needs upstream's `data.smooth`).
`LR.use` reordering, `mode = "merged"`, `subsetCommunication`'s `signaling` / `pairLR.use`
filters, `rankNet`'s `mode = "comparison"`, `datatype != "RNA"`, and `projectData` also stay
in R for now. `raw.use = FALSE` is Rust-backed: it is the same kernel reading a different
matrix, and what is missing is only the projection that produces the matrix. Rust-backed and identical,
messages and return classes included: `aggregateNet`'s filtered branch, `filterCommunication`,
`subsetCommunication`'s eight DEG thresholds together with its `netP` aggregation (which returns
a *tibble*, and is reproduced as one), and `rankNet(mode = "single")` -- see
`docs/SEMANTICS.md`.

Three more bugs, all in code that looked right:

22. **`aggregateNet`'s `net$weight` is an `LDOUBLE` sum across interactions**, and the
    accumulator must stay 80-bit for the whole loop. An accumulator that narrows back to
    `f64` after each add is an `f64` sum wearing a disguise, and it is 1 ulp off at 8 terms
    (2.27e-13 on 1230) -- invisible at 3 or 4 terms, which is why the shorter fixtures
    passed.
23. **`reshape2::melt` on a 3-D array is not in the array's own order.** It builds the label
    frame with `expand.grid`, which varies its *first* argument fastest, so melted rows run
    interaction-major, then target, then source. That order is `subsetCommunication`'s
    return value.
24. **`intersect(want, colnames(net))` sees column *presence*, not non-`NA` values.** An
    all-`NA` `evidence` column is in the output; a table with no `evidence` column at all is
    not. Deriving the columns from the surviving rows instead collapses them to the five
    melt columns, and deriving them from cell values drops `evidence` entirely.

`computeAveExpr` (all three `type` values, R's `intersect` order, the `match.arg` error)
and `subsetDB` (four configurations, the `non_protein` flip, the unknown-key error) are also
Rust-backed and identical. R-isms 22–28 are in `docs/SEMANTICS.md` alongside 12–21.

**Fixture discipline (learned the hard way in 1e).** `gen_expr_golden.R` now dumps the
416 × 6 expression matrix to `tests/fixtures/expr_matrix.tsv` as `%a` hex floats instead
of leaving the Rust test to re-derive it from `set.seed`/`runif`. Two rounds of "every
name resolves, every number is plausible, every value is wrong" came from that
re-derivation. Rules now enforced:

* the fixture is **ground truth in a file**, self-describing (`groups<TAB>n`,
  `genes<TAB>m`, names, a `values` marker, one value line per gene) — never inferred from
  line counts, because R's connection buffering reordered the two halves once;
* a hand-transcribed mirror of a generator literal (`LR_DF`) is checked by a test against a
  record the generator emits for the purpose;
* where a record's key was ambiguous (one name for a pair of cofactors), the generator
  records **both** names;
* `expr_matrix.tsv` is read, not regenerated; `MersenneTwister` is still checked against
  `set.seed(4242)` separately, as the join between the two.

Twenty-one substantive bugs were found *only* by differential testing against R and real
x87, all of them in code that looked correct and produced plausible, uniformly-distributed,
wrong answers:

1. `sample.int`'s loop passes the **shrinking** population to `R_unif_index`; using the
   original `n` yields duplicates.
2. `ceil(log2(dn))` is `bit_length(ceil(dn) - 1)`, not `bit_length(dn)` — off by one at
   every power of two, which shifts the whole stream.
3. R's `Int32` is `unsigned int`, so MT19937's shifts are **logical**; the signed
   variant yields a uniform but entirely wrong stream.
4. `FixupSeeds` overwrites `i_seed[0]` (the index `mti`) with the constant `624`, so a
   port that keeps the LCG value reads out of bounds.
5. R's `mean` is a two-pass corrected mean in 80-bit; `f64` gives `0` where R gives
   `0.33365885416666669`.

Plus **six** in the quantile and database layers, all found the same way:

6. collapse's live path is `dquickselect` (which tests the *reduced* weight, avoiding
   `0 * Inf = NaN`), not the `FQUANTILE_ORDVEC` macro in the same file.
7. `triMean`'s inner `mean` takes `na.rm = FALSE` — the `na.rm` belongs to `fquantile`.
8. `extractGeneSubset`'s "is this a complex?" test, `x %in% symbols == "FALSE"`, selects
   names **absent** from `geneInfo$Symbol` (it only works because R coerces logical to
   character). Read literally, it drops every complex and subunit.
9. `order()` on the annotation *factor* is level order; on the character vector it is
   alphabetical. Ordering after the `as.character` conversion gives a different `nLR1`.
10. `unlist()` on a multi-row data frame is column-major, and the empty cells are
    load-bearing: the human DB contains `IL12AB` = `IL12A, <empty>, IL12B`.
11. `c(agonist, antagonist, co_A, co_I)` is column-wise, not row-wise.
12. A plain gene's expression is **copied verbatim** by `computeExpr_LR`; only complexes go
    through `exp(mean(log(x)))`. They differ (`exp(log(3)) == 3.0000000000000004`), so
    single genes need their own path.
13. The **antagonist** Hill term is `Kh^n/(Kh^n + x^n)`, not the agonist's `x^n/(Kh^n+x^n)`.
    They coincide only at `x == Kh`.
14. In a **cells × genes** R matrix the *gene* index strides (`m[cell, g]` at
    `g * n_cells + cell`), the opposite of the genes × groups matrices used elsewhere in the
    port. Both layouts are correct for their own shape; confusing them reads a transposed
    matrix.
15. `aggregate()` is **group-major** (index `g * n_cols + j`) and the `computeExprGroup_*`
    paths then `t()` it. Picking the wrong stride returns the right multiset of numbers in
    the wrong order.
16. `unlist()` on a data-frame column slice is **column-major** even when the *intent* reads
    as row-major.
17. `aggregate()` **drops** empty factor levels rather than emitting `NA`. `r-core` returns a
    dense buffer with `NaN` instead — a documented divergence, unreachable from CellChat.

18. `P2` and `P3` are each `crossprod` of a 1 x k row, so `P2[a,b] = p2[a]*p2[b]` and
    `P3[a,b] = p3[a]*p3[b]`. Writing `p2[a]*p3[b]` agrees whenever one of them is 1 — i.e.
    for every pair without both an agonist and an antagonist.
19. `P1*P2*P3*P4*P.spatial` associates strictly left to right, and the roundings are
    observable.
20. `P.spatial` is multiplied by `adj.contact` on **every** iteration past `nLR1`, so
    `Pnull` uses `P0 * adj` there. Idempotent only because `adj.contact` is 0/1.
21. `Prob` can exceed 1: the co-agonist is a product of `1 + h` over all subunits, and the
    rank-1 outer product squares it. Upstream behaviour, and changing it would change every
    published CellChat result.

R-isms 12–21 are written up with sources and consequences in `docs/SEMANTICS.md`.

The pattern is consistent enough to be worth stating: **every one of these produces output of
the right length, the right set, and plausible numbers, differing only in order or in the
last ulp.** A set-based comparison would have passed all of them. Each was caught only by a
bit-level differential test against the pinned upstream.

Day-by-day for the remainder:

* **W1a** workspace, `r-core` skeleton, CI (fmt/clippy/`cargo test`/miri-free), MSRV pin.
* **W1b** `rng`: MT19937 + `unif_rand()` + `R_unif_index`; property test:
  `sample.int(n, n)` byte-identical to R for all `n ∈ [1, 5000]` and 10⁴ random `(n, seed)`.
* **W1c** `stats`: collapse type-7 quantile. **Port collapse's arithmetic, not R's.**
  Test: 10⁶ vectors, `identical()` to `collapse::fquantile`.
  Also R type-1 (`trimmedMean`), `median`, `geometricMean`, `thresholdedMean`.
* **W1d** `db`: parse `CellChatDB` (exported to Arrow/Parquet via `arrow`), flatten
  complex/cofactor subunit tables into `Vec<Range<u32>>` over a gene-symbol index.
* **W1e** `expr`: `computeExpr_LR`, `_complex`, `_coreceptor`, `_agonist`, `_antagonist`,
  `geometricMean` — unit parity each.
* **W1f** `aggregate`: L0 implementation, then L1 (zero short-circuit).
* **W2a** `prob` + `kernel`: the `i`/`b` loops, exact R-isms of §4.
* **W2b** `r-bindings` + R shim; first end-to-end parity run on the smoke fixture.
* **W2c–W2d** run the full fixture matrix; triage every mismatch into a test.

**Gate 1:** 100 % of fixtures at rung R1–R8 `identical()`. Any deviation below rung is a
blocking bug with a new regression test.

### Phase 2 — Drop-in shim (1 week)

* `CellChat:::.onLoad` `dyn.load`s the cdylib, generates wrappers from extendr metadata.
* `computeCommunProb(object, ...)` keeps the **exact upstream signature and defaults**,
  extracts slots, calls Rust, writes `object@net$prob/$pval`, `object@images$distance`,
  `object@options$parameter`, `object@options$run.time`.
* `options$parameter` must contain the *same keys in the same order*.
* Message parity: replicate upstream `cat()`/`print()` text under `verbose = TRUE`
  (gated by `getOption("cellchatrs.verbose", TRUE)`).
* Error parity: reproduce the `droplevels` check, the `scale.distance` bound check, and the
  `contact.range`/`contact.knn.k` check with identical messages.
* Interop test: `identical(computeCommunProb_rust(o), computeCommunProb_upstream(o))`
  on every fixture, using `all.equal` at tolerance 0 *and* `identical()`.

**Gate 2:** `identical()` on the full S4 object for all RNA fixtures.

### Phase 3 — Benchmark study (1 week)

Design in §9.2. Deliverables: raw timings CSV, fitted scaling model, speedup table with
CIs, Amdahl decomposition (kernel vs `identifyOverExpressedGenes` vs I/O), and 6 paper-grade
figures. Also profile with `perf` to confirm the model (memory-bound vs compute-bound).

**Memory headroom gate (locked host, §14.8):** 32 GB caps the grid. Before quoting any
`nC ≥ 50 000` number, verify that `data.use` (`nGenes·nC·8 B`) + `avgB` (~34 MB) +
per-pair `Pboot` working set fits with ≥20 % spare, measured with `/usr/bin/time -v`.
Cap the grid at 50 000 cells if it does not; report the cap rather than an OOM.

### Phase 4 — Spatial, comparison, downstream (3 weeks)

* `computeRegionDistance` — replaces `BiocNeighbors::queryKNN(AnnoyParam)` with an exact
  k-d tree (**locked**, §14.4). This changes results, because Annoy is approximate.
  Deliverable: a *divergence measurement* on the real visium fixture — the fraction of
  `queryKNN` neighbour sets that differ, and the resulting spread in
  `net$weight` / `net$count` — published alongside the speedup. Framed as a correction
  (removing an approximation), not as a regression.
* `rankNet(mode = "comparison")`, `rankNetPairwise`, `mergeCellChat`, comparison analysis.
  (`computeCommunProbPathway`, `aggregateNet`, `filterCommunication`,
  `subsetCommunication`'s DEG/`netP` branches and `rankNet(mode = "single")` have since
  landed; this list is the remaining work.)
* `identifyOverExpressedGenes` (Wilcoxon DE) — worth porting; it is the *other* hot spot
  for large data, and `presto`'s exact algorithm must be replicated.

### Phase 5 — Packaging, CI, docs (1 week)

* `R CMD check` clean; `pkgdown` site; vignette showing the measured speedup.
* **GitHub release first** (locked, §14.7), with the release notes stating the R-ism 11
  compatibility decision explicitly, so it is not mistaken for a bug.
* CRAN/Bioconductor submission checklist, prepared in parallel and submitted after
  CellChat maintainer sign-off. CellChat is GPL-3; the port must be GPL-3 too.
* Reproducible artifact: `renv.lock` + `Cargo.lock` + `UPSTREAM.md` SHA + figshare DOIs.

---

## 9. Testing strategy (the rigour requirement)

### 9.1 Fixture matrix (the oracle)

Axes and levels:

| Axis | Levels |
|---|---|
| `type` | `triMean` (default), `truncatedMean`, `thresholdedMean`, `median`; partial-match forms `"tri"`, `"med"` |
| `nC` | 50, 200, 1 000, 5 000, 21 557 (real), 50 000 |
| `K` | 1, 2, 3, 8, 25 (real), 40; plus one cell group with 1 cell; plus one with 0 cells (error path) |
| `nLR` | 1, 5, 50, 500, 2 239 (full), plus hand-built edge cases |
| `nboot` | 1, 2, 10, 100 |
| `population.size` | `FALSE`, `TRUE` |
| `raw.use` | `TRUE`, `FALSE` (needs `data.smooth`) |
| `datatype` | `RNA`, `"spatial"` (with `distance.use` `TRUE`/`FALSE`, `contact.dependent` `TRUE`/`FALSE`) |
| LR structure | single gene; 2–5 subunit complex; agonist only; antagonist only; co-A only; co-I only; all combined; complex with a missing subunit; cofactor not in DB; gene absent from `data.signaling` |
| `Kh`, `n` | (0.5, 1) default; (0.5, 0.5); (2, 2); (1e-3, 1); (1e3, 1) |
| `trim` | 0, 0.1, 0.25 (boundary), 0.3 (error path) |
| data | all-zero; all-constant; single nonzero cell; `NA`/`NaN`/`Inf` injected; integer counts; log-normalised; 99 % sparse; duplicated rownames |

Total ≈ 200–400 configurations. Cartesian explosion is avoided by
**pairwise coverage** (each axis fully crossed at its "typical" point, then varied
one-at-a-time), which gives far more defect coverage per R-oracle minute.

### 9.2 Test layers

| Layer | Tool | What it proves | Count |
|---|---|---|---|
| L0 unit | `cargo test` | each primitive matches R within tolerance | ~200 |
| L0 property | `proptest` | invariants: `triMean` monotone-equivariant under sorting, invariant under input permutation; Hill ∈ [0,1]; `Pval ∈ {k/nboot}` | ~50 |
| L1 RNG | golden | `sample.int` bit-identical to R | 10⁵ cases |
| L1 stats | golden | `fquantile` bit-identical to collapse | 10⁶ vectors |
| L2 kernel | oracle matrix | `Prob`/`Pval` bit-identical end-to-end | ~300 |
| L3 differential | `R` harness | `aggregateNet`, `rankNet`, `subsetCommunication`, `netP` identical | ~100 |
| L4 statistical | R | with a *different* RNG stream, p-value distributions are exchangeable (two-sample KS, and **sign-flip agreement rate** on calls at `thresh=0.05`) | 20 |
| L5 metamorphic | `proptest` | invariant under: cell permutation; cluster relabelling (values permute identically, dimnames follow); gene permutation; scaling `data` by a constant (prob is scale-invariant, since `data/max(data)`); duplicating all cells (prob unchanged) | ~15 |
| L6 regression | golden | known bug fixtures stay fixed | growing |
| L7 end-to-end | R script | run the upstream tutorial verbatim; all `net*` slots and all `subsetCommunication` frames identical | 4 tutorials |

**Invariants asserted on every kernel run** (cheap, catch a lot):

* `all(Prob >= 0 & Prob <= max(Prob))`, `all(is.finite(Prob))` unless inputs are non-finite;
* `all(Pval >= 0 & Pval <= 1)` and `all(Pval * nboot == round(Pval * nboot))`;
* `Pval[Prob == 0] == 1`;
* `dimnames(Prob) == dimnames(Pval)`, `nLR` last-dim size matches `nrow(LRsig)`;
* `Prob` unchanged when `nboot` changes (the observed/point estimate is independent of `nboot`) — **a very strong end-to-end check**;
* `dimnames` order equals `levels(object@idents)`.

**Fuzzing.** A `proptest`-driven config generator (not a fixed list) that emits random
`LRsig` data frames with random subunit/cofactor structures, random sparsity, and random
NaN patterns, and asserts the invariants above plus agreement with R on a 1 % sample
(running R for 100 % of fuzz cases is prohibitive; 1 % sampling with a 5-minute budget is
ample and is itself part of the CI cost model).

### 9.3 CI matrix

| Job | Matrix |
|---|---|
| `rust` | `stable`, `nightly` (miri on `stats` only), × ubuntu/macos/windows |
| `r-parity` | ubuntu + R 4.3/4.4/4.5, × {smoke fixture matrix (≈20 configs), full matrix nightly} |
| `bench` | nightly, self-hosted, pinned CPU (`taskset`), publishes JSON |
| `msrv` | oldest supported R + `extendr` min version |

Nightly full-matrix, PR-time smoke subset (target: PR CI < 6 min, nightly < 40 min).

---

## 10. Benchmark methodology

### 10.1 Dataset suite

* **Synthetic (controlled).** Negative-binomial counts, `nC ∈ {1e3, 5e3, 1e4, 2e4, 5e4}`,
  `K ∈ {8, 16, 25, 40}`, sparsity ∈ {2 %, 7 %, 20 %}, `nLR` ∈ {250, 1000, 2239},
  `nboot ∈ {10, 100}`. Purpose: fit the scaling model with orthogonal factors.
* **Real.** wound (21 557 × 25), humanSkin (7 563 × 12 × 2 conditions),
  visium mouse cortex (spatial), embryonic E13/E14 (comparison).
* **Stress.** Downsampled 50k- and 100k-cell sets from GEO, to reach the regime users
  complain about.

### 10.2 Protocol

1. Pin CPUs: `taskset -c 0-7` (8 physical cores); SMT off. Set
   `RAYON_NUM_THREADS`/`CELLCHATRS_THREADS` explicitly. Record `lscpu`, governor,
   `R.version.string`, commit SHAs.
2. Discard one warm-up run (page cache cold for the 209 MB matrix).
3. ≥5 timed repetitions; report **median** and bootstrap 95 % CI of the ratio.
4. Report **peak RSS** (`/usr/bin/time -v`) — a 2× memory regression is a failure.
5. Report **CPU time** as well as wall time, to separate parallel efficiency from work done.
6. Speedup = `median(t_R) / median(t_rust)` on identical inputs, with the CI propagated.
   Report **paired** speedups (same run, alternating) to cancel drift.
7. Cold-cache variants reported separately (`drop_caches` if permitted, else a 4 GB
   cache-polluting read loop).
8. End-to-end speedup = whole `createCellChat → … → computeCommunProb` pipeline,
   including `identifyOverExpressedGenes` and `subsetData`, so the Amdahl fraction is honest.

### 10.3 Reporting

`t.csv` (raw), `summary.csv` (median/CI/speedup/RSS), `scaling.json` (fitted model),
six figures: (1) speedup vs `nC` (log-log, one line per `nGenes`); (2) speedup vs `nLR`;
(3) speedup vs `nboot`; (4) R cost breakdown stacked bars; (5) strong-scaling efficiency
(1…8 threads); (6) Rust↔R agreement scatter (parity plot, showing bit-identity as a
diagonal of exact zeros).

---

## 11. Risk register

| ID | Risk | Likelihood | Impact | Mitigation |
|---|---|---|---|---|
| **R1** | `collapse` quantile arithmetic not reproducible | Medium | High (drops us to rung "≤2 ulp") | Collapse source is open; already measured the discrepancy is ≤4.4e-16 and localised to type-7 interpolation; port the C++ path directly and prove with 10⁶ vectors. Fallback rung is acceptable and publishable. |
| **R2** | `NMF`/`reticulate`/`BiocNeighbors` block a clean CellChat install | **High — mitigated** | Low | Already worked around: sourcing `R/modeling.R` into an environment plus a minimal S4 class is enough to run upstream `computeCommunProb` (it produced every number in §2.5). Only Phase 2 needs the real package installed. |
| **R3** | Annoy is approximate → spatial results diverge from upstream | **Certain — accepted** | Low | Decided: exact k-d tree (§14.4). Phase 4 must *quantify* the divergence on the real visium fixture and frame it as a correction. Residual risk is only that the measured divergence proves large enough to need per-dataset review. |
| **R4** | End-to-end speedup diluted by `identifyOverExpressedGenes` | High | Medium (only for the headline number) | Measure and publish the decomposition; port the DE step in Phase 4. Never claim end-to-end speedup without this. |
| **R5** | `log`/`exp` libm differences across platforms | Medium | Low | Parity CI on 3 OSes; if it bites, use a crate reproducing glibc's libm, or accept ≤1 ulp. |
| **R6** | S4 slot/attribute mismatch in the shim | Medium | Medium | `identical()` on the whole object is the gate; `slotNames`/`attributes` diff in the harness output. |
| **R7** | Memory blow-up at 100k cells | Medium | Medium | Streaming aggregation over gene tiles; never materialise `Pboot` for all pairs; enforce a configurable cap with a clear error. |
| **R8** | Oversubscription when combined with `future` | Medium | Low | Detect `nbrOfWorkers()`; scale the rayon pool accordingly. |
| **R9** | Upstream diverges (v3 / SpatialCellChat) | High over time | Medium | Pin a commit; keep the oracle suite parameterised by version; mirror `SpatialCellChat` separately. |
| **R10** | Slower than hoped because R's cost is not where modelled | Low | High | Phase 0's cost attribution exists precisely to catch this, before writing Rust. |

---

## 12. Milestones and deliverables

| ID | Deliverable | Acceptance |
|---|---|---|
| M0 | Reproducible R baseline + oracle manifest + benchmark harness | CV ≤ 2 %; oracle complete |
| M1 | `r-core` with RNG, stats, expr, aggregate, kernel | Unit + property tests green |
| M2 | `r-bindings` + R shim | `identical()` on all RNA fixtures |
| M3 | Benchmark report + scaling model + figures | CIs, RSS, Amdahl decomposition |
| M4 | Spatial + downstream numerics | Oracle-clean on Tier-3 fixtures |
| M5 | `R CMD check` clean, docs, CI, reproducible artifact | — |
| M6 | Manuscript draft | — |

---

## 13. Publication framing

The paper should claim exactly what the evidence supports:

1. **A drop-in, API-compatible accelerator** for CellChat: bit-identical outputs at
   `identical()`-level tolerance, verified on a 200–400 configuration matrix and 4 upstream
   tutorials.
2. **A measured, CI-backed speedup** with a fitted scaling model, and an explicit Amdahl
   decomposition showing what is *not* accelerated.
3. **An honest correction**: upstream recommends downsampling "to further speed up";
   we show it becomes unnecessary, and quantify the result on 50k–100k cells.
4. **A note on spatial**: exact KNN replaces Annoy, removing an approximation.

Avoid claiming a specific multiplier until M3 — state the model as a prediction in M1/M2
and only publish measured numbers in M3.

---

## 14. Decisions taken

All questions resolved on 2026-09-30. This section is now normative: the rest of the plan
is written to these decisions, and changing one here invalidates the sections that
reference it.

| # | Decision | Consequence for the plan |
|---|---|---|
| 1 | **Scope = numeric core + R shim.** Port `modeling.R` in full plus the numerics in `utilities.R` / `analysis.R`; leave `visualization.R` and `app.R` in R. | §1 tiers are final. Effort ≈ 8 weeks. The paper's claim is "same numbers, faster, plus the numerics we also accelerated", not "CellChat in Rust". |
| 2 | **Parity target = bit-identical.** `identical()` on the whole S4 object is the acceptance gate for every fixture. | §4 rungs R1–R11 are all **Exact**; the "1 ulp" and "Tol" rungs exist only as documented contingency, not as a target. §6.1 **L5 (f32 fast path) is removed** — it forfeits R6. Reassociating the aggregation is forbidden; operation order is part of the contract. |
| 3 | **R-ism 11 = replicate the error.** `computeExpr_complex` must abort with the identical message when a complex subunit is absent from `data.signaling`. | Strict parity, and the bug is reported upstream separately. The `r-core` error type needs a `MissingSubunit` variant whose R-side rendering is byte-compared against upstream's message. |
| 4 | **Spatial KNN = exact k-d tree.** Annoy is approximate; we replace it. | S1 stays a **documented divergence**, and Phase 4 must additionally *measure* the divergence: fraction of neighbour sets differing from Annoy on the real visium fixture. This is reported as a correction, not a regression. |
| 5 | **Binding = `extendr`.** Proven working end-to-end on this host. | §5 architecture is fixed. No Python/reticulate dependency. |
| 6 | **Upstream pinned to `main` @ `75253cd0c9e68410e6e721a6d3a0419a1d7e358f`** (v2.2.0.9001, 2026-03-04). | See `UPSTREAM.md`. The oracle manifest and a `test_pin` check record the SHA; re-pinning is a reviewed act. SpatialCellChat / v3 is deferred. |
| 7 | **Distribution = GitHub first, then CRAN.** | Phase 5 ships a GitHub release; CRAN submission follows CellChat maintainer sign-off. The paper is prepared in parallel, not blocked on CRAN. |
| 8 | **Benchmarks on this host, pinned.** 8c/16t, 32 GB, AVX2. | The 100k-cell row is bounded by RAM: `nGenes × nC × 8 B` for `data.use` plus a 34 MB `avgB` plus the `Pboot` working set must fit in 32 GB. Phase 3 must verify headroom explicitly before quoting any 100k number, and cap the grid at 50k if it does not. |

### Consequences that are easy to miss

* Decision 2 (bit-identical) makes the **`Option<f64>`-style guard rails illegal** in the
  hot path: no fast paths that change association, no `f32`, no FMA contraction that
  reorders a reduction. Build with `target-feature=-fma` where the compiler would
  otherwise contract, and add a test that asserts no FMA in the aggregation (compare
  against a `-C target-cpu` build).
* Decision 2 + R10 (R's `LONG_DOUBLE` accumulation) means the 4-element mean in
  `triMean` needs a real 64-bit-mantissa accumulator in Rust, not an `f64` one. Budget
  it; it is a known ~30-line piece, but it is not optional.
* Decision 3 turns a "fix" into a compatibility feature, which reviewers often
  misread. The Phase 2 release notes must state it explicitly.
### Phase 2g — `computeRegionDistance`, and a gate that was measuring nothing

Two things landed together, and the second is the more important.

**`computeRegionDistance` and `computeCellDistance` are ported.** The arithmetic, the neighbour
query and the R-side marshalling each have their own gate, because each has a different oracle:

| piece | oracle | gate |
|---|---|---|
| `mean(x, trim, na.rm)`, `fdist`, the exact k-d tree | R itself, and exhaustive search for the tree | `src/rust/crates/r-core/tests/spatial_parity.rs` (21 tests) |
| `computeRegionDistance`'s arithmetic | upstream's body with the neighbour *query* substituted | `src/rust/crates/r-core/tests/region_parity.rs` (10 tests, 13 fixtures) |
| the R shims | upstream for `computeCellDistance`; the corpus for `computeRegionDistance` | `computeCellDistance` / `computeRegionDistance` blocks of `check_identical.R` |

Upstream's `computeRegionDistance` **cannot be run here at all**: it evaluates
`BiocNeighbors::queryKNN(..., AnnoyParam())` unconditionally and the package is not installable
without network. Annoy is randomised and approximate, so upstream is not an oracle for a neighbour
query regardless. The corpus therefore lifts upstream's body verbatim with two expressions replaced
by exhaustive search, and **each fixture records its own nearest-neighbour margin `d1/d2`** -- the
substitution is only sound when the runner-up decisively loses to the nearest neighbour, and
`the_recorded_margins_hold_up` re-checks the recorded number in Rust so a later layout edit fails
in the suite and not only in the generator. The 13 fixtures span 2.4x to 590x.

Six upstream behaviours that read as bugs turned out to be behaviours, and are documented in
`docs/SEMANTICS.md` with the fixture that pins each. The two that would most plausibly be
"corrected" by a future reader:

* **`d.spatial` is symmetrised outside the `if (do.symmetric)` guard**, so it is symmetric either
  way while `adj.contact` is not.
* **A declared level with no cells gets the wrong row name.** Upstream's loop counter picks the
  output slot and a *compacted* `level.use` picks the level, so A's row is written into slot 1 and
  labelled with the first level.

**The gate was vacuous for every `computeCommunProb` comparison.** `computeCommunProb` gated its
Rust path on `nzchar(Sys.getenv("CELLCHATRS_FALLBACK", "0"))`, and `"0"` is four characters long,
so the predicate was always `TRUE` and every call was delegated to upstream. The gate compares the
shim against upstream, so a fully green run had been comparing upstream with itself, and
`parity.json` recorded those quantities at a passing rung.

The predicate is now an explicit affirmative, `check_identical.R` **stops** if the variable is set
before comparing anything, and `tests/parity/check_rust_path.R` proves the kernel is reached from
outside the package by pointing `options$db` at a nonexistent directory and requiring the
Rust-only error. With the Rust path genuinely live, the previously-hidden `lr_col()` bug appeared
immediately -- a zero-length LR vector is a Rust panic, not a missing optional -- and is fixed.

`computeCommunProb`'s spatial branch is now Rust-backed end to end, which is the point rather than
a completeness tick: falling back would have resolved the neighbours through `AnnoyParam`, the
approximation this port exists to remove.

### Phase 2h — the spatial branch of `computeCommunProb` is gated

`tests/parity/check_identical.R` now has a ten-case spatial block. The reference is upstream's whole
body with `computeRegionDistance` replaced **in its own environment** by the exact-neighbour oracle
from `tests/parity/exact_neighbour_oracle.R`, which is the file `gen_region_golden.R` also sources.
Sharing one copy is what makes the comparison honest: two copies that drifted would be different
functions with the same name, and the gate would compare the port against an accident.

The oracle file is upstream's text with two expressions replaced by exhaustive search -- the
`findKNN` call and the `queryKNN` call -- and nothing else. Everything the block actually tests is
therefore upstream's own code:

```
sp_distance_use       prob=TRUE  pval=TRUE  parameter=TRUE  msgs=TRUE  (64 positive)
sp_no_distance_use    prob=TRUE  pval=TRUE  parameter=TRUE  msgs=TRUE  (64 positive)
sp_scale_too_small    error=TRUE   identical stop(), message embeds format(1/d.min, digits = 2)
sp_forced_contact     prob=TRUE  pval=TRUE  parameter=TRUE  msgs=TRUE  (56 positive)
sp_all_contact        prob=TRUE  pval=TRUE  parameter=TRUE  msgs=TRUE  (56 positive)
sp_mixed              prob=TRUE  pval=TRUE  parameter=TRUE  msgs=TRUE  (64 positive)
sp_contact_false      prob=TRUE  pval=TRUE  parameter=TRUE  msgs=TRUE  (64 positive)
sp_k_min_too_large    prob=TRUE  pval=TRUE  parameter=TRUE  msgs=TRUE  (0 positive)
sp_one_sample         prob=TRUE  pval=TRUE  parameter=TRUE  msgs=TRUE  (64 positive)
sp_knn_k              prob=TRUE  pval=TRUE  parameter=TRUE  msgs=TRUE  (64 positive)

IDENTICAL: 0 failing comparisons out of 359, over 9 configurations
```

`sp_mixed` is the case that earns its place: it is the only fixture that reaches the fourth `nLR1`
branch, `nLR1 <- max(which(annotation %in% diffusible))`, where the contact-dependent and diffusible
pairs get *different* `P.spatial`. The three others cover `contact.dependent.forced`, the
contact-only branch, and the `contact.dependent = FALSE` fallback, and each has a distinct `cat`
text that is compared.

Building this found three bugs, all of which had been invisible:

1. **`options$parameter` recorded the RNA branch's constants unconditionally**, so every spatial run
   reported `distance.use = NULL` and `interaction.range = NULL` while having just used `TRUE` and
   `100`. It now records the values in force.
2. **A zero-length LR vector is a Rust panic, not a missing optional.** `as.character(NULL)` is
   `character(0)`, and a column-subset `LR.use` produces one. `lr_col()` normalises every LR
   metadata vector to `nLR` entries with `NA` as `""`.
3. **Lazy argument evaluation** meant the lengths the kernel reported were the lengths at *force*
   time, not at call time. See `docs/SEMANTICS.md`.

### Phase 2i — `rankNet(mode = "comparison")`

`mode = "comparison"` is now Rust-backed and gated, closing the largest pure-numeric gap left in
`R/analysis.R`. `rankNetPairwise` and `mergeCellChat`'s object plumbing remain.

**What moved to Rust** (`src/rust/crates/r-core/src/ranknet.rs::comparison_flow`): the per-comparison
`apply(prob, 3, sum)`, the `-1/log` transform, the **pooled** degenerate reassignment, the union of
pathway vocabularies, the per-comparison threshold cut and axis filters, and the raw relative
ratios. Nine fixtures, bit for bit.

**What deliberately stayed in R**, and why:

* `as.numeric(format(., digits = 1))`. `format`'s `digits` is *vector-dependent*: `format(1.2,
  digits = 1)` is `"1"`, `format(c(0.04, 1.2), digits = 1)` is `c("0.04", "1.20")`. The rounded
  value of an element is not a function of that element, and the `order()` below sorts on the
  formatted values. `signif(x, 1)` -- the obvious reimplementation -- gets the first case right and
  the rest wrong. R's `format` is exact and is one line away.
* the `order()` itself, for the same reason plus stability.
* the data-frame assembly: the row-name-indexed `df[[i]][pair.name[[i]], 2] <- ...` that leaves an
  absent pathway at 0, `factor()` of `name`, `rev(levels(group))`, and `rbind`'s uniquification
  (observed to append a bare `1`, no separator).

**The plot.** `identical()` on `gg.obj` cannot hold: `geom_bar`'s `position` holds a quosure whose
environment is the one the code ran in, and upstream's is sourced into `cellchatrs_upstream_env()`
while the shim's runs in the package namespace. The gate compares `ggplot_build(gg)$data` -- the
rendered geometry -- and it matches. PLAN.md 14.1 keeps visualization in R, and this is the
concrete reason the claim is about geometry rather than object identity.

```
rc_two_disjoint           value=TRUE  built=TRUE  (10 x 5, rows X,B,C,D)
rc_pooled_degenerate      value=TRUE  built=TRUE  (10 x 5, rows X,four,three,one)
rc_three_comparisons      value=TRUE  built=TRUE  (27 x 6, rows h,g,h,i,a)
rc_four_comparisons       value=TRUE  built=TRUE  (20 x 7, rows e,a,d,c)
rc_count_measure          value=TRUE  built=TRUE  (10 x 5, rows d,e,a,c)
rc_thresh_empties_second  error=TRUE   identical error text
rc_filters                value=TRUE  built=TRUE  (6 x 5)
rc_order_ties             value=TRUE  built=TRUE  (12 x 5, rows v,u,z,y)
rc_all_degenerate         value=TRUE  built=TRUE  (4 x 5)

IDENTICAL: 0 failing comparisons out of 359, over 9 configurations
```

**Two oracles, and where they disagree.** The corpus records the error upstream's *lifted body*
produces; the gate's oracle is upstream's *exported `rankNet`*. For `filters` and `all_degenerate`
the two disagree -- the harness stops, the function does not -- so the gate asserts that the two
sides behave the *same* and prints the unreproduced harness record rather than either failing or
dropping the case. The Rust suite asserts the harness record verbatim, since the harness is the
right oracle for the kernel.

Four bugs found on the way, all of which had been invisible:

1. **`rankNet` routed every comparison-mode call to upstream.** `needs_upstream` tested
   `any(comparison != comparison[1])`, and `comparison` *defaults* to `c(1, 2)` -- so the condition
   was TRUE by default and upstream won unconditionally, for `mode = "comparison"` as well as
   `single`.
2. **`pSum < 0` flags every pathway whose total exceeds 1,** not only the one that totals exactly
   1. `-1/log(x)` is negative for all `x > 1`, and a total above 1 is ordinary for a `k x k` sum of
   probabilities. The source comment claimed the arm was unreachable.
3. **The shim was missing the zero-dropping loop,** which lives in the tail *after* the `do.stat`
   block and so is easy to overlook.
4. **`options$parameter` was recorded from the RNA branch's constants,** so a spatial run reported
   `distance.use = NULL` having just used `TRUE`.

### Phase 3a — the configuration matrix, and the two bugs it found immediately

The definition of done asks for a 200-400 configuration matrix with pairwise coverage across
thirteen axes. `check_identical.R` answers a different question ("does each function work"), so this
is a separate artifact: `tests/parity/matrix.R` (axes, candidate pool, greedy covering selection,
coverage report), `tests/parity/matrix_objects.R` (the object family), and
`tests/parity/check_matrix.R` (the runner and the JSON report).

**The matrix.** Fourteen axes -- type, nC, K, nLR, nboot, population.size, raw.use, datatype, LR
structure, Kh, n, trim, pathological data, and a `scale` axis -- with 1051 distinct axis pairs. A
deterministic greedy selection covers all of them in **131 configurations**, and the matrix is padded
to 200 with configurations that overlap the existing ones *least*, so the padding adds breadth
rather than duplicating rows. `matrix.R` runs standalone and prints the coverage; no RNG anywhere,
so the same axes give the same matrix.

Two details that are easy to get wrong and would have made the coverage claim wrong:

* **Pairs are keyed by `(axis, level)`, not by the two level strings.** Keying on the levels alone
  collapses distinct required pairs -- `paste0(trim, "|", n)` and `paste0(Kh, "|", type)` can both
  be `"0.5|1"` -- and then one configuration is credited with covering both. Keying on levels
  reported 616 pairs where there are 875.
* **The greedy gain is an intersection with what is *still uncovered*.** The obvious
  `setdiff(pairs, c(covered, wanted))` subtracts the whole required set, scores every candidate at
  zero, and stops on the first iteration -- leaving a "matrix" that was 100% padding and reported
  429/616 pairs as covered.

Coverage is recomputed from the **effective** axes -- what each configuration actually ran with,
after any adaptation. A builder that forces `K = 1` for a one-cell object records that, so the row
cannot be counted as covering `K = 2, 3, 4`.

**Two bugs, both invisible until now.**

1. **`data.use <- data/max(data)` was missing from the kernel entirely.** Every fixture in the suite
   drew from `runif(0.01, 1)`, and `runif` reaches `1` often enough that `max(data) == 1` in all of
   them -- the division was the identity, so its absence changed nothing and the differential gate
   stayed green. The step's purpose is scale-*invariance*, and without it the port's `Prob` moved by
   more than an order of magnitude on the same data divided by two. Fixed in
   `r_core::prob::scale_by_max`, with a new `scale` axis in the matrix and a *direct* invariance
   check in the runner -- pairwise coverage cannot express a relation between two rows, so it would
   not have caught a regression here.
2. **`F80` had no representation for `NaN` or `Inf`**, and `from_f64` mapped both to zero. Every
   reduction in the crate accumulates through `F80`, so `sum(c(1, NaN))` returned `1` -- not a value
   R can produce. Combined with R's `if (sum(P1_Pspatial) == 0)`, where `NaN == 0` is `NA` and
   `if (NA)` raises, the kernel was *returning an all-`NaN` network* where upstream refuses. `F80`
   now carries both specials in sentinel exponents, `r_max` implements R's `max` rather than Rust's
   `maxNum`, and `KernelError::IfNa` raises upstream's message.

One "bug" was a bug in the fixture, and finding that is also worth recording: `dup_rownames` was
applied before the gene subsetting, so R's `d$expr[genes, , drop = FALSE]` -- which indexes by
*name* -- raised "subscript out of bounds" inside the builder. It now applies to the final matrix,
where a duplicate row name can also *remove a subunit from the universe under its own name*, which
is how that fixture legitimately reaches R-ism 11.

**Two more, both found by the matrix and both in the same function.** `thresholdedMean` is the only
`FunMean` that branches on a statistic instead of a value, and the statistic is not what its name
suggests: `Matrix::nnzero` returns **`NA`** for any vector containing a missing value, rather than a
smaller count, so `percent` is `NA` and `if (percent < trim)` raises. The port counted non-zero
entries and skipped the missing ones -- a reasonable reading of the name and of the `na.rm = TRUE` in
the signature -- and a unit test had locked that reading in.

Fixing that exposed the second half. Upstream's `aggregate` at `modelling.R:115` runs `FunMean` over
**every** gene, and the port does too, but then reads back only the interacting genes and discards
the rest -- including the `NaN` that upstream would have raised on. So a raise on a gene that no
interaction references became an all-zero `Prob`. A `NaN` in the observed aggregate is now a raise for
`thresholdedMean` and only for `thresholdedMean`; for `triMean` and `median` it is a value, and
upstream's own `if (sum(P1_Pspatial) == 0)` is what turns it into the same message. The port also ran
the bootstrap before the observed aggregate, the reverse of upstream's order, so its error precedence
was wrong independently of the arithmetic.

**A report that silently dropped an axis.** The JSON writer emitted `names(cfg)[-1]`, and the first
element happened to be `type`. Every row in `matrix_report.json` was therefore missing its `type`, which
is what made two rows with the same tag indistinguishable and sent the investigation down a false
trail for a while. A report that drops one axis cannot be used to tell two rows apart, and it dropped
the one that mattered.

**The matrix passes in full.** 200 configurations, 1051 of 1051 axis pairs covered, scale-invariance
asserted, zero non-passing rows -- 158 `identical` and 42 `error-equal`, where "error-equal" means
upstream itself raises (missing subunit, all-zero matrix, duplicated row name) and the port raises the
same message. The runner now carries a `type` field in every row, which the JSON writer had been
dropping.

Gates re-verified after the `thresholdedMean` and error-precedence fixes: `cargo test --release` 272
tests across 22 binaries; installed-shim gate 359/359; `parity.json` 124/129 quantities at a passing
rung with 5 untested; `R CMD check` `Status: OK`.

### Phase 3b — `rankNetPairwise`: the port, and what a pass-through costs

`rankNetPairwise` and `mergeCellChat` were the two functions in the required surface that the shim
did not provide at all -- `CellChat:::rankNetPairwise` did not exist, and both calls went
straight to upstream. `rankNetPairwise` is now ported: the shim owns the function and the ordering
`order(pval, -prob)` per `(i, j)` group pair is Rust (`ranknet_pairwise_orders`, over the existing
`order_f64_multi`). Everything else in its body -- the `data.frame()`, its `row.names`, the nested
`list()`, `names(temp) <- colnames(prob)` -- stays in R, on the same reasoning that keeps
`format(x, digits = 1)` and `rankNet`'s row assembly there: it is language behaviour, not numerics.

`mergeCellChat` remains R throughout -- it lives in `CellChat_class.R` rather than `analysis.R`, and
it is S4 slot assembly with no arithmetic, so 14.1 covers it. It is listed as its own parity quantity
at the **pass-through** rung rather than folded into `rankNet_comparison`, so the report says which
half is done.

Its gate, `tests/parity/check_merge.R`, is deliberately not a comparison of `mergeCellChat` with
upstream's `mergeCellChat`: there is only one implementation, so that comparison would report
`IDENTICAL` while testing nothing. It asserts the merge's own contract instead, and then runs the
ported `aggregateNet` / `subsetCommunication` / `filterCommunication` / `rankNet` on a per-dataset
net lifted back out of a merged object. That second half is the part that can fail, and writing it
found a real defect: the port returned `interaction_name` as `character` where upstream returns a
**factor** -- `var.convert` factors every column that came from a dimname, and `interaction_name`
is the third dimension's. Every value in `netfiltered_golden.txt` still matched, because the corpus
recorded no column *classes* at all. The corpus records them now, `NetTable` carries
`interaction_levels`, and two Rust tests pin them, so `cargo test` catches a regression rather than
needing the R gate. Details in `docs/SEMANTICS.md`.

### Verification state

| gate | result |
|---|---|
| `cargo fmt --all --check` | clean |
| `cargo clippy --workspace --all-targets -- -D warnings` | 0 findings (was 325) |
| `cargo test --workspace --release` | 318 passing, 26 binaries, 0 failing |
| `check_identical.R` (the acceptance gate) | IDENTICAL, 0/418 failing, 9 configurations |
| `check_merge.R` | IDENTICAL, 0/17 failing |
| `check_cli.R` (the standalone binary) | IDENTICAL, 0/231 failing, 14 configurations |
| `stat_equiv.R` | EQUIVALENT, 0/130 failing, 120 independent seeds |
| `metamorphic.R` | 40/40 |
| `check_matrix.R` | 200 configurations, 1051/1051 axis pairs, scale invariant |
| `scripts/codegen_variants.sh` | default / `+fma` / `target-cpu=native` byte-identical |
| `R CMD check` | Status: OK |
| `parity.json` | 159/159 at a passing rung, 0 untested |

### Centrality: split port (degrees/strength/betweenness in Rust, solvers in igraph)

`analysis.netP_centrality` was the last untested quantity and is now gated. Upstream's
`computeCentralityLocal` computes eleven measures; four are iterative solvers
(`hub_score` deprecated in igraph 2.0.3, `authority_score`, `eigen_centrality`/ARPACK,
`page_rank`/PRPACK) whose output disagrees with itself across runs on identical input --
measured, and asserted by the gate so the justification cannot rot. Bit-parity with a
nondeterministic oracle is meaningless, so those four call igraph directly (identical by
construction). The other seven are pure functions, ported to `r-core::centrality` and verified
bit-for-bit against installed igraph 2.3.4 on a 77-case corpus: unweighted degrees, weighted
strengths (plain sequential `f64` in edge-ID order -- R's long-double `rowSums` agrees only
3802/4669), and betweenness (Dijkstra with dist-plus-one encoding, exact 2-way-heap tie rules,
epsilon comparisons at 1e-10, Brandes accumulation), plus igraph's exact NA/positivity messages
and the tiny-weights warning. `sna::flowbet`/`infocent` run for real on both sides (upstream's
own `tryCatch`, copied verbatim). `check_centrality.R`: 17/17, and the tutorial exercises the
whole chain end to end.

### The standalone CLI

The objective requires `r-core` to ship "as a standalone library/CLI so the numerics are usable without
R", and there was no CLI. `src/rust/crates/cellchatrs-cli` is one: a `cellchatrs` binary that reads a
self-describing input file and writes the `Prob`/`Pval` networks, the averaged expression and
`aggregateNet`'s two matrices, with no R in the process. Subcommands are `run`, `describe`, `mean` and
`version`.

`tests/parity/check_cli.R` gates it against **pinned upstream** over 14 configurations -- all four
`type.mean` values (two as `match.arg` prefixes), `population.size` both ways, `nboot` 1..9, `Kh` across
six orders of magnitude, `n` at 1 and 2, `raw.use` both ways, two seeds -- requiring `identical()` on
`Prob`, `Pval`, `computeAveExpr`'s means and `aggregateNet`'s matrices. 0 of 231 comparisons failing.

It found three defects in the CLI and one upstream asymmetry that had been mislabelled:

- `raw.use = FALSE` is a **matrix choice and nothing else** -- there is no library-size normalisation
  anywhere in `computeCommunProb` -- and the CLI had invented one, which saturated `P1` to 1 and made it
  raise where upstream succeeded;
- `computeAveExpr` has no `raw.use` parameter, so it and `computeCommunProb` read *different* matrices
  for the same object, and the CLI reported the kernel's aggregate under the name `ave_expr`;
- `Prob` is **not bounded by 1** -- 16 of 96 values exceed 1, the largest 2025.1 -- so a range invariant
  written for the gate failed against output the port reproduces bit for bit. `Prob` is a score;
- the fixture generator originally defined `aggregate`, which shadowed `base::aggregate` inside the
  sourced upstream and made every call fail with `unused argument (FUN = FunMean)`, while the run
  reported success.

The three untested quantities are all blocked on an external dependency rather than on work:
`analysis.netP_centrality` (igraph 2.3.4 is now installed locally and the port is the remaining work),
`prob.spatial_divergence` (BiocNeighbors/Annoy), `utilities.identify_over_expressed_genes_fast`
(presto).

Three things the verification work turned up that were not previously known, all recorded in
`docs/SEMANTICS.md`:

- **The suite is demonstrably sensitive to rounding.** A control experiment -- `F80::mul` rounded
  through an `f64`, i.e. a 53-bit accumulator -- fails 9 tests across 6 binaries. "311 tests pass" on
  its own says nothing about whether the tests can see rounding at all; this says they can.
- **`interaction_name` in the melted table is a `factor`, not a `character`.** A real defect: the
  corpus recorded no column *classes*, so every value comparison still passed. The corpus records them
  now, `NetTable` carries `interaction_levels`, and two Rust tests pin them.
- **`seed.use` changes `Pval` and not `Prob`.** Upstream permutes the group labels and not the data,
  so the null distribution depends on the seed while the point estimate cannot. The obvious version of
  a seed-stream test is therefore vacuous, and the control in `stat_equiv.R` exists because this
  gate's first version got it backwards.

**Three bugs, and the reason the first one survived so long.** The kernel misread R's `dim()` as a
list of slice lengths, so it returned `__error`, and the shim's first version *silently answered
from upstream* on any mismatch it could not interpret. Upstream's answer is identical to upstream's
by definition, so all eighteen differential cases passed for a function that was not ported at all.
The general lesson is now structural rather than a note: the fallback warns, and the gate counts
fallback warnings as failures. A gate that cannot tell a port from a pass-through is not a gate.

The second was the slice *order*: a `k x k x n` array flattens column-major, so a fixed `(i, j)`
slice is a strided gather with stride `k * k`, not `chunks_exact(n)`. That version passed a
one-slice test 200 times and failed 539 of 540 slices at `c(3, 3, 6)`.

The third was in `order_f64_multi`, which `rankNet` also uses: a missing value in a **second** key
compared as equal rather than sorting last within its group, and an all-missing **primary** key
stopped the comparison instead of falling through to the remaining keys. Both are wrong against R
(`order(c(1,1,2), c(NA,5,3))` is `2 1 3`; `order(c(NA,NA,NA), c(3,1,2))` is `2 3 1`) and both are
invisible without repeated keys, which is why the gate now has a case built around exactly that:
repeated `pval`, one `NaN` probability.

Also corrected: `check_matrix.R`'s JSON writer emitted `names(cfg)[-1]`, and the first element
happened to be `type` -- so every row of `matrix_report.json` was missing its `type`. A report that
drops one axis cannot be used to tell two rows apart.

Gates after this phase: `cargo test --release` 274 tests across 22 binaries; installed-shim gate
**377/377**; matrix 200 configurations, 1051/1051 axis pairs, zero non-passing rows;
`parity.json` **129/134** at a passing rung with 5 untested; `R CMD check` `Status: OK`.

### Next: the three reachable quantities, ranked

`analysis.mergeCellChat`, `utilities.identify_over_expressed_genes_dataset` and
`analysis.netP_centrality` are the remaining gaps that do not first need a package installed.

`identifyOverExpressedGenes(group.dataset = ...)` is worth correcting in the ledger first, because
it was filed as blocked and is not. `group.dataset` appears in **both** halves of upstream's
`if (do.fast)` -- the presto half at `utilities.R:429-484` and the Wilcoxon half at `:512-520`. Only
the first is unreachable without `presto`, so `do.fast = FALSE` with `group.dataset` set is testable
today. What the port is missing is the dataset-comparison *selection*: `cell.use1`/`cell.use2`
chosen by dataset rather than by group complement, the `pos.dataset` validation (a `cat()` followed
by a bare `stop()`, so the message is R's default and the failure is an error with no text), the
`group.DE.combined` variant that pools cells across groups, and the `markers.all$datasets` factor with
its `order(datasets, pvalues, -logFC)`. The `group.dataset = NULL` path is already gated.

`mergeCellChat` is 115 lines of S4 slot combination and `net` reindexing across a list of objects,
with no arithmetic to move -- the same category as `rankNetPairwise`'s data-frame construction, where
the right answer turned out to be "port the one `order`, keep the plumbing" rather than "port it all".

`netP_centrality` is the only one of the three with real numerics, and it is blocked on `igraph`:
`subsetCommunication(measure = "centrality")` and `rankNet`'s DEG-weighted centralities both route
through `igraph::betweenness`, `closeness` and `eigen_centrality` on a built graph. Reproducing
igraph's exact outputs for those three is a different proposition from reproducing CellChat's
arithmetic, and it should not be attempted by guessing.

The other two untested quantities are genuinely dependency-blocked and should stay that way until the
package is available: `prob.spatial_divergence` needs `BiocNeighbors`/`Annoy` to measure the
approximate query against the exact k-d tree on real spatial data, and
`utilities.identify_over_expressed_genes_fast` needs `presto`, which is a different algorithm rather
than a faster one.

### Phase 3c — `group.dataset`: a gap that was filed for the wrong reason

`identifyOverExpressedGenes(group.dataset = ...)` was recorded as blocked on `presto` because it
sits next to `do.fast = TRUE` in the ledger. It is not blocked: upstream honours `group.dataset` in
both halves of `if (do.fast)` -- `utilities.R:429-484` for presto and `:512-520` for the Wilcoxon
test -- and only the first needs presto. It is now ported and gated, and `parity.json` is at
**134/138** with 4 untested.

The port's shape came out of the code rather than a design: the branch changes only how
`cell.use1` / `cell.use2` are chosen per group, so `identify_over_expressed_selected` takes the
selection as a parameter and both the plain and dataset paths share the numerics. That is the second
time this turn's work has reduced to "one primitive, parameterised" -- the first was
`rankNetPairwise`'s single `order`.

Four things the branch needed that are worth keeping as rules:

* `toString` collapses **every** non-positive dataset into one level, so four datasets still give two
  levels, one named `"D2, D3, D4"`. The `datasets` factor and the row order follow from that name.
* A group with exactly one percentage-passing feature produces **no** markers: `apply(X, 1, FUN)` on a
  one-row matrix returns an unnamed scalar, `FC[features]` is `NA`, and `features.diff` is
  `character(0)`. Caught by the `idents.use = "g1"` cases, which a threshold sweep never samples.
* An empty `cell.use2` has to *fail* with `-markers.all$logFC : invalid argument to unary operator`,
  after upstream has added a `datasets` column to a 0 x 0 frame. The shim keeps a 0 x 0 frame on this
  path and the collapsed 1 x 1 one otherwise, because the shape is part of the contract.
* `rbind` is not associative. Upstream rbinds one group at a time and R's row-name uniquification
  falls back to a bare digit once `x.1` is taken, so `do.call(rbind, frames)` renames two of
  forty-one rows differently under `only.pos = FALSE` while every column matches.

A shim bug of a different kind is also fixed: the frame loop guarded on the kernel's
`n_before_only_pos` counter, which over-estimates once `idents.use` has filtered rows away, and on an
empty filtered result the run-detection produced an index of `c(1, 0)` and `data.frame` raised
`row names contain missing values`. The guard is `length(res$features) > 0L` now, which is what
upstream's per-group `if (nrow(gde) > 0)` means.

Gates: `cargo test --release` 276 tests across 22 binaries; installed-shim gate **418/418**;
`parity.json` **134/138**; `R CMD check` `Status: OK`. The configuration matrix re-run is in flight
and is the one gate not yet re-verified after these changes.

### Phase 4 — the benchmark, and the bug it found that no test could

The Rust speedup had never been measured. `docs/BENCHMARKS.md` carried a cost model projecting
"~4 s / ~22×" for human skin, and every other deliverable that depends on a number -- the paper's
claims, the vignette, the pkgdown site -- was waiting on it. All three of the authors' datasets are
present locally (`humanSkin.rda`, `wound.rda`, `visium.rds`), so this was never actually blocked.

`bench-runner/fixture.R` is the CellChat input pipeline factored out of `prof_real.R` (which keeps
its own copy, since it is the source of the recorded baselines and editing it would invalidate
them). `bench-runner/bench_real.R` times both sides in one process on one object, discards a
warm-up, takes five timed repeats, reports the median with a 10 000-resample bootstrap CI and
`VmHWM`, pins with `taskset`, and asserts `identical()` on the whole `net` before reporting a
speedup. `bench-runner/analyse.py` emits the tables, the Amdahl fit and the parity SVGs with no
dependencies, so a reader can regenerate the paper artifacts.

**The kernel was single-threaded.** The first thread sweep was flat -- 10.66 s at 1 thread, 10.63 s
at 8 -- and `grep -rn "par_iter\|into_par_iter\|rayon::" src/rust/crates/r-core/src/` returned nothing. The
bootstrap replicate loop was a plain `for`, and 92-95 % of runtime is bootstrap aggregation. A
`rayon` dependency, a `request_num_threads` API and a pool builder were all wired up and all unused.
Nothing in 418 differential comparisons, 276 Rust tests or 200 matrix configurations could see this:
every one of them measures *whether the answer is right*, and the answer was right. Only a runtime
measurement can see that the right answer is being computed one core at a time.

Parallelising it is safe by construction and the split is forced by the RNG: R draws all `nboot`
permutations from **one** MT19937 stream, so the draws stay sequential and only the aggregation goes
parallel, where each replicate reads a disjoint permutation and writes a disjoint slot and indexed
`collect()` preserves replicate order -- leaving the `F80` accumulation order within a replicate
untouched. Verified: 272 -> 276 Rust tests, shim gate still 418/418, matrix still 200/200, `R CMD
check` still OK, and `identical()` on both real fixtures.

**Measured, on the authors' own data, at the default `nboot = 100`:**

| fixture | nC | R upstream | Rust, 8 threads | speedup | peak RSS | parity |
|---|---:|---:|---:|---:|---:|---|
| human skin | 7 563 | 93.01 s [92.94, 93.80] | 2.19 s [2.08, 2.24] | **42.45×** [41.49, 45.12] | 1 251 MB | identical |
| mouse wound | 21 557 | 217.23 s [215.58, 222.00] | 8.35 s [8.25, 8.51] | **26.01×** [25.33, 26.89] | 2 424 MB | identical |

All **1 207 952** `Prob` values across the two fixtures are bit-identical to upstream: max absolute
difference 0, and zero deviating points in either log-log scatter. That is the acceptance gate
applied to the authors' data rather than to a fixture.

Strong scaling on human skin: 10.58 s / 5.72 s / 3.52 s / 2.16 s at 1 / 2 / 4 / 8 threads, i.e.
**4.91× at 61 % efficiency** on 8 physical cores, and 2.34 s at 16 -- *slower* than 8, because the
extra eight are SMT siblings. That is the concrete reason `parallel::detectCores()` cannot be used
on this host: it reports 16, and a table that counted those as cores would claim an efficiency the
hardware does not have.

**Two harness defects worth keeping as rules.** The runner first defaulted `CELLCHATRS_DB` to the
human export, so the mouse fixture ran against the human database: upstream read the mouse database
and succeeded, the kernel read the human one, found that the human subunits of `TGFbR1_R2` are not
in the mouse expression matrix, and returned `subscript out of bounds` -- the *correct* answer, and
extendr turns a returned `Err` into that R condition. It read as a kernel bug and cost more than the
benchmark did, so the runner now checks `MANIFEST.tsv` for both species and the pinned SHA and
refuses to start on a mismatch. And the sweep recorded `CELLCHATRS_THREADS`, which is only a
*request*; printing rayon's actual pool exposed rows labelled 1 through 16 that had all run at 8.

The projection is now retired. It assumed the streaming was already parallel when it was not, and
underestimated per-replicate cost by 2.6×, so it was a plausible model of a sequential kernel written
while the kernel was sequential -- internally consistent and wrong, which is exactly the failure the
objective's "do not publish a multiplier before it is measured" is about.

### Phase 5 — property tests, and the accumulator bug they found

`proptest` has been a declared dev-dependency since the crate was created and **no test used it**.
The objective lists "property tests via proptest" as a required layer, and the randomised
metamorphic and fuzz layers were also largely absent. `src/rust/crates/r-core/tests/properties.rs` now has 18
properties over the primitives everything else is built on: `order_f64` / `order_f64_multi` against a
transcription of R's documented rule, `r_max`'s `NaN` propagation, `r_mean`'s bounds, `nnzero`'s
`NA`-on-any-missing rule, all five `GroupMean` variants on constants and bounds, `F80` summation
and special-value propagation and commutativity, `sample.int` always being a permutation, MT19937
reproducibility, and `aggregate_1`'s level ordering.

It found a real bug immediately. `F80::mul` opened with

```rust
if self.m == 0 || other.m == 0 { return F80::ZERO; }
```

and **infinity is encoded with a zero mantissa**, so that line caught every infinite operand:
`prod(c(-Inf))` returned `0` where R returns `-Inf`, and `(-Inf) * 2` was `+0`. `r_prod` is reachable
only through the cofactor path in `computeExpr_coreceptor`, no fixture puts an infinity in a
cofactor column, and every parity test pins a finite case -- so 276 unit tests, 418 differential
comparisons and 200 matrix configurations all missed it.

The fix is R's table, which is **not** IEEE: an infinite product is `Inf` with the XOR of the signs,
including `Inf * -Inf` = `-Inf` (IEEE says `NaN`), and `NaN` arises only from zero times infinity.
Two related consequences fell out: a zero product must keep the XOR sign (`prod(c(5, -0))` is `-0`,
so multiplication is commutative but a signed zero is not), and `to_f64` was returning a bare `0.0`
for a zero mantissa, dropping every signed zero on the way out.

Six of the eighteen properties failed on their first run, and in every case the *property* was the
thing that was wrong, not the kernel: a column-major buffer indexed as row-major (the third distinct
layout mistake in this port, after the slice order and the `dim` vector), a bit comparison that flags
correct signed-zero arithmetic, an F80 tolerance bounded by the sum rather than by the magnitude of
the cancelling partials, a `same_group` range assumed to start at output position 0, and two attempts
to use a naive `f64` product as the oracle for a LONG_DOUBLE product -- which cannot work, since
`prod(c(8.98e-288, 5.77e-49, Inf))` is `Inf` in R and `NaN` in a `f64` fold. Each is now a comment
explaining the trap, because a property that fails for a reason nobody can reconstruct is worse than
no property.

Gates after: `cargo test --release` **298 tests across 23 binaries**; shim gate **418/418**; matrix
**200/200**, 1051/1051 pairs; real-data parity still exact on both fixtures (1 207 952 `Prob` values,
max |diff| 0); `parity.json` **141/145**; `R CMD check` `Status: OK`.

Still absent from the objective's required test layers, and the honest next list: the metamorphic
transforms that are not yet a layer (cell permutation, cluster relabelling, gene permutation, cell
duplication, nboot invariance of `Prob` -- only constant scaling is covered, via the matrix's `scale`
axis); randomised fuzzing of L-R *database* structures; statistical equivalence against an
independent RNG stream; end-to-end tutorial reproduction; and CI across stable/nightly Rust and
multiple R versions.

### Phase 6 — metamorphic invariance and L-R database fuzzing

Two more of the objective's named layers, both absent until now.

**`tests/parity/metamorphic.R`** applies each of the six named transformations and asserts two
separate things per transformation: the port still matches upstream *on the transformed object*, and
the relation holds relative to the untransformed run. 40 checks, all passing. Three of the six
relations are not what they look like and the obvious assertion is false:

* **Cell permutation** is not an invariance of `Pval` -- the bootstrap draws `sample.int(nC, nC)` over
  cell indices, so permuting the input moves `Pval` (measured: it does) while the observed `Prob`
  stays put. The transform also has to carry the labels with the columns, or it is not a permutation
  at all but a reassignment of expressions to different groups.
* **Cluster relabelling** permutes the rows *and* the columns, since `prob[i, j, l]` is "group `i`
  sends to group `j`". The relation is `prob[σ(i), σ(j), l] == prob_orig[i, j, l]`.
* **Constant scaling** is invariant only up to rounding unless the factor is a power of two:
  `max(f*x)` is not `f*max(x)` in floating point. Measured residue 0 at `f = 0.5` and `1e4`, `1.1e-16`
  at `f = 3`, `1.7e-21` at `f = 7.25`. The claim that *is* exact for every factor is
  `identical(port, upstream)` on the scaled object.

Each block carries a negative control, because an invariance that holds because the transform did
nothing is not a test. The objective's four recorded invariants -- dimnames consistency,
`Pval` in `{k/nboot}`, `Pval[Prob == 0] == 1`, `Prob` in `[0, 1]` -- are checked on every object the
file builds.

**`src/rust/crates/r-core/tests/db_fuzz.rs`** generates hostile L-R databases (empty subunit lists, names
colliding across the complex/cofactor/symbol tables, self-referential and overlapping complexes) and
asserts that `resolve_entity` and `compute_expr_lr` are total: a value or an `ExprError`, never a
panic, which is what makes upstream's `subscript out of bounds` an `Err` rather than a crash. The
resolution is compared against a `BTreeMap` reference written independently in the test.

The kernel passed all four properties. The properties were wrong four times, and each time the
generator had found a rule I had mis-stated: a complex with an all-empty subunit list contributes
nothing at all; the empty name is dropped by `gene[gene != ""]`; names outside `geneInfo$Symbol` are
dropped by `checkGeneSymbol`; and the split is on the **official symbol list**, not on membership of
the complex table. In the real database the complex names are not symbols and the plain gene names
are, which is why that is invisible by hand.

Gates after: `cargo test --release` **302 tests across 24 binaries** (22 new property/fuzz tests, and
the six shrunken regression cases proptest found are committed so they re-run every time); shim gate
**418/418**; metamorphic **40/40**; matrix **200/200**; real-data parity exact on both fixtures;
`parity.json` **150/154**; `R CMD check` `Status: OK`.

Still absent from the objective's required test layers: statistical equivalence against an
independent RNG stream, end-to-end tutorial reproduction, and CI across stable/nightly Rust and
multiple R versions.
