# cellchat-rs

A Rust re-implementation of the inference kernel of
[CellChat](https://github.com/jinworks/CellChat) (Jin et al., *Nature Communications* 2021;
*Nature Protocols* 2024), exposed to R as a **drop-in, output-identical accelerator**.

**Read [`PLAN.md`](PLAN.md) first** — it contains the source audit, the measured R
baseline, the performance model, the parity contract, the phased work plan, the benchmark
methodology, and the locked decisions in §14. This README is only the map.

## Locked decisions

| | |
|---|---|
| Scope | numeric core + R shim; `visualization.R` / `app.R` stay in R |
| Parity | **bit-identical** — `identical()` on the whole S4 object is the gate |
| Missing complex subunit | replicate upstream's `subscript out of bounds` (and report it upstream) |
| Spatial KNN | exact k-d tree, replacing approximate Annoy; divergence published |
| Binding | `extendr` |
| Upstream | `main` @ `75253cd0` (v2.2.0.9001) — see [`UPSTREAM.md`](UPSTREAM.md) |
| Distribution | GitHub first, then CRAN after maintainer sign-off |
| Benchmarks | this host (8c/16t, 32 GB), pinned with `taskset` |

## Why

`CellChat::computeCommunProb()` (`R/modeling.R:63-327`) dominates CellChat's runtime.
Measured on this machine against upstream R, at CellChat's default `nboot = 100`
(taskset-pinned, medians over repeats with bootstrap CIs; see `bench-runner/results/`):

| fixture | cells | genes | L-R pairs | upstream | Rust | speedup |
|---|---:|---:|---:|---:|---:|---:|
| human skin (Figshare 24470719) | 7 563 | 1 094 | 1 583 | 101.4 s | 2.41 s | 42.1x |
| mouse wound (Figshare 21896400) | 21 557 | 1 101 | 1 568 | 245.3 s | 9.41 s | 26.1x |

Full-pipeline numbers (what a user observes, I/O through centrality): 16.7x skin, 9.7x
wound. See `docs/BENCHMARKS.md` for the Amdahl decomposition and `paper/paper.md` for
the evidence-limited summary.

92–95 % of that time is the bootstrap aggregation
`aggregate(t(data), list(group), triMean)`, repeated `nboot` times — pure numerics with
independent work per gene per bootstrap, i.e. ideally parallel and vectorisable.

## Status

| | |
|---|---|
| Phase 0 — baseline, fixtures, harness | **done** (numbers above) |
| Phase 1 — `r-core` numerics + parity | **done** (159/159 parity quantities; see `parity.json`) |
| Phase 2 — R shim | **done** (drop-in overrides + 56 delegated names; tutorial identical end to end) |
| Phase 3 — standalone CLI | **done** (`crates/cellchatrs-cli`, 14 configs identical to upstream) |
| Phase 4 — benchmarks | **done** (kernel, Amdahl, spatial divergence, synthetic grid; see `docs/BENCHMARKS.md`) |
| Phase 5 — publication | in progress (`vignettes/`, `_pkgdown.yml`, `paper/paper.md` drafted; not yet released) |

The R↔Rust plumbing is **working and verified**: `R CMD INSTALL .` builds the Rust code
via `configure`, and `R CMD check` is clean.

## Layout

| path | what |
|---|---|
| `crates/r-core/` | the numerics; no R, no I/O, all pure functions (parity-testable) |
| `crates/cellchatrs/` | extendr `cdylib`; argument marshalling only |
| `R/` | the `cellchatrs` R package: `.onLoad`, and later the `computeCommunProb` override |
| `bench-runner/` | the Phase-0 R harness that produced the baselines |
| `tests/` | fixtures, parity, fuzz (populated in Phase 1) |
| `docs/` | `SEMANTICS.md` (parity contract), `BENCHMARKS.md` (results) |

## Build

```sh
cargo build --release          # workspace
cargo test                     # numerics
R CMD INSTALL .                # R package (runs cargo via ./configure)
```

```r
library(cellchatrs)
cellchatrs_threads()                        # live rayon pool size
Sys.setenv(CELLCHATRS_THREADS = "8")        # pin before load; benchmarks must pin
```

`parallel::detectCores()` is **not** trustworthy here (it reports 16 for
`logical = FALSE` on an 8c/16t host). Always pin with `CELLCHATRS_THREADS` or `taskset`
when benchmarking.

## Licence

GPL-3, matching CellChat.
