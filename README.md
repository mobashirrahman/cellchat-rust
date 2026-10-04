# cellchat-rust

A Rust re-implementation of the inference kernel of
[CellChat](https://github.com/jinworks/CellChat) (Jin et al., *Nature Communications* 2021;
*Nature Protocols* 2024), exposed to R as a **drop-in, output-identical accelerator**.
The R package it installs is named `cellchatrs`.

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
| Phase 1 — `r-core` numerics + parity | **done** (160/160 parity quantities; see `parity.json`) |
| Phase 2 — R shim | **done** (drop-in overrides + 56 delegated names; tutorial identical end to end) |
| Phase 3 — standalone CLI | **done** (`src/rust/crates/cellchatrs-cli`, 14 configs identical to upstream) |
| Phase 4 — benchmarks | **done** (kernel, Amdahl, spatial divergence, synthetic grid; see `docs/BENCHMARKS.md`) |
| Phase 5 — publication | in progress (`vignettes/`, `_pkgdown.yml`, `paper/paper.md` drafted; not yet released) |

The R↔Rust plumbing is **working and verified**: `R CMD INSTALL .` builds the Rust code
via `configure`, and `R CMD check` is clean.

## Layout

| path | what |
|---|---|
| `src/rust/crates/r-core/` | the numerics; no R, no I/O, all pure functions (parity-testable) |
| `src/rust/crates/cellchatrs/` | extendr `cdylib`; argument marshalling only |
| `src/rust/crates/cellchatrs-cli/` | standalone binary: the same numerics without R |
| `R/` | the `cellchatrs` R package: `.onLoad`, the `computeCommunProb` override, and 56 generated pass-throughs to upstream names |
| `bench-runner/` | the R harness that produced every published measurement |
| `tests/` | fixtures and golden corpora (`fixtures/`), differential and parity gates (`parity/`), LR-structure fuzzing (`fuzz/`) |
| `docs/` | `SEMANTICS.md` (parity contract), `BENCHMARKS.md` (results and protocol) |
| `paper/` | the evidence-limited manuscript draft |
| `parity.json` | machine-checkable ledger: 160 quantities, the rung each one reached |

## Verify the claims yourself

Nothing in this README has to be taken on trust: every number has a committed artifact and a
command that regenerates it.

**Prerequisites.** Building the R package currently requires Rust and Cargo 1.84 or newer on
Linux x86_64. The exact pinned CellChat source and database are bundled with the package. To
recreate the parity ledger, also check out upstream at its pinned commit and download the larger
Figshare fixtures; the helper downloads fixed file IDs and verifies their hashes:

```sh
git clone https://github.com/jinworks/CellChat ../CellChat
git -C ../CellChat checkout 75253cd0c9e68410e6e721a6d3a0419a1d7e358f
./scripts/fetch_ci_fixtures.sh parity
export CELLCHAT_SRC=../CellChat
```

**The parity ledger.** `parity.json` is generated, not written by hand:
`scripts/parity_report.py` *runs* the Rust suite and every R gate as subprocesses and records the
rung each quantity reached, so regenerating it re-derives the whole claim.

```sh
python3 scripts/parity_report.py        # writes parity.json; exits nonzero on a regression
```

**The individual gates.**

```sh
cargo test --workspace --release        # 331 unit / property / RNG-golden / x87-oracle tests
R CMD check --no-manual                 # packaging, docs, vignette
Rscript tests/parity/check_identical.R  # 427 whole-S4 identical() comparisons, 9 configurations
Rscript tests/parity/check_matrix.R     # 200 configurations, 1051 axis pairs
Rscript tests/parity/tutorial_repro.R   # end-to-end tutorial, 16 comparisons
Rscript tests/parity/check_cpu_gate.R   # the benchmark host-contention gate, 47 checks
python3 bench-runner/test_analyse.py    # the benchmark analysis step, 15 checks
```

**The benchmarks.** Every figure in `bench-runner/results/` came from a runner under the protocol
in `docs/BENCHMARKS.md` — warm-up discarded, five timed repeats, median with a bootstrap CI, CPUs
pinned, a memory-headroom gate, and a host-contention gate that refuses to measure a busy machine
and records the per-core load when it proceeds. Re-running any of them writes its own evidence
beside the numbers.

## Build

```sh
cargo build --release          # workspace
cargo test                     # numerics
R CMD INSTALL .                # R package (runs cargo via ./configure)
```

The R package build is offline: `src/rust/r-build/vendor.tar.xz` contains the pinned Rust source
dependencies and is SHA-256 checked before compilation. Regenerate it after changing
`src/rust/r-build/Cargo.lock` with `scripts/vendor_rust_deps.sh`.

```r
library(CellChat)
cellchatrs_threads()                        # live rayon pool size
Sys.setenv(CELLCHATRS_THREADS = "8")        # pin before load; benchmarks must pin
```

`parallel::detectCores()` is **not** trustworthy here (it reports 16 for
`logical = FALSE` on an 8c/16t host). Always pin with `CELLCHATRS_THREADS` or `taskset`
when benchmarking.

## Licence

GPL-3, matching CellChat.
