# CellChat 2.2.0

Not yet tagged; changes since `v0.1.0`.

- The R package is renamed from `cellchatrs` to `CellChat` and installs independently: it
  bundles the pinned upstream source and databases, no longer needs the original package,
  and is loaded in its place. The Rust crates and the CLI keep the `cellchatrs` name.
- The version follows the pinned upstream (v2.2.0.9001, with the fourth component dropped
  for CRAN) rather than this repository's own tag sequence.
- `CellChatDB.human`, `CellChatDB.mouse`, `CellChatDB.zebrafish`, `PPI.human` and
  `PPI.mouse` are exported and bound on load to the bundled upstream copies, read on first
  use. Bare names work in every install route, as they do under upstream's `LazyData`;
  previously a built install (including `remotes::install_github`) had no way to reach
  them except `system.file()`. `data(CellChatDB.human)` still warns in a built install.
- `parity.json` separates compatibility quantities from alternative-algorithm checks: 148
  quantities, 147 at a passing rung, plus 12 checks of the exact spatial path. The one
  untested quantity is `package.independent_dropin`, which needs the original installed
  beside this package.

# cellchatrs 0.1.0

First public release: a drop-in, output-identical accelerator for CellChat inference.

## What it is

- `computeCommunProb()` and the surrounding numerics reimplemented in Rust (`r-core`,
  zero R dependency), exposed through an R package that overrides the matching CellChat
  functions. Existing CellChat code runs unchanged.
- Verified by 160 machine-checkable parity quantities (`parity.json`): `identical()` on
  whole S4 objects across a 200-configuration matrix, an installed-shim differential gate,
  metamorphic and fuzz suites, an independent-RNG stream, a full tutorial reproduction,
  and a standalone CLI checked against upstream.
- Measured speedups at the tutorial default `nboot = 100` on the pinned host (8c/16t AMD
  Zen 2, taskset-pinned, five timed repeats after a discarded warm-up, medians with
  bootstrap CIs, each run recording the host's per-core busy fractions): kernel 42.1x
  (human skin, 7,563 cells) and 26.1x (mouse wound, 21,557 cells); full pipeline 16.7x and
  9.7x with stated Amdahl limits of 24.8x and 13.3x.

## Corrections the evidence forced

- Upstream's approximate spatial neighbour search (Annoy) carries a reproducible bias:
  42/64 inter-group distances differ (2.8% max), moving 19.6% of `Prob` entries (1.0% of
  peak) on the authors' visium data. The port uses an exact k-d tree.
- `Prob` scores are unbounded above (observed max 2025.1); they are scores, not
  probabilities, and the port reproduces them unclamped.
- Four centrality measures (`hub`, `authority`, `eigen`, `page_rank`) call igraph
  directly: they are iterative solvers whose output varies run to run. Degrees, strengths
  and betweenness are ported bit-exact.
- `identifyOverExpressedGenes(do.fast = TRUE)` delegates to presto (a normal
  approximation, not a faster Wilcoxon); both branches select the same features on
  fixtures.
- With the I/O floor measured in seconds, downsampling large datasets for runtime
  reasons is unnecessary.

## Scope limits (locked, see PLAN.md 14.1)

Visualization and the Shiny app stay upstream R code. Plotting functions are not
overridden; numeric slots feeding the plots are all verified.
