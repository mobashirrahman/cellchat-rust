# Semantics contract

The normative specification that `r-core` implements. `PLAN.md` §4 is the summary;
this file is the version CI checks against. Every row has a test.

## Parity target: **bit-identical** (locked — `PLAN.md` §14.2)

Every requirement below is **Exact**: `identical()` to upstream R, bit for bit, on the
whole S4 object. The `1 ulp` and `Tol(n)` rungs appear in `parity.json` only as a
contingency if a requirement turns out to be unreachable; they are not a target.

Two consequences that are easy to forget:

* **No FMA contraction and no reassociation** in the aggregation. Contraction would
  change rounding. Build with `-C target-feature=-fma` where the compiler would
  otherwise contract, and add a test asserting the aggregation is bit-identical to a
  `-C target-cpu=native` build.
* **R's `LONG_DOUBLE` accumulation must be emulated with a real 64-bit mantissa**
  (`u128` or `long double` FFI). An `f64` sum of 4 values is *not* sufficient (R10).

Upstream is pinned at `75253cd0c9e68410e6e721a6d3a0419a1d7e358f` (v2.2.0.9001) — see
`UPSTREAM.md`. "Upstream" means exactly that commit.

## Requirements

| # | Quantity | Rung | Upstream source | Notes |
|---|---|---|---|---|
| R1 | bootstrap permutations | Exact | `R/modeling.R:207` `sample.int(nC, nC)` | MT19937 + `R_unif_index`; high value, low cost |
| R2 | `triMean` | Exact | `R/modeling.R:878` | via `collapse::fquantile`, **not** R's `quantile` |
| R3 | `truncatedMean` | Exact | `R/modeling.R:70` -> `mean.default` | R's `trim` slicing, incl. the `trim >= 0.5` short-circuit to `median` |
| R4 | `median`, `thresholdedMean`, `geometricMean` | Exact | `R/modeling.R:863,889,491` | even-`n` `median` is `mean()` of the two middles; `log`/`exp` must match libm |
| R5 | Hill function, outer products, `P.spatial` | Exact | `R/modeling.R:222-255` | pure IEEE, same association order |
| R6 | `net$prob` | Exact | — | follows from R1–R5 |
| R7 | `net$pval` | Exact (integers) | `R/modeling.R:303` | counting; order-independent |
| R8 | `dimnames`, `options$parameter`, slot contents | Exact | `R/modeling.R:311-320` | `identical()` on the whole S4 object |
| R9 | NaN / NA handling in `fquantile` | Exact | — | NaN is missing when `na.rm = TRUE`; also the `n > 1e5` radix fast path (R11) must trigger on the same condition |
| R10 | `mean()` | Exact | `src/main/summary.c:476` (`real_mean`) | **Two-pass corrected mean in `LDOUBLE`**, see below |
| S1 | KNN in `computeRegionDistance` | **intentional divergence** | `R/modeling.R:1207` | Annoy is approximate; replaced by an exact k-d tree (`PLAN.md` §14.4). Divergence is measured on real spatial data in Phase 4 and published. |

### Why R2 says "collapse", not R

`triMean` is `mean(collapse::fquantile(x, probs = c(.25,.5,.5,.75), na.rm = TRUE))`.
R's own `quantile(type = 7)` uses a **1-based** `h = (n-1)*p + 1`; collapse uses a
**0-based** `h = (n-1)*Q` with the interpolation weight `h - floor(h)`. The `+ 1` is a
rounding opportunity collapse avoids.

Measured against R 4.3.3 + collapse 2.1.8 on 20 000 random vectors (ties, `NaN` present),
reproducibly (`tests/parity/gen_stats_golden.R`, recorded in `stats_golden.txt`):

| level | disagrees | max absolute difference |
|---|---:|---:|
| `fquantile` outputs (4 per vector) | **1.13 %** | 2.22e-16 (1 ulp) |
| final `triMean` values | **3.51 %** | 2.22e-16 (1 ulp) |

So porting R's formula would be *wrong* for 1 vector in 30. We port **collapse's**.
(`collapse_uses_radix_order` and `collapse_divergence_is_real` keep this honest: if a
future R/collapse update made the two agree, the tests fail rather than rotting.)

A second, sharper difference: collapse's live path is **`dquickselect`**
(`fnth_fmedian_fquantile.c:414`), not the `FQUANTILE_ORDVEC` macro in the same file.
`dquickselect` subtracts the integer part of `h` *before* testing `h <= 0.0`, so a zero
interpolation weight returns the order statistic without ever touching the next value.
The macro tests the un-reduced `h`, which would evaluate `0.0 * Inf = NaN`. Observable:
for `c(1, Inf, 2)` at `Q = 0.5` collapse returns `2.0`, reachable only through the
`dquickselect` guard. `fquantile_type7` implements the `dquickselect` semantics.

For `triMean` the `radixorder` fast path is unreachable: its condition is
`length(x) > 1e5 && length(probs) > log(length(x))`, and with `length(probs) == 4` the
second clause needs `length(x) < e^4 ≈ 54.6`. Pinned by a test.

### R10: R's `mean` is a two-pass corrected mean, not a sum

`src/main/summary.c:476`:

```c
LDOUBLE s = 0.0;
for (k) s += dx[k];              // pass 1, 80-bit
if (R_FINITE((double) s)) s /= n;
else { s = 0.; for (k) s += dx[k]/n; }   // overflow recovery
if (finite_s && R_FINITE((double) s)) {
    LDOUBLE t = 0.0;
    for (k) t += (dx[k] - s);    // pass 2, residuals
    s += t/n;
} else if (R_FINITE((double) s)) { /* second recovery form */ }
return ScalarReal((double) s);
```

Three details, each of which changes the last bit:

1. `R_FINITE` is applied to **`(double) s`** — after narrowing, not to the 80-bit value.
   A finite accumulator that overflows `f64` therefore takes the recovery branch.
2. The two recovery forms are *not* interchangeable: the first computes `dx[k] / n` in
   **`f64`** and then widens, the second `(dx[k] - s) / n` in 80-bit.
3. `f64` accumulation is genuinely not equivalent. `mean(c(1e16, 1, -1e16))` is
   `0.33365885416666669` in R and `0` with a naive `f64` chain — and it is *not* `1/3`
   either, because pass 1 lands on the `1` only in 64 bits.

`crates/r-core/src/longdouble.rs` implements the 80-bit type and is differentially
tested against real x87 `long double` (`tests/parity/gen_f80_ref.c`, 6 180 cases
including 180 chains up to 50 000 elements).

### R9b: `triMean`'s inner `mean` uses `na.rm = FALSE`

`triMean <- function(x, na.rm = TRUE) mean(collapse::fquantile(x, probs = ..., na.rm = na.rm))`
— the `na.rm` belongs to `fquantile`, **not** to `mean`. So a `NaN` among the four
quantiles propagates. Observable for `c(1, -Inf, 2)`, whose quartiles are
`(NaN, 1, 1, 1.5)`: upstream's trimean is `NaN`, while a `na.rm = TRUE` mean would
return `1.1666...`. `geometricMean` and `thresholdedMean` do pass `na.rm` to `mean`.

### R-isms in the database layer

Four that all produced *plausible, correctly-sized, wrong-ordered* output, so none of
them would show up in a length or set comparison:

1. **`extractGeneSubset`'s complex test is `x %in% symbols == "FALSE"`.** That works only
   because R coerces the logical to character before the `==`, so `as.character(FALSE) ==
   "FALSE"`. The branch therefore selects names **absent** from `geneInfo$Symbol` — the
   complex names. The shipped DB has zero literal `"FALSE"` symbols, so reading the test
   literally drops every complex and every subunit.
2. **`order()` on the annotation factor is level order; on the character vector it is
   alphabetical.** `subsetData` orders *before* converting back to character. Ordering
   after conversion silently produces a different table with a different `nLR1`.
3. **`unlist()` on a multi-row data frame is column-major.** The subunits reach the
   kernel through `unlist(complex_input[match(complex, rownames(...)), ])`, so the order
   is: every row's `subunit_1`, then every row's `subunit_2`, .... The empty cells are
   load-bearing because they carry the column positions — the human DB really does contain
   `IL12AB` = `IL12A, <empty>, IL12B`, and dropping the empty would move `IL12B` into the
   `subunit_2` slot and shift every later subunit.
4. **`c(agonist, antagonist, co_A_receptor, co_I_receptor)` is column-wise.** R appends
   every interaction's agonist, then every antagonist, then every co-activator, then every
   co-inhibitor. Reading the four fields per row interleaves them and reorders the
   cofactor subunits that follow.

Plus one that the tests caught in the *export* rather than the port: `stopifnot(cond,
"message")` evaluates the string as a further condition, so the message itself becomes the
failure. Use `if (!cond) stop(...)`.

### R1 findings: two things that are not what they look like

Both were found by differential testing, not by reading, and both would have produced a
stream that is *uniform and plausible* while being completely wrong — the worst kind of
bug for a parity port.

1. **R's `Int32` is `unsigned int`, not `int`.** `RNG.c` states it in a comment
   (`/* typedef unsigned int Int32; in Random.h */`). The right shifts in the MT19937
   regeneration and tempering are therefore *logical*, and R's stream is the standard
   MT19937. Porting it with signed `i32` shifts yields a perfectly uniform but
   entirely wrong stream. Settled empirically: all 2 × 10⁶ draws of
   `set.seed(1); runif(2e6)` lie exactly on the `k · 2⁻³²` grid and none equals
   `fixup`'s fallback value — both impossible with signed shifts.
2. **`FixupSeeds` clobbers `mti` with a constant.** `RNG_Init` fills 625 words with an
   LCG, but `FixupSeeds(MERSENNE_TWISTER, initial = TRUE)` then sets `i_seed[0]`
   (which *is* `mti`) to the literal `624`. A port that keeps the LCG value reads
   `mt[negative]` and yields a degenerate "permutation" that still contains every
   value. Caught by the property test `sample_int_is_a_permutation`, not by any
   golden value.

A third, smaller trap: `sample.int`'s without-replacement loop decrements the
population as it goes (`x[j] = x[--n]`), and the *shrinking* `n` is what feeds
`R_unif_index`. Sampling against the original `n` produces duplicates.

## R-isms that must be reproduced verbatim

1. **`P.spatial` is mutated inside the loop.** `R/modeling.R:232`:
   `if (i > nLR1) P.spatial <- P.spatial * adj.contact`.
   It fires on the first iteration past `nLR1`, and iteration `nLR1+1` tests the
   *pre*-multiplied matrix in the `sum(P1_Pspatial) == 0` early exit. Idempotent in
   practice, but the iteration index is observable.
2. **Missing cofactor rows inject `NA`.** `cofactor_input[missing, ]` is an all-`NA` row;
   `unlist()` then produces `NA`, and `v[NA]` inserts an `NA` element. `intersect()`
   removes it later. Benign, but the path must exist.
3. **An unresolvable coreceptor means "no modulation", not zero.**
   `computeExpr_coreceptor` returns `matrix(1)` when `intersect(...)` is empty.
4. **A single zero kills a geometric mean.** `exp(mean(log(x), na.rm = TRUE))` is `0`
   whenever any subunit is `0`, because `log 0 = -Inf`.
5. **`max()` spans the whole matrix, zeros included**, and warns (`-Inf`) when
   everything is non-finite.
6. **`Pval[Prob == 0] <- 1` runs after the loop**, overriding per-pair values.
7. **`Pboot - Pnull > 0` is strict.** Exact ties are not rejections.
8. **`aggregate()` emits groups in factor order**, and `levels(object@idents)` defines
   the output `dimnames`. `nlevels != length(unique(group))` must raise upstream's error.
9. **`type` is `match.arg`-matched**, so `"tri"`, `"med"`, … are valid.
10. **`table(groupboot)/nC` is recomputed per `(i, b)`** in upstream when
    `population.size = TRUE`; we hoist it and must get identical values.
11. **`computeExpr_complex` has no bounds check** (`R/modeling.R:537`):
    `data.use[RsubunitsV, ]` uses raw subunit names, unlike the coreceptor path which
    intersects first. One absent subunit ⇒ upstream aborts with
    `subscript out of bounds`. Confirmed empirically: 8 of 1614 human-skin L-R pairs had
    to be dropped to get the fixture to run at all.
    **Decision (locked):** replicate the error with the identical message, and file it
    upstream as a bug — see `PLAN.md` §14.3. `r-core` therefore needs a `MissingSubunit`
    error variant whose R-side rendering is byte-compared against upstream's message.

### R-isms found in the `computeExpr_*` family

Found by `crates/r-core/src/expr.rs` + `tests/expr_parity.rs` against
`tests/fixtures/expr_golden.txt`, all of which change a returned bit and none of which
are visible from the source alone.

12. **A plain gene is copied verbatim; a complex is geometric-averaged.**
    `computeExpr_LR` branches on `which(geneLR %in% rownames(data.use))` and assigns
    `dataLavg[index.singleL, ] <- dataL1avg` — no `geometricMean`. Only the complex path
    takes `exp(mean(log(x), na.rm = TRUE))`. These are not the same function:
    `exp(log(3)) == 3.0000000000000004`. Folding single genes through the log/exp round
    trip perturbs **every** plain-gene ligand and receptor in the network by ~1 ulp, so
    `r-core` has a distinct `EntityRows::Single` variant rather than one code path.
13. **The antagonist Hill term flips the numerator.** Agonist is
    `data.avg^n / (Kh^n + data.avg^n)`; antagonist is `Kh^n / (Kh^n + data.avg^n)`. The
    two agree only at `x == Kh`, so sharing one helper is silently wrong for every
    cofactor whose expression is not exactly `Kh`. The invariant that does hold is
    `h_agonist + h_antagonist == 1`, hence `ag + an == 2` for a single subunit.
14. **A cells × genes matrix has the *gene* index striding.** R stores
    `m[cell, gene]` at `gene * n_rows + cell`. The genes × groups matrices elsewhere in
    the port stride the other way, so the two layouts coexist and both are correct for
    their own shape. Getting `CellExpr` wrong reads a transposed matrix: every group mean
    silently becomes the mean of two unrelated cells.
15. **`aggregate()` returns group-major, so `t()` is a real transpose.** R's result is
    `n_groups` rows × `n_cols` columns column-major (index `g * n_cols + j`); the
    `computeExprGroup_*` paths then apply `t(data.avg[, -1])` to get subunits × groups
    (index `j * n_groups + g`). An implementation that picks the wrong stride returns the
    right *multiset* of numbers in the wrong order, which is why it survives a
    multiset-only check.
16. **`unlist()` on a `data.frame` column slice is column-major, not row-major.** The
    fixture's gene universe is `unique(c(unlist(complex_input[, subunit_cols]),
    unlist(cofactor_input[, cofactor_cols])))`, i.e. every `subunit_1` down all rows, then
    every `subunit_2`. Enumerating row-by-row keeps the gene *set* identical and permutes
    the *values*, so every name still resolves and every number is still plausible.
17. **`aggregate()` drops empty levels rather than emitting `NA`.**
    `aggregate.data.frame` ends with `y <- y[match(sort(unique(grp)), grp, 0L), ]`, so a
    level with no cells gets no row at all. `r-core` returns a dense `n_groups` buffer
    with `NaN` instead. **Documented divergence**, unreachable from CellChat because
    `ident` is derived from the cell labels, and pinned by
    `empty_groups_are_dense_nan_and_that_is_a_documented_divergence`.

### R-isms found in `computeCommunProb`

Found by `crates/r-core/src/prob.rs` + `tests/prob_parity.rs` against
`tests/fixtures/prob_golden.txt` (eight parameter configurations, `Prob` and `Pval`
bit-exact for all of them). 18–21 were each worth at least one wrong-bit bug in an
implementation that looked correct.

18. **`P2` and `P3` are rank-1 outer products, and mixing their indices is invisible
    almost everywhere.** `computeExpr_agonist`/`_antagonist` return a 1 x k row and the
    kernel wraps it in `Matrix::crossprod`, so `P2[a,b] = p2[a] * p2[b]` and
    `P3[a,b] = p3[a] * p3[b]`. Writing `p2[a] * p3[b]` agrees whenever `p2 == p3` or
    either is 1 — i.e. for every pair that does not have *both* an agonist and an
    antagonist. Separately, the bootstrap loop naturally wants one index hoisted out of
    the inner loop, and squaring that hoisted value is right only on the diagonal.
19. **`P1*P2*P3*P4*P.spatial` associates strictly left to right**, and every intermediate
    rounding is a rounding. The multiplication order is reproduced literally rather than
    regrouped (e.g. `((p1 * p2ab) * p3ab) * p4` * `psp`).
20. **`P.spatial` is mutated once per iteration past `nLR1`,** so the value used by `Pnull`
    at iteration `i > nLR1` is `P0 * adj.contact` and not `P0`. This is only idempotent
    because `adj.contact` is 0/1. It is hoisted to `p_scaled` in [`crate::prob::Kernel`],
    and `P1_Pspatial` in the `sum(...) == 0` early exit sees the *already mutated* matrix
    too, so both branches agree on `P_eff(i)`.
21. **`Prob` can legitimately exceed 1.** The co-agonist term is a *product* of `1 + h`
    over all cofactor subunits, each in `[1, 2]`, so a 10-subunit agonist reaches `2^10` and
    the rank-1 outer product squares that. Not a bug and not something to "fix": a
    normalisation here would change every published CellChat result. Pinned by the corpus.

Also confirmed by the corpus rather than assumed:

* `data.signaling` is **genes x cells**; `aggregate(t(data.use), list(group), ...)`
  transposes. Building it the other way round fails inside `aggregate.data.frame` with
  `arguments must have same length`.
* `sample.int` returns **1-based** indices.
* `Prob[,,i]` and `Pval[,,i]` are `n_groups x n_groups` **column-major** slices, so
  `unravel(i, k)` gives the interaction first and the group pair last.
* `dimnames(Prob)` is `list(levels(group), levels(group), rownames(pairLRsig))` and
  `dimnames(Pval) == dimnames(Prob)`.

### R-isms found in `aggregateNet` / `subsetCommunication`

Found by `crates/r-core/src/net.rs` + `tests/net_parity.rs` against
`tests/fixtures/net_golden.txt` (four fixtures, three thresholds, three filter
combinations).

22. **`net$weight` is an `LDOUBLE` sum *across interactions*, not across cells.**
    `apply(prob, c(1,2), sum)` calls `FUN` on `prob[a, b, ]` — a length-`nLR` vector — so
    the accumulation is interaction-major within a single output cell, in 80-bit. An
    accumulator that narrows to `f64` after each add is an `f64` sum in disguise and is
    1 ulp off at 8 terms (2.27e-13 on 1230); at 3 or 4 terms it agrees, so the short
    fixtures pass and the long one does not.
23. **`reshape2::melt` on a 3-D array is not in the array's own order.** It builds the label
    frame with `expand.grid(dimnames(x))`, and `expand.grid` varies its **first** argument
    fastest, so the melted rows run interaction-major, then target, then source — the
    opposite of R's column-major array traversal. `as.vector(data)` is column-major and the
    labels are permuted to match, so the two agree on values and disagree on order. That
    order is the return value.
24. **`intersect(want, colnames(net))` sees column *presence*, not non-`NA` values.** An
    all-`NA` `evidence` column is in the output; a table with no `evidence` column at all is
    not. This is also why upstream's
    `net[rowSums(is.na(net)) != ncol(net), ]` guard drops a row only when it is `NA`
    *everywhere*, and why the output column set must not be derived from the surviving rows.
25. **`subsetCommunication_internal` raises rather than returning an empty table**
    (`stop("No significant signaling interactions are inferred based on the input!")`) when
    the *threshold* removes everything — but the later `sources.use` / `targets.use` filters
    can empty it afterwards, which only warns. The two are different shapes and the Rust
    binding distinguishes them.

### R-isms found in `computeAveExpr` / `subsetDB` / `subsetData`

Found by `crates/r-core/src/de.rs` + `tests/de_parity.rs` against
`tests/fixtures/de_golden.txt`. 26–28.

26. **`intersect` returns its *first* argument's order, not the matrix's.**
    `features.use <- intersect(features, rownames(data.use))` means
    `computeAveExpr(features = c("G10","G1","G10","GZZ","G2"))` produces rows
    `G10, G1, G2`: the caller's order, the duplicate collapsed, the absent name dropped.
    Sorting, or using the matrix's order, gives a matrix with the right *set* of rows in
    the wrong order — invisible to any value comparison that lines rows up by name.
27. **`subsetDB` flips `non_protein` to `TRUE` when `search` names
    `"Non-protein Signaling"`.** The flip happens *before* the
    `annotation != "Non-protein Signaling"` exclusion, so passing all four annotations
    explicitly keeps all 3233 rows even with `non_protein = FALSE`. Missing the flip is a
    994-row difference on the pinned human DB.
28. **`match.arg`'s error message uses U+201C/U+201D**, because R formats it with
    `sQuote()`. The message is compared byte for byte, so straight quotes are a mismatch.
    Same for `subsetDB`'s unknown-key error, which has no `(got "...")` suffix.

### A structural note about the exported database

`tests/parity/export_db.R` writes `interactions.tsv` **already in annotation order**, so
`r-core` never reproduces `order()` on the annotation factor, and `subsetData`'s reordering
is a no-op on the Rust side. The pinned `.rda` is *not* in that order, so the reordering is
a real change there. The invariant that makes this safe — and that
`de_parity.rs::the_corpus_has_records_for_everything_the_test_reads` checks — is:

> `r-core`'s interaction order equals upstream's order **after** `subsetData`.

`de_parity.rs` also exercises `order_by_annotation` directly on a reversed input, so the
sort is tested even though the shipped export does not need it.

## `identifyOverExpressedGenes` (`do.fast = FALSE`)

Nine R-isms, each of which produced a wrong answer before it was found. Every one is
pinned by `crates/r-core/tests/wilcox_parity.rs` against a corpus generated from the
pinned upstream (`tests/parity/gen_wilcox_golden.R`) and re-checked end to end through the
installed shim by `tests/parity/check_identical.R`.

1. **The Wilcoxon branch and the presto branch have different filters.** The presto
   branch (`do.fast = TRUE`) ends in
   `dplyr::filter(genes.de, pvalues < thresh.p, logFC_abs >= thresh.fc, pct.max > thresh.pc*100)`.
   The Wilcoxon branch has **no** `pct.max` filter at all; it selects with
   `alpha.min > thresh.pc` (a *fraction*) and `FC > thresh.fc`, then filters only on
   `pvalues < thresh.p` and, when `only.pos`, `logFC > 0`. Porting the presto filter
   wholesale into the Wilcoxon path silently drops features for any `thresh.pc > 0`, and
   `identifyOverExpressedLigandsReceptors` is routinely called with `thresh.pc = 0.01`.
2. **The marker table is ordered by p-value, not by feature.** Per group,
   `gde <- gde[order(gde$pvalues, -gde$logFC), ]`, then `rbind`. So the rows of
   `features.info` — and hence `var.features[["features"]]`, and hence the LR table — are
   in `(cluster, p-value ascending, logFC descending)` order. Emitting features in
   ascending order agrees on the *set* and disagrees on the *order*.
3. **`features.info` row names are the feature names**, for two unrelated reasons: in the
   Wilcoxon branch the unnamed `data.alpha[features, , drop = FALSE]` matrix argument
   carries `rownames == features` and `data.frame()` adopts them; in the `do.DE = FALSE`
   branch `rowSums()` on a `Matrix` returns a *named* vector and `data.frame()` adopts the
   first named argument's names. (With a dense `data.use` the second one gives `1..nrow`,
   so this row naming is sparse-vs-dense dependent upstream.)
4. **The schema collapses to `0x1` when nothing passes.** `markers.all <- data.frame()` is
   only ever `rbind`-ed when a group has a row clearing `thresh.p`; if none does, the frame
   keeps zero columns, and the following `markers.all$features <- as.character(...)` *adds*
   a single zero-length column. A "marker table" of `0x1` with just `features` is a real,
   reachable object (`thresh.p = 0.001` on ordinary data) and `identical()` sees the
   difference. The kernel returns `n_before_only_pos` so the shim knows which shape to build.
5. **`p.adjust(..., n = nrow(X))`, not `nrow(data.use)`.** `X` is the full
   `object@data.signaling`, so passing a `features` argument shrinks the adjusted p-values
   by the subset ratio and changes which features clear `thresh.p`.
6. **`min.cells` is only used by the `do.DE = FALSE` branch** (and by presto). The
   Wilcoxon branch ignores it, despite it being in the same signature — so a high
   `min.cells` is a no-op there, and the corpus pins that (`mincells` spec, `min.cells = 1000`).
7. **C's `sign(0)` is `0`, Rust's `f64::signum(0.0)` is `1.0`.** `z <- (z - sign(z) * 0.5)/SIGMA`
   therefore gets a *spurious* `-0.5` numerator at `z == 0`. For an all-ties gene `SIGMA == 0`,
   so `z` must be `NaN` (R returns `NaN`) and the p-value must be `NaN`; with the spurious
   numerator `z` becomes `-inf` and the p-value `0` — i.e. **every all-ties gene reported as
   maximally significant**. `r_core::wilcox::r_sign` is C's `sign`, and
   `an_all_ties_gene_gives_a_nan_p_value_not_zero` is the regression test.
8. **`rank(..., na.last = "keep")` ranks NA last**, in input order:
   `rank(c(1, NA, 2))` is `c(1, 3, 2)`, not `c(1, NA, 2)`. (`wilcox.test` drops NAs before
   ranking, so this only shows up through `rank` itself.)
9. **`wilcox.test`'s `correct` argument is load-bearing** and not visible in the returned
   object: for `W = 8` with `n.x * n.y` even, `correct = TRUE` gives `0.3976...` and
   `correct = FALSE` gives `0.3412...`. Reading a `correct = FALSE` golden record and
   comparing it against a `correct = TRUE` implementation passes by coin-flip.

### Corpus-generator hazards found here

* **`out$correct` does not exist.** `names(wilcox.test(x, y))` is exactly
  `statistic parameter p.value null.value alternative method data.name`. Asking for
  `$correct` yields `NULL`, `as.integer(NULL)` is `integer(0)`, and `sprintf("%s", integer(0))`
  returns `character(0)` — so the whole record collapses to a **blank line** rather than to
  an error, and the corpus silently contains 6 error records and 72 empty ones.
* **`rep(x, lengths = ...)` is not a length specifier.** `rep(c("g1","g2","g3"), lengths = c(30,30,20))`
  returns a length-3 vector, not 80. (`rep(x, times = c(30,30,20))` is right.) This produced
  a 3-cell fixture that every downstream comparison then "failed" on.

## `computeCommunProbPathway`

Four R-isms, each pinned by `crates/r-core/tests/pathway_parity.rs` (9 fixtures) and
re-checked end to end through the installed shim by `tests/parity/check_identical.R`.

1. **Both `sum`s accumulate in LONG_DOUBLE**, and it is observable: the `bigmag` fixture's
   upstream per-L-R totals are `30000000000000004`, where an `f64` loop gives
   `30000000000000000`. Since `pathways.sig` is then *ordered by those totals*, a 1-ulp
   difference can permute the whole `netP$pathways` list.
2. **The two `apply` calls sum in different orders.** `apply(prob, 3, sum)` visits each L-R
   slice in `(c, r)` order. `apply(prob, c(1, 2), by, group, sum)` splits by *column* (i.e.
   by L-R) and each column's `k*k` values are then summed in the sub-matrix's own
   column-major order, which is `(r, c)`. Using one order for both is the obvious mistake,
   and it moves the pathway totals.
3. **The pathway axis is *last*.** `apply(...)` with results on the last axis produces a
   `k x k x nPathways` array whose `aperm(..., c(2, 3, 1))` moves the pathway from first to
   last — the layout `object@netP$prob` actually has. Reading `prob.pathways` as "a stack of
   pathway matrices" and laying it out pathway-major gives the transpose, with entirely
   plausible numbers.
4. **`sort(..., decreasing = TRUE)` is `sort.int`, whose default method for a plain double
   vector is `"shell"` and is not stable.** Two fixtures (`tied`, `tied3`) have bit-identical
   pathway totals and show that R breaks the tie in ascending index order, so the port uses a
   stable descending sort. That is a *measured* claim, not an assumption: if the pinned
   upstream's shell sort ever permuted a tie, `sort_desc_r_index` and its tests are where the
   divergence would show up. Ties need two pathways whose L-R sets sum to bit-identical
   doubles, which does not occur on real data.

### One upstream failure is preserved rather than fixed

```r
prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))
```

`apply(prob, c(1, 2), by, group, sum)` returns a **2-d matrix** when `group` has one level and
a **3-d array** otherwise, so a one-pathway `LRsig` — every interaction in the same pathway,
which is ordinary — makes `aperm` raise `'perm' is of wrong length 3 (!= 2)`. The shim routes
that case to upstream; the Rust core computes a value for it and deliberately does not claim
it. Pinned by `single_pathway_errors_like_upstream`.

## `aggregateNet`'s filtered branch is broken upstream, and the port reproduces that

`aggregateNet` has two branches. The unfiltered one (`is.null(sources.use) & ...`) is a plain
`apply` over the probability array and is already bit-exact. The **filtered** one is:

```r
df.net$source_target <- paste(df.net$source, df.net$target, sep = "|")
df.net2 <- df.net %>% group_by(source_target) %>% summarize(count = n(), .groups = "drop")
a <- stringr::str_split(df.net2$source_target, "|", simplify = TRUE)
df.net2$source <- factor(df.net2$source, levels = cells.level[cells.level %in% unique(df.net2$source)])
df.net2$target <- factor(df.net2$target, levels = ...)
count <- tapply(df.net2[["count"]], list(df.net2$source, df.net2$target)), sum)
```

`stringr::str_split(x, "|", simplify = TRUE)` takes a **regular expression**, and `|` is the
alternation metacharacter, so it matches the empty string at every position rather than the
literal pipe:

```r
stringr::str_split(c("g1|g2", "g10|g3"), "|", simplify = TRUE)
#      "" "g" "1" "|" "g" "2" "" ""
#      "" "g" "1" "0" "|" "g" "3" "" ""
```

The matrix has `nchar(key) + 1` columns -- `""`, one per **character**, `""` -- so `a[, 1]`
is the empty string for every row and `a[, 2]` is a single character, not a group name. No
cell group is named `""`, and consequently:

* `remove.isolate = TRUE` restricts the levels to `character(0)`, and `tapply` returns a
  **`0 x 0`** matrix;
* `remove.isolate = FALSE` keeps `cells.level`, every entry is `NA`, and the two `is.na`
  fills turn it into a **`k x k` matrix of zeros**.

Either way `net$count` and `net$weight` carry no information: the aggregation upstream
intended never reaches the result. Verified for `sources.use`/`targets.use`, `signaling`, and
`pairLR.use` filters alike, at four thresholds, and with level sets that are not
alphabetically ordered. A drop-in replacement that returned real numbers here would *differ*
from upstream, so `r_core::net::aggregate_net_filtered` reproduces the empty/zero result and
`docs/SEMANTICS.md` records why. Pinned by
`crates/r-core/tests/netfiltered_parity.rs` and re-checked through the installed shim.

The parts of the branch that *are* correct, and are pinned alongside it because they would
matter if the regex were fixed:

* `dplyr::group_by` on the joined string key orders the groups **byte-wise**, not by factor
  level and not in locale collation order. With levels `c("g1", "g10", "g2", "g3")` the key
  order is `g10|g1, g10|g10, g10|g2, g10|g3, g1|g1, g1|g10, ...`: `"g10|g1"` precedes
  `"g1|g1"` because `'0' < '|'`, and `"g1|g10"` precedes `"g1|g2"` because `'1' < '2'`. A port
  that groups by `(source, target)` passes on a 3-level fixture and fails at 10 groups.
* `group_by` leaves one row per pair, so the two `tapply`s sum **length-1** vectors: they are
  a scatter, not a reduction.

## `filterCommunication`

Ten R-isms. The numerics are in `crates/r-core/src/filter.rs`, pinned by
`crates/r-core/tests/filter_parity.rs` against a corpus generated from the pinned upstream
(`tests/parity/gen_filter_golden.R`, 13 cases / 70 records) and re-checked end to end
through the installed shim by `tests/parity/check_identical.R`. Four of these were found
only by the *installed-shim* gate, after the Rust test was already green — the pattern
repeated from `aggregateNet` and is the reason the shim gate is not optional.

1. **The message order is not the computation order.** Upstream raises the
   `min.samples` precondition *after* printing the `min.cells` line and the per-sample
   lines, because the check is written in the middle of the function where the data it
   needs is already in hand. A port that hoists the check to the top is *more* correct and
   prints less. The sequence is: `nonFilter.keep` → `min.cells` → the `stop()` → one line
   per sample → the cross-sample percentage. The core therefore does not raise this one at
   all: it returns `min_samples_ok = false` and lets the shim raise, because the core cannot
   produce the preceding output. The other failure — the all-zero network — *is* raised by
   the core, because nothing precedes it.
2. **`scales::percent(x, accuracy = .1)` reduces to**
   `paste0(format(round(100 * x, 1), nsmall = 1, scientific = FALSE), "%")`, and
   `format()` pads to `nsmall` so the result is `"0.0%"` and `"100.0%"`, never `"0%"` or
   `"100%"`. `NaN` formats as `"NaN%"`, which is what `(0 - 0) / 0` needs.
3. **`cat()` is called twice for the percentage line**, once ending in `"!\t"` and once
   ending in `"!\n"`. So the tab and the percentage are on the same *logical* line but in
   different write calls, and `capture.output` reports them as separate entries. Comparing
   per-line output requires joining first; comparing raw output does not.
4. **`toString()` separates with `", "` — one space, no trailing space.** Trivial, and it is
   the only place a group-name list appears.
5. **`split(x, cumsum(...))` reorders the group names.** `split()` factorises its `f`
   argument, so the groups come back in *sorted* factor order. With a flat vector of group
   indices per sample that silently permutes the names: upstream prints `g1, g3` and a
   `split()`-based shim printed `g3, g1`. The per-sample loop uses index arithmetic, and
   the corpus pins group levels whose sorted order differs from their level order
   (`g1, g2, g10`) so the two cannot be confused.
6. **The complex-subunit table is read row-major by the binding.** `as.character()` on a
   matrix is column-major, so handing `cx[, subunit_cols]` straight to a kernel that indexes
   `row * n_cols + col` transposes the table: every complex gets the *next* one's subunits,
   `computeExpr_complex`'s geometric mean is taken over the wrong genes, the consistency
   mask keeps a different set of pairs, and the network comes out with 21 significant
   entries where upstream has 5. It does not error, and it is invisible on a fixture with a
   single complex.
7. **`net$prob` is a 3-d array, not a `k x (k*nLR)` matrix.** Upstream's
   `matrix(c(net$prob), ...)` *always* produces a matrix, so upstream's own
   `dim(net$prob) <- dim(net$pval)` is a no-op that only ever worked because the matrix
   already had 3 dimensions' worth of length. Rebuilding it with `dimnames=` requires
   `array()`: a matrix with the right length and the wrong number of dimnames is rejected
   with `length of 'dimnames' [3] must match that of 'dim'`. The flat order is R's
   column-major with the L-R axis fastest, which is what `as.numeric()` already delivers.
8. **A group absent from a sample counts as zero cells in that sample,** and is excluded
   for `min.samples >= 2`. `droplevels()` on the per-sample `idents` is what makes it
   absent, so the sample's group set is the set actually observed there — not the global
   `levels(object@idents)`.
9. **`nonFilter.keep = TRUE` is a no-op.** The augmented `net` *is* built and *is* written to
   `net$cellExcludes`; it is then discarded, because the `if (length(cell.excludes) > 0)`
   block ends with a second `object@net <- net` and `net` has been shadowed by the
   pre-augmentation local from the `if/else`. The core records the augmented arrays *and*
   the flag, and the shim must then use the un-augmented ones, so the test can assert on
   what the core computed without the shim being able to observe it.
10. **With `min.samples = 1` the whole per-sample half is skipped,** including
    `num.interaction0` and `num.interaction1`, so only the `min.cells` line is printed and
    the cross-sample percentage is never computed. `min.samples` is compared against
    `length(levels(sample.info))`, i.e. the *sample* factor, not the number of samples
    present in the data.
11. **All-zero network, `min.samples >= 2`:** `1:length(LR.nonzero) == 1:0` is
    `integer(0) == c(1, 0)`, which is `logical(0)`; indexing `dim(net$prob)[, 1:0, LR.nonzero]`
    with it raises `subscript out of bounds`. Preserved verbatim. Note the asymmetry with
    case 1: this one is a real indexing failure deep in the body, so the core raises it.

## `subsetCommunication`'s DEG and `netP` branches

Eleven R-isms, pinned by `crates/r-core/tests/subset_parity.rs` against a corpus generated
from the pinned upstream's own `subsetCommunication_internal` (`tests/parity/gen_subset_golden.R`,
34 cases) and re-checked through the installed shim by `tests/parity/check_identical.R`. The
input table is hand-built rather than generated: there is no RNG in this branch, so the corpus
is reproducible from the generator alone and every failure it pins is deterministic.

1. **A threshold does not drop a row with a missing value -- it blanks it.** Upstream filters
   with `net[net$ligand.pvalues <= ligand.pvalues, , drop = FALSE]`. For a row whose value is
   `NA`, the index entry is `NA`, and indexing a data frame by a logical vector containing `NA`
   **inserts a row of `NA`s**. The row is removed later, by the
   `rowSums(is.na(net)) != ncol(net)` guard, which is reached unconditionally -- and a later
   threshold applied to an all-`NA` row yields `NA` again, so the row stays blank through the
   rest of the chain and is dropped exactly once. A port that writes `filter(|r| r.p <= t)`
   gets the same *answer* and is wrong for a different reason: it also drops rows that were
   already entirely `NA` on input when **no** threshold was supplied, which is a separate
   check. The corpus has a case for each.
2. **`datasets` and `sources.use` / `targets.use` are the exceptions.** `%in%` is `match()`,
   which returns `FALSE` for a missing value and never `NA`, so those three drop the row
   outright and insert nothing.
3. **The sign of the `logFC` *argument* picks the comparison**, not the data:
   `if (ligand.logFC >= 0) ... >= else ... <=`. So a negative threshold keeps the
   *downregulated* genes, and `ligand.logFC = 0` takes the `>=` branch. It also means a
   permissive threshold for that column does not exist: the column minimum is negative, and
   using it as the argument selects the `<=` branch and keeps one row.
4. **The order of the eight `if`s is observable.** Each is an independent
   `if (!is.null(...))`, and each has its own `stop()` when the column is absent -- so asking
   for two thresholds whose columns are both missing raises the *first* one's message.
5. **The two `stop()` wordings differ.** Eight thresholds say `before using the threshold
   'ligand.pvalues'`; `datasets` says `before selecting 'datasets'`. A corpus that asserted
   one string for all nine would pass on eight cases and be wrong.
6. **`group_by` sorts on the *pasted* key, byte-wise.** The key is
   `paste(source, target, sep = "sourceTotarget")`, so `"g10sourceTotargetg1"` precedes
   `"g1sourceTotargetg1"` -- the second byte is `'0'` against `'s'`. Sorting by
   `(source, target, pathway)` as separate strings gives `"g1"` before `"g10"`, i.e. a
   *different* row order that looks correct. A fixture with three groups whose level order and
   byte order disagree is what keeps the two apart.
7. **`paste` turns a missing cell into the string `"NA"`, not into `NA_character_`.** It
   coerces with `as.character`, and `as.character(NA)` is `"NA"`. So a row with a missing
   `source` groups under a real key and comes back with `source` as that string. It is the one
   place in this function where a missing value survives as data, and it is why the all-`NA`
   row drop has to run *before* the aggregation.
8. **`mean(pval)` and `sum(prob)` accumulate in LONG_DOUBLE**, and they propagate `NA` and
   `NaN` *differently*: a group containing an `NA` reports `NA`, a group containing only `NaN`
   reports `NaN`. Both are NaN bit patterns, so a single `is_nan` cannot tell them apart and
   the two have to be tracked as separate flags. The kernel returns `NA` as a missing cell and
   `NaN` as a cell reading `"NaN"` -- *not* as the corpus's `<NA>` / `<NaN>` escape, which is a
   file format and would collapse both to `NaN` on the way back.
9. **The `netP` result is a tibble, not a data frame.** It is the only branch that passes
   through `dplyr::group_by |> summarize`, and that is what makes it a `tbl_df`; every other
   branch ends at `BiocGenerics::as.data.frame` and stays a plain `data.frame`. The class is
   part of the return value, so `identical()` sees it, and reproducing it needs
   `tibble::as_tibble` rather than a hand-built class vector -- a tibble's `row.names` are a
   compact `c(NA, -n)` integer vector, not `1:n`.
10. **`datasets` survives the final column selection only alongside `ligand.logFC`.** The
    `datasets` branch of the final `if/else` requires *both* columns, so a table carrying
    `datasets` and no `ligand.logFC` silently loses `datasets`. It reads like a bug and is the
    code.
11. **An unrecognised `slot.name` keeps every column.** The final `if/else` compares with `==`,
    so a value that is neither `"net"` nor `"netP"` falls through both branches and returns the
    table untouched.

Two `extendr` constraints shaped the binding, and both cost real time:

* **`NA_real_` cannot be passed for an `f64` parameter.** `extendr` 0.9 raises
  `Error::MustNotBeNA` during *argument conversion*, before the Rust body runs, and its message
  is the fixed string `"Must not be NA."` -- naming neither the argument nor the function. A
  "threshold not supplied" sentinel of `NA_real_` therefore never reaches the kernel; it aborts
  the call. `NaN` converts cleanly, and every threshold upstream accepts is finite, so `NaN` is
  the sentinel. `opt()` in the binding and `.cellchatrs_opt()` in the shim both say so.
* **An `Err` from an `#[extendr]` function loses its message**, becoming the same
  `MustNotBeNA` text. Since the message *is* the contract for these seven error cases, the
  kernel's error travels back as a value under `__error` and the shim raises it with `stop()`,
  which is also closer to upstream than an `extendr` condition: upstream's own `stop()` is a
  simpleError with no call, and `conditionMessage()` is what the gate compares.

The `netP` path also emits two `dplyr`/`tidyselect` deprecation warnings
("Using an external vector in selections was deprecated", "Setting row names on a tibble is
deprecated"). Those belong to the *installed dependency*, pinned to a dplyr release rather than
to CellChat, so the shim does not reproduce them: emitting another package's
version-specific text would be wrong on the next dplyr release. The Rust test matches the
dependency warnings by **prefix** and fails on any warning it does not recognise, so a genuinely
new observable still cannot slip through.

## `pSum < 0` flags every pathway whose total exceeds 1

`rankNet`'s rescaling is `pSum <- -1/log(pSum)`, and the reassignment targets
`which(is.infinite(pSum) | pSum < 0)`. The `-Inf` arm is the familiar one -- a total of exactly 1
gives `log(1) == 0` -- and it is natural to read the `< 0` arm as defensive. It is not: `-1/log(x)`
is **negative for every `x > 1`**, so the flag covers the whole range of totals at or above 1, and
`x == 1` is only its boundary.

That is not a corner case. `pSum.original` is `apply(prob, 3, sum)`, a sum over a `k x k` slice of
probabilities, so four cells at 0.5 already total 2. A port that flags only the exact-1 case
undercounts the reassignment, produces a shorter `values.assign`, and orders the bars differently.

In comparison mode the flag is also **pooled across comparisons**: one `values.assign` from
`max(unlist(pSum))` over all of them, one `pSum.original.all`, one `position` from sorting that,
and each comparison's slice addressed by an offset into the pooled sequence. Running the
reassignment per comparison -- the obvious reading, and what `mode = "single"` degenerates to --
gives different assigned values whenever more than one comparison has a flagged pathway, and
therefore a different row order.

## `rankNet`'s information flow

Eleven R-isms, pinned by `crates/r-core/tests/ranknet_parity.rs` against a corpus built by
`tests/parity/gen_ranknet_golden.R`, which **lifts upstream's own expressions verbatim** rather
than reimplementing them, so the oracle is upstream's arithmetic and not a reading of it. There is
no RNG: every probability array is written out, which is what makes the degenerate cases reachable
on purpose. Re-checked through the installed shim by `tests/parity/check_identical.R`.

1. **`-1/log(x)` has two degenerate cases and they are not the same case.** `log(0)` is `-Inf`,
   so `-1/-Inf` is `+0` -- an ordinary-looking value that is **not** caught by
   `is.infinite(pSum) | pSum < 0`, and so is never reassigned. `log(1)` is `+0`, so `-1/0` is
   `-Inf`, which *is* caught. And `log(x)` for `x < 0` is `NaN`, which
   `pSum[is.na(pSum)] <- 0` has already turned into `0` *before* the flag test runs -- so a
   negative total also escapes reassignment. The net rule, which is not obvious from the code:
   **only a total above 1 is flagged.**
2. **The transform inverts the ordering.** A larger total gives a *smaller* scaled value, because
   `-1/log` is decreasing on `(0, 1)`. The bars in a `rankNet` plot are therefore ordered
   inversely to the totals, and the reassignment exists to lift the degenerate ones back above
   the rest.
3. **The reassignment preserves the original order of the degenerate entries.**
   `position <- sort(pSum.original[idx1], index.return = TRUE)$ix` is the permutation that sorts
   the originals, and `values.assign[match(1:length(idx1), position)]` *inverts* it: entry `i`
   gets `values.assign` at its own rank. Reading it as "the sorted values in order" reverses the
   assignment, which is invisible on a symmetric set and wrong on any other. The corpus's
   `weight_three_degenerate` has four flagged entries with totals `1, 2, 2.5, 2` precisely so the
   tie in the middle makes the reversal detectable.
4. **`max(pSum)` is taken *before* the reassignment**, so it is the largest value among the
   entries that were *not* flagged. Using the post-reassignment maximum is a self-referential
   mistake that is off by the reassignment itself and still looks plausible.
5. **`seq(a, b, length.out = 1)` is `a`**, not the midpoint. It shows up whenever exactly one
   entry is degenerate, which is the common case.
6. **When *every* pathway is degenerate, `max(pSum)` is `-Inf` and R's `seq()` refuses it**, with
   `'from' must be a finite number`. `max()` over a numeric vector excludes `NA` but not `-Inf`.
   A port that computes the sequence itself yields `NaN` and fills the plot with missing bars
   instead of failing, which is the kind of divergence nobody reports until a user sees a blank
   figure.
7. **`order()` on a double vector is radix and therefore stable.** `rankNet` sorts by
   `contribution` and then compares `df$name == pair.name[i]`, so with a stable sort tied
   pathways keep their array order; an unstable sort permutes them, which changes `df`'s row
   order and therefore the factor levels of `df$name` *and* the result's row names.
8. **The zero-dropping loop leaves gaps in the row names.** `df[-which(...), ]` keeps the
   original integer row names, so a result can have row names `1 2 5 6`. `identical()` sees
   them.
9. **`sources.use` and `targets.use` are both validated against `dimnames(prob)[[1]]`** -- the
   *source* axis. A target naming a group that exists only on the second axis is rejected, and
   one naming a first-axis group that is not a target is accepted and then matches nothing. Both
   reproduced: "fixing" the second axis would change which calls work.
10. **The two group filters index different axes of the column-major array**: the source axis is
    the outer index and the target axis the middle one. Getting them the wrong way round is
    invisible on a symmetric `prob` and is caught by a fixture whose two axes differ.
11. **`contribution.relative` is rounded to one *significant* digit** before it is used as a sort
    key or compared against `1 ± tol`:
    `as.numeric(format(df[[n-i+1]]$contribution/df[[1]]$contribution, digits = 1))`. `format`
    returns *character* and `as.numeric` parses it back, so what survives is the rounding: two
    ratios agreeing to three digits are **equal** afterwards, and a port keeping full precision
    orders rows differently. `Inf` survives the round trip (`as.numeric("Inf")` is `Inf`); only
    `NaN` and `NA` become `NA` and are then reset to 0. The same line's `format` pads the vector
    to a common width, so the *characters* depend on the neighbours -- which is why the corpus
    records the parsed values and not the text.

A fixture-construction trap worth recording, because it produced three convincing false failures:
`meta[3]` is `dimnames(prob)[[3]]`, the *pathway* axis, and the two group axes are `k` long and
live in their own record. Reading the pathway names as the group levels works for any fixture
where `k == n` and fails for every one where they differ -- so it passes the first case and
misreports all the rest.

## `raw.use`

`computeCommunProb(raw.use = FALSE)` reads `object@data.smooth` instead of
`object@data.signaling` and then takes **exactly the same path** — there is no second kernel,
because upstream's own `if (raw.use) ... else ...` is two lines that choose a matrix. What is
*not* implemented is `projectData`, which produces `data.smooth` in the first place, so the
caller supplies the slot. Pinned by 9 corpus cases through the shim and 6 more in
`filter_parity.rs`.

Three things the corpus established, none of them obvious:

1. **With `min.samples = NULL` the two settings give identical results.** The `min.cells` half
   *looks* like it reads the data — it computes per-group means — but the exclusion test is
   `length(which(idents == i)) <= min.cells`, a count of **cells**, so the matrix is not
   consulted at all on that path. A port that computed the exclusion from expression values would
   differ here, and a corpus with only multi-sample cases would not notice.
2. **The `data <- data/max(data)` scaling applies to the *selected* matrix.** The two matrices
   have different maxima, so scaling by the wrong one divides every group mean by a different
   constant. That is invisible while the cross-sample mask is insensitive to it and wrong as soon
   as a per-sample score crosses the threshold.
3. **Smoothing can empty `LR.nonzero`.** A `data.smooth` built by shrinking towards the group
   mean is *constant within each group*, so every cell of a group scores identically, the
   cross-sample comparison ties, and `LR.nonzero` comes out empty — which is upstream's
   `1:length(LR.nonzero) == 1:0` → `"subscript out of bounds"` crash. Two of the six
   `raw.use = FALSE` fixtures hit it. That is a property of real projected data, not of the
   fixture, and it is worth knowing before someone reports it as a port bug.

The `data.smooth` guard is `if ("data.smooth" %in% methods::slotNames(object) == FALSE)` — with
**no** `!`. `%in%` binds tighter than `==`, so it reads "(the slot is in `slotNames`) == FALSE",
i.e. the slot is absent. This shim carried a spurious `!` for a while, on the reading that the
condition was doubly negated; that inverted the test, so every *valid* `raw.use = FALSE` call
failed with the "`object@data.smooth` is missing" message. The corpus did not catch it, because
`slotNames` reports a class's **declared** slots rather than the assigned ones, and the fixture
class declares `data.smooth = "ANY"` — so removing the assignment is not the same as removing the
slot, and the guard's firing direction is unreachable from that fixture. Exercising it needs a
class that does not declare the slot; until then the guard is covered by reading, not by a test,
and `filter_parity.rs` says so rather than asserting something that passes for the wrong reason.

## `computeRegionDistance`

`crates/r-core/src/spatial.rs`, in two layers. The **foundations** -- R's `mean(x, trim, na.rm)`,
`collapse::fdist`, and the exact k-d tree that replaces `AnnoyParam` -- are pinned by
`crates/r-core/tests/spatial_parity.rs` against a corpus generated from R
(`tests/parity/gen_spatial_golden.R`, 26 trimmed-mean cases and 7 distance matrices).

The **function** is pinned by `crates/r-core/tests/region_parity.rs` against
`tests/fixtures/region_golden.txt` (13 cases), and the R-side marshalling by the
`computeCellDistance` and `computeRegionDistance` blocks of `tests/parity/check_identical.R`.

### The corpus is an exact-substitution oracle, and says so

Upstream's `computeRegionDistance` cannot be *run* here at all: it evaluates
`BiocNeighbors::queryKNN(..., AnnoyParam())` unconditionally and the package is not installable
without network. Annoy is randomised and approximate, so upstream is not an oracle for a neighbour
query in any case -- two runs can disagree.

So `gen_region_golden.R` lifts upstream's body **verbatim** with exactly two expressions replaced
by exhaustive search, and each fixture records its own nearest-neighbour margin `d1 / d2`. The
substitution is sound only when the runner-up is decisively farther than the nearest neighbour, and
`the_recorded_margins_hold_up` re-checks the recorded number in Rust, so a layout edit that quietly
destroys the property fails in the test suite and not only in the generator. The 13 fixtures range
from 2.4x to 590x.

Getting that margin right took three attempts, and the failures are the useful part. Separating two
groups by a **large** gap gives `d1/d2 = 0.995` -- the opposite of intended -- because a big gap
makes every cell of the other group nearly equidistant, so the nearest neighbour is the *least*
determined thing in the fixture. Widening the gap made it worse. Geometric ray spacing gives
`d2 - d1 = f - 1` for every cell, which makes the ratio predictable instead.

### Six behaviours that read as bugs and are not

1. **`d.spatial` is symmetrised unconditionally.** Upstream's
   `d.spatial <- (d.spatial + t(d.spatial))/2` sits *outside* the `if (do.symmetric)` block that
   guards the three `adj.*` symmetrisations. Read as a block it looks like an oversight. Moving it
   inside -- which is what the port did first -- produced a `d.spatial` differing from upstream in
   three of nine entries on the `not_symmetric` fixture. `do.symmetric = FALSE` therefore leaves the
   returned `adj.contact` genuinely asymmetric and `d.spatial` symmetric.
2. **A declared level with no cells gets the wrong row name.** Upstream sizes the arrays with
   `numCluster <- nlevels(group)` but filters the level *names* into `level.use`, then loops
   `for (i in 1:numCluster)` reading `level.use[i]`. The loop counter picks the slot and `level.use`
   picks the level, and they are different counters. With `levels = c("C","A","B")` and only A and B
   populated, **A's row lands in slot 1**, which `rownames(d.spatial) <- levels(group)` labels "C";
   slot 3 is never written. The port reproduces the mislabelling. Collapsing the two counters into
   "the present levels" -- the obvious reading -- yields the transpose, with every value correct.
3. **`contact.range = NULL` is an empty selection, not a skipped pair.** `dist - NULL` is
   `numeric(0)`, so `which(numeric(0) < tol)` selects nothing and `adj.contact` is 0 throughout.
   This is the *normal* case when `contact.knn.k` is given instead, so treating it as "skip this
   pair" left `adj.spatial` and `d.spatial` unwritten -- all `NaN` -- for every fixture using
   `contact.knn.k`.
4. **`contact.range <- 10000` does not make every `adj.contact` entry 1.** Upstream's comment says
   so, and it is true of the *range test*: every query then passes it. But the entry is
   `length(unique(qout$index[selected rows])) >= k.min`, a count of **distinct target cells**, and
   three source cells sharing one nearest target give 1, so with `k.min = 2` the entry is 0.
   Reading the comment as a guarantee is what makes this look like a port bug.
5. **`ratio = NULL` empties the selection.** `qout$distance * NULL` is `numeric(0)`, so
   `FunMean(numeric(0))` is `NaN` and both adjacency tests select nothing. That is the *default*
   signature of `computeRegionDistance`, and a length-1 `ratio` is recycled the way R's `ratio[k]`
   is -- so a scalar `ratio` with several samples is the ordinary case, not an edge case.
6. **`computeCellDistance` returns a `dist`, and the rewrap is not `as.dist`.** `collapse::fdist`
   returns a `dist` whose attributes differ from `stats::as.dist`'s in three places: no `call`
   attribute, `method = "euclidean"` set, and `Labels` left `NULL` when the input has no row names
   (where `as.dist` synthesises `as.character(seq_len(n))`). `computeCellDistance` also raises
   `fdist`'s own "at least 2 rows" error for a single cell, *after* the two-column check.

### `computeCellDistance`'s other R-side details

`ratio` is applied when non-`NULL`, so `ratio = numeric(0)` *is* multiplied in -- unlike
`computeRegionDistance`, where `NULL` empties the vector. The threshold fires only when **both**
`interaction.range` and `tol` are given, and `tol` is an **addend** rather than a tolerance:
`> interaction.range + tol`. `colname` assignment to `c("x_cent", "y_cent")` happens before anything
can fail, which is why a one-column matrix reports the column error rather than anything about
coordinates.

### Flat-vector marshalling is row-major, and R's is column-major

Both kernels take a flat coordinate vector and index it row-major. `as.numeric(matrix)` is
**column-major**, so passing it straight through hands the kernel the transpose. For `fdist` the
transpose is invisible -- the output is symmetric, so every value matches -- which is exactly why
it has to be written down: the same mistake in `computeRegionDistance` moves every cell to a
different place and the k-d tree answers confidently about the wrong layout. The shims pass
`as.numeric(t(coordinates))`, with the reasoning inline.

`extendr` also derives the R wrapper's parameter names from the Rust ones, so they are snake_case:
passing `interaction.range` to `spatial_region_distance` is an "unused argument" error raised from
*inside* the shim, so the traceback points at the kernel call rather than at the mismatch.

### R's `mean(x, trim, na.rm)`

Reproduced from `base::mean.default`'s body rather than from a description of a trimmed mean,
because the two differ in three places that matter:

1. **`na.rm = TRUE` removes `NaN` as well as `NA`**, since `is.na(NaN)` is `TRUE`. `NaN` is a
   *value* everywhere else in R, so a port that filters only `NA` returns `NaN` where R returns
   a mean. `mean(c(1, NaN, 3), na.rm = TRUE)` is `2`.
2. **The count trimmed from each end is `floor(n * trim)`,** so the kept count is
   `n - 2 * floor(n * trim)`. It is not `round(n * (1 - trim))`, and it is **not monotone in
   `n`**: at `trim = 0.1`, `n = 9` trims nothing, `n = 10` trims one, `n = 19` one, `n = 20` two.
   The corpus pins all four.
3. **`trim >= 0.5` returns `median`**, and on an even `n` that is the mean of the two middle
   values. `computeRegionDistance` uses `trim = 0.1`, so this branch is unreachable there; it is
   implemented because `mean` is a general utility.

Plus the `anyNA` early return, which is unreachable when `na.rm = TRUE` and gives `NA` -- not a
mean of the rest -- when it is `FALSE` *and* `trim > 0`. With `trim == 0` the whole trimming
block is skipped and a `NaN` is summed through instead, so the two settings differ. The
accumulation is LONG_DOUBLE.

### `collapse::fdist` returns a `dist`, not a matrix

`computeCellDistance`'s `d.spatial <- collapse::fdist(coordinates)` is a **`dist`**: the compact
triangle, with `Size` / `Diag` / `Upper` / `method` attributes. The two statements that follow it
therefore go through `dist`'s own methods, and `d.spatial[d.spatial > x] <- NaN` writes only into
the half a `dist` exposes. `fdist` also refuses a one-row matrix outright, and the diagonal is
`0` for a matrix but `NA_real_` for a data frame -- `computeCellDistance` passes a matrix.

The corpus records the full symmetric matrix (computed in plain R) *and* `as.vector(d)`. The
matrix is compared bit for bit, which pins the arithmetic; `as.vector(d)` is compared as a
**multiset**, because its traversal is not documented and `collapse`'s is not base R's -- for a
3-4-5 triangle `collapse` orders the pairs `(1,2), (2,3), (1,3)` where base orders
`(1,2), (1,3), (2,3)`. Pinning the traversal would pin a property of `as.vector` rather than of
the distances, and the `dist` wrapper is re-applied on the R side.

### The exact k-d tree

`Annoy` is a randomised approximate index, so the spatial branch is the one place where a
bit-identical port is *impossible by construction*: two runs of upstream can disagree, and
upstream is not a usable oracle for a neighbour query. The locked decision (`PLAN.md` §14.5) is
to replace it with an exact index and **measure** the divergence on real spatial data.

So the tree is verified a different way: against **exhaustive search**, on point sets built to
break a k-d tree rather than to exercise one. `knn_matches_brute_force_on_adversarial_point_sets`
covers a regular lattice (every neighbour tied at distance 1), an all-duplicate set (every
distance tied at 0), a collinear run (one axis with no spread), a circle (points at equal radius),
two tight clusters far apart (so the split axis matters), and uniform 3-d noise — at
`k` in {1, 2, 3, 5, 8, 17}, querying from points in the set (so distance 0 occurs) and from
off-lattice points. A tree that is wrong on uniform noise is usually right on uniform noise.

Ties break by **lower point index**, which makes `(distance, index)` a total order. That is not
cosmetic: `computeRegionDistance` counts `length(intersect(knn.i, unique(qout$index[...]))`, and
`unique` keeps the *first* occurrence, so with equidistant candidates the chosen set is stable
but the order is not, and the order is what the count sees.

Five bugs were found on the way, and each is a trap rather than a typo:

1. **An arena, not a `2v+1 / 2v+2` heap layout.** Filling a `Vec` by node index with
   `resize(node + 1, 0)` **truncates** whenever a deep node is written before a shallow one —
   which an explicit stack guarantees. The tree ends up with holes where whole subtrees were, and
   the search walks into them and returns a plausible point. An arena, appended to, cannot fail
   this way; `the_arena_has_one_node_per_point` is the cheap assertion that would have caught it.
2. **Near/far sign.** With `delta = point - query`, `delta > 0` puts the split plane to the
   *right* of the query, so the query is on the left and the **left** child is near. The more
   common `query - point` spelling gives the opposite, and the search then descends the far
   subtree first — still a neighbour, just the wrong one.
3. **Only the far subtree's bound is tightened by `delta^2`.** The near subtree contains the
   query's own coordinate and inherits its parent's bound *unchanged*. Tightening both is the
   symmetric-looking reading, and it makes the near bound too large, so the search prunes the
   subtree that holds the answer.
4. **Insert-then-evict must compare against the root.** `push` followed by
   `swap(0, len-1); pop()` removes `heap[len-1]` *after* the swap — the old root, which for
   `k == 1` is the old *best*. Every candidate after the first then replaced the answer wholesale
   and "nearest" was simply the last point visited.
5. **An empty range is not a leaf.** `mid = lo + (hi - lo) / 2` leaves `[mid+1, hi)` empty
   whenever `mid + 1 == hi`, and `hi - lo == 0` is not `== 1`, so a `build` signature returning a
   bare `usize` recurses forever. This one at least fails loudly — a stack overflow rather than a
   wrong answer.

All five return a *plausible* point rather than failing, which is why the suite's centre of
gravity is exhaustive search on adversarial inputs and not spot checks. `the_search_is_subquadratic`
is there because a quadratic tree would satisfy every correctness test here and still be a
regression: `computeCommunProb` runs a 1-NN query per cell per cell group, so `O(n log n)`
against `O(n^2)` is the entire reason for replacing Annoy.

## R's `order` with a missing value in a *second* key

`rankNetPairwise` is one line of arithmetic:

```r
temp[[j]] <- data[with(data, order(pval, -prob)), ]
```

`order` with two numeric keys and the default `method = "auto"` is a **stable radix sort**, so the
permutation is a total function of `(pval, -prob, index)`. That makes it a faithful candidate for
Rust, and the port's version of the function is otherwise a pass-through to upstream.

It was wrong for two reasons, and both were invisible without *repeated* keys.

**A missing value in a later key compared as equal, not as greater.** Measured against R:

```text
order(c(1,1,2),  c(NA,5,3))   -> 2 1 3
order(c(1,1,2),  c(5,NA,3))   -> 1 2 3
```

The `NA` in the second key sorts **last within its group** of equal first keys. The port's
comparator treated `partial_cmp` returning `None` as "equal", which is `1 2 3` -- the input order.
So for `order(pval, -prob)` with a repeated `pval` and one `NaN` probability, the `NaN` stayed where
it was and R moved it to the end of the tie group. A different row order, from a plausible-looking
implementation.

**An all-missing primary key stopped the comparison dead.** `order(c(NA,NA,NA), c(3,1,2))` is
`2 3 1` in R: `na.last = TRUE` moves the missing values to the end, and the *remaining keys still
order them*. The port returned `Equal` as soon as both primary keys were missing, so it produced
index order. Missing values now sort last in **every** key, and when both sides of a key are missing
the comparison falls through to the next key rather than declaring a tie.

`order_f64_multi` is shared with `rankNet`, and the whole rankNet suite (274 Rust tests, 377 shim
comparisons) passes unchanged after the fix, which is the evidence that the old behaviour was simply
wrong rather than differently right.

## Upstream's `rankNetPairwise` reads a `k x k` block and does not reconcile `LR.use`

Three shapes in the ported function turned out to be load-bearing, and each one is an error the port
has to raise **from its own code**:

* **`LR.use` shorter than `dim(prob)[3]`.** Upstream uses `pairLR.use` for the data-frame's columns
  and `prob[i, j, ]` for the values, and never checks that they agree. Its
  `data.frame(..., row.names = rownames(pairLR.use))` raises `row names supplied are of the wrong
  length`. Slicing `prob` down to `LR.use` first -- the obvious reading of "use this subset of
  interactions" -- turns that into a successful shorter result.
* **A non-square `prob`.** `numCluster` is `dim(prob)[1]` and both loops run `1:numCluster`, so
  upstream reads the leading `k x k` block; then `names(temp) <- colnames(prob)` has the wrong length
  and raises `'names' attribute [5] must be the same length as the vector [2]`. The block is taken
  on `prob`/`pval` only -- `rownames(prob)` and `colnames(prob)` stay the originals, which is exactly
  where upstream fails.
* **`k = 1` and `n = 1`.** Upstream writes `data[idx, ]` with no `drop = FALSE`, so a
  single-interaction network comes back as a vector rather than a data frame. The port's
  `data[idx, , drop = FALSE]` would return a one-row data frame, so the `drop` argument is chosen to
  reproduce upstream.

## A pass-through and a port are both "identical" to upstream

The shim's first version answered from upstream whenever the kernel's result did not fit -- a
`length(ord)` mismatch, a declined shape, an empty interaction list. Every differential test passed,
because upstream's answer is by definition identical to upstream's. The gate was reporting success
for a function that was not ported at all, and it had been doing so for a real bug: the kernel read
R's `dim()` as a list of slice lengths, summed it to `k + k + n` instead of `k * k * n`, returned
`__error`, and the shim quietly answered from upstream.

Three changes, and the third is the one that matters:

1. the kernel's shape handling is explicit and correct, so the fallback is not being taken;
2. every fallback goes through a helper that **warns**, naming the reason; and
3. the gate counts those warnings as failures, so a pass-through can never again be reported as a
   port. `check_identical.R` runs 18 `rankNetPairwise` cases and 13 of them are `identical` while 5
   are upstream's own errors reproduced by the port's own code -- with zero fallbacks.

The `chunk_exact`-versus-stride bug in the same kernel is worth recording for the same reason. A
`k x k x n` array flattens column-major, so the first index varies fastest and a fixed `(i, j)` slice
is a **strided** gather with stride `k * k`, not a contiguous run:

```r
a <- array(1:8, c(2, 2, 2)); as.vector(a)
# 1 2 3 4 5 6 7 8  ==  a[1,1,1] a[2,1,1] a[1,2,1] a[2,2,1] a[1,1,2] ...
```

`flat.chunks_exact(n)` -- the obvious reading -- returns the run that varies `i`, so every slice
after the first is transposed. That version passed a one-slice test 200 times and failed 539 of 540
slices at `c(3, 3, 6)`. A test with a single group pair cannot see a slice-ordering bug; the gate
now covers `c(2, 5, 3)`, `c(5, 3, 4)` and `c(3, 1, 2)` as well as the square shapes.

## `group.dataset` was never blocked on presto, and it was not a small change

`identifyOverExpressedGenes(group.dataset = ...)` sat in the ledger as "the dataset-comparison
branch", filed next to `do.fast = TRUE` and therefore next to presto. It is not blocked. Upstream
honours `group.dataset` in **both** halves of its `if (do.fast)`:

* `utilities.R:429-484` -- the presto branch, and
* `utilities.R:512-520` -- the Wilcoxon branch.

Only the first needs presto, so `do.fast = FALSE` with `group.dataset` set was testable all along.
The ledger entry was wrong, and "wrong because it looked like its neighbour" is a failure mode worth
recording.

What the branch actually changes is one thing -- how the two cell sets are chosen per group:

```r
if (is.null(group.dataset)) {
  cell.use1 <- which(labels == level.use[i])
  cell.use2 <- base::setdiff(1:length(labels), cell.use1)
} else if ((!is.null(group.dataset)) & (group.DE.combined == FALSE)) {
  cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
  cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
} else if ((!is.null(group.dataset)) & (group.DE.combined == TRUE)) {
  cell.use1 <- which(labels.dataset == pos.dataset)
  cell.use2 <- which(labels.dataset != pos.dataset)
}
```

The percentage filter, `mean.fxn`, `wilcox.test` and the Bonferroni multiplier are identical
whichever way the sets are chosen. So the kernel takes the **selection as a parameter**
(`identify_over_expressed_selected`) rather than growing a `group.dataset` flag through the whole
Wilcoxon path -- the same reasoning as `rankNetPairwise`, where one `order` was the only arithmetic.

Four details in the plumbing around it are load-bearing, and each was a bug before it was a rule.

**The `toString` collapse.** `labels.dataset[labels.dataset != pos.dataset] <- toString(setdiff(
unique(labels.dataset), pos.dataset))` puts *every* non-positive dataset into a single level. A
four-dataset comparison therefore has **two** levels, one of them named `"D2, D3, D4"` -- `toString`
joins with `", "`. Keeping the datasets separate gives a different `datasets` factor and a different
row order.

**A single candidate feature means no markers at all.** `apply(X, 1, FUN)` on a one-row matrix
returns an unnamed *scalar*:

```r
a  <- apply(matrix(1:4, nrow = 1), 1, f)   # 3, and names(a) is NULL
FC <- a - a2
FC["gene1"]                                 # NA
names(which(FC > 0))                        # character(0)
```

so `features.diff` is empty, `intersect(features, features.diff)` is empty, and the loop's `next`
skips the group. A port that keeps feature names attached to values would report a marker here --
and "the only group whose marker count is one" is exactly the case a threshold sweep does not
sample.

**The empty `cell.use2` case has to fail.** When every cell of a group is in the positive dataset,
`cell.use2` is empty, `pct.2` is `NaN`, every feature is dropped by `which(alpha.min > thresh.pc)`, and
`markers.all` is still the bare `data.frame()`. Upstream then adds a `datasets` column to that 0 x 0
frame and dies on

```text
Error in -markers.all$logFC : invalid argument to unary operator
```

Substituting the collapsed one-column marker frame returns an empty table where upstream returns
nothing at all. The shape is part of the contract, so the shim keeps a 0 x 0 frame on this path and
a 1 x 1 one otherwise.

**`rbind` is not associative, and upstream rbinds one group at a time.** `markers.all` starts as
`data.frame()` and grows by `rbind(markers.all, gde)`. R's row-name uniquification renames a
collision to `x.1` only while `x.1` is free, and otherwise appends a bare digit -- so `x` becomes
`x1`, not `x.1`. How many earlier passes have happened therefore changes the names, and
`do.call(rbind, frames)` (uniquify the whole concatenation once) produces a different set. With
`only.pos = FALSE`, where a feature can be a marker in several groups, the two differed on two of
forty-one row names while every column matched.

One more shim bug is worth naming because it is a shape bug rather than a numeric one: the frame
loop guarded on the kernel's `n_before_only_pos` counter, which is an over-estimate by design once
`idents.use` has filtered rows away. On an *empty* filtered result,
`starts <- which(c(TRUE, logical(0)))` is `1` and `ends <- c(numeric(0), 0)` is `0`, so the frame
indexed `res$features[c(1, 0)]` and `data.frame` answered `row names contain missing values`. The
guard is now `length(res$features) > 0L`, which is what upstream's per-group
`if (nrow(gde) > 0)` actually means.

The gate runs 40 `group.dataset` configurations -- two and four datasets, `pos` first and last,
`group.DE.combined` both ways, `only.pos` both ways, the three threshold shapes, `idents.use`, the
aligned-labels case where `cell.use2` is empty, and the bare `stop()` for an unusable `pos.dataset` --
all against upstream's own function.

## Metamorphic relations: three of the six obvious statements are false

The differential gate compares the port against upstream on *one* object per configuration. That
says "the port agrees here" and nothing about whether the agreement would survive a re-labelling of
the input, because a bug that is symmetric in two coordinates is invisible while both are held
fixed. `tests/parity/metamorphic.R` applies each transformation and asserts two separate things: the
port still matches upstream *on the transformed object*, and the stated relation holds relative to
the untransformed run. 40 checks, all passing.

Three relations are not what they look like, and the obvious assertion is false:

* **Cell permutation is not an invariance of `Pval`.** The bootstrap draws `sample.int(nC, nC)` over
  *cell indices*, so permuting the input permutes which cells each replicate samples and `Pval`
  moves -- measured, it does move. Only the observed `Prob` is invariant. And the transform has to
  carry the labels with the columns: permuting the matrix and keeping the labels fixed does not
  permute anything, it reassigns each cell's expression to a different group.
* **Cluster relabelling permutes the rows *and* the columns.** `prob[i, j, l]` means group `i` sends
  to group `j`, so renaming `g1 <-> g2` swaps row 1 with row 2 *and* column 1 with column 2. The
  relation is `prob[σ(i), σ(j), l] == prob_orig[i, j, l]`, and asserting equality in place asserts
  something false.
* **Constant scaling is invariant only up to rounding, unless the factor is a power of two.**
  `data.use <- data/max(data)` is designed to cancel the factor, but `max(f * x)` is not
  `f * max(x)` in floating point, so `x / max(f * x)` and `x / max(x)` differ in the last bits.
  Measured residue: exactly 0 at `f = 0.5` and `f = 1e4`, `1.1e-16` at `f = 3` (one ulp of 0.451), and
  `1.7e-21` at `f = 7.25`. What *is* exact, for every factor, is `identical(port, upstream)` on the
  scaled object -- and that is the claim worth making.

Each block also carries a negative control, because an invariance that holds because the transform did
nothing is not a test: doubling one gene must move `Prob`, and an order-preserving relabel must
leave the values in place while the dimnames change.

The objective's four recorded invariants are checked on every object the file builds: dimnames
consistency, `Pval` in `{k/nboot}`, `Pval[Prob == 0] == 1`, and `Prob` in `[0, 1]`.

## Fuzzing the database structure found no kernel bug and four wrong properties

`crates/r-core/tests/db_fuzz.rs` generates hostile L-R databases -- empty subunit lists, names
colliding across the complex, cofactor and symbol tables, self-referential complexes, overlapping
subunits -- and checks that `resolve_entity` and `compute_expr_lr` are **total**: every entry point
returns a value or an `ExprError`, never a panic. The resolution result is compared against a
`BTreeMap` reference written independently in the test, so a shared misreading cannot make the
property vacuous.

The kernel passed all four properties. The *properties* were wrong four times, and each time the
generator had found a real rule I had mis-stated:

* a complex with an entirely empty subunit list contributes **nothing**, not even its own name;
* the empty name is dropped by `gene[gene != ""]`;
* names not in `geneInfo$Symbol` are dropped by `checkGeneSymbol` before anything else happens; and
* the split is on the **official symbol list**, not on membership of the complex table. Upstream:

```r
complex <- geneSet[which(geneSet %in% geneIfo$Symbol == "FALSE")]
geneSet  <- intersect(geneSet, geneIfo$Symbol)
```

So a symbol is kept as itself and never expanded, and a non-symbol is looked up in the complex table
and contributes only its non-empty subunits. In the real database the complex names are not symbols
and the plain gene names are, which is exactly why this is invisible by hand -- and why three
attempts at the property, each using "is in the complex table" as the discriminator, were all wrong.

## A property test found a bug in the 80-bit accumulator that 276 unit tests did not

`F80::mul` began with

```rust
if self.m == 0 || other.m == 0 {
    return F80::ZERO;
}
```

and **infinity is encoded with a zero mantissa**. So that line caught every infinite operand, and
`prod(c(-Inf))` -- which R evaluates as `-Inf` -- came out as `0`. Every multiplication involving an
infinity was wrong: `F80::from_f64(-Inf).mul(F80::from_f64(2.0))` was `+0`.

Nothing caught it. `r_prod` is only reachable through the cofactor path in `computeExpr_coreceptor`,
no fixture in the repository puts an infinity in a cofactor column, and every parity test pins a
*finite* case. A property test that generates `Inf` does not care that no fixture has one.

The fix follows R's measured table, which is **not** IEEE 754:

```text
prod(Inf, Inf)  = Inf     prod(Inf, -Inf)  = -Inf
prod(-Inf, Inf)  = -Inf    prod(-Inf, -Inf) = Inf
prod(Inf, 0)     = NaN     prod(0,  Inf)    = NaN
prod(-Inf, 3)    = -Inf    prod(-Inf, 0)    = NaN
```

An infinite product is `Inf` carrying the **XOR of the two signs** -- including `Inf * -Inf`, which
IEEE makes `NaN` and R does not -- and only a zero times an infinity is `NaN`. Writing the "obviously
correct" IEEE rule here would have been a fresh divergence.

Two smaller consequences of the same encoder. A zero product must keep the XOR sign, because
`prod(c(-0, 5))` is `-0x0p+0` and `prod(c(5, -0))` is *also* `-0x0p+0` -- multiplication is
commutative but a signed zero is not. And `F80::to_f64` returned a bare `0.0` for a zero mantissa,
so every signed zero came back positive on the way out. (`-Inf * 0.0` is `NaN`, not `-0.0`; the bit
pattern is the only way to spell a signed zero.)

The moral is about test *kinds*, and it is the reason `crates/r-core/tests/properties.rs` exists.
A named case pins an answer; a generated case checks a rule. Both bugs here produced plausible
outputs on every input the repository happened to contain, and would have shipped.

## The configuration matrix, and its result

`tests/parity/matrix.R` (axes, candidate pool, greedy covering selection, coverage report),
`tests/parity/matrix_objects.R` (the object family) and `tests/parity/check_matrix.R` (the runner and
the JSON report) are a separate artifact from `check_identical.R`: the gate asks "does each function
work", the matrix asks "does it work in every corner at once". Fourteen axes, 200 configurations,
**1051 of 1051 axis pairs covered**, and the runner asserts scale-invariance directly, because
pairwise coverage cannot express a relation *between* two rows.

```
configurations : 200
status counts  : {"error-equal": 42, "identical": 158}
coverage       : 1051 / 1051  complete = True
scale invariant: True
non-passing    : 0
```

The 42 `error-equal` rows are configurations where upstream itself raises -- a missing subunit, an
all-zero matrix, a duplicated row name -- and the port raises with the same message. Matching the
message is the requirement; returning a value where upstream refuses is not parity, which is why the
runner treats `identical` and `error-equal` as the two passing rungs and `error-differs` as a
failure.

Three things about the harness had to be right for those numbers to mean anything.

**Upstream writes to stdout, not through `message()`.** Its two `cat`s, its `print`s and
`setTxtProgressBar`'s carriage returns all bypass the condition system, so the runner sinks the
connection and captures the text through a `message` handler. The sink is popped on the error path
too, because an unbalanced `sink()` silently swallows the rest of the script's output -- and it did,
once, in a way that made a 40-minute run look like it had produced nothing.

**A `NaN` on one side and a value on the other is not a value comparison.** The runner compares
verdicts and the three comparable fields, never the whole record, because `effective$scale` differs by
construction and `identical(a, b)` would then be false for the one thing the check is meant to assert.

**When a row disagrees, the report says why.** The runner repeats the surviving side in place, in
three call shapes, and records whether the input object was the one it was built from. All of that
costs nothing on a passing sweep, and on a failing one it is the difference between a hypothesis and
a measurement. It is also how the `type`-dropping bug below was finally pinned.

## A fixture whose maximum is 1 hides an entire missing step

`computeCommunProb` begins with one line:

```r
data.use <- data/max(data)
```

The port **omitted it** for the whole life of the kernel, and nothing noticed. Every fixture in the
suite drew from `runif(0.01, 1)`, and with a few hundred draws `runif` reaches `1` often enough
that `max(data) == 1` every time -- so the division was the identity and its absence was
unobservable. The differential gate, which compares against upstream's real `computeCommunProb`,
stayed green at 359/359.

What the step actually buys is **scale-invariance**: upstream's `Prob` depends only on
`data/max(data)`, so multiplying the expression matrix by any constant leaves it bit-for-bit
unchanged. Without the division the port's `Prob` moved by more than an order of magnitude when the
same data was divided by two -- and `Kh` is fixed, so the Hill terms
`n^Kh * P^Kh / (Kh + P^Kh)` are not scale-invariant.

Two things now prevent a repeat:

* the configuration matrix carries a **`scale` axis** (`0.5`, `1`, `2`, `3`), and
* `check_matrix.R` asserts the invariance *directly* -- one configuration at `scale = 1` and the
  same at `scale = 3` -- because pairwise coverage cannot express a relation *between* two rows.

The degenerate maxima are not special-cased; each falls out of the arithmetic:

| `max(data)` | `data/max(data)` | upstream |
|---|---|---|
| finite, positive | ordinary | `Prob`, scale-invariant |
| `Inf` | finite entries `0`, the `Inf` entry `NaN` | `Prob` all zero, **no error** |
| `0` (all-zero matrix) | all `NaN` | error, see below |
| `NaN` or `NA` | all `NaN` | error, see below |

## `Matrix::nnzero` returns `NA`, and upstream's `if` cannot continue past that

`thresholdedMean` is three lines:

```r
thresholdedMean <- function(x, trim = 0.1, na.rm = TRUE) {
  percent <- Matrix::nnzero(x)/length(x)
  if (percent < trim) return(0) else return(mean(x, na.rm = na.rm))
}
```

It is the only one of the four `FunMean` variants that branches on a statistic rather than on a
value, and the statistic is not what its name suggests. Measured against `Matrix::nnzero`:

| `x` | `nnzero(x)` | `percent` | `thresholdedMean` |
|---|---|---|---|
| `rep(0, 15)` | `0` | `0` | `0` |
| `as.numeric(1:15)` | `15` | `1` | `8` |
| `c(NaN, 1, 3, 5)` | **`NA`** | `NA` | **error** |
| `c(NA, 0, 0, 0)` | **`NA`** | `NA` | **error** |
| `rep(NaN, 15)` | **`NA`** | `NA` | **error** |

A single missing entry does not reduce the count -- it makes the whole answer `NA`. So `percent` is
`NA`, and R evaluates `if (NA < trim)`, which is

```
Error in if (percent < trim) return(0) else return(mean(x, na.rm = na.rm)) :
  missing value where TRUE/FALSE needed
```

The port read `nnzero` as "count the non-zero entries and skip the missing ones", which is a
perfectly reasonable reading of the name and of the `na.rm = TRUE` sitting right there in the
signature. It made `thresholded_mean` return `0` where upstream raises, and a unit test
(`thresholded_mean_treats_nan_as_missing`) asserted that `c(NaN, 1, 3, 5)` averages to `3.0` --
locking the wrong reading in so that a later fix would look like a regression.

## Upstream applies `FunMean` to every gene, and the port threw the error away

Fixing `nnzero` was necessary and not sufficient. Upstream aggregates the whole expression matrix:

```r
data.use.avg <- aggregate(t(data.use), list(group), FUN = FunMean)   # modelling.R:115
```

so `FunMean` runs on every gene, not on the genes the interaction list mentions. The port does the
same -- it aggregates all genes -- and then reads back only the interacting ones, so the offending
`NaN` was computed and discarded. Six of the 200 matrix configurations returned an all-zero `Prob`
where upstream raises.

A `NaN` in the observed aggregate therefore has to be treated as a raise, and only for
`thresholdedMean`, because that is the only `FUN` for which a `NaN` means "upstream stopped" rather
than "this group had no finite value". For `triMean` and `median` a `NaN` group mean is an ordinary
value, and upstream's later `if (sum(P1_Pspatial) == 0)` is what converts it into the same message.

The port also computed the bootstrap replicates *before* the observed aggregate, so its error
precedence was upstream's reversed: `modelling.R:115` is the observed aggregate and `:208` is the
bootstrap. The two now run in upstream's order.

## R's `max()` is not Rust's `f64::max`, and `sum()` had no way to carry `NaN`

Both surfaced from the same fixture, and both are the kind of difference that a suite built only
from finite data cannot see.

**`f64::max` is a `maxNum`.** `1.0_f64.max(f64::NAN)` is `1.0`; R's `max(c(1, NA))` is `NA`, because
R's `max` propagates `NaN` and only ignores `NA` when `na.rm` says so. `-Inf` is a *value*, not a
missing one: `max(c(-Inf, 1))` is `1`. So `r_max` short-circuits on `NaN` and treats `-Inf` as
ordinary.

**`F80` had no representation for `NaN` or `Inf` at all** -- `from_f64` mapped every non-finite
input to `ZERO`. Since *every* reduction in this crate accumulates through `F80`, `sum(c(1, NaN))`
came out as `1`, which is not a value R can produce. `F80` now carries the two special values in
sentinel exponents, as the x87 extended format does, and `add` follows IEEE 754 and R: any `NaN`
poisons the sum, `Inf + -Inf` is `NaN`, otherwise an infinite operand wins. Finite arithmetic is
unchanged byte-for-byte; `longdouble_special_parity.rs` pins both the specials and the finite path.

The consequence that made this observable is R's own:

```r
if (sum(P1_Pspatial) == 0) { ... }
```

`NaN == 0` is `NA` in R, and `if (NA)` raises **`missing value where TRUE/FALSE needed`**. The port
has no `if`, so the comparison was simply false and the kernel carried on and *returned* an
all-`NaN` network. That is strictly worse than an error: the caller cannot tell "no communication"
from "the arithmetic went missing". `KernelError::IfNa` now raises it, checked before the zero test
because in Rust the zero test is simply false for `NaN`.

## `format(x, digits = 1)` is not a per-element function, so it stayed in R

`rankNet(mode = "comparison")` rounds its relative contributions with

```r
contribution.relative[[i]] <- as.numeric(format(df[[ncomp - i + 1]]$contribution / df[[1]]$contribution, digits = 1))
```

and then **sorts on the result**. That makes R's pretty-printer part of the numeric surface: get it
wrong and the row order changes, not just a digit.

The trap is that `digits` is not per element. Measured:

| input | `format(., digits = 1)` | `as.numeric` |
|---|---|---|
| `1.2` | `"1"` | `1` |
| `c(0.04, 1.2)` | `c("0.04", "1.20")` | `c(0.04, **1.2**)` |
| `c(Inf, 1.2, 0.4, 0, 0)` | `c("Inf", "1.2", "0.4", "0.0", "0.0")` | `c(Inf, **1.2**, 0.4, 0, 0)` |
| `c(1, 1.004, 1.04, 0.996)` | `c("1", "1", "1", "1")` | `c(1, 1, 1, 1)` |

The same value rounds differently depending on its neighbours, because `format` picks a *common*
scale for the vector and then pads to it. `signif(x, 1)` -- the obvious reimplementation -- gives
`1` for `1.2` in every case, which is right for the first row of the table and wrong for the other
three.

So the kernel returns the **unformatted** ratio and the shim applies R's own `format`, which is
exact by construction and is where the value is consumed anyway. The same reasoning applies to the
`order()`: it sorts on the formatted values, so it is R's `order` too. The corpus records both
(`cmp_ratio` and `cmp_relative`) so the distinction is visible rather than assumed.

## The comparison branch's reassignment is pooled

In `mode = "single"` the degenerate reassignment is over one network. In `mode = "comparison"` it
is **pooled across all of them**: one `values.assign` from `max(unlist(pSum))`, one
`pSum.original.all` concatenating every comparison's flagged entries in comparison order, and one
`position` from sorting that. Each comparison's slice is then addressed by an *offset*:

```r
pSum[[i]][idx[[i]]] <- values.assign[match(length(unlist(idx[1:i-1])) + 1:length(unlist(idx[1:i])), position)]
```

Running it per comparison -- the obvious generalisation, and exactly what the single-network path
degenerates to -- gives different assigned values whenever more than one comparison has a flagged
pathway, and therefore a different row order. `pooled_degenerate` has flagged pathways in both
comparisons precisely so the two cannot be confused.

## The comparison branch's row names are uniquified with a bare suffix

`do.call(rbind, df)` over per-comparison frames that share row names produces `X`, `B`, `C`, `A`,
`D`, `X1`, `B1`, ... -- a bare `1`, with no separator. `make.unique` would give `X.1`. The row
names are part of the return value and `identical()` sees them, so the shim lets R's own `rbind`
produce them rather than reconstructing them.

Two smaller things in the same branch, both of which the shim initially got wrong:

* **The zero-dropping loop runs after the `rbind`,** in the tail shared by both modes. A pathway is
  dropped only if it contributes zero in *every* comparison, and it is dropped from all of them at
  once, changing `nrow` and leaving gaps in the row names.
* **`color.use` is reversed before `colors.text` is computed,** so the first entry of the reversed
  vector is what a `contribution.relative > 1 + tol` row is painted. Reversing afterwards swaps the
  two colours on every coloured row.

## `cellchatrs_upstream_env()` sources `visualization.R`, and the caller must have ggplot2 attached

`rankNet(mode = "comparison")` builds a grouped bar chart and calls `CellChat_theme_opts()`, which
is defined in `visualization.R`. Without that file the branch cannot even be evaluated, so the
reference environment sources it. It has `parent = globalenv()`, which means upstream's body
resolves its *unqualified* names -- `ggplot`, `theme`, `element_text` -- from the caller's search
path; a real CellChat session has `ggplot2` attached because CellChat imports it. The differential
gate has to attach it explicitly, or the *reference* fails with "could not find function" while the
port, calling `ggplot2::` explicitly, succeeds.

`identical()` on the returned `gg.obj` cannot hold across the two environments: `geom_bar`'s
`position` carries a *quosure* -- an expression plus its capture environment. The gate compares
`ggplot_build(gg)$data` instead, which is the rendered geometry, and that does match.

## A green differential gate can be vacuous, and this one was

`computeCommunProb` gated its Rust path on

```r
fallback <- nzchar(Sys.getenv("CELLCHATRS_FALLBACK", "0"))
```

which is **always `TRUE`**: the unset default is the string `"0"`, and `"0"` has four characters.
The escape hatch was permanently engaged, every `computeCommunProb` call was delegated to pinned
upstream, and the differential gate -- which compares the shim against upstream -- reported
`IDENTICAL: 0 failing comparisons out of 273` for a run that had compared upstream with itself.

Nothing in the output distinguished the two cases. The same run's `prob.*` and
`rshim.computeCommunProb` quantities were all measuring nothing, and `parity.json` recorded them at
a passing rung.

Two things now prevent it:

* the predicate is an explicit affirmative (`1`/`true`/`yes`/`on`), and
* `tests/parity/check_identical.R` **stops** if `CELLCHATRS_FALLBACK` is set, before any
  comparison runs, so the vacuous configuration cannot be produced silently.

`tests/parity/check_rust_path.R` is the complementary evidence, from outside the package: it points
`options$db` at a directory that does not exist and requires the Rust-only
`no CellChatDB export found` error. Upstream uses `object@DB` and never looks at a directory, so a
call that *succeeds* under that condition did not reach the kernel.

Fixing the predicate immediately exposed a second bug that the same permanent fallback had been
hiding: `as.character(pairLR.use$agonist)` is `character(0)` when the column is absent, and a
zero-length LR vector is a Rust panic (`lr_agonist has 0 entries but lr_ligand has 1`) rather than a
missing optional. `lr_col()` now normalises every LR metadata vector to `nLR` entries with `NA` as
`""` -- equivalent, because every use upstream makes of those columns is `!is.na(x) & x != ""`, and
required, because `NA_character_` cannot cross into a Rust `Vec<String>` at all.

## A fixture can be wrong in a way that reads as a kernel bug

Two of them, both found by the spatial gate and both worth writing down because the symptom points
somewhere else entirely.

**A hand-made `object@DB`.** Upstream resolves ligand/receptor pairs through `object@DB`; the Rust
kernel reads the database from a **directory on disk** (`cellchatrs_db_dir`). A spatial fixture
that supplies its own small `object@DB` therefore has the two sides expanding the L-R pair
differently, and reports `Prob[A,A]` of 0.95 against upstream's 0.41 -- a 2.3x "error" that is the
fixture's fault entirely. The fix is to build the spatial fixture from the *same* `data.signaling`,
`LRsig` and `DB` as the RNA configurations, which is what it now does.

**Order-dependent globals.** `check_identical.R` builds its fixtures from globals (`LRsig`,
`data.signaling`, `cell_group`) created by evaluating `gen_prob_golden.R`'s preamble. Later blocks
in the same file reassign `LRsig` in the *global* environment -- `computeCommunProbPathway` and
`subsetCommunication` each build their own small one -- so a block added later that reaches for
`LRsig` silently gets a 1-row, ligand-less frame from whichever block ran last. The symptom was
`nrow(pairLR.use) = 1, ligand = 0` and a kernel complaint naming neither the fixture nor the
clobber. The RNA fixture is now snapshotted under private names (`RNA_LRsig`, `RNA_data`,
`RNA_group`) immediately after the preamble.

A third, smaller one: `reinstall.sh` restores the previous install when `R CMD INSTALL` fails, and
a `... | tail -1` pipeline hides that exit status. A broken `R/modeling.R` therefore leaves a
*stale, working* package installed, and every subsequent measurement is of the previous build. The
symptom here was an error message that a fix had already addressed, which is the most expensive kind
of stale state. Install now parses the source first and its exit status is checked.

## Lazy evaluation hid a kernel error for a while

`computeCommunProb` passed the LR metadata vectors inline:

```r
res <- compute_commun_prob(..., lr_agonist = lr_col(pairLR.use, "agonist", nLR), ...)
```

R evaluates call arguments lazily, so that is a promise the callee forces at an arbitrary later
point, after the spatial branch has rebound `nLR1` and after `nLR` has been used for other things.
The kernel's refusal -- `lr_agonist has 1 entries but lr_ligand has 8` -- names neither the caller
nor the fixture, and the length it reports is the length at *force* time rather than at call time.
The vectors are now bound to locals first, which makes the values eager, named and inspectable, and
the lengths are checked in R with a message that says what was actually handed over.

## `computeCommunProb`'s spatial branch

`P.spatial` and `adj.contact` come from `computeRegionDistance`, which is the Rust kernel with the
exact k-d tree, so the branch no longer falls back to upstream. Upstream's own path would resolve
the neighbours through `AnnoyParam` -- the approximation this port exists to remove -- so falling
back was not merely incomplete, it was the wrong answer.

The R-side transformation stays in R rather than moving into the kernel because it is mostly
*messages*: two `print`s whose text embeds `Sys.time()`, four `cat`s, and a `stop()` whose message
embeds a computed `format(1/d.min, digits = 2)`. `nLR1` -- the index past which
`P.spatial * adj.contact` is applied -- depends on the `annotation` column's contents and is what
distinguishes contact-dependent from diffusible L-R pairs, so the four branches and their exact
texts are reproduced.

## Invariants asserted on every kernel run

Cheap, and they catch a lot:

* `all(Prob >= 0)` and `all(Prob <= max(Prob))`; all finite unless inputs are non-finite;
* `all(Pval >= 0 & Pval <= 1)` and `all(Pval * nboot == round(Pval * nboot))`;
* `Pval[Prob == 0] == 1`;
* `dimnames(Prob) == dimnames(Pval)`; last dimension `== nrow(LRsig)`;
* `dimnames(Prob)[[1]] == levels(object@idents)`;
* **`Prob` is invariant to `nboot`** — the observed network does not depend on the
  number of permutations. A very strong end-to-end check.

## Metamorphic relations

| transform | expected |
|---|---|
| permute cells | `Prob`, `Pval` unchanged (given the same seed) |
| relabel / reorder cell groups | entries permute identically; `dimnames` follow |
| permute genes | unchanged |
| scale `data` by a constant | unchanged (`data/max(data)` is scale-free) |
| duplicate every cell | unchanged (per-group means are duplicated but order statistics scale) |
| `nboot` change | `Prob` unchanged, `Pval` in `{k/nboot}` |

## FMA, reassociation, and whether this test suite would notice

The contract is that FMA contraction is off and nothing is reassociated, because both change
rounding. Three separate claims, three separate pieces of evidence — recorded here because the
first two are the ones that are easy to assert and hard to actually check.

**1. The suite is genuinely sensitive to rounding.** Measured, by a control experiment: replace
`F80::mul` with a round trip through `f64`, i.e. simulate a 53-bit accumulator, and run the suite.
**9 tests fail across 6 test binaries** — the `r-core` unit tests, the x87 oracle, the probability
goldens, the expression goldens and the database fuzz suite. One truncation in one function is
caught by five independent layers. Without this control, "311 tests pass" says nothing about
whether the tests can see rounding at all.

**2. Codegen variants do not move the numbers.** `scripts/codegen_variants.sh` builds the
parity-critical tests three times and diffs the results:

```text
default      -C target-feature=-fma      311/311
-fma off     -C target-feature=+fma      311/311   byte-identical
host native  -C target-cpu=native         311/311   byte-identical
```

This is expected rather than lucky. **Rust does not contract floating-point expressions**: rustc
does not set LLVM's `contract` fast-math flag on `fmul`/`fadd`, so there is no contraction for a
target feature to enable. The flag in `.cargo/config.toml` is belt-and-braces on top of a language
guarantee, and this file's job is to keep it that way rather than to credit it with a guarantee the
language already provides. What it does still protect against is a `-ffast-math`-style flag
appearing in `RUSTFLAGS` somewhere, which the same script would catch.

The in-crate probes are in `crates/r-core/tests/no_fma.rs`, and they are written so they cannot be
mistaken for flag checks: `a_bare_multiply_add_in_this_crate_rounds_twice` states an observable
invariant, and passes with `+fma` too, because the language — not the flag — is what forbids the
contraction.

**3. Reassociation is observable, and is therefore testable.** This one needed work, because the
first oracle was not sharp enough to prove anything.

R's `real_mean` pass 1 is repeated addition of a whole vector in 80-bit, so the order is part of
the specification. The `chain` rows in `tests/fixtures/f80_ref.txt` record R's left-to-right
accumulation with the inputs written out. Reassociating them changes the answer in only **2 of 156**
rows — because `F80` has 11 more mantissa bits than `f64` and absorbs most reordering. That is a
genuine and reassuring property, and it is also why those rows cannot prove the order is right: the
order is barely *visible* there, so matching proves little.

So the oracle grew a `chain_cancel` class: interleaved `+B, +s, -B, +s` over
`B ∈ {1e16, 1e17, 1e300, 1e-300, 2^53}` and `s ∈ {1, 0.5, 3, 1e-8, 7}`, lengths 8–48. Left to right,
each `s` is swallowed by the `B` in front of it; a tree reduction recovers values R's order never
saw. **88 of 275 rows (32%)** now distinguish the two orders, and
`adversarial_chains_pin_the_order_hard` asserts every one of the 275 matches x87 while requiring at
least a quarter to be order-sensitive. The rate is structural, not a function of length — it did not
move when the oracle was extended from 5 to 11 lengths per `(B, s)` pair. What it depends on is
whether `B` exceeds what 11 extra mantissa bits can absorb.

The same widening made the `F80` probe meaningful. `F80::mul` keeps 64 bits, so `F80(a)·F80(b) + F80(c)`
*is* the fused result to `f64` precision — and that is the point, not a hazard. The protection is
that the accumulation is performed in `F80` explicitly rather than left to whatever the compiler
fuses out of an `f64` expression, so the rounding is a documented 64-bit-mantissa rule instead of a
codegen accident. `the_f80_accumulator_behaves_like_x87_not_like_f64` asserts it matches the fused
value *and* differs from the two-rounding `f64` one; if it ever agreed with the latter, `F80::mul`
had begun truncating to 53 bits and `triMean` would be quietly wrong.

Oracle size after the change: 9,576 rows — 3,000 add, 3,000 sub, 3,000 mul, 120 prod chains, 455
sum chains (180 benign + 275 adversarial). Regeneration is byte-identical for the pre-existing
9,301 rows, so the fixture grew rather than moved.

## Three factors in the melted table, not two

`subsetCommunication` returns `source`, `target` **and** `interaction_name` as factors. Everything
else in the table -- `ligand`, `receptor`, `annotation`, `evidence` -- is character.

The mechanism is that upstream melts the `K x K x N` array and calls `var.convert` on the result,
which factors every column that came from a **dimname**. `source` and `target` come from the first
two dimensions; `interaction_name` comes from the third. The join columns arrive from `LR` with the
melt, are not dimnames, and stay character.

This was a real defect in the port, not a documentation error. The shim returned
`interaction_name` as character, and **every value comparison in the golden corpus still passed**:
row count, column names, column order, `prob`, `pval`, the `source|target` key order, and the level
order of `source` and `target`. The corpus recorded `col_order` and `classes` for nothing -- it had
no `classes` record at all. What caught it was `tests/parity/check_merge.R`, which compares whole
data frames and therefore sees `class()`. The comment in `R/modeling.R` that asserted "only `source`
and `target` are factors upstream ... everything else stays character" was confidently wrong, and it
had survived because nothing compared `class()`.

Three changes, so it cannot recur quietly:

1. `tests/fixtures/netfiltered_golden.txt` records `col_order`, `classes`, `interaction_levels` and
   `df_source_levels` per case, captured **before** the generator grafts `source_target` on -- that
   column belongs to `aggregateNet`'s filtered branch, not to the returned table.
2. `NetTable` carries `interaction_levels`, and two new Rust tests pin the levels and the class of
   every column, so `cargo test` catches a regression rather than needing the R gate.
3. The shim factors `interaction_name` over those levels.

The levels are the **whole** `Prob` third dimnames, not the subset that survived the threshold. An
interaction filtered out of the table still contributes an unused level, so the factor is not
`droplevels(rownames(df))` and reproducing the unused levels is part of the contract.

## `mergeCellChat` is gated, and the gate is not a tautology

`mergeCellChat` is not reimplemented. It lives in `CellChat_class.R` (not `analysis.R`, as the
objective's file list has it), it contains no arithmetic, and the 14.1 scope decision leaves S4 slot
assembly in R. `tests/parity/check_merge.R` gates it, and the gate is deliberately **not** a
comparison of `mergeCellChat` with upstream's `mergeCellChat` -- there is only one implementation, so
that comparison would report `IDENTICAL` while testing nothing. It is the failure mode recorded above
for the `rankNetPairwise` fallback, and it is worth refusing by name.

The gate does two things instead.

**The merge's own contract**, which has to be written down somewhere because there is no second
implementation to compare against: `data.joint` is the `cbind` in argument order; `genes.use` is the
intersection **in A's order**; `meta.use` is the intersection of `meta` column names; `idents` on a
merged object is a *list* of per-dataset factors plus a `$joint` factor, and `$joint`'s levels are
`union(levels(a), levels(b))` -- first-appearance order, not sorted; `add.names` propagates to
`meta$datasets` and to the names of `net`, `netP`, `LR` and `idents`; and both `stop()` messages are
reproduced verbatim, double exclamation mark included.

**What consumes a merged object.** A merged object feeds the comparison analyses, not
`computeCommunProb` -- `@idents` is a list, and the kernel needs a factor, so re-running the kernel
on a merged object is not something upstream supports either. What the port has to get right is
whether the Rust `aggregateNet` / `subsetCommunication` / `filterCommunication` / `rankNet` accept
the net structures a merge carries. Each dataset's net is lifted back out of the merged object and
compared against upstream on the same net. That is non-vacuous: a merge that renumbered a group,
dropped a dimname or nested `net` at the wrong depth changes these results. And the test asserts the
level order is load-bearing by reversing it and requiring the answer to change.

Three upstream behaviours found while writing it, each pinned because it is silent:

- **`data.signaling` can come out empty.** `data.joint` is intersected to `genes.use`, but the
  subsetting of `data.signaling` is against `gene.signaling.joint`, which is the **union** of the two
  objects' gene names. So every row the intersection kept that only one dataset had is dropped again,
  and when the gene sets are disjoint the result is a `data.signaling` with **zero rows** and no
  colnames, from a merge that reported success.
- **`aggregateNet`'s two branches write different things.** Called with no filter arguments it leaves
  `prob` as the `K x K x N` array and adds `count` and `weight` as `K x K` matrices. Given
  `sources.use`, `targets.use`, `signaling` or `pairLR.use` it melts the table with `dplyr` and
  `tapply`s over the `source|target` character key. Neither writes a table into the object --
  `subsetCommunication` *returns* the table. So "run `aggregateNet`, then `rankNet`" cannot work, and
  `rankNet` reads the array directly.
- **`rankNet`'s `slot.name` defaults to `"netP"`.** On a single-dataset view of a merged object
  `@netP` is empty, so `rankNet` stops with "No inferred communications for the input!" regardless of
  the threshold. That message means "wrong slot", not "empty network". Relatedly,
  `subsetCommunication` tests `pval >= thresh` and `rankNet` tests `pval > thresh`; the same argument
  name, two conventions, and with `nboot = 5` the attainable p-values are {0, 0.2, 0.4, 0.6, 0.8, 1}
  so most "reasonable" thresholds are either vacuous or total.

## `seed.use` changes `Pval` and not `Prob`, and the reason is an asymmetry in the bootstrap

`computeCommunProb` bootstraps like this:

```r
set.seed(seed.use)
permutation <- replicate(nboot, sample.int(nC, size = nC))
data.use.avg.boot <- lapply(1:nboot, function(nE) {
  groupboot <- group[permutation[, nE]]
  aggregate(t(data.use), list(groupboot), FUN = FunMean)     # <- data.use is NOT permuted
})
```

**The labels are permuted and the data is not.** Each draw pairs cell `i` with a *different* cell's
group label, so the per-group means really do differ from draw to draw and the null distribution
`Pboot` depends on `seed.use`. But `Prob` is `Pnull`, computed from the *unpermuted* `data.use.avg`,
and is therefore **independent of `seed.use` for every configuration** — not as a consequence of
`population.size`, which is what an earlier draft of this file claimed, and which was wrong.

Measured on the standard fixture over 120 seeds: `Prob` is identical at all 120, while the `Pval`
vector takes 120 distinct values. So the seed has exactly one job, and it is the one that decides
whether a `Prob` is called significant.

This is worth writing down because it is easy to get backwards, and because it makes the obvious
version of a seed-stream test vacuous. `tests/parity/stat_equiv.R` was first written to watch `Prob`
for seed sensitivity, reported the seed as inert, and the control it carries is what caught it. The
gate now watches `Pval` where the seed acts, and **asserts** the `Prob` invariance rather than being
confused by it — a port that let the permutation leak into `Pnull` would fail on that assertion.

## The independent stream, and why Park-Miller rather than a 64-bit xorshift

`tests/parity/stat_equiv.R` draws its 120 seeds from a generator the test implements itself, salted
by the fixture's shape, so that neither the port nor upstream can influence the sequence. Without
that, a port and a reference sharing an RNG could satisfy the gate by construction — which is the
failure mode the repo already records for the `rankNetPairwise` fallback.

The first choice was a 64-bit xorshift, and it cannot be written in R at all: `bitwShiftL` is defined
on R *integers*, so shifting a word with the high bit set returns `NA` with "NAs introduced by
coercion to integer range", which propagates into `set.seed` as "supplied seed is not a valid
integer". Park-Miller (`s <- (s * 48271) %% 2147483647`) is exact in doubles — the product stays below
`2^47` — and reproducible forever. It is a weak generator, and that does not matter: the requirement
is that the seed sequence be *outside both implementations' control* and exactly reproducible, not
that it be statistically excellent. Rejecting a perfect stream in favour of a merely adequate one
would buy nothing.

One more trap worth recording: `<<-` inside a second `local()` block does **not** reach a variable
created in the first. `<<-` walks the enclosing scopes, and a `local()` frame is not in scope for a
sibling, so the two halves have to be created in one `local()` and exported with `assign`.

The gate reports `EQUIVALENT: 0 failing checks out of 130`, with all 120 seeds producing distinct
null distributions on both sides.

## The standalone CLI, and the four bugs its parity gate found

The objective requires `r-core` to ship "as a standalone library/CLI so the numerics are usable
without R". `crates/cellchatrs-cli` is that CLI: a `cellchatrs` binary that reads a self-describing
input file and writes the `Prob`/`Pval` networks, the averaged expression and `aggregateNet`'s two
matrices, with no R in the process at all.

`tests/parity/check_cli.R` runs the binary against **pinned upstream**, not against the shim, over 14
configurations — all four `type.mean` values (two supplied as `match.arg` prefixes), `population.size`
both ways, `nboot` from 1 to 9, `Kh` across six orders of magnitude, `n` at 1 and 2, `raw.use` both
ways, and two seeds. `Prob`, `Pval`, `computeAveExpr`'s means and `aggregateNet`'s two matrices are all
required `identical()`. The reference is a sourced copy of upstream's `modeling.R`, not this package's
shim, because the shim shares the kernel and comparing it against itself would prove nothing.
**IDENTICAL: 0 failing comparisons out of 231.**

The wire format uses `%a` hex floats and is written **gene-major, one line per gene**, with the reader
indexing `data[gene * n_cells + cell]`. Both are deliberate: hex makes a mis-transcribed digit a
*different double* rather than a near-miss, so a parity failure points at the kernel rather than at the
file; and the explicit transposition is the one `expr_parity.rs` records going wrong silently, where the
gene *set* stays correct while the *values* permute.

Four defects, three of them in the CLI itself. The gate is worth the trouble.

**`raw.use = FALSE` is a matrix choice and nothing else.** Upstream is

```r
if (raw.use) data <- as.matrix(object@data.signaling)
else          data <- as.matrix(object@data.smooth)
...
data.use <- data/max(data)          # modelling.R:111, unconditional
```

The CLI divided by the row mean when `raw.use` was `FALSE`, on the reasonable-sounding guess that
"raw versus projected" implied a normalisation. It does not, and there is no library-size
normalisation anywhere in `computeCommunProb`. Upstream produced a full `Prob` array; the CLI raised
`missing value where TRUE/FALSE needed`, because `data.smooth` is row-normalised, so `data / 1` fed the
Hill function values orders of magnitude above `Kh` and `P1` saturated to 1 everywhere. The R shim never
had this bug — it copies upstream's two lines verbatim — so the gate that caught it exists only for the
CLI, which is the argument for having one.

**`computeAveExpr` and `computeCommunProb` read different matrices.** `computeAveExpr` has no
`raw.use` parameter and defaults to `slot.name = "data.signaling"`. So for a `raw.use = FALSE` object
the two functions read `data.smooth` and `data.signaling` respectively, and that is upstream. The CLI
reported the *kernel's* observed aggregate under the name `ave_expr`, which agrees with
`computeAveExpr` for thirteen of the fourteen configurations and differs only for `raw.use = FALSE` —
where the gate said `computeAveExpr means MISMATCH 96 values`. An intermediate "fix" that made
`compute_ave_expr` honour `raw.use` was also wrong, and the gate rejected it too. The output now
carries both blocks, `ave_expr` and `kernel_ave`, correctly named.

**`Prob` is not bounded by 1.** An invariant check written here — `all(Prob >= 0 & Prob <= 1)` — failed
against an output the port reproduces bit for bit: **16 of 96 values exceed 1, the largest 2025.1**.
`Prob` is `P1 * P2 * P3 * P4 * P.spatial`; while `P1 = dataLR^n / (Kh^n + dataLR^n)` is bounded by 1,
the agonist and antagonist terms are `crossprod(matrix(x, nrow = 1))` — plain sums of expression values
— and `P4` is a group-proportion outer product, so any of them can exceed 1 and the product can. `Prob` is
a *score*, not a probability, and upstream never clamps it. Non-negativity *is* worth asserting: a
negative `Prob` would mean an expression path had a sign error, and `prob > 0` downstream would silently
drop the interaction. The objective's own recorded invariants do not include a range check on `Prob`, and
that turns out to be right.

**A generator that shadows a name upstream uses stops being a reference.** The fixture generator defines
`aggregate_branch` — it was originally called `aggregate` — and every `computeCommunProb` call failed
with `unused argument (FUN = FunMean)`. The sourced upstream environment has `parent = globalenv()`,
which is load-bearing for `Matrix::rowSums` (see the comment on `cellchatrs_upstream_env`), and it means
a function the *generator* defines in the global environment shadows the same name inside upstream's own
body. Upstream's bootstrap calls `aggregate(t(data.use), list(groupboot), FUN = FunMean)` and was reaching
this script's function instead of `base::aggregate`. Every golden recorded that error and the run reported
success. The rule for any script that sources upstream into `globalenv()` is the one the shim already
follows: never define a name upstream uses.

## `computeAveExpr` cannot be asked for `thresholdedMean`

Its `match.arg` offers three choices — `triMean`, `truncatedMean`, `median` — while
`computeCommunProb` offers four. `computeAveExpr(object, type = "thresholdedMean")` raises

```text
'arg' should be one of "triMean", "truncatedMean", "median"
```

so the `match.arg` error text is not shared between the two functions even though both declare
`type.mean`. The generator records the `ave_expr` block only for the types that can answer for it, rather
than working around the asymmetry, and the comment in `pipeline.rs` says why.

## `filterCommunication` with `NA` in `prob`: three adjustments, not one

The tutorial's real object carries `NA_real_` in `net$prob`, and upstream tolerates it end to end.
The port refused it in two places and miscounted it in a third, and the three had to be fixed
together because each one masks the next: fixing only the binding conversion moved the failure from
argument conversion to the counts, and fixing only the counts moved it to the messages.

**1. The binding rejected `NA` before the body ran.** `extendr` 0.9's generated `Vec<f64>`
conversion raises `Error::MustNotBeNA` during argument conversion -- naming neither the argument
nor the function -- so `filter_communication(prob = as.numeric(net$prob), ...)` died on any `NA`.
It now takes `Robj` and reinterprets the `REALSXP` payload with `as_real_slice()`, which performs
no missingness check and preserves the `NA_real_` bits exactly.

**2. The three `sum(net$prob > 0)` counts come back `NA`.** R's `sum()` uses `na.rm = FALSE` unless
told otherwise, so any `NA` (or `NaN` -- `NaN > 0` is `NA` in R, not `FALSE` as in IEEE) in the
buffer makes the count `NA`, and the `FilterResult` fields are `Option<usize>` for exactly that
reason. They cross back as `NA_integer_`, and R-side arithmetic on them propagates `NA` the way
upstream's does. The pre-fix code counted `v > 0.0` over everything, which is correct if and only
if the buffer is finite -- and no fixture covered anything else, because no fixture's `prob` held
an `NA`.

**3. `NA`-summing slices are dropped from `lr_nonzero`.** `which(apply(net$prob, 3, sum) != 0)`
evaluates to `NA` for such a slice and `which()` drops `NA`s. The pre-fix code accumulated every
slice through `F80`, so an `NA`-containing slice produced `F80::NAN`, compared `!= 0.0`, and was
*included* -- the exact opposite of R. `NaN`-containing slices stay included (`NaN != 0` is TRUE in
both languages), which `slice_sum` handles by checking the `NA_real_` bits *before* converting and
letting every other non-finite value through the existing `F80` path.

**4. The messages say "NA", not "NaN%".** `scales::percent(NA)` is `NA` (verified: `percent(NaN)` is
`NA` too), and the caller pastes it, so upstream prints e.g. "NA interactions are removed!".
`cellchatrs_percent` mapped missingness to `"NaN%"`, matching neither input flavor. Both now map
to `"NA"`.

None of this changed any existing fixture: the corpus has no `NA`/`NaN` in `prob`, so every count
stays `Some` and every slice stays included, which the full suite confirms. The tutorial's stage 5
is the end-to-end proof, and `filter::na_tests` pins the semantics in Rust.

A second `NA`, in a different argument of the same call, surfaced in the same run and is fixed in
the same change for the opposite reason: `DB$geneInfo$Symbol` holds exactly one `NA` in the pinned
human DB (1 of 26,827, audited), and `symbols: Vec<String>` rejected it the same way. Upstream never
notices it -- `extractGeneSubset` only asks `%in%`/`intersect()` with non-`NA` needles, and R matches
an `NA` table entry against nothing but an `NA` needle -- so `NA` entries are dropped at the binding
via `CanBeNA::is_na`, which distinguishes the `NA_string_` sentinel from a literal `"NA"` without
unsafe code. Inserting a placeholder (`""`, `"NA"`) would have been wrong: it could match a lookup
R would not match. (`Database::is_symbol` already documented this equivalence; the binding is what
could not express it.)

`data` keeps the loud rejection on purpose: `NA` in the expression matrix is not modeled anywhere
downstream (R's `max()` returns `NA` where the binding folds with `f64::max`, which *ignores* NaN
payloads), so refusing loudly beats diverging silently. The asymmetry is the contract: `prob` NA is
specified end-to-end, `data` NA is refused at the door.

## `subsetCommunication(slot.name = "netP")` emitted one row per input row

Upstream aggregates to one row per `(source, target, pathway)` group (`summarize(prob = sum(prob))`,
`summarize(pval = mean(pval))`); the port's `aggregate_netp` walked every *row* index in sorted
order and emitted a row for each, so a table with duplicate group keys came back unaggregated -- 1042
rows for 576 distinct keys on the tutorial's real `netP`, where the pathway interaction of two
L-R pairs between the same groups is routine. Fixed by emitting the first index per distinct key;
`subset::netp_tests` pins a duplicate-key table and the sorted-group order. No corpus case has a
duplicate key, so the suite stayed green throughout -- which is itself the lesson, recorded here
rather than as an excuse: the corpus proves the aggregation arithmetic, not that aggregation
happens, and only a table with a repeated group can tell those apart.

## Centrality: what is ported, what is delegated, and why the split is where it is

Upstream's `computeCentralityLocal` computes eleven measures per pathway network. Seven are pure
functions of the input and live in `r-core::centrality`; four stay in R and call igraph directly.
The split is not a judgment about difficulty -- it is forced by measurement, and the measurement is
gated (`tests/parity/check_centrality.R` asserts it) so the justification cannot rot while the code
changes around it.

**`hub_score`, `authority_score`, `eigen_centrality` and `page_rank` disagree with themselves.**
Run each twice on identical input with no seed set and at least hub/eigen/authority come back
different -- HITS power iteration, ARPACK and PRPACK all draw start vectors or stop on tolerances
that are not reproducible across runs. (`page_rank` happened to agree on the probe graphs; it stays
delegated anyway, because PRPACK is version-sensitive iterative code and a port would trade
identical-by-construction for a reimplementation that breaks the next time igraph retunes a
tolerance.) Bit-parity with a nondeterministic oracle is meaningless: a golden would pin one draw,
and a passing test would prove the port equals *a* draw rather than *the* computation. So those four
are the same package, the same function, the same input -- identical by construction -- and the gate
seeds both sides identically before comparing the assembled list, so any remaining difference is in
what the shim hands igraph or in the Rust measures, which is what is under test.

**Strength is plain sequential `f64` accumulation in edge-ID order -- not R's `rowSums`.**
`strength_all` iterates edges `0..E-1` accumulating `res[from] += weight`. Edge IDs follow
column-major creation (`j` outer, `i` inner in `igraph_i_weighted_adjacency_directed`), so per
vertex this is ascending-index order -- but with `f64` `+=`, not R's long-double `sum()`. Measured
4669/4669 against installed igraph 2.3.4 for sequential-`f64`, versus 3802/4669 for `rowSums`. The
difference is precision, not order, and it is exactly the kind of difference that looks like nothing
in a `digits = 6` printout and fails `identical()` every time. Unweighted degrees are integer counts
of strictly-positive entries (`rowSums(net > 0)` -- a negative entry *creates an edge* via
`M != 0.0` but is *not counted*, and the corpus caught `!= 0.0` meaning the wrong thing).

**Betweenness is Dijkstra with dist-plus-one encoding, exact 2-way-heap tie rules, epsilon
comparisons at 1e-10, and Brandes accumulation -- all replicated operation by operation.**
The load-bearing details, each verified rather than assumed: distances stored plus one (`0.0` means
infinity); the heap negates keys and `shift_up` swaps on `>=` (ties rise) while `sink` swaps only
on strict `<` preferring left (ties stay) -- any other tie rule changes pop order, which changes
rounding of every accumulated score without changing any parent *set*; `igraph_cmp_epsilon` with
its zero/subnormal/overflow branches taken in order; path counts (`nrgeo`) as `f64`; loops excluded
from traversal but included once in strength; normalisation factor exactly `1.0`, skipped as an
exact no-op. The golden (`centrality_golden.txt`, 77 cases: 60 adversarial random graphs, reciprocal
edge cases, 2 error cases, 6 real tutorial slices) is compared by parsed bits with a NaN wildcard,
because R's `%.17g` shortest-spelling and Rust's `{:.17e}` spell the same doubles differently.

**The installed igraph differs from main-branch source, so empirical behavior rules.**
The messages in 2.3.4 (`Weight vector must be positive. Invalid value`, `...must not contain NaN
values...`, `Some weights are smaller than epsilon...`, each with a `\nSource:
centrality/betweenness.c:NNN` suffix that is part of `conditionMessage()`) match neither the
main-branch strings nor their check order: NaN is checked *first* regardless of position (a vector
with both negative and NaN entries reports NaN), positivity second, and the epsilon warning fires
at `minweight <= 1e-10` (warns at exactly 1e-10, silent at 1.0000001e-10). The port implements the
measured behavior and the corpus pins the texts byte for byte; the `Source:` suffix is reproduced
because the gates compare `conditionMessage()`.
