# Golden corpus for `filterCommunication`, from the pinned upstream commit.
#
# The function has two independent halves, and the corpus needs both:
#
#   1. `min.cells`: zero `net$prob` on the rows/columns of cell groups with too few cells.
#      Always active. Emits a `cat()` whose percentage comes from `scales::percent(x,
#      accuracy = .1)`, so the *formatting* is part of the contract.
#   2. `min.samples >= 2`: recompute per-sample group means, score every L-R pair per
#      sample as an outer product, binarise, count the samples in which each
#      (source, target) is present, and keep only those present in >= `min.samples`.
#      Needs >= 2 samples; errors above the sample count.
#
# R-isms this stresses, all documented in src/rust/crates/r-core/src/filter.rs:
#   * `apply(score.LR[, , jj, ], c(1, 2), sum)` sums over the *sample* axis, in sample
#     order, and R's `sum` widens to LONG_DOUBLE -- though the values are 0/1, so the
#     interesting part is the axis, not the precision;
#   * `droplevels(group.use)` inside the sample loop changes the *group order* per sample,
#     so the zero-fill back to `levels(group)` has to be per sample, not once;
#   * `if ("data.smooth" %in% methods::slotNames(object) == FALSE)` -- `%in%` binds tighter
#     than `==`, so this is `(x %in% y) == FALSE`, which is what was meant;
#   * `cell.excludes.sample <- setdiff(cell.excludes.sample, cell.excludes)` at the end, so a
#     group excluded for having too few cells *overall* is not also treated as "rare";
#   * the two `cat()`s use `cat(..., '\t')` and `cat(paste0(..., '\n'))`, so the separator is
#     a literal tab and the messages are concatenated without a newline between them.
#
# Usage:  R_LIBS=.rlib R --vanilla -f tests/parity/gen_filter_golden.R
suppressWarnings(suppressMessages({library(Matrix); library(collapse); library(dplyr)}))

CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
env <- new.env()
for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R")) {
  sys.source(file.path(CC, "R", f), envir = env, keep.source = FALSE)
}

## The slots `filterCommunication` reads: `net`, `idents`, `meta$samples`,
## `options$parameter$raw.use` / `$type.mean` / `$trim`, `data.signaling`, and
## `DB$interaction` / `$complex` / `$geneInfo`.
setClass("MiniFilter", representation(
  net = "list", idents = "factor", meta = "list", options = "list",
  data.signaling = "ANY", DB = "list", data.smooth = "ANY"
))

fmt <- function(v) paste(sprintf("%.17g", as.numeric(v)), collapse = ",")
fmti <- function(v) paste(as.integer(v), collapse = ",")
## Write to a scratch name and rename at the end. Evaluating a generator's *preamble* is
## the tempting way to reuse its fixture definitions -- and the preamble contains the
## `file(..., "wt")` line, which truncates the committed corpus to zero bytes. That has
## happened here twice, once silently. A `.part` file makes the mistake harmless.
q <- file("tests/fixtures/filter_golden.txt.part", "wt")

## ------------------------------------------------------------------ the object factory
## A small CellChatDB-like table. Complexes and cofactors are what `computeExpr_LR` needs to
## split a multi-subunit receptor, so they are non-trivial on purpose.
mk_db <- function() {
  ## Mirrors the *real* CellChatDB schema, which differs from the obvious one in three ways
  ## that all matter here:
  ##   * `complex` has columns `subunit_1 .. subunit_5` (not a `subunit` column), and its
  ##     **row names** are the complex names -- `extractGeneSubset` does
  ##     `complex_input[match(complex, rownames(complex_input), nomatch = 0), ]`;
  ##   * `geneInfo` has a `Symbol` column (not `gene`);
  ##   * `interaction$receptor` holds the *complex name* for a multi-subunit receptor, e.g.
  ##     "TGFbR1_R2" / "FZD8_LRP5" -- not a subunit and not a separate complex column.
  complex <- data.frame(
    subunit_1 = c("R1a", "R2", "R3", "R4", "R3", "R4"),
    subunit_2 = c("R1b", "", "", "", "", ""),
    subunit_3 = rep("", 6), subunit_4 = rep("", 6), subunit_5 = rep("", 6),
    stringsAsFactors = FALSE
  )
  rownames(complex) <- c("R1", "R2", "R3", "R4", "R5", "R6")
  interaction <- data.frame(
    interaction_name = c("L1_R1", "L2_R2", "L3_R3", "L4_R4", "L5_R5", "L6_R6"),
    ligand = c("L1", "L2", "L3", "L4", "L5", "L6"),
    receptor = c("R1", "R2", "R3", "R4", "R5", "R6"),
    stringsAsFactors = FALSE
  )
  cofactor <- data.frame(
    cofactor_1 = character(0), cofactor_2 = character(0), stringsAsFactors = FALSE
  )
  ## `Symbol` lists the *subunits* and the ligands, and deliberately **not** the complex names
  ## R1/R5/R6. `extractGeneSubset` has
  ##   `complex <- geneSet[which(geneSet %in% geneIfo$Symbol == "FALSE")]`
  ## and `%in%` binds tighter than `==`, so the `==` coerces a logical to character and the
  ## test really asks "is this name *absent* from `Symbol`?". For a well-formed database that
  ## selects exactly the complex names, so the bug is benign -- but only because the complex
  ## names are absent. Put them in `Symbol` and complexes silently stop expanding.
  geneInfo <- data.frame(
    Symbol = c(paste0("L", 1:6), "R1a", "R1b", paste0("R", 2:4), paste0("X", 1:4)),
    isComplex = c(rep("FALSE", 11), rep("TRUE", 4)),
    stringsAsFactors = FALSE
  )
  ## 6 ligands + R1a, R1b + R2, R3, R4 = 11 single genes, plus 4 cofactors = 15. The
  ## complex names R1, R5 and R6 are absent on purpose (see the note above).
  stopifnot(nrow(geneInfo) == 15L, nrow(complex) == 6L,
            !any(c("R1", "R5", "R6") %in% geneInfo$Symbol))
  stopifnot(identical(rownames(complex), as.character(interaction$receptor)))
  ## The non-empty subunits, as a *set*: R1a, R1b (complex R1), R2, R3, R4, plus R3 (R5)
  ## and R4 (R6). Order is not asserted because `unlist()` on a data.frame goes column by
  ## column, not row by row.
  stopifnot(setequal(unlist(complex)[nzchar(unlist(complex))],
                    c("R1a", "R1b", "R2", "R3", "R4")))
  list(interaction = interaction, complex = complex, cofactor = cofactor, geneInfo = geneInfo)
}
mk <- function(net_prob, net_pval, levels, lr, idents, samples, data, db, type_mean = "triMean",
               trim = 0.1, raw_use = TRUE, with_smooth = FALSE, smooth = NULL) {
  dim(net_prob) <- c(length(levels), length(levels), length(lr))
  dimnames(net_prob) <- list(levels, levels, lr)
  dim(net_pval) <- dim(net_prob)
  dimnames(net_pval) <- dimnames(net_prob)
  o <- new("MiniFilter",
           net = list(prob = net_prob, pval = net_pval, prob.dim = dim(net_prob),
                      dimnames = dimnames(net_prob)),
           idents = idents, meta = list(samples = samples),
           options = list(parameter = list(raw.use = raw_use, type.mean = type_mean,
                                           trim = trim), mode = "single"),
           data.signaling = data, DB = db)
  ## `smooth` defaults to `data`, which makes `raw.use = FALSE` a no-op and the case worthless:
  ## the two settings have to read *different* matrices, or the branch is never exercised.
  if (with_smooth) o@data.smooth <- if (is.null(smooth)) data else smooth
  o
}

## Deterministic integer-ish expression so the outer products and the group means are exactly
## reproducible and free of RNG.
mk_data <- function(genes, cells, seed) {
  set.seed(seed)
  m <- matrix(rpois(length(genes) * cells, 3) + 1, nrow = length(genes),
              dimnames = list(genes, NULL))
  m
}

## ------------------------------------------------------------------- the case list
cases <- list()
genes <- c(paste0("L", 1:6), "R1a", "R1b", paste0("R", 2:4), paste0("X", 1:4))
stopifnot(length(genes) == 15L)
DB <- mk_db()
## The L-R names in `net$prob`'s third dimnames must be the DB's `interaction_name`s, because
## upstream does `idx <- match(LR.nonzero, interaction_input$interaction_name)` and then
## `interaction_input$ligand[idx]` -- an unmatched name is `NA`, and the failure surfaces
## much later as "subscript out of bounds" inside `computeExpr_LR`, not as a match error.
LR_ALL <- as.character(DB$interaction$interaction_name)

## 1. Single sample, `min.samples = NULL`: only the `min.cells` half runs. 3 groups of 10, 10
##    and 4 cells, so group 3 is excluded at `min.cells = 10`.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 10), rep("g2", 10), rep("g3", 4)), levels = lev)
  samples <- factor(rep("s1", 24), levels = "s1")
  prob <- array(rep(c(0.5, 0.2, 0, 0.3, 0.1, 0.4), each = 9), dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "one_sample_mincells", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 24, min_cells = 10L, min_samples = NULL,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes)
}

## 2. Same, with `min.cells` excluding nothing (3 groups of 8, `min.cells = 3`): the
##    `else` branch, so no percentage is printed at all.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 8), rep("g2", 8), rep("g3", 8)), levels = lev)
  samples <- factor(rep("s1", 24), levels = "s1")
  prob <- array(rep(c(0.5, 0.2, 0, 0.3, 0.1, 0.4), each = 9), dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "one_sample_keep", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 24, min_cells = 3L, min_samples = NULL,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes)
}

## 3. Two samples, `min.samples = 2`: the cross-sample half runs, with every group present in
##    both samples.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 8), rep("g2", 8), rep("g3", 8)), levels = lev)
  samples <- factor(c(rep("s1", 12), rep("s2", 12)), levels = c("s1", "s2"))
  prob <- array(rep(c(0.5, 0.2, 0, 0.3, 0.1, 0.4), each = 9), dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "two_samples_ms2", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 24, min_cells = 3L, min_samples = 2L,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes)
}

## 4. Three samples, `min.samples = 2`, and one L-R pair absent from a whole sample: the
##    interesting case, where the consistency filter must actually remove something.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 8), rep("g2", 8), rep("g3", 8)), levels = lev)
  samples <- factor(c(rep("s1", 8), rep("s2", 8), rep("s3", 8)), levels = c("s1", "s2", "s3"))
  prob <- array(rep(c(0.5, 0.2, 0, 0.3, 0.1, 0.4), each = 9), dim = c(3, 3, 6))
  ## Make the 4th L-R pair zero everywhere, so it drops out of `LR.nonzero` and never reaches
  ## the per-sample scoring.
  prob[, , 4] <- 0
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "three_samples_ms2", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 24, min_cells = 3L, min_samples = 2L,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes)
}

## 5. As (4) but `rare.keep = TRUE`, with a group that is present overall but *rare* -- only
##    a couple of cells in one sample -- so `rare.keep` has something to preserve.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ## g3 has 2 cells in s1 and 8 in s2/s3: present overall (10 >= min.cells) but rare in s1.
  ident <- factor(c(rep("g1", 8), rep("g2", 8), rep("g3", 10)), levels = lev)
  samples <- factor(c(rep("s1", 10), rep("s2", 8), rep("s3", 8)), levels = c("s1", "s2", "s3"))
  prob <- array(rep(c(0.5, 0.2, 0, 0.3, 0.1, 0.4), each = 9), dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "rare_keep", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 26, min_cells = 3L, min_samples = 2L,
    rare_keep = TRUE, nonfilter_keep = FALSE, genes = genes)
  cases[[length(cases) + 1]] <- list(
    name = "rare_noKeep", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 26, min_cells = 3L, min_samples = 2L,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes)
}

## 6. `nonFilter.keep = TRUE`, which stores the unfiltered arrays in two extra slots. The
##    corpus records those too, since `identical()` on the whole object sees them.
{
  lev <- c("g1", "g2"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 5), rep("g2", 3)), levels = lev)
  samples <- factor(rep("s1", 8), levels = "s1")
  prob <- array(rep(c(0.5, 0.2, 0, 0.3, 0.1, 0.4), each = 4), dim = c(2, 2, 6))
  pval <- array(0.01, dim = c(2, 2, 6))
  cases[[length(cases) + 1]] <- list(
    name = "nonfilter_keep", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 8, min_cells = 4L, min_samples = NULL,
    rare_keep = FALSE, nonfilter_keep = TRUE, genes = genes)
}

## 7. All-zero `net$prob`: `num.interaction0 == 0`, so the percentage divides by zero. R's
##    `scales::percent(NaN, accuracy = .1)` is part of the contract.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 4), rep("g2", 4), rep("g3", 2)), levels = lev)
  samples <- factor(c(rep("s1", 5), rep("s2", 5)), levels = c("s1", "s2"))
  prob <- array(0, dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "allzero", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 10, min_cells = 3L, min_samples = 2L,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes)
}

## 8. `min.samples` above the sample count: upstream's own `stop()`, byte for byte.
{
  lev <- c("g1", "g2"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 4), rep("g2", 4)), levels = lev)
  samples <- factor(c(rep("s1", 4), rep("s2", 4)), levels = c("s1", "s2"))
  prob <- array(rep(c(0.5, 0.2, 0, 0.3, 0.1, 0.4), each = 4), dim = c(2, 2, 6))
  pval <- array(0.01, dim = c(2, 2, 6))
  cases[[length(cases) + 1]] <- list(
    name = "too_many_samples", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 8, min_cells = 3L, min_samples = 5L,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes)
}

## 9. All three `type.mean` choices, over two samples, because `FunMean` is selected by a
##    `switch` and the per-sample group means differ between them.
for (tm in c("triMean", "truncatedMean", "thresholdedMean", "median")) {
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 6), rep("g2", 6), rep("g3", 6)), levels = lev)
  samples <- factor(c(rep("s1", 9), rep("s2", 9)), levels = c("s1", "s2"))
  prob <- array(rep(c(0.5, 0.2, 0, 0.3, 0.1, 0.4), each = 9), dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = paste0("mean_", tm), lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 18, min_cells = 3L, min_samples = 2L,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes, type_mean = tm)
}

## 11. `raw.use = FALSE`: the per-sample means come from `data.smooth` instead of
##     `data.signaling`, so the `min.cells` zero-fill and the consistency mask can differ even
##     when the network structure is identical. The smooth matrix here is the group means
##     broadcast back to cells -- which is roughly what `projectData` produces -- which lifts
##     every zero above the threshold and so *excludes nothing* where the raw data excluded a
##     group. The two settings therefore disagree, which is the only way the branch is tested.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 10), rep("g2", 10), rep("g3", 4)), levels = lev)
  samples <- factor(rep("s1", 24), levels = "s1")
  ## A `prob` that varies *within* each slice. `rep(..., each = 9)` gives a constant slice, and
  ## a constant slice makes the row/column zero-fill invisible in the recorded counts -- so the
  ## case would pass whether or not `data.smooth` were consulted.
  prob <- array(seq_len(3 * 3 * 6) / 54, dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "raw_use_false_mincells", lev = lev, lr = lr, ident = ident, samples = samples,
    ## `min.cells = 5`, not 10: the test is `ncell.g <= min.cells`, so 10 against 10-cell groups
    ## excludes *every* group and the whole network is zeroed -- which tests the arithmetic and
    ## nothing else. At 5 only g3 (4 cells) goes, so the zero-fill is visible as a partial one.
    prob = prob, pval = pval, ncell = 24, min_cells = 5L, min_samples = NULL,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes,
    raw_use = FALSE, smooth = "group_means")
}

## 12. `raw.use = FALSE` with a *non-default* `type.mean`, so the branch is exercised with the
##     per-sample mean itself changed rather than only the zero-fill. `median` is the most
##     different of the four from `triMean`.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 8), rep("g2", 8), rep("g3", 8)), levels = lev)
  samples <- factor(c(rep("s1", 12), rep("s2", 12)), levels = c("s1", "s2"))
  prob <- array(seq_len(3 * 3 * 6) / 54, dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "raw_use_false_median", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 24, min_cells = 3L, min_samples = 2L,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes, type_mean = "median",
    raw_use = FALSE, smooth = "group_means")
}

## 13. `raw.use = FALSE` and `rare.keep = TRUE`: the "keep every L-R with any signal in any
##     sample" branch, on the smooth data.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 8), rep("g2", 8), rep("g3", 8)), levels = lev)
  samples <- factor(c(rep("s1", 12), rep("s2", 12)), levels = c("s1", "s2"))
  prob <- array(seq_len(3 * 3 * 6) / 54, dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "raw_use_false_rare_keep", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 24, min_cells = 3L, min_samples = 2L,
    rare_keep = TRUE, nonfilter_keep = FALSE, genes = genes,
    raw_use = FALSE, smooth = "group_means")
}

## 14. `raw.use = FALSE` with the slot *absent*. Upstream writes
##     `if (!("data.smooth" %in% methods::slotNames(object)) == FALSE)` -- `%in%` binds tighter
##     than `!` and tighter than `==`, so the test is `(!in_set) == FALSE`, i.e. `!in_set`, and
##     the error fires when the slot **is** present. The condition is doubly-negated and the
##     message says the opposite of what the test computes; both reproduced.
{
  lev <- c("g1", "g2"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 8), rep("g2", 8)), levels = lev)
  samples <- factor(rep("s1", 16), levels = "s1")
  prob <- array(rep(c(0.5, 0.2, 0.3, 0.1), each = 4), dim = c(2, 2, 6))
  pval <- array(0.01, dim = c(2, 2, 6))
  cases[[length(cases) + 1]] <- list(
    name = "raw_use_false_slot_absent", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 16, min_cells = 3L, min_samples = NULL,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes,
    raw_use = FALSE, smooth = "group_means", drop_smooth = TRUE)
}

## 15. `raw.use = FALSE` with `nonFilter.keep = TRUE`, which is the only configuration in which
##     the min.cells zero-fill is *observable* at all: upstream writes the augmented
##     `net$cellExcludes` and the zeroed `net$prob`, and only then overwrites `object@net` with
##     the pre-augmentation local when `nonFilter.keep` is FALSE. So this case is what pins the
##     `raw.use` branch's arithmetic, and the messages pin the rest.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 10), rep("g2", 10), rep("g3", 4)), levels = lev)
  samples <- factor(rep("s1", 24), levels = "s1")
  prob <- array(seq_len(3 * 3 * 6) / 54, dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "raw_use_false_nonfilter_keep", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 24, min_cells = 5L, min_samples = NULL,
    rare_keep = FALSE, nonfilter_keep = TRUE, genes = genes,
    raw_use = FALSE, smooth = "group_means")
}

## 16. `raw.use = FALSE` with `min.samples = 2` and a smooth matrix that is *not* a function of
##     the group labels. The cases above use shrinkage towards the group mean, which is
##     constant within a group; that makes every cell of a group score identically, the
##     cross-sample comparison ties, and `LR.nonzero` comes out empty -- upstream's `1:0` crash.
##     So both multi-sample `raw.use` cases in this corpus *error*, and neither of them compares a
##     successful result. This one supplies an independent draw instead, which keeps the
##     per-cell variation the sample half needs while still being a matrix the `raw.use` flag has
##     to actually read.
{
  lev <- c("g1", "g2", "g3"); lr <- LR_ALL
  ident <- factor(c(rep("g1", 8), rep("g2", 8), rep("g3", 8)), levels = lev)
  samples <- factor(c(rep("s1", 12), rep("s2", 12)), levels = c("s1", "s2"))
  prob <- array(seq_len(3 * 3 * 6) / 54, dim = c(3, 3, 6))
  pval <- array(0.01, dim = c(3, 3, 6))
  cases[[length(cases) + 1]] <- list(
    name = "raw_use_false_two_samples", lev = lev, lr = lr, ident = ident, samples = samples,
    prob = prob, pval = pval, ncell = 24, min_cells = 3L, min_samples = 2L,
    rare_keep = FALSE, nonfilter_keep = FALSE, genes = genes,
    raw_use = FALSE, smooth = "alt_draw")
}

## ------------------------------------------------------------------- build the smooth data
## A "recipe" rather than a literal matrix, resolved **once, here**, for every case: the
## generator has two independent `for (cs in cases)` loops (the golden dump and the inputs dump)
## and a value computed inside one of them is invisible to the other. The first version built
## the smooth matrix in the golden loop only, so the inputs dump wrote the recipe *string* into
## `sprintf("%a", ...)`, which produced a bare "NAs introduced by coercion" and a `data.smooth`
## slot holding strings.
for (i in seq_along(cases)) {
  if (identical(cases[[i]]$smooth, "alt_draw")) {
    ## An independent draw of the same shape. Deliberately *not* a function of the group
    ## labels: a group-mean-based matrix is constant within a group, which makes the
    ## cross-sample comparison tie and empties `LR.nonzero`.
    cases[[i]]$smooth <- mk_data(cases[[i]]$genes, cases[[i]]$ncell,
                                 seed = 777000 + i)
  } else if (identical(cases[[i]]$smooth, "group_means")) {
    d <- mk_data(cases[[i]]$genes, cases[[i]]$ncell, seed = 20240440 + length(cases))
    grp <- as.integer(cases[[i]]$ident)
    ug <- sort(unique(grp))
    ## `data` is genes x cells, the group labels are per **cell**, so the selection is on the
    ## columns -- and the per-*gene* mean of that block is `rowMeans`, not `colMeans`. Both
    ## halves of that are the same class of mistake: a transposed read of which axis is which.
    ## `colMeans` on a 5 x 10 block returns 10 values, and `vapply` then reports a length
    ## mismatch that points at the fixture rather than at the axis.
    stopifnot(length(grp) == ncol(d),
              length(ug) == length(unique(grp)),
              nrow(d) == length(cases[[i]]$genes),
              ncol(d) == cases[[i]]$ncell)
    ## A plain loop rather than `vapply` + `t()`. The transposition is the whole point of the
    ## expression and `t(vapply(...))` gets it wrong in a way that is invisible until a later
    ## line multiplies the result by another matrix and reports "non-conformable arrays" --
    ## which reads as a shape bug three lines from the cause and is not one.
    gm <- matrix(0, nrow(d), length(ug))
    for (j in seq_along(ug)) gm[, j] <- rowMeans(d[, grp == ug[j], drop = FALSE])
    ## genes on the rows, groups on the columns -- the transpose of the `vapply` result. With
    ## the names the other way round, R does not merely mislabel: assigning a 3-element
    ## dimnames to a 15-row matrix *reshapes* it, and the reshape is silent.
    dimnames(gm) <- list(rownames(d), as.character(ug))
    ## 0.5 * group mean + 0.5 * raw, rather than the group mean outright. The pure group mean is
    ## constant within a group, so every cell of a group scores identically, the cross-sample
    ## comparison ties, and `LR.nonzero` comes out empty -- which is a genuine upstream crash
    ## (`1:0`) but makes every `min.samples >= 2` case collapse onto the error path and compare
    ## nothing. Shrinkage keeps the dropout reduction while retaining per-cell variation.
    ## `gm[, grp]`, not `gm[grp, ]`: `gm` is genes x groups and the broadcast is across the
    ## **columns**, one group mean per cell. Row-indexing yields `ncell x ngroups`, which is the
    ## wrong shape in a way that only shows up once something is done with the result -- and with
    ## no arithmetic attached it produced a plausible-looking matrix that was never the matrix
    ## any of the two code paths read.
    cases[[i]]$smooth <- 0.5 * gm[, grp, drop = FALSE] + 0.5 * d
  }
}

## ------------------------------------------------------------------- run and record
for (cs in cases) {
  data <- mk_data(cs$genes, cs$ncell, seed = 20240440 + length(cases))
  ## A "recipe" rather than a literal matrix, so the smooth data is derived from the same `data`
  ## and the same cell labels: the two matrices then agree cell-for-cell, and any difference in
  ## the result is attributable to the smoothing and not to a different fixture.
  cs$smooth <- if (identical(cs$smooth, "group_means")) {
    grp <- as.integer(cs$ident)
    ug <- sort(unique(grp))
    ## `unname()` because `colMeans` on a matrix with no column names returns a *named* vector,
    ## and a named vector dropped into a matrix assigns names that then disagree with the row
    ## labels -- which `array()` and `dim<-` then reject.
    ## `data` is genes x cells, the group labels are per **cell**, so the selection is on the
    ## columns -- and the per-*gene* mean of that block is `rowMeans`, not `colMeans`. Both
    ## halves of that are the same class of mistake: a transposed read of which axis is which.
    ## `colMeans` on a 5 x 10 block returns 10 values, and `vapply` then reports a length
    ## mismatch that points at the fixture rather than at the axis.
    gm <- t(vapply(ug, function(g) unname(rowMeans(data[, grp == g, drop = FALSE])),
                   numeric(nrow(data))))
    dimnames(gm) <- list(as.character(ug), rownames(data))
    gm[grp, , drop = FALSE]
  } else cs$smooth
  tm <- if (is.null(cs$type_mean)) "triMean" else cs$type_mean
  ru <- isTRUE(cs$raw_use)
  ## `with_smooth` is on whenever `raw_use` is FALSE: upstream reads `data.smooth` in that case,
  ## and a missing slot is a *different* error (the one the `raw_use_missing_smooth` case pins).
  ## `drop_smooth` is the inverse of `with_smooth`: it builds the object *without* the slot, so
  ## the mis-negated `!in %in% slotNames == FALSE` guard fires.
  want_smooth <- !ru && !isTRUE(cs$drop_smooth)
  mkargs <- list(type_mean = tm, raw_use = ru, with_smooth = want_smooth, smooth = cs$smooth)
  o <- do.call(mk, c(list(cs$prob, cs$pval, cs$lev, cs$lr, cs$ident, cs$samples, data, DB),
                     mkargs))
  o2 <- do.call(mk, c(list(cs$prob, cs$pval, cs$lev, cs$lr, cs$ident, cs$samples, data, DB),
                      mkargs))
  ## `cat()` output is captured, not printed: the messages are part of the contract.
  msgs <- character(0)
  out <- withCallingHandlers(
    tryCatch(
      get("filterCommunication", envir = env)(
        o, min.cells = cs$min_cells, min.samples = cs$min_samples,
        rare.keep = cs$rare_keep, nonFilter.keep = cs$nonfilter_keep),
      error = function(e) structure(conditionMessage(e), class = "rerr")),
    message = function(m) {
      msgs <<- c(msgs, conditionMessage(m)); invokeRestart("muffleMessage")
    })
  ## `min.samples` is `NULL` in half the cases, so it is written as a string, not `%d`.
  cat(sprintf("case\t%s\t%d\t%s\t%d\t%s\t%s\n", cs$name, cs$min_cells,
              if (is.null(cs$min_samples)) "NA" else as.character(cs$min_samples),
              length(cs$samples), cs$rare_keep, cs$nonfilter_keep), file = q)
  if (inherits(out, "rerr")) {
    cat(sprintf("filter\t%s\tERROR\t%s\n", cs$name, out), file = q)
    next
  }
  ## The `cat()`s go to stdout, so they are not captured by the message handler. Re-run with
  ## `capture.output` to record them exactly.
  o3 <- do.call(mk, c(list(cs$prob, cs$pval, cs$lev, cs$lr, cs$ident, cs$samples, data, DB),
                      mkargs))
  txt <- capture.output(
    invisible(get("filterCommunication", envir = env)(
      o3, min.cells = cs$min_cells, min.samples = cs$min_samples,
      rare.keep = cs$rare_keep, nonFilter.keep = cs$nonfilter_keep)))
  for (ln in txt) cat(sprintf("msg\t%s\t%s\n", cs$name, gsub("\t", "<TAB>", ln)), file = q)
  ## The interaction counts, as upstream computes them:
  ##   num.interaction0 <- sum(net$prob > 0)   before any filtering
  ##   num.interaction1 <- sum(net$prob > 0)   after the min.cells zero-fill
  ##   num.interaction2 <- sum(net$prob > 0)   after the cross-sample mask
  ##
  ## 1 and 2 are only *observable* when `nonFilter.keep = TRUE`. With the default `FALSE`,
  ## upstream writes the augmented `net$cellExcludes` and then overwrites `object@net` with the
  ## pre-augmentation local, so the returned `prob` is the unfiltered array and all three counts
  ## are equal -- by design, and the reason the `nonFilter.keep` case exists at all. The
  ## `raw.use` fixtures therefore discriminate through the `cat()` output (the percentages are
  ## computed from these counts) and through the `nonFilter.keep = TRUE` case, not here.
  ##
  ## The first version wrote `sum(cs$prob > 0) - 0L` for the second, which is just the first
  ## again, so the count meant to show the zero-fill never showed it.
  cat(sprintf("filter\t%s\tn_inter0=%d\tn_inter1=%d\tn_inter2=%d\n", cs$name,
              sum(o@net$prob > 0),
              sum(out@net$prob > 0),
              sum(out@net$prob > 0)), file = q)
  cat(sprintf("filter\t%s\tprob=%s\n", cs$name, fmt(out@net$prob)), file = q)
  cat(sprintf("filter\t%s\tdim=%s\n", cs$name,
              paste(dim(out@net$prob), collapse = "x")), file = q)
  cat(sprintf("filter\t%s\tdimnames=%s\n", cs$name,
              paste(unlist(lapply(dimnames(out@net$prob), paste, collapse = ",")),
                    collapse = "|")), file = q)
  if (cs$nonfilter_keep) {
    cat(sprintf("filter\t%s\tprob_nonFilter=%s\n", cs$name, fmt(out@net$prob.nonFilter)), file = q)
    cat(sprintf("filter\t%s\tpval_nonFilter=%s\n", cs$name, fmt(out@net$pval.nonFilter)), file = q)
  }
  ## The net's other slots must be untouched; `identical()` on the whole object sees them.
  cat(sprintf("filter\t%s\tpval=%s\n", cs$name, fmt(out@net$pval)), file = q)
  cat(sprintf("filter\t%s\tslot_names=%s\n", cs$name,
              paste(sort(names(out@net)), collapse = ",")), file = q)
}
close(q)
file.rename("tests/fixtures/filter_golden.txt.part", "tests/fixtures/filter_golden.txt")
## ------------------------------------------------------------------ dump the inputs
## One `|`-separated metadata line per case, then `prob`, `pval`, `data` and the per-cell
## label rows, each as `%a` hex or comma-separated text. The Rust test needs the *same*
## inputs, and a fixture that needs an RNG gets written out rather than re-derived.
##
## The metadata is a single line because a tagged-record-per-field format needed a
## field-index parser on both sides, and R's `strsplit` drops a trailing empty field --
## which silently shifts every later field when the last value is empty.
qi <- file("tests/fixtures/filter_inputs.tsv.part", "wt")
for (cs in cases) {
  data <- mk_data(cs$genes, cs$ncell, seed = 20240440 + length(cases))
  cat(sprintf("case\t%s\t%s\n", cs$name, paste(c(
    paste(cs$lev, collapse = ","),
    paste(cs$lr, collapse = ","),
    as.character(cs$ncell),
    as.character(cs$min_cells),
    if (is.null(cs$min_samples)) "NA" else as.character(cs$min_samples),
    cs$rare_keep,
    cs$nonfilter_keep,
    if (is.null(cs$type_mean)) "triMean" else cs$type_mean,
    "0.1",
    ## `raw.use`, appended last so the earlier field indices are untouched. **Defaults to
    ## TRUE** when the case does not set it, which is what upstream's own default is. Without the
    ## default every pre-existing case silently became `raw.use = FALSE` and the whole corpus
    ## started reading a `data.smooth` slot that only two of the cases have.
    if (!is.null(cs$raw_use) && !isTRUE(cs$raw_use)) "FALSE" else "TRUE",
    ## The `data.smooth` state, three values rather than two. `dropped` is the case where the
    ## object has **no** `data.smooth` slot at all even though a matrix was constructed for the
    ## record: upstream's guard
    ##   `if (!("data.smooth" %in% methods::slotNames(object)) == FALSE)`
    ## does not fire in that state, because `%in%` binds tighter than both `!` and `==`, so the
    ## test is `(!in_set) == FALSE`, i.e. "the slot **is** present". The condition is doubly
    ## negated and the message says the opposite of what it computes; both reproduced, and the
    ## absent state is a third distinct input rather than a variant of the present one.
    if (is.null(cs$smooth)) "-" else if (isTRUE(cs$drop_smooth)) "dropped" else "yes"
  ), collapse = "|")), file = qi)
  cat(sprintf("idents\t%s\t%s\n", cs$name, paste(as.integer(cs$ident) - 1L, collapse = ",")), file = qi)
  cat(sprintf("samples\t%s\t%s\n", cs$name, paste(as.integer(cs$samples) - 1L, collapse = ",")), file = qi)
  cat(sprintf("genes\t%s\t%s\n", cs$name, paste(cs$genes, collapse = ",")), file = qi)
  cat("prob\n", file = qi); writeLines(paste(sprintf("%a", as.numeric(cs$prob)), collapse = " "), qi)
  cat("pval\n", file = qi); writeLines(paste(sprintf("%a", as.numeric(cs$pval)), collapse = " "), qi)
  cat("data\n", file = qi); writeLines(paste(sprintf("%a", as.numeric(data)), collapse = " "), qi)
  ## `smooth` goes *after* `data`, so the Rust reader's fixed offsets hold and the record set
  ## stays in the same order on both sides. Written straight after the `case` line it shifted
  ## every later offset by two, and the failure surfaced as "expected idents, found smooth" in a
  ## completely different test.
  if (!is.null(cs$smooth)) {
    cat(sprintf("smooth\t%s\n", cs$name), file = qi)
    writeLines(paste(sprintf("%a", as.numeric(cs$smooth)), collapse = " "), qi)
  }
  ## The DB, once: the same for every case, so it goes in a single record.
  if (cs$name == cases[[1]]$name) {
    cat(sprintf("db_complexes\t%s\n", paste(rownames(DB$complex), collapse = ",")), file = qi)
    for (cn in colnames(DB$complex)) {
      cat(sprintf("db_subunit\t%s\t%s\n", cn, paste(DB$complex[[cn]], collapse = ",")), file = qi)
    }
    cat(sprintf("db_symbols\t%s\n", paste(DB$geneInfo$Symbol, collapse = ",")), file = qi)
    cat(sprintf("db_ligand\t%s\n", paste(DB$interaction$ligand, collapse = ",")), file = qi)
    cat(sprintf("db_receptor\t%s\n", paste(DB$interaction$receptor, collapse = ",")), file = qi)
  }
}
close(qi)
file.rename("tests/fixtures/filter_inputs.tsv.part", "tests/fixtures/filter_inputs.tsv")
