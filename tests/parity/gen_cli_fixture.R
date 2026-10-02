# Write the standalone CLI's input files, and record what pinned upstream computes from them.
#
# The CLI exists so `r-core`'s numerics are usable without R. That is only worth claiming if the CLI
# produces the *same* numbers R does, so for each configuration below this writes one input file and
# two things happen with it: pinned upstream's `computeCommunProb` reads the same data in R, and
# `tests/parity/check_cli.R` runs the `cellchatrs` binary and requires `identical()` on the full
# arrays. The R side here does not go through the shim at all -- the point is that the input file is a
# complete, self-contained description, not a projection of an S4 object.
#
# Two format decisions, both load-bearing:
#
#   * **Values are `%a` hex floats.** Every `f64` round-trips exactly, and a mis-transcribed digit
#     becomes a *different double* rather than a near-miss, so a parity failure points at the kernel
#     rather than at the file. With `%.17g` a lost digit is still a valid number and the failure
#     looks like a 1-ulp rounding bug.
#   * **The matrix is written gene-major, one line per gene.** R's column-major flat buffer for a
#     genes-by-cells matrix has, for each cell, all genes consecutive; writing one line per gene and
#     having the reader index `data[gene * n_cells + cell]` is the explicit form of that
#     transposition. `crates/r-core/tests/expr_parity.rs` records what happens when the inversion is
#     implicit instead: the gene *set* stays correct while the *values* permute, every name still
#     resolves, and every number still looks plausible.
#
# Usage:  R_LIBS=.rlib R --vanilla -f tests/parity/gen_cli_fixture.R
suppressWarnings(suppressMessages({
  library(methods); library(Matrix); library(collapse); library(dplyr)
}))

CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)

## The same L-R set and expression matrix the rest of the suite uses, so a failure in the CLI gate is
## about the CLI and not about different data.
src <- readLines("tests/parity/gen_prob_golden.R")
stop_at <- grep("^q <- file", src)[1]
eval(parse(text = paste(src[seq_len(stop_at - 1)], collapse = "\n")))

n_genes <- nrow(data.signaling)
n_cells <- ncol(data.signaling)
groups <- levels(cell_group)
n_groups <- length(groups)

## Trimmed to a size the goldens stay small enough to read, while keeping every structure that makes
## the kernel interesting: complexes, a cofactor, and an agonist/antagonist pair.
LRsig <- LRsig[seq_len(min(nrow(LRsig), 6)), , drop = FALSE]
n_lr <- nrow(LRsig)

## The configurations, which are the axes that change the *arithmetic* rather than only the shape.
##
## One fixture would prove the CLI can compute something. These prove it computes the same thing on
## the paths that are actually different in the kernel:
##
##   * all four `type.mean` values, including the two that read `trim` -- `truncatedMean` and
##     `thresholdedMean` are the only ones for which a wrong `trim` is visible at all;
##   * `population.size` both ways, which switches `P3.boot`/`P4.boot` between `matrix(1)` and the
##     permuted group proportions;
##   * `nboot` from 1 to 9, which exercises the exact-integer p-value counting and the degenerate
##     `1/nboot` case at `nboot = 1`;
##   * `Kh` six orders of magnitude apart, which moves the Hill coefficient into and out of the
##     regime where `Kh^n` dominates `dataLR^n`;
##   * `n` at 1 and 2, the exponent on `dataLR`;
##   * `raw.use` both ways, which switches the `data / max(data)` step for the `data / rowsum` one;
##   * and two seeds, so the bootstrap null differs and `Pval` is not trivially reproducible.
##
## `type` is given as a *prefix* for two of them on purpose. `match.arg` accepts a unique prefix and
## records the matched spelling, so `tri` and `thr` exercise that path -- and if the CLI recorded the
## typed string rather than the matched one, the golden's `type_mean` would disagree with it.
CONFIGS <- list(
  list(name = "tri_ps0",      type = "triMean",         trim = 0.1,  ps = FALSE, raw = TRUE,
       nboot = 5L, seed = 20240431L, Kh = 0.5, n = 1),
  list(name = "tri_ps1",      type = "triMean",         trim = 0.1,  ps = TRUE,  raw = TRUE,
       nboot = 5L, seed = 20240431L, Kh = 0.5, n = 1),
  list(name = "matcharg_tri", type = "tri",             trim = 0.1,  ps = FALSE, raw = TRUE,
       nboot = 3L, seed = 7L,      Kh = 0.5, n = 1),
  list(name = "matcharg_thr", type = "thr",             trim = 0.25, ps = FALSE, raw = TRUE,
       nboot = 4L, seed = 11L,     Kh = 0.5, n = 1),
  list(name = "trunc_t25",    type = "truncatedMean",   trim = 0.25, ps = FALSE, raw = TRUE,
       nboot = 4L, seed = 13L,     Kh = 0.5, n = 1),
  list(name = "thresh_t30",   type = "thresholdedMean", trim = 0.3,  ps = TRUE,  raw = TRUE,
       nboot = 4L, seed = 3L,      Kh = 0.5, n = 1),
  list(name = "median",       type = "median",          trim = 0.1,  ps = FALSE, raw = TRUE,
       nboot = 6L, seed = 4L,      Kh = 0.5, n = 1),
  list(name = "nboot1",       type = "triMean",         trim = 0.1,  ps = FALSE, raw = TRUE,
       nboot = 1L, seed = 5L,      Kh = 0.5, n = 1),
  list(name = "nboot9",       type = "triMean",         trim = 0.1,  ps = FALSE, raw = TRUE,
       nboot = 9L, seed = 6L,      Kh = 0.5, n = 1),
  list(name = "kh_tiny",      type = "triMean",         trim = 0.1,  ps = FALSE, raw = TRUE,
       nboot = 4L, seed = 8L,      Kh = 1e-3, n = 1),
  list(name = "kh_huge",      type = "triMean",         trim = 0.1,  ps = FALSE, raw = TRUE,
       nboot = 4L, seed = 9L,      Kh = 1e3,  n = 1),
  list(name = "n2",           type = "triMean",         trim = 0.1,  ps = FALSE, raw = TRUE,
       nboot = 5L, seed = 10L,     Kh = 0.5, n = 2),
  list(name = "raw_use_false", type = "triMean",        trim = 0.1,  ps = FALSE, raw = FALSE,
       nboot = 5L, seed = 12L,     Kh = 0.5, n = 1),
  list(name = "seed_b",       type = "triMean",         trim = 0.1,  ps = TRUE,  raw = TRUE,
       nboot = 5L, seed = 99999L,  Kh = 0.5, n = 1)
)

## ------------------------------------------------------------------ the shared fixture
## `MiniCellChat` is the same minimal S4 class the rest of the parity suite uses:
## `computeCommunProb` reads `@data.signaling`, `@LR$LRsig`, `@DB`, `@idents` and
## `@options$datatype`, and nothing else.
## `data.smooth` is in the class, and is *populated*, because `raw.use = FALSE` reads it:
## `if (raw.use) data <- as.matrix(object@data.signaling) else data <- as.matrix(object@data.smooth)`.
## Declaring the slot without filling it is not enough -- `as.matrix(NULL)` raised "'data' must be of a
## vector type, was 'NULL'", and `computeAveExpr` separately needs the slot to *exist* for
## `methods::slot` to validate the name. Both failures are properties of the fixture, not the port, and
## both are avoidable by supplying the matrix upstream's own smoothing step would.
##
## The values are the signaling matrix divided by its row sums, which is a plausible "projected"
## matrix and -- more to the point -- differs from `data.signaling`, so the `raw.use` axis is genuinely
## exercised rather than being a no-op comparison.
setClass("MiniCellChat", representation(
  data.signaling = "Matrix",
  data.smooth = "ANY",
  LR = "list", LRsig = "data.frame", DB = "ANY",
  idents = "factor", options = "list", net = "ANY"
))
obj <- new("MiniCellChat",
           data.signaling = data.signaling,
           data.smooth = {
             ## A row-sum-normalised copy, so `raw.use = FALSE` is a different matrix and not the
             ## same numbers relabelled. `drop = FALSE` because a one-gene fixture would otherwise
             ## collapse to a vector and `as.matrix` would take the other branch.
             z <- as.matrix(data.signaling)
             z / pmax(rowSums(z) / ncol(z), .Machine$double.xmin)
           },
           LR = list(LRsig = LRsig), LRsig = LRsig,
           DB = list(complex = DB$complex, cofactor = DB$cofactor),
           idents = cell_group,
           options = list(datatype = "RNA", mode = "single", db = normalizePath(DBDIR)))

## Upstream is sourced into a fresh environment, never reached through this package. The `globalenv()`
## parent is load-bearing and `cellchatrs_upstream_env()` explains why: a namespace lookup chain stops
## before the search path, so upstream's own `rowSums` would find `base::rowSums` rather than
## `Matrix`'s method for an `lgCMatrix`.
env <- new.env(parent = globalenv())
for (p in c("Matrix", "collapse", "dplyr", "future", "rlang")) {
  suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
}
invisible(utils::capture.output(suppressWarnings(suppressMessages(
  for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R", "visualization.R")) {
    sys.source(file.path(CC, "R", f), envir = env, keep.source = FALSE)
  }))))

## ------------------------------------------------------------------ one pair per configuration
## The matrix, the group assignment and the L-R table are identical across configurations; only the
## `type`/`trim`/`population.size`/`nboot`/`seed`/`Kh`/`n` block changes. So the expensive part --
## writing 24 x 60 hex floats -- is done once and reused, and a difference between two goldens is
## unambiguously a difference in the kernel's arithmetic for that configuration.
hex_lines <- function(m) vapply(seq_len(nrow(m)), function(g) {
  paste(sprintf("%a", as.numeric(m[g, ])), collapse = " ")
}, "")
gene_lines <- hex_lines(data.signaling)

cell_names <- colnames(data.signaling)
group_line <- paste(as.integer(cell_group) - 1L, collapse = " ")

write_input <- function(path, cfg) {
  con <- file(path, "wt")
  cat("version\t1\n", file = con)
  cat(sprintf("groups\t%d\n", n_groups), file = con)
  cat(sprintf("genes\t%d\n", n_genes), file = con)
  cat(sprintf("cells\t%d\n", n_cells), file = con)
  cat(sprintf("lr\t%d\n", n_lr), file = con)
  writeLines(groups, con)
  writeLines(rownames(data.signaling), con)
  cat("values\n", file = con)
  writeLines(gene_lines, con)
  if (!cfg$raw) {
    ## Written only for `raw.use = FALSE`. The format *rejects* the block when `raw_use` is `TRUE` and
    ## *requires* it when `FALSE`, so a fixture cannot carry a matrix the kernel will not read, nor
    ## omit one it will. That is stricter than it needs to be for a file format and is the point: a
    ## format that quietly accepts a matrix nobody reads stops describing the computation.
    cat("values_smooth\n", file = con)
    writeLines(hex_lines(as.matrix(obj@data.smooth)), con)
  }
  writeLines(cell_names, con)
  cat("cell_group\n", file = con)
  cat(sprintf("%s\n", group_line), file = con)
  cat("lr_columns\tinteraction_name\tligand\treceptor\tagonist\tantagonist\tco_A_receptor\tco_I_receptor\n",
      file = con)
  for (i in seq_len(n_lr)) {
    ## `""` for absent, never omitted: the reader treats a short row as an error precisely so a
    ## missing trailing field cannot be read as an empty string.
    writeLines(paste(c(rownames(LRsig)[i], LRsig$ligand[i], LRsig$receptor[i],
                       LRsig$agonist[i], LRsig$antagonist[i],
                       LRsig$co_A_receptor[i], LRsig$co_I_receptor[i]),
                     collapse = "\t"), con)
  }
  cat(sprintf("type\t%s\n", cfg$type), file = con)
  cat(sprintf("trim\t%.17g\n", cfg$trim), file = con)
  cat(sprintf("population_size\t%s\n", cfg$ps), file = con)
  cat(sprintf("raw_use\t%s\n", cfg$raw), file = con)
  cat(sprintf("nboot\t%d\n", cfg$nboot), file = con)
  cat(sprintf("seed\t%d\n", cfg$seed), file = con)
  cat(sprintf("Kh\t%.17g\n", cfg$Kh), file = con)
  cat(sprintf("n\t%.17g\n", cfg$n), file = con)
  close(con)
}

## `aggregateNet`'s unfiltered branch, written out rather than called, because upstream's function
## mutates and returns the *object* and the golden wants the two matrices on their own. The four
## statements are `modelling.R:404-411` verbatim, including the `is.na` guards, which matter whenever
## a group pair has no interactions at all.
##
## Named `aggregate_branch`, **not** `aggregate`, and that is not a style preference. The sourced
## upstream environment has `parent = globalenv()` -- which `cellchatrs_upstream_env()` explains is
## load-bearing for `Matrix::rowSums` -- so a function this script defines in the global environment
## shadows the same name *for upstream's own body*. Defining `aggregate` here made every
## `computeCommunProb` call fail with "unused argument (FUN = FunMean)": upstream's bootstrap does
## `aggregate(t(data.use), list(groupboot), FUN = FunMean)`, and it was reaching this four-argument
## function instead of `base::aggregate`. Every golden recorded that error and the run "succeeded".
## The general rule for any script that sources upstream into `globalenv()` is the same one the shim
## follows: never define a name upstream uses, or the reference stops being the reference.
aggregate_branch <- function(prob, pval, k, thresh) {
  pv <- pval
  pv[prob == 0] <- 1
  z <- prob
  z[pv >= thresh] <- 0
  count <- apply(z > 0, c(1, 2), sum)
  weight <- apply(z, c(1, 2), sum)
  weight[is.na(weight)] <- 0
  count[is.na(count)] <- 0
  list(count = count, weight = weight)
}

## Remove the previous run's pairs first. A regeneration that fails halfway would otherwise leave a
## directory where some goldens are new and some are stale, and a gate looping over the directory
## would compare the CLI against a mixture of two different kernels. Deleted by name, not by pattern
## that could catch an unrelated fixture.
unlink(sprintf("tests/fixtures/cli_input_%s.tsv", vapply(CONFIGS, `[[`, "", "name")))
unlink(sprintf("tests/fixtures/cli_golden_%s.txt", vapply(CONFIGS, `[[`, "", "name")))

n_errored <- 0L
for (cfg in CONFIGS) {
  inp <- sprintf("tests/fixtures/cli_input_%s.tsv", cfg$name)
  gold <- sprintf("tests/fixtures/cli_golden_%s.txt", cfg$name)
  write_input(inp, cfg)

  ## `computeAveExpr` is recorded per group so the CLI's `ave_expr` block is gated too, and it is the
  ## only place `trim` is visible for the two trimmed means.
  res <- tryCatch(
    suppressWarnings(suppressMessages(get("computeCommunProb", envir = env)(
      obj, type = cfg$type, trim = cfg$trim, population.size = cfg$ps, raw.use = cfg$raw,
      nboot = cfg$nboot, seed.use = cfg$seed, Kh = cfg$Kh, n = cfg$n))),
    error = function(e) structure(conditionMessage(e), class = "cli_error"))

  ## An error is part of the contract, so the golden records it and the gate requires the CLI to
  ## fail with the same message. The check has to come *before* any header line that reads
  ## `res@options$parameter`: an earlier version wrote the header first and died with "no applicable
  ## method for `@` applied to an object of class \"cli_error\"", which names neither the
  ## configuration nor the message upstream actually raised.
  if (inherits(res, "cli_error")) {
    qg <- file(gold, "wt")
    cat("version\t1\n", file = qg)
    cat(sprintf("config\t%s\n", cfg$name), file = qg)
    cat("error\t", file = qg)
    cat(as.character(res), "\n", file = qg)
    close(qg)
    n_errored <- n_errored + 1L
    next
  }

  qg <- file(gold, "wt")
  cat("version\t1\n", file = qg)
  cat(sprintf("config\t%s\n", cfg$name), file = qg)
  ## The *matched* spelling, not the typed one. `match.arg("tri")` records `triMean`, and the CLI
  ## has to agree or this line is the check.
  cat(sprintf("type\t%s\n", res@options$parameter$type.mean), file = qg)
  cat(sprintf("nboot\t%d\n", cfg$nboot), file = qg)
  cat(sprintf("groups\t%d\n", n_groups), file = qg)
  cat(sprintf("lr\t%d\n", n_lr), file = qg)
  cat(sprintf("dim\t%d,%d,%d\n", n_groups, n_groups, n_lr), file = qg)

  p <- array(as.numeric(res@net$prob), dim = c(n_groups, n_groups, n_lr))
  v <- array(as.numeric(res@net$pval), dim = c(n_groups, n_groups, n_lr))
  ## `%.17g` rather than `%a` here: the golden is a *human-inspectable* record of what R produced, and
  ## 17 significant digits is the shortest width that round-trips a double uniquely. The hex spelling
  ## is for the wire, in the input file, where a transcription error has to be loud.
  cat("prob\n", file = qg); writeLines(sprintf("%.17g", p), qg)
  cat("pval\n", file = qg); writeLines(sprintf("%.17g", v), qg)
  ## `computeAveExpr`'s own `match.arg` offers only `triMean`, `truncatedMean` and `median` -- three
  ## choices, not four. `thresholdedMean` is a `computeCommunProb` option that `computeAveExpr` cannot
  ## be asked for, and asking produces upstream's error verbatim:
  ## 'arg' should be one of "triMean", "truncatedMean", "median".
  ## That is a genuine asymmetry in the two functions and is recorded rather than worked around, so the
  ## `ave_expr` block is present only for the types it can answer for.
  ave <- tryCatch(
    as.numeric(get("computeAveExpr", envir = env)(obj, type = cfg$type, trim = cfg$trim)),
    error = function(e) NULL)
  if (!is.null(ave)) {
    cat("ave_expr\n", file = qg)
    writeLines(sprintf("%.17g", ave), qg)
  }
  ag <- aggregate_branch(p, v, n_groups, 0.05)
  cat("aggregate_count\n", file = qg)
  writeLines(apply(ag$count, 1, function(r) paste(sprintf("%.17g", r), collapse = "\t")), qg)
  cat("aggregate_weight\n", file = qg)
  writeLines(apply(ag$weight, 1, function(r) paste(sprintf("%.17g", r), collapse = "\t")), qg)
  cat("dimnames_source\n", file = qg); writeLines(paste(groups, collapse = "\t"), qg)
  cat("dimnames_interaction\n", file = qg)
  writeLines(paste(dimnames(res@net$prob)[[3]], collapse = "\t"), qg)
  close(qg)
}

cat(sprintf("wrote %d input/golden pairs (%d of which record an upstream error)\n",
            length(CONFIGS), n_errored), file = stderr())
