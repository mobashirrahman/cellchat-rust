# End-to-end reproduction of the upstream CellChat tutorial's computational pipeline.
#
# This is requirement 3's "end-to-end tutorial reproduction": the vignette
# `tutorial/CellChat-vignette.Rmd` in the pinned upstream, run step by step on real data,
# twice -- once with every function from the pinned upstream, once with every function from
# the installed shim -- requiring `identical()` on every non-visual slot after each stage.
#
# Two deliberate departures from a literal transcription, both stated:
#
#   * No plots. Scope decision 14.1 keeps visualization.R and app.R in R; a ggplot object
#     carries an environment and a call, so `identical()` on two plots is not a well-defined
#     acceptance test. The numeric slots feeding the plots are all compared, which is what
#     determines what the plots would draw.
#   * The tutorial's own input file (`data_humanSkin_CellChat.rda`, an author's laptop path) is
#     not available here. The input is rebuilt from the same authors' Figshare object this repo
#     already benchmarks (`humanSkin.rda`, 7,563 cells): log1p CPM over the signalling genes and
#     the L-R table filtered to resolvable pairs, exactly as `bench-runner/fixture.R` builds it.
#     The *steps* are the tutorial's, in the tutorial's order, with the tutorial's arguments.
#
# Pipeline (vignette order): createCellChat -> subsetData -> identifyOverExpressedGenes ->
# identifyOverExpressedInteractions -> computeCommunProb(triMean defaults, nboot = 100) ->
# filterCommunication(min.cells = 10) -> computeCommunProbPathway -> aggregateNet ->
# subsetCommunication spot checks -> netAnalysis_computeCentrality(slot netP).
#
# Timing note: the upstream side is the slow one on purpose. nboot = 100 is the tutorial's own
# default, and the point of reproducing it verbatim is that the shim earns its result against the
# exact computation a user runs -- not a downsampled one. Expect ~2 min for the reference kernel.
#
# Usage: R_LIBS=.rlib R --vanilla -f tests/parity/tutorial_repro.R
suppressWarnings(suppressMessages({
  library(methods); library(Matrix); library(collapse); library(dplyr)
}))
suppressPackageStartupMessages(library(Matrix))
suppressWarnings(suppressMessages(library(cellchatrs)))

CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")
SKIN <- Sys.getenv("TUTORIAL_DATA", "/scratch/mdra00001/tmp/opencode/data/humanSkin.rda")

if (tolower(Sys.getenv("CELLCHATRS_FALLBACK", "0")) %in% c("1", "true", "yes", "on")) {
  stop("CELLCHATRS_FALLBACK is set: the shim side would answer from upstream and this gate ",
       "would compare upstream against itself. Unset it.", call. = FALSE)
}

## igraph must be present for the centrality stage. If it is absent the script stops before any
## comparison rather than comparing two objects that both skipped the stage -- a skipped stage that
## still reports IDENTICAL is the failure mode this file exists to prevent.
if (!requireNamespace("igraph", quietly = TRUE)) {
  stop("igraph is not installed, so netAnalysis_computeCentrality cannot run on either side. ",
       "Install it rather than skipping the stage.", call. = FALSE)
}
suppressPackageStartupMessages(library(igraph))

fails <- 0L
n_cmp <- 0L
cmp <- function(label, ok, detail = "") {
  n_cmp <<- n_cmp + 1L
  cat(sprintf("%-58s %s%s\n", label, if (ok) "ok" else "MISMATCH",
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!ok) fails <<- fails + 1L
  invisible(ok)
}
q <- function(x) {
  invisible(utils::capture.output(value <- suppressWarnings(suppressMessages(x))))
  if (is.null(value)) stop("q(): NULL return; a comparison against it would be vacuous")
  value
}

## ---------------------------------------------------------------- upstream reference env
## Sourced once, into its own environment with parent = globalenv() (the lookup-chain rationale is
## on `cellchatrs_upstream_env`). Everything on the reference side is reached as `UP$fn` so there is
## no ambiguity about which implementation produced a value.
## `CellChat_class.R` goes into `globalenv()`, and ONLY there -- not into `UP`. This is load-bearing
## and the failure for getting it wrong is far from the cause.
##
## `methods::new("CellChat")` and `is(x, "Seurat")` resolve the class through
## `getClassDef(Class, where = topenv(parent.frame()))`, i.e. the *top* of the calling chain, not the
## calling environment. A `setClass` executed inside a private `new.env()` registers the definition
## under a synthetic timestamp package name ("Created a package name ... when none found"), and the
## first `new()`/`is()` then tries to *load* that package: `unable to find required package
## '2026-10-01 ...'`, raised inside the kernel call. Sourcing the class file into `globalenv()`
## registers it where every lookup chain ends, so one definition serves the reference env, the
## shim's cached env, and this script's own `slot()` calls. (`tests/parity/check_merge.R` does the
## same; `setClass` on an already-defined class returns quietly, so this is idempotent.)
## The five function files go into `UP` as before -- none of them defines the class.
UP <- new.env(parent = globalenv())
for (p in c("Matrix", "collapse", "dplyr", "future", "rlang")) {
  suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
}
for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R", "visualization.R")) {
  suppressWarnings(suppressMessages(
    sys.source(file.path(CC, "R", f), envir = UP, keep.source = FALSE)))
}
suppressWarnings(suppressMessages(
  sys.source(file.path(CC, "R", "CellChat_class.R"), envir = globalenv(),
             keep.source = FALSE)))
stopifnot(!is.null(methods::getClass("CellChat", where = globalenv())))

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)

## ---------------------------------------------------------------- input, as the vignette builds it
raw <- new.env(); load(SKIN, envir = raw); skin <- raw[[ls(raw)[1]]]
counts <- skin$data
grp <- if (!is.null(skin$meta)) factor(skin$meta$labels) else factor(skin$labels)
cat(sprintf("input: %d genes x %d cells, %d groups\n", nrow(counts), ncol(counts), nlevels(grp)))

## Log1p CPM over the full matrix first, exactly as the tutorial preprocesses, then hand both sides
## the same normalised matrix so the comparison starts below any preprocessing choice.
lib <- Matrix::colSums(counts)
cpm <- counts
cpm@x <- cpm@x / rep(lib, diff(cpm@p)) * 1e4
norm <- log1p(as.matrix(cpm))
meta <- data.frame(labels = grp, row.names = colnames(counts))

## The shim does not override createCellChat: it takes no numeric arguments, so there is nothing to
## accelerate and the reference constructor is the constructor. It is called unqualified -- it lives
## in `globalenv()` now (see above), not in `UP`, so `UP$createCellChat` is NULL and calling it fails
## with the unhelpful "attempt to apply non-function". Both sides start from the same object by
## construction, and the first compared step is subsetData.
obj_up <- createCellChat(object = norm, meta = meta, group.by = "labels")
obj_rs <- obj_up
cmp("createCellChat identical by construction", identical(obj_up, obj_rs))

## Tutorial-verbatim: the vignette assigns the database next (`cellchat@DB <- CellChatDB.human`),
## because `createCellChat` leaves `@DB` empty and `subsetData` selects from it -- without this both
## sides fail identically inside `dplyr::select` on a NULL table, which would compare two errors
## rather than two analyses. Same assignment on both sides, so still identical by construction.
obj_up@DB <- DB
obj_rs@DB <- DB
cmp("DB assignment identical by construction", identical(obj_up, obj_rs))

## Each stage runs the reference on obj_up and the shim on obj_rs, then compares every non-visual
## slot. The slot list is explicit: `images` and `dr` are visualisation products and are out of scope
## per 14.1, `data`/`data.signaling`/`data.smooth` are inputs, and `options$run.time` is wall-clock.
SLOTS <- c("data.signaling", "idents", "meta", "net", "netP", "LR", "DB", "var.features", "options")

stage <- function(label, up_call, rs_call, slots = SLOTS, exclude_options = c("run.time")) {
  t0 <- Sys.time()
  a <- q(up_call(obj_up))
  b <- q(rs_call(obj_rs))
  dt <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
  obj_up <<- a
  obj_rs <<- b
  ok_all <- TRUE
  for (s in slots) {
    sa <- tryCatch(methods::slot(a, s), error = function(e) paste("ERR", conditionMessage(e)))
    sb <- tryCatch(methods::slot(b, s), error = function(e) paste("ERR", conditionMessage(e)))
    if (s == "options" && is.list(sa) && is.list(sb)) {
      sa <- sa[setdiff(names(sa), exclude_options)]
      sb <- sb[setdiff(names(sb), exclude_options)]
    }
    if (!identical(sa, sb)) {
      ok_all <- FALSE
      da <- utils::capture.output(str(utils::head(sa, 2)))
      cat(sprintf("    slot %s differs: up %s\n", s, paste(da, collapse = " | ")))
    }
  }
  cmp(sprintf("%s [%.0fs]", label, dt), ok_all)
}

## ---------------------------------------------------------------- the vignette, step by step
stage("1 subsetData",
      function(o) UP$subsetData(o),
      function(o) cellchatrs::subsetData(o))

stage("2 identifyOverExpressedGenes (tutorial defaults: presto fast path)",
      function(o) UP$identifyOverExpressedGenes(o),
      function(o) cellchatrs::identifyOverExpressedGenes(o))

stage("3 identifyOverExpressedInteractions",
      function(o) UP$identifyOverExpressedInteractions(o),
      function(o) cellchatrs::identifyOverExpressedInteractions(o))

stage("4 computeCommunProb triMean, tutorial defaults incl. nboot = 100",
      function(o) UP$computeCommunProb(o, type = "triMean"),
      function(o) cellchatrs::computeCommunProb(o, type = "triMean"))

stage("5 filterCommunication min.cells = 10",
      function(o) UP$filterCommunication(o, min.cells = 10),
      function(o) cellchatrs::filterCommunication(o, min.cells = 10))

stage("6 computeCommunProbPathway",
      function(o) UP$computeCommunProbPathway(o),
      function(o) cellchatrs::computeCommunProbPathway(o))

stage("7 aggregateNet",
      function(o) UP$aggregateNet(o),
      function(o) cellchatrs::aggregateNet(o))

## ---------------------------------------------------------------- 8 downstream spot checks
## `subsetCommunication` returns data frames, not objects, so it is compared directly on both slot
## names -- the interaction table and the pathway table, which is what the tutorial queries next.
for (sn in c("net", "netP")) {
  a <- q(UP$subsetCommunication(obj_up, slot.name = sn))
  b <- q(cellchatrs::subsetCommunication(obj_rs, slot.name = sn))
  ok <- identical(a, b)
  cmp(sprintf("8 subsetCommunication slot %s (%d rows)", sn, NROW(a)), ok)
  if (!ok) {
    cat(sprintf("    nrow up=%d rs=%d ncol up=%d rs=%d\n", NROW(a), NROW(b), NCOL(a), NCOL(b)))
    cat(sprintf("    cols up: %s\n", paste(colnames(a), collapse = ",")))
    cat(sprintf("    cols rs: %s\n", paste(colnames(b), collapse = ",")))
    for (col in intersect(colnames(a), colnames(b))) {
      ea <- a[[col]]; eb <- b[[col]]
      same <- identical(ea, eb)
      if (!same) {
        cat(sprintf("    col %s: class up=%s rs=%s", col,
                    paste(class(ea), collapse = "/"), paste(class(eb), collapse = "/")))
        if (is.factor(ea) || is.factor(eb)) {
          cat(sprintf(" levels up=[%s] rs=[%s]",
                      paste(utils::head(levels(ea), 8), collapse = ","),
                      paste(utils::head(levels(eb), 8), collapse = ",")))
        }
        if (is.numeric(ea) && is.numeric(eb) && length(ea) == length(eb)) {
          d <- abs(ea - eb); d <- d[!(is.na(ea) & is.na(eb))]
          if (length(d)) cat(sprintf(" max|diff|=%.3g n_differ=%d", max(d), sum(d > 0)))
        }
        cat("\n")
      }
    }
    only_up <- setdiff(rownames(a), rownames(b))
    only_rs <- setdiff(rownames(b), rownames(a))
    if (length(only_up) || length(only_rs)) {
      cat(sprintf("    rownames only-up=%d only-rs=%d\n", length(only_up), length(only_rs)))
    }
  }
}

## ---------------------------------------------------------------- 9 centrality
## igraph's ARPACK start vector is random, so both sides are seeded identically first. Same package
## version, same seed, same input: any remaining difference would be in what the shim hands igraph,
## which is what is under test. `set.seed` here seeds igraph's RNG, not the kernel's -- the kernel's
## own seed.use is an explicit argument, unchanged.
set.seed(20240501)
a9 <- q(UP$netAnalysis_computeCentrality(obj_up, slot.name = "netP"))
set.seed(20240501)
b9 <- q(cellchatrs::netAnalysis_computeCentrality(obj_rs, slot.name = "netP"))
obj_up <- a9
obj_rs <- b9
cmp("9 netAnalysis_computeCentrality netP (same igraph seed both sides)",
    identical(a9@netP, b9@netP))
if (!identical(a9@netP, b9@netP)) {
  na <- names(a9@netP$centr); nb <- names(b9@netP$centr)
  cat(sprintf("    pathways: up %d vs rs %d\n", length(na), length(nb)))
}

## ---------------------------------------------------------------- 10 whole-object verdict
## Everything above compared stage by stage so a failure points at its stage. This is the acceptance
## gate from the objective: `identical()` on the whole S4 object, minus the two fields that are
## wall-clock or visualisation by definition.
strip <- function(o) {
  o@options$run.time <- NULL
  o@images <- list()
  o@dr <- list()
  o
}
cmp("10 whole object identical() modulo run.time/images/dr",
    identical(strip(obj_up), strip(obj_rs)))

## The objective's recorded invariants, on the tutorial's own output.
p <- as.vector(obj_rs@net$pval); qv <- as.vector(obj_rs@net$prob)
nb <- obj_rs@options$parameter$nboot
cmp("Pval in {k/nboot}", all(abs(p * nb - round(p * nb)) < 1e-9))
cmp("Pval[Prob==0]==1", all(p[qv == 0] == 1))
cmp("dimnames consistent",
    identical(dimnames(obj_rs@net$prob)[[1]], levels(obj_rs@idents)) &&
    identical(dimnames(obj_rs@net$prob)[[2]], levels(obj_rs@idents)))

cat(sprintf("\n%s: %d failing comparisons out of %d\n",
            if (fails == 0L) "TUTORIAL-IDENTICAL" else "MISMATCH", fails, n_cmp))
quit(status = if (fails == 0L) 0L else 1L)
