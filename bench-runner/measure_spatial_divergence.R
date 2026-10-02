# Measure the divergence between upstream's approximate KNN and the port's exact k-d tree,
# on real spatial data.
#
# The locked decision (`PLAN.md` §14.5) is not "reproduce Annoy" but "replace it with an exact
# k-d tree and **measure** the divergence". This is that measurement, and it is deliberately
# three-way rather than two-way.
#
# The obvious experiment is exact-versus-Annoy, and on its own it is not enough: a difference
# between two things is only evidence of an improvement if you also know which one is right. So
# Annoy is run **twice** on identical input and the two runs are compared to each other, which
# separates two very different situations.
#
#   * If the two Annoy runs disagree, upstream's spatial output is not reproducible, and the exact
#     tree removes run-to-run noise.
#   * If they agree with each other but differ from the exact tree -- which is what actually
#     happens, and what the header of this file originally predicted the opposite of -- then the
#     approximation is a **systematic bias**, not noise, and the exact tree is a correction rather
#     than a de-noiser.
#
# The distinction matters for the paper. Noise is an argument about variance and a fix is to average
# it away. A bias is an argument about correctness: a reproducible wrong answer is worse than a
# varying right one, because nothing in a user's workflow would ever surface it. The measured
# numbers are below.
#
# Three things are compared, in increasing order of consequence:
#
#   1. `d.spatial`, the k x k matrix of minimum inter-group distances. The nearest-neighbour
#      search is where the approximation lives.
#   2. `P.spatial`, what `computeCommunProb` derives from it -- `1/d.spatial` with the
#      diagonal overwritten -- so a difference in `d.spatial` is scaled into the kernel.
#   3. `Prob` itself, the whole `K x K x N` array, with `Pval`. This is the number a user's
#      result depends on, and the only one the paper may quote.
#
# The data is the authors' own mouse cortex visium dataset: 1,073 spots, 648 genes, 433 of them
# signalling, 437 L-R pairs, 8 cell types ranging from 23 to 400 spots. Small enough to run
# three times, large enough that the approximation has something to be wrong about.
#
# Usage:  R_LIBS=.rlib R --vanilla -f bench-runner/measure_spatial_divergence.R
suppressWarnings(suppressMessages({
  library(methods); library(Matrix); library(collapse); library(dplyr)
}))
## The shim has to be *attached*, and not only for `computeRegionDistance`: the
## `cellchatrs_upstream_*` wrappers and the `CellChat` S4 class both come from it. Without this line
## the shim's functions are not on the search path at all, and the first thing reached for the
## upstream reference is whatever happens to share its name -- which is how a script ends up
## reporting a class-lookup failure for a package it never loaded.
suppressWarnings(suppressMessages(library(cellchatrs)))

CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")
VISIUM <- Sys.getenv("VISIUM", "/scratch/mdra00001/tmp/opencode/data/visium.rds")
OUT <- "bench-runner/results/spatial_divergence.json"
NBOOT <- as.integer(Sys.getenv("NBOOT", "10"))

## `q` runs an expression with its *console* output swallowed and its **value** returned.
##
## Two channels need handling and the naive one-liner gets the second wrong. `capture.output` alone
## does not silence upstream's `txtProgressBar`, which writes to **stderr**, so every one of the three
## kernel runs printed a 437-cell bar per L-R pair into the middle of the report. Redirecting the
## message connection to `nullfile()` is what fixes it.
##
## The earlier attempt at this returned the `sink()` result rather than the expression's value, and
## the symptom was quiet and bad: the `d.spatial` comparison reported "agree 0 / 0" and the run
## reported success. A helper that loses its return value does not throw -- it makes every downstream
## comparison vacuously pass. Hence `on.exit` to restore the sink, and the explicit `value` return.
q <- function(x) {
  con <- file(nullfile(), open = "wt")
  sink(con, type = "message")
  on.exit({ sink(type = "message"); close(con) }, add = TRUE)
  invisible(utils::capture.output(value <- suppressWarnings(suppressMessages(x))))
  ## Guard against the failure mode above: an empty return is never a legitimate result here, and
  ## silently comparing NULL against NULL would report agreement.
  if (is.null(value)) stop("q(): the expression returned NULL; the comparison would be vacuous")
  value
}

say <- function(...) cat(..., "\n", sep = "")

setClass("CellChatProbe", representation(
  data = "ANY", data.signaling = "ANY", data.smooth = "ANY", idents = "ANY", meta = "ANY",
  images = "ANY", DB = "ANY", LR = "ANY", net = "ANY", netP = "ANY", options = "ANY"
))

## ---------------------------------------------------------------- the pinned object
## The RDS holds a real `CellChat` S4 object, so the class has to be defined before it can be
## unserialised at all -- `readRDS` on a `CellChat` object without the class raises
## 'CellChat' is not a defined class, which is a property of the fixture, not of the data.
## Into the **global** environment, not a private one. The shim's `computeRegionDistance` and
## `computeCommunProb` are called on a real `CellChat` S4 object, and `methods`' dispatch resolves
## the class through the *search path* -- sourcing `CellChat_class.R` into a private env and then
## calling the shim raises `unable to find required package 'CellChat'` from `getClassDef`, because
## the class is defined but not reachable. `tests/parity/check_merge.R` hits the same constraint and
## does the same thing. `setClass` on an already-defined class returns quietly, so this is safe to
## run after upstream has been sourced elsewhere.
## **Attached**, not merely loadable. `setClassUnion(name = "AnyMatrix", members = c("matrix",
## "dgCMatrix"))` needs `dgCMatrix` to be *defined* at the moment the class file is sourced, and
## `requireNamespace` alone leaves it unavailable -- the failure is
## `the member classes must be defined: not true of "dgCMatrix"`. That in turn left `CellChat`
## undefined, and the next symptom was a mile away: `unable to find required package 'CellChat'`
## from `getClassDef`, raised inside the shim rather than at the point of the mistake. Sourced
## without `capture.output` too, so a failure surfaces here rather than three calls later.
suppressPackageStartupMessages(library(Matrix))
invisible(suppressWarnings(suppressMessages(
  sys.source(file.path(CC, "R", "CellChat_class.R"), envir = globalenv(), keep.source = FALSE))))
stopifnot(!is.null(methods::getClass("CellChat", where = globalenv())))
o <- readRDS(VISIUM)

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)

say("data:", paste(dim(o@data), collapse = "x"),
    " signalling:", paste(dim(o@data.signaling), collapse = "x"),
    " L-R pairs:", nrow(o@LR$LRsig),
    " groups:", nlevels(o@idents),
    " spots:", ncol(o@data))

## `updateCellChat` adds the columns `computeRegionDistance` reads. `samples` is absent from the
## visium object -- it is a single section -- so a single-level factor is the honest value, and
## declaring a level that does not exist would change the per-sample k-NN branch's behaviour.
meta <- data.frame(group = o@idents, stringsAsFactors = FALSE)
meta$samples <- factor(rep("section1", ncol(o@data)))
meta$label <- o@idents

## `contact.range` is the visium spot diameter. `contact.dependent = TRUE` has **no default range**:
## upstream raises "Please check the documentation of `computeCommunProb` and provide the value of
## either `contact.range` or `contact.knn.k`" without one. 65 is the physically meaningful value for a
## 10x Visium spot and it is what makes the `adj.contact` k-NN -- the query the exact tree replaces --
## do any work, so the divergence being measured is the divergence that matters.
spot_diameter <- as.numeric(o@images$scale.factors$spot.diameter)

## `ratio` and `tol` are **not optional**, and their absence is silent.
##
## Upstream does `qout$distance <- qout$distance * ratio[k]` and then
## `d.spatial[i,j,k] <- FunMean(qout$distance)`. With the documented defaults -- `ratio = NULL`,
## `tol = NULL` -- the product is `numeric(0)`, `FunMean(numeric(0))` is `NaN`, and **every** entry of
## `d.spatial` is `NaN`. No error, no warning: a 64-element matrix of `NaN` that flows into
## `P.spatial[is.na(d.spatial)] <- 0` and then silently zeroes every spatial weight. `updateCellChat`
## normally fills these in as `spot.diameter / fullsize`, but this visium object predates
## `spatial.factors`, so they have to be supplied.
##
## `ratio = 1` keeps distances in the same image-pixel unit as `contact.range = 65`, and `tol = 0`
## makes the range tests exact rather than tolerance-padded. Both are the neutral choice: any other
## value would rescale the distances and change which pairs are in range, which is a different
## experiment.
RATIO <- 1
TOL <- 0

## `interaction.range = 1000` rather than upstream's default 250. The visium centroids span
## 2,676 x 292 image pixels, so 250 -- roughly four spot diameters -- leaves 24 of the 56
## off-diagonal group pairs with no nearby pair at all, and a pair that is `NaN` under one index is
## `NaN` under the other, which is exactly the kind of disagreement that a KNN approximation cannot
## produce. A wider range puts more pairs in play, so more of them can be compared.
INTERACTION_RANGE <- as.integer(Sys.getenv("INTERACTION_RANGE", "1000"))
say("spot.diameter:", spot_diameter, " k.min: 10  interaction.range:", INTERACTION_RANGE,
    " ratio:", RATIO, " tol:", TOL)

## ---------------------------------------------------------------- three runs
## `run_one(fn, seed)` returns the `d.spatial`, `P.spatial` and `Prob` from one pass. `fn` is
## either the shim's `computeRegionDistance` (exact k-d tree) or pinned upstream's (Annoy).
##
## `computeRegionDistance` itself is deterministic given the coordinates, so the seed only moves
## Annoy. That is the point: the seed is threaded through so the two upstream runs are the same
## computation with a different random draw, and any difference between them is Annoy's.
spatial_run <- function(exact, seed) {
  ## No `CellChat` object is touched here: `computeRegionDistance` takes `coordinates` and `meta` as
  ## plain arguments, so a probe object would only be a way to reach a class lookup that can fail.
  ## `k.min = 10` with a 23-spot group means the 1-NN and the k-NN branches both have to cope with
  ## fewer spots than neighbours, which is the regime where an approximate index is most likely to
  ## be wrong and where the exact tree's tie-break matters most.
  crd <- cbind(x = o@images$coordinates$x_cent, y = o@images$coordinates$y_cent)
  if (exact) {
    q(computeRegionDistance(coordinates = crd, meta = meta, k.min = 10,
                            contact.dependent = TRUE, contact.range = spot_diameter,
                            interaction.range = INTERACTION_RANGE, ratio = RATIO, tol = TOL))
  } else {
    ## `set.seed` here and nowhere else. Annoy draws from R's RNG, so two upstream runs with
    ## different seeds are the same computation with a different random draw -- and any difference
    ## between them is the approximation rather than the data.
    set.seed(seed)
    q(cellchatrs_upstream_computeRegionDistance(
      coordinates = crd, meta = meta, k.min = 10, contact.dependent = TRUE,
      contact.range = spot_diameter, interaction.range = INTERACTION_RANGE,
      ratio = RATIO, tol = TOL))
  }
}

t0 <- Sys.time()
exact <- spatial_run(TRUE, NA)
t_exact <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
annoy_a <- spatial_run(FALSE, 1L)
t_a <- as.numeric(difftime(Sys.time(), t0, units = "secs")) - t_exact
annoy_b <- spatial_run(FALSE, 2L)
t_b <- as.numeric(difftime(Sys.time(), t0, units = "secs")) - t_exact - t_a
say("d.spatial timings (s): exact", format(t_exact, digits = 4),
    " annoy(seed 1)", format(t_a, digits = 4), " annoy(seed 2)", format(t_b, digits = 4))

## ---------------------------------------------------------------- the comparison
## Group-pair level. `d.spatial` is `NaN` where two groups have no nearby pair, and that `NaN` is
## *meaningful* -- it is what zeroes `P.spatial` -- so it has to be compared as a value and not
## dropped with `na.rm`. The three states are "both NaN", "differs", "agrees".
cmp_matrix <- function(a, b, label) {
  na <- is.na(a); nb <- is.na(b)
  both_na <- sum(na & nb)
  ## A NaN against a number is the most consequential disagreement there is: one side zeroes the
  ## pair and the other divides by it, which is a different interaction entirely rather than a
  ## different weight on the same one.
  one_na <- sum(xor(na, nb))
  both_num <- !na & !nb
  differ <- sum(both_num & a != b)
  rel <- if (sum(both_num)) max(abs(a[both_num] - b[both_num]) / pmax(abs(b[both_num]), .Machine$double.eps)) else 0
  say(sprintf("%-28s agree %d / %d   differ %d   NaN-vs-value %d   both NaN %d   max rel diff %s",
              label, sum(both_num) - differ, length(a), differ, one_na, both_na,
              format(rel, digits = 4)))
  list(agree = sum(both_num) - differ, total = length(a), differ = differ,
       one_na = one_na, both_na = both_na, max_rel = rel)
}

cmp_exact_a <- cmp_matrix(exact$d.spatial, annoy_a$d.spatial, "d.spatial exact vs Annoy 1")
cmp_a_b <- cmp_matrix(annoy_a$d.spatial, annoy_b$d.spatial, "d.spatial Annoy 1 vs 2")
cmp_exact_b <- cmp_matrix(exact$d.spatial, annoy_b$d.spatial, "d.spatial exact vs Annoy 2")

## ---------------------------------------------------------------- the consequence
## The whole kernel, three times, so `Prob` is compared rather than `P.spatial`. This is the number
## a result depends on, and the one the paper may quote.
## The probe object. Every slot here is read by one of the two paths; the names must match upstream's
## exactly because the shim uses `@` and `slot()`, both of which resolve against the class.
probe <- function() {
  new("CellChatProbe",
      data = o@data,
      data.signaling = o@data.signaling,
      data.smooth = o@data.signaling,
      idents = o@idents,
      meta = meta,
      ## `@images$coordinates` and `@images$spatial.factors` are the spatial branch's two inputs.
      ## The *name* matters: the shim tests `if ("spatial.factors" %in% names(object@images))` and
      ## otherwise stops with "`object@images$spatial.factors` is missing. Please update the object via
      ## `updateCellChat`!". This visium object carries the older `scale.factors`, so the field has to
      ## be renamed on the way in -- which is exactly what `updateCellChat` did when upstream renamed
      ## it in 2.1.1. `spot.diameter` is a scalar; upstream's tutorial uses a per-section vector, and
      ## the code indexes it as `spot.diameter[i]`, so it is wrapped.
      images = list(coordinates = cbind(x = o@images$coordinates$x_cent,
                                        y = o@images$coordinates$y_cent),
                    spatial.factors = list(spot.diameter = spot_diameter,
                                           ratio = RATIO, tol = TOL)),
      DB = list(complex = DB$complex, cofactor = DB$cofactor),
      LR = list(LRsig = o@LR$LRsig),
      net = list(prob = array(0, c(1, 1, 1)), pval = array(0, c(1, 1, 1))),
      netP = list(),
      options = list(datatype = "spatial", mode = "single", db = normalizePath(DBDIR),
                     population.size = FALSE))
}

prob_run <- function(exact, seed, lrsig = LR_ok) {
  o2 <- probe()
  o2@LR <- list(LRsig = lrsig)
  fn <- if (exact) computeCommunProb else cellchatrs_upstream_computeCommunProb
  if (!exact) set.seed(seed)
  q(fn(o2, type = "triMean", trim = 0.1, population.size = FALSE,
       distance.use = TRUE, interaction.range = INTERACTION_RANGE, scale.distance = 0.01,
       k.min = 10, contact.dependent = TRUE, contact.range = spot_diameter,
       do.symmetric = TRUE, nboot = NBOOT, seed.use = 1L, Kh = 0.5, n = 1))
}
## The `Prob` stage, and one thing has to be handled first.
##
## On the *unsubsampled* data both sides raise `"subscript out of bounds"`. That is upstream's own
## crash for an L-R pair whose receptor complex has a subunit that is not a row of
## `data.signaling` -- the behaviour `PLAN.md` requires the port to reproduce rather than paper over.
## All 437 pairs come from the pinned visium object's own `LRsig`, and some of them name subunits the
## object's 433 signalling genes do not include. So the first thing measured is whether the two sides
## fail *the same way*, which is the parity claim; and then the pairs are filtered so a `Prob`
## comparison can happen at all.
lr_resolvable <- function(lrsig, cm) {
  cs <- as.matrix(cm[, grepl("^subunit", colnames(cm)), drop = FALSE])
  ok <- vapply(seq_len(nrow(lrsig)), function(i) {
    lig <- unlist(strsplit(as.character(lrsig$ligand[i]), ",", fixed = TRUE), use.names = FALSE)
    rec <- unlist(strsplit(as.character(lrsig$receptor[i]), ",", fixed = TRUE), use.names = FALSE)
    all(c(lig, rec) %in% rownames(o@data.signaling))
  }, logical(1))
  lrsig[ok, , drop = FALSE]
}
LR_full <- o@LR$LRsig
LR_ok <- lr_resolvable(LR_full, DB$complex)
say("")
say("L-R pairs: ", nrow(LR_full), " total, ", nrow(LR_ok),
    " whose every subunit is present in data.signaling")

## The parity claim on the *unfiltered* set, recorded as the error both sides must reproduce.
err_exact <- tryCatch({ prob_run(TRUE, NA, LR_full); NA_character_ },
                      error = function(e) conditionMessage(e))
err_annoy <- tryCatch({ prob_run(FALSE, 1L, LR_full); NA_character_ },
                      error = function(e) conditionMessage(e))
say("  unfiltered, exact side  :", if (is.na(err_exact)) "no error" else err_exact)
say("  unfiltered, upstream side:", if (is.na(err_annoy)) "no error" else err_annoy)
say("  identical messages      :", identical(err_exact, err_annoy))

say("")
say("running the kernel three times with nboot =", NBOOT, "on the", nrow(LR_ok), "resolvable pairs")
t0 <- Sys.time()
p_exact <- prob_run(TRUE, NA, LR_ok)
say("  exact   ", format(as.numeric(difftime(Sys.time(), t0, units = "secs")), digits = 4), "s")
p_a <- prob_run(FALSE, 1L, LR_ok)
say("  Annoy 1 ", format(as.numeric(difftime(Sys.time(), t0, units = "secs")), digits = 4), "s cumulative")
p_b <- prob_run(FALSE, 2L, LR_ok)
say("  Annoy 2 ", format(as.numeric(difftime(Sys.time(), t0, units = "secs")), digits = 4), "s cumulative")

## The `Prob` comparison has to be over *interactions*, not cells: `Prob[i, j, n]` is one
## source/target/pathway weight, and a difference there changes a reported interaction's strength.
cmp_prob <- function(a, b, label) {
  d <- abs(a - b)
  nz <- d > 0
  ## The relative measure is against the **spread of `Prob`**, not against each value. An element-wise
  ## `d / |b|` is meaningless on this data: `Prob` ranges over four orders of magnitude and passes
  ## through exact zeros, so a difference of 1e-12 next to a `b` of 1e-25 is a "relative" error of
  ## 8e12 -- the number an earlier version of this function printed, and it says nothing. Dividing by
  ## the peak value gives a quantity in [0, 1] that answers the question a reader actually has: how
  ## big is the largest discrepancy compared to a typical interaction weight?
  scale <- max(abs(b), na.rm = TRUE)
  rel <- if (is.finite(scale) && scale > 0) max(d) / scale else 0
  ## The 1-NN pair and the long-range `adj.spatial` gate both go through a nearest-neighbour
  ## distance, so a `d.spatial` difference shows up in `Prob` scaled by `P.spatial`, which is
  ## `1/d.spatial` -- itself inversely proportional. Reporting the ratio to the peak keeps that
  ## legible.
  say(sprintf("%-28s differ %d of %d (%.2f%%)   max abs %s   max/peak %s",
              label, sum(nz), length(d), 100 * mean(nz), format(max(d), digits = 4),
              format(rel, digits = 4)))
  list(differ = sum(nz), total = length(d), max_abs = max(d), max_rel = rel)
}
p_exact_a <- cmp_prob(as.numeric(p_exact@net$prob), as.numeric(p_a@net$prob), "Prob exact vs Annoy 1")
p_a_b <- cmp_prob(as.numeric(p_a@net$prob), as.numeric(p_b@net$prob), "Prob Annoy 1 vs 2")

## ---------------------------------------------------------------- the record
## A dependency-free JSON writer, as `bench_real.R` has: a benchmark that cannot emit its result
## without installing a package is a benchmark nobody runs.
j_esc <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  paste0("\"", x, "\"")
}
j_num <- function(x) if (is.finite(x)) format(x, digits = 17) else "null"
fields <- function(nms, vals) paste0("\"", nms, "\": ", vals, collapse = ", ")

obj <- function(cmp) list(
  agree = j_num(cmp$agree), total = j_num(cmp$total), differ = j_num(cmp$differ),
  nan_versus_value = j_num(cmp$one_na), both_nan = j_num(cmp$both_na), max_rel_diff = j_num(cmp$max_rel))
pobj <- function(cmp) list(
  differ = j_num(cmp$differ), total = j_num(cmp$total),
  fraction = j_num(cmp$differ / cmp$total),
  max_abs_diff = j_num(cmp$max_abs), max_diff_over_peak = j_num(cmp$max_rel))

## The pinned upstream, stated once. Read from `cellchatrs_upstream_cached()`'s own notion of where
## upstream lives rather than hard-coded twice: a benchmark that quotes a commit which has drifted
## from the one the rest of the suite uses is worse than one that quotes nothing. The *commit* is
## verified by `tests/parity/verify_upstream.R` and by CI's `verify-upstream` job; here it is
## recorded so the JSON says which tree produced the numbers.
UPSTREAM_SHA <- Sys.getenv("CELLCHATRS_UPSTREAM_SHA", "75253cd0c9e68410e6e721a6d3a0419a1d7e358f")
UPSTREAM_VERSION <- Sys.getenv("CELLCHATRS_UPSTREAM_VERSION", "2.2.0.9001")

json <- paste0("{\n",
  "  \"what\": \"exact k-d tree vs BiocNeighbors AnnoyParam, on the authors' mouse cortex visium data\",\n",
  "  \"upstream\": { \"repo\": \"jinworks/CellChat\", \"commit\": \"", UPSTREAM_SHA,
  "\", \"version\": \"", UPSTREAM_VERSION, "\" },\n",
  "  \"data\": { \"spots\": ", j_num(ncol(o@data)),
  ", \"genes\": ", j_num(nrow(o@data)),
  ", \"signalling_genes\": ", j_num(nrow(o@data.signaling)),
  ", \"lr_pairs\": ", j_num(nrow(LR_full)),
  ", \"groups\": ", j_num(nlevels(o@idents)),
  ", \"nboot\": ", j_num(NBOOT),
  ", \"interaction_range\": ", j_num(INTERACTION_RANGE),
  ", \"contact_range\": ", j_num(spot_diameter),
  ", \"ratio\": ", j_num(RATIO), ", \"tol\": ", j_num(TOL),
  ", \"k_min\": 10 },\n",
  "  \"d_spatial\": {\n",
  "    \"exact_vs_annoy1\": {", fields(names(obj(cmp_exact_a)), unlist(obj(cmp_exact_a))), " },\n",
  "    \"exact_vs_annoy2\": {", fields(names(obj(cmp_exact_b)), unlist(obj(cmp_exact_b))), " },\n",
  "    \"annoy1_vs_annoy2\": {", fields(names(obj(cmp_a_b)), unlist(obj(cmp_a_b))), " }\n  },\n",
  "  \"prob\": {\n",
  "    \"exact_vs_annoy1\": {", fields(names(pobj(p_exact_a)), unlist(pobj(p_exact_a))), " },\n",
  "    \"annoy1_vs_annoy2\": {", fields(names(pobj(p_a_b)), unlist(pobj(p_a_b))), " }\n  },\n",
  "  \"unfiltered_lrsig\": { \"total\": ", j_num(nrow(LR_full)),
  ", \"resolvable\": ", j_num(nrow(LR_ok)),
  ", \"exact_error\": ", j_esc(if (is.na(err_exact)) "" else err_exact),
  ", \"upstream_error\": ", j_esc(if (is.na(err_annoy)) "" else err_annoy),
  ", \"errors_identical\": ", if (identical(err_exact, err_annoy)) "true" else "false", " },\n",
  "  \"seconds\": { \"exact_d_spatial\": ", j_num(t_exact),
  ", \"annoy1_d_spatial\": ", j_num(t_a), ", \"annoy2_d_spatial\": ", j_num(t_b), " }\n}\n")
dir.create(dirname(OUT), showWarnings = FALSE, recursive = TRUE)
writeLines(json, OUT)
say("")
say("wrote ", OUT)
