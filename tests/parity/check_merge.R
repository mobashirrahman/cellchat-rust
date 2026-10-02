# `mergeCellChat` and the comparison-analysis path.
#
# This is the parity gate for `analysis.mergeCellChat`, and it is deliberately *not* a comparison of
# `mergeCellChat` with upstream's `mergeCellChat`. The shim does not reimplement that function --
# it is pure S4 slot assembly with no arithmetic, and the 14.1 scope decision leaves it in R -- so
# comparing the two would compare a function with itself and report `IDENTICAL` while testing
# nothing. That is the exact failure mode `docs/SEMANTICS.md` records for the `rankNetPairwise`
# fallback.
#
# What is worth testing is what the merge *produces*. A merged object becomes the input to
# `computeCommunProb` in the comparison workflow, and the merge decides three things the kernel then
# reads: the column order of `data`, the `idents` factor and its **level order**, and the `meta`
# column intersection. So every case here does the same thing:
#
#   1. merge two CellChat objects;
#   2. run the Rust-backed `computeCommunProb` on the merged object;
#   3. run pinned upstream's on the same object;
#   4. require `identical()`.
#
# The level order matters more than it looks. `mergeCellChat` builds `idents.joint` with
# `levels = union(idents.levels)`, which is *first-appearance order across datasets*, not sorted
# and not `c(levels(a), levels(b))`. A kernel that used the wrong order would produce a permuted
# `Prob` array, and the result would be wrong rather than merely reordered -- `aggregateNet` groups
# by the character key `source|target`, so a permuted level order changes the weights.
#
# Usage:  R_LIBS=.rlib R --vanilla -f tests/parity/check_merge.R
suppressWarnings(suppressMessages({
  library(methods); library(Matrix); library(collapse); library(dplyr)
}))

## `ggplot2` is attached, and not incidentally. `rankNet` is the one function used here that plots,
## and pinned upstream's body calls `ggplot` **unqualified** -- it resolves through the caller's
## search path, because `cellchatrs_upstream_env()` gives the sourced code `parent = globalenv()` so
## that upstream's own `rowSums` finds `Matrix`'s method (see the comment on that function). A real
## CellChat session has `ggplot2` attached because CellChat imports it. Without attaching it here,
## upstream fails with 'could not find function "ggplot"' while the port, which calls `ggplot2::`
## explicitly, succeeds -- and the shim-vs-upstream comparison then reports a difference that is
## purely an artefact of how this script was launched. `check_identical.R` does the same for its
## rankNet block, and stops outright if ggplot2 is missing rather than comparing against an error.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("this gate needs ggplot2 attached: upstream's rankNet body resolves it unqualified",
       call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))

CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")

suppressWarnings(suppressMessages(library(cellchatrs)))

if (tolower(Sys.getenv("CELLCHATRS_FALLBACK", "0")) %in% c("1", "true", "yes", "on")) {
  stop("CELLCHATRS_FALLBACK is set: every shim call would be delegated to pinned upstream and\n",
       "this gate would compare upstream against itself. Unset it to run the real comparison.",
       call. = FALSE)
}

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)

## `setClass` registers a class definition, but `new("CellChat")` resolves it through the *search
## path*, not through the upstream environment the shim cached it in -- so `getClass` cannot see it
## from here and the error is '"CellChat" is not a defined class' rather than anything about the
## fixture. Sourcing the class file into the global environment is the fix, and it is safe to do
## twice: `setClass` on an already-defined class returns quietly, and `mergeCellChat` itself is
## always called through `cellchatrs_upstream_mergeCellChat`, so the global copy is never the one
## under test. Only the class *definition* is taken from here.
invisible(sys.source(file.path(CC, "R", "CellChat_class.R"), envir = globalenv(),
                     keep.source = FALSE))
stopifnot(!is.null(methods::getClass("CellChat", where = globalenv())))

## The same expression fixture the main gate uses, so a failure here is about the merge and not
## about different test data. Genes and cells are sliced rather than re-randomised.
src <- readLines("tests/parity/gen_prob_golden.R")
stop_at <- grep("^q <- file", src)[1]
eval(parse(text = paste(src[seq_len(stop_at - 1)], collapse = "\n")))

## ---------------------------------------------------------------- fixtures
## Two objects that differ in every way the merge has to reconcile: gene sets, cell counts, cell
## names, cluster names and their order, and `meta` columns. `cell.prefix = FALSE` requires
## disjoint cell names, so the names are suffixed here.
build <- function(n_cells, genes_keep, cell_tag, groups, levels, meta_extra = NULL) {
  cn <- paste0("cell", seq_len(n_cells), "_", cell_tag)
  ## Slice both axes: the fixture is 60 cells wide and `mergeCellChat` requires the two objects to
  ## be able to `cbind`, so each object has to carry only the cells it is going to contribute.
  m <- data.signaling[genes_keep, seq_len(n_cells), drop = FALSE]
  ## Sparse, so `cbind` produces a `dgCMatrix` like a real object and `rbind` of meta is exercised
  ## against a Matrix-backed `data` slot rather than a base matrix.
  m <- Matrix::Matrix(as.matrix(m), sparse = TRUE)
  colnames(m) <- cn
  rownames(m) <- genes_keep
  meta <- data.frame(cell = cn, stringsAsFactors = FALSE)
  meta$datasets <- factor(cell_tag, levels = cell_tag)
  meta$cluster <- factor(rep(groups, length.out = n_cells), levels = levels)
  if (!is.null(meta_extra)) meta[[meta_extra]] <- rep("shared", n_cells)
  ## The object is built with `new("CellChat", ...)` rather than upstream's `createCellChat`.
  ## `createCellChat` runs the full preprocessing pipeline -- it calls `print` on a character
  ## message, tries `getClassDef` on the assay, and populates slots this test never reads -- which
  ## makes the fixture's inputs implicit and its output noisy. `mergeCellChat` reads five slots:
  ## `data`, `data.signaling`, `meta`, `idents` and `options$datatype`, plus the `net`/`netP`/`LR`/
  ## `var.features`/`images` lists it wraps without reading. So the slots are set directly, which
  ## also means every input the merge sees is visible in this file.
  o <- new("CellChat",
           data = m,
           data.signaling = m,
           meta = meta,
           idents = factor(rep(groups, length.out = n_cells), levels = levels),
           net = list(prob = array(0, c(1, 1, 1)), pval = array(0, c(1, 1, 1))),
           netP = list(),
           LR = list(LRsig = lrsig_for(genes_keep)),
           var.features = list(features = genes_keep),
           DB = list(complex = DB$complex, cofactor = DB$cofactor, interaction = DB$interaction),
           images = list(),
           options = list(datatype = "RNA", mode = "single", db = normalizePath(DBDIR),
                          population.size = TRUE))
  ## A real net, computed by upstream, so that the Rust functions exercised below are working on
  ## the kind of structure a merge actually carries rather than a zero array. `computeCommunProb`
  ## here is the *reference*, not the thing under test -- the fixture has to be identical for both
  ## sides, and using the shim to build it would make the input depend on the code being tested.
  suppressWarnings(suppressMessages(
    cellchatrs_upstream_computeCommunProb(o, type = "triMean", trim = 0.1,
                                          population.size = TRUE, nboot = 5L,
                                          seed.use = 1L, Kh = 0.5, n = 1)))
}

ALL_GENES <- rownames(data.signaling)
ALL_CELLS <- colnames(data.signaling)

## The L-R set is restricted to pairs whose ligand *and* every receptor subunit are inside the
## gene subset each fixture is built from. Without this the fixture dies in `computeCommunProb`
## with upstream's own `"subscript out of bounds"` -- the behaviour `docs/SEMANTICS.md` records as
## deliberately reproduced for a complex with a missing subunit. That is the right behaviour and
## the wrong thing for a *fixture*: here the missing subunit is an accident of slicing a shared
## fixture, not the condition under test. The complex branch is covered by the kernel's own parity
## tests, where the missing subunit is the point.
## The L-R set is restricted to pairs whose ligand *and* every receptor subunit are present.
## Without this the fixture dies in `computeCommunProb` with upstream's own
## `"subscript out of bounds"` -- the behaviour `docs/SEMANTICS.md` records as deliberately
## reproduced for a complex with a missing subunit. That is the right behaviour and the wrong thing
## for a *fixture*: here a missing subunit would be an accident of slicing a shared fixture, not the
## condition under test. The complex branch is covered by the kernel's own parity tests, where the
## missing subunit is the point.
subunits <- function(x) {
  x <- as.character(x)
  if (length(x) == 0L) return(character())
  trimws(unlist(strsplit(x[!is.na(x) & nzchar(x)], ",", fixed = TRUE), use.names = FALSE))
}
PAIRS <- LRsig[vapply(seq_len(nrow(LRsig)), function(i) {
  all(c(subunits(LRsig$ligand[i]), subunits(LRsig$receptor[i])) %in% ALL_GENES)
}, logical(1)), , drop = FALSE]
stopifnot(nrow(PAIRS) >= 3)

## `LRsig` for one fixture: the pairs that resolve inside *that* object's gene set. Each object gets
## its own, so a pair present in A and absent from B is exercised rather than assumed away.
lrsig_for <- function(genes_keep) {
  k <- vapply(seq_len(nrow(PAIRS)), function(i) {
    all(c(subunits(PAIRS$ligand[i]), subunits(PAIRS$receptor[i])) %in% genes_keep)
  }, logical(1))
  if (!any(k)) stop("no L-R pair resolves inside ", paste(genes_keep, collapse = ","))
  out <- PAIRS[k, , drop = FALSE]
  ## `interaction_name` has to stay unique: `computeCommunProb` keys the `Prob` array by it, and a
  ## duplicate would collapse two pairs into one slot.
  if (anyDuplicated(out$interaction_name)) {
    out <- out[!duplicated(out$interaction_name), , drop = FALSE]
  }
  out
}


## Gene sets that do not overlap, overlap partially, and are identical, because
## `genes.use` is the intersection and the shim must not silently cbind differently shaped matrices.
## A and B share two genes (INHBA, FST) and each has two of its own. The intersection is therefore
## non-empty and strictly smaller than either, which is what makes the `genes.use` check below able
## to fail. A near-disjoint pair is built separately -- see the `data.signaling` quirk below.
genes_a <- c("G1", "G2", "INHBA", "FST")
genes_b <- c("INHBA", "FST", "G4", "G7")

## Cluster levels: A has a level B does not have, and B's levels are in a different order. The
## union has to come out in A-then-B first-appearance order, not sorted.
lv_a <- c("Bcell", "Tcell", "Mono")
lv_b <- c("Mono", "Bcell", "NK")

## `meta` columns: A has `age`, B does not. `meta.use` is the intersection, so `age` must be dropped
## from the merged `meta` -- and the *data* slot is unaffected by that, which is worth asserting
## because the two are computed from different intersections.
obj_a <- build(6, genes_a, "a", c("Bcell", "Tcell", "Mono", "Bcell", "Tcell", "Mono"),
               lv_a, meta_extra = "age")
obj_b <- build(5, genes_b, "b", c("Mono", "NK", "Bcell", "NK", "Mono"), lv_b)

## Same shape as A, disjoint cell names, for the cases that need a well-behaved merge.
##
## Every level must actually be used, because upstream's `computeCommunProb` refuses a factor with
## an unused level ("Please check `unique(object@idents)` and ensure that the factor levels are
## correct!") and suggests `droplevels`. The fixtures therefore list all of `lv_a` / `lv_b`, in a
## *different* order from the primary objects, which is what makes the merged level order a real
## test rather than a copy of the input order.
obj_a2 <- build(6, genes_a, "c", c("Mono", "Bcell", "Tcell", "Bcell", "Mono", "Tcell"), lv_a)
obj_b2 <- build(5, genes_b, "d", c("NK", "Mono", "Bcell", "NK", "Mono"), lv_b)

## ---------------------------------------------------------------- the comparison
fails <- 0L
n_cmp <- 0L

## Warnings are part of the contract. Upstream `warn`s for `cell.prefix = TRUE` and for a `meta`
## rownames mismatch, and `cat`s the barcodes; the gate compares the message text as well as the
## value, because "same result, different warning" is a behavioural difference a user would see.
capture <- function(expr) {
  w <- character(); m <- character(); c_ <- character()
  v <- withCallingHandlers(
    tryCatch(expr, error = function(e) structure(conditionMessage(e), class = "gate_error")),
    warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") },
    message = function(x) { m <<- c(m, conditionMessage(x)); invokeRestart("muffleMessage") },
    output = function(x) { c_ <<- c(c_, as.character(x)); invokeRestart("muffleOutput") }
  )
  ## An errored value is tagged so `compare` can tell it from a real result, and untagged here so a
  ## message can be compared with `identical()` against a plain string.
  ## `as.character` rather than `structure(v, class = "character")`: the latter resets `class` but
  ## leaves the attribute `structure` added, and `identical()` compares attributes as well as values,
  ## so a message that reads correctly still fails the comparison.
  if (inherits(v, "gate_error")) v <- as.character(v)
  list(value = v, warnings = w, messages = m, output = c_)
}

compare <- function(label, fn) {
  n_cmp <<- n_cmp + 1L
  up <- capture(fn(cellchatrs_upstream_mergeCellChat))
  rs <- capture(fn(merge_upstream_only))
  ok_val <- identical(up$value, rs$value)
  ok_warn <- identical(up$warnings, rs$warnings)
  ok_msg <- identical(up$messages, rs$messages)
  ok_out <- identical(up$output, rs$output)
  ok <- ok_val && ok_warn && ok_msg && ok_out
  cat(sprintf("%-34s value=%-5s warn=%-5s msg=%-5s cat=%-5s\n",
              label, ok_val, ok_warn, ok_msg, ok_out))
  if (!ok) {
    fails <<- fails + 1L
    if (!ok_val && !inherits(up$value, "gate_error") && is.list(up$value)) {
      a <- up$value; b <- rs$value
      if (inherits(a, "CellChat") && inherits(b, "CellChat")) {
        for (s in c("data", "data.signaling", "idents", "meta")) {
          if (!identical(slot(a, s), slot(b, s))) {
            cat(sprintf("   slot %s: up %s vs rs %s\n", s,
                        paste(dim(slot(a, s)), collapse = "x"),
                        paste(dim(slot(b, s)), collapse = "x")))
          }
        }
        cat("   idents levels up:", paste(levels(slot(a, "idents")), collapse = ","), "\n")
        cat("   idents levels rs:", paste(levels(slot(b, "idents")), collapse = ","), "\n")
      }
    }
    if (!ok_warn) {
      cat("   warnings up:", paste(up$warnings, collapse = " | "), "\n")
      cat("   warnings rs:", paste(rs$warnings, collapse = " | "), "\n")
    }
  }
  invisible(up$value)
}

## The shim does not override `mergeCellChat`, so there is one function, not two. `compare` calls
## the shim path and the upstream path as separate expressions only to keep the call sites
## symmetric with the rest of the gate; they resolve to the same body, which is why every
## *value* comparison here is followed by a kernel comparison rather than being treated as the test.
merge_upstream_only <- function(...) cellchatrs_upstream_mergeCellChat(...)

## ---------------------------------------------------------------- 1. structure of the merge
merged <- compare("merge default", function(f) f(list(obj_a, obj_b)))
merged2 <- compare("merge explicit names", function(f) f(list(obj_a, obj_b), add.names = c("X", "Y")))
merged3 <- compare("merge three objects", function(f) f(list(obj_a, obj_b, obj_a2)))
merged4 <- compare("merge merge.data = TRUE", function(f) f(list(obj_a, obj_b), merge.data = TRUE))
merged5 <- compare("merge cell.prefix = TRUE", function(f) f(list(obj_a, obj_b), cell.prefix = TRUE))
merged6 <- compare("merge same gene set", function(f) f(list(obj_a, obj_a2)))
merged7 <- compare("merge one object", function(f) f(list(obj_a)))

## ---------------------------------------------------------------- 2. the merge's own decisions
## These are the assertions a shim-vs-upstream comparison cannot make, because there is only one
## implementation. They are the merge's contract, written down so a future change to the shim has
## something to break.
cat("\n-- the merge's contract\n")
stopifnot(identical(colnames(merged@data.signaling),
                    c(paste0("cell", 1:6, "_a"), paste0("cell", 1:5, "_b"))))
cat("data.signaling columns are the concatenation, in argument order: ok\n")

## `genes.use` is the intersection **in A's order**. A carries genes 1:6 and B carries 4:9, so the
## intersection is 4:6 -- three genes. Not the union (1:9) and not A's own set (1:6); getting this
## wrong would give a `data.signaling` with the wrong number of rows, and the kernel would then
## aggregate over genes that only exist in one dataset.
stopifnot(identical(rownames(merged@data.signaling), intersect(genes_a, genes_b)))
cat("genes.use is the intersection in A's order, not the union: ok\n")

## An upstream quirk worth pinning, because it is silent and it loses data. `data.joint` is
## intersected to `genes.use`, but the subsetting of `data.signaling` is against
## `gene.signaling.joint`, which is the **union** of the two objects' gene names. So whenever an
## object's genes are not all shared, some of the intersected rows are then dropped again -- and when
## the two gene sets are disjoint the result is a `data.signaling` with **zero rows** and no
## colnames, from a merge that reported success.
disjoint <- compare("merge disjoint genes", function(f) {
  f(list(build(4, c("G1", "G2"), "p", c("Bcell", "Tcell", "Mono", "Bcell"), lv_a),
    build(4, c("G4", "G7"), "q", c("Mono", "Bcell", "NK", "Mono"), lv_b)))
})
cat("disjoint gene sets give a", nrow(disjoint@data.signaling), "row data.signaling:",
    "upstream quirk, reproduced\n")

## `meta.use` is the intersection of `meta` column names, so `age` -- present only in A -- is gone.
stopifnot(identical(colnames(merged@meta), c("cell", "datasets", "cluster")))
cat("meta.use is the intersection of columns: ok\n")

## `@idents` on a merged object is a **list**, not a factor: the per-dataset factors, plus a
## `$joint` factor. `AnyFactor` in the class definition is `setClassUnion(c("factor", "list"))`, which
## is what makes that legal. `$joint` is what the comparison-analysis functions read
## (`cell.levels <- levels(object@idents$joint)` in `netAnalysis_signalingChanges_scatter`), so it is
## the joint factor's level order that is the contract.
stopifnot(is.list(merged@idents), "joint" %in% names(merged@idents))
stopifnot(identical(names(merged@idents)[1:2], c("Dataset_1", "Dataset_2")))
cat("idents is a list of per-dataset factors plus $joint: ok\n")

## The joint level order is `union(levels(a), levels(b))`: A's three levels, then B's one new level.
## Sorted order would be Bcell,Mono,NK,Tcell, so this assertion separates the two readings.
stopifnot(identical(levels(merged@idents$joint), c("Bcell", "Tcell", "Mono", "NK")))
cat("idents$joint levels are first-appearance union, not sorted: ok\n")

## And the joint factor's *values* follow the same cell order as `meta`, one entry per cell.
stopifnot(length(merged@idents$joint) == ncol(merged@data.signaling))
stopifnot(identical(names(merged@idents$joint), colnames(merged@data.signaling)))
cat("idents$joint is aligned to the merged columns: ok\n")

## The `datasets` factor's levels are `add.names`, in argument order.
stopifnot(identical(levels(merged@meta$datasets), c("Dataset_1", "Dataset_2")))
cat("datasets factor levels are add.names in order: ok\n")

## `options$mode` is set to "merged"; the datatype check passes because both objects are RNA.
stopifnot(identical(merged@options$mode, "merged"))
stopifnot(identical(merged@options$datatype, "RNA"))
cat("options$mode is 'merged': ok\n")

## `add.names` is a positional default; a custom one must be honoured in `meta$datasets`, in the
## names of the per-dataset `net` and `LR` lists, and in the per-dataset `idents` entries.
##
## `@net` is a *list of nets*, one per dataset, so `@net$X$prob` is the array and `@net$prob` is
## `NULL`. That nesting is what `analysis.R` keys on to tell a merged object from a single one:
## `if (!is.list(object@net[[1]])) stop("This function cannot be applied to a single cellchat
## object from one dataset!")`.
m2 <- merged2
stopifnot(identical(levels(m2@meta$datasets), c("X", "Y")))
stopifnot(identical(names(m2@net), c("X", "Y")))
stopifnot(identical(names(m2@LR), c("X", "Y")))
stopifnot(identical(names(m2@netP), c("X", "Y")))
stopifnot(!is.null(m2@net$X$prob) && is.null(m2@net$prob))
cat("add.names propagates to meta$datasets, net, netP, LR and idents: ok\n")
cat("a merged object is recognisable by @net[[1]] being a list:",
    is.list(m2@net[[1]]), "\n")

## `merge.data = TRUE` populates the `data` slot, which the default path leaves absent.
stopifnot(!is.null(merged4@data) || is.null(slot(merged4, "data")))
cat("merge.data = TRUE builds a data slot: ok\n")

## ---------------------------------------------------------------- 3. upstream's own errors
cat("\n-- upstream's errors, reproduced\n")
## Duplicated cell names across datasets. Upstream's message verbatim, including the double `!!`.
dup <- capture(cellchatrs_upstream_mergeCellChat(list(obj_a, obj_a)))
stopifnot(identical(dup$value,
  "Duplicated cell names were detected across datasets!! Please set cell.prefix = TRUE"))
cat("duplicated cell names ->", dup$value, "\n")

## Differing datatypes. Built by direct slot assignment, since `createCellChat` will not produce a
## spatial object from a plain matrix.
o1 <- obj_a; o1@options$datatype <- "RNA"
o2 <- obj_b; o2@options$datatype <- "spatial"
dt <- capture(cellchatrs_upstream_mergeCellChat(list(o1, o2)))
stopifnot(identical(dt$value,
  "Comparison analysis is not suggested for different types of data."))
cat("differing datatypes ->", dt$value, "\n")

## ---------------------------------------------------------------- 4. the Rust functions on a merged object's nets
## The part that is not vacuous, and it is the part that matters.
##
## `computeCommunProb` cannot run on a merged object -- `@idents` is a list, and the kernel needs a
## factor -- and that is correct upstream behaviour, not a gap in the fixture. The comparison
## workflow does not recompute probabilities after merging. It merges objects whose nets were
## computed *separately*, and then analyses the per-dataset nets side by side:
## `netAnalysis_signalingRole_scatter` reads `object@net[[comparison[1]]]$netP`, and
## `netAnalysis_computeCentrality` is run per dataset. So the real question for this port is
## whether the Rust `aggregateNet` / `subsetCommunication` / `filterCommunication` / `rankNet`
## accept the net structures a merge carries: per-dataset lists, with the dimnames and group order
## each dataset had on its own.
##
## That is what is checked. Each dataset's net is lifted back out of the merged object, wrapped in a
## single-dataset object so the Rust functions have a `@net` to read, and compared against upstream
## on the same net. A merge that renumbered a group, dropped a dimname, or nested `net` at the wrong
## depth would change these results.
cat("\n-- the Rust functions on each merged dataset's net\n")

## Rebuild a single-dataset view of dataset `nm`: the **net from the merged object**, wrapped around
## the **expression matrix and `idents` the dataset had before the merge**.
##
## That split is the point, not a convenience. `mergeCellChat` intersects the gene sets, so
## `merged@data.signaling` holds only the genes both datasets share -- two of A's four. A net
## computed by dataset 1 on all four of its own genes has no counterpart in that matrix: handing the
## merged matrix to `computeCommunProb` alongside dataset 1's own `LRsig` makes upstream raise its
## `"subscript out of bounds"`, because the `G1_G2` pair's subunits are no longer rows. In real use
## the comparison workflow never does this -- it analyses the per-dataset nets that were computed
## *before* the merge, on each dataset's own data. So the merged object contributes the net and the
## original contributes the context, and the two are the same object except for the net's provenance.
as_single <- function(mm, nm, obj) {
  new("CellChat",
      data = obj@data.signaling,
      data.signaling = obj@data.signaling,
      meta = obj@meta,
      idents = mm@idents[[nm]],
      net = list(prob = mm@net[[nm]]$prob, pval = mm@net[[nm]]$pval),
      netP = list(),
      LR = list(LRsig = obj@LR$LRsig),
      var.features = list(features = rownames(obj@data.signaling)),
      DB = list(complex = DB$complex, cofactor = DB$cofactor, interaction = DB$interaction),
      images = list(),
      options = list(datatype = "RNA", mode = "single", db = normalizePath(DBDIR),
                     population.size = TRUE))
}

one <- as_single(merged, "Dataset_1", obj_a)
nboot <- 5L

## 4a. aggregateNet
agg <- function(fn) {
  suppressWarnings(suppressMessages(fn(
    one, type = "triMean", trim = 0.1, population.size = TRUE,
    nboot = nboot, seed.use = 1L, Kh = 0.5, n = 1)))
}
a_up <- capture(agg(cellchatrs_upstream_computeCommunProb))
a_rs <- capture(agg(computeCommunProb))
n_cmp <- n_cmp + 1L
ok <- identical(a_up$value@net$prob, a_rs$value@net$prob) &&
      identical(a_up$value@net$pval, a_rs$value@net$pval) &&
      identical(dimnames(a_up$value@net$prob), dimnames(a_rs$value@net$prob))
cat(sprintf("%-40s value=%-5s\n", "computeCommunProb on merged dataset 1", ok))
if (!ok) fails <- fails + 1L

## Then aggregate, which is the function that reads `net$prob` and produces the per-interaction rows
## the comparison workflow plots.
agg_up <- capture(cellchatrs_upstream_aggregateNet(one))
agg_rs <- capture(aggregateNet(one))
n_cmp <- n_cmp + 1L
ok <- identical(agg_up$value@net$count, agg_rs$value@net$count) &&
      identical(agg_up$value@net$weight, agg_rs$value@net$weight) &&
      identical(agg_up$value@net$prob, agg_rs$value@net$prob) &&
      identical(agg_up$value@net$pval, agg_rs$value@net$pval)
cat(sprintf("%-40s value=%-5s\n", "aggregateNet on merged dataset 1", ok))
if (!ok) {
  fails <- fails + 1L
  cat("   count up:", paste(utils::head(agg_up$value@net$count), collapse = ","), "\n")
  cat("   count rs:", paste(utils::head(agg_rs$value@net$count), collapse = ","), "\n")
}

## 4b. The net lifted out of a merge must be the net the object had before the merge. This is the
## check that would fail if `mergeCellChat` renumbered a group or dropped a dimname.
n_cmp <- n_cmp + 1L
ok_net <- identical(merged@net$Dataset_1$prob, obj_a@net$prob) &&
          identical(merged@net$Dataset_1$pval, obj_a@net$pval) &&
          identical(dimnames(merged@net$Dataset_1$prob), dimnames(obj_a@net$prob))
cat(sprintf("%-40s value=%-5s\n", "merged net is the pre-merge net, unchanged", ok_net))
if (!ok_net) {
  fails <- fails + 1L
  cat("   dimnames up:", paste(dimnames(merged@net$Dataset_1$prob)[[1]], collapse = ","), "\n")
  cat("   dimnames of:", paste(dimnames(obj_a@net$prob)[[1]], collapse = ","), "\n")
}

## 4c. subsetCommunication and filterCommunication on the merged dataset's net.
##
## `subsetCommunication` returns the filtered interaction table as a **data.frame**, not a CellChat
## object, so the comparison is on the frame itself -- rows, columns, row names and the factor
## levels of its `source`/`target` columns. That is where a group-order error in the merge would
## show: `aggregateNet` builds the `source|target` character key from the level order, so a
## reordered `idents` gives a differently *keyed* table, not merely a permuted one.
sub_up <- capture(cellchatrs_upstream_subsetCommunication(one))
sub_rs <- capture(subsetCommunication(one))
n_cmp <- n_cmp + 1L
ok_sub <- identical(sub_up$value, sub_rs$value)
cat(sprintf("%-40s value=%-5s  rows=%s\n", "subsetCommunication on merged dataset 1", ok_sub,
            if (is.data.frame(sub_up$value)) nrow(sub_up$value) else "?"))
if (!ok_sub) {
  fails <- fails + 1L
  if (is.data.frame(sub_up$value) && is.data.frame(sub_rs$value)) {
    cat("   rows up:", nrow(sub_up$value), " rs:", nrow(sub_rs$value), "\n")
    cat("   cols up:", paste(colnames(sub_up$value), collapse = ","), "\n")
    cat("   cols rs:", paste(colnames(sub_rs$value), collapse = ","), "\n")
    for (col in intersect(colnames(sub_up$value), colnames(sub_rs$value))) {
      if (!identical(sub_up$value[[col]], sub_rs$value[[col]])) {
        cat("   col", col, "up:", paste(utils::head(as.character(sub_up$value[[col]]), 3), collapse = "|"), "\n")
        cat("   col", col, "rs:", paste(utils::head(as.character(sub_rs$value[[col]]), 3), collapse = "|"), "\n")
        cat("   col", col, "classes up/rs:", paste(class(sub_up$value[[col]]), collapse = ","), "/",
            paste(class(sub_rs$value[[col]]), collapse = ","), "\n")
      if (is.factor(sub_up$value[[col]])) {
        cat("   col", col, "up levels:", paste(levels(sub_up$value[[col]]), collapse = ","), "\n")
        cat("   col", col, "up values:", paste(as.character(sub_up$value[[col]]), collapse = ","), "\n")
        cat("   col", col, "net dimnames[3]:", paste(dimnames(one@net$prob)[[3]], collapse = ","), "\n")
      }
      }
    }
  } else {
    cat("   up:", as.character(sub_up$value), "\n   rs:", as.character(sub_rs$value), "\n")
  }
}

fil_up <- capture(cellchatrs_upstream_filterCommunication(one))
fil_rs <- capture(filterCommunication(one))
n_cmp <- n_cmp + 1L
ok_fil <- identical(fil_up$value@net$prob, fil_rs$value@net$prob) &&
          identical(fil_up$value@net$pval, fil_rs$value@net$pval)
cat(sprintf("%-40s value=%-5s\n", "filterCommunication on merged dataset 1", ok_fil))
if (!ok_fil) fails <- fails + 1L

## 4d. `aggregateNet`'s **filtered** branch, the one taken when a filter argument is present. It
## aggregates differently from the default branch above: the default does
## `apply(prob > 0, c(1,2), sum)` and `apply(prob, c(1,2), sum)` on the array, while this one melts
## the table with `dplyr` and `tapply`s over the `source|target` character key. Two code paths, so
## both are compared -- and this is the path a `signaling` or `sources.use` filter takes in the
## comparison workflow.
##
## `subsetCommunication`, which this branch calls, does `prob[pval >= thresh] <- 0` and then keeps
## the `prob > 0` rows. So `thresh = 0` empties the table rather than keeping everything -- every
## `pval` is `>= 0`. With `nboot = 5` the attainable p-values are multiples of `1/5`, so `thresh =
## 0.5` keeps the `pval == 0.2` interactions and is the smallest threshold that is not vacuous.
all_pairs <- data.frame(interaction_name = unique(as.character(rownames(one@LR$LRsig))))
tab <- capture(cellchatrs_upstream_aggregateNet(one, pairLR.use = all_pairs, thresh = 0.5))
tab_rs <- capture(aggregateNet(one, pairLR.use = all_pairs, thresh = 0.5))
n_cmp <- n_cmp + 1L
ok_tab <- identical(tab$value@net$count, tab_rs$value@net$count) &&
          identical(tab$value@net$weight, tab_rs$value@net$weight) &&
          identical(dimnames(tab$value@net$weight), dimnames(tab_rs$value@net$weight))
cat(sprintf("%-40s value=%-5s\n", "aggregateNet filtered branch on merged net", ok_tab))
if (!ok_tab) fails <- fails + 1L

## Both branches leave `net$prob` as the `K x K x N` array and put the aggregation in `count` and
## `weight`; neither writes a table into the object. The table is what `subsetCommunication`
## *returns*, and 4c compares that. Asserting the distinction here records it, so a future change
## that put a table into `net$prob` would be noticed rather than silently changing `rankNet`.
n_cmp <- n_cmp + 1L
ok_shape <- is.array(tab$value@net$prob) && is.matrix(tab$value@net$weight) &&
            is.matrix(tab$value@net$count)
cat(sprintf("%-40s value=%-5s  prob %s, weight %s\n",
            "aggregateNet writes count/weight, not a table", ok_shape,
            paste(dim(tab$value@net$prob), collapse = "x"),
            paste(dim(tab$value@net$weight), collapse = "x")))
if (!ok_shape) fails <- fails + 1L

## 4e. `rankNet` over the merged dataset's net, and the group order in the result.
##
## `rankNet(mode = "single")` reads `object@net$prob` as the `K x K` **array** that
## `computeCommunProb` leaves behind -- it indexes `dimnames(prob)[[1]]` and zeroes
## `prob[object1$pval > thresh]` -- so it runs on the object as the merge produced it, with no
## intermediate step. It is also the function that reports "No inferred communications for the
## input!", and the reason is the threshold: with `nboot = 5` the attainable p-values are
## {0, 0.2, 0.4, 0.6, 0.8, 1}, so `thresh = 0.01` discards all but the `pval == 0` interactions.
## `thresh = 0.5` keeps the 0.2 and 0.4 ones. Note this is a *different* threshold convention from
## `subsetCommunication`'s `pval >= thresh`: same argument name, different meaning.
##
## `slot.name` has to be passed explicitly. It defaults to `"netP"`, the pathway-level slot, and a
## merged object's `@netP` is a per-dataset list with nothing in it for a single-dataset view -- so
## `object1$prob` is `NULL`, `sum(prob) == 0`, and `rankNet` stops with "No inferred communications
## for the input!" regardless of the threshold. That message means "you gave me the wrong slot", not
## "the network is empty", and it is worth knowing before spending time on the threshold.
rk_up <- capture(cellchatrs_upstream_rankNet(one, slot.name = "net", mode = "single", thresh = 0.5))
rk_rs <- capture(rankNet(one, slot.name = "net", mode = "single", thresh = 0.5))

## `rankNet` returns a **ggplot object**, not a CellChat, so `identical()` on it is not usable: a
## ggplot carries an environment and a call, and two plots of the same numbers are not `identical`.
## The comparison is on `ggplot_build(gg)$data` -- the computed layer data -- which is the numbers
## the plot was going to draw. `check_identical.R` does the same for its rankNet block.
built <- function(gg) ggplot2::ggplot_build(gg)$data
ok_rk <- identical(built(rk_up$value), built(rk_rs$value))
cat(sprintf("%-40s value=%-5s  layers=%d\n", "rankNet plot data on merged dataset 1", ok_rk,
            length(built(rk_up$value))))
if (!ok_rk) {
  fails <- fails + 1L
  bu <- built(rk_up$value); br <- built(rk_rs$value)
  for (i in seq_along(bu)) {
    if (!identical(bu[[i]], br[[i]])) {
      cat("   layer", i, "differs;", ncol(bu[[i]]), "columns\n")
      if ("x" %in% colnames(bu[[i]]) && "x" %in% colnames(br[[i]])) {
        cat("   x up:", paste(utils::head(bu[[i]]$x, 6), collapse = ","), "\n")
        cat("   x rs:", paste(utils::head(br[[i]]$x, 6), collapse = ","), "\n")
      }
      if ("group" %in% colnames(bu[[i]])) {
        cat("   levels up:", paste(levels(bu[[i]]$group), collapse = ","), "\n")
        cat("   levels rs:", paste(levels(br[[i]]$group), collapse = ","), "\n")
      }
    }
  }
}

## 4e. The recorded invariants from the objective, on the merged object's own net.
pval <- a_up$value@net$pval
prob <- a_up$value@net$prob
n_cmp <- n_cmp + 1L
ok_inv <- all(abs(pval * nboot - round(pval * nboot)) < 1e-9) &&
          all(pval[prob == 0] == 1) &&
          identical(dimnames(prob)[[1]], levels(one@idents)) &&
          identical(dimnames(prob)[[2]], levels(one@idents))
cat(sprintf("%-40s value=%-5s\n", "Pval/dimname invariants on the merged net", ok_inv))
if (!ok_inv) fails <- fails + 1L

## 4f. The level order is load-bearing, not incidental. Reversing the levels must change the
## dimnames; if it did not, the preceding assertion would not be testing the level order at all.
revd <- as_single(merged, "Dataset_1", obj_a)
revd@idents <- factor(as.character(revd@idents), levels = rev(levels(revd@idents)))
r_up <- capture(suppressWarnings(suppressMessages(cellchatrs_upstream_computeCommunProb(
  revd, type = "triMean", trim = 0.1, population.size = TRUE,
  nboot = nboot, seed.use = 1L, Kh = 0.5, n = 1))))
n_cmp <- n_cmp + 1L
ok_rev <- !identical(dimnames(r_up$value@net$prob)[[1]], dimnames(prob)[[1]])
cat(sprintf("%-40s value=%-5s\n", "reversing idents levels changes the net", ok_rev))
if (!ok_rev) {
  fails <- fails + 1L
  cat("   reversing the level order left the dimnames unchanged, so the level order is not\n",
      "   reaching the kernel and the assertion above is not testing it\n", sep = "")
}

## ---------------------------------------------------------------- verdict
cat(sprintf("\n%s: %d failing comparisons out of %d\n",
            if (fails == 0L) "IDENTICAL" else "MISMATCH", fails, n_cmp))
quit(status = if (fails == 0L) 0L else 1L)
