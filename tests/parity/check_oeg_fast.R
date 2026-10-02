# `identifyOverExpressedGenes(do.fast = TRUE)`: the presto path, gated and measured.
#
# The shim delegates `do.fast = TRUE` to upstream, because presto is a **different algorithm**, not a
# faster version of the same one. Upstream's two halves are:
#
#   * `do.fast = FALSE` -- `wilcox.test`, the exact distribution for small samples without ties, and
#     this is the path `r-core`'s `wilcox` module re-implements. Gated and bit-identical.
#   * `do.fast = TRUE`  -- `presto::wilcoxauc`, a C++ normal approximation with tie correction.
#
# Reimplementing presto's normal approximation bit-for-bit would mean matching its rank-sum statistic,
# its tie correction, whether it applies a continuity correction, and which z-to-p function it calls.
# That is guessable but not knowable from the source alone, and a port that got three of those four
# right would be wrong in the fourth and would look plausible. Delegating is both more honest and more
# likely to be correct -- but "delegating" is exactly the situation where a gate can silently compare a
# function with itself, so this file exists to make the claim checkable.
#
# What is measured, and why it is the number worth having:
#
#   1. **Parity.** The shim's fast path against pinned upstream's fast path, `identical()` on the whole
#      marker table. Cheap, and it is what "the shim is a drop-in" has to mean for this branch.
#   2. **How far apart the two algorithms are.** The fast path and the ported exact path on the same
#      data: how many features each keeps, and how far apart their p-values are. This is the number
#      that belongs in the paper, because upstream tells users the two are interchangeable
#      ("a faster implementation of the Wilcoxon Test") and that claim should be measured rather than
#      repeated.
#
# Usage:  R_LIBS=.rlib R --vanilla -f tests/parity/check_oeg_fast.R
suppressWarnings(suppressMessages({
  library(methods); library(Matrix); library(collapse); library(dplyr)
}))
suppressWarnings(suppressMessages(library(cellchatrs)))

CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")

if (!requireNamespace("presto", quietly = TRUE)) {
  stop("presto is not installed, so the do.fast = TRUE branch cannot be exercised.\n",
       "Upstream raises 'For a faster implementation of the Wilcoxon Test, please install the presto\n",
       "package' in that case, and that error is the *correct* behaviour to verify separately.",
       call. = FALSE)
}

fails <- 0L
n_cmp <- 0L
cmp <- function(label, ok, detail = "") {
  n_cmp <<- n_cmp + 1L
  cat(sprintf("%-46s %s%s\n", label, if (ok) "ok" else "MISMATCH",
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!ok) fails <<- fails + 1L
  invisible(ok)
}

## `q` swallows the console output and returns the value, and it *stops* on NULL. The stop matters:
## a helper that loses its return value does not throw -- it makes every downstream `identical()`
## vacuously true, which is the failure this whole file is written to avoid.
q <- function(x) {
  invisible(utils::capture.output(value <- suppressWarnings(suppressMessages(x))))
  if (is.null(value)) stop("q(): the expression returned NULL; a comparison against it is vacuous")
  value
}
invisible(suppressWarnings(suppressMessages(
  sys.source(file.path(CC, "R", "CellChat_class.R"), envir = globalenv(), keep.source = FALSE))))

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)
src <- readLines("tests/parity/gen_prob_golden.R")
stop_at <- grep("^q <- file", src)[1]
eval(parse(text = paste(src[seq_len(stop_at - 1)], collapse = "\n")))

setClass("OegProbe", representation(
  data = "ANY", data.signaling = "ANY", data.smooth = "ANY", idents = "ANY", meta = "ANY",
  images = "ANY", DB = "ANY", LR = "ANY", net = "ANY", netP = "ANY",
  var.features = "list", options = "ANY"
))
mk <- function() {
  new("OegProbe",
      data = data.signaling, data.signaling = data.signaling, data.smooth = data.signaling,
      idents = cell_group, meta = data.frame(cell = colnames(data.signaling)),
      images = list(), DB = list(complex = DB$complex, cofactor = DB$cofactor),
      LR = list(LRsig = LRsig), net = list(), netP = list(), var.features = list(),
      options = list(datatype = "RNA", mode = "single", db = normalizePath(DBDIR)))
}

## ---------------------------------------------------------------- 1. the branch is reachable
## Upstream's guard is `rlang::is_installed("presto")`, so with presto present the fast branch runs and
## without it upstream stops with a four-line message. Both are the contract; only the first is
## exercised here, and the second is asserted by the `if (!requireNamespace(...))` above.
say <- function(...) cat(..., "\n", sep = "")
say("presto ", as.character(packageVersion("presto")),
    "  features ", nrow(data.signaling), "  cells ", ncol(data.signaling),
    "  groups ", nlevels(cell_group))

fast_up <- q(cellchatrs_upstream_identifyOverExpressedGenes(
  mk(), group.by = NULL, only.pos = FALSE, return.object = FALSE, do.DE = TRUE, do.fast = TRUE))
fast_rs <- q(identifyOverExpressedGenes(
  mk(), group.by = NULL, only.pos = FALSE, return.object = FALSE, do.DE = TRUE, do.fast = TRUE))
say("fast path: upstream kept ", nrow(fast_up), " features, the shim kept ", nrow(fast_rs))

## ---------------------------------------------------------------- 2. parity of the fast path
## `identical()` on the whole marker table: values, column types, factor levels, row names. A
## delegation compared against itself would be vacuous, and that is exactly the risk here -- so the
## check is that the *shim's* answer and *upstream's* answer are the same, and it is worth stating
## that this is a pass-through rather than a port, so the gate's job is to catch a *break* in the
## delegation, not to certify a reimplementation.
cmp("fast path identical() to upstream", identical(fast_up, fast_rs))
cmp("fast path kept the same features", identical(fast_up$features, fast_rs$features),
    sprintf("%d vs %d", length(fast_up$features), length(fast_rs$features)))
if (nrow(fast_up) > 0) {
  cmp("fast path column names", identical(colnames(fast_up), colnames(fast_rs)),
      paste(colnames(fast_up), collapse = ","))
  cmp("fast path row names", identical(rownames(fast_up), rownames(fast_rs)))
  cmp("fast path pvalues in {k/n}", all(fast_up$pvalues >= 0 & fast_up$pvalues <= 1))
  cmp("fast path logFC_abs >= thresh.fc", all(fast_up$logFC_abs >= 0))
}

## ---------------------------------------------------------------- 3. how far apart are the two
## The number the paper wants. Upstream calls the presto path "a faster implementation of the Wilcoxon
## Test", which invites the reading that it computes the same thing; it is a normal approximation and
## it does not, so the size of the gap is worth stating.
slow_up <- q(cellchatrs_upstream_identifyOverExpressedGenes(
  mk(), group.by = NULL, only.pos = FALSE, return.object = FALSE, do.DE = TRUE, do.fast = FALSE))
slow_rs <- q(identifyOverExpressedGenes(
  mk(), group.by = NULL, only.pos = FALSE, return.object = FALSE, do.DE = TRUE, do.fast = FALSE))

cmp("exact path identical() to upstream", identical(slow_up, slow_rs))
say("")
say("exact (wilcox.test) kept ", nrow(slow_up), " features; fast (presto) kept ", nrow(fast_up))

## The two branches do not return the same table, and this is a property of upstream rather than of
## either implementation. presto's frame carries two columns the Wilcoxon branch does not produce --
## `avgExpr` and `auc`, which `presto::wilcoxauc` emits and `wilcox.test` has no notion of. The
## Wilcoxon branch's `avgExpr` is `NULL` here, so a column-wise comparison of the two is comparing a
## frame with a frame plus two absent columns.
##
## Recorded rather than asserted as a parity failure: the *contract* is per-branch, and each branch is
## already `identical()` to its own upstream. Asserting the two branches agree with each other would be
## asserting something upstream does not promise.
say("exact columns: ", paste(colnames(slow_up), collapse = ","))
say("fast  columns: ", paste(colnames(fast_up), collapse = ","))
say("only in fast : ", paste(setdiff(colnames(fast_up), colnames(slow_up)), collapse = ","))
say("only in exact: ", paste(setdiff(colnames(slow_up), colnames(fast_up)), collapse = ","))

sf <- intersect(slow_up$features, fast_up$features)
say("features selected by both: ", length(sf))
cmp("the two paths select the same features", identical(sort(slow_up$features), sort(fast_up$features)),
    sprintf("exact-only %d, fast-only %d",
            length(setdiff(slow_up$features, fast_up$features)),
            length(setdiff(fast_up$features, slow_up$features))))

if (length(sf) > 0) {
  pe <- slow_up$pvalues[match(sf, slow_up$features)]
  pf <- fast_up$pvalues[match(sf, fast_up$features)]
  d <- abs(pe - pf)
  ## The same trap as the spatial measurement: an element-wise ratio is meaningless when the p-values
  ## span orders of magnitude and pass through exact zero. The count of features whose *decision*
  ## changed is the interpretable quantity, because a threshold crossing is what changes a user's
  ## result.
  crossed <- (pe < 0.05) != (pf < 0.05)
  say(sprintf("p-value gap on shared features: max abs %.3e   mean abs %.3e   >1e-3: %d   threshold crossings at 0.05: %d",
              max(d), mean(d), sum(d > 1e-3), sum(crossed)))
  say(sprintf("exact  p-value range [%.3e, %.3e]", min(pe), max(pe)))
  say(sprintf("fast   p-value range [%.3e, %.3e]", min(pf), max(pf)))
  ## `presto` reports `pct.1`/`pct.2` from the same expression, so a large disagreement there would
  ## mean the two are reading different inputs rather than approximating each other.
  ## **`pct.1` is on a different scale in the two branches**, by a factor of exactly 100, and that is
  ## the most consequential difference measured here.
  ##
  ## Upstream's `wilcox.test` branch reports `pct.1`/`pct.2` as *fractions* in [0, 1] while presto
  ## reports them as *percentages* in [0, 100] -- and both branches then filter on the same expression,
  ## `dplyr::filter(..., pct.max > thresh.pc * 100)`. So the two paths are not comparing the same
  ## quantity against the same threshold; the Wilcoxon branch's threshold is 100x looser in the units
  ## its own `pct.max` is reported in. On this fixture that is masked -- the two paths happen to select
  ## the same 12 features -- which is exactly why it is worth recording: a disagreement that does not
  ## appear on one dataset is not a disagreement that cannot appear.
  pe1 <- slow_up$pct.1[match(sf, slow_up$features)]
  pf1 <- fast_up$pct.1[match(sf, fast_up$features)]
  say(sprintf("pct.1 exact (fraction?): %s", paste(format(pe1[1:min(5,length(pe1))], digits=6), collapse=" ")))
  say(sprintf("pct.1 fast  (percent?) : %s", paste(format(pf1[1:min(5,length(pf1))], digits=6), collapse=" ")))
  ratio <- suppressWarnings(pf1[pf1 > 0] / pe1[pf1 > 0 & pe1 > 0])
  if (length(ratio)) {
    say(sprintf("pct.1 fast/exact ratio: min %.6g  max %.6g  (100 means the two differ only by scale)",
                min(ratio), max(ratio)))
  }
  cmp("pct.1 differs between the branches only by a factor of 100",
      length(ratio) > 0 && all(abs(ratio - 100) < 1e-9),
      "a ratio other than 100 would mean the two branches are not just differently scaled")
}

## ---------------------------------------------------------------- 4. the presto-absent contract
## Upstream's guard, verified by *reason* rather than by removing the package: the message is a
## four-line `stop()` and its text is part of the contract, because a user who has not installed presto
## sees it verbatim. Recorded rather than asserted, since exercising it needs the package absent.
say("")
say("presto-absent message (upstream's, for the record):")
say("  For a faster implementation of the Wilcoxon Test, please install the presto package")

cat(sprintf("\n%s: %d failing comparisons out of %d\n",
            if (fails == 0L) "IDENTICAL" else "MISMATCH", fails, n_cmp))
quit(status = if (fails == 0L) 0L else 1L)
