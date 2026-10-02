## Decisive check that `computeCommunProb` reaches the Rust kernel and not pinned upstream.
##
## A green differential gate cannot distinguish "the port matches" from "the port was never
## called" -- that is exactly how the `nzchar(Sys.getenv(..., "0"))` bug survived a fully green
## 273/273 run. The two paths are made distinguishable from the outside: the Rust path resolves
## its database through `cellchatrs_db_dir`, which *requires a readable directory on disk*, while
## the upstream path uses `object@DB` and needs no directory at all. Point `options$db` at a
## directory that does not exist and the Rust path must fail with the db error; if it succeeds, the
## call went upstream.
suppressWarnings(suppressMessages({
  library(Matrix); library(collapse); library(dplyr); library(cellchatrs)
}))
setClass("PChat", representation(
  data.signaling = "matrix", data.smooth = "matrix", idents = "factor", meta = "data.frame",
  images = "list", LR = "list", net = "list", options = "list", DB = "list"))

ncell <- 12L; ngene <- 4L
expr <- matrix(abs(seq_len(ncell * ngene) / 37) + 0.1, nrow = ngene, ncol = ncell)
rownames(expr) <- paste0("g", seq_len(ngene)); colnames(expr) <- paste0("c", seq_len(ncell))
levs <- c("A", "B")
group <- factor(rep(levs, each = ncell / 2), levels = levs)
LRsig <- data.frame(ligand = "g1", receptor = "g2", annotation = "Secreted Signaling",
                    interaction_name = "LR1", pathway_name = "pA", stringsAsFactors = FALSE)
mk <- function(dbval) {
  new("PChat", data.signaling = expr, data.smooth = expr, idents = group,
      meta = data.frame(samples = factor(rep("s1", ncell), levels = "s1"),
                        row.names = colnames(expr)),
      images = list(), LR = list(LRsig = LRsig), net = list(),
      options = list(mode = "single", datatype = "RNA", db = dbval),
      DB = list(interaction = data.frame(interaction_name_2 = "LR1", ligand = "g1"),
                complex = data.frame(interaction_name_2 = character(0), cofactor = character(0)),
                cofactor = data.frame(interaction_name_2 = character(0),
                                      co_A_receptor = character(0), co_I_receptor = character(0)),
                geneInfo = data.frame(symbol = paste0("g", seq_len(ngene)))))
}
## `capture.output` returns the *text*, so the assignment has to happen in this frame and the
## object returned explicitly. `capture.output(...)` alone yields a character vector, and
## `r@net$prob` then fails on a `character` -- which looks like a kernel bug and is not.
quiet <- function(o, ...) {
  suppressWarnings(suppressMessages({
    invisible(capture.output(res <- computeCommunProb(o, nboot = 2, ...)))
    res
  }))
}

## 1. The db path exists -> the Rust path runs to completion.
o <- mk(getOption("cellchatrs.probe.db", Sys.getenv("CELLCHATRS_DB", "")))
if (nzchar(o@options$db) && dir.exists(o@options$db)) {
  r <- try(quiet(o), silent = TRUE)
  cat("with a real db dir:",
      if (inherits(r, "try-error")) paste("ERROR:", conditionMessage(attr(r, "condition")))
      else sprintf("ok, %d positive of %d", sum(r@net$prob > 0), length(r@net$prob)), "\n")
} else {
  cat("with a real db dir: SKIPPED (set CELLCHATRS_DB)\n")
}

## 2. The db path does not exist. The Rust path must refuse; upstream never looks at a directory.
o_bad <- mk("/nonexistent/cellchatrs-db-probe")
r2 <- try(quiet(o_bad), silent = TRUE)
if (inherits(r2, "try-error")) {
  cat("with a bogus db dir: ERROR ->", conditionMessage(attr(r2, "condition")), "\n")
  cat("  => reached the Rust kernel:", grepl("CellChatDB export", conditionMessage(attr(r2, "condition"))), "\n")
} else {
  cat("with a bogus db dir: SUCCEEDED => the call did NOT go through cellchatrs_db_dir\n")
}

## 3. The fallback predicate itself, which is the bug that hid all of this.
cat("fallback predicate is FALSE by default:",
    !(tolower(Sys.getenv("CELLCHATRS_FALLBACK", "0")) %in% c("1", "true", "yes", "on")), "\n")
cat("fallback predicate is TRUE when set to 1:",
    tolower(Sys.getenv("CELLCHATRS_FALLBACK", "1")) %in% c("1", "true", "yes", "on"), "\n")
