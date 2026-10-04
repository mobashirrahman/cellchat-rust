## Decisive check that `computeCommunProb` reaches the Rust kernel and not pinned upstream.
##
## Trace the actual native binding. A custom object DB works even without disk exports,
## so successful output alone cannot prove the Rust kernel was called.
suppressWarnings(suppressMessages({
  library(Matrix); library(collapse); library(dplyr); library(CellChat)
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

## Trace the registered Rust entry point as well as requiring a kernel error. The deliberately
## unknown interaction reaches the extension with a valid export directory, then fails in kernel
## database resolution. This makes the Rust and upstream paths observably different.
stopifnot(!tolower(Sys.getenv("CELLCHATRS_FALLBACK", "0")) %in% c("1", "true", "yes", "on"))
ns <- asNamespace("CellChat")
stopifnot(is.function(get("compute_commun_prob", envir = ns, inherits = FALSE)))
.cellchatrs_native_probe_count <- 0L
trace("compute_commun_prob", where = ns,
      tracer = quote(assign(".cellchatrs_native_probe_count",
                           .GlobalEnv$.cellchatrs_native_probe_count + 1L,
                           envir = .GlobalEnv)), print = FALSE)

Sys.unsetenv("CELLCHATRS_DB")
o_bad <- mk(tempfile("nonexistent-db-"))
r2 <- try(quiet(o_bad), silent = TRUE)
stopifnot(!inherits(r2, "try-error"))
stopifnot(identical(.cellchatrs_native_probe_count, 1L))
stopifnot(is.list(r2@net), identical(dim(r2@net$prob), c(2L, 2L, 1L)))
untrace("compute_commun_prob", where = ns)
rm(.cellchatrs_native_probe_count, envir = .GlobalEnv)
cat("Rust path: reached the native kernel and returned a probability array\n")
