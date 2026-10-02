## The CellChat input pipeline, factored out so the R baseline and the Rust measurement are
## provably the same computation.
##
## This is `bench-runner/prof_real.R`'s preamble, moved verbatim. That file is the source of the
## recorded R baselines (89.7 s human skin, 206.5 s mouse wound) and must not be edited in a way
## that would invalidate them, so the shared code lives here and `prof_real.R` keeps its own copy
## for exactly that reason. `bench_real.R` asserts that the object built here has the same
## dimensions, nnz and `nLR` as `prof_real.R` reports, so a drift between the two is a *failing
## check* rather than a silently incomparable speedup.
##
## Emulates what `updateCellChat()` leaves in a CellChat object: log1p CPM over the signalling
## genes, `data.signaling` restricted to genes that any L-R pair needs (directly or as complex
## subunits), and `LRsig` filtered to pairs whose every subunit is present.
suppressMessages({library(collapse); library(Matrix); library(dplyr); library(cellchatrs)})
CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
DB <- local({
  e <- new.env()
  load(file.path(CC, "data", paste0("CellChatDB.", Sys.getenv("SPECIES", "human"), ".rda")), envir = e)
  get(ls(e)[1L], e)
})
upstream_env <- function() {
  e <- new.env()
  source(file.path(CC, "R", "modeling.R"), local = e)
  e
}

## A probe class rather than `CellChat` itself: `computeCommunProb` touches `@data.signaling`,
## `@data.smooth`, `@idents`, `@LR$LRsig`, `@options$datatype` and `@net`, and instantiating a
## real CellChat would drag in the databases and the plotting surface for no gain. `prof_real.R`
## defines the identical class, and the two must stay identical or the comparison is void.
setClass("CellChatProbe", representation(
  data = "ANY", data.signaling = "ANY", data.smooth = "ANY", idents = "ANY", meta = "ANY",
  images = "ANY", DB = "ANY", LR = "ANY", net = "ANY", netP = "ANY", options = "ANY"
))

extractSig <- function(DB) {
  ii <- DB$interaction
  cs <- as.matrix(DB$complex[, grepl("^subunit", colnames(DB$complex)), drop = FALSE])
  cof <- unique(c(ii$agonist, ii$antagonist, ii$co_A_receptor, ii$co_I_receptor))
  cof <- cof[!is.na(cof) & cof != ""]
  cfm <- DB$cofactor[match(cof, rownames(DB$cofactor), nomatch = 0L), , drop = FALSE]
  cm <- as.matrix(cfm[, grepl("^cofactor", colnames(cfm)), drop = FALSE])
  unique(c(ii$ligand, ii$receptor,
           unique(as.character(cs)[cs != ""]), unique(as.character(cm)[cm != ""])))
}

## `build_fixture(path, max_lr = 0L)` -> list(object = <CellChatProbe>, info = list(...)).
build_fixture <- function(path, max_lr = 0L, seed = 1L) {
  raw <- new.env()
  load(path, envir = raw)
  o <- raw[[ls(raw)[1L]]]
  counts <- o$data
  grp <- if (!is.null(o$meta)) factor(o$meta$labels) else factor(o$labels)
  nC <- ncol(counts)
  genes.sig <- sort(intersect(extractSig(DB), rownames(counts)))
  X <- as.matrix(counts[genes.sig, , drop = FALSE])
  lib <- colSums(counts)
  X <- X / rep(lib, each = nrow(X)) * 1e4
  X <- log1p(X)
  subOK <- function(g) {
    if (g %in% rownames(counts)) return(TRUE)
    if (!(g %in% rownames(DB$complex))) return(FALSE)
    s <- as.character(DB$complex[g, ]); s <- s[!is.na(s) & s != ""]
    length(s) > 0L && all(s %in% rownames(counts))
  }
  LR <- DB$interaction[DB$interaction$annotation %in%
    c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact"), ]
  keep <- vapply(LR$ligand, subOK, TRUE) & vapply(LR$receptor, subOK, TRUE)
  LR <- LR[keep, , drop = FALSE]
  if (max_lr > 0L) LR <- LR[seq_len(min(max_lr, nrow(LR))), , drop = FALSE]
  M <- as(Matrix(X, sparse = FALSE), "dgCMatrix")
  obj <- methods::new("CellChatProbe", data = M, data.signaling = M, data.smooth = NULL,
    idents = grp, meta = data.frame(samples = factor(rep("s1", nC))), images = list(),
    DB = DB, LR = list(LRsig = LR), net = list(), netP = list(),
    options = list(datatype = "RNA", mode = "single"))
  list(
    object = obj,
    info = list(
      fixture = basename(path), n_genes_total = nrow(counts), n_cells = nC,
      n_groups = nlevels(grp), n_genes_signaling = length(genes.sig), n_lr = nrow(LR),
      nnz_pct = 100 * Matrix::nnzero(as.matrix(counts[genes.sig, ])) / length(counts[genes.sig, ])
    )
  )
}

## Peak resident set of *this* R process, in MB. `gc()`'s own accounting is R's heap only and
## understates the cost of a dgCMatrix-to-dense conversion by more than an order of magnitude on
## the wound fixture, so it is not what a user would see in `ps`.
peak_rss_mb <- function() {
  st <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  hit <- grep("^VmHWM:", st, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", hit)) / 1024
}

## Available memory, in MB. The objective forbids quoting a result at or above 50 000 cells
## without an explicit headroom check, and the check has to read the *kernel's* number:
## MemFree excludes reclaimable page cache and reads near zero on a healthy machine.
mem_available_mb <- function() {
  st <- tryCatch(readLines("/proc/meminfo"), error = function(e) character(0))
  hit <- grep("^MemAvailable:", st, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", hit)) / 1024
}

## The dense buffer the kernel streams once per bootstrap replicate is the thing that actually
## decides whether a run fits: `n_genes_signaling * n_cells * 8` bytes, plus the sparse input.
required_mb <- function(info) {
  (info$n_genes_signaling * info$n_cells * 8) / 1024^2
}

check_headroom <- function(info, factor = 1.6) {
  need <- required_mb(info) * factor
  have <- mem_available_mb()
  ok <- is.na(have) || have >= need
  list(required_mb = need, available_mb = have, ok = ok,
      message = if (ok) "ok" else
        sprintf("NOT QUOTABLE: needs ~%.0f MB (x%.1f) but only %.0f MB is available",
                need, factor, have))
}
