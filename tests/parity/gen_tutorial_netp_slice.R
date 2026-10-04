# Dump real pathway slices for the centrality corpus (one-off fixture generation).
suppressWarnings(suppressMessages({
  library(methods); library(Matrix); library(collapse); library(dplyr)
}))
suppressWarnings(suppressMessages(library(CellChat)))
CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")
SKIN <- Sys.getenv("TUTORIAL_DATA", "data/humanSkin.rda")
q <- function(x) { invisible(utils::capture.output(value <- suppressWarnings(suppressMessages(x)))); value }
UP <- new.env(parent = globalenv())
for (p in c("Matrix", "collapse", "dplyr", "future", "rlang")) suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R", "visualization.R")) suppressWarnings(suppressMessages(sys.source(file.path(CC, "R", f), envir = UP, keep.source = FALSE)))
suppressWarnings(suppressMessages(sys.source(file.path(CC, "R", "CellChat_class.R"), envir = globalenv(), keep.source = FALSE)))
E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E); DB <- get(ls(E)[1], E)
raw <- new.env(); load(SKIN, envir = raw); skin <- raw[[ls(raw)[1]]]
counts <- skin$data
grp <- if (!is.null(skin$meta)) factor(skin$meta$labels) else factor(skin$labels)
lib <- Matrix::colSums(counts); cpm <- counts
cpm@x <- cpm@x / rep(lib, diff(cpm@p)) * 1e4
norm <- log1p(as.matrix(cpm))
meta <- data.frame(labels = grp, row.names = colnames(counts))
obj <- createCellChat(object = norm, meta = meta, group.by = "labels")
obj@DB <- DB
obj <- q(UP$subsetData(obj))
obj <- q(CellChat::identifyOverExpressedGenes(obj))
obj <- q(CellChat::identifyOverExpressedInteractions(obj))
obj <- q(CellChat::computeCommunProb(obj, type = "triMean"))
obj <- q(CellChat::filterCommunication(obj, min.cells = 10))
obj <- q(CellChat::computeCommunProbPathway(obj))
np <- obj@netP
con <- file("tests/fixtures/tutorial_netp_slice.tsv", "wt")
k <- dim(np$prob)[1]; npath <- min(6L, dim(np$prob)[3])
cat(sprintf("k\t%d\nnpath\t%d\n", k, npath), file = con)
flat <- vapply(seq_len(npath), function(p) paste(sprintf("%.17g", as.numeric(np$prob[, , p])), collapse = " "), "")
writeLines(paste(flat, collapse = " "), con)
close(con)
cat("dumped", npath, "pathway slices, k =", k, "\n")
