# Golden corpus for r-core's prob.rs: runs upstream `computeCommunProb` (sourced verbatim
# from the pinned commit, not reimplemented) over a small fixture matrix and records
# Prob/Pval for several parameter configurations.
suppressWarnings(suppressMessages({library(Matrix); library(collapse); library(future); library(dplyr)}))

CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
env <- new.env()
source(file.path(CC, "R", "modeling.R"), local = env)

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)

## ------------------------------------------------------------------ minimal S4 object
## `computeCommunProb` only ever touches `@data.signaling`, `@LR$LRsig`, `@DB$complex`,
## `@DB$cofactor`, `@idents`, `@options$datatype` and (for non-RNA) `@images`/`@meta`.
## A minimal class is enough and keeps the oracle independent of the rest of the package,
## which needs deps (FNN, BiocNeighbors, igraph) that the numeric surface does not.
setClass("MiniCellChat", representation(
  data.signaling = "Matrix",
  LR             = "list",
  LRsig          = "data.frame",
  DB             = "ANY",
  idents         = "factor",
  options        = "list",
  net            = "ANY"
))

## ------------------------------------------------------------------ the fixture data
set.seed(20240101)
NC <- 60L
NG <- 4L
## Group labels deliberately *not* in first-appearance order: the first cell is a member of
## "g3" so that anything indexing groups by appearance rather than by factor level is
## immediately wrong. Upstream errors on this only if a level is empty, which none is.
cell_group <- factor(
  c(rep("g1", 20), rep("g2", 15), rep("g3", 15), rep("g4", 10))[sample.int(NC)],
  levels = c("g1", "g2", "g3", "g4")
)
stopifnot(nlevels(cell_group) == length(unique(cell_group)))

## Genes: a few plain ones plus the subunits of every complex/cofactor the L-R table names.
LR <- data.frame(
  ligand       = c("G1", "Activin AB", "IL12AB", "G3", "G4", "INHBA", "IL23A", "G5"),
  receptor     = c("G2", "IL12B", "G6", "Activin AB", "G7", "FST", "G8", "IL12A"),
  agonist      = c("", "", "TGFb agonist", "", "", "TGFb antagonist", "", ""),
  antagonist   = c("", "ACTIVIN antagonist", "", "NODAL agonist", "", "", "", ""),
  co_A_receptor = c("", "", "", "TGFb agonist", "", "", "", ""),
  co_I_receptor = c("", "", "TGFb inhibition receptor", "", "", "", "", ""),
  annotation   = c("Secreted Signaling", "Secreted Signaling", "Secreted Signaling",
                    "ECM-Receptor", "ECM-Receptor", "Non-protein Signaling",
                    "Secreted Signaling", "Secreted Signaling"),
  stringsAsFactors = FALSE
)
LRsig <- LR
LRsig$ligand_antagonist <- ifelse(LR$antagonist == "", LR$ligand,
                                  paste0(LR$ligand, " - ", LR$antagonist))
LRsig$interaction_name <- paste(LR$ligand, LR$receptor, sep = "_")
rownames(LRsig) <- LRsig$interaction_name

complex_subunits_of <- function(n) {
  if (n %in% rownames(DB$complex)) {
    as.character(DB$complex[n, grepl("^subunit", colnames(DB$complex))])
  } else {
    NA_character_
  }
}
cofactor_subunits_of <- function(n) {
  if (n %in% rownames(DB$cofactor)) {
    as.character(DB$cofactor[n, grepl("cofactor", colnames(DB$cofactor))])
  } else {
    NA_character_
  }
}
subunits <- unique(c(
  unlist(lapply(c(LR$ligand, LR$receptor), complex_subunits_of)),
  unlist(lapply(c(LR$agonist, LR$antagonist, LR$co_A_receptor, LR$co_I_receptor),
                cofactor_subunits_of))
))
subunits <- subunits[!is.na(subunits) & subunits != ""]
## Plain (non-complex) ligand/receptor names are looked up in the *matrix*, so they have
## to be in the gene universe even when they are nobody's subunit. `IL23A` is exactly
## that case, and omitting it made the whole fixture die with upstream's
## `subscript out of bounds` (R-ism 11) and no hint that the fault was the fixture.
plain_lr <- setdiff(c(LR$ligand, LR$receptor), rownames(DB$complex))
genes <- unique(c(paste0("G", 1:8), plain_lr, subunits))
NGENES <- length(genes)

## `object@data.signaling` is **genes x cells** upstream: `nC <- ncol(data.use)` counts
## cells, and `aggregate(t(data.use), list(group), ...)` transposes to cells x genes first.
## Getting this the other way round fails with `arguments must have same length` from
## `aggregate.data.frame`, which is a much clearer error than it deserves to be.
expr <- matrix(runif(NGENES * NC, 0.01, 1), nrow = NGENES, ncol = NC,
               dimnames = list(genes, paste0("c", 1:NC)))
## Engineered structure so the branches are actually reached rather than drowned in noise:
##   - an all-zero gene (L1 zero short-circuit in the trimean)
##   - a gene expressed in a single group only
##   - a gene constant across cells
##   - a zero ligand/receptor pair, so `sum(P1_Pspatial) == 0` fires
expr["G1", ] <- 0
expr["G2", cell_group != "g2"] <- 0
expr["G3", ] <- 1
expr["G4", ] <- 0.5
if ("IL12A" %in% genes) expr["IL12A", ] <- 0
if ("IL12B" %in% genes) expr["IL12B", ] <- 0
data.signaling <- as(expr, "dgCMatrix")

cat(sprintf("fixture: %d cells x %d genes, %d groups, %d L-R pairs, %d subunits\n",
            NC, NGENES, NG, nrow(LRsig), length(subunits)), file = stderr())
cat(sprintf("genes: %s\n", paste(genes, collapse = ",")), file = stderr())

## Pre-flight: every name the kernel will resolve must be present. Without this the fixture
## dies later with upstream's `subscript out of bounds` (R-ism 11) and gives no hint that
## the fault is the fixture rather than the port. Cofactor names themselves are looked up
## in `DB$cofactor` by `computeExpr_coreceptor`, not in the matrix, so they are exempt.
stopifnot(all(subunits %in% genes), all(plain_lr %in% genes))

mk <- function(datatype = "RNA") {
  o <- new("MiniCellChat",
           data.signaling = data.signaling,
           LR = list(LRsig = LRsig),
           LRsig = LRsig,
           DB = list(complex = DB$complex, cofactor = DB$cofactor),
           idents = cell_group,
           options = list(datatype = datatype))
  o
}

## ------------------------------------------------------------------ configs
## The matrix is deliberately spread over every branch that changes a bit: the four
## `type.mean` values, the population.size switch, the Hill parameters, and nboot.
configs <- list(
  list(name = "tri_ps0",    type = "triMean",          trim = 0.1, pop = FALSE, nboot = 5,  seed = 1L, Kh = 0.5, n = 1),
  list(name = "tri_ps1",    type = "triMean",          trim = 0.1, pop = TRUE,  nboot = 5,  seed = 1L, Kh = 0.5, n = 1),
  list(name = "trim_ps0",   type = "truncatedMean",    trim = 0.2, pop = FALSE, nboot = 7,  seed = 2L, Kh = 0.5, n = 1),
  list(name = "thresh_ps1", type = "thresholdedMean",  trim = 0.3, pop = TRUE,  nboot = 4,  seed = 3L, Kh = 0.5, n = 1),
  list(name = "median_ps0", type = "median",           trim = 0.1, pop = FALSE, nboot = 6,  seed = 4L, Kh = 0.5, n = 1),
  list(name = "hill2",      type = "triMean",          trim = 0.1, pop = FALSE, nboot = 5,  seed = 1L, Kh = 0.5, n = 2),
  list(name = "kh1e3",      type = "triMean",          trim = 0.1, pop = TRUE,  nboot = 8,  seed = 7L, Kh = 1e3, n = 1),
  list(name = "nboot1",     type = "triMean",          trim = 0.1, pop = FALSE, nboot = 1,  seed = 5L, Kh = 0.5, n = 1)
)

q <- file("tests/fixtures/prob_golden.txt", "wt")
fmt <- function(v) paste(format(as.numeric(v), digits = 17), collapse = ",")
## `format()` right-pads to a common width, so the reader must trim; `sprintf("%.17g")`
## does not. Use the same helper on both sides to keep the fixture unambiguous.
fmt2 <- function(v) paste(sprintf("%.17g", as.numeric(v)), collapse = ",")

## ------------------------------------------------------------------ inputs
## The Rust test cannot re-derive the fixture (same lesson as `expr_matrix.tsv`: the
## expression matrix is dumped, not reconstructed), and the L-R table is hand-built, so it
## is recorded here and mirrored in Rust under an explicit equality test.
qf <- file("tests/fixtures/prob_inputs.tsv", "wt")
cat(sprintf("cells\t%d\ngenes\t%d\ngroups\t%d\n", NC, NGENES, NG), file = qf)
writeLines(genes, qf)
cat("groups_labels\n", file = qf)
writeLines(levels(cell_group), qf)
cat("cell_group\n", file = qf)
writeLines(as.character(cell_group), qf)
cat("data\n", file = qf)
## `t(data.use)`: cells x genes, as the kernel consumes it. `%a` is exact hex.
writeLines(apply(as.matrix(t(as.matrix(data.signaling) / max(as.matrix(data.signaling)))), 1,
                 function(r) paste(sprintf("%a", r), collapse = " ")), qf)
cat("lr\n", file = qf)
for (i in seq_len(nrow(LRsig))) {
  cat(sprintf("%s\t%s\t%s\t%s\t%s\t%s\t%s\n", LRsig$ligand[i], LRsig$receptor[i],
              LRsig$agonist[i], LRsig$antagonist[i], LRsig$co_A_receptor[i],
              LRsig$co_I_receptor[i], rownames(LRsig)[i]), file = qf)
}
close(qf)

## The observed pass, recorded so the Rust side's `avg` can be checked against R's rather
## than assumed to be right (it is the same `aggregate` call as the bootstrap, so a bug
## there would otherwise cancel out on both sides).
q2 <- file("tests/fixtures/prob_avg.txt", "wt")
for (cfg in configs) {
  fu <- switch(cfg$type, triMean = function(x) env$triMean(x),
               truncatedMean = function(x) mean(x, trim = cfg$trim, na.rm = TRUE),
               thresholdedMean = function(x) env$thresholdedMean(x, trim = cfg$trim, na.rm = TRUE),
               median = function(x) median(x, na.rm = TRUE))
  du <- as.matrix(data.signaling)
  du <- du / max(du)
  au <- t(stats::aggregate(t(du), list(cell_group), FUN = fu)[, -1])
  cat(sprintf("avg\t%s\t%s\n", cfg$name, fmt2(au)), file = q2)
}
close(q2)

for (cfg in configs) {
  obj <- mk()
  res <- tryCatch(
    suppressWarnings(suppressMessages(env$computeCommunProb(
      obj,
      type = cfg$type, trim = cfg$trim, population.size = cfg$pop,
      nboot = cfg$nboot, seed.use = cfg$seed, Kh = cfg$Kh, n = cfg$n
    ))),
    error = function(e) e
  )
  if (inherits(res, "condition")) {
    cat(sprintf("config\t%s\tERROR\t%s\n", cfg$name, conditionMessage(res)), file = q)
    next
  }
  p <- res@net$prob; v <- res@net$pval
  cat(sprintf("config\t%s\t%dx%dx%d\n", cfg$name, dim(p)[1], dim(p)[2], dim(p)[3]), file = q)
  cat(sprintf("dimnames\t%s\t%s\t%s\n", cfg$name,
              paste(dimnames(p)[[1]], collapse = ","),
              paste(dimnames(p)[[3]], collapse = ",")), file = q)
  cat(sprintf("prob\t%s\t%s\n", cfg$name, fmt(p)), file = q)
  cat(sprintf("pval\t%s\t%s\n", cfg$name, fmt(v)), file = q)
  ## `options$parameter` is part of the drop-in contract (objective 2), so record it.
  cat(sprintf("param\t%s\t%s\n", cfg$name,
              paste(names(res@options$parameter), vapply(res@options$parameter, function(z)
                paste(format(z, digits = 15), collapse = "|"), ""), collapse = "\t")), file = q)
}
close(q)
cat("wrote tests/fixtures/prob_golden.txt\n", file = stderr())
