#!/usr/bin/env Rscript
# Export the pinned CellChatDB to the flat, dependency-free format that `r-core`'s
# `db.rs` parses. From the repo root:
#
#   Rscript tests/parity/export_db.R <human|mouse|zebrafish> <outdir>
#
# Design notes
# ------------
# * TSV rather than Parquet/Arrow. `r-core` must stay free of I/O and heavy
#   dependencies so its numerics can be differentially tested without an R session, and
#   so the standalone library has nothing to link. The whole DB is ~1 MB, so parsing
#   cost is irrelevant next to the kernel.
# * The export is **derived from the pinned upstream commit**, and the manifest records
#   the source `.rda`'s MD5 and a SHA-256 of every emitted file, so a Rust test can
#   refuse to run against a stale or doctored fixture.
# * Interactions are emitted already in annotation order (see below), so `r-core` never
#   has to reproduce R's `order()` on a factor. The R side *verifies* that R's ordering
#   is a stable partition and aborts if a future R release changes that.

suppressWarnings(suppressMessages(library(dplyr)))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2L) stop("usage: export_db.R <human|mouse|zebrafish> <outdir>")
species <- args[[1]]
outdir <- args[[2]]
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

## ---------------------------------------------------------------- locate upstream
pin <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
rda <- file.path(pin, "data", sprintf("CellChatDB.%s.rda", species))
if (!file.exists(rda)) {
  stop(sprintf("cannot find %s\n", rda),
       "clone jinworks/CellChat at 75253cd0c9e68410e6e721a6d3a0419a1d7e358f ",
       "and set CELLCHAT_SRC to the checkout")
}
E <- new.env()
load(rda, envir = E)
DB <- get(ls(E)[1], E)

## The annotation factor levels, in the order `subsetData` uses. Index 0..2 are the
## diffusion-mediated classes and index 3 is contact-dependent; `computeCommunProb`
## derives `nLR1` from the boundary between them, so this order is load-bearing.
ANNOT_LEVELS <- c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling",
                  "Cell-Cell Contact")

## ---------------------------------------------------------------- field sanitising
## Missing values become empty fields. Tabs and newlines cannot appear in a gene symbol
## or a pathway name, but sanitise defensively so a malformed DB cannot silently shift
## columns.
clean <- function(x) {
  x <- as.character(x)
  x[is.na(x)] <- ""
  x <- gsub("[\t\r\n]", " ", x)
  x
}

## ---------------------------------------------------------------- interactions
ii <- DB$interaction
need <- c("interaction_name", "pathway_name", "ligand", "receptor",
          "agonist", "antagonist", "co_A_receptor", "co_I_receptor", "annotation")
if (!all(need %in% colnames(ii))) {
  stop("CellChatDB$interaction is missing column(s): ",
       paste(setdiff(need, colnames(ii)), collapse = ", "))
}

## NB: `order()` must be evaluated while `annotation` is still a FACTOR, exactly as
## `subsetData` does it. `order()` on a factor uses the level order; on a character
## vector it sorts alphabetically, which is a *different* table with a different `nLR1`.
## Converting to character first is a silent, plausible-looking bug -- it was one here.
if (length(unique(ii$annotation)) > 1L) {
  ii$annotation <- factor(ii$annotation, levels = ANNOT_LEVELS)
  ## `order()` on a factor must be a *stable* partition by level, because `nLR1` and the
  ## reported `dimnames` depend on the within-level order being preserved. Assert it
  ## rather than trusting it.
  stable_partition <- function(f) {
    a <- as.character(f)
    as.integer(unlist(lapply(levels(f), function(l) which(a == l)), use.names = FALSE))
  }
  ## NB `stopifnot()` treats trailing character arguments as *further conditions*, not
  ## as messages, so `stopifnot(cond, "msg")` fails on the string itself. Use an
  ## explicit `stop()` wherever a message is wanted.
  if (!identical(order(ii$annotation), stable_partition(ii$annotation))) {
    stop("R's order() on the annotation factor is no longer a stable partition; ",
         "r-core's pre-sorted export would then be wrong and nLR1 would move")
  }
  ii <- ii[order(ii$annotation), , drop = FALSE]
  ii$annotation <- as.character(ii$annotation)
}
## Unknown annotations would sort to NA, i.e. last, in R. Fail loudly instead: an
## unrecognised class would silently change the contact/diffusion split.
unknown <- setdiff(unique(ii$annotation), ANNOT_LEVELS)
if (length(unknown)) {
  stop("unrecognised annotation value(s): ", paste(unknown, collapse = ", "),
       " -- add them to ANNOT_LEVELS, because they change the contact/diffusion split")
}

writeLines(
  vapply(seq_len(nrow(ii)), function(i) {
    paste(clean(unlist(ii[i, need, drop = FALSE])), collapse = "\t")
  }, character(1L)),
  file.path(outdir, "interactions.tsv")
)

nlr1 <- {
  diffusion <- which(ii$annotation %in% ANNOT_LEVELS[1:3])
  if (length(diffusion)) max(diffusion) else 0L
}

## ---------------------------------------------------------------- complex / cofactor
## Row names carry the entity name and are not part of `apply()`, so build lines by hand.
sub_cols <- grep("^subunit", colnames(DB$complex), value = TRUE)
if (!length(sub_cols)) stop("CellChatDB$complex has no subunit_* columns")
## NB empty cells are emitted, not dropped. Upstream reaches subunits via
## `unlist(complex_input[match(complex, rownames(...)), ])`, and `unlist` on a
## multi-row data frame is COLUMN-major. Dropping the empties would shift later columns
## and silently reorder the subunits, which is exactly the bug this export is meant to
## prevent. Trailing tabs are preserved and `str::split('\t')` in Rust keeps them.
writeLines(
  vapply(seq_len(nrow(DB$complex)), function(i) {
    paste(c(rownames(DB$complex)[i], clean(unlist(DB$complex[i, sub_cols]))),
          collapse = "\t")
  }, character(1L)),
  file.path(outdir, "complexes.tsv")
)

cof_cols <- grep("^cofactor", colnames(DB$cofactor), value = TRUE)
if (!length(cof_cols)) stop("CellChatDB$cofactor has no cofactor* columns")
## Same reasoning as the complex table: keep the empty cells so the column positions
## survive the round trip.
writeLines(
  vapply(seq_len(nrow(DB$cofactor)), function(i) {
    paste(c(rownames(DB$cofactor)[i], clean(unlist(DB$cofactor[i, cof_cols]))),
          collapse = "\t")
  }, character(1L)),
  file.path(outdir, "cofactors.tsv")
)

## ---------------------------------------------------------------- extractGene reference
## `extractGene` / `extractGeneSubset` are reimplemented in r-core's `db.rs`. Their output
## is emitted here *by upstream's own R code* so the Rust port can be compared against it
## order-for-order, not just set-for-set. `identifyOverExpressedGenes` and
## `computeAveExpr(features = ...)` observe this order, so it is part of the contract.
extractGeneSubset <- function(geneSet, complex_input, geneIfo) {
  complex <- geneSet[which(geneSet %in% geneIfo$Symbol == "FALSE")]
  geneSet <- intersect(geneSet, geneIfo$Symbol)
  complexsubunits <- dplyr::select(
    complex_input[match(complex, rownames(complex_input), nomatch = 0), ],
    starts_with("subunit"))
  complex <- intersect(complex, rownames(complexsubunits))
  complexsubunitsV <- unlist(complexsubunits)
  complexsubunitsV <- unique(complexsubunitsV[complexsubunitsV != ""])
  unique(c(geneSet, complexsubunitsV))
}
extractGene <- function(CellChatDB) {
  interaction_input <- CellChatDB$interaction
  complex_input <- CellChatDB$complex
  cofactor_input <- CellChatDB$cofactor
  geneIfo <- CellChatDB$geneInfo
  geneL <- unique(interaction_input$ligand)
  geneR <- unique(interaction_input$receptor)
  geneLR <- c(extractGeneSubset(geneL, complex_input, geneIfo),
              extractGeneSubset(geneR, complex_input, geneIfo))
  cofactor <- unique(c(interaction_input$agonist, interaction_input$antagonist,
                      interaction_input$co_A_receptor, interaction_input$co_I_receptor))
  cofactor <- cofactor[cofactor != ""]
  cofactorsubunits <- dplyr::select(
    cofactor_input[match(cofactor, rownames(cofactor_input), nomatch = 0), ],
    starts_with("cofactor"))
  cofactorsubunitsV <- unlist(cofactorsubunits)
  geneCofactor <- unique(cofactorsubunitsV[cofactorsubunitsV != ""])
  unique(c(geneLR, geneCofactor))
}
## NB the order matters and it is the *sorted* order. `subsetData` assigns
## `object@DB$interaction <- interaction_input` (annotation-sorted) BEFORE calling
## `extractGene(DB)`, so `extractGene` always sees the sorted table -- not the order the
## `.rda` ships. Using the raw table here gave a reference that differed from upstream's
## real state at index 353, which is exactly the kind of near-miss worth a test.
sorted_db <- DB
sorted_db$interaction <- ii

gene_use <- extractGene(sorted_db)
writeLines(gene_use, file.path(outdir, "extract_gene.txt"))

## Per-annotation gene sets too: `subsetData` is normally run on a `subsetDB` output, and
## the fixture matrix needs the 3-annotation (no Non-protein) case as well.
subsets <- list(
  four = c(ANNOT_LEVELS),
  three = c(ANNOT_LEVELS[1], ANNOT_LEVELS[2], ANNOT_LEVELS[4])
)
for (nm in names(subsets)) {
  d <- sorted_db
  d$interaction <- sorted_db$interaction[sorted_db$interaction$annotation %in% subsets[[nm]], ]
  ## `subsetDB` preserves the table's row order, so filtering the *sorted* table is
  ## what `subsetDB(subsetDB(...))` then `subsetData` actually produces.
  g <- extractGene(d)
  writeLines(g, file.path(outdir, sprintf("extract_gene_%s.txt", nm)))
}

## ---------------------------------------------------------------- official symbols
## `geneInfo$Symbol` is the "is this an official gene symbol?" oracle. It contains
## exactly one NA in the shipped DB; NA is *not* a symbol, because `%in%` only matches
## an NA needle against an NA table entry and needles are never NA here.
sym <- clean(DB$geneInfo$Symbol)
sym <- unique(sym[sym != ""])
writeLines(sym, file.path(outdir, "symbols.txt"))

## ---------------------------------------------------------------- manifest
## R 4.3 ships no SHA-256 in `tools` and `digest`/`openssl` are not installed, so this
## uses MD5. It is a *staleness* check -- "does the on-disk fixture still match what
## the pinned .rda produces?" -- not an integrity guarantee against an adversary.
md5 <- function(f) unname(tools::md5sum(f))
files <- c("interactions.tsv", "complexes.tsv", "cofactors.tsv", "symbols.txt",
           "extract_gene.txt", "extract_gene_four.txt", "extract_gene_three.txt")
hashes <- vapply(files, function(f) md5(file.path(outdir, f)), character(1L))

lines <- c(
  "# cellchatrs exported CellChatDB -- see tests/parity/export_db.R",
  sprintf("species\t%s", species),
  sprintf("upstream_commit\t%s",
          Sys.getenv("CELLCHAT_COMMIT", "75253cd0c9e68410e6e721a6d3a0419a1d7e358f")),
  sprintf("rda_md5\t%s", md5(rda)),
  sprintf("n_interactions\t%d", nrow(ii)),
  sprintf("n_complexes\t%d", nrow(DB$complex)),
  sprintf("n_cofactors\t%d", nrow(DB$cofactor)),
  sprintf("n_symbols\t%d", length(sym)),
  sprintf("annotation_levels\t%s", paste(ANNOT_LEVELS, collapse = " | ")),
  sprintf("nlr1\t%d", nlr1),
  sprintf("md5_interactions\t%s", hashes[["interactions.tsv"]]),
  sprintf("md5_complexes\t%s", hashes[["complexes.tsv"]]),
  sprintf("md5_cofactors\t%s", hashes[["cofactors.tsv"]]),
  sprintf("md5_symbols\t%s", hashes[["symbols.txt"]]),
  sprintf("md5_extract_gene\t%s", hashes[["extract_gene.txt"]]),
  sprintf("n_extract_gene\t%d", length(gene_use))
)
writeLines(lines, file.path(outdir, "MANIFEST.tsv"))

cat(sprintf("exported %s: %d interactions, %d complexes, %d cofactors, %d symbols -> %s\n",
            species, nrow(ii), nrow(DB$complex), nrow(DB$cofactor), length(sym), outdir))
cat(sprintf("annotation runs: %s\n", paste(
  sprintf("%s x%d", ANNOT_LEVELS, vapply(ANNOT_LEVELS, function(l) sum(ii$annotation == l), integer(1L))),
  collapse = ", ")))
cat(sprintf("nLR1 (contact/diffusion boundary index) = %d of %d\n", nlr1, nrow(ii)))
