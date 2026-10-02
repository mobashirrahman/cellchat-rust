# Golden corpus for r-core's de.rs: evaluates upstream `computeAveExpr`, `subsetDB` and
# `subsetData` (sourced verbatim from the pinned commit).
#
# `subsetData` and `subsetDB` mutate S4 / CellChatDB objects, so the oracle calls them on
# minimal ones and records the *resulting* structure -- the reordering, the kept rows, the
# gene list -- rather than the object.
suppressWarnings(suppressMessages({library(collapse); library(dplyr)}))

CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
env <- new.env()
for (f in c("modeling.R", "utilities.R", "database.R")) {
  sys.source(file.path(CC, "R", f), envir = env, keep.source = FALSE)
}

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)

setClass("MiniAve", representation(
  data = "matrix", data.signaling = "ANY", DB = "ANY", idents = "factor",
  meta = "list", var.features = "list", options = "list"
))

fmt <- function(v) paste(sprintf("%.17g", as.numeric(v)), collapse = ",")
q <- file("tests/fixtures/de_golden.txt", "wt")

## ------------------------------------------------------------------ fixture data
set.seed(20240303)
NGENE <- 24L
NCELL <- 60L
genes <- paste0("G", seq_len(NGENE))
m <- matrix(runif(NGENE * NCELL, 0, 5), nrow = NGENE, dimnames = list(genes, NULL))
m[1, ] <- 0                                  # all zero
m[2, seq_len(30)] <- 0                       # zero in the first half
m[3, ] <- 2                                  # constant
m[4, seq(1, NCELL, by = 3)] <- 0             # sparse, so medians and trimeans differ
m[5, seq_len(NCELL)] <- NA_real_              # a fully missing row
cell_group <- factor(
  c(rep("g1", 25), rep("g2", 20), rep("g3", 15))[sample.int(NCELL)],
  levels = c("g1", "g2", "g3")
)
cat(sprintf("fixture: %d genes x %d cells, %d groups\n", NGENE, NCELL, nlevels(cell_group)),
    file = stderr())

mk <- function() {
  new("MiniAve", data = m, DB = DB, idents = cell_group, meta = list(),
      var.features = list(), options = list(datatype = "RNA"))
}

## ------------------------------------------------------------------ computeAveExpr
## `data.use` is genes x cells; the function transposes, aggregates and transposes back.
ave_specs <- list(
  list(name = "tri",    type = "triMean",        trim = NULL, features = NULL),
  list(name = "trim01", type = "truncatedMean",  trim = 0.1,  features = NULL),
  list(name = "trim05", type = "truncatedMean",  trim = 0.5,  features = NULL),
  list(name = "median", type = "median",         trim = NULL, features = NULL),
  ## `intersect(features, rownames)` keeps the *caller's* order, and duplicates collapse.
  list(name = "featsel", type = "triMean",       trim = NULL,
       features = c("G10", "G1", "G10", "GZZ", "G2")),
  ## `match.arg` partial matching: "tri" is a unique prefix of "triMean".
  list(name = "matcharg", type = "tri",          trim = NULL, features = NULL),
  ## A group with a single cell, and a group with a single expressed cell.
  list(name = "median_missing", type = "median", trim = NULL, features = c("G5", "G1", "G4"))
)

for (spec in ave_specs) {
  o <- mk()
  out <- tryCatch(
    suppressWarnings(suppressMessages(env$computeAveExpr(
      o, features = spec$features, type = spec$type, trim = spec$trim,
      slot.name = "data"))),
    error = function(e) structure(conditionMessage(e), class = "rerr"))
  if (inherits(out, "rerr")) {
    cat(sprintf("ave\t%s\tERROR\t%s\n", spec$name, out), file = q)
    next
  }
  cat(sprintf("ave\t%s\t%sx%d\n", spec$name, nrow(out), ncol(out)), file = q)
  cat(sprintf("ave_rows\t%s\t%s\n", spec$name, paste(rownames(out), collapse = ",")), file = q)
  cat(sprintf("ave_cols\t%s\t%s\n", spec$name, paste(colnames(out), collapse = ",")), file = q)
  cat(sprintf("ave_vals\t%s\t%s\n", spec$name, fmt(out)), file = q)
}

## `match.arg` must reject an ambiguous prefix with upstream's message.
amb <- tryCatch(env$computeAveExpr(mk(), type = "t"), error = conditionMessage)
cat(sprintf("ave_argerr\t%s\n", amb), file = q)

## ------------------------------------------------------------------ subsetData
## Needs `extractGene` from utilities.R, and a `data.signaling` slot.
sub_specs <- list(
  list(name = "default", features = NULL),
  list(name = "explicit", features = c("G7", "G1", "G7", "NOPE")),
  list(name = "empty", features = character(0))
)
## The pre-`subsetData` annotation shuffling is part of what we want to pin, so record the
## DB's annotation order before and after.
int0 <- DB$interaction
cat(sprintf("db_int_rows\t%d\n", nrow(int0)), file = q)
cat(sprintf("db_int_ann0\t%s\n", paste(int0$annotation, collapse = ",")), file = q)

for (spec in sub_specs) {
  o <- mk()
  o2 <- tryCatch(env$subsetData(o, features = spec$features),
                 error = function(e) structure(conditionMessage(e), class = "rerr"))
  if (inherits(o2, "rerr")) {
    cat(sprintf("subset_data\t%s\tERROR\t%s\n", spec$name, o2), file = q)
    next
  }
  cat(sprintf("subset_data\t%s\t%dx%d\n", spec$name,
              nrow(o2@data.signaling), ncol(o2@data.signaling)), file = q)
  cat(sprintf("subset_data_rows\t%s\t%s\n", spec$name,
              paste(rownames(o2@data.signaling), collapse = ",")), file = q)
  cat(sprintf("subset_db_ann\t%s\t%s\n", spec$name,
              paste(o2@DB$interaction$annotation, collapse = ",")), file = q)
  cat(sprintf("subset_db_rows\t%d\n", nrow(o2@DB$interaction)), file = q)
}

## ------------------------------------------------------------------ subsetDB
subdb_specs <- list(
  list(name = "default", search = NULL, non_protein = FALSE),
  list(name = "nonprotein", search = NULL, non_protein = TRUE),
  list(name = "explicit_all", search = c("Secreted Signaling", "ECM-Receptor",
                                         "Cell-Cell Contact", "Non-protein Signaling"),
       non_protein = FALSE),
  list(name = "contact_only", search = "Cell-Cell Contact", non_protein = FALSE),
  list(name = "empty_search", search = character(0), non_protein = FALSE)
)
for (spec in subdb_specs) {
  d <- DB
  out <- tryCatch(
    suppressMessages(env$subsetDB(d, search = spec$search, key = "annotation",
                                  non_protein = spec$non_protein)),
    error = function(e) structure(conditionMessage(e), class = "rerr"))
  if (inherits(out, "rerr")) {
    cat(sprintf("subset_db\t%s\tERROR\t%s\n", spec$name, out), file = q)
    next
  }
  cat(sprintf("subset_db\t%s\t%d\n", spec$name, nrow(out$interaction)), file = q)
  cat(sprintf("subset_db_annset\t%s\t%s\n", spec$name,
              paste(sort(unique(out$interaction$annotation)), collapse = ",")), file = q)
  cat(sprintf("subset_db_first\t%s\t%s\n", spec$name,
              paste(head(out$interaction$interaction_name, 3), collapse = ",")), file = q)
}
## A bad key must give upstream's message.
badkey <- tryCatch(env$subsetDB(DB, search = "Secreted Signaling", key = "nope"),
                   error = conditionMessage)
cat(sprintf("subset_db_keyerr\t%s\n", badkey), file = q)

close(q)
cat("wrote tests/fixtures/de_golden.txt\n", file = stderr())
