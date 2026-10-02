# Golden corpus for r-core's wilcox.rs: evaluates upstream `stats::wilcox.test`,
# `stats::p.adjust`, `stats::rank` and `identifyOverExpressedGenes(do.fast = FALSE)` from
# the pinned commit.
suppressWarnings(suppressMessages({library(Matrix); library(collapse); library(dplyr)}))

CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
env <- new.env()
for (f in c("modeling.R", "utilities.R", "database.R")) {
  sys.source(file.path(CC, "R", f), envir = env, keep.source = FALSE)
}

setClass("MiniDE", representation(
  data = "matrix", data.signaling = "ANY", DB = "ANY", idents = "factor",
  meta = "list", var.features = "list", options = "list"
))

fmt <- function(v) paste(sprintf("%.17g", as.numeric(v)), collapse = ",")
q <- file("tests/fixtures/wilcox_golden.txt", "wt")

## ------------------------------------------------------------------ wilcox.test
## Small, awkward pairs: ties, zeroes, one- and two-sided, both branches.
set.seed(20240404)
w_cases <- list(
  list(x = c(1, 2, 3, 4),           y = c(5, 6, 7, 8)),          # exact, no ties
  list(x = c(1, 2, 3, 4, 5),        y = c(2, 3, 4, 5, 6)),       # ties -> normal
  list(x = c(0, 0, 0, 1),           y = c(0, 1, 1, 1)),          # heavy ties + zeroes
  list(x = c(1, 1),                 y = c(2, 2)),                 # minimal, all ties
  list(x = c(5),                    y = c(1, 2, 3, 4)),           # n.x = 1
  list(x = c(1, 1, 2, 2, 3, 3, 4, 4), y = c(1, 2, 2, 3, 3, 4, 4, 5)),
  list(x = rnorm(30),               y = rnorm(30)),                # exact-ish, no ties
  list(x = rnorm(60),               y = rnorm(60)),                # normal branch, big
  list(x = c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10), y = c(11, 12, 13, 14, 15)),
  ## a perfectly balanced split: W - n.x*n.y/2 == 0, so sign(z) == 0 and the continuity
  ## correction is zero. A discontinuity that a naive port gets wrong.
  list(x = c(1, 2, 3, 4),           y = c(1, 2, 3, 4)),
  list(x = c(2, 2, 2),              y = c(1, 1, 1)),
  ## NA handling: x drops NAs, y drops NAs, and they are dropped *separately*.
  list(x = c(1, 2, NA, 4, 5),       y = c(3, NA, 5, 6, 7)),
  list(x = c(NA_real_, NA_real_),   y = c(1, 2, 3))
)
for (i in seq_along(w_cases)) {
  cs <- w_cases[[i]]
  for (alt in c("two.sided", "greater", "less")) {
    for (corr in c(TRUE, FALSE)) {
      out <- tryCatch(
        suppressWarnings(stats::wilcox.test(cs$x, cs$y, alternative = alt, correct = corr)),
        error = function(e) structure(conditionMessage(e), class = "rerr"))
      if (inherits(out, "rerr")) {
        cat(sprintf("wilcox\t%d\t%s\t%s\tERROR\t%s\n", i, alt, corr, out), file = q)
        next
      }
      ## NB an `htest` object has no `$correct`: `names(out)` is exactly
      ## `statistic parameter p.value null.value alternative method data.name`. Asking for
      ## `$correct` gives NULL, `as.integer(NULL)` is `integer(0)`, and `sprintf("%s", ...)`
      ## then returns `character(0)` -- so the whole record collapses to a blank line
      ## rather than to an error. The branch is already visible in `$method`, which says
      ## "exact test" or "with continuity correction", so that is what is recorded.
      stopifnot(is.null(out$correct))
      cat(sprintf("wilcox\t%d\t%s\t%s\t%s\t%s\t%s\t%s\n", i, alt, corr,
                  out$method, format(out$statistic, digits = 17),
                  format(out$p.value, digits = 17), out$alternative, out$null.value),
          file = q)
    }
  }
}

## ------------------------------------------------------------------ rank and p.adjust
cat(sprintf("rank\t%s\t%s\n",
            paste(sapply(list(c(1, 2, 3), c(1, 1, 2), c(5, 5, 5, 5), c(3, 1, 2),
                              c(1, NA, 2), c(2, 1, 1, 3, 3, 3)),
                         function(v) paste(rank(v), collapse = "/")), collapse = ","),
            ""), file = q)
for (nm in c("bonferroni", "holm", "none")) {
  p <- c(0.001, 0.008, 0.039, 0.041, 0.042, 0.06, 0.074, 0.205, 0.212, 0.216)
  cat(sprintf("padjust\t%s\t%s\t%s\n", nm, length(p),
              paste(sprintf("%.17g", p.adjust(p, method = nm)), collapse = ",")), file = q)
}

## ------------------------------------------------------------------ identifyOverExpressedGenes
## A genes x cells matrix with the structure the function branches on: all-zero rows, rows
## expressed in one group only, constant rows, and a fully missing row.
set.seed(20240405)
NG <- 30L
NC <- 80L
genes <- paste0("G", seq_len(NG))
m <- matrix(round(rpois(NG * NC, 3), 3), nrow = NG, dimnames = list(genes, NULL))
m[1, ] <- 0                          # all zero
m[2, seq_len(NC / 2)] <- 0          # zero in the first half only
m[3, ] <- 1                          # constant
grp <- factor(c(rep("g1", 30), rep("g2", 30), rep("g3", 20))[sample.int(NC)],
              levels = c("g1", "g2", "g3"))
cat(sprintf("fixture: %d genes x %d cells, %d groups\n", NG, NC, nlevels(grp)), file = stderr())

## Dump `X` and the group vector. The Rust test used to re-derive `rpois(3)` from R's
## Poisson RNG, which is a fourth re-derivation of a fixture in a test -- exactly the
## failure mode `tests/parity/README.md` warns about, and it failed: the inversion loop's
## `unif_rand()` accounting differs from the Knuth product method. Ground truth in a file.
qm <- file("tests/fixtures/wilcox_matrix.tsv", "wt")
cat(sprintf("genes\t%d\ncells\t%d\n", NG, NC), file = qm)
writeLines(genes, qm)
cat("groups\n", file = qm)
writeLines(levels(grp), qm)
cat("cell_group\n", file = qm)
writeLines(as.character(grp), qm)
cat("data\n", file = qm)
writeLines(apply(X = m, MARGIN = 1, FUN = function(r) paste(sprintf("%a", r), collapse = " ")),
           con = qm)
close(qm)

mk_de <- function() {
  o <- new("MiniDE", data = m, idents = grp, meta = list(), var.features = list(),
           options = list(datatype = "RNA", mode = "single"))
  o@data.signaling <- as(o@data, "dgCMatrix")
  o
}

de_specs <- list(
  list(name = "default",  thresh_pc = 0,   thresh_fc = 0, thresh_p = 0.05,
       only_pos = TRUE,  do_DE = TRUE,  min_cells = 10),
  ## `thresh.pc` is compared against a *percentage* in the final filter and a *fraction* in
  ## the selection step, so a non-zero value is the interesting case.
  list(name = "pc01",     thresh_pc = 0.1, thresh_fc = 0, thresh_p = 0.05,
       only_pos = TRUE,  do_DE = TRUE,  min_cells = 10),
  list(name = "fc005",    thresh_pc = 0,   thresh_fc = 0.5, thresh_p = 0.1,
       only_pos = TRUE,  do_DE = TRUE,  min_cells = 10),
  list(name = "twosided", thresh_pc = 0,   thresh_fc = 0, thresh_p = 0.05,
       only_pos = FALSE, do_DE = TRUE,  min_cells = 10),
  list(name = "p001",     thresh_pc = 0,   thresh_fc = 0, thresh_p = 0.001,
       only_pos = TRUE,  do_DE = TRUE,  min_cells = 10),
  ## `do.DE = FALSE` skips the Wilcoxon entirely and selects on the percentage only.
  list(name = "node",     thresh_pc = 0.2, thresh_fc = 0, thresh_p = 1,
       only_pos = TRUE,  do_DE = FALSE, min_cells = 10),
  ## A high min.cells that every group fails.
  list(name = "mincells", thresh_pc = 0,   thresh_fc = 0, thresh_p = 0.05,
       only_pos = TRUE,  do_DE = TRUE,  min_cells = 1000)
)

for (spec in de_specs) {
  o <- mk_de()
  out <- tryCatch(
    suppressWarnings(suppressMessages(env$identifyOverExpressedGenes(
      o, group.by = NULL, idents.use = NULL, invert = FALSE,
      features.name = "features", only.pos = spec$only_pos, features = NULL,
      return.object = TRUE, thresh.pc = spec$thresh_pc, thresh.fc = spec$thresh_fc,
      thresh.p = spec$thresh_p, do.DE = spec$do_DE, do.fast = FALSE,
      min.cells = spec$min_cells))),
    error = function(e) structure(conditionMessage(e), class = "rerr"))
  if (inherits(out, "rerr")) {
    cat(sprintf("de\t%s\tERROR\t%s\n", spec$name, out), file = q)
    next
  }
  info <- out@var.features[["features.info"]]
  cat(sprintf("de\t%s\t%dx%d\n", spec$name, nrow(info), ncol(info)), file = q)
  cat(sprintf("de_cols\t%s\t%s\n", spec$name, paste(colnames(info), collapse = ",")), file = q)
  ## The row names are part of the object and `identical()` checks them. They are the
  ## *feature names*, in both branches, and for two different reasons:
  ##  * `do.DE = TRUE`:  `data.frame(..., data.alpha[features, , drop = FALSE], ...)` -- the
  ##    unnamed matrix argument carries `rownames == features` and `data.frame()` adopts them.
  ##  * `do.DE = FALSE`: `data.frame(features = <unnamed chr>, nCells = rowSums(data.use > 0))`
  ##    -- `rowSums` on a `Matrix` returns a *named* vector and `data.frame()` adopts the first
  ##    named argument's names. (With a dense `data.use` it would be 1..nrow instead, so this
  ##    row naming is sparse-vs-dense dependent upstream.)
  cat(sprintf("de_rownames\t%s\t%s\n", spec$name,
              paste(rownames(info), collapse = ",")), file = q)
  ## The marker list itself, in the order upstream produces, plus the `features` vector.
  cat(sprintf("de_features\t%s\t%s\n", spec$name,
              paste(out@var.features[["features"]], collapse = ",")), file = q)
  for (cc in colnames(info)) {
    v <- info[[cc]]
    cat(sprintf("de_col\t%s\t%s\t%s\n", spec$name, cc,
                paste(ifelse(is.na(v), "NA", format(v, digits = 17)), collapse = ",")), file = q)
  }
  ## One Wilcoxon p-value per gene per group, for *every* gene, so the fixture covers the
  ## genes that are filtered out as well as the ones that survive.
  lev <- levels(grp)
  for (i in seq_along(lev)) {
    cells <- which(as.character(grp) == lev[i])
    rest <- setdiff(seq_len(NC), cells)
    pv <- vapply(seq_len(NG), function(g) {
      stats::wilcox.test(m[g, cells], m[g, rest])$p.value
    }, numeric(1))
    cat(sprintf("de_pval\t%s\t%s\t%s\n", spec$name, lev[i], paste(sprintf("%.17g", pv), collapse = ",")),
        file = q)
  }
  ## And the intermediate quantities the Rust side computes directly.
  for (i in seq_along(lev)) {
    cells <- which(as.character(grp) == lev[i])
    rest <- setdiff(seq_len(NC), cells)
    pct1 <- round(rowSums(m[, cells, drop = FALSE] > 0) / length(cells), 3)
    pct2 <- round(rowSums(m[, rest, drop = FALSE] > 0) / length(rest), 3)
    mean.fxn <- function(x) log(x = mean(x = expm1(x = x)) + 1)
    d1 <- apply(m[, cells, drop = FALSE], 1, mean.fxn)
    d2 <- apply(m[, rest, drop = FALSE], 1, mean.fxn)
    cat(sprintf("de_pct1\t%s\t%s\t%s\n", spec$name, lev[i], paste(sprintf("%.17g", pct1), collapse = ",")), file = q)
    cat(sprintf("de_pct2\t%s\t%s\t%s\n", spec$name, lev[i], paste(sprintf("%.17g", pct2), collapse = ",")), file = q)
    cat(sprintf("de_fc\t%s\t%s\t%s\n", spec$name, lev[i], paste(sprintf("%.17g", d1 - d2), collapse = ",")), file = q)
  }
}

close(q)
cat("wrote tests/fixtures/wilcox_golden.txt\n", file = stderr())
