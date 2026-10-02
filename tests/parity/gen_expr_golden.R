# Golden corpus for r-core's expr.rs: evaluates the *upstream* computeExpr_* functions
# from R/modeling.R (sourced verbatim, not reimplemented) on a set of fixtures.
suppressWarnings(suppressMessages({library(Matrix); library(dplyr); library(collapse)}))

CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
env <- new.env()
source(file.path(CC, "R", "modeling.R"), local = env)

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)
complex_input <- DB$complex
cofactor_input <- DB$cofactor

# Fixture data: genes x K cell groups, with values engineered to exercise the branches.
#   - plain genes (0.1 .. 0.9)
#   - all-zero rows (the L1 zero short-circuit in the aggregation)
#   - a row with a single non-zero (produces a zero trimean)
#   - subunit rows for a real complex
set.seed(4242)
K <- 6L
## The gene set must include **every** complex and cofactor subunit, otherwise the
## `intersect(..., rownames(data.use))` inside computeExpr_coreceptor / _agonist /
## _antagonist silently degrades to the all-ones no-op and the fixture stops testing
## anything. That happened on the first attempt and looked like a pass.
all_subunits <- unique(c(
  unlist(complex_input[, grepl("^subunit", colnames(complex_input))]),
  unlist(cofactor_input[, grepl("^cofactor", colnames(cofactor_input))])
))
all_subunits <- all_subunits[!is.na(all_subunits) & all_subunits != ""]
genes <- unique(c(paste0("G", 1:12), all_subunits))
X <- matrix(runif(length(genes) * K, 0.01, 1), nrow = length(genes), ncol = K,
            dimnames = list(genes, paste0("c", 1:K)))
cat(sprintf("fixture: %d genes x %d groups (%d DB subunits)\n",
            nrow(X), K, length(all_subunits)), file = stderr())
X[1, ] <- 0                                   # all zero
X[2, c(1, 2, 3, 4, 5)] <- 0                  # >= 75% zero in every group
X[3, ] <- c(0.9, 0, 0, 0, 0, 0)              # sparse
X[4, ] <- 1                                   # constant

fmt <- function(v) paste(format(as.numeric(v), digits = 17), collapse = ",")
q <- file("tests/fixtures/expr_golden.txt", "wt")

## ---------------------------------------------------------------- the fixture matrix
## Dump `X` itself, in exact hex, rather than making the Rust test re-derive it from
## `set.seed`/`runif`. The re-derivation was a bug factory: it forced the test to
## reproduce R's `unlist(data.frame[, cols])` column-major order and the `unique()`
## first-appearance order, and when the Rust side got either wrong the fixture still had
## the right gene *set*, so every name resolved and every value looked plausible while
## being the value of a different gene. Ground truth in a file cannot rot that way.
## `%a` is C99 hex float: exactly representable, so the round trip is bit-exact.
## One connection for the whole file: opening a fresh `file(..., "at")` per call lets R
## flush them out of order (it happened, and the reader cannot tell a reordered file from
## a corrupted one). The section markers make the layout self-describing so the Rust
## reader never has to infer it.
xf <- file("tests/fixtures/expr_matrix.tsv", "wt")
cat(sprintf("groups\t%d\ngenes\t%d\n", K, length(genes)), file = xf)
writeLines(genes, xf)
cat("values\n", file = xf)
writeLines(apply(X, 1, function(r) paste(sprintf("%a", r), collapse = " ")), xf)
close(xf)

## ---------------------------------------------------------------- computeExpr_LR
## Names chosen to hit every branch of `index.singleL` / `index.complexL`.
lr_cases <- list(
  list("single",              c("G1", "G2")),
  list("complex",             c("Activin AB", "IL12AB")),
  list("mixed",               c("G1", "Activin AB", "IL12AB", "G4")),
  list("missing_subunit",     c("G1", "IL12AB")),   # fine here; both subunits present
  list("absent_not_complex",  c("G1", "NOT_A_COMPLEX")),
  list("repeat",             c("G1", "G1", "Activin AB", "Activin AB")),
  list("single_row",          "G5")
)
for (case in lr_cases) {
  nm <- case[[1]]; geneset <- case[[2]]
  out <- tryCatch(env$computeExpr_LR(geneset, X, complex_input),
                  error = function(e) structure(conditionMessage(e), class = "rerr"))
  if (inherits(out, "rerr")) {
    cat(sprintf("lr\t%s\tERROR\t%s\n", nm, out), file = q)
  } else {
    cat(sprintf("lr\t%s\t%s\t%s\n", nm, paste(dim(out), collapse = "x"), fmt(out)), file = q)
  }
}

## A complex whose subunit is absent from the matrix must reproduce upstream's
## "subscript out of bounds" (R-ism 11), not a graceful no-op.
Xm <- X[!(rownames(X) %in% c("IL12B")), , drop = FALSE]
out <- tryCatch(env$computeExpr_LR(c("IL12AB"), Xm, complex_input),
                error = function(e) conditionMessage(e))
cat(sprintf("lr\tmissing_subunit_row\t%s\t%s\n",
            if (is.character(out)) "ERROR" else "ok", out), file = q)

## ---------------------------------------------------------------- computeExpr_coreceptor
## NB these must be *cofactor row names* (e.g. "TGFb agonist"), not subunit names.
## Passing a subunit gives a missing row and therefore the all-ones no-op, which looks
## like a passing test but exercises nothing.
lr_df <- data.frame(
  co_A_receptor    = c("", "ACTIVIN antagonist", "TGFb antagonist", "", "NOT_A_COFACTOR"),
  co_I_receptor    = c("", "TGFb inhibition receptor", "", "", ""),
  stringsAsFactors = FALSE
)
## Emit the table itself so the Rust test can assert its hand-typed mirror still agrees.
for (col in c("co_A_receptor", "co_I_receptor")) {
  cat(sprintf("lr_df\t%s\t%s\n", col, paste(lr_df[[col]], collapse = "\t")), file = q)
}
for (ty in c("A", "I")) {
  out <- tryCatch(env$computeExpr_coreceptor(cofactor_input, X, lr_df, type = ty),
                  error = function(e) structure(conditionMessage(e), class = "rerr"))
  if (inherits(out, "rerr")) {
    cat(sprintf("coreceptor\t%s\tERROR\t%s\n", ty, out), file = q)
  } else {
    cat(sprintf("coreceptor\t%s\t%s\t%s\n", ty, paste(dim(out), collapse = "x"), fmt(out)),
        file = q)
  }
}
## Subset with a single row, as the bootstrap loop calls it.
for (i in c(2L, 3L, 5L)) {
  one <- lr_df[i, , drop = FALSE]
  out <- tryCatch(env$computeExpr_coreceptor(cofactor_input, X, one, type = "A"),
                  error = function(e) structure(conditionMessage(e), class = "rerr"))
  cat(sprintf("coreceptor_one\t%d\t%s\t%s\n", i,
              if (inherits(out, "rerr")) "ERROR" else fmt(out), ""), file = q)
}

## ---------------------------------------------------------------- agonist / antagonist
ag_cases <- list(
  list("ACTIVIN antagonist", "ACTIVIN inhibition receptor"),  # 1 vs 1 subunit
  list("TGFb agonist", "TGFb antagonist"),        # 2 vs 3 subunits
  list("TGFb inhibition receptor", "WNT antagonist"),  # 8 and 1 subunits
  list("NOT_A_COFACTOR", ""),                     # missing cofactor row -> all ones
  list("IL12AB", ""),                             # a *complex* name in the cofactor slot
  list("NODAL agonist", "")                       # multi-subunit, single field
)
for (case in ag_cases) {
  for (Kh_n in list(c(0.5, 1), c(0.5, 0.5), c(2, 2), c(1e-3, 1), c(1e3, 1))) {
    Kh <- Kh_n[[1]]; n <- Kh_n[[2]]
    for (idx in 1) {
      d <- data.frame(agonist = case[[1]], antagonist = case[[2]], stringsAsFactors = FALSE)
      a <- tryCatch(env$computeExpr_agonist(X, d, cofactor_input, idx, Kh, n),
                    error = function(e) structure(conditionMessage(e), class = "rerr"))
      an <- tryCatch(env$computeExpr_antagonist(X, d, cofactor_input, idx, Kh, n),
                     error = function(e) structure(conditionMessage(e), class = "rerr"))
      ## Both cofactor names are recorded: `computeExpr_agonist` reads
      ## `pairLRsig$agonist[index]` and the antagonist twin reads
      ## `pairLRsig$antagonist[index]`, so a record keyed on one name alone cannot be
      ## checked by a reader holding only that name.
      cat(sprintf("agonist\t%s\t%s\t%g\t%g\t%s\t%s\n", case[[1]], case[[2]], Kh, n,
                  if (inherits(a, "rerr")) "ERROR" else fmt(a),
                  if (inherits(an, "rerr")) "ERROR" else fmt(an)), file = q)
    }
  }
}

## ---------------------------------------------------------------- geometricMean
gm_cases <- list(c(1, 4), c(0, 1, 4), c(2, 8), c(1e-300, 1e300), c(1, 1, 1, 1, 1),
                 c(1, NA), c(NA, NA), c(1, Inf, 2), c(1, -Inf, 2), c(0.1, 0.2, 0.3))
for (i in seq_along(gm_cases)) {
  v <- gm_cases[[i]]
  o <- tryCatch(env$geometricMean(v), error = function(e) structure(conditionMessage(e), class = "rerr"))
  cat(sprintf("geomean\t%d\t%s\n", i, if (inherits(o, "rerr")) "ERROR" else format(o, digits = 17)),
      file = q)
}

## ---------------------------------------------------------------- group-mean variants
## `computeExprGroup_agonist` / `_antagonist` are exported but never called by
## CellChat itself; the fixture is here so they cannot be silently skipped.
grp <- factor(rep(paste0("g", 1:3), each = 2), levels = paste0("g", 1:3))
d <- data.frame(agonist = "TGFb agonist", antagonist = "TGFb antagonist",
                stringsAsFactors = FALSE)
cat(sprintf("group_agonist\t%s\t%s\t%s\n", d$agonist[1], d$antagonist[1],
            tryCatch(fmt(env$computeExprGroup_agonist(X, d, cofactor_input, grp, 1L,
                                                      0.5, env$triMean, 1)),
                     error = function(e) paste("ERROR", conditionMessage(e)))), file = q)
cat(sprintf("group_antagonist\t%s\t%s\t%s\n", d$agonist[1], d$antagonist[1],
            tryCatch(fmt(env$computeExprGroup_antagonist(X, d, cofactor_input, grp, 1L,
                                                         0.5, env$triMean, 1)),
                     error = function(e) paste("ERROR", conditionMessage(e)))), file = q)

close(q)
cat("wrote tests/fixtures/expr_golden.txt\n")
