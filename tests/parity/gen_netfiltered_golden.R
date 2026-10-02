# Golden corpus for r-core's `aggregateNet` filtered branch
# (`net::aggregate_net_filtered`), plus the unfiltered branch for contrast.
#
# The point of the filtered branch is that the grouping key is `paste(source, target, sep =
# "|")` -- a *string* -- and `dplyr::group_by` orders by that string. So the corpus
# deliberately includes:
#   * a group name that is a prefix of another ("g1" vs "g10"), where string order and
#     factor-level order disagree;
#   * group levels that are *not* alphabetically ordered, so `cells.level` order and sorted
#     order disagree;
#   * a subset that leaves sources and targets with *different* level sets, so the result is
#     not square and the two `dimnames` differ;
#   * `remove.isolate = FALSE`, which restores the full level set;
#   * an empty cell, which is `NA` until the two `is.na` fills.
#
# Usage:  R_LIBS=.rlib R --vanilla -f tests/parity/gen_netfiltered_golden.R
suppressWarnings(suppressMessages({library(Matrix); library(collapse); library(dplyr)}))

CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
env <- new.env()
for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R")) {
  sys.source(file.path(CC, "R", f), envir = env, keep.source = FALSE)
}

## Enough of a `CellChat` for `subsetCommunication` and `aggregateNet`: a `net` slot and the
## `idents` whose levels define `cells.level`.
setClass("MiniAgg", representation(net = "list", idents = "factor", LR = "list",
                                   LRsig = "data.frame", options = "list"))

fmt <- function(v) paste(sprintf("%.17g", as.numeric(v)), collapse = ",")
q <- file("tests/fixtures/netfiltered_golden.txt", "wt")

mk <- function(prob, pval, levels, lr) {
  k <- dim(prob)[1]
  dimnames(prob) <- list(levels, levels, lr)
  dimnames(pval) <- list(levels, levels, lr)
  ## Two pathways, so the `signaling` filter has something to keep and something to drop.
  ## With a single pathway a `signaling` filter either keeps everything (so the case is not
  ## a filtered case at all -- `aggregateNet`'s `is.null(...) & ...` guard sends it down the
  ## unfiltered branch) or removes everything (so it errors).
  ## The full `LRsig` column set. `subsetCommunication`'s `signaling` branch `dplyr::select`s
  ## on `pathway_name` *and* on the annotation/evidence columns, so a two-column table fails
  ## with "Can't select columns that don't exist" before reaching the filter.
  LRsig <- data.frame(
    interaction_name = lr,
    interaction_name_2 = lr,
    pathway_name = rep(c("PWA", "PWB"), length.out = length(lr)),
    ligand = paste0("L", seq_along(lr)),
    receptor = paste0("R", seq_along(lr)),
    annotation = rep("Secreted Signaling", length(lr)),
    evidence = NA_character_,
    stringsAsFactors = FALSE
  )
  new("MiniAgg",
      net = list(prob = prob, pval = pval, prob.dim = dim(prob),
                 dimnames = list(list(source = levels, target = levels),
                                 list(source = levels, target = levels),
                                 list(interaction_name = lr))),
      idents = factor(levels, levels = levels),
      LR = list(LRsig = LRsig), LRsig = LRsig,
      options = list(datatype = "RNA", mode = "single"))
}

## --------------------------------------------------------------------------- the corpus
cases <- list()

## 1. The plain case: 3 levels, a `sources.use` filter that keeps everything, 2 thresholds.
{
  k <- 3L; nlr <- 6L
  lev <- c("g1", "g2", "g3")
  lr <- paste0("LR", seq_len(nlr))
  set.seed(20240431)
  prob <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  pval <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  cases[[length(cases) + 1]] <- list(name = "plain3", k = k, nlr = nlr, lev = lev, lr = lr,
                                    prob = prob, pval = pval, thresh = 0.5,
                                    sources = NULL, targets = NULL, signaling = NULL,
                                    pairLR = NULL, remove = TRUE)
}

## 2. **"g1" is a prefix of "g10".** `paste(source, target, sep = "|")` sorts `"g1|g2"`
##    before `"g10|g1"` as a string, but factor order puts `g10` before `g1` only if the
##    levels say so -- here they do not, so the two orders agree. The disagreement shows up
##    in `df.net2`'s row order, which is what `df.net2$prob <- df.net3$prob` and the
##    per-group summation order depend on.
{
  lev <- c("g1", "g10", "g2", "g3")
  k <- length(lev); nlr <- 5L
  lr <- paste0("LR", seq_len(nlr))
  set.seed(20240432)
  prob <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  pval <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  cases[[length(cases) + 1]] <- list(name = "prefix4", k = k, nlr = nlr, lev = lev, lr = lr,
                                    prob = prob, pval = pval, thresh = 0.5,
                                    sources = NULL, targets = NULL, signaling = NULL,
                                    pairLR = NULL, remove = TRUE)
}

## 3. **Level order is not alphabetical**, so `cells.level` order and `sort()` order
##    disagree: the `dimnames` must follow `cells.level`, not the sorted key.
{
  lev <- c("g3", "g1", "g10", "g2")
  k <- length(lev); nlr <- 4L
  lr <- paste0("LR", seq_len(nlr))
  set.seed(20240433)
  prob <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  pval <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  cases[[length(cases) + 1]] <- list(name = "unsorted4", k = k, nlr = nlr, lev = lev, lr = lr,
                                    prob = prob, pval = pval, thresh = 0.5,
                                    sources = NULL, targets = NULL, signaling = NULL,
                                    pairLR = NULL, remove = TRUE)
}

## 4. **Asymmetric subsets**: `sources.use` and `targets.use` keep different level sets, so
##    `remove.isolate = TRUE` gives a non-square result with different `dimnames` per axis.
{
  lev <- c("g1", "g2", "g3", "g4")
  k <- length(lev); nlr <- 5L
  lr <- paste0("LR", seq_len(nlr))
  set.seed(20240434)
  prob <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  pval <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  cases[[length(cases) + 1]] <- list(name = "asym4", k = k, nlr = nlr, lev = lev, lr = lr,
                                    prob = prob, pval = pval, thresh = 0.5,
                                    sources = c("g1", "g2"), targets = c("g2", "g3", "g4"),
                                    signaling = NULL, pairLR = NULL, remove = TRUE)
  cases[[length(cases) + 1]] <- list(name = "asym4_keep", k = k, nlr = nlr, lev = lev, lr = lr,
                                    prob = prob, pval = pval, thresh = 0.5,
                                    sources = c("g1", "g2"), targets = c("g2", "g3", "g4"),
                                    signaling = NULL, pairLR = NULL, remove = FALSE)
}

## 5. **A `signaling` filter**, which is a *different* code path inside
##    `subsetCommunication` (`grepl` on the pathway column) and can leave whole cells empty.
{
  lev <- c("g1", "g2", "g3")
  k <- 3L; nlr <- 4L
  lr <- paste0("LR", seq_len(nlr))
  set.seed(20240435)
  prob <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  pval <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  cases[[length(cases) + 1]] <- list(name = "signaling3", k = k, nlr = nlr, lev = lev, lr = lr,
                                    prob = prob, pval = pval, thresh = 0.5,
                                    sources = NULL, targets = NULL, signaling = "PWA",
                                    pairLR = NULL, remove = TRUE)
  cases[[length(cases) + 1]] <- list(name = "signaling3_b", k = k, nlr = nlr, lev = lev, lr = lr,
                                    prob = prob, pval = pval, thresh = 0.5,
                                    sources = NULL, targets = NULL, signaling = "PWB",
                                    pairLR = NULL, remove = TRUE)
}

## 6. **A `pairLR.use` subset**: one interaction only, so most cells are empty and every
##    `weight` that survives is a single interaction's probability.
{
  lev <- c("g1", "g2", "g3")
  k <- 3L; nlr <- 6L
  lr <- paste0("LR", seq_len(nlr))
  set.seed(20240436)
  prob <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  pval <- array(round(runif(k * k * nlr), 3), dim = c(k, k, nlr))
  cases[[length(cases) + 1]] <- list(name = "pairlr3", k = k, nlr = nlr, lev = lev, lr = lr,
                                    prob = prob, pval = pval, thresh = 0.5,
                                    sources = NULL, targets = NULL, signaling = NULL,
                                    ## `pairLR.use` must be a one-column data.frame named
                                    ## `interaction_name` or `pathway_name`; a bare character
                                    ## vector is rejected with "pairLR.use should be a data
                                    ## frame with a signle column named either ...".
                                    pairLR = data.frame(interaction_name = c("LR2", "LR5"),
                                                        stringsAsFactors = FALSE),
                                    remove = TRUE)
}

## 7. Thresholds at both ends, including one that keeps everything (`thresh = 1`).
{
  lev <- c("g1", "g2")
  k <- 2L; nlr <- 3L
  lr <- paste0("LR", seq_len(nlr))
  prob <- array(c(0.5, 0.25, 0.75, 0.125, 0.625, 0.375, 0.875, 0.0625, 0.3125, 0.4375, 0.5625,
                  0.1875), dim = c(k, k, nlr))
  pval <- array(c(0.01, 0.02, 0.03, 0.04, 0.05, 0.06, 0.07, 0.08, 0.09, 0.10, 0.50, 0.99),
                dim = c(k, k, nlr))
  for (thr in c(0, 0.05, 0.5, 1)) {
    cases[[length(cases) + 1]] <- list(
      name = paste0("thr", gsub("0\\.", "", format(thr))), k = k, nlr = nlr, lev = lev, lr = lr,
      prob = prob, pval = pval, thresh = thr, sources = NULL, targets = NULL,
      signaling = NULL, pairLR = NULL, remove = TRUE
    )
  }
}

## ------------------------------------------------------------------- run and record
for (cs in cases) {
  o <- mk(cs$prob, cs$pval, cs$lev, cs$lr)
  LRsig <- o@LR$LRsig
  run <- function(f) suppressWarnings(suppressMessages(
    f(o, sources.use = cs$sources, targets.use = cs$targets, signaling = cs$signaling,
      pairLR.use = cs$pairLR, remove.isolate = cs$remove, thresh = cs$thresh,
      return.object = TRUE)))
  out <- tryCatch(run(get("aggregateNet", envir = env)), error = function(e) structure(
    conditionMessage(e), class = "rerr"))
  cat(sprintf("case\t%s\t%d\t%d\n", cs$name, cs$k, cs$nlr), file = q)
  if (inherits(out, "rerr")) {
    cat(sprintf("agg\t%s\tERROR\t%s\n", cs$name, out), file = q)
    next
  }
  cat(sprintf("agg\t%s\tdim=%s\n", cs$name,
              paste(dim(out@net$count), collapse = "x")), file = q)
  cat(sprintf("agg\t%s\tcount=%s\n", cs$name, fmt(out@net$count)), file = q)
  cat(sprintf("agg\t%s\tweight=%s\n", cs$name, fmt(out@net$weight)), file = q)
  cat(sprintf("agg\t%s\tsource_levels=%s\n", cs$name,
              paste(rownames(out@net$count), collapse = ",")), file = q)
  cat(sprintf("agg\t%s\ttarget_levels=%s\n", cs$name,
              paste(colnames(out@net$count), collapse = ",")), file = q)
  cat(sprintf("agg\t%s\trn_source=%s\n", cs$name,
              paste(rownames(out@net$count), collapse = ",")), file = q)
  cat(sprintf("agg\t%s\trn_target=%s\n", cs$name,
              paste(colnames(out@net$count), collapse = ",")), file = q)
  ## The intermediate `df.net2` order, which is what the string-keyed `group_by` produces.
  ## `subsetCommunication` itself raises when a filter leaves nothing; the aggregateNet case
  ## then has no `df.net2` to record. Guarded so one such case does not abort the corpus.
  df <- tryCatch(suppressWarnings(suppressMessages(get("subsetCommunication", envir = env)(
    object = o, slot.name = "net", sources.use = cs$sources, targets.use = cs$targets,
    signaling = cs$signaling, pairLR.use = cs$pairLR, thresh = cs$thresh))),
    error = function(e) NULL)
  if (!is.null(df)) {
    ## The column set and classes are recorded *before* `source_target` is grafted on. That column
    ## belongs to `aggregateNet`'s filtered branch, not to `subsetCommunication`'s return value, so
    ## recording it here would pin a column the function does not return.
    cat(sprintf("agg\t%s\tcol_order=%s\n", cs$name,
                paste(colnames(df), collapse = ",")), file = q)
    cat(sprintf("agg\t%s\tclasses=%s\n", cs$name,
                paste(vapply(df, function(x) class(x)[1], ""), collapse = ",")), file = q)
    if (is.factor(df$interaction_name)) {
      cat(sprintf("agg\t%s\tinteraction_levels=%s\n", cs$name,
                  paste(levels(df$interaction_name), collapse = ",")), file = q)
    }
    if (is.factor(df$source)) {
      cat(sprintf("agg\t%s\tdf_source_levels=%s\n", cs$name,
                  paste(levels(df$source), collapse = ",")), file = q)
    }
    df$source_target <- paste(df$source, df$target, sep = "|")
    d2 <- df %>% group_by(source_target) %>%
      summarize(count = n(), prob = sum(prob), .groups = "drop")
    cat(sprintf("agg\t%s\tdf_order=%s\n", cs$name,
                paste(d2$source_target, collapse = ",")), file = q)
    cat(sprintf("agg\t%s\tdf_count=%s\n", cs$name, fmt(d2$count)), file = q)
    cat(sprintf("agg\t%s\tdf_prob=%s\n", cs$name, fmt(d2$prob)), file = q)

  }
}
close(q)

## ------------------------------------------------------------------ dump the inputs
## The filtered branch's *output* is empty or all zeros whatever the inputs, so a test that
## only checked the result would pass on any implementation at all. The inputs are dumped so
## the Rust side can also pin the part of the branch that is correct -- the byte-wise key
## order of `df.net2` -- and so it can prove the zeros are the *upstream* zeros.
qi <- file("tests/fixtures/netfiltered_inputs.tsv", "wt")
## All the metadata for a case is on **one** line, `|`-separated, so a reader is a single
## `strsplit` rather than a scan for tagged records. (An earlier version wrote one tagged
## record per field and needed a field-index parser on both the R and Rust sides; two
## hand-rolled parsers, two chances to be subtly wrong, for no benefit.)
enc <- function(v) if (is.null(v)) "none" else if (is.data.frame(v)) "df" else "vec"
encv <- function(v) {
  if (is.null(v)) "" else if (is.data.frame(v)) paste(v[[1]], collapse = ",") else paste(v, collapse = ",")
}
for (cs in cases) {
  meta <- c(
    cs$k, cs$nlr, paste(cs$lev, collapse = ","), paste(cs$lr, collapse = ","),
    sprintf("%.17g", cs$thresh), cs$remove,
    enc(cs$sources), encv(cs$sources), enc(cs$targets), encv(cs$targets),
    enc(cs$signaling), encv(cs$signaling), enc(cs$pairLR), encv(cs$pairLR)
  )
  cat(sprintf("case\t%s\t%s\n", cs$name, paste(meta, collapse = "|")), file = qi)
  cat("prob\n", file = qi)
  writeLines(paste(sprintf("%a", as.numeric(cs$prob)), collapse = " "), qi)
  cat("pval\n", file = qi)
  writeLines(paste(sprintf("%a", as.numeric(cs$pval)), collapse = " "), qi)
}
close(qi)
cat("wrote tests/fixtures/netfiltered_golden.txt and netfiltered_inputs.tsv\n", file = stderr())
