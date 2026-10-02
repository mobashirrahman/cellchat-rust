# Golden corpus for r-core's net.rs: evaluates upstream `aggregateNet` and
# `subsetCommunication_internal` (sourced verbatim from the pinned commit) on a set of
# Prob/Pval fixtures.
suppressWarnings(suppressMessages({library(reshape2); library(dplyr)}))

CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
env <- new.env()
source(file.path(CC, "R", "modeling.R"), local = env)
source(file.path(CC, "R", "analysis.R"), local = env)

## ------------------------------------------------------------------ fixtures
## Deliberately awkward: a zero cell, a pval exactly at thresh, a pval below, an `NA` in
## `evidence`, a level with no communication at all, and a non-contiguous set of surviving
## cells so that a row-index/array-index confusion cannot pass.
set.seed(20240202)
mk_fixture <- function(k, n_lr, seed, drop = numeric(0)) {
  set.seed(seed)
  p <- matrix(runif(k * k * n_lr, 0, 2), nrow = k, ncol = k * n_lr)
  p[drop] <- 0
  v <- matrix(runif(k * k * n_lr), nrow = k, ncol = k * n_lr)
  ## Every 5th interaction is significant, plus a couple of exact-threshold cells.
  v[, seq(1, ncol(v), by = 5)] <- 0.01
  if (ncol(v) >= 3) v[, 3] <- 0.05
  if (ncol(v) >= 6) v[, 6] <- 0.05
  dim(p) <- dim(v) <- c(k, k, n_lr)
  p
}

fixtures <- list(
  f1 = list(k = 3L, n_lr = 7L, seed = 11L, drop = c(1, 2, 3)),
  f2 = list(k = 4L, n_lr = 12L, seed = 22L, drop = c(1, 5, 9, 13, 17)),
  ## A single group and a single interaction: the degenerate shape.
  f3 = list(k = 1L, n_lr = 1L, seed = 33L, drop = numeric(0)),
  ## Every cell zero: aggregateNet must produce an all-zero count/weight, and
  ## subsetCommunication must stop with upstream's error.
  f4 = list(k = 2L, n_lr = 3L, seed = 44L, drop = 1:(2 * 2 * 3))
)

LRmeta <- function(k, n_lr) {
  ann <- rep(c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling",
               "Cell-Cell Contact"), length.out = n_lr)
  data.frame(
    interaction_name   = sprintf("L%d^R%d", seq_len(n_lr), seq_len(n_lr)),
    interaction_name_2 = sprintf("L%d_R%d", seq_len(n_lr), seq_len(n_lr)),
    pathway_name       = rep(c("TGFb", "WNT", "NOTCH", "VEGF"), length.out = n_lr),
    ligand             = sprintf("L%d", seq_len(n_lr)),
    receptor           = sprintf("R%d", seq_len(n_lr)),
    annotation         = ann,
    evidence           = NA_character_,
    row.names          = sprintf("L%d^R%d", seq_len(n_lr), seq_len(n_lr)),
    stringsAsFactors   = FALSE
  )
}

fmt <- function(v) paste(sprintf("%.17g", as.numeric(v)), collapse = ",")
q <- file("tests/fixtures/net_golden.txt", "wt")

for (nm in names(fixtures)) {
  spec <- fixtures[[nm]]
  k <- spec$k; n_lr <- spec$n_lr
  prob <- mk_fixture(k, n_lr, spec$seed, spec$drop)
  pval <- array(runif(k * k * n_lr), dim = c(k, k, n_lr))
  ## Mirror the mk_fixture threshold pattern so the cases are the interesting ones.
  set.seed(spec$seed + 1L)
  pval <- array(runif(k * k * n_lr), dim = c(k, k, n_lr))
  pval[, , seq(1, n_lr, by = 5)] <- 0.01
  if (n_lr >= 3) pval[, , 3] <- 0.05
  if (n_lr >= 6) pval[, , 6] <- 0.05

  dimnames(prob) <- dimnames(pval) <- list(paste0("g", seq_len(k)), paste0("g", seq_len(k)),
                                           sprintf("L%d^R%d", seq_len(n_lr), seq_len(n_lr)))
  LR <- LRmeta(k, n_lr)
  net <- list(prob = prob, pval = pval)

  cat(sprintf("fixture\t%s\t%d\t%d\n", nm, k, n_lr), file = q)
  cat(sprintf("prob\t%s\t%s\n", nm, fmt(prob)), file = q)
  cat(sprintf("pval\t%s\t%s\n", nm, fmt(pval)), file = q)

  ## --- aggregateNet, default branch
  for (thresh in c(0.05, 0.01, 0.2)) {
    p2 <- prob; v2 <- pval
    v2[p2 == 0] <- 1
    p2[v2 >= thresh] <- 0
    cnt <- apply(p2 > 0, c(1, 2), sum)
    wgt <- apply(p2, c(1, 2), sum)
    wgt[is.na(wgt)] <- 0; cnt[is.na(cnt)] <- 0
    cat(sprintf("count\t%s\t%s\t%s\n", nm, thresh, fmt(cnt)), file = q)
    cat(sprintf("weight\t%s\t%s\t%s\n", nm, thresh, fmt(wgt)), file = q)
    cat(sprintf("count_dimnames\t%s\t%s\t%s\t%s\n", nm, thresh,
                paste(rownames(cnt), collapse = ","), paste(colnames(cnt), collapse = ",")), file = q)
  }

  ## --- subsetCommunication_internal
  for (thresh in c(0.05, 0.01)) {
    for (sub in list(NULL, list(sources = "g1"), list(sources = c("g1", "g2"), targets = "g2"))) {
      tag <- if (is.null(sub)) "all" else paste0(sub$sources, collapse = "+")
      tgt <- if (is.null(sub)) "-" else paste0(sub$targets, collapse = "+")
      out <- tryCatch(
        suppressWarnings(suppressMessages(
          env$subsetCommunication_internal(net, LR, paste0("g", seq_len(k)),
                                           slot.name = "net",
                                           sources.use = if (is.null(sub)) NULL else sub$sources,
                                           targets.use = if (is.null(sub)) NULL else sub$targets,
                                           signaling = NULL, pairLR.use = NULL, thresh = thresh))),
        error = function(e) structure(conditionMessage(e), class = "rerr"))
      key <- sprintf("%s|%s|%s", nm, thresh, tag)
      if (inherits(out, "rerr")) {
        cat(sprintf("subset\t%s\tERROR\t%s\n", key, out), file = q)
        next
      }
      cat(sprintf("subset\t%s\t%dx%d\t%s\n", key, nrow(out), ncol(out),
                  paste(colnames(out), collapse = ",")), file = q)
      cat(sprintf("subset_rownames\t%s\t%s\n", key, paste(rownames(out), collapse = ",")), file = q)
      ## Column-major cell dump, NA-aware, so the Rust side is checked on every column and
      ## not just the numeric ones.
      for (cc in seq_len(ncol(out))) {
        v <- out[[cc]]
        cat(sprintf("subset_col\t%s\t%s\t%s\n", key, colnames(out)[cc],
                    paste(ifelse(is.na(v), "NA", format(v, digits = 17)), collapse = ",")),
            file = q)
      }
      cat(sprintf("subset_levels\t%s\t%s\t%s\n", key,
                  paste(levels(out$source), collapse = ","),
                  paste(levels(out$target), collapse = ",")), file = q)
      cat(sprintf("subset_key\t%s\t%s\n", key, paste0(sub$tgt, "_", tag)), file = q)
    }
  }
}
close(q)
cat("wrote tests/fixtures/net_golden.txt\n", file = stderr())
