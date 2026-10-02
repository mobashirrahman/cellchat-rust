# Golden corpus for r-core's pathway.rs: `computeCommunProbPathway` from the pinned commit.
#
# The fixture is built to stress the three R-isms documented in `crates/r-core/src/pathway.rs`:
#   * both `sum`s accumulate in LONG_DOUBLE, so the totals differ from an f64 loop;
#   * the two `apply` calls sum in *different* orders ((c, r) for `LR.sig`, L-R-then-(r, c)
#     for the pathway totals);
#   * `sort(..., decreasing = TRUE)` is R's shell sort, whose tie order is not obviously
#     ascending -- so the corpus includes pathways with *bit-identical* totals.
#
# Usage:  R_LIBS=.rlib R --vanilla -f tests/parity/gen_pathway_golden.R
suppressWarnings(suppressMessages(library(Matrix)))

CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
env <- new.env()
for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R")) {
  sys.source(file.path(CC, "R", f), envir = env, keep.source = FALSE)
}

## A small class carrying only the `net` slot the function reads.
## `netP` is a slot of its own upstream: the `object` branch writes
## `object@net$LRs` and `object@netP$pathways` / `$prob`.
setClass("MiniNet", representation(net = "list", netP = "list"))

q <- file("tests/fixtures/pathway_golden.txt", "wt")

## `dimnames` is built here rather than taken from the fixture, so every spec has it even
## when the `prob` array was constructed without names.
mk_net <- function(prob, pval, levels, lr) {
  prob <- array(prob, dim = dim(prob), dimnames = list(levels, levels, lr))
  pval <- array(pval, dim = dim(pval))
  new("MiniNet", net = list(
    prob = prob, pval = pval, prob.dim = dim(prob),
    dimnames = list(list(source = levels, target = levels),
                    list(source = levels, target = levels),
                    list(interaction_name = lr))
  ))
}

## --------------------------------------------------------------------------- the corpus
## `specs` describes networks to run through the pinned upstream. `prob`/`pval` are given
## as flat vectors in R's column-major order, so nothing is re-derived downstream.
fmt <- function(v) paste(sprintf("%.17g", as.numeric(v)), collapse = ",")

specs <- list()

## 1. Uniform random-ish, no zeros: the ordinary case.
{
  set.seed(20240421)
  k <- 4L; nlr <- 7L
  lev <- paste0("g", seq_len(k))
  lr <- paste0("LR", seq_len(nlr))
  prob <- array(runif(k * k * nlr), dim = c(k, k, nlr), dimnames = list(lev, lev, lr))
  pval <- array(runif(k * k * nlr), dim = c(k, k, nlr))
  specs[[length(specs) + 1]] <- list(name = "rand4", prob = prob, pval = pval, lev = lev, lr = lr,
                                    thresh = 0.5)
}

## 2. Several p-values above `thresh`, so `prob[pval > thresh] <- 0` bites and whole L-R
##    pairs drop out of `LR.sig`.
{
  k <- 3L; nlr <- 5L
  lev <- paste0("g", seq_len(k))
  lr <- paste0("LR", seq_len(nlr))
  prob <- array(seq_len(k * k * nlr) / 7, dim = c(k, k, nlr))
  pval <- array(rep(c(0.1, 0.9, 0.5, 0.05, 0.99), length.out = k * k * nlr), dim = c(k, k, nlr))
  specs[[length(specs) + 1]] <- list(name = "thresh3", prob = prob, pval = pval, lev = lev, lr = lr,
                                    thresh = 0.5)
}

## 3. An all-zero network: every total is 0, so `LR.sig` and `pathways.sig` are both empty
##    and the result is a 0-length pathway list. The degenerate case a port tends to
##    special-case wrongly.
{
  k <- 3L; nlr <- 4L
  lev <- paste0("g", seq_len(k))
  lr <- paste0("LR", seq_len(nlr))
  specs[[length(specs) + 1]] <- list(
    name = "allzero", prob = array(0, dim = c(k, k, nlr)), pval = array(1, dim = c(k, k, nlr)),
    lev = lev, lr = lr, thresh = 0.5
  )
}

## 4. **Tied pathway totals.** Two pathways get *exactly* the same L-R values, so their
##    totals are bit-identical and `sort(decreasing = TRUE)` has to break the tie. This is
##    the case that decides whether a stable descending sort reproduces R.
##    Pathways are attached by making LR1/LR2 and LR3/LR4 identical.
{
  k <- 3L
  lev <- paste0("g", seq_len(k))
  lr <- c("LR1", "LR2", "LR3", "LR4")
  base <- array(c(0.5, 0.25, 0.125, 1.5, 0.75, 0.375, 0.0625, 0.03125, 0.015625), dim = c(k, k, 2))
  prob <- array(0, dim = c(k, k, 4L))
  prob[, , 1:2] <- base   # pathway A: LR1, LR2
  prob[, , 3:4] <- base   # pathway B: LR3, LR4 -- identical, so tied
  pval <- array(0.1, dim = c(k, k, 4L))
  specs[[length(specs) + 1]] <- list(name = "tied", prob = prob, pval = pval, lev = lev, lr = lr,
                                    thresh = 0.5)
}

## 5. Three-way tie, and a tie *between a zero and a non-zero* total (only the non-zero
##    survives `!= 0`).
{
  k <- 2L
  lev <- paste0("g", seq_len(k))
  lr <- paste0("LR", 1:6)
  prob <- array(0, dim = c(k, k, 6L))
  prob[, , 1] <- 1.0; prob[, , 3] <- 1.0; prob[, , 5] <- 1.0   # pathways A, B, C all equal
  prob[, , 2] <- 2.0; prob[, , 4] <- 2.0                        # pathways D, E tied
  pval <- array(0.1, dim = c(k, k, 6L))
  specs[[length(specs) + 1]] <- list(name = "tied3", prob = prob, pval = pval, lev = lev, lr = lr,
                                    thresh = 0.5)
}

## 6. Unequal-magnitude values, so a LONG_DOUBLE sum and an f64 sum genuinely disagree and
##    the corpus can tell the two apart.
{
  k <- 3L; nlr <- 4L
  lev <- paste0("g", seq_len(k))
  lr <- paste0("LR", seq_len(nlr))
  prob <- array(0, dim = c(k, k, nlr))
  for (l in seq_len(nlr)) {
    for (r in seq_len(k)) {
      for (c in seq_len(k)) {
        ## 1e16 + 1 - 1e16 is exactly 0 in f64 but 1.0 with an 80-bit accumulator; the
        ## classic LONG_DOUBLE demonstration.
        prob[r, c, l] <- 1e16 + 1 - 1e16 * (r - 1) * (c - 1) * 0
        prob[r, c, l] <- 1e16 + 1
      }
      prob[r, , l] <- c(1e16, 1.0, -1e16 + 1e16)[seq_len(k)]
    }
  }
  pval <- array(0.1, dim = c(k, k, nlr))
  specs[[length(specs) + 1]] <- list(name = "bigmag", prob = prob, pval = pval, lev = lev, lr = lr,
                                    thresh = 0.5)
}

## 7. A threshold that zeroes everything (thresh = 0 with pval = 0) and one that zeroes
##    nothing (thresh = 0 with pval > 0). `pval > thresh` is strict, so pval == thresh
##    survives -- a boundary worth pinning.
{
  k <- 2L; nlr <- 3L
  lev <- paste0("g", seq_len(k))
  lr <- paste0("LR", seq_len(nlr))
  prob <- array(seq_len(k * k * nlr) / 3, dim = c(k, k, nlr))
  pval <- array(0, dim = c(k, k, nlr))
  specs[[length(specs) + 1]] <- list(name = "thr0", prob = prob, pval = pval, lev = lev, lr = lr,
                                    thresh = 0)
  pval2 <- array(0.25, dim = c(k, k, nlr))
  specs[[length(specs) + 1]] <- list(name = "thr025", prob = prob, pval = pval2, lev = lev, lr = lr,
                                     thresh = 0.25)
}

## 8. A single group (k = 1) and a single L-R: the minimum non-degenerate shape.
{
  specs[[length(specs) + 1]] <- list(
    name = "k1", prob = array(0.5, dim = c(1, 1, 1)), pval = array(0.1, dim = c(1, 1, 1)),
    lev = "g1", lr = "LR1", thresh = 0.5
  )
}

## ------------------------------------------------------------------- run and record
for (sp in specs) {
  o <- mk_net(sp$prob, sp$pval, sp$lev, sp$lr)
  ## `pairLR.use` is a data.frame with a `pathway_name` column. Give every *pair* of
  ## consecutive L-Rs its own pathway, so the corpus exercises several pathways per
  ## network and the ordering.
  nlr <- dim(sp$prob)[3]
  ## Two consecutive L-Rs per pathway. Written out rather than
  ## `rep(paste0("PW", ceiling(seq_len(nlr)/2)), each = 2, length.out = nlr)`: `rep` applies
  ## `each` to the *already expanded* vector, so that expression yields all "PW1" for an
  ## even `nlr` and collapses the network to a single pathway -- the same
  ## `rep(..., lengths=/each= + length.out)` trap as elsewhere in this repo.
  pw <- vapply(seq_len(nlr), function(i) paste0("PW", ((i - 1L) %/% 2L) + 1L), character(1))
  LRsig <- data.frame(
    interaction_name = sp$lr,
    pathway_name = pw, stringsAsFactors = FALSE
  )
  out <- tryCatch(
    suppressWarnings(suppressMessages(get("computeCommunProbPathway", envir = env)(
      object = o, net = NULL, pairLR.use = LRsig, thresh = sp$thresh))),
    error = function(e) structure(conditionMessage(e), class = "rerr")
  )
  cat(sprintf("spec\t%s\t%d\t%d\n", sp$name, dim(sp$prob)[1], dim(sp$prob)[3]), file = q)
  if (inherits(out, "rerr")) {
    cat(sprintf("pathway\t%s\tERROR\t%s\n", sp$name, out), file = q)
    next
  }
  ## `object` is not NULL, so upstream returns the S4 object: the interesting values are
  ## `out@net$LRs`, `out@netP$pathways` and `out@netP$prob`. The `object = NULL` branch
  ## returns a bare `list(pathways, prob)` and is recorded below as well, so both return
  ## shapes are pinned.
  cat(sprintf("pathway\t%s\tpathways=%s\n", sp$name,
              paste(out@netP$pathways, collapse = ",")), file = q)
  cat(sprintf("pathway\t%s\tlr_sig=%s\n", sp$name, paste(out@net$LRs, collapse = ",")), file = q)
  cat(sprintf("pathway\t%s\tdim=%s\n", sp$name, paste(dim(out@netP$prob), collapse = "x")), file = q)
  ## Flattened in R's column-major order, which is what the Rust side compares against.
  cat(sprintf("pathway\t%s\tprob=%s\n", sp$name, fmt(out@netP$prob)), file = q)
  cat(sprintf("pathway\t%s\tdimnames=%s\n", sp$name,
              paste(unlist(lapply(dimnames(out@netP$prob), paste, collapse = "|")), collapse = ",")),
      file = q)
  bare <- tryCatch(
    suppressWarnings(suppressMessages(get("computeCommunProbPathway", envir = env)(
      object = NULL, net = mk_net(sp$prob, sp$pval, sp$lev, sp$lr)@net, pairLR.use = LRsig,
      thresh = sp$thresh))),
    error = function(e) structure(conditionMessage(e), class = "rerr"))
  if (inherits(bare, "rerr")) {
    cat(sprintf("pathway\t%s\tbare=ERROR\t%s\n", sp$name, bare), file = q)
  } else {
    cat(sprintf("pathway\t%s\tbare=%s\n", sp$name, paste(bare$pathways, collapse = ",")), file = q)
    cat(sprintf("pathway\t%s\tbare_prob=%s\n", sp$name, fmt(bare$prob)), file = q)
  }
  ## The pre-`aperm` intermediate and both `apply` sums, so the two different summation
  ## orders are pinned separately rather than only through their difference.
  probz <- sp$prob
  probz[sp$pval > sp$thresh] <- 0
  cat(sprintf("pathway\t%s\tlr_sums=%s\n", sp$name, fmt(apply(probz, 3, sum))), file = q)
  pwp <- aperm(apply(probz, c(1, 2), by, factor(pw, levels = unique(pw)), sum), c(2, 3, 1))
  cat(sprintf("pathway\t%s\tpw_sums=%s\n", sp$name, fmt(apply(pwp, 3, sum))), file = q)
  ## The `k x k x nPathways` intermediate *before* `aperm(..., c(2, 3, 1))`, so the
  ## per-(source, target) pathway values are pinned independently of the permutation and
  ## of the descending-total reordering. Without it a mismatch in `netP$prob` can only be
  ## localised to "the sum or the permutation", not to which.
  ## NB `pwp` here is already `aperm(apply(...), c(2, 3, 1))`, i.e. the **`k x k x
  ## nPathways`** array -- the same axis order as `netP$prob`, only without the
  ## decreasing-total reordering. The pre-`aperm` array has the pathway *first*, and
  ## recording that one instead is an easy mistake: `apply` assembles
  ## `dim = c(dn[MARGIN], nResults)` and the results land on the last axis, so the pathway
  ## is not where the code's variable names suggest.
  cat(sprintf("pathway\t%s\tpwp_aperm=%s\n", sp$name, fmt(pwp)), file = q)
  cat(sprintf("pathway\t%s\tpwp_dim=%s\n", sp$name, paste(dim(pwp), collapse = "x")), file = q)
}
close(q)

## ------------------------------------------------------------------ dump the inputs
## The Rust test needs the *same* `prob`/`pval`, and this repo's rule is mechanical: if a
## fixture needs an RNG, the generator writes it out. `rand4` is `runif()` after
## `set.seed(20240421)`, and re-deriving R's `runif` in Rust would be the same mistake the
## wilcox fixture already made twice. `%a` hex, so the bits survive the round trip exactly.
qi <- file("tests/fixtures/pathway_inputs.tsv", "wt")
for (sp in specs) {
  k <- dim(sp$prob)[1]
  nlr <- dim(sp$prob)[3]
  lev <- sp$lev
  lr <- sp$lr
  cat(sprintf("spec\t%s\t%d\t%d\t%s\t%s\n", sp$name, k, nlr,
              paste(lev, collapse = ","), paste(lr, collapse = ",")), file = qi)
  pw <- vapply(seq_len(nlr), function(i) paste0("PW", ((i - 1L) %/% 2L) + 1L), character(1))
  cat(sprintf("pathways\t%s\t%s\n", sp$name, paste(unique(pw), collapse = ",")), file = qi)
  cat(sprintf("thresh\t%s\t%.17g\n", sp$name, sp$thresh), file = qi)
  cat("prob\n", file = qi)
  writeLines(paste(sprintf("%a", as.numeric(sp$prob)), collapse = " "), qi)
  cat("pval\n", file = qi)
  writeLines(paste(sprintf("%a", as.numeric(sp$pval)), collapse = " "), qi)
}
close(qi)
cat("wrote tests/fixtures/pathway_golden.txt and pathway_inputs.tsv\n", file = stderr())
