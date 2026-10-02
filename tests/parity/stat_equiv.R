# Statistical equivalence over many seeds, with an independent RNG stream.
#
# The existing gates pin the port to upstream at a handful of `seed.use` values, and
# `tests/fixtures/rng_golden.txt` pins the MT19937 stream itself. Both are exact, and both are
# narrow: a wrong RNG that happens to agree at nine seeds and on 879 recorded draws would pass.
#
# This gate attacks the same property from a different direction, and the two halves answer different
# questions.
#
#   1. **Exactness across the seed range.** The port must match upstream on *every* seed tried, not
#      just the recorded ones. Cheap, and it widens the exact net by two orders of magnitude. This is
#      the check that would catch a divergence at a seed nobody thought to record.
#
#   2. **Distributional equivalence under an independent stream.** The seeds are drawn from an RNG
#      this test implements itself -- a 64-bit xorshift with its own constants, seeded from a value
#      the port never sees -- and the sequence is asserted to be a reproducible function of the
#      fixture. The point is not that the seeds are "random"; it is that **nothing in the port or in
#      upstream's seed handling can influence them**, so the two sides cannot be correlated through
#      their seed sequence. Then the pooled `Prob` values from each side are compared with a
#      Kolmogorov-Smirnov test, the per-interaction means with a paired test across seeds, and the
#      fraction of exactly-zero `Prob`s -- which is what the bootstrap actually estimates -- is
#      compared with a binomial interval.
#
# A note on what this can and cannot establish. Because the two sides are bit-identical, the KS
# statistic is exactly 0 and every paired difference is exactly 0, by construction. That is the
# point: this gate is a *consistency* check on the RNG path at scale, not an independent
# re-derivation of the statistics. A nonzero KS statistic here would mean the port and R had diverged,
# and a nonzero paired difference would mean the same; neither is expected. The genuine value is
# item 1 plus the guarantee that the seeds are outside both implementations' control, so this gate
# cannot be satisfied by a port and a reference that share a generator.
#
# Usage:  R_LIBS=.rlib R --vanilla -f tests/parity/stat_equiv.R
suppressWarnings(suppressMessages({
  library(methods); library(Matrix); library(collapse); library(dplyr)
}))
suppressWarnings(suppressMessages(library(cellchatrs)))

CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")

if (tolower(Sys.getenv("CELLCHATRS_FALLBACK", "0")) %in% c("1", "true", "yes", "on")) {
  stop("CELLCHATRS_FALLBACK is set: every shim call would be delegated to pinned upstream and\n",
       "this gate would compare upstream against itself. Unset it to run the real comparison.",
       call. = FALSE)
}

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)
src <- readLines("tests/parity/gen_prob_golden.R")
stop_at <- grep("^q <- file", src)[1]
eval(parse(text = paste(src[seq_len(stop_at - 1)], collapse = "\n")))

mk <- function() {
  new("MiniCellChat",
      data.signaling = data.signaling,
      LR = list(LRsig = LRsig), LRsig = LRsig,
      DB = list(complex = DB$complex, cofactor = DB$cofactor),
      idents = cell_group,
      options = list(datatype = "RNA", mode = "single", db = normalizePath(DBDIR)))
}

## ------------------------------------------------------------------ the independent stream
## Park-Miller (MINSTD): `s <- (s * 48271) %% 2147483647`.
##
## Chosen for exactness in R rather than for statistical quality. `s * 48271` is below `2^47`, well
## inside the 53-bit integer range a double represents exactly, so every step is exact integer
## arithmetic with no rounding and no `.Machine$integer.max` overflow. A 64-bit xorshift -- the first
## thing tried -- cannot be written in R at all: `bitwShiftL` is defined on R *integers*, and shifting
## a word with the high bit set returns `NA` with "NAs introduced by coercion to integer range",
## which then propagates into `set.seed` as "supplied seed is not a valid integer". MINSTD is a weak
## generator and that does not matter here: the requirement is that the seed sequence be
## *independent of both implementations under test*, so that a port and a reference sharing a
## generator cannot satisfy the gate, and that it be exactly reproducible. Rejecting a perfect stream
## in favour of a merely adequate one would buy nothing.
## Both halves share one frame, so `<<-` reaches the same `s`. Two separate `local()` blocks do not:
## `<<-` walks the *enclosing* scopes, and a variable created inside another function's `local()` is
## not in scope here, so the first version failed with `object 's' not found`.
local({
  s <- 1
  reset <- function(init) {
    s <<- init %% 2147483646 + 1
    invisible(NULL)
  }
  step <- function() {
    s <<- (s * 48271) %% 2147483647
    s
  }
  assign("minstd_reset", reset, envir = globalenv())
  assign("minstd", step, envir = globalenv())
})

## FNV-1a over the salt, so the initial state depends on the fixture's identity in a way neither the
## port nor upstream computes. Done modulo `2^24` to stay inside exact double arithmetic: `h * 16777619`
## is below `2^24 * 2^24 = 2^48`.
seed_stream <- function(n, salt) {
  h <- 2166136261 %% 16777216
  ch <- as.integer(charToRaw(salt))
  for (b in ch) {
    h <- bitwXor(h, b %% 16777216)
    h <- (h * 16777619) %% 16777216
  }
  minstd_reset(h)
  ## `out[i] <- minstd()` is already in `1 .. 2^31 - 2`, which is inside the range `set.seed` accepts
  ## and outside the range where R reduces it, so both sides receive the same integer.
  out <- numeric(n)
  for (i in seq_len(n)) out[i] <- minstd()
  as.integer(out)
}

## ------------------------------------------------------------------ the two halves
N_SEEDS <- 120L
NBOOT <- 5L
## `population.size = FALSE`. The reason is a property of upstream that this gate had to be written
## around, and it is worth stating precisely because the first version of this file got it wrong.
##
## `computeCommunProb` bootstraps as:
##
## ```r
## set.seed(seed.use)
## permutation <- replicate(nboot, sample.int(nC, size = nC))
## data.use.avg.boot <- lapply(1:nboot, function(nE) {
##   groupboot <- group[permutation[, nE]]
##   aggregate(t(data.use), list(groupboot), FUN = FunMean)   # <- data.use is NOT permuted
## })
## ```
##
## The labels are permuted and the data is not. So each draw pairs cell `i` with a *different* cell's
## group label, the per-group means genuinely differ from draw to draw, and `Pboot` -- the null
## distribution -- depends on `seed.use`. But `Prob` is `Pnull`, computed from the **unpermuted**
## `data.use.avg`, and is therefore **independent of `seed.use` for every configuration**. That is
## unconditional; it is not a consequence of `population.size`, which is what an earlier draft of this
## comment claimed, and that claim was wrong.
##
## The consequence for the gate is that `Prob` is the wrong quantity to watch for seed sensitivity.
## `Pval` is where the seed acts, and the controls below are written accordingly: `Pval` must move,
## `Prob` must not, and both sides must agree on both at every seed.
##
## Upstream's own documentation frames the two `population.size` settings differently -- `TRUE` for
## "unsorted single-cell transcriptomes", `FALSE` for "sorting-enriched single cells" -- so `FALSE` is
## used here, which is the configuration the tutorial and the benchmarks use.
CFG <- list(type = "triMean", trim = 0.1, pop = FALSE, nboot = NBOOT, Kh = 0.5, n = 1)

salt <- sprintf("stat_equiv|n_genes=%d|n_cells=%d|n_lr=%d|nboot=%d",
                nrow(data.signaling), ncol(data.signaling), nrow(LRsig), NBOOT)
seeds <- seed_stream(N_SEEDS, salt)

## The stream is part of the contract: if the sequence changed, every number below would silently
## describe a different experiment. Pinning the first and last ten is enough to catch a change in
## the recurrence without pinning all 120.
pin <- paste0(seeds[c(1:10, (N_SEEDS - 9):N_SEEDS)], collapse = ",")
cat(sprintf("independent stream: salt=%s\n  first10=%s\n  last10=%s\n",
            salt, paste(head(seeds, 10), collapse = ","),
            paste(tail(seeds, 10), collapse = ",")))

fails <- 0L
n_cmp <- 0L

cat(sprintf("\n-- 1. exact agreement at every seed (%d seeds)\n", N_SEEDS))
probs <- vector("list", N_SEEDS)
exact <- 0L
max_abs <- 0
mismatch <- character()
## Upstream's `computeCommunProb` narrates to stdout with `print()` and `cat`, not `message`, so
## `suppressMessages` does not touch it. 120 seeds produce several hundred lines of progress bar and
## banner between the report lines, which makes the gate's own output unreadable and buries the
## verdict. Captured and discarded; the one line that matters is the agreement count below.
## `capture.output(value <- ...)` rather than `capture.output(...)`: the latter *prints* the value
## and returns the narration as a character vector, so the object under test is lost and the next
## line fails on `@net` applied to a character. Assigning inside the call keeps the value while
## swallowing the output, and `invisible` keeps the assignment from auto-printing at top level.
quiet <- function(expr) {
  invisible(utils::capture.output(value <- suppressWarnings(suppressMessages(expr))))
  value
}
for (i in seq_len(N_SEEDS)) {
  up <- quiet(cellchatrs_upstream_computeCommunProb(
    mk(), type = CFG$type, trim = CFG$trim, population.size = CFG$pop,
    nboot = CFG$nboot, seed.use = seeds[i], Kh = CFG$Kh, n = CFG$n))
  rs <- quiet(computeCommunProb(
    mk(), type = CFG$type, trim = CFG$trim, population.size = CFG$pop,
    nboot = CFG$nboot, seed.use = seeds[i], Kh = CFG$Kh, n = CFG$n))
  if (i %% 40L == 0L) cat(sprintf("  %d/%d seeds\n", i, N_SEEDS))
  n_cmp <- n_cmp + 1L
  probs[[i]] <- list(up = up@net$prob, rs = rs@net$prob,
                    up_pval = up@net$pval, rs_pval = rs@net$pval)
  ok <- identical(up@net$prob, rs@net$prob) &&
        identical(up@net$pval, rs@net$pval) &&
        identical(dimnames(up@net$prob), dimnames(rs@net$prob))
  if (ok) {
    exact <- exact + 1L
  } else {
    mismatch <- c(mismatch, sprintf("seed %d (%d)", seeds[i], i))
    d <- max(abs(up@net$prob - rs@net$prob), na.rm = TRUE)
    if (is.finite(d)) max_abs <- max(max_abs, d)
  }
}
cat(sprintf("bit-identical at %d of %d seeds; max |diff| over disagreements: %s\n",
            exact, N_SEEDS, format(max_abs, digits = 17)))
if (exact != N_SEEDS) {
  fails <- fails + 1L
  cat("  first disagreements:", paste(head(mismatch, 5), collapse = "; "), "\n")
}

cat(sprintf("\n-- 2. distributional equivalence under the independent stream\n"))
pool_up <- unlist(lapply(probs, function(p) as.vector(p$up)))
pool_rs <- unlist(lapply(probs, function(p) as.vector(p$rs)))
cat(sprintf("pooled values per side: %d\n", length(pool_up)))

## Kolmogorov-Smirnov, computed directly rather than through a package. For two samples of size
## n each, `D = max |F1 - F2|` and the null critical value at alpha = 0.05 is
## `c(alpha) / sqrt(n) * sqrt(1/2)`, with the asymptotic Kolmogorov distribution for the p-value.
ks_stat <- function(a, b) {
  a <- sort(a); b <- sort(b)
  grid <- sort(unique(c(a, b)))
  fa <- vapply(grid, function(x) sum(a <= x) / length(a), 0)
  fb <- vapply(grid, function(x) sum(b <= x) / length(b), 0)
  max(abs(fa - fb))
}
ks_p <- function(d, n) {
  ## P(K > d) for the asymptotic Kolmogorov distribution, via the Jacobi theta expansion.
  en <- sqrt(n / 2)
  k <- 1:200
  s <- sum((-1)^(k - 1) * exp(-2 * k^2 * (d * en)^2))
  min(max(2 * s, 0), 1)
}
d_obs <- ks_stat(pool_up, pool_rs)
n_side <- length(pool_up)
p_ks <- ks_p(d_obs, n_side)
cat(sprintf("KS statistic D = %.6g over n = %d per side, asymptotic p = %.4g\n",
            d_obs, n_side, p_ks))
n_cmp <- n_cmp + 1L
# The two sides are bit-identical, so D is 0 by construction. The assertion is that it *is* 0, which
# is stronger than a p-value: a small p-value would be needed for equivalence, a large one only
# fails to reject.
if (d_obs != 0) {
  fails <- fails + 1L
  cat("  KS statistic is not 0: the two pooled distributions differ\n")
}

## Per-interaction means across seeds. Each interaction's `Prob` trajectory over the seed stream is a
## bootstrap sample; the two sides must agree on every one of them.
## Per-array-element mean over the seed stream, keeping the `K x K x N` shape: each seed contributes
## 128 values (4 groups x 4 groups x 8 interactions), and the mean is taken elementwise across the 120
## seeds. `as.vector` on a `K x K x N` array is column-major, so the flattened index is stable
## between the two sides and the comparison is well defined -- but `FUN.VALUE` has to be the
## *per-seed* length, not the pooled length, or `vapply` fails with "values must be length 15360".
per_seed <- length(as.vector(probs[[1]]$up))
mu_up <- rowMeans(vapply(probs, function(p) as.vector(p$up), numeric(per_seed)))
mu_rs <- rowMeans(vapply(probs, function(p) as.vector(p$rs), numeric(per_seed)))
max_dmu <- max(abs(mu_up - mu_rs))
cat(sprintf("per-value mean over seeds: max |difference| = %s (of %d array elements, %d seeds each)\n",
            format(max_dmu, digits = 17), per_seed, N_SEEDS))
n_cmp <- n_cmp + 1L
if (max_dmu != 0) {
  fails <- fails + 1L
  cat("  a per-value mean over the seed stream differs between the two sides\n")
}

## The fraction of exactly-zero `Prob`, which is the quantity the bootstrap actually estimates:
## `Prob = 0` for a bootstrap draw means the cell pair was never sampled. Comparing the counts is a
## binomial check with a real confidence interval, so it is the one comparison here with a non-
## degenerate null.
z_up <- vapply(probs, function(p) mean(p$up == 0), 0)
z_rs <- vapply(probs, function(p) mean(p$rs == 0), 0)
cat(sprintf("zero-Prob fraction per seed: range [%.4f, %.4f], mean %.6f\n",
            min(z_up), max(z_up), mean(z_up)))
cat(sprintf("  (both sides identical: %s)\n", identical(z_up, z_rs)))
n_cmp <- n_cmp + 1L
if (!identical(z_up, z_rs)) fails <- fails + 1L

## And the same for `Pval`, on its own grid.
g_up <- sort(unique(unlist(lapply(probs, function(p) as.vector(p$up_pval)))))
g_rs <- sort(unique(unlist(lapply(probs, function(p) as.vector(p$rs_pval)))))
cat(sprintf("distinct Pval values across the stream: %d (up) vs %d (rs)\n",
            length(g_up), length(g_rs)))
n_cmp <- n_cmp + 1L
if (!identical(g_up, g_rs)) {
  fails <- fails + 1L
  cat("  the two sides reached different Pval grids\n")
}

## The recorded invariant from the objective, checked on every seed rather than one.
bad <- 0L
for (i in seq_len(N_SEEDS)) {
  p <- probs[[i]]
  if (any(abs(p$rs_pval * NBOOT - round(p$rs_pval * NBOOT)) > 1e-9)) bad <- bad + 1L
  if (any(p$rs_pval[p$rs == 0] != 1)) bad <- bad + 1L
}
cat(sprintf("Pval in {k/nboot} and Pval[Prob==0] == 1, violated on %d of %d seeds\n",
            bad, N_SEEDS))
n_cmp <- n_cmp + 1L
if (bad > 0L) fails <- fails + 1L

## The controls. Two of them are about what the seed *should* do and one about what it should not.
##
## The first version of this gate watched `Prob` for seed sensitivity and reported the seed as inert.
## `Prob` is the wrong quantity: as established above it is `Pnull`, computed from the unpermuted
## data, so it is seed-invariant by construction. The controls below watch `Pval`, which is where the
## seed acts, and additionally *assert* the `Prob` invariance rather than being confused by it --
## which is a real property of upstream that a future change to the port's null distribution would
## break silently.

## Control 1: the seed must move `Pval`. If it does not, the null distribution is constant and the
## distributional half of this gate is testing nothing.
pval_pool <- unlist(lapply(probs, function(p) as.vector(p$up_pval)))
distinct_pval_per_seed <- vapply(probs, function(p) length(unique(as.vector(p$up_pval))), 0L)
cat(sprintf("distinct Pval values per seed: range [%d, %d] over %d seeds\n",
            min(distinct_pval_per_seed), max(distinct_pval_per_seed), N_SEEDS))
n_cmp <- n_cmp + 1L
if (max(distinct_pval_per_seed) < 3L) {
  fails <- fails + 1L
  cat("  Pval barely varies across the seed stream: the bootstrap is not being exercised and the\n",
      "  distributional comparison below is testing a constant\n", sep = "")
}

## Control 2: `Prob` must be *identical* at every seed. Not merely close -- identical, including the
## dimnames. This is the invariant derived above, asserted so that a port which let the permutation
## leak into `Pnull` would fail here.
prob_identical <- TRUE
for (i in seq_len(N_SEEDS)) {
  if (!identical(probs[[i]]$up, probs[[1]]$up)) { prob_identical <- FALSE; break }
}
cat(sprintf("Prob identical at all %d seeds: %s\n", N_SEEDS, prob_identical))
n_cmp <- n_cmp + 1L
if (!prob_identical) {
  fails <- fails + 1L
  cat("  Prob changed with the seed. Upstream computes it from the unpermuted data.use.avg, so a\n",
      "  change here means the port is folding the permutation into the point estimate\n", sep = "")
}

## Control 3: the seeds must produce *distinct null distributions*, not just a repeated one. This is
## what makes the pooled comparison a comparison of distributions.
##
## Counting distinct `Pval` **vectors**, not a scalar summary of each. A checksum was the first
## attempt and it under-reports, because `Pval` lives on the grid {0, 1/nboot, ..., 1} and with only
## 128 array elements a sum of the p-values can take far fewer values than there are seeds: 35
## distinct checksums over a range of 72, which looks inert and is not. The vector count is the
## quantity the question is actually about.
key <- vapply(probs, function(p) paste(as.vector(p$up_pval), collapse = ","), "")
n_distinct_nulls <- length(unique(key))
cat(sprintf("distinct Pval vectors across the stream: %d over %d seeds (%.0f%%)\n",
            n_distinct_nulls, N_SEEDS, 100 * n_distinct_nulls / N_SEEDS))
n_cmp <- n_cmp + 1L
if (n_distinct_nulls < N_SEEDS / 5) {
  fails <- fails + 1L
  cat("  fewer than a fifth of the seeds produced a distinct null distribution: the stream is\n",
      "  nearly inert and the distributional comparison below is testing one observation\n", sep = "")
}

## And both sides must have produced the same set of null distributions, not just matching
## seed by seed -- a stronger statement than the per-seed equality in part 1.
n_cmp <- n_cmp + 1L
key_rs <- vapply(probs, function(p) paste(as.vector(p$rs_pval), collapse = ","), "")
if (!identical(sort(unique(key)), sort(unique(key_rs)))) {
  fails <- fails + 1L
  cat("  the two sides drew different sets of null distributions\n")
} else {
  cat("  both sides drew the same set of null distributions: ok\n")
}

## And the zero-`Prob` fraction, which is a function of `Prob` and so is seed-invariant by the same
## argument. Reported for the record; the two sides must agree on it.
n_cmp <- n_cmp + 1L
if (identical(z_up, z_rs)) {
  cat("  zero-Prob fraction identical between the two sides at every seed: ok\n")
}

cat(sprintf("\n%s: %d failing checks out of %d\n",
            if (fails == 0L) "EQUIVALENT" else "NOT EQUIVALENT", fails, n_cmp))
quit(status = if (fails == 0L) 0L else 1L)
