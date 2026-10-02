#!/usr/bin/env Rscript
# Golden corpus for the statistics port, produced from the *pinned* R + collapse.
#
#   Rscript tests/parity/gen_stats_golden.R > tests/fixtures/stats_golden.txt
#
# Two things are under test and neither may be assumed:
#
#  1. `triMean` goes through `collapse::fquantile(type = 7)`, which is NOT bit-identical
#     to R's `quantile(type = 7)` (measured: 7.4% of cases, up to 2 ulp). The corpus
#     therefore records `collapse`'s answer, not R's, and also records R's so the test
#     can prove the two really do differ -- otherwise a port that implemented the wrong
#     one would still pass.
#  2. R's `mean` accumulates in LONG_DOUBLE (64-bit mantissa). Whether a plain f64
#     accumulation is equivalent is an *empirical* question; the corpus records enough
#     to answer it rather than leaving it to argument.
#
# Requires: collapse 2.1.8 (the version this project validated against).

.libPaths(".rlib")
suppressMessages(library(collapse))
stopifnot(identical(RNGkind()[1], "Mersenne-Twister"))

## Vector store: the Rust test cannot call R, so the *inputs* have to be shipped
## alongside the expected outputs. Written as little-endian f64 in one flat blob, with
## a TSV index of (id, kind, offset, length, params...).
VBLOB  <- file("tests/fixtures/stats_vectors.bin", "wb")
VIDX   <- file("tests/fixtures/stats_vectors.idx", "wt")
vcount <- 0L
store_vec <- function(id, x) {
  vcount <<- vcount + 1L
  off <- seek(VBLOB, where = 0, origin = "end")
  writeBin(as.numeric(x), VBLOB, size = 8, endian = "little", useBytes = TRUE)
  len <- length(x)
  cat(sprintf("%s\t%s\t%s\t%s\n", id, as.integer(off), as.integer(len), vcount), file = VIDX)
  invisible(NULL)
}

triMean <- function(x, na.rm = TRUE)
  mean(collapse::fquantile(x, probs = c(0.25, 0.50, 0.50, 0.75), na.rm = na.rm))
geometricMean <- function(x, na.rm = TRUE) exp(mean(log(x), na.rm = na.rm))
thresholdedMean <- function(x, trim = 0.1, na.rm = TRUE) {
  percent <- Matrix::nnzero(x) / length(x)
  if (percent < trim) return(0) else return(mean(x, na.rm = na.rm))
}
f <- function(v) paste(format(v, digits = 17), collapse = ",")

set.seed(20240304)  # fixed: the corpus must be reproducible

# ---------------------------------------------------------------------------
# 1. triMean over many vector shapes.
# ---------------------------------------------------------------------------
gens <- list(
  uniform01   = function(n) runif(n),
  counts      = function(n) rpois(n, 3),
  rounded     = function(n) round(runif(n) * 5),
  binary      = function(n) as.numeric(rbinom(n, 1, 0.3)),
  normal      = function(n) rnorm(n),
  skewed      = function(n) rexp(n),
  heavy_zero  = function(n) { v <- rpois(n, 0.4); v },
  lognorm     = function(n) log1p(rgamma(n, 1.5)),
  allzero     = function(n) rep(0, n),
  constant    = function(n) rep(2.5, n),
  two_valued  = function(n) sample(c(0, 1), n, TRUE)
)

for (gname in names(gens)) {
  g <- gens[[gname]]
  for (n in c(1, 2, 3, 4, 5, 7, 8, 16, 17, 31, 32, 33, 63, 100, 255, 256, 257,
              862, 863, 1000, 2000, 5000, 10000)) {
    for (rep in 1:3) {
      x <- g(n)
      q <- as.numeric(collapse::fquantile(x, probs = c(.25, .5, .5, .75), na.rm = TRUE))
      cat(sprintf("trimean\t%s\t%d\t%d\t%s\t%s\n", gname, n, rep, f(q), format(triMean(x), digits = 17)))
      store_vec(sprintf("trimean:%s:%d:%d", gname, n, rep), x)
    }
  }
}

# ---------------------------------------------------------------------------
# 1b. Quantify the collapse-vs-R divergence at BOTH levels, so the plan's claim is
#     measured rather than remembered. Element level is what matters for a port that
#     reimplements fquantile; trimean level is what a user would observe.
# ---------------------------------------------------------------------------
elem_dis <- 0; elem_tot <- 0; elem_max <- 0
tri_dis  <- 0; tri_tot  <- 0; tri_max  <- 0
set.seed(7)
for (i in 1:20000) {
  n <- sample.int(5000, 1)
  x <- switch(sample(1:5, 1), round(runif(n), 2), sample(c(0, 1), n, TRUE),
              rpois(n, 3) / 3, rnorm(n), sample(0:9, n, TRUE))
  if (i %% 5 == 0) x[sample.int(n, 1)] <- NaN
  keep <- !is.na(x)
  ca <- as.numeric(collapse::fquantile(x, probs = c(.25,.5,.5,.75), na.rm = TRUE))
  ra <- as.numeric(stats::quantile(x[keep], probs = c(.25,.5,.5,.75), type = 7))
  elem_tot <- elem_tot + length(ca)
  d <- abs(ca - ra)
  elem_dis <- elem_dis + sum(d > 0); elem_max <- max(elem_max, d)
  a <- mean(ca); b <- mean(ra)
  tri_tot <- tri_tot + 1
  if (!identical(a, b)) { tri_dis <- tri_dis + 1; tri_max <- max(tri_max, abs(a - b)) }
}
cat(sprintf("#divergence\tvectors=%d\telem_disagree=%d\telem_frac=%.4f\telem_max=%.6e\ttri_disagree=%d\ttri_frac=%.4f\ttri_max=%.6e\n",
            tri_tot, elem_dis, elem_dis / elem_tot, elem_max, tri_dis, tri_dis / tri_tot, tri_max))

# ---------------------------------------------------------------------------
# 2. Long-vector triMean: exercises collapse's radix-order path and large n.
# ---------------------------------------------------------------------------
for (n in c(50000, 100001, 200000)) {
  x <- c(rpois(n - 5, 0.3), rep(0, 5))
  cat(sprintf("trimean_big\t%d\t%s\n", n, format(triMean(x), digits = 17)))
  store_vec(sprintf("trimean_big:%d", n), x)
}

# ---------------------------------------------------------------------------
# 3. geometricMean: log-space mean, including the zero -> 0 rule.
# ---------------------------------------------------------------------------
for (n in c(1, 2, 3, 4, 5, 8, 16, 64, 256, 1024)) {
  for (rep in 1:3) {
    x <- runif(n, 0.01, 10)
    cat(sprintf("geomean\t%d\t%d\t%s\n", n, rep, format(geometricMean(x), digits = 17)))
    cat(sprintf("geomean_zero\t%d\t%d\t%s\n", n, rep,
                format(geometricMean(c(x, 0)), digits = 17)))
    store_vec(sprintf("geomean:%d:%d", n, rep), x)
  }
}
set.seed(11)
for (i in 1:500) {
  n <- sample.int(600, 1)
  x <- runif(n, 0.001, 5)
  if (i %% 4 == 0) x[sample.int(n, 1)] <- 0
  if (i %% 7 == 0) x[sample.int(n, 1)] <- NaN
  cat(sprintf("geomean_rand\t%d\t%s\n", i, format(geometricMean(x), digits = 17)))
  store_vec(sprintf("geomean_rand:%d", i), x)
}

# ---------------------------------------------------------------------------
# 4. thresholdedMean across the trim boundary.
# ---------------------------------------------------------------------------
for (n in c(1, 10, 100, 1000)) {
  fracs <- c(0, 0.1, 0.5, 1)
  nnzs <- sort(unique(c(0, 1, n - 1, n, as.integer(floor(fracs * n)))))
  nnzs <- nnzs[nnzs >= 0 & nnzs <= n]
  for (nnz in nnzs) {
    x <- c(rep(1, nnz), rep(0, n - nnz))
    if (length(x) != n) next
    ## Avoid a degenerate all-ones vector (mean = 1 exactly) except when nnz is 0.
    if (nnz > 0 && nnz == n) x <- x + 0.25
    for (trim in c(0, 0.1, 0.25, 0.5, 1)) {
      cat(sprintf("threshmean\t%d\t%d\t%.2f\t%s\n", n, nnz, trim,
                  format(thresholdedMean(x, trim = trim), digits = 17)))
      store_vec(sprintf("threshmean:%d:%d:%.2f", n, nnz, trim), x)
    }
  }
}

# ---------------------------------------------------------------------------
# 5. NaN/NA handling of every summary (R9).
# ---------------------------------------------------------------------------
specials <- list(
  with_nan      = c(1, NaN, 2, 3, 4),
  all_nan       = c(NaN, NaN, NaN),
  nan_heavy     = c(rep(NaN, 9), 5),
  with_na       = c(1, NA, 2, 3, 4),
  single_nan    = c(NaN),
  with_inf      = c(1, Inf, 2),
  with_neinf    = c(1, -Inf, 2),
  zeros         = c(0, 0, 0, 0)
)
for (nm in names(specials)) {
  x <- specials[[nm]]
  cat(sprintf("special_trimean\t%s\t%s\n", nm,
              format(triMean(x), digits = 17)))
  cat(sprintf("special_geomean\t%s\t%s\n", nm,
              format(geometricMean(x), digits = 17)))
  cat(sprintf("special_median\t%s\t%s\n", nm,
              format(median(x, na.rm = TRUE), digits = 17)))
  store_vec(sprintf("special:%s", nm), x)
}

close(VBLOB); close(VIDX)
cat(sprintf("#vectors\t%d\n", vcount))

# ---------------------------------------------------------------------------
# 6. Fingerprints, so a corpus/consumer mismatch is loud.
# ---------------------------------------------------------------------------
cat(sprintf("#env\tR=%s\tcollapse=%s\n", getRversion(),
            as.character(packageVersion("collapse"))))
cat(sprintf("#r_mean_kind\tLONG_DOUBLE\tfirst=%.17e\n", mean(c(1e16, 1, -1e16))))
