#!/usr/bin/env Rscript
# Large randomised differential corpus for the RNG port.
#
#   Rscript tests/parity/gen_rng_stress.R > tests/fixtures/rng_stress.txt
#
# The point is breadth: `rng_parity.rs` proves a handful of documented values, this
# proves there is no shape of (seed, n) where the streams diverge. The lock target is
# "absolutely the same results", so a handful of seeds is not evidence.
set.seed(20240304)   # fixed so the corpus is reproducible
out <- stdout()

# (a) permutations across many population sizes, including every power of two and
#     the sizes where bits_needed() changes behaviour.
ns <- c(1:40, 2^(1:17), 2^(1:17) + 1L, 2^(1:17) - 1L, 1000L, 21557L, 65536L, 100000L)
for (n in unique(ns)) {
  set.seed(1000L + n)
  v <- sample.int(n, n)
  cat(sprintf("perm\t%d\t%s\n", n, paste(v, collapse = ",")), file = out)
}

# (b) a pseudorandom sweep of (seed, n) pairs
for (i in 1:400) {
  sd <- sample.int(2^31 - 1L, 1L)
  n  <- sample.int(20000L, 1L) + 1L
  set.seed(sd)
  v <- sample.int(n, n)
  cat(sprintf("perm\t%d\t%d\t%s\n", sd, n, paste(v, collapse = ",")), file = out)
}

# (c) unif_rand streams for the same seed sweep (first 16 draws only, to keep the
#     file small; rng_parity.rs covers long streams for named seeds)
for (i in 1:400) {
  sd <- sample.int(2^31 - 1L, 1L)
  set.seed(sd)
  cat(sprintf("unif\t%d\t%s\n", sd,
              paste(format(runif(16), digits = 17), collapse = ",")), file = out)
}
