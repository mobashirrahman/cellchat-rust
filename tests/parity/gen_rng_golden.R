#!/usr/bin/env Rscript
# Generate the RNG golden file consumed by src/rust/crates/r-core/tests/rng_parity.rs.
#
# Ground truth comes from R itself, so the constants in the Rust tests are never
# hand-transcribed. Re-run whenever the pinned R version changes:
#
#   Rscript tests/parity/gen_rng_golden.R > tests/fixtures/rng_golden.txt

# --- unif_rand stream -------------------------------------------------------
# runif() consumes exactly one unif_rand() per value for the default normal kind
# (Runif_kind 3 leaves unif_rand untouched: the Box-Muller variant used does not
# cache, so consecutive runif() calls map 1:1 onto unif_rand()).
for (s in c(1L, 2L, 42L, 1234567L, -1L, .Machine$integer.max)) {
  set.seed(s)
  cat(sprintf("unif_rand\tseed=%d\t%s\n", s,
              paste(format(runif(24), digits = 17), collapse = ",")))
}

# --- sample.int streams -----------------------------------------------------
# Shape used by computeCommunProb: sample.int(nC, nC), no replacement, nC >= 2.
for (s in c(1L, 2L, 7L, 99L)) {
  for (n in c(2L, 3L, 10L, 64L, 257L, 1000L, 21557L)) {
    set.seed(s)
    v <- sample.int(n, n)
    cat(sprintf("sample_int\tseed=%d\tn=%d\t%s\n", s, n,
                paste(v, collapse = ",")))
  }
}

# --- sample.int partial (k < n) --------------------------------------------
for (s in c(1L, 3L)) {
  for (cfg in list(c(10L, 3L), c(100L, 40L), c(1000L, 999L), c(5000L, 2500L))) {
    set.seed(s)
    v <- sample.int(cfg[1], cfg[2])
    cat(sprintf("sample_int_k\tseed=%d\tn=%d\tk=%d\t%s\n", s, cfg[1], cfg[2],
                paste(v, collapse = ",")))
  }
}

# --- sample.int WITH replacement (sample()'s default path) ------------------
for (s in 1:2) {
  set.seed(s)
  v <- sample.int(100, 20, replace = TRUE)
  cat(sprintf("sample_int_rep\tseed=%d\t%s\n", s, paste(v, collapse = ",")))
}

# --- R environment fingerprint ---------------------------------------------
cat(sprintf("#R\t%s\t%s\tkind=%s\tsample.kind=%s\n",
            R.version.string, getRversion(),
            RNGkind()[1], RNGkind()[3]))
