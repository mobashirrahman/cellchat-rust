## Bit-identity on the authors' own data, as a *check* rather than a measurement.
##
## `bench-runner/bench_real.R` does this at the default `nboot = 100` while timing both sides; this
## script does only the comparison, at a smaller `nboot`, so it is cheap enough to be a CI step.
## `nboot` does not affect parity -- the permutations are a deterministic function of
## `set.seed(seed.use)` and the aggregation is per-replicate -- and `NB=100` re-runs it at the
## setting the published numbers use.
##
## Exits non-zero on any mismatch, so `parity_report.py` can use it as a pin and a green run means
## something.
suppressWarnings(suppressMessages({
  library(collapse); library(Matrix); library(dplyr); library(CellChat)
}))
`%||%` <- function(a, b) if (is.null(a)) b else a

here <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
fixture_path <- file.path(here, "fixture.R")

## `fixture.R` loads the database **at source time**, reading `Sys.getenv("SPECIES")`. Sourcing it
## once above this loop -- the obvious arrangement -- silently gives every case after the first the
## *first* case's database, and the symptom is upstream itself crashing in `geometricMean` on a
## "missing subunit" that is not missing. This has now bitten twice, once here and once in
## `bench_real.R` (there via a mismatched `CELLCHATRS_DB`), so the species is set and the file
## re-sourced per case rather than once up front.
load_fixture_for <- function(species) {
  Sys.setenv(SPECIES = species)
  source(fixture_path, local = FALSE)
}

NB <- as.integer(Sys.getenv("NB", "20"))
DATA <- Sys.getenv("BENCH_DATA", "data")

quiet <- function(e) {
  tf <- tempfile(); sink(tf)
  out <- withCallingHandlers(e, message = function(m) invokeRestart("muffleMessage"),
                             warning = function(w) invokeRestart("muffleWarning"))
  sink(); unlink(tf); out
}

## Species, fixture, and the export each one must be compared against. Checking the export's
## MANIFEST is not optional: a human database against a mouse fixture makes the *kernel* fail
## correctly and the run look broken.
CASES <- list(
  list(species = "human", file = "humanSkin.rda", db = "db_human"),
  list(species = "mouse", file = "wound.rda",      db = "db_mouse")
)

failures <- 0L
checked <- 0L
missing <- character()
for (case in CASES) {
  load_fixture_for(case$species)
  Sys.setenv(CELLCHATRS_DB = file.path(getwd(), "tests", "fixtures", case$db))
  path <- file.path(DATA, case$file)
  if (!file.exists(path)) {
    missing <- c(missing, path)
    cat(sprintf("parity %-14s MISSING (no fixture at %s)\n", case$file, path))
    next
  }
  checked <- checked + 1L
  built <- quiet(build_fixture(path, max_lr = 0L))
  info <- built$info
  cpp <- get("computeCommunProb", envir = quiet(upstream_env()))
  a <- quiet(cpp(built$object, nboot = NB, seed.use = 1L))
  b <- quiet(computeCommunProb(built$object, nboot = NB, seed.use = 1L))
  n <- length(as.vector(a@net$prob))
  ## The whole `net`, which is what the S4 object stores: prob, pval, and any other member the
  ## kernel fills in. A `prob`-only comparison would miss a pval divergence.
  ok <- identical(a@net, b@net)
  maxd <- suppressWarnings(max(abs(a@net$prob - b@net$prob)))
  cat(sprintf("parity %-14s nC=%-6d K=%-3d nLR=%-5d values=%-8d net=%-5s max|diff|=%g\n",
              case$file, info$n_cells, info$n_groups, info$n_lr, n, ok, maxd))
  if (!ok || maxd != 0) {
    failures <- failures + 1L
    for (nm in union(names(a@net), names(b@net))) {
      if (identical(a@net[[nm]], b@net[[nm]])) next
      cat("   differs:", nm, "\n")
      print(all.equal(a@net[[nm]], b@net[[nm]]))
    }
  }
}
if (length(missing) > 0L || checked != length(CASES)) {
  cat(sprintf("parity: required real-data fixtures missing (%d/%d checked):\n  %s\n",
              checked, length(CASES), paste(missing, collapse = "\n  ")))
  quit(status = 1L)
}
if (failures > 0L) {
  cat(sprintf("parity: %d fixture(s) diverged from upstream\n", failures))
  quit(status = 1L)
}
cat("parity: all available real-data fixtures are identical to upstream\n")
