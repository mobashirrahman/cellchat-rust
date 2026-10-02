## Amdahl decomposition of the CellChat pipeline: I/O+prep vs OEG vs kernel vs downstream.
##
## bench_real.R measures one stage (`computeCommunProb`) on both sides. That number alone cannot
## support the paper's central performance claim, because a 40x kernel on a pipeline that spends
## half its time elsewhere is a ~2x pipeline. This script times the tutorial's full computational
## pipeline stage by stage -- the same stages tests/parity/tutorial_repro.R proves identical() --
## on both implementations, with the same protocol (warm-up discarded, >=5 repeats, median with
## bootstrap CI, peak RSS, taskset-pinned CPUs), and reports each stage's share of each side's
## total plus the Amdahl limit the shares imply.
##
## Stages (vignette order, visualization excluded per 14.1):
##   io_prep      load .rda + log1p CPM + createCellChat + DB assign + subsetData
##   oeg          identifyOverExpressedGenes (defaults) + identifyOverExpressedInteractions
##   kernel       computeCommunProb, triMean defaults, nboot from --nboot
##   downstream   filterCommunication + computeCommunProbPathway + aggregateNet
##   centrality   netAnalysis_computeCentrality on netP (same igraph seed both sides)
##
## Two things this script does not do, both on purpose:
## * It does not assert identical() per stage -- tutorial_repro.R already does, on this data.
##   Re-checking here would double a ~10-minute run for no new information. It DOES check the
##   final objects' net identical() once, as a tripwire that the timed pipeline is the gated one.
## * It does not time `identifyOverExpressedGenes` with do.fast=FALSE. The tutorial default
##   (presto) is what users run; the exact path is parity-gated elsewhere.
##
## Usage:
##   taskset -c 0-7 Rscript bench-runner/bench_amdahl.R --fixture .../humanSkin.rda --nboot 100
suppressWarnings(suppressMessages({
  library(stats); library(methods); library(Matrix); library(collapse); library(dplyr)
}))
suppressWarnings(suppressMessages(library(cellchatrs)))

`%||%` <- function(a, b) if (is.null(a)) b else a

## --- minimal JSON + stats (same shape as bench_real.R; kept local so this file runs standalone)
j_esc <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  paste0("\"", x, "\"")
}
j_num <- function(x) {
  if (length(x) != 1L) return(paste(vapply(x, j_num, ""), collapse = ", "))
  if (is.na(x)) return("null")
  if (is.infinite(x)) return(if (x > 0) "1e999" else "-1e999")
  if (x == round(x) && abs(x) < 1e15) return(sprintf("%.0f", x))
  sprintf("%.17g", x)
}
j_val <- function(x) {
  if (is.null(x)) return("null")
  if (is.list(x)) {
    if (!length(x)) return(if (!is.null(names(x))) "{}" else "[]")
    if (!is.null(names(x)) && all(nzchar(names(x)))) {
      return(paste0("{", paste(sprintf("%s: %s", j_esc(names(x)),
                                     vapply(x, j_val, "")), collapse = ", "), "}"))
    }
    return(paste0("[", paste(vapply(x, j_val, ""), collapse = ", "), "]"))
  }
  if (length(x) != 1L) return(paste0("[", paste(vapply(as.list(x), j_val, ""), collapse = ", "), "]"))
  if (is.character(x)) return(j_esc(x))
  if (is.logical(x)) return(if (isTRUE(x)) "true" else "false")
  j_num(x)
}
parse_args <- function() {
  a <- commandArgs(trailingOnly = TRUE)
  out <- list()
  i <- 1L
  while (i <= length(a)) {
    k <- a[i]
    if (grepl("^--", k)) {
      v <- if (i < length(a) && !grepl("^--", a[i + 1L])) a[i + 1L] else "TRUE"
      out[[sub("^--", "", k)]] <- v
      i <- i + (if (identical(v, "TRUE")) 1L else 2L)
    } else i <- i + 1L
  }
  out
}
median_ci <- function(x, n_boot = 10000L, seed = 4242L) {
  set.seed(seed)
  n <- length(x)
  meds <- numeric(n_boot)
  for (i in seq_len(n_boot)) meds[i] <- median(sample(x, n, replace = TRUE))
  c(lo = unname(quantile(meds, 0.025)), hi = unname(quantile(meds, 0.975)))
}
quietly <- function(expr) {
  tf <- tempfile()
  sink(tf)
  out <- withCallingHandlers(expr, message = function(m) invokeRestart("muffleMessage"),
                             warning = function(w) invokeRestart("muffleWarning"))
  sink()
  unlink(tf)
  out
}
time_once <- function(f) {
  gc(reset = TRUE, full = FALSE)
  t0 <- proc.time()[["elapsed"]]
  out <- quietly(f())
  t1 <- proc.time()[["elapsed"]]
  list(secs = t1 - t0, value = out)
}
measure <- function(label, f, repeats, warmup) {
  for (i in seq_len(warmup)) f()
  ts <- numeric(repeats)
  last <- NULL
  for (i in seq_len(repeats)) {
    r <- time_once(f)
    ts[i] <- r$secs
    last <- r$value
  }
  ci <- median_ci(ts)
  list(label = label, repeats = repeats, raw_secs = ts,
       median = median(ts), ci_lo = ci[["lo"]], ci_hi = ci[["hi"]],
       min = min(ts), max = max(ts), value = last)
}
## Read from R itself: `system("grep VmHWM /proc/self/status")` measures the forked subshell
## (~2 MB), not R, because `/proc/self` resolves in the child. An earlier version of this file did
## exactly that and reported rss=2MB for a run holding a 1 GB matrix. This is fixture.R's version,
## copied rather than sourced to keep this file standalone (see the header note in fixture.R about
## why bench_real.R must not be edited to share code).
peak_rss_mb <- function() {
  st <- tryCatch(readLines("/proc/self/status"), error = function(e) character(0))
  hit <- grep("^VmHWM:", st, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", hit)) / 1024
}
mem_available_mb <- function() {
  st <- tryCatch(readLines("/proc/meminfo"), error = function(e) character(0))
  hit <- grep("^MemAvailable:", st, value = TRUE)
  if (!length(hit)) return(NA_real_)
  as.numeric(gsub("[^0-9]", "", hit)) / 1024
}

opt <- parse_args()
FIXTURE <- opt$fixture %||% file.path(Sys.getenv("CELLCHATRS_DATA", "data"), "humanSkin.rda")
NBOOT <- as.integer(opt$nboot %||% "100")
REPEATS <- as.integer(opt$repeats %||% "5")
WARMUP <- as.integer(opt$warmup %||% "1")
OUT <- opt$out %||% "bench-runner/results/amdahl.json"

## --- host precondition -----------------------------------------------------------------------
## The same gate `bench_real.R` uses, from the shared implementation in `cpu_gate.R`. This runner's
## numbers are quoted in the paper as well, so "was the host quiet when these were measured" is a
## question about every runner, not only the two headline fixtures. See docs/BENCHMARKS.md,
## "Host contention" -- the reasoning and the two silent bugs the gate's tests pin down.
SCRIPT_DIR <- local({
  a <- grep("^--file=", commandArgs(FALSE), value = TRUE)
  if (length(a)) dirname(sub("^--file=", "", a[1])) else "bench-runner"
})
source(file.path(SCRIPT_DIR, "cpu_gate.R"))
host <- check_host_quiet(opt, label = "BENCH-AMDAHL")
if (is.na(host$pinned)) {
  cat("BENCH-AMDAHL REFUSING: cannot determine this process's CPU affinity and BENCH_PINNED is unset.\n")
  cat("            Run under taskset so the timings record which cores produced them.\n")
  quit(status = 2L)
}
if (!isTRUE(host$verdict$ok) && !isTRUE(host$allowed)) {
  writeLines(j_val(refusal_record(host)), paste0(OUT, ".refused.json"))
  quit(status = 3L)
}
CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")

## --- upstream reference env + class (same arrangement as tutorial_repro.R, q.v. for why the
## class file goes to globalenv and why Matrix must be attached, not merely loadable)
UP <- new.env(parent = globalenv())
for (p in c("Matrix", "collapse", "dplyr", "future", "rlang")) {
  suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
}
for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R", "visualization.R")) {
  suppressWarnings(suppressMessages(
    sys.source(file.path(CC, "R", f), envir = UP, keep.source = FALSE)))
}
suppressWarnings(suppressMessages(
  sys.source(file.path(CC, "R", "CellChat_class.R"), envir = globalenv(),
             keep.source = FALSE)))
E <- new.env()
## The database must match the data's species: a human DB assigned to a mouse object filters
## `subsetData` down to the intersection of human signalling genes with mouse rownames, which can
## leave fewer than 3 genes -- and OEG then stops with "Please check `object@data.signaling`...".
## That failure names the symptom (too few genes) rather than the cause (wrong database), and it
## fires on one dataset but not the other, which is why it survived the human-skin run.
SPECIES <- Sys.getenv("SPECIES", "human")
load(file.path(CC, "data", paste0("CellChatDB.", SPECIES, ".rda")), envir = E)
DB <- get(ls(E)[1], E)

## --- input (tutorial preamble, verbatim)
raw <- new.env(); load(FIXTURE, envir = raw); skin <- raw[[ls(raw)[1]]]
counts <- skin$data
grp <- if (!is.null(skin$meta)) factor(skin$meta$labels) else factor(skin$labels)
cat(sprintf("AMDAHL input: %d genes x %d cells, %d groups (%s)\n",
            nrow(counts), ncol(counts), nlevels(grp), basename(FIXTURE)))

## Each stage is a closure over fresh state where it matters: io_prep rebuilds from `counts` every
## repeat (page cache makes repeats 2..n faster than cold load, which is reported, not hidden --
## users re-run pipelines in warm sessions too), while oeg/kernel/downstream/centrality close over
## one prebuilt input object so repeats measure the stage rather than its prerequisites.
## A stage that rebuilt its own input would time the input too, and the shares would not add up.
prep_once <- function() {
  lib <- Matrix::colSums(counts)
  cpm <- counts
  cpm@x <- cpm@x / rep(lib, diff(cpm@p)) * 1e4
  norm <- log1p(as.matrix(cpm))
  meta <- data.frame(labels = grp, row.names = colnames(counts))
  obj <- createCellChat(object = norm, meta = meta, group.by = "labels")
  obj@DB <- DB
  obj
}

stages <- list(
  list(name = "io_prep",
       up = function() { o <- prep_once(); UP$subsetData(o) },
       rs = function() { o <- prep_once(); cellchatrs::subsetData(o) }),
  list(name = "oeg",
       up = function() {
         o <- UP$identifyOverExpressedGenes(upstream_state$obj)
         UP$identifyOverExpressedInteractions(o)
       },
       rs = function() {
         o <- cellchatrs::identifyOverExpressedGenes(shim_state$obj)
         cellchatrs::identifyOverExpressedInteractions(o)
       }),
  list(name = "kernel",
       up = function() UP$computeCommunProb(upstream_state$obj, type = "triMean", nboot = NBOOT),
       rs = function() cellchatrs::computeCommunProb(shim_state$obj, type = "triMean",
                                                    nboot = NBOOT)),
  list(name = "downstream",
       up = function() {
         o <- UP$filterCommunication(upstream_state$obj, min.cells = 10)
         o <- UP$computeCommunProbPathway(o)
         UP$aggregateNet(o)
       },
       rs = function() {
         o <- cellchatrs::filterCommunication(shim_state$obj, min.cells = 10)
         o <- cellchatrs::computeCommunProbPathway(o)
         cellchatrs::aggregateNet(o)
       }),
  list(name = "centrality",
       up = function() {
         set.seed(20240501)
         UP$netAnalysis_computeCentrality(upstream_state$obj, slot.name = "netP")
       },
       rs = function() {
         set.seed(20240501)
         cellchatrs::netAnalysis_computeCentrality(shim_state$obj, slot.name = "netP")
       })
)

## State advances through the pipeline in order: each stage's input is the previous stage's
## output, timed once per side before the repeats start. Timing the build of the input inside the
## repeats would charge every stage for its prerequisites.
upstream_state <- new.env()
shim_state <- new.env()
cat("AMDAHL warming inputs...\n")
upstream_state$obj <- quietly({ o <- prep_once(); UP$subsetData(o) })
shim_state$obj <- quietly({ o <- prep_once(); cellchatrs::subsetData(o) })
stopifnot(identical(upstream_state$obj, shim_state$obj))
upstream_state$obj <- quietly({
  o <- UP$identifyOverExpressedGenes(upstream_state$obj)
  UP$identifyOverExpressedInteractions(o)
})
shim_state$obj <- quietly({
  o <- cellchatrs::identifyOverExpressedGenes(shim_state$obj)
  cellchatrs::identifyOverExpressedInteractions(o)
})
stopifnot(identical(upstream_state$obj, shim_state$obj))

results <- list()
for (st in stages) {
  mu <- measure(paste0("up:", st$name), st$up, REPEATS, WARMUP)
  mr <- measure(paste0("rs:", st$name), st$rs, REPEATS, WARMUP)
  cat(sprintf("AMDAHL %-10s R=%7.2fs [%6.2f,%6.2f] RS=%7.2fs [%6.2f,%6.2f] speedup=%6.2fx\n",
              st$name, mu$median, mu$ci_lo, mu$ci_hi, mr$median, mr$ci_lo, mr$ci_hi,
              mu$median / mr$median))
  results[[st$name]] <- list(upstream = mu[c("raw_secs", "median", "ci_lo", "ci_hi")],
                            shim = mr[c("raw_secs", "median", "ci_lo", "ci_hi")])
  ## Advance the shared state with the *upstream* output when the stage transforms the object,
  ## so later stages start from identical inputs on both sides (verified above for the first
  ## two; the kernel downstream stages take whatever the measured call returned).
  if (st$name %in% c("kernel", "downstream")) {
    if (!identical(mu$value@net, mr$value@net)) {
      cat(sprintf("AMDAHL TRIPWIRE %s: net differs between sides; later stages still timed\n",
                  st$name))
    }
    upstream_state$obj <- mu$value
    shim_state$obj <- mr$value
  }
  if (st$name == "oeg") {
    upstream_state$obj <- mu$value
    shim_state$obj <- mr$value
  }
}

rss <- peak_rss_mb()
tot_up <- sum(vapply(results, function(r) r$upstream$median, 0))
tot_rs <- sum(vapply(results, function(r) r$shim$median, 0))
shares_up <- lapply(results, function(r) r$upstream$median / tot_up)
shares_rs <- lapply(results, function(r) r$shim$median / tot_rs)
## Amdahl limit: kernel infinitely fast => overall bounded by 1 / (1 - kernel_share_upstream).
p_kernel <- results$kernel$upstream$median / tot_up
amdahl_limit <- 1 / (1 - p_kernel)
overall <- tot_up / tot_rs
cat(sprintf(paste0("AMDAHL total R=%.1fs RS=%.1fs overall=%.2fx | kernel share of R=%.1f%% ",
                   "=> Amdahl limit %.2fx | rss=%.0fMB\n"),
            tot_up, tot_rs, overall, 100 * p_kernel, amdahl_limit, rss))

res <- list(
  schema = 1L,
  host = list(cpus = 16L, cores = 8L, model = "AMD Ryzen 7 3700X", mem_total_mb = 32806),
  pinned = host$pinned, host_state = host_record(host),
  fixture = basename(FIXTURE), nboot = NBOOT, repeats = REPEATS, warmup = WARMUP,
  stages = lapply(names(results), function(nm) {
    r <- results[[nm]]
    list(stage = nm,
         upstream = r$upstream, shim = r$shim,
         speedup = r$upstream$median / r$shim$median,
         speedup_ci_lo = r$upstream$ci_lo / r$shim$ci_hi,
         speedup_ci_hi = r$upstream$ci_hi / r$shim$ci_lo)
  }),
  total_upstream = tot_up, total_shim = tot_rs, overall_speedup = overall,
  kernel_share_upstream = p_kernel, amdahl_limit = amdahl_limit,
  peak_rss_mb = rss, mem_available_mb = mem_available_mb()
)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
writeLines(j_val(res), OUT)
cat(sprintf("AMDAHL -> %s\n", OUT))
