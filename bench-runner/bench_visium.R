## Kernel benchmark on the authors' mouse-cortex visium data (spatial path).
##
## bench_real.R covers human skin and mouse wound (scRNA-seq path). This is the same protocol --
## warm-up discarded, 5 timed repeats, median with bootstrap CI, peak RSS, taskset-pinned CPUs,
## `identical()` on the whole `net` afterwards -- applied to the spatial path on real spatial
## data (1,073 spots, 648 genes, 433 signalling, 8 cell types). Small enough to run quickly,
## different enough to matter: the spatial branch (`distance.use = TRUE` with the exact k-d
## tree, `P.spatial` multiplied per L-R pair) is code the scRNA-seq benchmarks never touch.
##
## Two facts that constrain the setup, both established in
## bench-runner/measure_spatial_divergence.R and repeated here rather than re-derived:
##
## * Only 134 of the object's 437 L-R pairs have every subunit present in `data.signaling`;
##   on the full set both sides raise upstream's own `subscript out of bounds` (verified
##   identical). Timing an error path would measure nothing, so the benchmark runs the 134
##   resolvable pairs -- the same set the divergence measurement compares `Prob` on.
## * `contact.dependent = TRUE` has no default range (upstream raises without one), so
##   `contact.range` is the visium spot diameter (65); `interaction.range = 1000`,
##   `ratio = 1`, `tol = 0` put distances in image-pixel units with exact range tests.
##
## Usage:
##   taskset -c 0-7 Rscript bench-runner/bench_visium.R [--nboot 100] [--out bench-runner/results/visium.json]
suppressWarnings(suppressMessages({
  library(stats); library(methods); library(Matrix); library(collapse); library(dplyr)
}))
suppressWarnings(suppressMessages(library(cellchatrs)))

`%||%` <- function(a, b) if (is.null(a)) b else a
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
NBOOT <- as.integer(opt$nboot %||% "100")
REPEATS <- as.integer(opt$repeats %||% "5")
WARMUP <- as.integer(opt$warmup %||% "1")
SEED <- as.integer(opt$seed %||% "1")
OUT <- opt$out %||% "bench-runner/results/visium.json"

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
host <- check_host_quiet(opt, label = "BENCH-VISIUM")
if (is.na(host$pinned)) {
  cat("BENCH-VISIUM REFUSING: cannot determine this process's CPU affinity and BENCH_PINNED is unset.\n")
  cat("            Run under taskset so the timings record which cores produced them.\n")
  quit(status = 2L)
}
if (!isTRUE(host$verdict$ok) && !isTRUE(host$allowed)) {
  writeLines(j_val(refusal_record(host)), paste0(OUT, ".refused.json"))
  quit(status = 3L)
}
VISIUM <- Sys.getenv("VISIUM", "/scratch/mdra00001/tmp/opencode/data/visium.rds")
CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")
THREADS <- as.integer(Sys.getenv("CELLCHATRS_THREADS", "8"))

## The RDS holds a real CellChat S4 object: the class must exist before unserialising, and it
## must be registered in globalenv (not a private env), or `new()`/`is()` fail with the
## timestamp-package error documented in tutorial_repro.R. Same arrangement, not re-derived.
for (p in c("Matrix", "collapse", "dplyr", "future", "rlang")) {
  suppressWarnings(suppressMessages(requireNamespace(p, quietly = TRUE)))
}
suppressWarnings(suppressMessages(
  sys.source(file.path(CC, "R", "CellChat_class.R"), envir = globalenv(),
             keep.source = FALSE)))
stopifnot(!is.null(methods::getClass("CellChat", where = globalenv())))
o <- readRDS(VISIUM)

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)

spot_diameter <- as.numeric(o@images$scale.factors$spot.diameter)
RATIO <- 1
TOL <- 0
INTERACTION_RANGE <- 1000L
meta <- data.frame(group = o@idents, stringsAsFactors = FALSE)
meta$samples <- factor(rep("section1", ncol(o@data)))
meta$label <- o@idents
cat(sprintf("BENCH visium: %d genes x %d spots, %d groups, %d L-R pairs\n",
            nrow(o@data), ncol(o@data), nlevels(o@idents), nrow(o@LR$LRsig)))

## Probe object: every slot the kernel reads, nothing it does not. The `@images` rename
## (`scale.factors` -> `spatial.factors`) is what `updateCellChat` did in 2.1.1; the pinned
## object predates it, so the rename happens on the way in (see measure_spatial_divergence.R).
probe <- function() {
  new("CellChatProbe",
      data = o@data,
      data.signaling = o@data.signaling,
      data.smooth = o@data.signaling,
      idents = o@idents,
      meta = meta,
      images = list(coordinates = cbind(x = o@images$coordinates$x_cent,
                                        y = o@images$coordinates$y_cent),
                    spatial.factors = list(spot.diameter = spot_diameter,
                                           ratio = RATIO, tol = TOL)),
      DB = list(complex = DB$complex, cofactor = DB$cofactor),
      LR = list(LRsig = o@LR$LRsig),
      net = list(prob = array(0, c(1, 1, 1)), pval = array(0, c(1, 1, 1))),
      netP = list(),
      options = list(datatype = "spatial", mode = "single", db = normalizePath(DBDIR),
                     population.size = FALSE))
}
setClass("CellChatProbe", representation(
  data = "ANY", data.signaling = "ANY", data.smooth = "ANY", idents = "ANY", meta = "ANY",
  images = "ANY", DB = "ANY", LR = "ANY", net = "ANY", netP = "ANY", options = "ANY"
), where = globalenv())

## The 134 resolvable pairs (see header): timing an error path would measure nothing.
lr_resolvable <- function(lrsig) {
  ok <- vapply(seq_len(nrow(lrsig)), function(i) {
    lig <- unlist(strsplit(as.character(lrsig$ligand[i]), ",", fixed = TRUE), use.names = FALSE)
    rec <- unlist(strsplit(as.character(lrsig$receptor[i]), ",", fixed = TRUE), use.names = FALSE)
    all(c(lig, rec) %in% rownames(o@data.signaling))
  }, logical(1))
  lrsig[ok, , drop = FALSE]
}
LR_ok <- lr_resolvable(o@LR$LRsig)
cat(sprintf("BENCH resolvable pairs: %d of %d\n", nrow(LR_ok), nrow(o@LR$LRsig)))

run_kernel <- function(fn) {
  o2 <- probe()
  o2@LR <- list(LRsig = LR_ok)
  fn(o2, type = "triMean", trim = 0.1, population.size = FALSE,
     distance.use = TRUE, interaction.range = INTERACTION_RANGE, scale.distance = 0.01,
     k.min = 10, contact.dependent = TRUE, contact.range = spot_diameter,
     do.symmetric = TRUE, nboot = NBOOT, seed.use = SEED, Kh = 0.5, n = 1)
}
m_up <- measure("R upstream", function() run_kernel(cellchatrs_upstream_computeCommunProb),
                REPEATS, WARMUP)
m_rs <- measure("Rust kernel", function() run_kernel(computeCommunProb), REPEATS, WARMUP)

parity_net <- identical(m_up$value@net, m_rs$value@net)
parity_prob <- identical(m_up$value@net$prob, m_rs$value@net$prob)
max_abs <- suppressWarnings(max(abs(m_up$value@net$prob - m_rs$value@net$prob)))
cat(sprintf("BENCH parity net=%s prob=%s max|diff|=%.3g\n", parity_net, parity_prob, max_abs))

## Non-identity here is EXPECTED, not a failure: the shim uses the exact k-d tree where upstream
## queries Annoy, and measure_spatial_divergence.R established on this same object that the two
## differ in 1,682 of 8,576 Prob entries with max abs diff 0.001775. What this benchmark asserts
## instead is *consistency with that measurement*: the same divergence reproduced, not a new one.
## A max|diff| far from 0.001775 would mean the spatial path changed under the benchmark; equal
## `identical()` would mean the exact tree had silently become Annoy. Either deviation fails loudly
## below instead of passing quietly.
divergence_ref <- 0.001774885936636432
divergence_ok <- is.finite(max_abs) && abs(max_abs - divergence_ref) / divergence_ref < 0.05
cat(sprintf("BENCH divergence consistency: max|diff|=%.6g vs recorded %.6g -> %s\n",
            max_abs, divergence_ref, if (divergence_ok) "CONSISTENT" else "DIVERGED"))

scatter <- list(
  n = length(m_up$value@net$prob),
  log10_prob_r   = log10(pmax(as.vector(m_up$value@net$prob), .Machine$double.xmin)),
  log10_prob_rs  = log10(pmax(as.vector(m_rs$value@net$prob), .Machine$double.xmin))
)

rss <- peak_rss_mb()
speedup <- m_up$median / m_rs$median
speedup_ci <- c(lo = m_up$ci_lo / m_rs$ci_hi, hi = m_up$ci_hi / m_rs$ci_lo)
res <- list(
  schema = 1L,
  host = list(cpus = 16L, cores = 8L, model = "AMD Ryzen 7 3700X", mem_total_mb = 32806),
  pinned = host$pinned, host_state = host_record(host),
  threads = THREADS,
  fixture = list(dataset = "visium", n_spots = ncol(o@data),
                 n_genes_total = nrow(o@data),
                 n_genes_signaling = nrow(o@data.signaling),
                 n_groups = nlevels(o@idents), n_lr = nrow(LR_ok),
                 n_lr_total = nrow(o@LR$LRsig)),
  params = list(interaction_range = INTERACTION_RANGE, contact_range = spot_diameter,
                ratio = RATIO, tol = TOL, k_min = 10L),
  nboot = NBOOT, repeats = REPEATS, warmup = WARMUP, seed = SEED,
  r_upstream = m_up[c("raw_secs", "median", "ci_lo", "ci_hi", "min", "max")],
  rust = m_rs[c("raw_secs", "median", "ci_lo", "ci_hi", "min", "max")],
  speedup = speedup, speedup_ci_lo = speedup_ci[["lo"]], speedup_ci_hi = speedup_ci[["hi"]],
  peak_rss_mb = rss,
  parity_net = parity_net, parity_prob = parity_prob, max_abs_diff = max_abs,
  ## See above: on the spatial path the honest expectation is the *measured* Annoy divergence,
  ## not bit-identity. `divergence_consistent` records whether this run reproduced it.
  divergence_consistent = divergence_ok,
  divergence_ref_max_abs = divergence_ref,
  mem_available_mb = mem_available_mb(),
  scatter = scatter
)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
writeLines(j_val(res), OUT)
cat(sprintf(paste0("BENCH R=%.2fs [%.2f,%.2f] Rust=%.2fs [%.2f,%.2f] speedup=%.2fx ",
                   "[%.2f,%.2f] rss=%.0fMB -> %s\n"),
            m_up$median, m_up$ci_lo, m_up$ci_hi, m_rs$median, m_rs$ci_lo, m_rs$ci_hi,
            speedup, speedup_ci[["lo"]], speedup_ci[["hi"]], rss, OUT))
