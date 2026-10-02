## Controlled synthetic grid: kernel time and parity as a function of cell count.
##
## The author datasets fix nC at 7,563 and 21,557. A scaling claim needs more than two points,
## and the downsampling correction ("subset cells to fit runtime") needs evidence ABOVE the sizes
## the authors ship -- otherwise it is an extrapolation, not a measurement. This script builds
## synthetic datasets on a grid of cell counts with everything else held fixed (K groups, nLR
## pairs, nGenes signaling genes, nboot), times upstream R against the shim on each point with the
## standard protocol, and checks `identical()` per point.
##
## Design constraints, each load-bearing:
## * Real L-R pairs from the bundled human export (single-gene ligands/receptors, like the
##   vignette picker), so complex/cofactor resolution runs for real. Invented pair names would
##   resolve against nothing and time an error path.
## * Counts are negative-binomial-ish (rpois with group-specific logFC on a DE subset), so the
##   matrix has the sparsity structure of scRNA-seq rather than uniform noise. Uniform noise
##   makes every gene uninformative and the bootstrap trivially fast -- a benchmark that flatters
##   by being easy.
## * The memory-headroom gate from the objective is enforced IN CODE, not in prose: no grid
##   point at or above 50,000 cells runs unless the host reports enough available memory for the
##   dense buffer times a safety factor. Quoting a 50k number from a run that swapped would be
##   measuring the page cache, not the kernel.
## * Parity per grid point, not just at the end: bit-identity must not degrade with size (e.g.
##   via an accumulation order that only shows at scale), so each point asserts it separately.
##
## Usage:
##   taskset -c 0-7 Rscript bench-runner/bench_synth.R [--grid 2000,5000,10000,20000,50000]
##     [--nboot 100] [--out bench-runner/results/synth_grid.json]
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
GRID <- as.integer(strsplit(opt$grid %||% "2000,5000,10000,20000,50000", ",", fixed = TRUE)[[1]])
NBOOT <- as.integer(opt$nboot %||% "100")
REPEATS <- as.integer(opt$repeats %||% "5")
WARMUP <- as.integer(opt$warmup %||% "1")
SEED <- as.integer(opt$seed %||% "7")
OUT <- opt$out %||% "bench-runner/results/synth_grid.json"

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
host <- check_host_quiet(opt, label = "BENCH-SYNTH")
if (is.na(host$pinned)) {
  cat("BENCH-SYNTH REFUSING: cannot determine this process's CPU affinity and BENCH_PINNED is unset.\n")
  cat("            Run under taskset so the timings record which cores produced them.\n")
  quit(status = 2L)
}
if (!isTRUE(host$verdict$ok) && !isTRUE(host$allowed)) {
  writeLines(j_val(refusal_record(host)), paste0(OUT, ".refused.json"))
  quit(status = 3L)
}
K_GROUPS <- 6L
N_LR <- 200L
N_GENES <- 320L
HEADROOM_FACTOR <- 3.0

CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")

## Upstream reference env (function bodies only; class globally, per tutorial_repro.R).
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

## --- synthetic data with scRNA-seq structure -------------------------------------------
## 200 real pairs (single-gene ends, non-empty annotation) from the bundled export; the gene
## universe is their subunits plus filler. Group-specific logFC on 60 genes gives the bootstrap
## something to find; the rest is background. Counts are Poisson around group means -- sparse,
## overdispersed across groups, nothing like uniform noise.
E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DBFULL <- get(ls(E)[1], E)
## A ligand/receptor symbol that coincides with a *complex* name is resolved as the complex, whose
## subunits are then required to be rows of the matrix. On real data they are; here only the picked
## genes exist. Excluding such pairs is not dumbing the fixture down -- it is the same condition
## the kernel enforces everywhere ("subscript out of bounds" otherwise), applied at selection time
## rather than discovered mid-run.
CXNAMES <- rownames(DBFULL$complex)
dbix <- system.file("db", "human", "interactions.tsv", package = "cellchatrs")
if (!nzchar(dbix)) stop("bundled human DB export missing; reinstall the package", call. = FALSE)
dblines <- strsplit(readLines(dbix), "\t", fixed = TRUE)
is_simple <- vapply(dblines, function(f) {
  length(f) >= 9L && nzchar(f[3]) && nzchar(f[4]) &&
    !grepl("[,;_]", f[3]) && !grepl("[,;_]", f[4]) && nzchar(f[9]) &&
    !(f[3] %in% CXNAMES) && !(f[4] %in% CXNAMES)
}, TRUE)
cand <- which(is_simple)
if (length(cand) < N_LR) stop("only ", length(cand), " simple pairs; need ", N_LR, call. = FALSE)
set.seed(SEED)
take <- sort(sample(cand, N_LR))
mkrow <- function(i) {
  f <- dblines[[i]]
  data.frame(interaction_name = f[1], pathway_name = f[2], ligand = f[3], receptor = f[4],
             agonist = "", antagonist = "", co_A_receptor = "", co_I_receptor = "",
             annotation = f[9], stringsAsFactors = FALSE)
}
LRsig <- do.call(rbind, lapply(take, mkrow))
rownames(LRsig) <- LRsig$interaction_name
sig_genes <- sort(unique(c(LRsig$ligand, LRsig$receptor)))
genes <- c(sig_genes, paste0("BG", seq_len(N_GENES - length(sig_genes))))
cat(sprintf("SYNTH gene universe: %d genes (%d signaling), %d pairs\n",
            length(genes), length(sig_genes), nrow(LRsig)))

build_object <- function(nC, make) {
  ## make = "up" or "rs": the two sides need structurally identical but independent objects,
  ## because computeCommunProb mutates slots in place downstream (filterCommunication writes
  ## net$prob). Sharing one object between sides would time the second side on the first side's
  ## output. A shared builder with a side flag keeps them identical by construction.
  set.seed(SEED + nC)
  grp <- factor(rep(seq_len(K_GROUPS), length.out = nC),
                levels = as.character(seq_len(K_GROUPS)))
  base <- rpois(length(genes) * nC, lambda = 2)
  X <- matrix(base, nrow = length(genes),
              dimnames = list(genes, paste0("c", seq_len(nC))))
  ## DE subset: 60 genes up 4x in their own group (cycled), so group means genuinely differ and
  ## the kernel aggregates something other than background. Without this the bootstrap compares
  ## identical distributions and the timing measures the fast path only.
  de <- genes[seq(1, min(360L, length(genes)), by = 6)]
  for (g in seq_along(de)) {
    grp_idx <- ((g - 1L) %% K_GROUPS) + 1L
    cells <- which(as.integer(grp) == grp_idx)
    X[de[g], cells] <- rpois(length(cells), lambda = 8)
  }
  lib <- colSums(X) + 1
  C <- X / rep(lib, each = nrow(X)) * 1e4
  N <- log1p(as.matrix(C))
  meta <- data.frame(labels = grp, row.names = colnames(N))
  obj <- createCellChat(object = N, meta = meta, group.by = "labels")
  ## Tutorial-verbatim: the vignette assigns the database next, because createCellChat leaves
  ## `@DB` empty and subsetData selects from it.
  obj@DB <- DBFULL
  obj <- UP$subsetData(obj)
  obj@LR$LRsig <- LRsig
  obj
}

results <- list()
for (nC in GRID) {
  ## The headroom gate, enforced in code: dense buffer (genes x cells x 8 bytes) times the safety
  ## factor must fit in MemAvailable, or the point is refused. At or above 50,000 cells this gate
  ## is what the objective requires before any number may be quoted; below it the same check runs
  ## and is reported, so a small-grid number from a loaded host is not silently trusted either.
  need_mb <- length(genes) * nC * 8 / 1024^2 * HEADROOM_FACTOR
  avail_mb <- mem_available_mb()
  ok <- is.finite(avail_mb) && avail_mb > need_mb
  cat(sprintf("SYNTH nC=%d need=%.0fMB avail=%.0fMB headroom=%s\n",
              nC, need_mb, avail_mb, if (ok) "ok" else "REFUSED"))
  if (!ok) {
    results[[as.character(nC)]] <- list(n_cells = nC, refused = TRUE,
                                        need_mb = need_mb, avail_mb = avail_mb)
    next
  }
  obj_up <- build_object(nC, "up")
  obj_rs <- build_object(nC, "rs")
  stopifnot(identical(obj_up, obj_rs))
  nlr_here <- nrow(obj_up@LR$LRsig)
  mu <- measure("R upstream", function() UP$computeCommunProb(
    obj_up, type = "triMean", nboot = NBOOT, seed.use = SEED), REPEATS, WARMUP)
  mr <- measure("Rust kernel", function() cellchatrs::computeCommunProb(
    obj_rs, type = "triMean", nboot = NBOOT, seed.use = SEED), REPEATS, WARMUP)
  par_net <- identical(mu$value@net, mr$value@net)
  par_prob <- identical(mu$value@net$prob, mr$value@net$prob)
  max_abs <- suppressWarnings(max(abs(mu$value@net$prob - mr$value@net$prob)))
  cat(sprintf(paste0("SYNTH nC=%d K=%d nLR=%d R=%.2fs [%.2f,%.2f] RS=%.2fs [%.2f,%.2f] ",
                     "speedup=%.2fx parity_net=%s parity_prob=%s max|diff|=%.3g\n"),
              nC, K_GROUPS, nlr_here,
              mu$median, mu$ci_lo, mu$ci_hi, mr$median, mr$ci_lo, mr$ci_hi,
              mu$median / mr$median, par_net, par_prob, max_abs))
  results[[as.character(nC)]] <- list(
    n_cells = nC, n_groups = K_GROUPS, n_lr = nlr_here, n_genes = length(genes),
    nboot = NBOOT, refused = FALSE, need_mb = need_mb, avail_mb = avail_mb,
    upstream = mu[c("raw_secs", "median", "ci_lo", "ci_hi", "min", "max")],
    shim = mr[c("raw_secs", "median", "ci_lo", "ci_hi", "min", "max")],
    speedup = mu$median / mr$median,
    speedup_ci_lo = mu$ci_lo / mr$ci_hi, speedup_ci_hi = mu$ci_hi / mr$ci_lo,
    parity_net = par_net, parity_prob = par_prob, max_abs_diff = max_abs)
  rm(obj_up, obj_rs, mu, mr); gc(full = FALSE)
}

rss <- peak_rss_mb()
res <- list(
  schema = 1L,
  host = list(cpus = 16L, cores = 8L, model = "AMD Ryzen 7 3700X", mem_total_mb = 32806),
  pinned = host$pinned, host_state = host_record(host), seed = SEED,
  design = list(k_groups = K_GROUPS, n_lr = N_LR, n_genes = N_GENES, nboot = NBOOT,
                repeats = REPEATS, warmup = WARMUP, headroom_factor = HEADROOM_FACTOR),
  grid = results, peak_rss_mb = rss, mem_available_mb = mem_available_mb()
)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
writeLines(j_val(res), OUT)
cat(sprintf("SYNTH -> %s (rss=%.0fMB)\n", OUT, rss))
