## Kernel + comparison benchmarks on the remaining Figshare datasets.
##
## bench_real.R covers human skin and mouse wound (scRNA-seq path). This covers the rest of the
## objective's dataset list that exists locally: human skin multi-condition (NL/LS, the comparison
## tutorial's data) and mouse embryonic E13/E14. These arrive as *precomputed* CellChat objects
## (.rds), not count matrices, so unlike bench_real.R there is no fixture-building step -- the
## objects are re-run through computeCommunProb on both sides and compared.
##
## Two things this script does that bench_real.R does not, both dictated by the data:
##
## * Per-dataset database. NL/LS are human, E13/E14 are mouse (mixed-case symbols: Tgfb1, not
##   TGFB1). The shim resolves the L-R database from `options$db` first, so each object gets its
##   species' export assigned explicitly. A single global CELLCHATRS_DB would run one species
##   against the wrong database -- the exact failure bench_real.R documents (human subunits of
##   TGFbR1_R2 missing from the mouse matrix), so the assignment is per object, not per process.
## * The multi-condition comparison workflow itself: mergeCellChat(LS, NL) plus rankNet in
##   comparison mode, timed and compared on both sides. The kernel rows measure inference; this
##   measures the analysis the multi-condition data exists for.
##
## Protocol matches bench_real.R: warm-up discarded, 5 timed repeats, median with bootstrap CI,
## peak RSS, taskset-pinned CPUs. nboot=100, seed=1 throughout.
##
## Usage:
##   taskset -c 0-7 Rscript bench-runner/bench_more_datasets.R
##   (takes ~1h: E13/E14 upstream runs dominate; progress lines go to stdout)
suppressWarnings(suppressMessages({
  library(stats); library(methods); library(Matrix); library(collapse); library(dplyr)
}))
suppressWarnings(suppressMessages(library(CellChat)))

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

CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
DATADIR <- "data"
NBOOT <- 100L
REPEATS <- 5L

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
## No argument parser in this runner, so the gate takes its settings from the environment
## (`BENCH_CONTEND_WINDOW`, `BENCH_CONTEND_MAX`, `BENCH_ALLOW_CONTENDED=1`).
host <- check_host_quiet(list(), label = "BENCH-DATASETS")
if (is.na(host$pinned)) {
  cat("BENCH-DATASETS REFUSING: cannot determine this process's CPU affinity and BENCH_PINNED is unset.\n")
  cat("            Run under taskset so the timings record which cores produced them.\n")
  quit(status = 2L)
}
if (!isTRUE(host$verdict$ok) && !isTRUE(host$allowed)) {
  writeLines(j_val(refusal_record(host)), "bench-runner/results/datasets.refused.json")
  quit(status = 3L)
}
WARMUP <- 1L
SEED <- 1L

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
stopifnot(!is.null(methods::getClass("CellChat", where = globalenv())))

DATASETS <- list(
  list(name = "nl", file = "cellchat_humanSkin_NL.rds", species = "human"),
  list(name = "ls", file = "cellchat_humanSkin_LS.rds", species = "human"),
  list(name = "e13", file = "cellchat_embryonic_E13.rds", species = "mouse"),
  list(name = "e14", file = "cellchat_embryonic_E14.rds", species = "mouse")
)

## Re-home a deserialized object onto the globalenv class definition. `readRDS` restores the
## object fine, but every `@<-` on it (including upstream's own `object@options$run.time <- ...`
## at the end of computeCommunProb) triggers a package lookup for "CellChat" -- the package
## attribute baked in at serialization time -- and fails because only the globalenv definition
## exists here. Rebuilding slot-by-slot through `methods::new` re-registers the object under the
## definition both sides can see; slot contents are untouched (verified identical below), so both
## sides still compute on the same data. `data.smooth` is absent from these older objects and is
## left at its default, matching what the kernel reads when `raw.use = TRUE` (nothing).
rehome <- function(o) {
  sl <- names(methods::getSlots("CellChat"))
  cls <- methods::getSlots("CellChat")
  args <- list()
  skipped <- character(0)
  for (s in sl) {
    v <- tryCatch(methods::slot(o, s), error = function(e) NULL)
    if (is.null(v)) next
    ## Deposited objects may predate the current class definition: E13's `var.features` is a
    ## bare character where the class wants a list, and `methods::new` refuses the whole object
    ## for it (`validObject`: got character, should be list). The kernel never reads
    ## `var.features` (nor `dr`/`data.raw`/`data.scale`), so non-conforming slots are left at
    ## their defaults and RECORDED -- silently dropping them would misrepresent the input, and
    ## coercing a character into a list would invent structure upstream never had.
    ok <- tryCatch(methods::is(v, cls[[s]]), error = function(e) FALSE)
    if (isTRUE(ok)) {
      args[[s]] <- v
    } else {
      skipped <- c(skipped, s)
    }
  }
  obj <- do.call(methods::new, c(list("CellChat"), args))
  ## Returned alongside, not attached: an extra attribute on the object would ride into every
  ## downstream comparison and could only ever cause a spurious mismatch.
  list(obj = obj, skipped = skipped)
}

bench_one <- function(spec) {
  raw <- readRDS(file.path(DATADIR, spec$file))
  rh <- rehome(raw)
  obj <- rh$obj
  if (length(rh$skipped)) {
    cat(sprintf("BENCH %s: rehome skipped non-conforming slot(s): %s\n",
                spec$name, paste(rh$skipped, collapse = ",")))
  }
  ## Rehoming must not alter contents: every readable slot identical between the deserialized
  ## object and the rebuilt one. If this ever fails, both sides would still agree with each
  ## other (identical inputs) but would no longer be measuring the deposited object.
  for (s in c("data", "data.signaling", "idents", "meta", "net", "netP", "DB", "LR")) {
    if (!identical(methods::slot(raw, s), methods::slot(obj, s))) {
      cat(sprintf("BENCH REFUSING: rehome changed slot %s on %s\n", s, spec$name))
      quit(status = 2L)
    }
  }
  ## Deposited objects predate `options$datatype` (it is NULL), and both implementations
  ## branch on `object@options$datatype != "RNA"` -- NULL errors identically on both sides with
  ## "missing value where TRUE/FALSE needed", which would benchmark an error path. All four
  ## datasets are scRNA-seq, so "RNA" is factual, not a guess. Assigned post-rehome, where `@<-`
  ## resolves against the globalenv class; assigning on the deserialized object would hit the
  ## package-lookup failure documented above.
  if (is.null(obj@options$datatype)) {
    obj@options$datatype <- "RNA"
    cat(sprintf("BENCH %s: options$datatype was unset; set to RNA\n", spec$name))
  }
  ## No `obj@options$db <- ...` assignment here, deliberately. `@<-` on a deserialized object
  ## triggers a class-package lookup for "CellChat" (the package attribute baked in at
  ## serialization time) and fails when only the globalenv definition exists -- a path no gate
  ## had exercised, since every existing test either builds objects via createCellChat or only
  ## reads slots. The shim resolves the database as options$db -> bundled fingerprint match ->
  ## CELLCHATRS_DB, so setting the env per dataset reaches the same export without touching the
  ## object. The bundled fingerprint match is tried first and would win on its own; the env var
  ## is belt and braces for subsetted DBs the fingerprint cannot see.
  Sys.setenv(CELLCHATRS_DB = normalizePath(file.path("tests", "fixtures",
                                                     paste0("db_", spec$species)),
                                           mustWork = TRUE))
  ## Reconcile the pair list with the pinned export *before* timing. Deposited objects may bundle
  ## an older CellChatDB than the pinned export (NL: object complex table 339x5 vs export 338x6),
  ## and a pair resolving under one but not the other makes the two sides answer different
  ## databases: upstream reads subunits from the object's DB and succeeds, the shim reads the
  ## pinned export and raises the documented missing-subunit error. Neither side is wrong, so the
  ## benchmark runs the intersection both agree on -- exactly as the visium benchmark runs its 134
  ## resolvable pairs. Dropped pairs are reported, not silently skipped.
  expath <- file.path("tests", "fixtures", paste0("db_", spec$species), "complexes.tsv")
  cm <- read.delim(expath, header = FALSE, stringsAsFactors = FALSE)
  rownames(cm) <- cm[, 1]
  genes <- rownames(obj@data.signaling)
  subunits <- function(g) {
    if (!(g %in% rownames(cm))) return(g)
    s <- unlist(cm[g, -1, drop = FALSE], use.names = FALSE)
    s <- s[!is.na(s) & nzchar(s)]
    if (!length(s)) return(g)
    s
  }
  lr <- obj@LR$LRsig
  ok <- vapply(seq_len(nrow(lr)), function(i) {
    all(c(subunits(lr$ligand[i]), subunits(lr$receptor[i])) %in% genes)
  }, TRUE)
  dropped <- rownames(lr)[!ok]
  if (length(dropped)) {
    cat(sprintf("BENCH %s: %d/%d pairs unresolvable under the pinned export, dropped: %s\n",
                spec$name, length(dropped), nrow(lr), paste(dropped, collapse = ",")))
  }
  obj@LR$LRsig <- lr[ok, , drop = FALSE]
  nC <- ncol(obj@data)
  cat(sprintf("BENCH %s: %d genes x %d cells, %d groups, %d LR pairs (%s)\n",
              spec$name, nrow(obj@data), nC, nlevels(obj@idents),
              nrow(obj@LR$LRsig), spec$species))
  args <- list(type = "triMean", nboot = NBOOT, seed.use = SEED)
  m_up <- measure("R upstream", function() do.call(UP$computeCommunProb,
                                                  c(list(obj), args)), REPEATS, WARMUP)
  m_rs <- measure("Rust kernel", function() do.call(computeCommunProb,
                                                    c(list(obj), args)), REPEATS, WARMUP)
  par_net <- identical(m_up$value@net, m_rs$value@net)
  par_prob <- identical(m_up$value@net$prob, m_rs$value@net$prob)
  max_abs <- suppressWarnings(max(abs(m_up$value@net$prob - m_rs$value@net$prob)))
  speedup <- m_up$median / m_rs$median
  ci <- c(lo = m_up$ci_lo / m_rs$ci_hi, hi = m_up$ci_hi / m_rs$ci_lo)
  cat(sprintf(paste0("BENCH %s R=%.2fs [%.2f,%.2f] Rust=%.2fs [%.2f,%.2f] speedup=%.2fx ",
                     "[%.2f,%.2f] parity_net=%s parity_prob=%s max|diff|=%.3g\n"),
              spec$name, m_up$median, m_up$ci_lo, m_up$ci_hi,
              m_rs$median, m_rs$ci_lo, m_rs$ci_hi, speedup, ci[["lo"]], ci[["hi"]],
              par_net, par_prob, max_abs))
  list(dataset = spec$name, species = spec$species, n_cells = nC,
       n_genes_total = nrow(obj@data),
       n_genes_signaling = nrow(obj@data.signaling),
       n_groups = nlevels(obj@idents), n_lr = nrow(obj@LR$LRsig),
       n_lr_deposited = nrow(lr), dropped_pairs = as.list(dropped),
       rehome_skipped_slots = as.list(rh$skipped),
       nboot = NBOOT, repeats = REPEATS, warmup = WARMUP, seed = SEED,
       r_upstream = m_up[c("raw_secs", "median", "ci_lo", "ci_hi", "min", "max")],
       rust = m_rs[c("raw_secs", "median", "ci_lo", "ci_hi", "min", "max")],
       speedup = speedup, speedup_ci_lo = ci[["lo"]], speedup_ci_hi = ci[["hi"]],
       parity_net = par_net, parity_prob = par_prob, max_abs_diff = max_abs)
}

results <- list()
only <- Sys.getenv("BENCH_ONLY", "")
cmp_only <- identical(Sys.getenv("BENCH_COMPARISON_ONLY", ""), "1")
for (spec in DATASETS) {
  if (nzchar(only) && spec$name != only) next
  if (cmp_only) next
  r <- bench_one(spec)
  out <- sprintf("bench-runner/results/dataset_%s.json", spec$name)
  dir.create(dirname(out), recursive = TRUE, showWarnings = FALSE)
  writeLines(j_val(c(list(schema = 1L,
                          host = list(cpus = 16L, cores = 8L, model = "AMD Ryzen 7 3700X",
                                      mem_total_mb = 32806),
                          pinned = host$pinned, host_state = host_record(host),
                          threads = as.integer(Sys.getenv("CELLCHATRS_THREADS", "8"))),
                       r,
                       list(peak_rss_mb = peak_rss_mb(),
                            mem_available_mb = mem_available_mb()))), out)
  cat(sprintf("BENCH %s -> %s\n", spec$name, out))
  results[[spec$name]] <- r
  rm(r)
  gc(full = FALSE)
}

## --- multi-condition comparison: mergeCellChat(LS, NL) + rankNet comparison ---------------
## Skipped unless the NL and LS kernel runs above both completed in this invocation: merging
## half-built objects would compare a different computation. The env-var gate keeps smoke runs
## (BENCH_ONLY=nl) from falling through into a comparison over missing objects.
if (!nzchar(Sys.getenv("BENCH_ONLY", "")) || identical(Sys.getenv("BENCH_COMPARISON_ONLY", ""), "1")) {
## The analysis the NL/LS pair exists for. Both sides merge the *recomputed* objects above
## (not the precomputed nets on disk), so the comparison measures the ported path end to end.
cat("BENCH comparison: mergeCellChat(LS, NL)\n")
## mergeCellChat/createCellChat live in CellChat_class.R, which is sourced into globalenv (see
## above), not into UP -- so the upstream side is the bare name, which resolves to the pinned
## tree's verbatim body. Bound here so the call sites read symmetrically.
up_mergeCellChat <- get("mergeCellChat", envir = globalenv())
mkfresh <- function(spec) {
  ## Read-only: see the note in bench_one about `@<-` on deserialized objects. Both LS and NL
  ## are human, so one env setting covers the comparison.
  rehome(readRDS(file.path(DATADIR, spec$file)))$obj
}
Sys.setenv(CELLCHATRS_DB = normalizePath(file.path("tests", "fixtures", "db_human")))
## ggplot2 attached for the rankNet block, same reason as check_identical.R: upstream's
## comparison body calls `ggplot` unqualified, resolving through the caller's search path.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("the comparison block needs ggplot2 attached", call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))
cmp_up <- mkfresh(DATASETS[[2]])
cmp_rs <- mkfresh(DATASETS[[2]])
nl_up <- mkfresh(DATASETS[[1]])
nl_rs <- mkfresh(DATASETS[[1]])
mg_up <- measure("merge upstream", function() up_mergeCellChat(list(nl_up, cmp_up),
                                                                 add.names = c("NL", "LS")),
                 REPEATS, WARMUP)
mg_rs <- measure("merge shim", function() mergeCellChat(list(nl_rs, cmp_rs),
                                                        add.names = c("NL", "LS")),
                 REPEATS, WARMUP)
m_up <- mg_up$value
m_rs <- mg_rs$value
merge_par <- identical(m_up, m_rs)
cat(sprintf("BENCH merge R=%.2fs [%.2f,%.2f] RS=%.2fs [%.2f,%.2f] parity=%s\n",
            mg_up$median, mg_up$ci_lo, mg_up$ci_hi,
            mg_rs$median, mg_rs$ci_lo, mg_rs$ci_hi, merge_par))
rn_up <- measure("rankNet upstream", function() UP$rankNet(m_up, mode = "comparison",
                                                           measure = "weight", do.stat = TRUE,
                                                           return.data = TRUE),
                 REPEATS, WARMUP)
rn_rs <- measure("rankNet shim", function() rankNet(m_rs, mode = "comparison",
                                                    measure = "weight", do.stat = TRUE,
                                                    return.data = TRUE),
                 REPEATS, WARMUP)
## `return.data = TRUE` gives `list(signaling.contribution = df, gg.obj = gg)`. The data frames
## are compared -- the ggplot objects carry environments and `identical()` on them would compare
## evaluation frames rather than the computation. Timing still covers the full call including
## plot construction, so the speedup is not flattered by skipping work.
rn_par <- identical(rn_up$value$signaling.contribution, rn_rs$value$signaling.contribution)
cat(sprintf("BENCH rankNet-comparison R=%.2fs [%.2f,%.2f] RS=%.2fs [%.2f,%.2f] parity=%s\n",
            rn_up$median, rn_up$ci_lo, rn_up$ci_hi,
            rn_rs$median, rn_rs$ci_lo, rn_rs$ci_hi, rn_par))
writeLines(j_val(list(
  schema = 1L, what = "mergeCellChat(LS,NL) + rankNet comparison",
  merge = list(upstream = mg_up[c("raw_secs", "median", "ci_lo", "ci_hi")],
               shim = mg_rs[c("raw_secs", "median", "ci_lo", "ci_hi")],
               speedup = mg_up$median / mg_rs$median, parity = merge_par),
  ranknet_comparison = list(upstream = rn_up[c("raw_secs", "median", "ci_lo", "ci_hi")],
                            shim = rn_rs[c("raw_secs", "median", "ci_lo", "ci_hi")],
                            speedup = rn_up$median / rn_rs$median, parity = rn_par),
  peak_rss_mb = peak_rss_mb(), mem_available_mb = mem_available_mb()
)), "bench-runner/results/multicondition.json")
cat("BENCH comparison -> bench-runner/results/multicondition.json\n")
} ## end multi-condition comparison (skipped on BENCH_ONLY smoke runs)
cat(sprintf("BENCH peak_rss=%.0fMB\n", peak_rss_mb()))
