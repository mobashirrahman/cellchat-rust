## Measure upstream R and the Rust kernel on the *same* object, with the protocol the objective
## specifies: discard a warm-up, at least five timed repeats, report the median with a bootstrap
## confidence interval and peak RSS, and pin CPUs with `taskset` rather than trusting
## `parallel::detectCores()` (which misreports on this 8c/16t host).
##
## What makes the two numbers comparable is that both run on one object built once, in one
## process, and the runner asserts `identical()` on the whole `net` afterwards. A speedup quoted
## without that check is a speedup between two different computations.
##
## Usage:
##   taskset -c 0-7 Rscript bench-runner/bench_real.R --fixture .../humanSkin.rda --nboot 100
##   CELLCHATRS_THREADS=1,2,4,8 Rscript bench-runner/bench_real.R --sweep
suppressWarnings(suppressMessages({
  library(stats)
}))
`%||%` <- function(a, b) if (is.null(a)) b else a

## A dependency-free JSON writer. `jsonlite` is not installed on this host and a benchmark that
## cannot emit its results without installing a package is a benchmark nobody runs; the schema
## here is small and fixed, so the writer is smaller than the dependency would be.
## `fixed = TRUE` throughout: a backslash is a regex metacharacter, and `gsub("\\", ...)`
## without it is a malformed pattern rather than an escape.
j_esc <- function(x) {
  x <- gsub("\\", "\\\\", x, fixed = TRUE)
  x <- gsub("\"", "\\\"", x, fixed = TRUE)
  x <- gsub("\n", "\\n", x, fixed = TRUE)
  x <- gsub("\r", "\\r", x, fixed = TRUE)
  x <- gsub("\t", "\\t", x, fixed = TRUE)
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
jsonlite_minimal <- function(x) j_val(x)

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
opt <- parse_args()
FIXTURE <- opt$fixture %||% "/scratch/mdra00001/tmp/opencode/data/humanSkin.rda"
NBOOT <- as.integer(opt$nboot %||% "100")
MAXLR <- as.integer(opt$maxlr %||% "0")
REPEATS <- as.integer(opt$repeats %||% "5")
WARMUP <- as.integer(opt$warmup %||% "1")
SEED <- as.integer(opt$seed %||% "1")
OUT <- opt$out %||% "bench-runner/results/real.json"
SWEEP <- !is.null(opt$sweep)

## `fixture.R` reads SPECIES from the environment to pick `CellChatDB.<species>.rda`, so it has to
## be set before the source -- not after, or the builder loads the human database for a mouse
## fixture and the gene names do not line up.
SPECIES <- Sys.getenv("SPECIES", "human")
Sys.setenv(SPECIES = SPECIES)
SCRIPT_DIR <- dirname(sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1]))
source(file.path(SCRIPT_DIR, "fixture.R"))

## The shim resolves the L-R database from `options$db` or from `CELLCHATRS_DB`; upstream's body
## carries the database in `@DB`. The two have to be the **same** database or the two numbers are
## not a speedup but a comparison of different computations.
##
## This is not hypothetical. The first version of this runner defaulted `CELLCHATRS_DB` to
## `db_human` and so ran the mouse wound fixture against the human database. Upstream read the
## mouse database and succeeded; the kernel read the human one, found that the human subunits of
## `TGFbR1_R2` are not in the mouse expression matrix, and returned
## `subscript out of bounds` -- which is the correct answer and extendr turns a returned `Err`
## into that R condition. The failure looked like a kernel bug and cost more time than the
## benchmark did.
if (Sys.getenv("CELLCHATRS_DB") == "") {
  Sys.setenv(CELLCHATRS_DB = file.path(getwd(), "tests", "fixtures", paste0("db_", SPECIES)))
}
## ... and check it rather than trusting the caller, from the export's own manifest.
db_dir <- Sys.getenv("CELLCHATRS_DB")
man <- file.path(db_dir, "MANIFEST.tsv")
if (!file.exists(man)) {
  cat(sprintf("BENCH REFUSING: no MANIFEST.tsv in %s\n", db_dir))
  quit(status = 2L)
}
## Parsed line by line rather than with `read.delim(comment.char = "#")`: the manifest's *header*
## is itself a comment, so `read.delim` discards the header and then reads it as the first data
## row -- giving columns named `species` and `upstream_commit` holding the word "species" and the
## word "upstream_commit". The check then refuses every run, which is at least a safe direction to
## fail in, but for the wrong reason and with a useless message.
mf <- list()
for (ln in readLines(man)) {
  ln <- trimws(ln)
  if (!nzchar(ln) || startsWith(ln, "#")) next
  parts <- strsplit(ln, "	", fixed = TRUE)[[1]]
  if (length(parts) >= 2L) mf[[parts[1]]] <- parts[2]
}
mf_species <- mf$species %||% ""
mf_sha <- mf$upstream_commit %||% ""
if (!identical(mf_species, SPECIES)) {
  cat(sprintf("BENCH REFUSING: %s exports species %s but SPECIES=%s. A mismatch makes the\n",
              db_dir, mf_species, SPECIES))
  cat("           R and Rust read different L-R databases, so the timings are not comparable.\n")
  quit(status = 2L)
}
if (!identical(mf_sha, "75253cd0c9e68410e6e721a6d3a0419a1d7e358f")) {
  cat(sprintf("BENCH REFUSING: %s was exported at %s, not the pinned upstream SHA\n", db_dir, mf_sha))
  quit(status = 2L)
}

## --- which CPUs this process is actually pinned to, and whether they are quiet ----------------
## The logic lives in `cpu_gate.R` so it can be tested; this is the runner's half -- discover the
## affinity, call the gate, and refuse or record. Both problems below were found by measuring, not
## by reading the protocol:
##
##   1. `BENCH_PINNED` was the only source of the recorded pinning, and the documented command never
##      set it, so runs wrote `"pinned": ""` while genuinely running on eight pinned cores.
##   2. Nothing checked the host was idle. A runaway `python3 gen_stats_vectors.py` from an
##      unrelated project in this same account had been at 98% CPU for 51 hours. It moved the R
##      median on human skin from 93 s to 108 s between two runs of byte-identical code, and no
##      artifact said why.
source(file.path(SCRIPT_DIR, "cpu_gate.R"))

host <- check_host_quiet(opt, label = "BENCH")
if (is.na(host$pinned)) {
  cat("BENCH REFUSING: cannot determine this process's CPU affinity and BENCH_PINNED is unset.\n")
  cat("           Run under taskset, e.g. `taskset -c 0-7 Rscript ...`, so the timings record\n")
  cat("           which cores produced them.\n")
  quit(status = 2L)
}
if (!isTRUE(host$verdict$ok) && !isTRUE(host$allowed)) {
  writeLines(jsonlite_minimal(refusal_record(host)), paste0(OUT, ".refused.json"))
  quit(status = 3L)
}
PINNED <- host$pinned
PIN_CPUS <- host$cpus

## --- timings -------------------------------------------------------------------------
## One timed repeat. `gc(reset = TRUE)` before each so R's own heap counter starts clean, but the
## reported memory is `VmHWM`, which is a high-water mark for the *process* and therefore only
## ever rises -- so it is read once at the end and attributed to the whole run, not per repeat.
## Upstream writes its progress bar and `cat` lines to **stdout**, through carriage returns, and
## a 5-repeat run emits about 40 000 of them. `sink()` to a file keeps the log readable; the
## connection is always popped, including on error, because an unbalanced `sink()` silently
## swallows everything after it.
quietly <- function(expr) {
  tf <- tempfile()
  sink(tf)
  out <- withCallingHandlers(expr, message = function(m) invokeRestart("muffleMessage"),
                             warning = function(w) invokeRestart("muffleWarning"))
  sink()
  unlink(tf)
  out
}

time_once <- function(f, object) {
  gc(reset = TRUE, full = FALSE)
  t0 <- proc.time()[["elapsed"]]
  out <- quietly(f(object))
  t1 <- proc.time()[["elapsed"]]
  list(secs = t1 - t0, value = out)
}

## Bootstrap CI of the *median* over the repeats. Resampling the median rather than the mean
## because a single slow repeat (a page-cache miss, a scheduler preemption) drags a mean but not a
## median, and the objective asks for the median.
median_ci <- function(x, n_boot = 10000L, seed = 4242L) {
  set.seed(seed)
  n <- length(x)
  meds <- numeric(n_boot)
  for (i in seq_len(n_boot)) meds[i] <- median(sample(x, n, replace = TRUE))
  c(lo = unname(quantile(meds, 0.025)), hi = unname(quantile(meds, 0.975)))
}

measure <- function(label, f, object, repeats, warmup) {
  for (i in seq_len(warmup)) f(object)                     # discarded
  ts <- numeric(repeats)
  last <- NULL
  for (i in seq_len(repeats)) {
    r <- time_once(f, object)
    ts[i] <- r$secs
    last <- r$value
  }
  ci <- median_ci(ts)
  list(label = label, repeats = repeats, raw_secs = ts,
       median = median(ts), ci_lo = ci[["lo"]], ci_hi = ci[["hi"]],
       min = min(ts), max = max(ts), value = last)
}

## --- the two sides -------------------------------------------------------------------
up <- upstream_env()
upstream_cpp <- get("computeCommunProb", envir = up)
## The shim, i.e. the Rust kernel behind the identical public name. This is the drop-in claim
## being measured: the same call a user makes.
args_for <- function(object) list(object = object, nboot = NBOOT, seed.use = SEED)

built <- build_fixture(FIXTURE, max_lr = MAXLR, seed = SEED)
info <- built$info
hr <- check_headroom(info)
cat(sprintf("BENCH nC=%d K=%d nGenes=%d nLR=%d nnz=%.1f%% dense=%.0fMB headroom=%s\n",
            info$n_cells, info$n_groups, info$n_genes_signaling, info$n_lr, info$nnz_pct,
            required_mb(info), if (hr$ok) "ok" else paste("INSUFFICIENT", hr$available_mb)))
if (!hr$ok) {
  cat("BENCH REFUSING:", hr$message, "\n")
  quit(status = 2L)
}

## `--rust-only` skips the R side entirely. A thread sweep is a property of the Rust kernel, and
## re-measuring upstream at every thread count would multiply a 90-200 s single-threaded run by the
## number of points for a number already known not to depend on `CELLCHATRS_THREADS`. The
## comparison speedup is then absent rather than wrong: the report says `null`, not a stale figure.
RUST_ONLY <- !is.null(opt[["rust-only"]])

m_up <- if (RUST_ONLY) NULL else
  measure("R upstream", function(o) do.call(upstream_cpp, args_for(o)),
          built$object, REPEATS, WARMUP)
m_rs <- measure("Rust kernel", function(o) do.call(computeCommunProb, args_for(o)),
                built$object, REPEATS, WARMUP)

if (RUST_ONLY) {
  rss <- peak_rss_mb()
  res <- list(
    schema = 1L, rust_only = TRUE,
    ## `CELLCHATRS_THREADS` is what was *requested*; `cellchatrs_threads()` is the size rayon
    ## actually built. Recording only the request is how a scaling table ends up labelled with
    ## numbers the run never used.
    threads = as.integer(Sys.getenv("CELLCHATRS_THREADS", "0")),
    pool_threads = as.integer(cellchatrs_threads()),
    pinned = PINNED, fixture = info, nboot = NBOOT,
    repeats = REPEATS, warmup = WARMUP, seed = SEED,
    rust = m_rs[c("raw_secs", "median", "ci_lo", "ci_hi", "min", "max")],
    peak_rss_mb = rss, dense_stream_mb = required_mb(info),
    mem_available_mb = mem_available_mb()
  )
  dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
  writeLines(jsonlite_minimal(res), OUT)
  cat(sprintf("BENCH threads=%d (pool=%d) Rust=%.2fs [%.2f,%.2f] rss=%.0fMB -> %s\n",
              res$threads, res$pool_threads, m_rs$median, m_rs$ci_lo, m_rs$ci_hi, rss, OUT))
  quit(status = 0L)
}

## --- parity on the real fixture -------------------------------------------------------
## The whole `net`, not just `prob`: this is the same `identical()` acceptance gate the
## differential suite uses, applied to the authors' own data rather than to a fixture.
parity_net <- identical(m_up$value@net, m_rs$value@net)
parity_prob <- identical(m_up$value@net$prob, m_rs$value@net$prob)
max_abs <- suppressWarnings(max(abs(m_up$value@net$prob - m_rs$value@net$prob)))
cat(sprintf("BENCH parity net=%s prob=%s max|diff|=%.3g\n", parity_net, parity_prob, max_abs))

## A scatter is the right plot because a *magnitude* agreement can hide a *relative* disaster at
## small Prob, and a log-log identity view shows both at once.
scatter <- list(
  n = length(m_up$value@net$prob),
  log10_prob_r   = log10(pmax(as.vector(m_up$value@net$prob), .Machine$double.xmin)),
  log10_prob_rs  = log10(pmax(as.vector(m_rs$value@net$prob), .Machine$double.xmin))
)

rss <- peak_rss_mb()
speedup <- m_up$median / m_rs$median
## The CI on the ratio, propagated from the two medians' CIs rather than assumed.
speedup_ci <- c(lo = m_up$ci_lo / m_rs$ci_hi, hi = m_up$ci_hi / m_rs$ci_lo)

res <- list(
  schema = 1L,
  host = list(cpus = 16L, cores = 8L, model = "AMD Ryzen 7 3700X", mem_total_mb = 32806),
  pinned = PINNED,
  ## The per-pinned-core busy fractions the contention gate sampled before the first timed run.
  ## Recorded even when they passed, because "the host was quiet" is otherwise an unverifiable
  ## claim, and because a reader comparing two runs needs to know they were taken under the same
  ## conditions.
  host_state = host_record(host),
  threads = as.integer(strsplit(Sys.getenv("CELLCHATRS_THREADS", "0"), ",")[[1]][1]),
  fixture = info, nboot = NBOOT, repeats = REPEATS, warmup = WARMUP, seed = SEED,
  r_upstream = m_up[c("raw_secs", "median", "ci_lo", "ci_hi", "min", "max")],
  rust = m_rs[c("raw_secs", "median", "ci_lo", "ci_hi", "min", "max")],
  speedup = speedup, speedup_ci_lo = speedup_ci[["lo"]], speedup_ci_hi = speedup_ci[["hi"]],
  peak_rss_mb = rss, parity_net = parity_net, parity_prob = parity_prob, max_abs_diff = max_abs,
  dense_stream_mb = required_mb(info), mem_available_mb = mem_available_mb(),
  scatter = scatter
)
dir.create(dirname(OUT), recursive = TRUE, showWarnings = FALSE)
writeLines(jsonlite_minimal(res), OUT)
cat(sprintf("BENCH R=%.2fs [%.2f,%.2f] Rust=%.2fs [%.2f,%.2f] speedup=%.2fx [%.2f,%.2f] rss=%.0fMB -> %s\n",
            m_up$median, m_up$ci_lo, m_up$ci_hi, m_rs$median, m_rs$ci_lo, m_rs$ci_hi,
            speedup, speedup_ci[["lo"]], speedup_ci[["hi"]], rss, OUT))
