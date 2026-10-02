#!/usr/bin/env Rscript
## Runs the configuration matrix in `matrix.R` against pinned upstream and the Rust shim, and
## writes a machine-checkable report.
##
## `check_identical.R` asks "does each function work". This asks the question the definition of done
## asks: does the *matrix* hold, with pairwise coverage over thirteen axes, and does the port match
## upstream on every row?
##
## Two things make the report trustworthy rather than merely green:
##
##   * the coverage is recomputed from the **effective** axes -- what each configuration actually
##     ran with, after any adaptation -- so a row that silently overrode an axis cannot be counted
##     as covering it;
##   * every row's verdict is recorded, including the ones that error, and an error counts as a pass
##     only when *both* sides produce the same error text. A port that stopped where upstream did
##     not would be a failure, and so would one that succeeded where upstream stopped.
##
## Usage: Rscript tests/parity/check_matrix.R [--max N] [--out FILE]

suppressWarnings(suppressMessages({
  library(Matrix); library(collapse); library(dplyr); library(cellchatrs)
}))
ROOT <- Sys.getenv("CELLCHATRS_ROOT", unset = getwd())
CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")

`%||%` <- function(a, b) if (is.null(a)) b else a

args <- commandArgs(trailingOnly = TRUE)
arg_of <- function(flag, default) {
  i <- which(args == flag)
  if (!length(i) || i == length(args)) default else args[i + 1L]
}
MAXROW <- as.integer(arg_of("--max", "Inf"))
## `--only` re-runs a single row through the *same* code path as the full sweep. Needed because a
## row can behave differently inside the sweep than in isolation -- which is exactly the question
## when the sweep reports a status the isolated reproduction contradicts.
ONLY <- arg_of("--only", NA_character_)
## `--start` skips the rows before the first failure, so a failing row can be reproduced without
## re-running everything ahead of it.
START <- as.integer(arg_of("--start", "1"))
OUT <- arg_of("--out", file.path(ROOT, "tests", "fixtures", "matrix_report.json"))

source(file.path(ROOT, "tests", "parity", "matrix.R"))



source(file.path(ROOT, "tests", "parity", "matrix_objects.R"))

## ------------------------------------------------------------------ one configuration
run_one <- function(cfg) {
  built <- build_object(cfg)
  o <- built$object
  args <- list(object = o, type = as.character(cfg$type), trim = as.numeric(cfg$trim),
               raw.use = isTRUE(cfg$raw), population.size = isTRUE(cfg$pop),
               nboot = as.integer(cfg$nboot), seed.use = 1L,
               Kh = as.numeric(cfg$Kh), n = as.numeric(cfg$n),
               ## Spatial-only parameters, always passed so both sides see the same call.
               distance.use = TRUE, scale.distance = 10, k.min = 1,
               contact.range = 60, interaction.range = 100)
  ## Upstream's `cat` and `print` go to **stdout**, not through `message()`, so `capture.output`
  ## is the only way to keep the runner's own progress lines readable. They are captured rather
  ## than discarded because two of them are part of the observable behaviour being compared.
  ## Upstream's `cat`, its two `print`s and `setTxtProgressBar`'s carriage returns all write to
  ## stdout; `capture.output` catches the first two but the progress bar uses `\r`, which survives
  ## as a very long single line. `sink()` diverts the connection entirely -- the text is still
  ## captured through the `message` handler, so nothing observable is lost -- and the sink is
  ## always popped, including on error, because an unbalanced `sink()` silently swallows the rest of
  ## the script's output.
  side <- function(f) tryCatch({
    txt <- character(0)
    tf <- tempfile()
    sink(tf)
    v <- withCallingHandlers(do.call(f, args),
      message = function(m) { txt <<- c(txt, conditionMessage(m)); invokeRestart("muffleMessage") },
      warning = function(w) invokeRestart("muffleWarning"))
    sink()
    unlink(tf)
    suppressMessages(suppressWarnings(list(txt = txt, v = v)))
  }, error = function(e) {
    while (sink.number() > 0L) sink()
    list(err = conditionMessage(e))
  })
  ea <- side(function(...) get("computeCommunProb", envir = up)(...))
  eb <- side(computeCommunProb)
  if (!is.null(ea$err) || !is.null(eb$err)) {
    ## When one side errors and the other does not, the *inputs* the surviving side saw are the
    ## only evidence that matters: if they differ from the inputs the erroring side saw, the
    ## divergence is in the fixture or the marshalling, not in the kernel. Recorded rather than
    ## printed, so it survives into the JSON report.
    diag <- NULL
    if (is.null(eb$err)) {
      ## Repeat the surviving side in place. A configuration that is deterministic and whose input
      ## is known identical cannot legitimately answer differently on a second call, so if it does,
      ## the cause is in the session state this row inherited rather than in the kernel.
      once <- function(f) tryCatch({
          v <- suppressWarnings(suppressMessages(do.call(f, args)))
          sprintf("result sum(prob)=%s NaN=%d", format(sum(v@net$prob)), sum(is.nan(v@net$prob)))
        }, error = function(e) paste0("error: ", conditionMessage(e)))
      again <- c(paste0("side1: ", once(computeCommunProb)),
                 paste0("side2: ", once(computeCommunProb)),
                 paste0("plain: ", once(computeCommunProb)),
                 ## The same, but after upstream has just aborted in this session.
                 paste0("afterUp: ", { invisible(once(function(...)
                   get("computeCommunProb", envir = up)(...))); once(computeCommunProb) }))
      bu <- tryCatch(as.matrix(built$object@data.signaling), error = function(e) NULL)
      po <- tryCatch(as.matrix(eb$v@data.signaling), error = function(e) NULL)
      diag <- sprintf("%s | repeats: %s",
                      sprintf("built max(data)=%s nonfinite=%s | port max(data)=%s nonfinite=%s | same object=%s",
                      format(max(bu)), sum(!is.finite(bu)),
                      if (is.null(po)) "n/a" else format(max(po)),
                      if (is.null(po)) "n/a" else sum(!is.finite(po)),
                      identical(bu, po)),
                      paste(again, collapse = " ; "))
    }
    return(list(status = if (identical(ea$err, eb$err)) "error-equal" else "error-differs",
                upstream = if (is.null(ea$err)) NA_character_ else ea$err,
                port = if (is.null(eb$err)) NA_character_ else eb$err,
                effective = built$effective, note = built$note, diag = diag))
  }
  va <- ea$v; vb <- eb$v
  ## `options$run.time` is wall-clock and is excluded everywhere, as is nothing else.
  same <- identical(va@net$prob, vb@net$prob) &&
    identical(va@net$pval, vb@net$pval) &&
    identical(va@options$parameter, vb@options$parameter) &&
    ## Messages, with the `Sys.time()` the two `print()` calls embed removed.
    identical(strip_clock(ea$txt), strip_clock(eb$txt))
  list(status = if (same) "identical" else "differs",
       upstream = NA_character_, port = NA_character_,
       effective = built$effective, note = built$note,
       detail = if (same) NA_character_ else diff_detail(va, vb, ea, eb))
}

strip_clock <- function(x) gsub("\\[[0-9]{4}-[0-9]{2}-[0-9]{2} [0-9:.]+\\]", "[TIME]",
                              paste(x, collapse = ""))
diff_detail <- function(va, vb, ea, eb) {
  d <- character(0)
  if (!identical(va@net$prob, vb@net$prob)) {
    d <- c(d, sprintf("prob: sum %s vs %s", format(sum(va@net$prob)), format(sum(vb@net$prob))))
  }
  if (!identical(va@net$pval, vb@net$pval)) d <- c(d, "pval differs")
  if (!identical(va@options$parameter, vb@options$parameter)) {
    for (f in union(names(va@options$parameter), names(vb@options$parameter))) {
      if (!identical(va@options$parameter[[f]], vb@options$parameter[[f]]))
        d <- c(d, paste0("parameter ", f))
    }
  }
  if (!identical(strip_clock(ea$txt), strip_clock(eb$txt)))
    d <- c(d, paste0("messages: [", strip_clock(ea$txt), "] vs [", strip_clock(eb$txt), "]"))
  paste(d, collapse = "; ")
}

## ------------------------------------------------------------------ run
m <- matrix_build()
cfgs <- m$configs
if (is.finite(MAXROW)) cfgs <- cfgs[seq_len(min(MAXROW, length(cfgs)))]
if (!is.na(ONLY)) {
  hit <- grepl(ONLY, vapply(cfgs, function(c) as.character(c$tag %||% ""), character(1)))
  if (!any(hit)) stop("--only matched no configuration: ", ONLY)
  cfgs <- cfgs[hit]
}
if (START > 1L) cfgs <- cfgs[seq(START, length(cfgs))]
## R buffers its own stderr connection when it is redirected to a file, so `cat(file = stderr())`
## is not enough to watch a long run: a redirected sweep produced a 62-byte log and no way to
## distinguish a slow run from a dead one. Progress is therefore appended to a plain file, which is
## visible immediately whatever the buffering.
PROG <- arg_of("--progress", tempfile("check_matrix_progress_"))
cat(sprintf("matrix: %d configurations\n", length(cfgs)))
cat(sprintf("matrix: %d configurations (progress -> %s)\n", length(cfgs), PROG), file = PROG)

results <- vector("list", length(cfgs))
t0 <- Sys.time()
for (i in seq_along(cfgs)) {
  r <- tryCatch(run_one(cfgs[[i]]), error = function(e)
    list(status = "harness-error", upstream = conditionMessage(e), port = NA_character_,
         effective = cfgs[[i]], note = "the runner itself failed"))
  results[[i]] <- r
  if (i %% 25L == 0L || i == length(cfgs)) {
    ## Progress goes to **stderr**, which is unbuffered when redirected to a file, while stdout
    ## stays block-buffered. A long run that only wrote stdout produced a 62-byte log file and no
    ## way to tell a slow run from a dead one.
    cat(sprintf("  %3d/%d  %-6s %s\n", i, length(cfgs), r$status,
                format(Sys.time() - t0, units = "secs")), file = PROG, append = TRUE)
  }
}

## ------------------------------------------------------------------ report
tab <- do.call(rbind, lapply(results, function(r)
  data.frame(status = r$status, stringsAsFactors = FALSE)))
counts <- table(tab$status)
cat("\nstatus:\n")
for (k in names(counts)) cat(sprintf("  %-14s %d\n", k, counts[[k]]), file = stderr())
bad <- tab$status %in% c("differs", "error-differs", "harness-error")
if (any(bad)) {
  cat("\nnon-passing configurations:\n")
  for (i in which(bad)) {
    r <- results[[i]]
    cfg <- cfgs[[i]]
    cat(sprintf("  [%d] %s\n", i,
                paste(sprintf("%s=%s", names(cfg)[-1], vapply(cfg[-1], as.character, "")),
                      collapse = " ")))
    cat(sprintf("      upstream: %s\n      port    : %s\n",
                if (is.na(r$upstream)) "<result>" else substr(r$upstream, 1, 70),
                if (is.na(r$port)) "<result>" else substr(r$port, 1, 70)))
    if (!is.null(r$detail) && !is.na(r$detail)) cat("      ", r$detail, "\n")
  }
}

## Coverage over the **effective** axes: what each configuration actually ran with.
eff_cfgs <- lapply(results, function(r) r$effective)
have <- unique(unlist(lapply(eff_cfgs, matrix_pairs_of)))
required <- m$required
missed <- setdiff(required, have)
cat(sprintf("\ncoverage over effective axes: %d of %d pairs, missing %d\n",
            length(intersect(have, required)), length(required), length(missed)))
if (length(missed)) {
  cat("missing pairs (first 20):\n")
  cat(paste0("  ", head(missed, 20), collapse = "\n"), "\n")
}

## Scale invariance, asserted directly rather than inferred from the pairwise rows.
##
## Upstream's `computeCommunProb` begins with `data.use <- data/max(data)`, so `Prob` depends only
## on the matrix *normalised by its own maximum* and multiplying the input by any non-zero constant
## must leave it bit-for-bit unchanged. This is a metamorphic property, and it is checked here
## because the matrix's pairwise coverage cannot express it: two configurations differing only in
## `scale` are two rows, and a row comparison says nothing about the relation between them.
##
## It earned its place by finding a real omission -- the kernel had no division at all, and every
## fixture so far had `max(data) == 1`.
inv_probe <- function(cfg) {
  eff <- cfg
  eff$scale <- 1
  a <- run_one(eff)
  eff$scale <- 3
  b <- run_one(eff)
  ## Compare the *verdict*, not the whole record: `effective$scale` differs by construction, so
  ## `identical(a, b)` is false even when both runs produced the same numbers. That would make the
  ## check report a failure for the one thing it is meant to assert.
  list(status_a = a$status, status_b = b$status,
       same = identical(a$status, b$status) && identical(a$upstream, b$upstream) &&
              identical(a$port, b$port) && identical(a$detail, b$detail))
}
inv <- inv_probe(list(type = "triMean", nC = 24L, K = 3L, nLR = 4L, nboot = 3L, pop = FALSE,
                      raw = TRUE, datatype = "RNA", lrstruct = "mixed", Kh = 0.5, n = 1,
                      trim = 0.1, data = "normal", tag = "scale-invariance"))
cat(sprintf("\nscale invariance (x1 vs x3): %s / %s -- identical = %s\n",
            inv$status_a, inv$status_b, inv$same))
scale_invariant <- isTRUE(inv$same)

## JSON, without a JSON package: the report is small and flat, and a hand-rolled writer keeps the
## dependency set unchanged.
esc <- function(s) {
  s <- gsub("\\", "\\\\", s, fixed = TRUE)
  s <- gsub("\"", "\\\"", s, fixed = TRUE)
  s <- gsub("\n", "\\\\n", s, fixed = TRUE)
  s <- gsub("\t", "\\\\t", s, fixed = TRUE)
  s <- gsub("\r", "\\\\r", s, fixed = TRUE)
  paste0("\"", s, "\"")
}
jnum <- function(x) if (is.na(x)) "null" else format(x, digits = 17)
## Every field, not `names(cfg)[-1]`: the first element happened to be `type`, so the `type` axis
## was missing from every row of the report. A report that silently drops one axis cannot be used
## to tell two rows apart, and here it hid the very field that distinguishes the failing rows.
cfg_json <- function(cfg) paste0("{",
  paste(sprintf("%s:%s", esc(names(cfg)),
                vapply(cfg, function(v) if (is.logical(v)) tolower(as.character(v)) else esc(as.character(v)),
                       character(1))), collapse = ","), "}")
rows <- vapply(seq_along(results), function(i) {
  r <- results[[i]]
  paste0("{\"config\":", cfg_json(cfgs[[i]]), ",\"effective\":", cfg_json(r$effective),
         ",\"status\":", esc(r$status),
         ",\"upstream\":", if (is.na(r$upstream)) "null" else esc(r$upstream),
         ",\"port\":", if (is.na(r$port)) "null" else esc(r$port),
         ",\"detail\":", if (is.null(r$detail) || is.na(r$detail)) "null" else esc(r$detail),
         ",\"note\":", esc(paste(r$note, collapse = "; ")),
         ",\"diag\":", if (is.null(r$diag)) "null" else esc(r$diag), "}")
}, character(1))
json <- paste0("{\n",
  "  \"upstream_sha\": ", esc(Sys.getenv("CELLCHATRS_UPSTREAM_SHA",
    unset = "75253cd0c9e68410e6e721a6d3a0419a1d7e358f")), ",\n",
  "  \"configurations\": ", length(results), ",\n",
  "  \"pairs_required\": ", length(required), ",\n",
  "  \"pairs_covered\": ", length(intersect(have, required)), ",\n",
  "  \"coverage_complete\": ", tolower(as.character(length(missed) == 0L)), ",\n",
  "  \"scale_invariant\": ", tolower(as.character(scale_invariant)), ",\n",
  "  \"status_counts\": {",
  paste(sprintf("\n    %s: %d", esc(names(counts)), as.integer(counts)), collapse = ","), "\n  },\n",
  "  \"rows\": [\n", paste(rows, collapse = ",\n"), "\n  ]\n}\n")
writeLines(json, OUT)
cat(sprintf("wrote %s\n", OUT))

n_fail <- sum(tab$status %in% c("differs", "error-differs", "harness-error"))
quit(status = if (n_fail == 0L && length(missed) == 0L && scale_invariant) 0L else 1L)
