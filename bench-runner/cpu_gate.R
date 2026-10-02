## CPU-affinity discovery and host-contention gate for the benchmark runners.
##
## Sourced by `bench_real.R`, `bench_amdahl.R` and `bench_synth.R`, and tested by
## `tests/parity/check_cpu_gate.R`. It lives in its own file because the two bugs it exists to
## prevent were both invisible from inside a 400-line runner: a gate nobody can call is a gate
## nobody can test.
##
## `cpu_busy_over()` takes its reader and its sleep as arguments. That is not for mocking's sake --
## it is what lets the precision bug below be reproduced deterministically. A test cannot reliably
## saturate a real core and then assert a measurement of 50%; with an injected reader it can.

`%||%` <- function(a, b) if (is.null(a)) b else a

## --- affinity ---------------------------------------------------------------------------------
## Read from the kernel rather than from an environment variable. `BENCH_PINNED` was the original
## source and the documented command never set it, so published runs recorded `"pinned": ""` while
## genuinely running on eight pinned cores. A timing table whose pinning is unknown cannot be
## reproduced, and the kernel's answer is the one that is actually true.
affinity_of_pid <- function(pid = Sys.getpid()) {
  out <- suppressWarnings(try(system2("taskset", c("-pc", pid), stdout = TRUE, stderr = FALSE),
                              silent = TRUE))
  if (inherits(out, "try-error") || !length(out)) return(NA_character_)
  m <- regmatches(out, regexpr(": [0-9,\\-]+", out))
  if (!length(m)) return(NA_character_)
  trimws(sub("^: ", "", m))
}

## "0-3,8,10-11" -> 0 1 2 3 8 10 11
expand_cpu_list <- function(spec) {
  out <- integer(0)
  for (part in strsplit(spec, ",", fixed = TRUE)[[1]]) {
    if (!nzchar(part)) next
    if (grepl("-", part, fixed = TRUE)) {
      b <- as.integer(strsplit(part, "-", fixed = TRUE)[[1]])
      if (length(b) != 2L || anyNA(b)) stop("malformed CPU range: ", part)
      out <- c(out, seq.int(b[1], b[2]))
    } else {
      out <- c(out, as.integer(part))
    }
  }
  sort(unique(out))
}

## --- contention ------------------------------------------------------------------------------
## /proc/stat per-CPU fields, in order: user, nice, system, idle, iowait, irq, softirq, steal,
## guest, guest_nice. Being 1-indexed in R, busy is x[1]+x[2]+x[3]+x[6]+x[7] and idle is x[4]+x[5].
## iowait counts as idle: a benchmark blocked on its own page cache is not being competed with.
## steal is in neither, because on this host it is always zero and folding it in would make a
## containerised run look busy for reasons it cannot control.
##
## The CPU ids are parsed out of the line prefixes and kept. They are what callers index by, and
## they have to be: `/proc/stat` lists cpu0 first, so cpu *n* sits at vector position *n + 1*, and
## indexing a CPU number directly (`v[0:7]`) drops the 0 and shifts every other reading by one --
## silently measuring cpus 1..7 and never measuring cpu 0 at all.
parse_cpu_stat <- function(lines) {
  v <- lapply(strsplit(lines, " +"), function(p) as.numeric(p[-1]))
  ids <- as.integer(sub("^cpu", "", sub(" .*$", "", lines)))
  list(cpu = ids,
       busy = vapply(v, function(x) x[1] + x[2] + x[3] + x[6] + x[7], 0),
       idle = vapply(v, function(x) x[4] + x[5], 0))
}

cpu_ticks <- function() parse_cpu_stat(grep("^cpu[0-9]", readLines("/proc/stat"), value = TRUE))

## Busy fraction per CPU across `window` seconds, named by CPU id.
##
## Difference the raw ticks, *then* divide. Differencing the pre-divided fractions instead is the
## obvious mistake and it fails silently in the permissive direction: the counters are of order 1e9,
## a double resolves roughly 1e-7 of that as rounding, and a 1.5 s window at 100 Hz advances them
## by only about 150 ticks -- so the difference of two fractions is pure rounding noise and returns
## 0.0 on every core. A saturated host then reads as a perfectly idle one.
##
## A core whose counters did not advance at all over the window (offline, or outside this cgroup's
## accounting) gets NA rather than 0. Reporting it as idle would be a false clean bill, so the
## caller refuses instead.
cpu_busy_over <- function(window, reader = cpu_ticks, sleep = Sys.sleep) {
  a <- reader()
  sleep(window)
  b <- reader()
  db <- b$busy - a$busy
  di <- b$idle - a$idle
  tot <- db + di
  setNames(ifelse(tot <= 0, NA_real_, db / tot), as.character(b$cpu))
}

## The verdict for a set of pinned cores. Returns a list; the caller decides what to do about it,
## so that this stays a pure function and the message formatting lives with the runner.
contention_verdict <- function(busy_frac, cpus, max_busy = 0.10) {
  pct <- 100 * busy_frac[as.character(cpus)]
  names(pct) <- cpus
  unknown <- cpus[is.na(pct)]
  over <- pct[!is.na(pct) & pct > 100 * max_busy]
  list(per_cpu_busy_pct = pct, unknown = unknown, over = over,
       threshold_pct = 100 * max_busy,
       ok = length(unknown) == 0L && length(over) == 0L)
}

## --- the runner's half -------------------------------------------------------------------------
## One function so that all five runners enforce the same protocol rather than each carrying its own
## copy of it. It reports and returns; it does not quit, because the caller has to know the outcome
## in order to decide whether to record a refusal artifact and with what path.
##
## `label` is the prefix for the printed lines ("BENCH", "BENCH-SYNTH", ...), so that the output of a
## suite of runners in one log is attributable.
check_host_quiet <- function(opt = list(), label = "BENCH", speak = cat) {
  pinned <- affinity_of_pid()
  if (is.na(pinned) || !nzchar(pinned)) pinned <- Sys.getenv("BENCH_PINNED", "")
  if (!nzchar(pinned)) {
    return(list(pinned = NA_character_, cpus = integer(0), reason = "unpinned",
                window_secs = NA_real_, max_busy = NA_real_, allowed = FALSE,
                verdict = list(ok = FALSE, per_cpu_busy_pct = numeric(0), unknown = integer(0),
                               over = numeric(0), threshold_pct = NA_real_)))
  }
  cpus <- expand_cpu_list(pinned)
  window <- as.numeric(opt[["contend-window"]] %||% Sys.getenv("BENCH_CONTEND_WINDOW", "1.5"))
  max_busy <- as.numeric(opt[["contend-max"]] %||% "0.10")
  ## `bench_more_datasets.R` has no argument parser, so the override has to be reachable without
  ## one. `BENCH_ALLOW_CONTENDED=1` is therefore accepted everywhere, not just by that runner.
  allowed <- !is.null(opt[["allow-contended"]]) ||
    identical(Sys.getenv("BENCH_ALLOW_CONTENDED", ""), "1")
  v <- contention_verdict(cpu_busy_over(window), cpus, max_busy)
  speak(sprintf("%s affinity=%s window=%.1fs per-pinned-cpu busy%%=%s\n", label, pinned, window,
                paste(sprintf("%d:%s", cpus,
                              ifelse(is.na(v$per_cpu_busy_pct), "?",
                                     sprintf("%.1f", v$per_cpu_busy_pct))), collapse = " ")))
  if (length(v$unknown)) {
    speak(sprintf("%s REFUSING: no CPU accounting for pinned core(s) %s over %.1fs.\n",
                  label, paste(v$unknown, collapse = ","), window))
    speak("           A core the kernel is not accounting for cannot be shown to be idle.\n")
    return(list(pinned = pinned, cpus = cpus, reason = "unaccounted", window_secs = window,
                max_busy = max_busy, allowed = allowed, verdict = v))
  }
  if (length(v$over)) {
    speak(sprintf("%s CONTENDED: %s (threshold %.0f%% over %.1fs)\n", label,
                  paste(sprintf("cpu %d at %.1f%%", as.integer(names(v$over)), v$over),
                        collapse = ", "), v$threshold_pct, window))
    speak("           Another process is using the pinned cores, so these timings would not be\n")
    speak("           reproducible and a rerun would not agree with them. Find it with\n")
    speak("           `ps -eo pcpu,pid,args --sort=-pcpu`, wait for it, or pin the run to idle cores.\n")
    speak("           `--allow-contended` records the contention and proceeds anyway.\n")
    if (!allowed)
      return(list(pinned = pinned, cpus = cpus, reason = "contended", window_secs = window,
                  max_busy = max_busy, allowed = allowed, verdict = v))
    speak("           --allow-contended given: proceeding, and recording the contention above.\n")
  }
  list(pinned = pinned, cpus = cpus, reason = NA_character_, window_secs = window,
       max_busy = max_busy, allowed = allowed, verdict = v)
}

## The block every runner embeds in its result JSON, so a reader can tell under what conditions the
## timings were taken without having to trust that the run happened to be quiet.
host_record <- function(host) {
  list(window_secs = host$window_secs, threshold_pct = 100 * host$max_busy,
       per_cpu_busy_pct = as.list(host$verdict$per_cpu_busy_pct),
       contended = !isTRUE(host$verdict$ok), override = isTRUE(host$allowed))
}

## The refusal artifact, written when a runner declines to measure. Returns the JSON text so each
## runner can hand it to whatever writer it has -- the five runners deliberately do not share a
## JSON writer, because each is meant to run standalone with no dependencies at all.
refusal_record <- function(host) {
  list(schema = 1L, refused = host$reason, pinned = host$pinned,
       window_secs = host$window_secs, threshold_pct = 100 * host$max_busy,
       per_cpu_busy_pct = as.list(host$verdict$per_cpu_busy_pct),
       offending = as.list(host$verdict$over))
}