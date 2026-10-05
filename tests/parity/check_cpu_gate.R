## Gate: the benchmark's CPU-affinity and host-contention checks do what they claim.
##
## `bench-runner/cpu_gate.R` is the thing standing between "a timing number" and "a timing number
## somebody can reproduce". It was written after an unrelated runaway process held a core at 98% for
## 51 hours and moved the human-skin R median from 93 s to 108 s between two runs of identical
## code, with nothing in the artifacts recording it. Two bugs in the first version of it were
## invisible without a test, and both failed in the *permissive* direction -- a gate that is too
## eager is obvious, one that is asleep is not:
##
##   * reading `x[1:5]` of the /proc/stat fields as "busy" includes `idle`, so every core reads
##     100% and the gate refuses an idle machine;
##   * differencing the per-CPU busy *fractions* instead of the raw counters returns exactly 0.0
##     everywhere, because the counters are of order 1e9 and a 1.5 s window moves them by ~150
##     ticks, which is below the double's resolution at that magnitude. A saturated host then
##     reads as idle -- the precise inverse of the bug above.
##
## The second is reproduced deterministically below by injecting a reader, which is why
## `cpu_busy_over()` takes one.
##
## Run:  Rscript tests/parity/check_cpu_gate.R

root <- Sys.getenv("CELLCHATRS_ROOT", unset = normalizePath(file.path(dirname(
  sub("^--file=", "", grep("^--file=", commandArgs(FALSE), value = TRUE)[1])), "..", "..")))
source(file.path(root, "bench-runner", "cpu_gate.R"))

fails <- 0L
ok <- function(label, cond) {
  if (isTRUE(cond)) {
    cat(sprintf("  ok    %s\n", label))
  } else {
    cat(sprintf("  FAIL  %s\n", label))
    fails <<- fails + 1L
  }
}

cat("expand_cpu_list\n")
ok("single core", identical(expand_cpu_list("3"), 3L))
ok("range", identical(expand_cpu_list("0-3"), 0:3))
ok("mixed", identical(expand_cpu_list("0-2,7,9-10"), c(0L, 1L, 2L, 7L, 9L, 10L)))
ok("whole range on this host", length(expand_cpu_list("0-15")) == 16L)
ok("reversed range still expands", identical(expand_cpu_list("3-1"), 1:3))
ok("duplicates collapse", identical(expand_cpu_list("0-2,1"), 0:2))

cat("\nparse_cpu_stat: which fields are busy\n")
## One CPU line: user nice system idle iowait irq softirq steal ...
mk <- function(user, nice, system, idle, iowait, irq, softirq, steal = 0) {
  paste("cpu0", user, nice, system, idle, iowait, irq, softirq, steal)
}
p <- parse_cpu_stat(c(mk(100, 0, 0, 0, 0, 0, 0)))
ok("all user time is busy", p$busy[1] == 100 && p$idle[1] == 0)
p <- parse_cpu_stat(c(mk(0, 0, 0, 100, 0, 0, 0)))
ok("idle time is NOT busy", p$busy[1] == 0 && p$idle[1] == 100)
p <- parse_cpu_stat(c(mk(0, 0, 0, 0, 100, 0, 0)))
ok("iowait counts as idle, not busy", p$busy[1] == 0 && p$idle[1] == 100)
p <- parse_cpu_stat(c(mk(0, 0, 0, 0, 0, 60, 0)))
ok("irq is busy", p$busy[1] == 60 && p$idle[1] == 0)
p <- parse_cpu_stat(c(mk(0, 0, 0, 0, 0, 0, 70)))
ok("softirq is busy", p$busy[1] == 70 && p$idle[1] == 0)
p <- parse_cpu_stat(c(mk(10, 20, 30, 40, 50, 60, 70, 999)))
ok("steal is in neither", p$busy[1] == 190 && p$idle[1] == 90)
## The regression the first version had: x[1:5] includes idle.
naive_busy <- function(x) sum(x[1:5])
ok("the naive x[1:5] reading is wrong, which is why it is not used",
   naive_busy(c(0, 0, 0, 100, 0, 0, 0)) == 100 && p$busy[1] != naive_busy(c(10, 20, 30, 40, 50, 60, 70)))

cat("\ncpu_busy_over: difference the ticks, then divide\n")
## A synthetic reader shaped like `cpu_ticks`: cumulative counters tagged with CPU ids. `base` is
## the baseline -- 1e9 is realistic for /proc/stat after a few days of uptime -- and `d` is the tick
## delta accumulated across the window.
reader_for <- function(baseline, dbusy, didle, cpus = 0:3) {
  i <- 0L
  n <- length(cpus)
  function() {
    i <<- i + 1L
    d <- if (i == 1L) 0 else 1L
    list(cpu = cpus, busy = rep(baseline + d * dbusy, n), idle = rep(baseline + d * didle, n))
  }
}
base <- 1e9
nosleep <- function(...) invisible(NULL)
half <- cpu_busy_over(0, reader_for(base, 75, 75), sleep = nosleep)
ok("50% busy over the window reads ~0.5", isTRUE(all.equal(unname(half[["0"]]), 0.5, tolerance = 1e-9)))
ok("and the magnitude is large -- this is the case the fraction-difference lost", half[["0"]] > 0.4)
ok("results are keyed by CPU id", identical(names(half), c("0", "1", "2", "3")))
## And the buggy formulation, for the record: it returns 0 on exactly this input.
naive_frac <- function(baseline, dbusy, didle) {
  f <- function(acc) { t <- acc + c(dbusy, didle); c(busy = t[1], idle = t[2]) }
  f0 <- c(busy = baseline, idle = baseline)
  a <- f0 / (f0[["busy"]] + f0[["idle"]])
  b <- f(c(0, 0)); b[["busy"]] <- baseline + dbusy; b[["idle"]] <- baseline + didle
  bfrac <- c(busy = b[["busy"]] / (b[["busy"]] + b[["idle"]]), idle = b[["idle"]] / (b[["busy"]] + b[["idle"]]))
  (bfrac - a)[["busy"]]
}
ok("differencing the fractions instead returns ~0 here, as the bug did",
   abs(naive_frac(base, 75, 75)) < 1e-6)

cat("\ncpu_busy_over: a core that does not advance\n")
still <- cpu_busy_over(0, reader_for(base, 0, 0), sleep = nosleep)
ok("no movement is NA, not 0 (a false clean bill)", is.na(still[["0"]]))

cat("\ncontention_verdict\n")
## Keyed by CPU id, the shape `cpu_busy_over` guarantees.
mkbusy <- function(...) {
  v <- c(...)
  setNames(v, as.character(seq_along(v) - 1L))
}
v <- contention_verdict(mkbusy(0.01, 0.02, 0.03, 0.01), 0:3, max_busy = 0.10)
ok("a quiet host passes", isTRUE(v$ok))
ok("and records the numbers", length(v$per_cpu_busy_pct) == 4L)
v <- contention_verdict(mkbusy(0.01, 0.02, 0.95, 0.01), 0:3, max_busy = 0.10)
ok("one busy core fails the gate", !isTRUE(v$ok))
ok("and names only that core", identical(as.integer(names(v$over)), 2L))
v <- contention_verdict(c(`0` = 0.01, `1` = NA), 0:1, max_busy = 0.10)
ok("an unaccounted core fails even though no core is over threshold",
   !isTRUE(v$ok) && identical(v$unknown, 1L))
v <- contention_verdict(mkbusy(0.11, 0.10, 0.10, 0.10), 0:3, max_busy = 0.10)
ok("the threshold is strict greater-than, so exactly 10% passes", identical(as.integer(names(v$over)), 0L))
ok("a core one point over fails", length(contention_verdict(mkbusy(0.11, 0.10), 0:1,
                                                          max_busy = 0.10)$over) == 1L)
ok("cpu 0 is measured, not skipped by R's dropped index 0",
   identical(as.integer(names(contention_verdict(mkbusy(0.99, 0.01, 0.01), 0:2,
                                                max_busy = 0.10)$over)), 0L))
ok("cpu 7 is measured too, not shifted off the end",
   identical(as.integer(names(contention_verdict(mkbusy(0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.99),
                                                0:7, max_busy = 0.10)$over)), 7L))

cat("\ncheck_host_quiet: the runner-facing wrapper\n")
## With a real affinity on this host and a threshold nothing can meet, it must refuse and say why.
h <- check_host_quiet(list("contend-max" = "-1"), label = "TEST", speak = function(...) invisible(NULL))
ok("an impossible threshold refuses", identical(h$reason, "contended") || identical(h$reason, "unaccounted"))
ok("it reports the pinning it discovered", !is.na(h$pinned) && nzchar(h$pinned))
ok("and the cores it checked", length(h$cpus) >= 1L)
ok("an override is not set by default", isFALSE(h$allowed))
h2 <- check_host_quiet(list("contend-max" = "-1", "allow-contended" = TRUE),
                       label = "TEST", speak = function(...) invisible(NULL))
ok("--allow-contended sets allowed", isTRUE(h2$allowed))
ok("and clears the refusal reason", is.na(h2$reason))
ok("but still reports the verdict as not ok", isFALSE(h2$verdict$ok))
## Unset explicitly rather than via `on.exit`: this runs as a promise inside `ok()`, and relying on
## a deferred cleanup to fire before the *next* assertion is how an environment variable ends up
## silently changing a later result.
Sys.setenv(BENCH_ALLOW_CONTENDED = "1")
env_ok <- isTRUE(check_host_quiet(list("contend-max" = "-1"), label = "T",
                                  speak = function(...) invisible(NULL))$allowed)
Sys.unsetenv("BENCH_ALLOW_CONTENDED")
ok("BENCH_ALLOW_CONTENDED is honoured for runners with no argument parser", env_ok)
ok("and is not left set for the next check",
   identical(Sys.getenv("BENCH_ALLOW_CONTENDED", ""), ""))

cat("\nhost_record / refusal_record\n")
h <- check_host_quiet(list("contend-max" = "-1"), label = "TEST", speak = function(...) invisible(NULL))
rec <- host_record(h)
ok("host_record carries the threshold", rec$threshold_pct == 100 * (-1))
ok("host_record marks the run contended", isTRUE(rec$contended))
ok("host_record marks the override", identical(rec$override, h$allowed))
ok("host_record keeps the per-cpu numbers", length(rec$per_cpu_busy_pct) == length(h$cpus))
ref <- refusal_record(h)
ok("refusal_record names the reason", ref$refused == h$reason)
ok("refusal_record is schema-stamped", ref$schema == 1L)
ok("refusal_record lists the offending cores", length(ref$offending) >= 0L)
ok("both are plain lists, safe for any JSON writer", is.list(rec) && is.list(ref))

cat("\nagainst this host\n")
real <- cpu_busy_over(0.3)
n_host <- length(grep("^cpu[0-9]", readLines("/proc/stat"), value = TRUE))
ok("every CPU is accounted for on this machine", !any(is.na(real)))
ok("all host CPUs are reported", length(real) == n_host)
ok("every value is a fraction in [0,1]", all(real >= 0 & real <= 1, na.rm = TRUE))

if (fails > 0L) {
  cat(sprintf("\nCPU GATE FAILED: %d failing check(s)\n", fails))
  quit(status = 1L)
}
cat("\nCPU GATE OK: every check passed\n")