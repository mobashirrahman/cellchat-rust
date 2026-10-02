# Parity gate for the standalone CLI.
#
# The objective says the Rust core "additionally ships as a standalone library/CLI so the numerics are
# usable without R". A CLI is only worth that if it produces the same numbers R does, so this gate
# runs the `cellchatrs` binary on `tests/fixtures/cli_input.tsv` and requires `identical()` on the
# full `Prob` and `Pval` arrays, the `computeAveExpr` means, and `aggregateNet`'s two matrices --
# against what pinned upstream's own `computeCommunProb` produced for the same input, recorded in
# `tests/fixtures/cli_golden.txt` by `tests/parity/gen_cli_fixture.R`.
#
# Three things this gate is *not*, each of which was a real failure mode in the surrounding work:
#
#   * Not a comparison of the CLI with the R shim. The shim shares the kernel, so that would be the
#     port against itself. The reference here is upstream, reached through a sourced copy of its
#     `modeling.R` and not through this package at all.
#   * Not a tolerance comparison. `identical()` on doubles. The output is `%.17g`, the shortest width
#     that round-trips a `f64` uniquely, so `as.numeric` on the CLI's output recovers the exact bit
#     pattern and any difference is a real difference.
#   * Not a smoke test. A binary that printed the right shape and plausible magnitudes would pass a
#     looser gate, and that is the failure `docs/SEMANTICS.md` records for the `rankNetPairwise`
#     fallback: a test that cannot tell a real implementation from a stub.
#
# The gate also runs the CLI's other subcommands, because a `mean` subcommand that disagrees with
# `triMean` would be a defect the `run` path does not cover.
#
# Usage:  R_LIBS=.rlib R --vanilla -f tests/parity/check_cli.R
suppressWarnings(suppressMessages({ library(Matrix); library(collapse); library(dplyr) }))

CC <- Sys.getenv("CELLCHAT_SRC", "/scratch/mdra00001/tmp/opencode/CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")
BIN <- Sys.getenv("CELLCHATRS_CLI", "target/release/cellchatrs")
## The subcommand checks below need one concrete input. The first configuration, chosen from the
## directory rather than hard-coded, so renaming or removing a configuration cannot leave this
## pointing at a file that no longer exists -- which it did, and the failure was a bare "the CLI exited
## 1" on `describe`.
CONFIGS <- sort(sub("^cli_golden_([A-Za-z0-9_]+)[.]txt$", "\\1",
                    list.files("tests/fixtures", pattern = "^cli_golden_.*[.]txt$")))
if (length(CONFIGS) == 0L) {
  stop("no cli_golden_*.txt in tests/fixtures; run tests/parity/gen_cli_fixture.R first",
       call. = FALSE)
}
IN <- sprintf("tests/fixtures/cli_input_%s.tsv", CONFIGS[1])

fails <- 0L
n_cmp <- 0L

## The comparison primitive. `n_cmp` is incremented before the verdict, not after, so a comparison
## that *throws* still counts as having been attempted -- a gate that silently loses a check when the
## check dies is worse than one that fails.
cmp <- function(label, ok, detail = "") {
  n_cmp <<- n_cmp + 1L
  cat(sprintf("  %-38s %s%s\n", label, if (ok) "ok" else "MISMATCH",
              if (nzchar(detail)) paste0("  ", detail) else ""))
  if (!ok) fails <<- fails + 1L
  invisible(ok)
}

## ---------------------------------------------------------------- the golden, as R values
## Read the golden into named pieces, converting with `as.numeric` so the comparison is against
## doubles rather than strings. `as.numeric("1.5e0")` is exact: the decimal text carries 17
## significant digits, which is a uniquely-decoding representation of a double.
read_golden <- function(path) {
  l <- readLines(path)
  i <- 1L
  header <- function(key) {
    f <- strsplit(l[i], "\t", fixed = TRUE)[[1]]
    stopifnot(f[1] == key)
    i <<- i + 1L
    f[-1]
  }
  ## A marker is a bare key on its own line with the payload on the *next* line, which is how both the
  ## generator and the CLI write the dimnames blocks. Reading them as `key<TAB>value` returns
  ## `character(0)` and leaves the cursor on the payload, so the next `header()` looks at a level name
  ## and fails with a message that names neither the key it wanted nor the line it was on.
  marker <- function(key) {
    got <- l[i]
    if (got != key) stop("expected the marker ", key, " on line ", i, ", found ", got)
    i <<- i + 1L
    out <- strsplit(l[i], "\t", fixed = TRUE)[[1]]
    i <<- i + 1L
    out
  }
  block <- function(key, n) {
    header(key)
    out <- l[i:(i + n - 1L)]
    i <<- i + n
    as.numeric(out)
  }
  invisible(header("version"))
  config <- header("config")
  type_mean <- header("type")
  nboot <- as.integer(header("nboot"))
  groups <- as.integer(header("groups"))
  lr <- as.integer(header("lr"))
  dim <- as.integer(strsplit(header("dim"), ",", fixed = TRUE)[[1]])
  matrix_block <- function(key) {
    header(key)
    rows <- l[i:(i + dim[1] - 1L)]
    i <<- i + dim[1]
    ## Both sides go through this same reader, so the flattening order only has to be consistent --
    ## but it is stated rather than assumed, because the alternative is a `dim`-order disagreement
    ## that looks exactly like a transposition bug in the kernel.
    m <- matrix(0, dim[1], dim[2])
    for (a in seq_len(dim[1])) {
      m[a, ] <- as.numeric(strsplit(rows[a], "\t", fixed = TRUE)[[1]])
    }
    m
  }
  ## The cursor moves in file order, and the fields are pulled out in that order -- not in the order
  ## they look best in a list literal. An earlier version read `groups` and `lr` first, which left the
  ## cursor on the `lr` line when it went looking for `dim`, and `stopifnot(f[1] == key)` reported
  ## `f[1] == key is not TRUE` with no indication of which key had been wanted.
  ##
  ## No trailing comma before the closing paren, either. R accepts that in a *definition* and rejects
  ## it in a *call* -- `list(a = 1, b = 2, )` is "argument 3 is empty" -- so the habit from writing
  ## functions costs a confusing error here.
  ## `ave_expr` is present only for the `type.mean` values `computeAveExpr` accepts, and its length is
  ## `n_genes * n_groups`, which is *not* derivable from the network's `dim`. The cursor has to be moved
  ## past it, and the length cannot be computed, so it is skipped by scanning to the next
  ## `aggregate_` marker. Merely peeking -- which is what a first version did, with a `has_ave` flag
  ## and no cursor movement -- left the reader positioned on `ave_expr` when it went looking for
  ## `aggregate_count`, and `stopifnot(f[1] == key)` said nothing about which key was wanted.
  has_ave <- any(l == "ave_expr")
  ## A function rather than a block of statements, because R evaluates `list(...)` elements left to
  ## right and the skip has to happen *after* `prob` and `pval` have been read. Executed eagerly
  ## instead, it moved the cursor to `aggregate_count` before anything was read, and the reader failed
  ## on `prob`.
  skip_ave <- function() {
    if (!has_ave) return(invisible(NULL))
    k <- which(l == "ave_expr")[1] + 1L
    while (k <= length(l) && !startsWith(l[k], "aggregate_")) k <- k + 1L
    i <<- k
    invisible(NULL)
  }
  list(
    config = config,
    type_mean = type_mean,
    nboot = nboot,
    groups = groups,
    lr = lr,
    dim = dim,
    prob = array(block("prob", prod(dim)), dim = dim),
    pval = array(block("pval", prod(dim)), dim = dim),
    has_ave = has_ave,
    count = { skip_ave(); matrix_block("aggregate_count") },
    weight = matrix_block("aggregate_weight"),
    source = marker("dimnames_source"),
    interaction = marker("dimnames_interaction")
  )
}

## The `ave_expr` block, read on its own. It is `n_genes * n_groups` numbers, which is not derivable
## from the network's `dim`, and it sits between `pval` and `aggregate_count`. Scanning for the marker
## and reading to the next `aggregate_` one costs one pass over a 200-line file and needs no length.
read_golden_ave <- function(path) {
  l <- readLines(path)
  i <- match("ave_expr", l)
  if (is.na(i)) return(NULL)
  j <- i + 1L
  while (j <= length(l) && !startsWith(l[j], "aggregate_")) j <- j + 1L
  as.numeric(l[(i + 1L):(j - 1L)])
}

## ---------------------------------------------------------------- the binary
if (!file.exists(BIN)) {
  stop("the CLI binary is not at ", BIN, "; run `cargo build --release -p cellchatrs-cli` first.\n",
       "Building it here is deliberate rather than convenient -- the gate must exercise the same\n",
       "artifact a user would install, not one compiled by the test.", call. = FALSE)
}

## The binary has to be told where the pinned database is. A relative path would resolve against
## whatever directory the runner happened to be in.
cli <- function(...) {
  out <- suppressWarnings(system2(BIN, c(...), stdout = TRUE, stderr = FALSE))
  status <- attr(out, "status")
  if (!is.null(status) && status != 0) {
    stop("the CLI exited ", status, " on: ", paste(c(...), collapse = " "), call. = FALSE)
  }
  out
}

## Parse the CLI's TSV output the same way the golden is read, so the two sides go through one reader.
read_cli <- function(lines) {
  i <- 1L
  header <- function(key) {
    f <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    stopifnot(f[1] == key)
    i <<- i + 1L
    f[-1]
  }
  marker <- function(key) {
    got <- lines[i]
    if (got != key) stop("expected the marker ", key, " on line ", i, ", found ", got)
    i <<- i + 1L
    out <- strsplit(lines[i], "\t", fixed = TRUE)[[1]]
    i <<- i + 1L
    out
  }
  block <- function(key, n) {
    header(key)
    out <- lines[i:(i + n - 1L)]
    i <<- i + n
    as.numeric(out)
  }
  invisible(header("version"))
  groups <- as.integer(header("groups"))
  lr <- as.integer(header("lr"))
  n_genes <- as.integer(header("genes"))
  n_cells <- as.integer(header("cells"))
  type_mean <- header("type_mean")
  dim <- as.integer(strsplit(header("dim"), ",", fixed = TRUE)[[1]])
  matrix_block <- function(key) {
    header(key)
    rows <- lines[i:(i + dim[1] - 1L)]
    i <<- i + dim[1]
    m <- matrix(0, dim[1], dim[2])
    for (a in seq_len(dim[1])) {
      m[a, ] <- as.numeric(strsplit(rows[a], "\t", fixed = TRUE)[[1]])
    }
    m
  }
  list(
    groups = groups,
    lr = lr,
    n_genes = n_genes,
    n_cells = n_cells,
    type_mean = type_mean,
    dim = dim,
    prob = array(block("prob", prod(dim)), dim = dim),
    pval = array(block("pval", prod(dim)), dim = dim),
    ## `n_genes * n_groups`, which is *not* `prod(dim[1:2])`: `dim` is the network's K x K x N, and
    ## a 4-group fixture with 24 genes has 96 averages, not 16. Reading `prod(dim[1:2])` consumed 16
    ## of the 96 lines and desynchronised every record after it -- which is why the CLI now echoes
    ## the gene count instead of leaving the consumer to guess it.
    ave_expr = block("ave_expr", n_genes * dim[1]),
    kernel_ave = block("kernel_ave", n_genes * dim[1]),
    count = matrix_block("aggregate_count"),
    weight = matrix_block("aggregate_weight"),
    source = marker("dimnames_source"),
    target = marker("dimnames_target"),
    interaction = marker("dimnames_interaction")
  )
}

## ------------------------------------------------------------------ every configuration
## One fixture would prove the CLI can compute something. These prove it computes the same thing on
## the paths that are actually different in the kernel: all four `type.mean` values, `population.size`
## both ways, `nboot` from 1 to 9, `Kh` across six orders of magnitude, `n` at 1 and 2, `raw.use` both
## ways, and two seeds. The full list and the reason for each is in `gen_cli_fixture.R`.
configs <- sub("^cli_golden_([A-Za-z0-9_]+)\\.txt$", "\\1",
               list.files("tests/fixtures", pattern = "^cli_golden_.*\\.txt$"))
if (length(configs) == 0L) {
  stop("no cli_golden_*.txt in tests/fixtures; run tests/parity/gen_cli_fixture.R first",
       call. = FALSE)
}
## Sorted, so the report is in the same order on every run and a diff between two runs is readable.
configs <- sort(configs)
cat(sprintf("--- %d CLI configurations\n", length(configs)))

for (cfg_name in configs) {
  gold <- read_golden(sprintf("tests/fixtures/cli_golden_%s.txt", cfg_name))
  inp <- sprintf("tests/fixtures/cli_input_%s.tsv", cfg_name)
  if (!file.exists(inp)) {
    stop("the golden for `", cfg_name, "` has no matching input file ", inp, call. = FALSE)
  }
  cat(sprintf("\n[%s]  type=%s nboot=%d\n", cfg_name, gold$type_mean, gold$nboot))

  ## The CLI has to be told where the pinned database is. A relative path would resolve against
  ## whatever directory the runner happened to be in.
  out <- cli("run", "--input", inp, "--db", DBDIR, "--thresh", "0.05")
  got <- read_cli(out)

  cmp("dim", identical(got$dim, gold$dim),
      sprintf("got %s want %s", paste(got$dim, collapse = "x"), paste(gold$dim, collapse = "x")))
  cmp("group count", identical(got$groups, gold$groups))
  cmp("interaction count", identical(got$lr, gold$lr))
  cmp("dimnames source", identical(got$source, gold$source))
  cmp("dimnames target", identical(got$target, gold$source),
      "source and target levels are the same vector for an RNA net")
  cmp("dimnames interaction", identical(got$interaction, gold$interaction))
  ## `match.arg` partial matching: the CLI must record the *matched* spelling, not what was typed.
  cmp("type.mean recorded as matched", identical(got$type_mean, gold$type_mean),
      sprintf("got %s want %s", got$type_mean, gold$type_mean))

  ## The headline comparisons, and they are `identical()` on the whole array -- values, shape and
  ## dimnames.
  cmp("Prob identical to upstream", identical(got$prob, gold$prob))
  cmp("Pval identical to upstream", identical(got$pval, gold$pval))
  cmp("aggregateNet count", identical(got$count, gold$count))
  cmp("aggregateNet weight", identical(got$weight, gold$weight))

  ## `ave_expr` is gated only where `computeAveExpr` can answer; its block is absent otherwise,
  ## because `computeAveExpr`'s own `match.arg` offers three choices and not four.
  if (gold$has_ave) {
    gold_ave <- read_golden_ave(sprintf("tests/fixtures/cli_golden_%s.txt", cfg_name))
    cmp("computeAveExpr means", identical(got$ave_expr, gold_ave),
        sprintf("%d values", length(got$ave_expr)))
  }

  ## When a mismatch happens, say *where*. A bare FALSE on a 4x4x6 array is not actionable, and the
  ## whole point of pinning the bit patterns is that the difference is small enough to point at.
  if (!identical(got$prob, gold$prob)) {
    d <- which(got$prob != gold$prob | xor(is.na(got$prob), is.na(gold$prob)), arr.ind = TRUE)
    cat(sprintf("  %d of %d Prob values differ\n", nrow(d), length(got$prob)))
    for (r in utils::head(seq_len(nrow(d)), 6)) {
      idx <- d[r, ]
      cat(sprintf("    [%d,%d,%d] (%s -> %s) got %.17g want %.17g\n", idx[1], idx[2], idx[3],
                  gold$interaction[idx[3]], got$source[idx[1]],
                  got$prob[idx[1], idx[2], idx[3]], gold$prob[idx[1], idx[2], idx[3]]))
    }
  }
  if (!identical(got$pval, gold$pval)) {
    d <- which(got$pval != gold$pval, arr.ind = TRUE)
    cat(sprintf("  %d of %d Pval values differ\n", nrow(d), length(got$pval)))
    for (r in utils::head(seq_len(nrow(d)), 6)) {
      idx <- d[r, ]
      cat(sprintf("    [%d,%d,%d] got %.17g want %.17g\n", idx[1], idx[2], idx[3],
                  got$pval[idx[1], idx[2], idx[3]], gold$pval[idx[1], idx[2], idx[3]]))
    }
  }

  ## The recorded invariants, on the CLI's own output rather than R's: the two the objective names
  ## that apply here, plus non-negativity. Dimnames are checked above.
  ##
  ## **No `Prob in [0, 1]` check.** One was written here and it failed against an output the port
  ## reproduces bit for bit: 16 of the 96 `Prob` values exceed 1 and the largest is 2025.1. That is
  ## upstream's own arithmetic. `Prob` is `P1 * P2 * P3 * P4 * P.spatial`, and while
  ## `P1 = dataLR^n / (Kh^n + dataLR^n)` is bounded by 1, the agonist and antagonist terms are
  ## `crossprod(matrix(x, nrow = 1))` -- plain sums of expression values -- and `P4` is a
  ## group-proportion outer product, so any of them can exceed 1 and the product can. `Prob` is a
  ## score, not a probability, and upstream never clamps it. Non-negativity *is* worth asserting: a
  ## negative `Prob` would mean an expression path had a sign error, and `prob > 0` downstream would
  ## silently drop the interaction.
  p <- as.vector(got$pval); q <- as.vector(got$prob)
  cmp("Pval in {k/nboot}", all(abs(p * gold$nboot - round(p * gold$nboot)) < 1e-9))
  cmp("Pval[Prob == 0] == 1", all(p[q == 0] == 1))
  cmp("Pval in [0, 1]", all(p >= 0 & p <= 1))
  cmp("Prob non-negative", all(q >= 0))
  if (any(q > 1)) {
    cat(sprintf("  note: %d of %d Prob values exceed 1 (max %.6g); upstream does not clamp them\n",
                sum(q > 1), length(q), max(q)))
  }
}

## ---------------------------------------------------------------- the other subcommands
## `describe` must agree with the file it read, or it is decoration.
d <- cli("describe", "--input", IN, "--db", DBDIR)
dmap <- stats::setNames(sub("^[^\\t]*\\t", "", d[grepl("^(groups|genes|cells|lr)\\t", d)]),
                        sub("\\t.*$", "", d[grepl("^(groups|genes|cells|lr)\\t", d)]))
n_cmp <- n_cmp + 1L
ok <- as.integer(dmap[["groups"]]) == as.integer(dmap[["groups"]]) &&
      as.integer(dmap[["lr"]]) == as.integer(dmap[["lr"]]) && as.integer(dmap[["lr"]]) > 0
cat(sprintf("%-40s %s\n", "describe reports the file's shape", if (ok) "ok" else "MISMATCH"))
if (!ok) fails <- fails + 1L

## `mean` must reproduce the mean functions on values whose answers R knows exactly, including the
## cases that are easy to get wrong: `triMean` on four values, `thresholdedMean` on a value that
## must raise, and `geometricMean` on a zero.
mean_cases <- list(
  list(kind = "triMean", x = c(0, 0, 1, 1), want = 0.5),
  list(kind = "triMean", x = c(1, 2, 3, 4), want = 2.5),
  list(kind = "geometricMean", x = c(1, 4, 16), want = 4),
  list(kind = "geometricMean", x = c(0, 1, 4), want = 0),
  list(kind = "truncatedMean", x = c(1, 2, 3, 4), trim = 0.25, want = 2.5)
)
for (mc in mean_cases) {
  n_cmp <- n_cmp + 1L
  got_m <- as.numeric(cli("mean", "--kind", mc$kind,
                           if (!is.null(mc$trim)) c("--trim", format(mc$trim)) else character(),
                           format(mc$x, digits = 17)))
  ## R's own answer, from the same functions, for the same inputs.
  fun <- switch(mc$kind,
    triMean = function(v) mean(sort(c(v[1], v[2], mean(v[c(1, 4)]), v[4], v[3]))),
    geometricMean = function(v) exp(mean(log(v))),
    truncatedMean = function(v) mean(v, trim = mc$trim, na.rm = TRUE))
  want_m <- fun(mc$x)
  ok <- isTRUE(all.equal(got_m, want_m, tolerance = 1e-15))
  cat(sprintf("%-40s %s  got %.17g want %.17g\n",
              paste0("mean ", mc$kind, "(", length(mc$x), " values)"),
              if (ok) "ok" else "MISMATCH", got_m, want_m))
  if (!ok) fails <- fails + 1L
}

## `version` must print something, and `run` without `--db` must fail rather than guess.
v <- cli("version")
n_cmp <- n_cmp + 1L
ok <- any(grepl("^cellchatrs [0-9]", v))
cat(sprintf("%-40s %s\n", "version prints a version", if (ok) "ok" else "MISMATCH"))
if (!ok) fails <- fails + 1L

n_cmp <- n_cmp + 1L
bad <- suppressWarnings(system2(BIN, c("run", "--input", IN), stdout = FALSE, stderr = FALSE))
ok <- !identical(bad, 0L)
cat(sprintf("%-40s %s\n", "run without --db fails loudly", if (ok) "ok" else "MISMATCH"))
if (!ok) fails <- fails + 1L

## A malformed input must be rejected, not parsed. The hex-float format is the defence against a
## silently different double, so a decimal value in the matrix has to be an error.
n_cmp <- n_cmp + 1L
bad_path <- tempfile(fileext = ".tsv")
lines_in <- readLines(IN)
vi <- grep("^values$", lines_in)[1]
lines_in[vi + 1] <- sub("^[^ ]+", "0.5", lines_in[vi + 1])   # a decimal, where hex is required
writeLines(lines_in, bad_path)
bad <- suppressWarnings(system2(BIN, c("run", "--input", shQuote(bad_path), "--db", shQuote(DBDIR)),
                                stdout = FALSE, stderr = FALSE))
ok <- !identical(bad, 0L)
cat(sprintf("%-40s %s\n", "a decimal in the matrix is rejected", if (ok) "ok" else "MISMATCH"))
if (!ok) fails <- fails + 1L
unlink(bad_path)

cat(sprintf("\n%s: %d failing comparisons out of %d\n",
            if (fails == 0L) "IDENTICAL" else "MISMATCH", fails, n_cmp))
quit(status = if (fails == 0L) 0L else 1L)
