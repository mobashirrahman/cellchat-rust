#!/usr/bin/env Rscript
## Golden corpus for `computeRegionDistance`'s foundations: R's `mean(x, trim, na.rm)` and
## `collapse::fdist`.
##
## The k-d tree is deliberately **not** in this corpus. Annoy is approximate, so upstream is not
## a usable oracle for a neighbour query -- two runs of upstream can disagree. The exact tree is
## verified against exhaustive search instead (see `knn_matches_brute_force` in
## `src/rust/crates/r-core/tests/spatial_parity.rs`), and the *divergence* from Annoy is a separate
## measurement on real spatial data. What upstream *is* a usable oracle for is the arithmetic
## that consumes a neighbour's index and distance, and that is what this corpus pins.
##
## No RNG: every vector is written out, so the corpus is reproducible from the generator alone.
suppressWarnings(suppressMessages(library(collapse)))

## R 4.4 added `%||%`; this generator also runs on older R, where the idiom is spelled out.
`%||%` <- function(a, b) if (is.null(a)) b else a

esc <- function(v) {
  if (is.character(v)) return(gsub("\t", "<TAB>", v, fixed = TRUE))
  if (is.na(v) && !is.nan(v)) return("<NA>")
  if (is.nan(v)) return("<NaN>")
  if (is.infinite(v)) return(if (v > 0) "Inf" else "-Inf")
  format(as.numeric(v), digits = 17, scientific = TRUE, trim = TRUE)
}
## `vapply` over a zero-length vector needs a FUN.VALUE it can infer, and an explicit `""` is
## ambiguous -- R then tries to use `v` as the template and errors with "If v is left empty, x
## needs to be a matrix with at least 2 rows", which says nothing about the actual problem. The
## corpus has a deliberately empty input (`mean(numeric(0))`), so the empty case is real.
vec <- function(tag, name, v) {
  body_text <- if (length(v) == 0L) "-" else paste(vapply(v, esc, ""), collapse = "\t")
  sprintf("%s\t%s\t%s", tag, name, body_text)
}

## ------------------------------------------------------------------ the trimmed-mean cases
## `mean.default`'s three behaviours, each on an input that isolates it.
##
## 1. `na.rm = TRUE` drops `NaN` as well as `NA`, because `is.na(NaN)` is `TRUE`. A port that
##    filters only `NA` returns `NaN` here.
## 2. The trim count is `floor(n * trim)` from *each* end, so the kept count is
##    `n - 2*floor(n*trim)`. `n = 9` trims nothing, `n = 10` trims one, `n = 19` one and
##    `n = 20` two: not monotone in `n`, and different from `round(n * (1 - trim))`.
## 3. `trim >= 0.5` returns `median()`, and on an even `n` that is the mean of the two middle
##    values -- with `trim = 0.5` on an even vector, `(x[n/2] + x[n/2+1]) / 2`.
tm <- list(
  ## n below 10: `floor(n * 0.1) == 0`, so no trim and `mean` is the plain mean.
  list(name = "n1",  x = 1,                          trim = 0.1, na_rm = TRUE),
  list(name = "n2",  x = c(1, 2),                    trim = 0.1, na_rm = TRUE),
  list(name = "n9",  x = as.numeric(1:9),            trim = 0.1, na_rm = TRUE),
  ## n = 10: `floor(1) == 1` from each end, 8 kept.
  list(name = "n10", x = as.numeric(1:10),           trim = 0.1, na_rm = TRUE),
  list(name = "n11", x = as.numeric(1:11),           trim = 0.1, na_rm = TRUE),
  list(name = "n19", x = as.numeric(1:19),           trim = 0.1, na_rm = TRUE),
  ## n = 20: `floor(2) == 2` from each end, 16 kept.
  list(name = "n20", x = as.numeric(1:20),           trim = 0.1, na_rm = TRUE),
  list(name = "n100", x = as.numeric(seq_len(100)),  trim = 0.1, na_rm = TRUE),
  list(name = "n1000", x = as.numeric(seq_len(1000)), trim = 0.1, na_rm = TRUE),
  ## Unsorted input, and ties at the trim boundary: the *multiset* kept is the same whichever
  ## equal element the partial sort picks, so the mean must be too.
  list(name = "unsorted", x = c(9, 1, 8, 2, 7, 3, 6, 4, 5, 10), trim = 0.1, na_rm = TRUE),
  list(name = "ties", x = c(2, 2, 2, 2, 2, 2, 2, 2, 2, 9),  trim = 0.1, na_rm = TRUE),
  ## `NaN` with `na.rm = TRUE`: dropped, so the mean is over the rest.
  list(name = "nan_removed", x = c(1, NaN, 3, 4, 5, 6, 7, 8, 9, 10), trim = 0.1, na_rm = TRUE),
  ## `NA` with `na.rm = TRUE`: dropped. The count is taken *after* the removal, so the trim
  ## applies to the shortened vector -- 9 values, so `floor(0.9) == 0` and nothing is trimmed.
  list(name = "na_removed", x = c(1, NA, 3, 4, 5, 6, 7, 8, 9, 10), trim = 0.1, na_rm = TRUE),
  list(name = "na_nan_removed",
       x = c(1, NA, NaN, 4, 5, 6, 7, 8, 9, 10, 11, 12), trim = 0.1, na_rm = TRUE),
  ## `na.rm = FALSE` with an `NA` anywhere: `anyNA(x)` fires and the result is `NA`, *not* a
  ## mean of the non-missing values.
  list(name = "na_kept", x = c(1, NA, 3, 4, 5), trim = 0.1, na_rm = FALSE),
  list(name = "nan_kept", x = c(1, NaN, 3, 4, 5), trim = 0.1, na_rm = FALSE),
  ## All missing: `mean(numeric(0))` is `NaN`.
  list(name = "all_na", x = c(NA_real_, NA_real_), trim = 0.1, na_rm = TRUE),
  list(name = "empty", x = numeric(0),               trim = 0.1, na_rm = TRUE),
  ## `trim = 0` skips the whole trimming block, so the `anyNA` early return does not fire and
  ## `na.rm = FALSE` gives a plain mean -- `NaN` included, so `NaN`.
  list(name = "trim0_nan", x = c(1, NaN, 3),        trim = 0,   na_rm = FALSE),
  list(name = "trim0_na", x = c(1, NA, 3),          trim = 0,   na_rm = FALSE),
  ## The `trim >= 0.5` branch: `median`, including the even-length two-middle-value average.
  list(name = "trim_half_even", x = as.numeric(1:10), trim = 0.5, na_rm = TRUE),
  list(name = "trim_half_odd",  x = as.numeric(1:9),  trim = 0.5, na_rm = TRUE),
  list(name = "trim_over_half", x = as.numeric(1:10), trim = 0.75, na_rm = TRUE),
  ## `Inf` is not missing: it goes into the sum and the result is `Inf` (or `NaN` when the
  ## trimmed window is `+Inf - Inf`, which is what a symmetric `Inf` pair gives).
  list(name = "inf_kept", x = c(1, 2, Inf, 4, 5, 6, 7, 8, 9, 10), trim = 0.1, na_rm = TRUE),
  ## Negative values and a fractional vector, so the LONG_DOUBLE accumulation is exercised on
  ## something whose `f64` sum is not exact.
  list(name = "negative", x = c(-1e8, 1, -1, 1e8, 3, -3, 1e-8, -1e-8, 7, -7),
       trim = 0.1, na_rm = TRUE),
  list(name = "fractional", x = c(0.1, 0.2, 0.3, 0.4, 1/3, 1/7, 1e-300, 1e300, 0.5, 0.25),
       trim = 0.1, na_rm = TRUE)
)

## ------------------------------------------------------------------ the fdist cases
## `collapse::fdist` on a matrix. Written out rather than generated, and deliberately
## degenerate: a regular lattice (visium-like), duplicate points (the diagonal-adjacent
## duplicates a real slide produces), a single point, and two points.
## `fdist` rejects a single-row matrix outright -- `collapse` needs at least two rows to form a
## pair. Recorded as an error case rather than dropped, because "one coordinate" is exactly the
## degenerate input a caller reaches for and the message is the only thing that distinguishes
## "too few points" from "wrong shape".
fd <- list(
  list(name = "one", x = matrix(c(0, 0), nrow = 1, byrow = TRUE), expect_error = TRUE),
  list(name = "two", x = matrix(c(0, 0, 3, 4), nrow = 2, byrow = TRUE)),
  ## 3-4-5 triangle: the diagonal is `0`, and the off-diagonals are exactly 5 and 5.
  list(name = "right_triangle", x = matrix(c(0, 0, 3, 0, 0, 4), nrow = 3, byrow = TRUE)),
  ## A regular 3x3 lattice with unit spacing, so the distances are exact integers and the
  ## diagonal neighbours are all at distance 1 -- a tie-heavy case.
  list(name = "lattice3", x = as.matrix(expand.grid(x = 0:2, y = 0:2))),
  ## Two identical points: distance exactly `0`, and `fdist` does not turn that into `NA`.
  list(name = "duplicate", x = matrix(c(1, 1, 1, 1, 2, 2, 2, 2), nrow = 4, byrow = TRUE)),
  ## A collinear run, where the k-d tree's widest-axis rule has only one axis to choose from.
  list(name = "collinear", x = matrix(as.numeric(0:9), ncol = 1)),
  ## Irrational coordinates, so the sums are not exact and the last ulp matters.
  list(name = "irrational",
       x = matrix(c(sqrt(2), sqrt(3), sqrt(5), sqrt(7), pi, exp(1), 1/sqrt(2), 1/sqrt(3)),
                  ncol = 2, byrow = TRUE))
)

## ------------------------------------------------------------------ run and record
out <- file.path(Sys.getenv("CELLCHATRS_ROOT", unset = getwd()), "tests/fixtures/spatial_golden.txt")
q <- file(paste0(out, ".part"), "wt")
writeLines(sprintf("n_trim\t%d", length(tm)), q)
for (c in tm) {
  v <- suppressWarnings(mean(c$x, trim = c$trim, na.rm = c$na_rm))
  writeLines(sprintf("tm\t%s\t%s\t%s", c$name,
    format(as.numeric(c$trim), digits = 17, scientific = TRUE, trim = TRUE),
    if (c$na_rm) "TRUE" else "FALSE"), q)
  writeLines(vec("tm_x", c$name, as.numeric(c$x)), q)
  writeLines(sprintf("tm_mean\t%s\t%s", c$name, esc(v)), q)
}
writeLines(sprintf("n_fd\t%d", length(fd)), q)
for (c in fd) {
  writeLines(sprintf("fd\t%s\t%d\t%d", c$name, nrow(c$x), ncol(c$x)), q)
  if (isTRUE(c$expect_error)) {
    e <- tryCatch({ fdist(c$x); NULL }, error = function(e) conditionMessage(e))
    writeLines(sprintf("fd_error\t%s\t%s", c$name, e %||% "<no error>"), q)
    next
  }
  d <- suppressWarnings(fdist(c$x))
  ## `collapse::fdist` returns a **`dist`**, not a matrix: the lower triangle *without* the
  ## diagonal, in column-major order, with `Size` / `Diag` / `Upper` / `method` attributes. So
  ## `computeCellDistance`'s `d.spatial` is a `dist`, and the two statements that follow it --
  ## `d.spatial * ratio` and `d.spatial[d.spatial > x] <- NaN` -- go through `dist`'s own
  ## methods. The second one in particular writes only into the half a `dist` exposes, so the
  ## NaNs it introduces are invisible to anything that reads the result as a matrix. The
  ## corpus records `as.vector(d)`, which is the observable.
  stopifnot(inherits(d, "dist"), isFALSE(attr(d, "Diag")), isFALSE(attr(d, "Upper")),
            identical(attr(d, "method"), "euclidean"), attr(d, "Size") == nrow(c$x))
  ## Two records, because a `dist` has no `[` method and the traversal `as.vector()` uses is
  ## neither documented nor base R's.
  ##
  ## * `fd_dist` is upstream's exact observable: `as.vector(d)`, the compact triangle in
  ##   `collapse`'s own order. The test compares it as a **multiset** against the port's
  ##   triangle, which is the strongest statement available without pinning a traversal that
  ##   belongs to `as.vector` rather than to the arithmetic.
  ## * `fd_row` is the full symmetric matrix, computed here in plain R from the coordinates.
  ##   The test compares that bit for bit, which is what pins the *arithmetic*: the accumulation
  ##   order and the `sqrt`, where an `f64` sum in the wrong order is 1 ulp off on some inputs.
  ##
  ## The two together are as strong as pinning the traversal, and neither depends on a detail of
  ## `as.vector` that could change between `collapse` releases. The `dist` wrapper itself is an
  ## R-side concern, and `computeCellDistance`'s shim re-wraps the matrix.
  writeLines(vec("fd_dist", c$name, as.vector(d)), q)
  full <- matrix(0, nrow(c$x), nrow(c$x))
  for (i in seq_len(nrow(c$x))) for (j in seq_len(nrow(c$x))) {
    full[i, j] <- sqrt(sum((c$x[i, ] - c$x[j, ])^2))
  }
  writeLines(sprintf("fd_n\t%s\t%d", c$name, nrow(c$x)), q)
  ## The coordinates as well as the resulting matrix. They are *not* interchangeable: a reader
  ## that reconstructs coordinates from the distance rows is feeding the distances back in as if
  ## they were positions, and a two-point fixture then asks for the distance between (0, 5) and
  ## (5, 0) instead of between (0, 0) and (3, 4) -- which is sqrt(50), not 5. That is a
  ## plausible-looking number, so the two records are kept separate on purpose.
  ## `as.vector(t(x))`, **not** `as.numeric(x)`. `collapse::fdist` treats a matrix's *rows* as
  ## the points, so the coordinates have to be flattened row-major; `as.numeric()` on a matrix is
  ## column-major and transposes the point set. For a 2x2 that is a plausible-looking wrong
  ## answer -- the 3-4-5 triangle becomes the distance between (0, 3) and (0, 4), which is 1.
  writeLines(vec("fd_x", c$name, as.vector(t(c$x))), q)
  for (i in seq_len(nrow(c$x))) writeLines(vec("fd_row", sprintf("%s_%d", c$name, i), full[i, ]), q)
}
close(q)
if (!file.rename(paste0(out, ".part"), out)) stop("could not install ", out)
cat("wrote", length(tm), "trimmed-mean cases and", length(fd), "fdist cases to", out, "\n")
