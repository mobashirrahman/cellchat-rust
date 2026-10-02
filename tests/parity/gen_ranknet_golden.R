#!/usr/bin/env Rscript
## Golden corpus for `rankNet`'s numeric core: the per-pathway information flow, the
## `-1/log` rescaling with its degenerate-entry reassignment, `order()`, and the
## one-significant-digit relative contributions.
##
## `rankNet` is mostly ggplot, and the plot stays in R (PLAN.md 14.1). What this corpus pins is
## the part that decides *which value goes in which row*, which is the part `identical()` can
## check. The upstream expressions are therefore **lifted verbatim** into small local functions
## rather than reimplemented, so the oracle is upstream's arithmetic and not my reading of it.
##
## There is no RNG: every probability array is written out literally, which is what makes the
## degenerate cases (`pSum` exactly 1, exactly 0, negative, `NaN`) reachable on purpose rather
## than by luck.
##
## Output is tab-separated, NA escaped as <NA> and NaN as <NaN>.
ROOT <- Sys.getenv("CELLCHATRS_ROOT", unset = getwd())
if (!dir.exists(file.path(ROOT, "R"))) stop("run from the repository root, or set CELLCHATRS_ROOT")
suppressWarnings(suppressMessages(library(collapse)))

esc <- function(v) {
  if (is.character(v)) return(gsub("\t", "<TAB>", v, fixed = TRUE))
  if (is.logical(v)) return(ifelse(v, "TRUE", "FALSE"))
  if (is.na(v) && !is.nan(v)) return("<NA>")
  if (is.nan(v)) return("<NaN>")
  format(as.numeric(v), digits = 17, scientific = TRUE, trim = TRUE)
}
vec_rec <- function(tag, name, v) sprintf("%s\t%s\t%s", tag, name,
  paste(vapply(v, esc, ""), collapse = "\t"))

## ------------------------------------------------------------------ upstream, verbatim
## `rankNet`'s `mode = "single"` body, up to (but not including) the ggplot. Returned as a list
## so the corpus records every intermediate, not just the answer: `pSum` before and after the
## transform, which entries were flagged, the reassignment, and the post-sort frame.
up_single <- function(prob, pval, thresh, measure, sources.use = NULL, targets.use = NULL) {
  prob[as.vector(pval) > thresh] <- 0
  if (measure == "count") {
    prob <- 1 * (prob > 0)
  }
  if (!is.null(sources.use)) {
    if (all(sources.use %in% dimnames(prob)[[1]])) {
      sources.use <- match(sources.use, dimnames(prob)[[1]])
    } else {
      stop("The input `sources.use` should be cell group names or a numerical vector!")
    }
    idx.t <- setdiff(seq_len(nrow(prob)), sources.use)
    prob[idx.t, , ] <- 0
  }
  if (!is.null(targets.use)) {
    if (all(targets.use %in% dimnames(prob)[[1]])) {
      targets.use <- match(targets.use, dimnames(prob)[[1]])
    } else {
      stop("The input `targets.use` should be cell group names or a numerical vector!")
    }
    idx.t <- setdiff(seq_len(nrow(prob)), targets.use)
    prob[, idx.t, ] <- 0
  }
  if (sum(prob) == 0) stop("No inferred communications for the input!")

  pSum <- apply(prob, 3, sum)
  pSum.original <- pSum
  if (measure == "weight") {
    pSum <- -1 / log(pSum)
    pSum[is.na(pSum)] <- 0
    idx1 <- which(is.infinite(pSum) | pSum < 0)
    values.assign <- seq(max(pSum) * 1.1, max(pSum) * 1.5, length.out = length(idx1))
    position <- sort(pSum.original[idx1], index.return = TRUE)$ix
    pSum[idx1] <- values.assign[match(seq_along(idx1), position)]
  } else if (measure == "count") {
    pSum <- pSum.original
  }
  ## `idx1`, `values.assign` and `position` exist only inside the `measure == "weight"`
  ## branch. Returning them unconditionally makes the *count* cases fail with "object 'idx1'
  ## not found", which is an artefact of the harness rather than a property of `rankNet` --
  ## and it is exactly the kind of thing that gets recorded as an upstream bug by mistake.
  if (measure == "weight") {
    list(original = pSum.original, scaled = pSum, flagged = idx1,
         values.assign = values.assign, position = position)
  } else {
    list(original = pSum.original, scaled = pSum, flagged = integer(0),
         values.assign = numeric(0), position = integer(0))
  }
}

## `as.numeric(format(x, digits = 1))` then the `is.na` reset, verbatim.
up_relative <- function(num, den) {
  r <- as.numeric(format(num / den, digits = 1))
  r[is.na(r)] <- 0
  r
}

## `idx <- with(df, order(df$contribution))` -- R's `order` on a double vector, radix and so
## stable. The corpus records the permutation itself, which is the thing a port must match.
up_order <- function(x) as.integer(order(x))

## ------------------------------------------------------------------ upstream, verbatim (comparison)
## `rankNet`'s `mode = "comparison"` body, up to the ggplot. The expressions are copied from the
## pinned source rather than paraphrased, so the corpus records upstream's arithmetic.
##
## Returned per comparison plus the pooled intermediates, because in comparison mode the
## degenerate reassignment is **pooled across comparisons**: one `values.assign`, one
## `pSum.original.all`, one `position`. That is the detail the Rust port is most likely to get
## wrong, so the corpus records it.
up_comparison <- function(prob.list, pval.list, names.list, thresh, measure,
                          sources.use = NULL, targets.use = NULL) {
  ncomp <- length(prob.list)
  pSum <- list(); pSum.original <- list(); pair.name <- list(); idx <- list()
  pSum.original.all <- c(); object.names.comparison <- c()
  for (i in seq_len(ncomp)) {
    prob <- prob.list[[i]]; prob[as.vector(pval.list[[i]]) > thresh] <- 0
    if (measure == "count") prob <- 1 * (prob > 0)
    if (!is.null(sources.use)) {
      if (all(sources.use %in% dimnames(prob)[[1]])) {
        sources.use <- match(sources.use, dimnames(prob)[[1]])
      } else stop("The input `sources.use` should be cell group names or a numerical vector!")
      idx.t <- setdiff(seq_len(nrow(prob)), sources.use)
      prob[idx.t, , ] <- 0
    }
    if (!is.null(targets.use)) {
      if (all(targets.use %in% dimnames(prob)[[1]])) {
        targets.use <- match(targets.use, dimnames(prob)[[1]])
      } else stop("The input `targets.use` should be cell group names or a numerical vector!")
      idx.t <- setdiff(seq_len(nrow(prob)), targets.use)
      prob[, idx.t, ] <- 0
    }
    if (sum(prob) == 0) stop("No inferred communications for the input!")
    pSum.original[[i]] <- apply(prob, 3, sum)
    if (measure == "weight") {
      pSum[[i]] <- -1 / log(pSum.original[[i]])
      pSum[[i]][is.na(pSum[[i]])] <- 0
      idx[[i]] <- which(is.infinite(pSum[[i]]) | pSum[[i]] < 0)
      pSum.original.all <- c(pSum.original.all, pSum.original[[i]][idx[[i]]])
    } else if (measure == "count") {
      pSum[[i]] <- pSum.original[[i]]
    }
    pair.name[[i]] <- names(pSum.original[[i]])
    object.names.comparison <- c(object.names.comparison, sprintf("d%d", i))
  }
  if (measure == "weight") {
    values.assign <- seq(max(unlist(pSum)) * 1.1, max(unlist(pSum)) * 1.5,
                         length.out = length(unlist(idx)))
    position <- sort(pSum.original.all, index.return = TRUE)$ix
    for (i in seq_len(ncomp)) {
      if (i == 1) {
        pSum[[i]][idx[[i]]] <- values.assign[match(seq_along(idx[[i]]), position)]
      } else {
        pSum[[i]][idx[[i]]] <- values.assign[
          match(length(unlist(idx[seq_len(i - 1)])) + seq_along(unlist(idx[seq_len(i)])), position)]
      }
    }
  } else {
    values.assign <- numeric(0); position <- integer(0)
  }
  pair.name.all <- as.character(unique(unlist(pair.name)))
  ## `df[[i]] <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0,
  ## group = ..., row.names = pair.name.all)` then the **row-name-indexed** assignment. Pathways
  ## missing from a comparison keep their 0, which is what makes a union across datasets work.
  df <- list()
  for (i in seq_len(ncomp)) {
    df[[i]] <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0,
                          group = object.names.comparison[i], row.names = pair.name.all,
                          stringsAsFactors = FALSE)
    df[[i]][pair.name[[i]], 3] <- pSum[[i]]
    df[[i]][pair.name[[i]], 2] <- pSum.original[[i]]
  }
  ## Recorded twice on purpose: `ratio.raw` is the arithmetic the kernel returns, and
  ## `contribution.relative` is what `format(..., digits = 1)` makes of it. They are not
  ## interchangeable -- `format` is vector-dependent, so the rounded value of one element depends
  ## on its neighbours -- and only the first is a kernel quantity.
  ratio.raw <- list()
  contribution.relative <- list()
  for (i in seq_len(ncomp - 1)) {
    ratio.raw[[i]] <- df[[ncomp - i + 1]]$contribution / df[[1]]$contribution
    ratio.raw[[i]][is.na(ratio.raw[[i]])] <- 0
    contribution.relative[[i]] <- as.numeric(format(ratio.raw[[i]], digits = 1))
    contribution.relative[[i]][is.na(contribution.relative[[i]])] <- 0
  }
  names(contribution.relative) <- paste0("contribution.relative.", seq_along(contribution.relative))
  for (i in seq_len(ncomp)) for (j in seq_along(contribution.relative)) {
    df[[i]][[names(contribution.relative)[j]]] <- contribution.relative[[j]]
  }
  data2 <- df[[ncomp]]$contribution
  c1 <- -contribution.relative[[1]]
  keys <- list(c1)
  if (ncomp == 3) keys <- list(c1, -contribution.relative[[2]], df[[1]]$contribution, -data2)
  if (ncomp >= 4) {
    keys <- list(c1, -contribution.relative[[2]], -contribution.relative[[3]],
                 df[[1]]$contribution, -data2)
  }
  if (ncomp == 2) keys <- list(c1, df[[1]]$contribution, -data2)
  ord <- do.call(order, keys)
  list(original = pSum.original, scaled = pSum, flagged = idx,
       values.assign = values.assign, position = position,
       pair.name.all = pair.name.all, ord = as.integer(ord),
       relative = contribution.relative, ratio.raw = ratio.raw,
       contrib = lapply(df, `[[`, "contribution"),
       scaled_all = lapply(df, `[[`, "contribution.scaled"),
       names.all = lapply(df, function(d) as.character(d$name)[ord]),
       object.names = object.names.comparison)
}

## ------------------------------------------------------------------ the fixtures
## `k = 4` groups, `n` pathways. Every array is written out so the corpus is reproducible.
LEV <- c("g1", "g2", "g3", "g10")
mk <- function(prob_rows, names, k = 4) {
  ## The extent is the number of *pathways*, not the length of the first slice: `prob_rows[[1]]`
  ## is one k x k matrix flattened, and using it as the extent is a shape error that only shows
  ## up when k^2 happens to differ from n.
  stopifnot(length(names) == length(prob_rows),
            all(vapply(prob_rows, length, 0L) == k * k))
  arr <- array(0, dim = c(k, k, length(prob_rows)),
               dimnames = list(LEV, LEV, names))
  for (c in seq_along(names)) {
    m <- matrix(prob_rows[[c]], k, k)
    arr[, , c] <- m
  }
  arr
}
## A p-value array with the same shape; `pv` is a function of the index so the `> thresh` cut
## and the `Pval[Prob==0]==1` invariant are both satisfiable by construction.
mkp <- function(k, n, f) array(vapply(seq_len(k * k * n), f, 0), dim = c(k, k, n))

cases <- list()

## -- the degenerate `-1/log` cases, one per kind ------------------------------------------
## A pathway whose total is exactly 1 gives `log(1) == 0` and `-1/0 == -Inf`, which *is* flagged.
## One whose total is exactly 0 gives `log(0) == -Inf` and `-1/-Inf == +0`, which is **not**
## flagged and looks like an ordinary value. One with a negative total gives `NaN`, which the
## `is.na` reset turns into 0 *before* the flag test, so it is never reassigned either. All
## three in one array, because whether a flag survives depends on the others.
cases$weight_degenerate <- list(
  names = c("one", "zero", "neg", "two"),
  prob = mk(list(c(0.25, 0, 0, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(-0.5, 0, 0, 0, -0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0.5, 0.5, 0, 0, 0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), c("one", "zero", "neg", "two")),
  pval = mkp(4, 4, function(i) 0.01), thresh = 0.05, measure = "weight")

## Three flagged entries, so `values.assign` has length 3 and the rank inversion is visible.
## Totals: 1, 1.5 and 3 (only the first is degenerate), so `max` is finite and the sequence is
## well defined.
cases$weight_three_degenerate <- list(
  names = c("d1", "d2", "d3", "ok1", "ok2"),
  prob = mk(list(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0.5, 0.5, 0.5, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(1, 0.5, 0.5, 0, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0.1, 0.1, 0, 0, 0.1, 0.1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
          c("d1", "d2", "d3", "ok1", "ok2")),
  pval = mkp(4, 5, function(i) 0.01), thresh = 0.05, measure = "weight")

## Every pathway degenerate: `max(pSum)` is then finite only if at least one is not `-Inf`.
## Here all totals are 1, so `pSum` is all `-Inf`, `max` is `-Inf`, and
## `seq(-Inf, -Inf, length.out = n)` is `NaN` -- so the reassignment writes `NaN`. Reached only
## by having *no* non-degenerate pathway, which is why it needs its own case.
cases$weight_all_degenerate <- list(
  names = c("a", "b", "c"),
  prob = mk(list(c(0.5, 0.5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0.25, 0.25, 0.25, 0, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), c("a", "b", "c")),
  pval = mkp(4, 3, function(i) 0.01), thresh = 0.05, measure = "weight")

## -- `thresh` cutting, and the count measure ---------------------------------------------
cases$thresh_cuts <- list(
  names = c("p1", "p2", "p3"),
  prob = mk(list(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
               c(0.1, 0, 0, 0, 0.2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0, 0, 0, 0, 0, 0, 0, 0, 3, 0, 0, 0, 0, 0, 0, 0)), c("p1", "p2", "p3")),
  pval = mkp(4, 3, function(i) if (i <= 16) c(0.01, 0.04, 0.06, 0.2, 0.01, 0.01, 0.3, 0.01,
                                            0.06, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01, 0.01)[i] else 0.01),
  thresh = 0.05, measure = "weight")

cases$count_measure <- list(
  names = c("p1", "p2", "p3"),
  prob = mk(list(c(1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0.25, 0.25, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), c("p1", "p2", "p3")),
  pval = mkp(4, 3, function(i) 0.01), thresh = 0.05, measure = "count")

## -- group filters ------------------------------------------------------------------------
cases$sources_filter <- list(
  names = c("p1", "p2"),
  prob = mk(list(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
               c(1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), c("p1", "p2")),
  pval = mkp(4, 2, function(i) 0.01), thresh = 0.05, measure = "weight",
  sources.use = c("g1", "g10"))

## `targets.use` is validated against `dimnames(prob)[[1]]` -- the *source* axis -- so a name
## that exists only on the target axis is rejected. `g1` exists on both, so this one is
## accepted; the rejected variant is its own case below.
cases$targets_filter <- list(
  names = c("p1", "p2"),
  prob = mk(list(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
               c(1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), c("p1", "p2")),
  pval = mkp(4, 2, function(i) 0.01), thresh = 0.05, measure = "weight",
  targets.use = c("g1", "g2"))

cases$both_filters <- list(
  names = c("p1", "p2", "p3"),
  prob = mk(list(c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16),
               c(1, 1, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(0, 0, 0, 0, 0, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0)), c("p1", "p2", "p3")),
  pval = mkp(4, 3, function(i) 0.01), thresh = 0.05, measure = "weight",
  sources.use = "g1", targets.use = "g2")

## An all-zero network: `sum(prob) == 0` stops *before* the per-pathway sums, so this is an
## error and not a vector of zeros.
cases$all_zero <- list(
  names = c("p1", "p2"),
  prob = mk(list(rep(0, 16), rep(0, 16)), c("p1", "p2")),
  pval = mkp(4, 2, function(i) 0.01), thresh = 0.05, measure = "weight")

## The `thresh` cut can empty the network even when the raw array is not all zero.
cases$thresh_empties <- list(
  names = c("p1", "p2"),
  prob = mk(list(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)), c("p1", "p2")),
  pval = mkp(4, 2, function(i) 0.5), thresh = 0.05, measure = "weight")

## -- `order()` ----------------------------------------------------------------------------
## Ties are the point: a stable radix sort keeps the array order of equal contributions, and an
## unstable sort permutes them. The group names include `g10` so a byte-wise or level-order sort
## is distinguishable from the numeric one -- though here the key is a number, so that is a
## check on `order_f64` rather than on string ordering.
## `order()` is only a test of *stability* if the key has both ties and distinct values, so the
## pathways differ in how many cells are non-zero: 1, 3, 1, 2, 3. Under `measure = "count"` the
## totals are those counts, giving two tied pairs (a with c, b with e) around a distinct d.
cases$order_ties <- list(
  names = c("a", "b", "c", "d", "e"),
  prob = mk(list(c(1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(7, 7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0),
               c(2, 2, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0)),
          c("a", "b", "c", "d", "e")),
  pval = mkp(4, 5, function(i) 0.01), thresh = 0.05, measure = "count")

## -- the one-significant-digit relative contribution ---------------------------------------
## `as.numeric(format(num/den, digits = 1))`: ratios that agree to three digits collapse, and
## the rounding is *significant* digits, so 0.0432 -> 0.04 while 1234 -> 1e+03.
rel_cases <- list(
  list(name = "collapse", num = c(1, 1.004, 1.04, 0.996), den = c(1, 1, 1, 1)),
  list(name = "significant", num = c(0.0432, 0.00432, 1234, 98765), den = c(1, 1, 1, 1)),
  list(name = "signs", num = c(-2, 2, -0.5, 0.5), den = c(1, 1, 1, 1)),
  list(name = "zeros", num = c(0, 1e-12, 1e12, 0), den = c(1, 1, 1, 1e-12)),
  ## A zero denominator gives `Inf` or `NaN`, `as.numeric(format(NaN))` is `NA`, and the
  ## `is.na` reset turns it into 0. A port that keeps `NaN` orders the row differently.
  list(name = "zero_denominator", num = c(1, 2, 3), den = c(1, 0, 0)),
  ## `format` pads the vector to a common width, so the *characters* differ between neighbours
  ## while the parsed values do not. Recorded to prove only the rounding survives.
  list(name = "padded", num = c(0.5, 50, 5000), den = c(1, 1, 1))
)

## A 4x4 matrix with `a` at [1,1] and `b` at [2,2], flattened column-major for `mk`. Writing the
## comparison fixtures as nested `rep(list(rep(c(...))))` literals made them unreadable and
## impossible to bracket-balance by eye; this says what each fixture actually is.
diag2 <- function(a, b) {
  m <- matrix(0, 4, 4); m[1, 1] <- a; m[2, 2] <- b; as.vector(m)
}
## `arr3` takes the *flattened* matrices back to back -- so `rep(diag2(a, b), n)` is `n`
## identical pathways -- and splits them, because `mk` wants a list of `k * k` vectors and
## making the caller write `list(v, v, v, v)` by hand is what made the fixtures unreadable.
arr3 <- function(vecs, names) {
  stopifnot(length(vecs) == 16L * length(names))
  mk(lapply(seq_along(names), function(i) as.numeric(vecs[((i - 1L) * 16L + 1L):(i * 16L)])),
     names)
}

## -- the comparison-mode fixtures -------------------------------------------------------------
## `k = 4` groups, each comparison its own pathway vocabulary so that `pair.name.all` is a real
## union and the row-name-indexed assignment has to leave absent pathways at 0.
cmp_cases <- list(
  ## Two datasets with different pathway vocabularies (A B C D vs B X C). `X` is absent from
  ## comparison 1, so its relative contribution is `contribution[X]/0` = `Inf` -- and `Inf`
  ## **survives** the `is.na` reset, which is why `X` sorts first. A port that treats "absent"
  ## as 0 orders the whole frame differently.
  list(name = "two_disjoint",
       names = list(c("A", "B", "C", "D"), c("B", "X", "C")),
       prob = list(
         arr3(rep(diag2(0.25, 0.25), 4), c("A", "B", "C", "D")),
         arr3(c(diag2(0.5, 0.5), diag2(0.1, 0.1), diag2(0.2, 0.2)), c("B", "X", "C"))),
       thresh = 0.05, measure = "weight"),
  ## **Degenerate entries in both comparisons at once**, which is the only way to see the
  ## pooling. Comparison 1 has a pathway whose total is exactly 1 (`-1/log(1)` = `-Inf`,
  ## flagged) and comparison 2 has one too, so the pooled `values.assign` spans both and
  ## `position` ranks across comparisons. Reassigning per comparison would give different
  ## assigned values here, and therefore a different `order`.
  list(name = "pooled_degenerate",
       names = list(c("one", "two", "three", "four"), c("one", "X", "three", "four")),
       ## Totals: c1 = 1, 2, 0.5, 0.25 and c2 = 1, 3, 0.5, 0.25, so **two** pathways per
       ## comparison are flagged and the pooled `values.assign` has four entries.
       ##
       ## Worth being precise about which, because the obvious answer is wrong. `-1/log(x)` is
       ## negative for every `x > 1`, so `which(is.infinite(pSum) | pSum < 0)` flags *any* pathway
       ## whose total reaches 1 -- the `log(1) == 0` case is just the boundary of a whole range,
       ## not a special case. Only totals below 1 stay positive. That matters here because a
       ## total above 1 is entirely ordinary: `pSum.original` sums a `k x k` slice of
       ## probabilities, and four cells at 0.5 already reach 2.
       ##
       ## The remaining two pathways per comparison are non-degenerate, so pooled `max(pSum)` is
       ## finite and `seq` does not stop. Without one this case is indistinguishable from
       ## `all_degenerate` and tests nothing about pooling.
       prob = list(
         arr3(c(diag2(0.5, 0.5), diag2(1.0, 1.0), diag2(0.25, 0.25), diag2(0.125, 0.125)),
              c("one", "two", "three", "four")),
         arr3(c(diag2(0.5, 0.5), diag2(1.5, 1.5), diag2(0.25, 0.25), diag2(0.125, 0.125)),
              c("one", "X", "three", "four"))),
       thresh = 0.05, measure = "weight"),
  ## Three comparisons: `relative.1` is `df[[3]]/df[[1]]` and `relative.2` is `df[[2]]/df[[1]]`
  ## -- the **last** over the first, not consecutive pairs -- so the sort has three keys before
  ## `contribution`. Dividing consecutive pairs produces a different `order`.
  list(name = "three_comparisons",
       names = list(c("a", "b", "c", "d"), c("a", "b", "e", "f"), c("a", "g", "h", "i")),
       prob = list(
         arr3(c(diag2(0.1, 0.1), diag2(0.2, 0.2), diag2(0.3, 0.3), diag2(0.4, 0.4)),
              c("a", "b", "c", "d")),
         arr3(c(diag2(0.2, 0.2), diag2(0.4, 0.4), diag2(0.6, 0.6), diag2(0.8, 0.8)),
              c("a", "b", "e", "f")),
         arr3(c(diag2(0.3, 0.3), diag2(0.6, 0.6), diag2(0.9, 0.9), diag2(1.2, 1.2)),
              c("a", "g", "h", "i"))),
       thresh = 0.05, measure = "weight"),
  ## Four comparisons, the widest `order(...)` upstream spells out before its `else` reaches for
  ## a `contribution.relative.4` that does not exist.
  list(name = "four_comparisons",
       names = list(c("a", "b"), c("a", "c"), c("a", "d"), c("a", "e")),
       prob = list(
         arr3(c(diag2(0.1, 0.1), diag2(0.2, 0.2)), c("a", "b")),
         arr3(c(diag2(0.2, 0.2), diag2(0.4, 0.4)), c("a", "c")),
         arr3(c(diag2(0.3, 0.3), diag2(0.6, 0.6)), c("a", "d")),
         arr3(c(diag2(0.4, 0.4), diag2(0.8, 0.8)), c("a", "e"))),
       thresh = 0.05, measure = "weight"),
  ## `measure = "count"`: no transform, no flagged entries, `values.assign` of length 0, and
  ## `contribution.scaled == contribution`. The `else` upstream takes is never reached, so the
  ## harness must not record `idx1` for these.
  list(name = "count_measure",
       names = list(c("a", "b", "c"), c("a", "d", "e")),
       prob = list(
         arr3(c(diag2(0.2, 0.3), diag2(0.4, 0.4), diag2(0.6, 0.6)), c("a", "b", "c")),
         arr3(c(diag2(0.9, 0.1), diag2(0.8, 0.2), diag2(0.7, 0.3)), c("a", "d", "e"))),
       thresh = 0.05, measure = "count"),
  ## The `thresh` cut is applied per comparison, and the `sum(prob) == 0` stop is checked
  ## **inside** the loop -- so a comparison that empties stops the whole call with the earlier
  ## comparisons' work already done and discarded.
  list(name = "thresh_empties_second",
       names = list(c("a", "b"), c("a", "b")),
       prob = list(
         arr3(c(diag2(0.5, 0.5), diag2(0.5, 0.5)), c("a", "b")),
         arr3(c(diag2(0.5, 0.5), diag2(0.5, 0.5)), c("a", "b"))),
       pval = list(mkp(4, 2, function(i) 0.01), mkp(4, 2, function(i) 0.9)),
       thresh = 0.05, measure = "weight"),
  ## The axis filters, applied per comparison. `targets.use` is validated against
  ## `dimnames(prob)[[1]]` upstream -- the *first* axis -- so a name that exists only as a target
  ## is accepted and then matches nothing. `g1` exists on both axes, `g2` only as a target.
  list(name = "filters",
       names = list(c("a", "b", "c"), c("a", "b", "c")),
       prob = list(
         arr3(c(diag2(0.5, 0.5), diag2(0.4, 0.4), diag2(0.3, 0.3)), c("a", "b", "c")),
         arr3(c(diag2(0.2, 0.2), diag2(0.6, 0.6), diag2(0.1, 0.1)), c("a", "b", "c"))),
       thresh = 0.05, measure = "weight",
       sources.use = "g1", targets.use = "g2"),
  ## Ties on every sort key, so the radix `order`'s **stability** decides the row order -- and the
  ## stable order is the order of `pair.name.all`, i.e. first appearance across comparisons.
  ## Ties on every sort key, so the radix `order`'s **stability** decides the row order -- and the
  ## stable order is the order of `pair.name.all`, i.e. first appearance across comparisons.
  ## Totals are 2 rather than 1: a total of exactly 1 makes `pSum` `-Inf`, every pathway gets
  ## flagged, `max(pSum)` is `-Inf` and `seq` stops -- so "all tied" has to be built with a
  ## *finite* `pSum` or the tie is never reached.
  list(name = "order_ties",
       names = list(c("z", "y", "x", "w"), c("z", "y", "v", "u")),
       prob = list(
         arr3(rep(diag2(1.0, 1.0), 4), c("z", "y", "x", "w")),
         arr3(rep(diag2(1.0, 1.0), 4), c("z", "y", "v", "u"))),
       thresh = 0.05, measure = "weight"),
  ## Every pathway degenerate in every comparison, so pooled `max(pSum)` is `-Inf` and `seq`
  ## stops with "'from' must be a finite number". One non-degenerate entry *anywhere* would keep
  ## `seq` finite, which is why the all-degenerate case is its own fixture.
  list(name = "all_degenerate",
       names = list(c("one", "two"), c("one", "two")),
       prob = list(
         arr3(c(diag2(0.5, 0.5), diag2(0.5, 0.5)), c("one", "two")),
         arr3(c(diag2(0.5, 0.5), diag2(0.5, 0.5)), c("one", "two"))),
       thresh = 0.05, measure = "weight")
)

## ------------------------------------------------------------------ run and record
out <- file.path(ROOT, "tests/fixtures/ranknet_golden.txt")
q <- file(paste0(out, ".part"), "wt")
writeLines(sprintf("n_cases\t%d", length(cases)), q)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  k <- 4L; n <- length(cs$names)
  writeLines(sprintf("case\t%s", nm), q)
  writeLines(sprintf("meta\t%s\t%s\t%s\t%s\t%d\t%d", nm, cs$measure,
    format(as.numeric(cs$thresh), digits = 17, scientific = TRUE, trim = TRUE),
    paste(cs$names, collapse = ","), k, n), q)
  ## The group levels, because `sources.use` / `targets.use` are cell-group names and the
  ## pathway names are a different vocabulary entirely -- resolving one against the other
  ## silently yields an empty index set and an all-zero flow.
  writeLines(sprintf("levels\t%s\t%s", nm, paste(LEV, collapse = ",")), q)
  writeLines(sprintf("sources\t%s\t%s", nm,
    if (is.null(cs$sources.use)) "-" else paste(cs$sources.use, collapse = ",")), q)
  writeLines(sprintf("targets\t%s\t%s", nm,
    if (is.null(cs$targets.use)) "-" else paste(cs$targets.use, collapse = ",")), q)
  writeLines(vec_rec("prob", nm, as.vector(cs$prob)), q)
  writeLines(vec_rec("pval", nm, as.vector(cs$pval)), q)
  res <- tryCatch(
    suppressWarnings(suppressMessages(up_single(cs$prob, cs$pval, cs$thresh, cs$measure,
                                                cs$sources.use, cs$targets.use))),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "rn_error"))
  if (inherits(res, "rn_error")) {
    writeLines(sprintf("error\t%s\t%s", nm, gsub("\t", "<TAB>", res$msg, fixed = TRUE)), q)
  } else {
    writeLines(vec_rec("original", nm, res$original), q)
    writeLines(vec_rec("scaled", nm, res$scaled), q)
    writeLines(sprintf("flagged\t%s\t%s", nm, paste(res$flagged, collapse = ",")), q)
    writeLines(sprintf("values_assign\t%s\t%s", nm,
      if (length(res$values.assign)) paste(vapply(res$values.assign, esc, ""), collapse = ",") else "-"), q)
    writeLines(sprintf("position\t%s\t%s", nm,
      if (length(res$position)) paste(res$position, collapse = ",") else "-"), q)
    ## `idx <- with(df, order(df$contribution))` -- the sort is on the *original* totals.
    writeLines(sprintf("order\t%s\t%s", nm, paste(up_order(res$original), collapse = ",")), q)
  }
}
writeLines(sprintf("n_cmp_cases\t%d", length(cmp_cases)), q)
for (cc in cmp_cases) {
  nm <- cc$name
  ncomp <- length(cc$names)
  pval.list <- if (!is.null(cc$pval)) cc$pval else
    lapply(cc$names, function(nms) mkp(4, length(nms), function(i) 0.01))
  writeLines(sprintf("cmp_case\t%s", nm), q)
  writeLines(sprintf("cmp_meta\t%s\t%s\t%s\t%d\t%s", nm, cc$measure,
    format(as.numeric(cc$thresh), digits = 17, scientific = TRUE, trim = TRUE),
    ncomp, paste(vapply(cc$names, paste, character(1), collapse = ","), collapse = ";")), q)
  writeLines(sprintf("cmp_levels\t%s\t%s", nm, paste(LEV, collapse = ",")), q)
  writeLines(sprintf("cmp_sources\t%s\t%s", nm,
    if (is.null(cc$sources.use)) "-" else cc$sources.use), q)
  writeLines(sprintf("cmp_targets\t%s\t%s", nm,
    if (is.null(cc$targets.use)) "-" else cc$targets.use), q)
  for (i in seq_len(ncomp)) {
    writeLines(vec_rec("cmp_prob", sprintf("%s#%d", nm, i), as.vector(cc$prob[[i]])), q)
    writeLines(vec_rec("cmp_pval", sprintf("%s#%d", nm, i), as.vector(pval.list[[i]])), q)
  }
  res <- tryCatch(
    suppressWarnings(suppressMessages(up_comparison(cc$prob, pval.list, cc$names, cc$thresh,
                                                    cc$measure, cc$sources.use, cc$targets.use))),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "rn_error"))
  if (inherits(res, "rn_error")) {
    writeLines(sprintf("cmp_error\t%s\t%s", nm, gsub("\t", "<TAB>", res$msg, fixed = TRUE)), q)
    next
  }
  writeLines(sprintf("cmp_names_all\t%s\t%s", nm, paste(res$pair.name.all, collapse = ",")), q)
  writeLines(sprintf("cmp_values_assign\t%s\t%s", nm,
    if (length(res$values.assign)) paste(vapply(res$values.assign, esc, ""), collapse = ",")
    else "-"), q)
  writeLines(sprintf("cmp_position\t%s\t%s", nm,
    if (length(res$position)) paste(res$position, collapse = ",") else "-"), q)
  writeLines(sprintf("cmp_order\t%s\t%s", nm, paste(res$ord, collapse = ",")), q)
  for (i in seq_len(ncomp)) {
    tag <- sprintf("%s#%d", nm, i)
    writeLines(vec_rec("cmp_original", tag, res$original[[i]]), q)
    writeLines(vec_rec("cmp_scaled", tag, res$scaled[[i]]), q)
    ## `idx[[i]]` only exists inside the `measure == "weight"` branch -- for `count`, upstream
    ## leaves `idx` as the empty `list()` it was initialised to. Recorded as absent rather than as
    ## an empty vector, because "no flagged entries" and "this list was never populated" are
    ## different facts and only the first is a claim about the data.
    writeLines(sprintf("cmp_flagged\t%s\t%s", tag,
      if (length(res$flagged) < i) "-" else
      if (length(res$flagged[[i]])) paste(res$flagged[[i]], collapse = ",") else "-"), q)
    writeLines(vec_rec("cmp_contrib", tag, res$contrib[[i]]), q)
    writeLines(vec_rec("cmp_scaled_all", tag, res$scaled_all[[i]]), q)
    writeLines(sprintf("cmp_names_ordered\t%s\t%s", tag,
      paste(res$names.all[[i]], collapse = ",")), q)
  }
  for (j in seq_along(res$relative)) {
    writeLines(vec_rec("cmp_ratio", sprintf("%s#%d", nm, j), res$ratio.raw[[j]]), q)
    writeLines(vec_rec("cmp_relative", sprintf("%s#%d", nm, j), res$relative[[j]]), q)
  }
  writeLines(sprintf("cmp_object_names\t%s\t%s", nm,
    paste(res$object.names, collapse = ",")), q)
}
for (rc in rel_cases) {
  writeLines(sprintf("relcase\t%s", rc$name), q)
  writeLines(vec_rec("rel_num", rc$name, rc$num), q)
  writeLines(vec_rec("rel_den", rc$name, rc$den), q)
  writeLines(vec_rec("relative", rc$name, suppressWarnings(up_relative(rc$num, rc$den))), q)
}
close(q)
if (!file.rename(paste0(out, ".part"), out)) stop("could not install ", out)
cat("wrote", length(cases), "flow cases,", length(cmp_cases), "comparison cases and",
    length(rel_cases), "relative cases to", out, "\n")
