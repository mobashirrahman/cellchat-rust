#!/usr/bin/env Rscript
## The configuration matrix for the kernel, and a machine-checkable coverage report.
##
## `check_identical.R` covers the *functions*: one or two configurations each, chosen to reach a
## branch. The definition of done asks for a **200-400 configuration matrix with pairwise coverage**
## across thirteen axes, which is a different thing: not "does each branch work" but "does every
## *pair* of settings ever co-occur, and is that true for enough settings that the matrix means
## something".
##
## This file builds the matrix. It does not run anything -- `check_matrix.R` does that. Keeping the
## two apart means the matrix can be regenerated, inspected and coverage-checked without spending
## the runtime, and the coverage claim can be checked independently of whether the ports pass.
##
## `build_matrix()` is **deterministic**: same axes in, same matrix out, no RNG. A matrix that
## changed shape between runs could not be cited as a coverage claim.

## ------------------------------------------------------------------ the axes
##
## Thirteen, matching the definition of done: type, nC, K, nLR, nboot, population.size, raw.use,
## datatype, LR structure, Kh, n, trim, and pathological data.
##
## Levels are ordered so that `config[i, axis]` cycling gives a sensible single-factor sweep, and
## so the *first* level of each axis is the "default" one -- which is what the function defaults
## say. `levels[[axis]][1]` is therefore the value a caller who changes nothing would get.
matrix_axes <- list(
  type     = c("triMean", "truncatedMean", "thresholdedMean", "median"),
  ## nC: cell counts. 12 is the smallest that gives every group >= 2 cells with K = 4; 60 is the
  ## size of the pre-existing hand fixture, so the matrix overlaps what is already gated.
  nC       = c(60L, 12L, 24L),
  K        = c(4L, 2L, 3L),
  ## nLR: how many of the L-R table's rows to keep. 1 is the degenerate "one pathway" case.
  nLR      = c(8L, 1L, 2L, 4L),
  ## nboot: 1 exercises the "every bootstrap is the observed data" path, and with `Prob`
  ## invariance to nboot it must give pval == 1 everywhere.
  nboot    = c(5L, 1L, 3L),
  pop      = c(FALSE, TRUE),
  raw      = c(TRUE, FALSE),
  datatype = c("RNA", "spatial"),
  ## LR structure. `missing_subunit` is the R-ism 11 crash: a complex whose subunit is absent from
  ## the matrix, which upstream answers with "subscript out of bounds" and the port must match.
  lrstruct = c("plain", "complex", "cofactor", "agonist", "antagonist", "mixed", "missing_subunit"),
  Kh       = c(0.5, 1e3),
  n        = c(1, 2),
  trim     = c(0.1, 0.2, 0.3),
  ## Pathological data. `single_cell` has one cell, which no group can have two of.
  data     = c("normal", "all_zero", "single_cell", "with_na", "with_nan", "with_inf",
               "dup_rownames"),
  ## A constant factor on the whole expression matrix. `computeCommunProb` starts with
  ## `data.use <- data/max(data)`, so **upstream's `Prob` is invariant to this** -- it is a
  ## metamorphic axis, not a parameter axis.
  ##
  ## It is here because its absence hid a real bug for the whole life of the port. Every earlier
  ## fixture drew from `runif(0.01, 1)`, which reaches `1` often enough that `max(data) == 1`
  ## every time, the division was the identity, and a kernel that had *omitted the division
  ## entirely* still matched upstream on every gated configuration. The scale axis makes the
  ## maximum vary, and `check_matrix.R` additionally asserts the invariance across configurations
  ## that differ only here -- so a future removal of the step fails loudly instead of silently
  ## depending on an RNG draw.
  scale    = c(1, 0.5, 2, 3)
)

## Short names, used for the fixture identifiers and the coverage report.
matrix_axis_names <- names(matrix_axes)

## ------------------------------------------------------------------ pairwise coverage
##
## A *pair* is an axis, a level of it, another axis, and a level of that -- so `(nC, 12, K, 3)` is
## one pair. Pairwise coverage means every such combination appears at least once. Counting the pairs
## gives the theoretical minimum for the matrix size, which is what makes the 200-400 band
## meaningful rather than arbitrary.
##
## The key spells out **both axis names**, not just the two level strings. Keying on the levels
## alone is wrong in a way that flatters the result: `paste0(trim, "|", n)` and
## `paste0(Kh, "|", type)` can both be `"0.5|1"`, and `unique()` then collapses two genuinely
## distinct required pairs into one, which both understates the requirement (875 pairs collapse to
## 616) and lets one configuration be credited with covering both. The axis names remove the
## collision; `pair_key` is the single place a key is built, so the two readers cannot disagree.
pair_key <- function(axis_a, level_a, axis_b, level_b) {
  paste0(axis_a, "=", level_a, "|", axis_b, "=", level_b)
}

## Every pair the matrix is required to contain, as a character vector.
matrix_all_pairs <- function() {
  ax <- matrix_axes
  keys <- character(0)
  for (i in seq_along(ax)) {
    for (j in seq_along(ax)) {
      if (j <= i) next
      for (u in ax[[i]]) {
        for (v in ax[[j]]) {
          keys <- c(keys, pair_key(matrix_axis_names[i], u, matrix_axis_names[j], v))
        }
      }
    }
  }
  unique(keys)
}

## Which pairs a single configuration covers.
matrix_pairs_of <- function(cfg) {
  ax <- matrix_axes
  keys <- character(0)
  for (i in seq_along(ax)) {
    for (j in seq_along(ax)) {
      if (j <= i) next
      ## `ax[[i]]` is the axis's *levels*; the name is what indexes a configuration.
      keys <- c(keys, pair_key(names(ax)[i], cfg[[names(ax)[i]]],
                               names(ax)[j], cfg[[names(ax)[j]]]))
    }
  }
  keys
}

## ------------------------------------------------------------------ the candidate pool
##
## Enough candidates that the greedy selection below has room to choose. Built by sweeping each
## axis while holding the rest at their defaults, then sweeping *pairs* of axes, then triples --
## which is what lets the selection reach coverage instead of stalling on pairs that no
## single-axis sweep reaches.
matrix_candidates <- function() {
  ax <- matrix_axes
  base <- lapply(ax, `[`, 1)
  names(base) <- matrix_axis_names
  out <- list()
  add <- function(cfg, tag) {
    cfg <- as.list(cfg)
    cfg[["tag"]] <- tag
    out[[length(out) + 1L]] <<- cfg
  }
  add(base, "base")
  ## Single-axis sweeps.
  for (a in matrix_axis_names) {
    for (u in ax[[a]]) {
      cfg <- base; cfg[[a]] <- u
      add(cfg, paste0("one:", a, "=", u))
    }
  }
  ## Pair sweeps, both orders.
  for (a in matrix_axis_names) {
    for (b in matrix_axis_names) {
      if (b == a) next
      for (u in ax[[a]]) {
        for (v in ax[[b]]) {
          cfg <- base; cfg[[a]] <- u; cfg[[b]] <- v
          add(cfg, paste0("two:", a, "=", u, ";", b, "=", v))
        }
      }
    }
  }
  ## Triple sweeps over the axes with the most levels, where a pair sweep leaves gaps.
  wide <- names(ax)[order(-vapply(ax, length, 0L))][seq_len(4L)]
  for (a in wide) for (b in wide) for (c in wide) {
    if (a == b || b == c || a == c) next
    for (u in ax[[a]]) for (v in ax[[b]]) for (w in ax[[c]]) {
      cfg <- base; cfg[[a]] <- u; cfg[[b]] <- v; cfg[[c]] <- w
      add(cfg, paste0("three:", a, b, c))
    }
  }
  ## A deterministic Latin-style sweep over every axis at once, so the pool is not made only of
  ## near-default configurations and the greedy selection has genuinely diverse material.
  N <- 400L
  for (t in seq_len(N)) {
    cfg <- base
    for (k in seq_along(ax)) {
      cfg[[matrix_axis_names[k]]] <- ax[[k]][[((t - 1L) + (k - 1L) * 7L) %% length(ax[[k]]) + 1L]]
    }
    add(cfg, paste0("sweep:", t))
  }
  out
}

## ------------------------------------------------------------------ greedy selection
##
## Pick configurations that each add the most uncovered pairs, until the target size is reached or
## nothing adds anything. Deterministic: ties are broken by candidate order, and the candidates are
## generated in a fixed order.
matrix_select <- function(candidates, target) {
  want <- matrix_all_pairs()
  covered <- character(0)
  chosen <- integer(0)
  remaining <- seq_along(candidates)
  repeat {
    if (length(chosen) >= target || !length(remaining)) break
    ## Gain is the number of *still-uncovered required* pairs this candidate would add. Written as
    ## an intersection with `setdiff(want, covered)`: the obvious `setdiff(pairs, c(covered, want))`
    ## subtracts the whole required set and therefore scores every candidate at zero, which stops
    ## the loop immediately and leaves a matrix of nothing but padding.
    need_now <- setdiff(want, covered)
    gains <- vapply(remaining, function(i) {
      length(intersect(matrix_pairs_of(candidates[[i]]), need_now))
    }, 0L)
    ## Ties by candidate index: `which.max` returns the first maximum, and `remaining` is in
    ## generation order, so the choice is reproducible.
    best <- remaining[which.max(gains)]
    if (gains[match(best, remaining)] == 0L) break
    chosen <- c(chosen, best)
    covered <- union(covered, matrix_pairs_of(candidates[[best]]))
    remaining <- setdiff(remaining, best)
  }
  list(index = chosen, covered = covered, required = want)
}

## ------------------------------------------------------------------ the matrix
##
## `target` is the lower end of the 200-400 band the definition of done asks for. The greedy
## selection stops as soon as every pair is covered, which may be well below it -- a matrix of 300
## configurations that covers 300 distinct pairs is not better than one of 90 that covers all of
## them. So after coverage is reached the matrix is *padded* deterministically, with candidates
## that are not already chosen, until `target` rows exist. The padding adds breadth (triples, and
## repeated axis-level pairings in different contexts) rather than duplicating rows, and the report
## records how many rows each phase contributed so the distinction stays visible.
matrix_build <- function(target = 200L, max_size = 400L) {
  cands <- matrix_candidates()
  sel <- matrix_select(cands, target)
  idx <- sel$index
  n_min <- length(idx)
  ## Pad with unused candidates, preferring ones that share the fewest pairs with what is already
  ## chosen -- the opposite of the coverage phase, so the padding is genuinely additional.
  if (n_min < target) {
    remaining <- setdiff(seq_along(cands), idx)
    overlap <- vapply(remaining, function(i) {
      length(intersect(matrix_pairs_of(cands[[i]]), sel$covered))
    }, 0L)
    remaining <- remaining[order(overlap, remaining)]
    need <- min(target, max_size) - n_min
    idx <- c(idx, remaining[seq_len(need)])
  }
  cfg <- lapply(idx, function(i) cands[[i]])
  ## Coverage of the final matrix, recomputed rather than inherited, because the padding changes it.
  covered <- unique(unlist(lapply(cfg, matrix_pairs_of)))
  list(
    configs = cfg,
    n_minimal = n_min,
    covered = covered,
    required = sel$required,
    complete = setequal(covered, sel$required)
  )
}

## ------------------------------------------------------------------ reporting
##
## Coverage per axis pair, as a data frame. This is the artifact the "pairwise coverage" claim
## rests on, so it is emitted rather than summarised: a claim of "covers everything" should be
## checkable by reading which pairs are missing, and there should be none.
matrix_coverage_table <- function(m) {
  ax <- matrix_axes
  rows <- list()
  for (i in seq_along(ax)) {
    for (j in seq_along(ax)) {
      if (j <= i) next
      need <- as.vector(outer(ax[[i]], ax[[j]], function(u, v)
        pair_key(matrix_axis_names[i], u, matrix_axis_names[j], v)))
      have <- unique(unlist(lapply(m$configs, function(cfg) {
        c(pair_key(matrix_axis_names[i], cfg[[matrix_axis_names[i]]],
                   matrix_axis_names[j], cfg[[matrix_axis_names[j]]]))
      })))
      ## `outer` builds the first axis fastest-varying; normalise both sides to a set.
      rows[[length(rows) + 1L]] <- data.frame(
        axis_a = matrix_axis_names[i], axis_b = matrix_axis_names[j],
        levels_a = length(ax[[i]]), levels_b = length(ax[[j]]),
        pairs_required = length(need), pairs_covered = sum(need %in% have),
        complete = all(need %in% have), stringsAsFactors = FALSE)
    }
  }
  do.call(rbind, rows)
}

if (identical(environment(), globalenv()) && !interactive()) {
  m <- matrix_build()
  cat(sprintf("matrix: %d configurations (%d minimal for coverage, %d padded)\n",
              length(m$configs), m$n_minimal, length(m$configs) - m$n_minimal))
  cat(sprintf("pairs: %d covered of %d required; complete = %s\n",
              length(m$covered), length(m$required), m$complete))
  cov <- matrix_coverage_table(m)
  cat(sprintf("axis pairs: %d, all complete = %s\n", nrow(cov), all(cov$complete)))
  if (!all(cov$complete)) print(cov[!cov$complete, ])
}
