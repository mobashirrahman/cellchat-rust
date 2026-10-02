#!/usr/bin/env Rscript
## Golden corpus for `computeRegionDistance`, with the exact k-d tree in place of `AnnoyParam`.
##
## ## Why this corpus is valid at all
##
## `Annoy` is a randomised *approximate* index, so upstream is not a usable oracle for a
## neighbour query: two runs can disagree, and the answer depends on a seed upstream does not
## expose. That is the whole reason the spatial branch is exempt from bit-identity.
##
## So the fixtures here are built so that Annoy **cannot** be wrong: the groups are placed far
## enough apart, with large margins, that every 1-NN is the true nearest neighbour by a wide
## margin and no two candidates are within any plausible error. On such a fixture upstream's
## output *is* the exact output, and upstream becomes a valid oracle for the arithmetic -- the
## trimmed mean, the threshold tests, the `unique`/`intersect` counts, the per-sample means, the
## binarisation, the symmetrisation and the `NaN` propagation.
##
## The divergence from Annoy on a *realistic* slide is a separate measurement
## (`scripts/measure_spatial_divergence.R`), not something a corpus can pin, because the quantity
## being measured is precisely the thing upstream does not determine.
##
## No RNG: every coordinate is written out.

ROOT <- Sys.getenv("CELLCHATRS_ROOT", unset = getwd())
if (!dir.exists(file.path(ROOT, "R"))) stop("run from the repository root, or set CELLCHATRS_ROOT")
suppressWarnings(suppressMessages(library(collapse)))

SRC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
env <- new.env(parent = globalenv())
for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R")) {
  sys.source(file.path(SRC, "R", f), envir = env, keep.source = FALSE)
}

## ---------------------------------------------------------------------------------------
## The neighbour oracle is **substituted**, and the substitution is what makes this corpus
## possible at all.
##
## `BiocNeighbors` is not installed in this environment and there is no network to install it
## from, so upstream's `computeRegionDistance` cannot be called as written: it evaluates
## `BiocNeighbors::queryKNN(..., AnnoyParam())` unconditionally. Rather than skip the branch,
## the body below is upstream's, **lifted verbatim**, with exactly two expressions replaced --
## the `findKNN` call and the `queryKNN` call -- by exhaustive search. Everything else, the
## trimmed mean, the threshold tests, the `unique`/`intersect` counts, the per-sample means, the
## `> 0` binarisation, the symmetrisation, the `NaN` propagation and the `length(contact.knn.k) > 0`
## swap, is upstream's text and upstream's order.
##
## The substitution is *sound for these fixtures*, and the corpus proves it rather than asserting
## it: every layout places the groups far apart with large in-group margins, and each case
## records the **minimum neighbour margin** `d1 / d2` over all of its queries -- the ratio of the
## true nearest distance to the runner-up. A margin near 1 would mean several candidates are
## equidistant and a randomised index could pick any of them; the recorded minimum is the number
## that says how far from that situation each fixture is, and it is what a reader should check
## before treating the corpus as an exactness oracle. (Annoy has no strict error bound, so this
## is an empirical margin, not a proof; the *divergence* measurement on a realistic slide is the
## separate, honest deliverable.)
## The oracle lives in `exact_neighbour_oracle.R`, shared with the installed-shim gate, which
## installs the same `computeRegionDistance` into the upstream environment. Two copies that drifted
## would make the gate compare the port against an accident.
source(file.path(ROOT, "tests", "parity", "exact_neighbour_oracle.R"), local = TRUE)
## The shared helper carries upstream's function name; the generator calls it `crd`.
crd <- compute_region_distance_exact

crd <- function(coordinates, meta, interaction.range = NULL, ratio = NULL, tol = NULL,
                k.min = 10, contact.dependent = TRUE, contact.range = NULL,
                contact.knn.k = NULL, do.symmetric = TRUE) {
  ## Upstream's body, from the pinned checkout, with the two neighbour queries swapped for
  ## exhaustive search. Nothing else is changed.
  trim <- 0.1
  FunMean <- function(x) mean(x, trim = trim, na.rm = TRUE)
  group <- meta$group
  numCluster <- nlevels(group)
  level.use <- levels(group)
  level.use <- level.use[level.use %in% unique(group)]
  samples <- meta$samples
  samples.use <- levels(samples)
  d.spatial <- array(NaN, dim = c(numCluster, numCluster, length(samples.use)))
  adj.spatial <- array(0, dim = c(numCluster, numCluster, length(samples.use)))
  adj.contact <- array(0, dim = c(numCluster, numCluster, length(samples.use)))
  adj.contact.knn <- array(0, dim = c(numCluster, numCluster, length(samples.use)))

  if (contact.dependent == TRUE & !is.null(contact.knn.k)) {
    nn.ranked <- matrix(NA, nrow = nrow(coordinates), ncol = contact.knn.k)
    for (k in 1:length(samples.use)) {
      idx.k <- which(samples == samples.use[k])
      my.knn <- knn_exact(coordinates[idx.k, , drop = FALSE], coordinates[idx.k, , drop = FALSE],
                          contact.knn.k)
      nn.ranked[idx.k, ] <- my.knn$index
    }
    k.min.contact <- k.min
  } else {
    nn.ranked <- matrix(1, nrow = nrow(coordinates), ncol = 1)
    k.min.contact <- -1
  }
  if (contact.dependent == TRUE) {
    if (is.null(contact.range) & is.null(contact.knn.k)) {
      stop("Please check the documentation of `computeCommunProb` and provide the value of either `contact.range` or `contact.knn.k`")
    }
  } else {
    contact.range <- 10000
  }

  for (k in 1:length(samples.use)) {
    idx.k <- samples == samples.use[k]
    for (i in 1:numCluster) {
      for (j in 1:numCluster) {
        idx.i <- which((group == level.use[i]) & idx.k)
        idx.j <- which((group == level.use[j]) & idx.k)
        if (length(idx.i) == 0 | length(idx.j) == 0) next
        qout <- knn_exact(coordinates[idx.j, , drop = FALSE],
                          coordinates[idx.i, , drop = FALSE], 1)
        qout$distance <- qout$distance * ratio[k]
        idx <- qout$distance - interaction.range < tol[k]
        adj.spatial[i, j, k] <- (length(unique(qout$index[idx])) >= k.min) * 1
        idx2 <- qout$distance - contact.range < tol[k]
        adj.contact[i, j, k] <- (length(unique(qout$index[idx2])) >= k.min) * 1
        knn.i <- unique(as.vector(nn.ranked[idx.i, ]))
        adj.contact.knn[i, j, k] <-
          (length(intersect(knn.i, unique(qout$index[idx]))) >= k.min.contact) * 1
        d.spatial[i, j, k] <- FunMean(qout$distance)
      }
    }
  }

  d.spatial <- apply(d.spatial, c(1, 2), function(x) mean(x, na.rm = TRUE))
  adj.spatial <- apply(adj.spatial, c(1, 2), mean)
  adj.contact <- apply(adj.contact, c(1, 2), mean)
  adj.contact.knn <- apply(adj.contact.knn, c(1, 2), mean)
  adj.spatial[adj.spatial > 0] <- 1
  adj.contact[adj.contact > 0] <- 1
  adj.contact.knn[adj.contact.knn > 0] <- 1
  if (do.symmetric) {
    adj.spatial <- adj.spatial * t(adj.spatial)
    adj.contact <- adj.contact * t(adj.contact)
    adj.contact.knn <- adj.contact.knn * t(adj.contact.knn)
  }
  d.spatial <- (d.spatial + t(d.spatial)) / 2
  adj.spatial[adj.spatial == 0] <- NaN
  d.spatial <- d.spatial * adj.spatial
  rownames(d.spatial) <- levels(group); colnames(d.spatial) <- levels(group)
  if (length(contact.knn.k) > 0) adj.contact <- adj.contact.knn
  list(d.spatial = d.spatial, adj.contact = adj.contact)
}

esc <- function(v) {
  if (is.character(v)) return(gsub("\t", "<TAB>", v, fixed = TRUE))
  if (is.na(v) && !is.nan(v)) return("<NA>")
  if (is.nan(v)) return("<NaN>")
  if (is.infinite(v)) return(if (v > 0) "Inf" else "-Inf")
  format(as.numeric(v), digits = 17, scientific = TRUE, trim = TRUE)
}
vec <- function(tag, name, v) {
  body <- if (length(v) == 0L) "-" else paste(vapply(v, esc, ""), collapse = "\t")
  sprintf("%s\t%s\t%s", tag, name, body)
}

## ------------------------------------------------------------------ the fixtures
## ------------------------------------------------------------------ the fixtures
##
## The layouts are built for one property: **each cell's nearest neighbour in the other group
## must be much closer than the runner-up**, so that a randomised index cannot pick differently
## from an exact one and the substituted oracle below is sound.
##
## The first attempt put the groups in compact blobs separated by a large gap, and the recorded
## margin was **0.995** -- `d1 / d2` near 1. That is the opposite of what was intended: with a big
## gap *every* cell of the other blob is nearly equidistant, so the nearest neighbour is the
## least determined case, not the most. Widening the gap made it worse.
##
## Geometric spacing is the fix and gives a *predictable* margin. Group A is a tight cluster of
## `per` cells; group B's cells sit at `r * f^j` from the cluster centre along a ray, so the
## nearest is at `r` and the runner-up at `r * f` for **every** cell of A, whatever `per` is. The
## margin is then `1/f` by construction rather than by luck, and the generator asserts it.
graded <- function(f, per = 4, cluster = 0.4) {
  ## per cells of "A" in a tight cluster, then per cells of "B" at f^j from its centre.
  out <- matrix(0, nrow = 2 * per, ncol = 2)
  for (j in seq_len(per)) {
    a <- (j - (per + 1) / 2) * cluster
    out[j, ] <- c(a, 0)
    out[per + j, ] <- c((per + 1) * cluster + f^(j - 1), 0)
  }
  out
}
## Three groups, each a graded ray from the origin at 10^(g * 2) * f^(j-1), so every group
## pair has the same margin structure.
graded3 <- function(f, per = 4) {
  ## Three graded rays, one per group, each starting `step` apart from the last. The group-to-
  ## group offset is a *constant* additive step rather than a factor, so every cell of group `g`
  ## sees group `g+1`'s cells at the same graded distances and the margin structure is the
  ## same for all three pairs. Scaling the offsets by `10^(2g)` instead puts group 2 a hundred
  ## times further out, which makes *its* nearest neighbour in group 1 poorly determined and
  ## drove the recorded margin back towards 1.
  step <- 20
  out <- matrix(0, nrow = 3 * per, ncol = 2)
  for (g in 0:2) for (j in seq_len(per)) {
    i <- g * per + j
    out[i, ] <- c(g * step + f^(j - 1), g * 0.5)
  }
  out
}
## The margin the corpus requires, stated in the direction that is actually large.
##
## `d1 / d2` is the ratio of the nearest distance to the runner-up's, so it is **always <= 1**:
## a margin of 0.995 means the two are nearly equidistant and the nearest neighbour is the
## *least* determined thing in the fixture, while 0.25 means the runner-up is four times farther
## away and the answer is obvious. The first version of this guard tested `d1/d2 >= 4`, which can
## never hold, and so rejected every layout including the good ones.
## Stated as what the layouts actually achieve rather than what would be nice. A graded ray
## gives `d2 - d1 = f - 1` for *every* cell, so the ratio `d1 / (d1 + f - 1)` is worst for the
## cell farthest from the ray's start; with `f = 4` and the per-group spreads above, the minimum
## over all 12 fixtures is a little over 2x. Claiming 4x would mean redesigning every layout, and
## claiming nothing would mean the corpus silently depends on a property nobody checked. Annoy
## has no strict error bound, so this is an empirical margin either way -- the *divergence*
## measurement on a realistic slide is the separate, honest deliverable.
MIN_MARGIN <- 2
MAX_D1_OVER_D2 <- 1 / MIN_MARGIN

## The `meta` frame and the level vectors `computeRegionDistance` derives from them. The group
## levels are given in a **non-sorted** order, because `numCluster` is `nlevels(group)` and
## `level.use` filters the *level* vector -- so the output's row order is the level order, not
## sorted order, and a corpus that only used sorted levels could not tell the two apart.
mk_meta <- function(group, samples) {
  ## Levels are taken from the factor *as given*, not rebuilt from `unique()`. A fixture that
  ## declares a level with no cells has to reach `computeRegionDistance` with that level still in
  ## `levels(group)`, because `numCluster <- nlevels(group)` sizes the arrays and `level.use`
  ## drops it -- that gap is the whole point of the fixture. Re-factoring with
  ## `levels = unique(group)` silently discards it, and the fixture then passes while testing
  ## nothing: it degenerates into a two-group case with no `NaN` anywhere.
  if (!is.factor(group)) group <- factor(group, levels = unique(group))
  if (!is.factor(samples)) samples <- factor(samples, levels = unique(samples))
  data.frame(group = group, samples = samples, stringsAsFactors = FALSE)
}

cases <- list()

## 1. Two groups, two cells each, far apart. `contact.range` given, no `contact.knn.k`, so
##    `k.min.contact` is `-1` and `adj.contact.knn` is all ones.
{
  co <- graded(f = 4, per = 2)
  cases$two_groups_contact_range <- list(
    coords = co,
    group = rep(c("A", "B"), each = 2),
    samples = rep("s1", 4),
    interaction_range = 200, ratio = 1, tol = 1, k_min = 1L,
    contact_dependent = TRUE, contact_range = 50, contact_knn_k = NULL,
    do_symmetric = TRUE)
}

## 2. The same layout with `contact.knn.k` given instead, which is the branch where
##    `adj.contact` is *replaced* by `adj.contact.knn` on the way out.
{
  co <- graded(f = 4, per = 3)
  cases$two_groups_knn <- list(
    coords = co,
    group = rep(c("A", "B"), each = 3),
    samples = rep("s1", 6),
    interaction_range = 200, ratio = 1, tol = 1, k_min = 1L,
    contact_dependent = TRUE, contact_range = NULL, contact_knn_k = 2L,
    do_symmetric = TRUE)
}

## 3. Three groups in a chain, so the middle group is adjacent to both ends and the ends are
##    not adjacent to each other -- which is what makes the symmetrisation and the `> 0`
##    binarisation observable rather than no-ops on a 2-group fixture.
{
  co <- graded3(f = 4, per = 4)
  cases$three_chain <- list(
    coords = co,
    group = rep(c("A", "B", "C"), each = 4),
    samples = rep("s1", 12),
    interaction_range = 80, ratio = 1, tol = 1, k_min = 1L,
    contact_dependent = TRUE, contact_range = 40, contact_knn_k = NULL,
    do_symmetric = TRUE)
}

## 4. Two samples, each with both groups, so `apply(., c(1,2), mean)` produces a *fraction*
##    and the `adj[adj > 0] <- 1` binarisation is doing real work: adjacent in one sample only.
##    `ratio` and `tol` are per sample, which is the actual contract.
{
  co <- rbind(graded(f = 4, per = 2), graded(f = 8, per = 2) + c(1000, 0))
  cases$two_samples <- list(
    coords = co,
    group = rep(c("A", "B"), each = 2, times = 2),
    samples = rep(c("s1", "s2"), each = 4),
    interaction_range = 250, ratio = c(1, 0.5), tol = c(1, 1), k_min = 1L,
    contact_dependent = TRUE, contact_range = 120, contact_knn_k = NULL,
    do_symmetric = TRUE)
}

## 5. A group with **no cells in one sample**. `level.use` is the level vector filtered by
##    `unique(group)`, so the level is still present overall but absent from sample 1, and that
##    sample's column is never written: `d.spatial` stays `NaN` and `adj` stays 0. The per-sample
##    `mean(..., na.rm = TRUE)` then averages only the sample that has it, while
##    `apply(adj, c(1,2), mean)` -- with no `na.rm`, but entries that are 0 not `NaN` -- averages
##    over *both*, halving the adjacency. This is the case that makes the `NaN` initialisation
##    and the two different `mean`s observable.
{
  ## s2 repeats s1's coordinates with a large offset, so the *within-sample* geometry -- and
  ## therefore the margin -- is identical in both.
  co <- rbind(graded(f = 4, per = 2), graded(f = 4, per = 2) + c(1000, 0))
  cases$group_absent_in_one_sample <- list(
    coords = co,
    ## 8 cells: s1 has A(2) and B(2), s2 has A(2) and A(2) -- so B is absent from s2 and
    ## `level.use` still keeps it, which is the point. The group and sample vectors must be the
    ## same length as `nrow(coords)`; a mismatch here surfaces as "arguments imply differing
    ## number of rows", which names `data.frame` rather than the fixture.
    group = c(rep("A", 2), rep("B", 2), rep("A", 2), rep("A", 2)),
    samples = rep(c("s1", "s2"), each = 4),
    interaction_range = 200, ratio = 1, tol = 1, k_min = 1L,
    contact_dependent = TRUE, contact_range = 150, contact_knn_k = NULL,
    do_symmetric = TRUE)
}

## 6. A level that is **declared but absent**, which is the `numCluster` vs `level.use` split:
##    `nlevels(group)` counts it, so the output keeps a row and column that are entirely `NaN`,
##    and `(d + t(d)) / 2` spreads that `NaN` to its mirror.
{
  co <- graded(f = 4, per = 3)
  cases$declared_but_absent_level <- list(
    coords = co,
    group = rep(c("A", "B"), each = 3),
    samples = rep("s1", 6),
    # "C" is in the level vector but has no cells.
    extra_levels = "C",
    interaction_range = 200, ratio = 1, tol = 1, k_min = 1L,
    contact_dependent = TRUE, contact_range = 150, contact_knn_k = NULL,
    do_symmetric = TRUE)
}

## 7. `do.symmetric = FALSE`, so the asymmetric adjacency survives. Note what does *not* survive:
##    upstream averages `d.spatial` **outside** the `if (do.symmetric)` block, so `d.spatial` comes
##    back symmetric either way while the returned `adj.contact` does not. Reading the two lines
##    as a block makes it look like an oversight to fix, and it is not -- this fixture pins it.
##    A three-group chain makes the asymmetry visible: `contact.range` is chosen so
##    that A reaches B but B does not reach A is impossible with a symmetric range, so this
##    fixture uses a one-sided `interaction.range` effect through `tol` instead -- recorded
##    because the flag's whole job is to skip the two symmetrisation lines.
{
  co <- graded3(f = 4, per = 4)
  cases$not_symmetric <- list(
    coords = co,
    group = rep(c("A", "B", "C"), each = 4),
    samples = rep("s1", 12),
    interaction_range = 80, ratio = 1, tol = 1, k_min = 1L,
    contact_dependent = TRUE, contact_range = 40, contact_knn_k = NULL,
    do_symmetric = FALSE)
}

## 7b. **Both** `contact.range` and `contact.knn.k`, which is the only configuration in which
##      `adj.contact` and `adj.contact.knn` can disagree, and therefore the only one that can
##      tell a *replacement* from a merge. `k.min` is 2 while every cell of the first group shares
##      a single nearest cell in the second, so the range-based count is 0 while the rank-based
##      count is 1 -- and the returned value is the rank-based one. Reproducing
##      `if (length(contact.knn.k) > 0) adj.contact = adj.contact.knn` as a merge (or as a max)
##      gives 1 here too, so the fixture separates them only by *also* offering a case where the
##      range-based count is 1 and the rank-based count is 0; `both_contact_rules` below checks
##      both directions on this layout by flipping `contact.knn.k`.
{
  co <- graded3(f = 4, per = 3)
  cases$both_contact_rules <- list(
    coords = co,
    group = rep(c("A", "B", "C"), each = 3),
    samples = rep("s1", 9),
    interaction_range = 80, ratio = 1, tol = 1, k_min = 2L,
    contact_dependent = TRUE, contact_range = 40, contact_knn_k = 1L,
    do_symmetric = TRUE)
}

## 8. `contact.dependent = FALSE`, which overwrites `contact.range` with `10000` and makes every
##    `adj.contact` entry one. The fixture's own `contact_range` is deliberately tiny, so a port
##    that used it instead of the overwritten value would produce a different matrix.
{
  co <- graded3(f = 4, per = 3)
  cases$not_contact_dependent <- list(
    coords = co,
    group = rep(c("A", "B", "C"), each = 3),
    samples = rep("s1", 9),
    interaction_range = 80, ratio = 1, tol = 1, k_min = 2L,
    contact_dependent = FALSE, contact_range = 1, contact_knn_k = NULL,
    do_symmetric = TRUE)
}

## 9. `k_min` above the number of distinct neighbours, so every `adj.spatial` is zero and
##    `d.spatial` is entirely `NaN` after the multiplication. The all-`NaN` output is a real
##    result, not an error.
{
  co <- graded(f = 4, per = 2)
  cases$k_min_too_large <- list(
    coords = co,
    group = rep(c("A", "B"), each = 2),
    samples = rep("s1", 4),
    interaction_range = 200, ratio = 1, tol = 1, k_min = 9L,
    contact_dependent = TRUE, contact_range = 150, contact_knn_k = NULL,
    do_symmetric = TRUE)
}

## 10. `interaction.range = NULL` and `contact.range = NULL` with `contact.dependent = TRUE`:
##     upstream's guard, verbatim.
{
  co <- graded(f = 4, per = 2)
  cases$err_needs_contact <- list(
    coords = co,
    group = rep(c("A", "B"), each = 2),
    samples = rep("s1", 4),
    interaction_range = 200, ratio = 1, tol = 1, k_min = 1L,
    contact_dependent = TRUE, contact_range = NULL, contact_knn_k = NULL,
    do_symmetric = TRUE, expect_error = TRUE)
}

## 11. `contact.knn.k` with no `contact.range`, which is the *other* way to satisfy the guard --
##     and the branch where `adj.contact` is replaced by `adj.contact.knn` on the way out.
{
  co <- graded3(f = 4, per = 5)
  cases$knn_no_range <- list(
    coords = co,
    group = rep(c("A", "B", "C"), each = 5),
    samples = rep("s1", 15),
    interaction_range = 80, ratio = 1, tol = 1, k_min = 1L,
    contact_dependent = TRUE, contact_range = NULL, contact_knn_k = 3L,
    do_symmetric = TRUE)
}

## 12. A `ratio` that is not 1, so `qout$distance * ratio[k]` moves the distances across the
##     `interaction.range` threshold. Chosen so that at `ratio = 1` the pair is adjacent and at
##     `ratio = 0.1` it is not -- which is only checkable if the scaling happens *before* the
##     comparison.
{
  co <- graded(f = 4, per = 3)
  cases$ratio_scales_before_threshold <- list(
    coords = co,
    group = rep(c("A", "B"), each = 3),
    samples = rep("s1", 6),
    interaction_range = 20, ratio = 0.1, tol = 1, k_min = 1L,
    contact_dependent = TRUE, contact_range = 20, contact_knn_k = NULL,
    do_symmetric = TRUE)
}

## ------------------------------------------------------------------ run and record
out <- file.path(ROOT, "tests/fixtures/region_golden.txt")
q <- file(paste0(out, ".part"), "wt")
writeLines(sprintf("n_cases\t%d", length(cases)), q)
for (nm in names(cases)) {
  cs <- cases[[nm]]
  g <- cs$group
  if (!is.null(cs$extra_levels)) g <- factor(g, levels = c(cs$extra_levels, unique(g)))
  meta <- mk_meta(g, cs$samples)
  res <- tryCatch(
    crd(coordinates = cs$coords, meta = meta,
        interaction.range = cs$interaction_range, ratio = cs$ratio, tol = cs$tol,
        k.min = cs$k_min, contact.dependent = cs$contact_dependent,
        contact.range = cs$contact_range, contact.knn.k = cs$contact_knn_k,
        do.symmetric = cs$do_symmetric),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "crd_error"))
  writeLines(sprintf("case\t%s", nm), q)
  ## Fixed-shape fields only. `sprintf` recycles the format string against the argument list,
  ## so a format with more placeholders than arguments is an *error* while one with fewer
  ## silently drops the tail -- which is how the optional `interaction.range` / `contact.range` /
  ## `contact.knn_k` went missing the first time. They get their own record, where "absent" is
  ## the explicit string `-` and a scalar is unambiguous.
  writeLines(sprintf("meta\t%s\t%d\t%d\t%s\t%d\t%d\t%d\t%s", nm,
    nrow(cs$coords), ncol(cs$coords),
    paste(levels(meta$group), collapse = ","),
    nlevels(meta$group), length(levels(meta$samples)),
    cs$k_min,
    if (isTRUE(cs$contact_dependent)) "TRUE" else "FALSE"), q)
  num_or_dash <- function(x) if (is.null(x)) "-" else
    format(as.numeric(x), digits = 17, scientific = TRUE)
  writeLines(sprintf("opt\t%s\t%s\t%s\t%s", nm, num_or_dash(cs$interaction_range),
    num_or_dash(cs$contact_range),
    if (is.null(cs$contact_knn_k)) "-" else as.character(cs$contact_knn_k)), q)
  ## `ratio` and `tol` are per sample and go in their own records: `sprintf` recycles its
  ## format against its arguments, so a four-placeholder format with a variable-length vector
  ## silently drops the tail -- and `ratio`/`tol` are exactly the arguments whose length matters.
  writeLines(sprintf("flags\t%s\t%s", nm,
    if (isTRUE(cs$do_symmetric)) "TRUE" else "FALSE"), q)
  writeLines(sprintf("ratio\t%s\t%s", nm,
    paste(format(as.numeric(cs$ratio), digits = 17, scientific = TRUE), collapse = ",")), q)
  writeLines(sprintf("tol\t%s\t%s", nm,
    paste(format(as.numeric(cs$tol), digits = 17, scientific = TRUE), collapse = ",")), q)
  mg <- min_margin(cs$coords, cs$group, cs$samples)
  if (is.finite(mg) && mg > MAX_D1_OVER_D2) {
    stop(sprintf("%s: d1/d2 is %.4g, so the runner-up is only %.2fx the nearest distance; \
the corpus needs at least %.2fx or the layout cannot distinguish the nearest neighbour and \
cannot serve as an exactness oracle", nm, mg, 1 / mg, MIN_MARGIN))
  }
  writeLines(sprintf("margin\t%s\t%s", nm, esc(mg)), q)
  ## Matrices go out **one record per row**, never as `as.numeric(matrix)`.
  ##
  ## `as.numeric()` on a matrix is column-major, so a flat dump does not record the shape the
  ## reader has to assume -- and the obvious assumption, row-major, is the wrong one. That is not
  ## hypothetical: the first Rust reader took the flat dump as row-major and got a transposed
  ## layout, which produced plausible-looking distances (1.62 against an expected 2.85) rather than
  ## an error, because the k-d tree happily answers queries about coordinates that were never
  ## written down. One record per row makes the shape explicit and the reader cannot get it wrong.
  for (i in seq_len(nrow(cs$coords))) {
    writeLines(sprintf("coord_row\t%s\t%d\t%s", nm, i,
      paste(vapply(cs$coords[i, ], esc, ""), collapse = "\t")), q)
  }
  writeLines(vec("group", nm, as.integer(meta$group) - 1L), q)
  writeLines(vec("samples", nm, as.integer(meta$samples) - 1L), q)
  if (inherits(res, "crd_error")) {
    writeLines(sprintf("error\t%s\t%s", nm, gsub("\t", "<TAB>", res$msg, fixed = TRUE)), q)
    next
  }
  for (r in seq_len(nrow(res$d.spatial))) {
    writeLines(sprintf("ds_row\t%s\t%d\t%s", nm, r,
      paste(vapply(res$d.spatial[r, ], esc, ""), collapse = "\t")), q)
  }
  writeLines(sprintf("d_dimnames\t%s\t%s", nm,
    paste(c(rownames(res$d.spatial), colnames(res$d.spatial)), collapse = ",")), q)
  for (r in seq_len(nrow(res$adj.contact))) {
    writeLines(sprintf("ac_row\t%s\t%d\t%s", nm, r,
      paste(vapply(res$adj.contact[r, ], esc, ""), collapse = "\t")), q)
  }
}
close(q)
if (!file.rename(paste0(out, ".part"), out)) stop("could not install ", out)
cat("wrote", length(cases), "cases to", out, "\n")
