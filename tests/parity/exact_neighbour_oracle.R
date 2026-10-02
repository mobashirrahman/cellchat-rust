## The exact-neighbour oracle for `computeRegionDistance`, shared by the corpus generator and the
## installed-shim gate.
##
## Upstream's `computeRegionDistance` cannot be called in this environment: it evaluates
## `BiocNeighbors::queryKNN(..., AnnoyParam())` unconditionally and `BiocNeighbors` is not
## installable without network. Annoy is also randomised and approximate, so upstream is not an
## oracle for a neighbour query even when it *can* run -- two runs can disagree.
##
## So this file holds upstream's `computeRegionDistance` **body**, lifted verbatim, with exactly two
## expressions replaced -- the `findKNN` call and the `queryKNN` call -- by exhaustive search.
## Everything else is upstream's text in upstream's order: the trimmed mean, the threshold tests,
## the `unique`/`intersect` counts, the per-sample means, the `> 0` binarisation, the
## symmetrisation, the `NaN` propagation, and the `length(contact.knn.k) > 0` swap.
##
## Two consumers, for two different jobs:
##
##   * `gen_region_golden.R` runs it over the layouts to produce the fixture corpus, recording each
##     layout's neighbour margin so the substitution can be *checked* rather than asserted.
##   * `check_identical.R` installs it into the pinned-upstream environment as
##     `computeRegionDistance`, so the spatial branch of upstream's `computeCommunProb` can run at
##     all and be compared against the Rust shim end to end.
##
## Sharing it is what makes the second use honest. Two copies that drifted would give the generator
## and the gate different functions with the same name, and the gate would be comparing the port
## against an accident.

knn_exact <- function(x, q, k) {
  ## Rows of `q` against rows of `x`, exhaustive. Returns a list with `index` (1-based, into `x`)
  ## and `distance`, shaped like `BiocNeighbors`' return value.
  d <- as.matrix(outer(seq_len(nrow(q)), seq_len(nrow(x)), Vectorize(function(a, b)
    sqrt(sum((q[a, ] - x[b, ])^2)))))
  ord <- t(apply(d, 1, order))
  idx <- matrix(0L, nrow(q), k)
  dst <- matrix(0, nrow(q), k)
  for (a in seq_len(nrow(q))) {
    o <- ord[a, seq_len(min(k, nrow(x)))]
    idx[a, seq_along(o)] <- o
    dst[a, seq_along(o)] <- d[a, o]
  }
  ## `order` is 1-based here already; keep the pairs and the runner-up ratio for the margin.
  attr(ord, "full") <- d
  list(index = idx, distance = dst, ord = ord)
}

## The minimum `d1 / d2` over every query in the layout, a measure of how unambiguous the
## nearest neighbour is. `Inf` when a group has a single cell (no runner-up).
min_margin <- function(coords, group, samples) {
  best <- Inf
  for (s in unique(samples)) {
    for (ga in unique(group)) for (gb in unique(group)) {
      if (ga == gb) next
      qi <- which(group == ga & samples == s)
      qj <- which(group == gb & samples == s)
      if (!length(qi) || !length(qj)) next
      for (i in qi) {
        dd <- sqrt(rowSums((matrix(coords[qj, ], ncol = ncol(coords)) -
                              matrix(coords[i, ], nrow = length(qj), ncol = ncol(coords),
                                     byrow = TRUE))^2))
        dd <- sort(dd)
        if (length(dd) >= 2 && dd[1] > 0) best <- min(best, dd[1] / dd[2])
      }
    }
  }
  best
}


compute_region_distance_exact <- function(coordinates, meta, interaction.range = NULL, ratio = NULL, tol = NULL,
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

