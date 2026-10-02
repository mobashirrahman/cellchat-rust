# Golden corpus for the deterministic half of `computeCentralityLocal`, from igraph itself.
#
# `crates/r-core/src/centrality.rs` reimplements igraph's strength summation order (edges in
# edge-ID order, plain sequential `f64`), its Dijkstra (dist-plus-one encoding, exact 2-way-heap
# tie rules, epsilon comparisons at 1e-10) and its Brandes accumulation, all read out of the C
# sources. This file is the evidence that the reimplementation is exact: every record below is what
# installed igraph 2.3.4 produced, and the Rust tests require bit-for-bit agreement.
#
# Two things this corpus deliberately does not do, both of which would make it weaker:
#
#   * It does not cover `hub_score` / `authority_score` / `eigen_centrality` / `page_rank`. Those
#     are iterative solvers whose output differs run to run on identical input (measured:
#     `hub_score` and `eigen_centrality` both disagree with themselves across two calls with no
#     seed set). A golden for a nondeterministic oracle pins nothing.
#   * It does not avoid ties, zeros, loops, or near-ties. The whole point of replicating the heap
#     tie rules and the epsilon comparison is that ties are where implementations diverge, so the
#     corpus is built to be full of them: duplicate points in `runif` discretisation, exact ties
#     by construction, magnitudes chosen so sequential and compensated summation disagree.
#
# Usage: R_LIBS=.rlib R --vanilla -f tests/parity/gen_centrality_golden.R
suppressWarnings(suppressMessages({ library(Matrix); library(igraph) }))

OUT <- "tests/fixtures/centrality_golden.txt"
con <- file(OUT, "wt")

fmt <- function(x) sprintf("%.17g", x)

## One case: the matrix, then everything igraph says about it. Values with full precision, so the
## Rust side compares bits rather than decimals.
emit <- function(name, m) {
  k <- nrow(m)
  cat(sprintf("case\t%s\t%d\n", name, k), file = con)
  cat("matrix\n", file = con)
  writeLines(apply(m, 1, function(r) paste(fmt(r), collapse = " ")), con)
  ## Error text, one record per line. igraph's weight messages arrive as two lines
  ## ("...Invalid value\nSource: centrality/betweenness.c:439") and the construction error as
  ## one; the `Source:` suffix is part of `conditionMessage()`, so it is recorded rather than
  ## stripped, and split across `error` / `error2` so every record stays single-line.
  write_error <- function(msg) {
    parts <- strsplit(msg, "\n", fixed = TRUE)[[1]]
    cat(sprintf("error\t%s\n", parts[1]), file = con)
    if (length(parts) > 1) {
      cat(sprintf("error2\t%s\n", paste(parts[-1], collapse = " ")), file = con)
    }
  }
  ## Construction first: `graph_from_adjacency_matrix` raises on any NA/NaN before any edge
  ## exists, so degrees, strength and betweenness are all unreached and the case records only
  ## the construction error. Upstream's `computeCentralityLocal` builds the graph before
  ## computing anything else, so the port must fail at the same point with the same text.
  G <- tryCatch(
    suppressWarnings(igraph::graph_from_adjacency_matrix(m, mode = "directed",
                                                        weighted = TRUE)),
    error = function(e) structure(conditionMessage(e), class = "centrality_error"))
  if (inherits(G, "centrality_error")) {
    write_error(G)
    return(invisible(NULL))
  }
  cat(sprintf("nedges\t%d\n", length(igraph::E(G))), file = con)
  ## Degrees and strength always succeed on a constructed graph (no reciprocal involved), so they
  ## are recorded unconditionally -- including for matrices whose betweenness raises below. Splitting
  ## construction errors (both entry points fail) from reciprocal errors (only betweenness fails)
  ## is what lets the Rust tests assert each entry point separately.
  for (nm in c("outdeg_unweighted", "indeg_unweighted")) {
    v <- if (nm == "outdeg_unweighted") as.numeric(rowSums(m > 0)) else as.numeric(colSums(m > 0))
    cat(sprintf("%s\n", nm), file = con)
    writeLines(paste(fmt(v), collapse = " "), con)
  }
  for (nm in c("outdeg", "indeg")) {
    v <- if (nm == "outdeg") igraph::strength(G, mode = "out") else igraph::strength(G, mode = "in")
    cat(sprintf("%s\n", nm), file = con)
    writeLines(paste(fmt(v), collapse = " "), con)
  }
  ## A warning is not a failure: the value is still computed and still recorded, and a bare
  ## `warning` marker records that it fired. `invokeRestart("muffleWarning")` only exists inside
  ## `withCallingHandlers`, not inside `tryCatch`'s `warning=` handler -- calling it there raises
  ## "no 'restart' 'muffleWarning' found", which an earlier version of this generator recorded as
  ## fifteen bogus `betweenness_error` entries with that text as the "message". A handler that
  ## turns warnings into errors is worse than no handler, because the corpus then pins the
  ## handler's failure instead of igraph's behaviour.
  warned <- FALSE
  btw <- withCallingHandlers(
    tryCatch({
      igraph::E(G)$weight <- 1 / igraph::E(G)$weight
      igraph::betweenness(G)
    }, error = function(e) structure(conditionMessage(e), class = "centrality_error")),
    warning = function(w) {
      warned <<- TRUE
      invokeRestart("muffleWarning")
    })
  if (inherits(btw, "centrality_error")) {
    cat("betweenness_error\n", file = con)
    write_error_lines(btw)
  } else {
    if (warned) {
      cat("warning\n", file = con)
    }
    cat("betweenness\n", file = con)
    writeLines(paste(fmt(btw), collapse = " "), con)
  }
  invisible(NULL)
}

## `write_error` writes `error`/`error2`; this writes a bare message already split the same way,
## for the `betweenness_error` marker which carries no message of its own.
write_error_lines <- function(msg) {
  parts <- strsplit(msg, "\n", fixed = TRUE)[[1]]
  cat(sprintf("error\t%s\n", parts[1]), file = con)
  if (length(parts) > 1) {
    cat(sprintf("error2\t%s\n", paste(parts[-1], collapse = " ")), file = con)
  }
}

## ---------------------------------------------------------------- randomized graphs
## Small, dense-ish, and adversarial by construction: a quarter of entries zeroed (zero-dropping),
## values spanning 1e-6..1e6 in signed pairs (summation-order sensitivity), exact ties (heap ties),
## and diagonals left nonzero half the time (loops counted once in strength, excluded from traversal).
## Non-negative throughout, like the probabilities this port actually sees (`Prob >= 0` is asserted
## on the tutorial's output). An earlier version drew `runif(-1, 1)`, so 60 of 75 cases died in
## igraph's reciprocal-positivity check and the corpus pinned the error text sixty times instead of
## any arithmetic. Negatives survive only in the dedicated error cases below, where the error *is*
## the assertion.
##
## Summation-order sensitivity without signs comes from magnitudes: values spanning 1e-6..1e6 in the
## same row make sequential, compensated and pairwise summation disagree, which is what pins the
## edge-ID order rather than the fact of addition.
set.seed(20240502)
ncases <- 0L
for (rep in 1:60) {
  k <- sample(2:8, 1)
  m <- matrix(runif(k * k) * 10^sample(-6:6, k * k, replace = TRUE), k, k)
  m[runif(k * k) < 0.25] <- 0
  if (rep %% 3 == 0) {
    ## Exact ties: duplicate a row's nonzero pattern into another row, so two vertices have
    ## identical out-neighbourhoods and the heap must break ties between them.
    i <- sample(k, 1); j <- sample(setdiff(seq_len(k), i), 1)
    m[j, ] <- m[i, ]
  }
  if (rep %% 4 == 0) {
    ## Magnitudes engineered so summation *order* is observable without signs: a 1e16 entry beside
    ## unit-scale ones, where sequential accumulation in different orders rounds differently and only
    ## the exact edge-ID order matches. (Signed cancellation would exercise the same property, but
    ## signed entries make `betweenness` raise on the reciprocal check instead of computing, so the
    ## arithmetic would go untested. One signed case lives in the error section below, where the
    ## raised message is the assertion.)
    m[] <- round(m, 1)
    m[sample(k * k, 1)] <- 1e16
  }
  emit(sprintf("rand%02d", rep), m)
  ncases <- ncases + 1L
}

## ---------------------------------------------------------------- reciprocal edge cases
## `betweenness` runs on `1/weight`, so tiny probabilities become huge distances and huge
## probabilities become tiny ones. These cases put reciprocals near the 1e-10 warning boundary,
## at zero (1/Inf), and at infinity (1/0 is impossible -- zeros are dropped -- but 1/tiny is not).
emit("recip_tiny", matrix(c(1e-11, 0.5, 0.3, 0.2), 2, 2, byrow = TRUE))
emit("recip_huge", matrix(c(1e10, 0.5, 0.3, 0.2), 2, 2, byrow = TRUE))
emit("recip_inf", matrix(c(Inf, 0.5, 0.3, 0.2), 2, 2, byrow = TRUE))
emit("allzero", matrix(0, 3, 3))
emit("single", matrix(0.5, 1, 1))
emit("two_disconnected", matrix(c(0.5, 0, 0, 0.25), 2, 2, byrow = TRUE))
emit("selfloop_only", matrix(c(0.7, 0, 0, 0), 2, 2, byrow = TRUE))
ncases <- ncases + 7L

## ---------------------------------------------------------------- error cases
## Whatever makes `graph_from_adjacency_matrix` raise must make the port raise the identical text,
## and whatever makes the reciprocal checks raise must do the same -- while degrees and strength on
## the same matrices still succeed, because construction succeeded. That split is why `emit`
## records measures and betweenness separately: a single early-return error record could not express
## "degrees fine, betweenness raises".
mna <- matrix(c(0.5, NA, 0.25, 0.1), 2, 2, byrow = TRUE)
emit("na_matrix", mna)
mnan <- matrix(c(0.5, NaN, 0.25, 0.1), 2, 2, byrow = TRUE)
emit("nan_matrix", mnan)
## Negative entries survive construction (only NA/NaN raise there) and die in the reciprocal
## positivity check -- so these pin both the error text and the fact that degrees/strength are
## unaffected by it.
mneg <- matrix(c(0.5, -0.1, 0.25, 0.1), 2, 2, byrow = TRUE)
emit("neg_weight", mneg)
mneg2 <- matrix(c(0.5, 0.2, 0.25, -1e-300), 2, 2, byrow = TRUE)
emit("neg_tiny_weight", mneg2)
ncases <- ncases + 4L

## ---------------------------------------------------------------- real pathway slices
## The tutorial's `netP$prob` slices, so the corpus includes the shapes and value ranges the
## comparison analysis actually feeds centrality (12 groups, sparse, probabilities with heavy
## tails). Regenerated from the tutorial object, not hand-made: hand-made fixtures test what the
## author imagined, and the tutorial tests what exists.
if (file.exists("tests/fixtures/tutorial_netp_slice.tsv")) {
  lines <- readLines("tests/fixtures/tutorial_netp_slice.tsv")
  k <- as.integer(strsplit(lines[1], "\t", fixed = TRUE)[[1]][2])
  npath <- as.integer(strsplit(lines[2], "\t", fixed = TRUE)[[1]][2])
  vals <- as.numeric(strsplit(lines[3], " ", fixed = TRUE)[[1]])
  for (p in seq_len(npath)) {
    m <- matrix(vals[((p - 1) * k * k + 1):(p * k * k)], k, k)
    emit(sprintf("tutorial_p%02d", p), m)
    ncases <- ncases + 1L
  }
}
close(con)
cat(sprintf("wrote %s (%d cases)\n", OUT, ncases), file = stderr())
