\name{netAnalysis_computeCentrality}
\alias{netAnalysis_computeCentrality}
\title{`netAnalysis_computeCentrality`, with the deterministic measures in Rust}
\description{
Upstream's `computeCentralityLocal` computes eleven measures per pathway network. Four are
iterative solvers whose output is not stable across runs -- `hub_score` (deprecated in igraph
2.0.3), `authority_score`, `eigen_centrality` (ARPACK) and `page_rank` (PRPACK) all disagree
with themselves on identical input, measured with no seed set. Bit-parity with a
nondeterministic oracle is meaningless, so those four stay in R and call igraph directly --
the same package, the same function, the same input, which is identical by construction rather
than by reimplementation.

The other seven are pure functions of the input and live in `r-core::centrality`, verified
bit-for-bit against installed igraph 2.3.4 on a 77-case corpus (`centrality_parity`):
unweighted in/out degrees, weighted strengths (plain sequential `f64` accumulation in edge-ID
order -- *not* R's long-double `rowSums`, which agrees only 3802/4669), and betweenness
(Dijkstra with dist-plus-one encoding, exact 2-way-heap tie rules, epsilon comparisons at
1e-10, Brandes accumulation). `flowbet`/`infocent` are `sna` calls in upstream's own
`tryCatch`, copied verbatim including the zeros fallback: `sna` is installed here, so both
sides run the real functions and agree on nonzero values; a host without `sna` gets zeros
from both sides through the same `tryCatch`.

The preamble (`pval[prob == 0] <- 1`, `prob[pval >= thresh] <- 0`, `net.name`, the
`pbapply`/`future.apply` dispatch, the `names(centr.all)` assignment and the slot writeback)
is upstream verbatim, because none of it is arithmetic the kernel owns. `igraph`'s
`hub_score` deprecation warning fires here exactly as often as upstream -- same function,
same lifecycle state -- so warning behaviour matches by construction rather than by silencing.
}
\usage{
netAnalysis_computeCentrality(object = NULL, slot.name = "netP", net = NULL, net.name = NULL, thresh = 0.05)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{slot.name}{Argument `slot.name` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{net}{Argument `net` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{net.name}{Argument `net.name` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{thresh}{Argument `thresh` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
