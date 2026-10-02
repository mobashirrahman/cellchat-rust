\name{computeCommunProbPathway}
\alias{computeCommunProbPathway}
\title{`computeCommunProbPathway`, backed by the Rust kernel}
\description{
The aggregation, the LONG_DOUBLE sums, the `!= 0` selection and the decreasing-total
reordering are all in `r-core`; only the two `return.object` shapes and the S4 slot
writes are here.

## One upstream failure is preserved, not fixed

```r
prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))
```

`apply(prob, c(1, 2), by, group, sum)` returns a **2-d matrix** when `group` has a single
level and a **3-d array** otherwise, so `aperm(x, c(2, 3, 1))` raises
`'perm' is of wrong length 3 (!= 2)` for a one-pathway `LRsig`. That is reachable from
ordinary use -- an `LRsig` where every interaction belongs to the same pathway -- and a
drop-in replacement has to fail identically, so this case routes to upstream. The Rust
core computes a value for it; it is the shim that must not claim it. Pinned by
`crates/r-core/tests/pathway_parity.rs::single_pathway_errors_like_upstream`.
}
\usage{
computeCommunProbPathway(object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{net}{Argument `net` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{pairLR.use}{Argument `pairLR.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{thresh}{Argument `thresh` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
