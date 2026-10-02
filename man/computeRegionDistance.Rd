\name{computeRegionDistance}
\alias{computeRegionDistance}
\title{`computeRegionDistance` from pinned upstream, with `BiocNeighbors::queryKNN(..., AnnoyParam())`}
\description{
replaced by an **exact** k-d tree.

Annoy is randomised and approximate, so this is the one place the port is deliberately not
bit-identical to upstream: it removes an approximation rather than introducing one. The
arithmetic is upstream's, in upstream's order, including six behaviours that read as bugs -- most
consequentially that `d.spatial` is symmetrised outside the `if (do.symmetric)` guard, and that a
declared level with no cells receives another level's row name. `docs/SEMANTICS.md` records each
with the fixture that pins it.

The R side owns `levels(factor)` and the sample order, because upstream sizes its arrays with
`nlevels(group)` while filtering the level *names*, and the caller has to be able to see that.
}
\usage{
computeRegionDistance(coordinates, meta, interaction.range = NULL, ratio = NULL, tol = NULL, k.min = 10, contact.dependent = TRUE, contact.range = NULL, contact.knn.k = NULL, do.symmetric = TRUE)
}
\arguments{
\item{coordinates}{Argument `coordinates` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{meta}{Argument `meta` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{interaction.range}{Argument `interaction.range` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{ratio}{Argument `ratio` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{tol}{Argument `tol` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{k.min}{Argument `k.min` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{contact.dependent}{Argument `contact.dependent` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{contact.range}{Argument `contact.range` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{contact.knn.k}{Argument `contact.knn.k` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{do.symmetric}{Argument `do.symmetric` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
