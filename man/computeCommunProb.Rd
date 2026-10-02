\name{computeCommunProb}
\alias{computeCommunProb}
\title{`computeCommunProb`, backed by the Rust kernel}
\description{
The whole numeric surface is in `r-core`: the group means, the `i x b` bootstrap
permutation test, the Hill functions, the p-value counting and the `P.spatial`
loop-carried mutation. See the file header for why the shim does so little, and
`docs/SEMANTICS.md` for the R-isms that are reproduced verbatim.

`datatype = "spatial"` routes to the pinned upstream body; see `NOT_PORTED`.

`raw.use = FALSE` is Rust-backed. It reads `object@data.smooth` rather than
`object@data.signaling` and then takes exactly the same path, because upstream's own
`if (raw.use) ... else ...` is two lines that choose a matrix and there is no separate kernel
afterwards. What this package does *not* implement is `projectData`, which is what produces
`data.smooth` in the first place -- so the caller has to supply the slot, and the corpus
exercises the path with an independently constructed one.
}
\usage{
computeCommunProb(object, type = c("triMean", "truncatedMean", "thresholdedMean", "median"), trim = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE, distance.use = TRUE, interaction.range = 250, scale.distance = 0.01, k.min = 10, contact.dependent = TRUE, contact.range = NULL, contact.knn.k = NULL, contact.dependent.forced = FALSE, do.symmetric = TRUE, nboot = 100, seed.use = 1L, Kh = 0.5, n = 1)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{type}{Argument `type` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{trim}{Argument `trim` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{LR.use}{Argument `LR.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{raw.use}{Argument `raw.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{population.size}{Argument `population.size` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{distance.use}{Argument `distance.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{interaction.range}{Argument `interaction.range` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{scale.distance}{Argument `scale.distance` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{k.min}{Argument `k.min` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{contact.dependent}{Argument `contact.dependent` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{contact.range}{Argument `contact.range` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{contact.knn.k}{Argument `contact.knn.k` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{contact.dependent.forced}{Argument `contact.dependent.forced` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{do.symmetric}{Argument `do.symmetric` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{nboot}{Argument `nboot` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{seed.use}{Argument `seed.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{Kh}{Argument `Kh` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{n}{Argument `n` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
