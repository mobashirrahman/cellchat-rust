\name{computeAveExpr}
\alias{computeAveExpr}
\title{`computeAveExpr`, backed by the Rust kernel}
\description{
Handles `type` in `c("triMean", "truncatedMean", "median")` with R's `match.arg`
partial matching, and `features` in R's `intersect` order.
}
\usage{
computeAveExpr(object, features = NULL, group.by = NULL, type = c("triMean", "truncatedMean", "median"), trim = NULL, slot.name = c("data.signaling", "data"), data.use = NULL)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{features}{Argument `features` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{group.by}{Argument `group.by` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{type}{Argument `type` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{trim}{Argument `trim` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{slot.name}{Argument `slot.name` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{data.use}{Argument `data.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
