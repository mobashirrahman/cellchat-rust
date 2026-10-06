\name{cellchatrs_computeRegionDistance_exact}
\alias{cellchatrs_computeRegionDistance_exact}
\title{Exact spatial distances, an explicit alternative to CellChat compatibility}
\description{
Exact spatial distances, an explicit alternative to CellChat compatibility.
}
\usage{
    cellchatrs_computeRegionDistance_exact(coordinates, meta, interaction.range = NULL,
    ratio = NULL, tol = NULL, k.min = 10, contact.dependent = TRUE, contact.range = NULL,
    contact.knn.k = NULL, do.symmetric = TRUE)
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
