\name{computeCellDistance}
\alias{computeCellDistance}
\title{`computeCellDistance` from pinned upstream, with `collapse::fdist` replaced by the Rust kernel}
\description{
Returns a `dist` object, as upstream does -- `collapse::fdist` does not return a matrix, and the
rewrap reproduces the three places `stats::as.dist` disagrees with it (no `call` attribute,
`method = "euclidean"` present, `Labels` `NULL` rather than synthesised). Also reproduces the
two-column check, `fdist`'s own "at least 2 rows" error, `ratio` being applied whenever it is
non-`NULL`, and the threshold that needs **both** `interaction.range` and `tol` with `tol` as an
addend. See `docs/SEMANTICS.md`.
}
\usage{
computeCellDistance(coordinates, interaction.range = NULL, ratio = NULL, tol = NULL)
}
\arguments{
\item{coordinates}{Argument `coordinates` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{interaction.range}{Argument `interaction.range` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{ratio}{Argument `ratio` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{tol}{Argument `tol` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
