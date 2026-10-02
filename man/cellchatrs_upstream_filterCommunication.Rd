\name{cellchatrs_upstream_filterCommunication}
\alias{cellchatrs_upstream_filterCommunication}
\title{The pinned upstream `filterCommunication`, verbatim}
\description{
The reference for the differential test and the fallback under
`CELLCHATRS_FALLBACK=1`. Wrapped rather than reached for directly because
`cellchatrs_upstream_cached` is internal, and `get("filterCommunication", envir = env)`
from a script would otherwise find *this package's* shim -- the two functions share a
name, and the whole point of the comparison is that they differ.
}
\usage{
cellchatrs_upstream_filterCommunication(object, ...)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{...}{Arguments passed through unchanged. See the corresponding function in the \pkg{CellChat} package.}
}
\keyword{internal}
