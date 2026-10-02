\name{buildSNN}
\alias{buildSNN}
\title{The pinned upstream `buildSNN`, verbatim}
\description{
Drop-in delegation: resolves the name for `library(cellchatrs)` users by calling the pinned
upstream body (see the header of this file for what delegation does and does not claim).
Upstream file: `R/analysis.R`. Needs Seurat at call time, exactly as upstream does.
}
\usage{
buildSNN(...)
}
\arguments{
\item{...}{Arguments passed through unchanged. See the corresponding function in the \pkg{CellChat} package.}
}
\keyword{internal}
