\name{netEmbedding}
\alias{netEmbedding}
\title{The pinned upstream `netEmbedding`, verbatim}
\description{
Drop-in delegation: resolves the name for `library(cellchatrs)` users by calling the pinned
upstream body (see the header of this file for what delegation does and does not claim).
Upstream file: `R/analysis.R`. Needs Seurat/IRLBA at call time, exactly as upstream does.
}
\usage{
netEmbedding(...)
}
\arguments{
\item{...}{Arguments passed through unchanged. See the corresponding function in the \pkg{CellChat} package.}
}
\keyword{internal}
