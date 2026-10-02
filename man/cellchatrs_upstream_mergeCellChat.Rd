\name{cellchatrs_upstream_mergeCellChat}
\alias{cellchatrs_upstream_mergeCellChat}
\title{The pinned upstream `mergeCellChat`, verbatim}
\description{
`mergeCellChat` is **not** reimplemented in Rust, and the gate below says so at a named rung
rather than leaving it untested. It is pure S4 slot assembly -- `rbind` of `meta`, `cbind` of
`data`, a `union` of factor levels, and two `stop()` calls -- with no arithmetic anywhere, so
there is nothing to move to `r-core` and the 14.1 scope decision covers it.

What is worth testing is what it *produces*: a merged object whose `data`, `idents` and `meta`
become the input to `computeCommunProb` in the comparison-analysis workflow, so
`tests/parity/check_merge.R` runs the Rust kernel on a merged object and requires
`identical()` to upstream on that same object. That is a real test, and unlike comparing this
function with itself it is not vacuous.
}
\usage{
cellchatrs_upstream_mergeCellChat(...)
}
\arguments{
\item{...}{Arguments passed through unchanged. See the corresponding function in the \pkg{CellChat} package.}
}
\keyword{internal}
