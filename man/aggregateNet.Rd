\name{aggregateNet}
\alias{aggregateNet}
\title{`aggregateNet`, backed by the Rust kernel}
\description{
Both branches are ported. The unfiltered one aggregates the probability array directly.

The filtered one (`sources.use` / `targets.use` / `signaling` / `pairLR.use`) is ported
because it turned out to be *broken upstream* at this commit, and reproducing a bug is
only possible if you know it is one. `stringr::str_split(key, "|", simplify = TRUE)`
takes a regular expression, `|` is the alternation metacharacter, so `a[, 1]` is `""` and
`a[, 2]` is the key's first character; the level restriction then finds no cell group and
`tapply` yields a `0 x 0` matrix under `remove.isolate = TRUE` and a `k x k` matrix of
zeros under `remove.isolate = FALSE`. `net$count` and `net$weight` therefore carry no
information for **any** filter. See `docs/SEMANTICS.md` and
`r_core::net::aggregate_net_filtered`.

The `message()` about isolate cell groups is emitted for `remove.isolate = TRUE` in the
filtered branch only, as upstream does.
}
\usage{
aggregateNet(object, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05, return.object = TRUE)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{sources.use}{Argument `sources.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{targets.use}{Argument `targets.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{signaling}{Argument `signaling` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{pairLR.use}{Argument `pairLR.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{remove.isolate}{Argument `remove.isolate` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{thresh}{Argument `thresh` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{return.object}{Argument `return.object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
