\name{subsetCommunication}
\alias{subsetCommunication}
\title{Two paths, split at the melt exactly as upstream splits it:}
\description{
* `net` is a 3-d array, no DEG threshold and `slot.name = "net"`: the melt, the `thresh`
cut, the L-R join and the final column selection.
* otherwise -- `net` is already a data frame, or any of the eight DEG thresholds is given,
or `slot.name = "netP"` -- the same melt, then the thresholds, the `netP` aggregation and
the final selection. The melt runs *unselected* in this path, because upstream's
`rowSums(is.na(net)) != ncol(net)` counts over the full column set.

`mode = "merged"` and `pairLR.use` / `signaling` still route to upstream: `signaling` is a
database query (`searchPair`) rather than a numeric kernel, and `pairLR.use` needs the
tryCatch-over-two-columns fallback.
}
\usage{
subsetCommunication(object = NULL, net = NULL, slot.name = "net", sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR.use = NULL, thresh = 0.05, datasets = NULL, ligand.pvalues = NULL, ligand.logFC = NULL, ligand.pct.1 = NULL, ligand.pct.2 = NULL, receptor.pvalues = NULL, receptor.logFC = NULL, receptor.pct.1 = NULL, receptor.pct.2 = NULL)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{net}{Argument `net` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{slot.name}{Argument `slot.name` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{sources.use}{Argument `sources.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{targets.use}{Argument `targets.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{signaling}{Argument `signaling` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{pairLR.use}{Argument `pairLR.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{thresh}{Argument `thresh` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{datasets}{Argument `datasets` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{ligand.pvalues}{Argument `ligand.pvalues` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{ligand.logFC}{Argument `ligand.logFC` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{ligand.pct.1}{Argument `ligand.pct.1` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{ligand.pct.2}{Argument `ligand.pct.2` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{receptor.pvalues}{Argument `receptor.pvalues` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{receptor.logFC}{Argument `receptor.logFC` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{receptor.pct.1}{Argument `receptor.pct.1` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{receptor.pct.2}{Argument `receptor.pct.2` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
