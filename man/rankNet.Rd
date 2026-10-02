\name{rankNet}
\alias{rankNet}
\title{Only `mode = "single"` is Rust-backed, and only for the part that decides *which value goes in}
\description{
which row*: the `thresh` cut, the `measure = "count"` binarisation, the group filters, the
per-pathway `apply(prob, 3, sum)`, the `-1/log` rescaling with its degenerate-entry
reassignment, and the `order()`. The ggplot stays here (`PLAN.md` 14.1) because the plot is
the deliverable in R and the numerics are what `identical()` can check.

`mode = "comparison"`, `signaling.type`, and the `signaling` / `pairLR` filters still route to
upstream. The comparison mode needs the Wilcoxon kernel for `do.stat`, which is ported but not
yet wired in here; routing it out is honest and the gap is in `NOT_PORTED`.
}
\usage{
rankNet(object, slot.name = "netP", measure = c("weight", "count"), mode = c("comparison", "single"), comparison = c(1, 2), color.use = NULL, stacked = FALSE, sources.use = NULL, targets.use = NULL, signaling = NULL, pairLR = NULL, signaling.type = NULL, do.stat = FALSE, paired.test = TRUE, cutoff.pvalue = 0.05, tol = 0.05, thresh = 0.05, show.raw = FALSE, return.data = FALSE, x.rotation = 90, title = NULL, bar.w = 0.75, font.size = 8, do.flip = TRUE, x.angle = NULL, y.angle = 0, x.hjust = 1, y.hjust = 1, axis.gap = FALSE, ylim = NULL, segments = NULL, tick_width = NULL, rel_heights = c(0.9, 0, 0.1))
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{slot.name}{Argument `slot.name` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{measure}{Argument `measure` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{mode}{Argument `mode` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{comparison}{Argument `comparison` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{color.use}{Argument `color.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{stacked}{Argument `stacked` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{sources.use}{Argument `sources.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{targets.use}{Argument `targets.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{signaling}{Argument `signaling` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{pairLR}{Argument `pairLR` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{signaling.type}{Argument `signaling.type` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{do.stat}{Argument `do.stat` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{paired.test}{Argument `paired.test` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{cutoff.pvalue}{Argument `cutoff.pvalue` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{tol}{Argument `tol` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{thresh}{Argument `thresh` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{show.raw}{Argument `show.raw` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{return.data}{Argument `return.data` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{x.rotation}{Argument `x.rotation` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{title}{Argument `title` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{bar.w}{Argument `bar.w` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{font.size}{Argument `font.size` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{do.flip}{Argument `do.flip` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{x.angle}{Argument `x.angle` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{y.angle}{Argument `y.angle` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{x.hjust}{Argument `x.hjust` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{y.hjust}{Argument `y.hjust` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{axis.gap}{Argument `axis.gap` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{ylim}{Argument `ylim` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{segments}{Argument `segments` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{tick_width}{Argument `tick_width` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{rel_heights}{Argument `rel_heights` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
