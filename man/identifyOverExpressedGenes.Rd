\name{identifyOverExpressedGenes}
\alias{identifyOverExpressedGenes}
\title{`identifyOverExpressedGenes`, backed by the Rust kernel}
\description{
Routes to Rust only where the numerics are a *pure function of the input matrix*:
the `do.DE = FALSE` branch (percentage and `min.cells` only) and the
`do.DE = TRUE, do.fast = FALSE` branch (the Wilcoxon path), with or without
`group.dataset` / `group.DE.combined`. Everything else -- `do.fast = TRUE` (presto, an
optional C++ dependency that may not be installed, and a *different algorithm* rather than a
faster one) -- goes to the pinned upstream body.

`group.dataset` was long filed as blocked on presto and is not: upstream honours it in **both**
halves of `if (do.fast)`, at `utilities.R:429-484` for presto and `:512-520` for the Wilcoxon
test. Only the first needs presto. What it changes is the *cell selection* --

```r
cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
```

or, with `group.DE.combined = TRUE`, both pooled across groups with no `labels` term at all --
so the kernel takes the selection as a parameter instead of growing a flag through the path.

The argument resolution (`group.by` / `idents.use` / `invert` / `features.use`) and the
two `stop()` calls stay in R, so the error messages and the `message()` for a non-factor
`group.by` are upstream's. `features.use <- intersect(features, row.names(X))` is
R's `intersect`, which keeps the *first* argument's order -- so a caller-supplied
`features` vector reorders the output relative to the rownames, and that has to happen
in R, not in Rust.

## The schema collapses when nothing passes

Upstream builds `markers.all <- data.frame()` and only `rbind`s a group's rows once
they survive `pvalues < thresh.p`. If no group does, the frame stays zero-column, and
the next line, `markers.all$features <- as.character(markers.all$features)`, *adds* a
single zero-length column: the "marker table" is `0x1` with only `features`, not `0x7`.
The kernel therefore returns `n_before_only_pos` so the shim knows which of the two
shapes to build; see `r_core::wilcox::GroupMarkers`.
}
\usage{
identifyOverExpressedGenes(object, group.by = NULL, idents.use = NULL, invert = FALSE, group.dataset = NULL, pos.dataset = NULL, group.DE.combined = FALSE, features.name = "features", only.pos = TRUE, features = NULL, return.object = TRUE, thresh.pc = 0, thresh.fc = 0, thresh.p = 0.05, do.DE = TRUE, do.fast = TRUE, min.cells = 10)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{group.by}{Argument `group.by` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{idents.use}{Argument `idents.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{invert}{Argument `invert` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{group.dataset}{Argument `group.dataset` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{pos.dataset}{Argument `pos.dataset` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{group.DE.combined}{Argument `group.DE.combined` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{features.name}{Argument `features.name` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{only.pos}{Argument `only.pos` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{features}{Argument `features` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{return.object}{Argument `return.object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{thresh.pc}{Argument `thresh.pc` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{thresh.fc}{Argument `thresh.fc` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{thresh.p}{Argument `thresh.p` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{do.DE}{Argument `do.DE` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{do.fast}{Argument `do.fast` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{min.cells}{Argument `min.cells` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
