\name{rankNetPairwise}
\alias{rankNetPairwise}
\title{`rankNetPairwise`, backed by the Rust kernel's ordering}
\description{
Upstream's body is one `order` and a great deal of R plumbing:

```r
for (i in 1:numCluster) for (j in 1:numCluster) \{
data <- data.frame(pathway_index = index, interaction_name = ..., ..., row.names = rownames(pairLR.use))
temp[[j]] <- data[with(data, order(pval, -prob)), ]
\}
object@net$pairwiseRank <- pairwiseLR
```

The **ordering** is the only arithmetic, and it is done in Rust: `order(pval, -prob)` with two
numeric keys is a stable radix sort, so the permutation is a total function of
`(pval, -prob, index)` and `order_f64_multi` is a faithful transcription of it. `na.last` is R's
default `TRUE`, which is why a missing `pval` sorts last rather than first.

Everything else stays in R, for the same reason `format` and `rankNet`'s row assembly do: the
`data.frame()` column order, `row.names = rownames(pairLR.use)`, the nested `list()` and the
`names(temp) <- colnames(prob)` are language behaviour rather than numerics, and reproducing
them in Rust would buy nothing but new ways to be wrong.

`pairLR.use` supplies only the data-frame's **columns**; the ordering is over the whole third
dimension of `prob`/`pval`. Upstream does not reconcile the two, and that is load-bearing: with
`LR.use` shorter than `dim(prob)[3]`, its `data.frame(..., row.names = rownames(pairLR.use))`
raises `row names supplied are of the wrong length`. Slicing `prob` down to `LR.use` first --
which reads like a faithful interpretation of "use this subset of interactions" -- silently
turns that error into a successful shorter result.

Upstream's loops are `for (i in 1:numCluster) for (j in 1:numCluster)` with
`numCluster <- dim(prob)[1]`, so it reads the leading `k x k` block and never touches the rest.
`rownames(prob)` and `colnames(prob)` are still used for the names, which is where a non-square
`prob` raises.
}
\usage{
rankNetPairwise(object, LR.use = NULL)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{LR.use}{Argument `LR.use` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
