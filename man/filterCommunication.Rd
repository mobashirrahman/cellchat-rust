\name{filterCommunication}
\alias{filterCommunication}
\title{`filterCommunication`, backed by the Rust kernel}
\description{
The numerics -- the `min.cells` zeroing, the per-sample group means, the outer products,
the binarisation and the cross-sample consistency mask -- are in `r-core`. The `cat()`
output stays here because it is `scales::percent` and `toString(levels(...))`, i.e. R
formatting, and because *where* the messages appear is part of the observable behaviour.

## Two upstream behaviours that are reproduced rather than repaired

* **`nonFilter.keep = TRUE` is a no-op.** `net <- object@net` is bound *before* the two
`object@net$prob.nonFilter <- ...` assignments, and the function ends with
`object@net <- net`, which throws them away. The `cat()` still fires, so the flag has an
observable effect and no other one. See `crates/r-core/tests/filter_parity.rs`.
* **An all-zero `net$prob` with `min.samples >= 2` is `"subscript out of bounds"`**, from
`for (jj in 1:length(LR.nonzero))` with `LR.nonzero` empty being `1:0` and then
`score.LR[, , 0, i] <- ...`. The kernel returns that as an error carrying upstream's
text, so it is re-raised unchanged.

The `raw.use` setting is read back out of `options$parameter`, not taken as an argument, so a
`filterCommunication` on a projected object filters the *projected* data. The corpus exercises
both settings, because the two take visibly different paths through the kernel: the per-sample
means come from a different matrix, so the zero-fill and the consistency mask can differ even
when the network structure is identical.
}
\usage{
filterCommunication(object, min.cells = 10, min.samples = NULL, rare.keep = FALSE, nonFilter.keep = FALSE)
}
\arguments{
\item{object}{Argument `object` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{min.cells}{Argument `min.cells` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{min.samples}{Argument `min.samples` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{rare.keep}{Argument `rare.keep` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{nonFilter.keep}{Argument `nonFilter.keep` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
