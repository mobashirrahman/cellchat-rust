\name{subsetDB}
\alias{subsetDB}
\title{`subsetDB`, backed by the Rust kernel}
\description{
Only `key = "annotation"` is routed here; any other key falls back, because the message
for an unknown key is part of the contract.
}
\usage{
subsetDB(CellChatDB, search = c(), key = "annotation", non_protein = FALSE)
}
\arguments{
\item{CellChatDB}{Argument `CellChatDB` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{search}{Argument `search` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{key}{Argument `key` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
\item{non_protein}{Argument `non_protein` of the corresponding \pkg{CellChat} function; see its documentation and \code{docs/SEMANTICS.md} in this package.}
}
\keyword{internal}
