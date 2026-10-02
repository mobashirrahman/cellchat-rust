# The single function upstream CellChat uses from BiocGenerics, and the only one needed.
#
# `subsetCommunication_internal` ends with
# `net <- BiocGenerics::as.data.frame(net, stringsAsFactors = FALSE)`. By that point `net`
# is already a data.frame, and BiocGenerics' method is base R's `as.data.frame.data.frame`,
# i.e. a no-op apart from dropping rownames. Returning the input unchanged is therefore
# bit-equivalent, and it lets the parity oracle run without a Bioconductor toolchain.
as.data.frame <- function(x, ...) {
  if (is.data.frame(x)) return(x)
  base::as.data.frame(x, ...)
}
