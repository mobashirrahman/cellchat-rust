## The upstream S4 class is registered when cellchatrs loads, so exported constructors and
## deserialized objects work without lazily sourcing CellChat into .GlobalEnv.
methods::setClassUnion("AnyMatrix", c("matrix", "dgCMatrix"))
methods::setClassUnion("AnyFactor", c("factor", "list"))

CellChat <- methods::setClass(
  "CellChat",
  slots = c(
    data.raw = "AnyMatrix",
    data = "AnyMatrix",
    data.signaling = "AnyMatrix",
    data.scale = "matrix",
    data.smooth = "AnyMatrix",
    images = "list",
    net = "list",
    netP = "list",
    meta = "data.frame",
    idents = "AnyFactor",
    DB = "list",
    LR = "list",
    var.features = "list",
    dr = "list",
    options = "list"
  )
)

methods::setMethod(
  "show",
  "CellChat",
  function(object) {
    if (object@options$mode == "single") {
      cat("An object of class", class(object), "created from a single dataset", "\n",
          nrow(object@data), "genes.\n", ncol(object@data), "cells. \n")
    } else if (object@options$mode == "merged") {
      cat("An object of class", class(object),
          "created from a merged object with multiple datasets", "\n",
          nrow(object@data.signaling), "signaling genes.\n",
          ncol(object@data.signaling), "cells. \n")
    }
    if (object@options$datatype == "RNA") {
      cat("CellChat analysis of single cell RNA-seq data! \n")
    } else {
      cat("CellChat analysis of", object@options$datatype,
          "data! The input spatial locations are \n")
      print(utils::head(object@images$coordinates))
    }
    invisible(NULL)
  }
)
