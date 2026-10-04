suppressPackageStartupMessages(library(CellChat))
local({
  reference <- get("cellchatrs_upstream_cached", asNamespace("CellChat"))()
  bytes <- function(x) serialize(x, NULL, version = 3L)
  quiet <- function(f) {v <- NULL; invisible(capture.output(v <- suppressMessages(f()))); v}
  x <- matrix(seq_len(48)/48, 4, dimnames=list(c("G1","G2","G3","G4"),paste0("c",1:12)))
  obj <- methods::new("CellChat", data=x, data.signaling=x,
      idents=factor(rep(c("A","B"),each=6)), options=list(mode="single",datatype="RNA"))
  supplied <- x[1:2,,drop=FALSE]
  for (fast in c(FALSE, TRUE)) {
    a <- quiet(function() reference$identifyOverExpressedGenes(obj, data.use=supplied,
          do.DE=FALSE,do.fast=fast,min.cells=1))
    b <- quiet(function() identifyOverExpressedGenes(obj, data.use=supplied,
          do.DE=FALSE,do.fast=fast,min.cells=1))
    stopifnot(identical(bytes(a), bytes(b)))
  }
  a <- quiet(function() reference$identifyOverExpressedGenes(obj,supplied,do.DE=FALSE,do.fast=FALSE,min.cells=1))
  b <- quiet(function() identifyOverExpressedGenes(obj,supplied,do.DE=FALSE,do.fast=FALSE,min.cells=1))
  stopifnot(identical(bytes(a),bytes(b)),
      identical(formals(identifyOverExpressedGenes),formals(reference$identifyOverExpressedGenes)))

  db <- reference$CellChatDB.human
  lr <- db$interaction["TGFB1_TGFBR1_TGFBR2",,drop=FALSE]
  obj@DB <- db
  obj@LR <- list(LRsig=lr)
  rng_call <- function(f, object) {
    set.seed(314159L)
    error <- tryCatch({quiet(function() f(object,nboot=2));NULL},error=conditionMessage)
    list(error=error,seed=.Random.seed,next.draw=runif(1))
  }
  # A missing complex subunit fails before set.seed and bootstrapping.
  stopifnot(identical(bytes(rng_call(reference$computeCommunProb,obj)),
                     bytes(rng_call(computeCommunProb,obj))))
  # A NaN Hill probability fails after bootstrapping, consuming the same draws.
  obj@LR$LRsig$ligand <- "G1"
  obj@LR$LRsig$receptor <- "G2"
  obj@LR$LRsig$co_A_receptor <- obj@LR$LRsig$co_I_receptor <- ""
  obj@data.signaling[,] <- 0
  stopifnot(identical(bytes(rng_call(reference$computeCommunProb,obj)),
                     bytes(rng_call(computeCommunProb,obj))))
})
