suppressPackageStartupMessages(library(CellChat))
local({
  reference <- get("cellchatrs_upstream_cached",asNamespace("CellChat"))()
  exports <- readLines(system.file("upstream/CellChat-75253cd0/NAMESPACE",package="CellChat"))
  exports <- sub("^export\\((.*)\\)$","\\1",grep("^export\\(",exports,value=TRUE))
  stopifnot(length(exports)==109L,all(exports %in% getNamespaceExports("CellChat")))
  for (nm in exports) {
    stopifnot(is.function(getExportedValue("CellChat",nm)),
              identical(formals(getExportedValue("CellChat",nm)),formals(reference[[nm]])))
  }
  for (nm in c("CellChatDB.human","CellChatDB.mouse","CellChatDB.zebrafish","PPI.human","PPI.mouse")) {
    e <- new.env(parent=emptyenv())
    utils::data(list=nm,package="CellChat",envir=e)
    stopifnot(exists(nm,e,inherits=FALSE),identical(serialize(e[[nm]],NULL),serialize(reference[[nm]],NULL)))
  }
  stopifnot(identical(attr(class(methods::new("CellChat")),"package"),"CellChat"))
  x <- matrix(seq_len(48)/48,4,dimnames=list(c("A","B","C","D"),paste0("c",1:12)))
  meta <- data.frame(labels=rep(c("G1","G2"),each=6),row.names=colnames(x))
  a <- suppressMessages(reference$createCellChat(x,meta=meta,group.by="labels"))
  b <- suppressMessages(createCellChat(x,meta=meta,group.by="labels"))
  stopifnot(methods::validObject(a),methods::validObject(b),identical(serialize(a,NULL),serialize(b,NULL)))
  stopifnot(identical(CellChat::scPalette(8),reference$scPalette(8)))
  # Construct a plot through a formerly absent public function without opening a device.
  a <- reference$compareInteractions(methods::new("CellChat",net=list(Control=list(count=diag(2)),Treatment=list(count=diag(2))),
             options=list(mode="merged"),meta=data.frame()))
  b <- compareInteractions(methods::new("CellChat",net=list(Control=list(count=diag(2)),Treatment=list(count=diag(2))),
             options=list(mode="merged"),meta=data.frame()))
  stopifnot(inherits(b,"ggplot"),identical(a$data,b$data))
})
