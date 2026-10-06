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
  ## The five datasets are byte-identical to the pinned upstream copies, and that is what this
  ## asserts. It used to go through `utils::data(list=nm, package="CellChat")`, which required
  ## shipping a second copy of all five in `data/` -- 2.9 MB, and with `inst/upstream` also
  ## present that is what puts the tarball over CRAN's 5 MB limit. The copy inside
  ## `inst/upstream/<sha>/data/` still ships, so load it from there directly and keep the
  ## identity assertion, which is the part with teeth. What is given up is only the LazyData
  ## routing, i.e. `data(CellChatDB.human)` after a tarball install.
  for (nm in c("CellChatDB.human","CellChatDB.mouse","CellChatDB.zebrafish","PPI.human","PPI.mouse")) {
    e <- new.env(parent=emptyenv())
    load(file.path(system.file("upstream/CellChat-75253cd0/data",package="CellChat"),paste0(nm,".rda")),envir=e)
    stopifnot(exists(nm,e,inherits=FALSE),identical(serialize(e[[nm]],NULL),serialize(reference[[nm]],NULL)))
    ## Upstream's LazyData makes the bare name resolve after `library(CellChat)`. `.onLoad` binds
    ## the same names to the bundled copy; assert the binding is exported, visible on the search
    ## path, and the same bytes -- without `data()`, which a built install cannot satisfy.
    stopifnot(nm %in% getNamespaceExports("CellChat"),
              identical(serialize(getExportedValue("CellChat",nm),NULL),serialize(e[[nm]],NULL)),
              identical(serialize(get(nm,envir=as.environment("package:CellChat"),inherits=FALSE),NULL),serialize(e[[nm]],NULL)))
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
