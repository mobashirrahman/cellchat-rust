#!/usr/bin/env Rscript
# Runs in a clean process with either the independently installed original or this package.
args <- commandArgs(trailingOnly=TRUE)
suppressPackageStartupMessages(library(CellChat))
out <- args[[1]]
dir.create(out, recursive=TRUE,showWarnings=FALSE)
raw <- function(x) serialize(x,NULL,version=3L)
write_raw <- function(name,x) writeBin(raw(x),file.path(out,paste0(name,".bin")))
quiet <- function(f) {v<-NULL;invisible(capture.output(v<-suppressMessages(f())));v}
canonical <- function(x) {if(methods::is(x,"CellChat"))x@options$run.time<-NULL;x}
exports <- sort(getNamespaceExports("CellChat"))
exports <- exports[vapply(exports,function(nm)is.function(getExportedValue("CellChat",nm)),logical(1))]
write_raw("api",setNames(lapply(exports,function(nm) if(is.function(getExportedValue("CellChat",nm)))
  formals(getExportedValue("CellChat",nm)) else NULL),exports))
for(nm in c("CellChatDB.human","CellChatDB.mouse","CellChatDB.zebrafish","PPI.human","PPI.mouse")){
 e<-new.env();utils::data(list=nm,package="CellChat",envir=e);write_raw(nm,e[[nm]])
}
set.seed(101)
x<-matrix(runif(72,.1,1),6,dimnames=list(c("G1","G2","G3","INHBA","INHBB","TGFBR1"),paste0("c",1:12)))
meta<-data.frame(labels=rep(c("A","B"),each=6),samples=factor(rep("s1",12)),row.names=colnames(x))
obj<-quiet(function()createCellChat(x,meta=meta,group.by="labels"))
stopifnot(methods::validObject(obj),identical(attr(class(obj),"package"),"CellChat"))
write_raw("constructor",obj)
saveRDS(obj,file.path(out,"original-object.rds"),version=3L)
if(length(args)>1){ imported<-readRDS(args[[2]]);stopifnot(methods::validObject(imported));write_raw("imported-object",imported) }
db<-CellChatDB.human
db$complex["Activin AB","subunit_1"]<-"G3"
obj@data.signaling<-x
obj@DB<-db
lr<-db$interaction[1,,drop=FALSE]
lr$ligand<-"Activin AB";lr$receptor<-"TGFBR1"
lr$agonist<-lr$antagonist<-lr$co_A_receptor<-lr$co_I_receptor<-""
obj@LR<-list(LRsig=lr)
computed<-quiet(function()computeCommunProb(obj,nboot=2))
write_raw("probability",canonical(computed));write_raw("probability-rng",.Random.seed)
write_raw("oeg-supplied",quiet(function()identifyOverExpressedGenes(obj,data.use=x[1:2,,drop=FALSE],do.DE=FALSE,do.fast=FALSE,min.cells=1)))
set.seed(41)
co<-matrix(runif(96),48,2)
spmeta<-data.frame(group=factor(rep(letters[1:4],each=12)),samples=factor(rep("s1",48)))
write_raw("spatial",computeRegionDistance(co,spmeta,ratio=1,tol=.1,k.min=1,
  interaction.range=10,contact.range=2,contact.knn.k=3))
set.seed(25)
snnx<-matrix(rnorm(90),6,dimnames=list(paste0("g",1:6),paste0("c",1:15)))
for(k in c(2,4,8))write_raw(paste0("snn-",k),buildSNN(snnx,k=k))
# The original native helper is an independent oracle for CSC structure, duplicates and pruning.
helper<-get("ComputeSNN",asNamespace("CellChat"))
nn<-matrix(c(1L,2L,3L,4L,2L,1L,4L,3L),4,2)
for(prune in c(-1,0,1/15,.4,1,2))write_raw(paste0("snn-prune-",prune),helper(nn,prune))
duplicates<-matrix(rep(1:4,2),4,2)
write_raw("snn-duplicates",helper(duplicates,0))
# Exercise a restored public rendering entry point on a deterministic graph.
set.seed(11)
net<-matrix(c(0,.2,.3,.4,0,.5,.6,.7,0),3,dimnames=list(LETTERS[1:3],LETTERS[1:3]))
grDevices::png(file.path(out,"circle.png"),width=640,height=480)
quiet(function()netVisual_circle(net,vertex.weight=c(1,2,3),weight.scale=TRUE))
grDevices::dev.off()
write_raw("palette",scPalette(8))
write_raw("class-provenance",class(obj))
