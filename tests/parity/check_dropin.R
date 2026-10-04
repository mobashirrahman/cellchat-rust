#!/usr/bin/env Rscript
# Compare independently installed packages in separate clean processes.
root<-Sys.getenv("CELLCHATRS_ROOT",getwd())
original<-Sys.getenv("CELLCHAT_ORACLE_LIB")
if(!nzchar(original)||!file.exists(file.path(original,"CellChat","DESCRIPTION")))
 stop("CELLCHAT_ORACLE_LIB must contain the independently installed pinned original package")
replacement<-Sys.getenv("CELLCHAT_REPLACEMENT_LIB",file.path(root,".rlib"))
out<-Sys.getenv("CELLCHAT_CONTRACT_OUT",tempfile("cellchat-contract-"))
dir.create(out,recursive=TRUE,showWarnings=FALSE)
run_side<-function(lib,side,import=NULL){
 dest<-file.path(out,side);dir.create(dest,showWarnings=FALSE)
 env<-c(paste0("R_LIBS=",paste(c(lib,Sys.getenv("R_LIBS")),collapse=.Platform$path.sep)))
 args<-c("--vanilla",shQuote(file.path(root,"tests/parity/package_worker.R")),shQuote(dest))
 if(!is.null(import))args<-c(args,shQuote(import))
 status<-system2(file.path(R.home("bin"),"Rscript"),args,env=env,
                 stdout=file.path(out,paste0(side,".log")),stderr=file.path(out,paste0(side,".log")))
 if(status!=0L)stop(side," worker failed: ",paste(tail(readLines(file.path(out,paste0(side,".log"))),12),collapse="\n"))
 dest
}
a<-run_side(original,"original")
b<-run_side(replacement,"replacement",file.path(a,"original-object.rds"))
api_a<-unserialize(readBin(file.path(a,"api.bin"),"raw",n=file.info(file.path(a,"api.bin"))$size))
api_b<-unserialize(readBin(file.path(b,"api.bin"),"raw",n=file.info(file.path(b,"api.bin"))$size))
stopifnot(length(api_a)==109L,all(names(api_a)%in%names(api_b)),identical(api_a,api_b[names(api_a)]))
files<-setdiff(list.files(a,pattern="\\.(bin|png)$"),"api.bin")
for(name in files){
 read<-function(path)readBin(path,"raw",n=file.info(path)$size)
 stopifnot(identical(read(file.path(a,name)),read(file.path(b,name))))
}
stopifnot(identical(readBin(file.path(a,"constructor.bin"),"raw",n=1e7),
                    readBin(file.path(b,"imported-object.bin"),"raw",n=1e7)))
cat("Independent package contract:",length(files)+2L,"byte comparisons passed;",length(api_a),"public signatures matched\n")
cat("Evidence:",out,"\n")
