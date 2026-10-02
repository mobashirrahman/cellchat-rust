.libPaths(".rlib")
suppressMessages({library(collapse); library(Matrix); library(dplyr)})
CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
methods::setClass("CellChatProbe", representation(data="ANY", data.signaling="ANY", data.smooth="ANY",
  idents="ANY", meta="ANY", images="ANY", DB="ANY", LR="ANY", net="ANY", netP="ANY", options="ANY"))
env <- new.env(); source(file.path(CC,"R","modeling.R"), local=env)
e <- new.env(); load(file.path(CC,"data","CellChatDB.human.rda"), envir=e); DB <- get(ls(e)[1], e)
NLR <- as.integer(Sys.getenv("NLR","2239")); NB <- as.integer(Sys.getenv("NB","10"))
nC <- as.integer(Sys.getenv("NC","5000")); K <- as.integer(Sys.getenv("K","25"))
ii <- DB$interaction
cx <- as.matrix(DB$complex[,1:5]); cs <- unique(as.character(cx)[cx!=""])
cof <- unique(c(ii$agonist, ii$antagonist, ii$co_A_receptor, ii$co_I_receptor))
cof <- cof[!is.na(cof) & cof != ""]
cm <- as.matrix(DB$cofactor[match(cof, rownames(DB$cofactor), nomatch=0), 1:16, drop=FALSE])
mv <- unique(as.character(cm)[cm != ""])
genes <- unique(c(stats::na.omit(unlist(DB$complex)), stats::na.omit(unlist(DB$cofactor)),
                ii$ligand, ii$receptor, cs, mv))
cat("n signaling genes:", length(genes), "\n")
set.seed(2); M <- rsparsematrix(length(genes), nC, 0.07, rand.x=function(z) rpois(z,3))
rownames(M) <- genes; colnames(M) <- paste0("c", 1:nC)
LR <- ii[ii$annotation %in% c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact") &
          ii$ligand %in% genes & ii$receptor %in% genes, ]
LR <- LR[seq_len(min(NLR, nrow(LR))), ]
cat("nLR:", nrow(LR), "\n"); gc(reset=TRUE, full=FALSE)
obj <- methods::new("CellChatProbe", data=M, data.signaling=M, data.smooth=NULL,
  idents=factor(sample(paste0("cl",1:K), nC, TRUE)),
  meta=data.frame(samples=factor(rep("s1", nC))), images=list(),
  DB=DB, LR=list(LRsig=LR), net=list(), netP=list(), options=list(datatype="RNA"))
cat("setup done; maxmem:", round(sum(gc()[,2]),1), "MB\n")

tA <- Sys.time(); X <- as.matrix(obj@data.signaling)
cat("as.matrix          :", round(as.numeric(difftime(Sys.time(),tA,units="secs")),3), "s\n")
tB <- Sys.time(); Xn <- X/max(X)
cat("max+divide         :", round(as.numeric(difftime(Sys.time(),tB,units="secs")),3), "s\n")
g <- obj@idents
tC <- Sys.time(); avg <- aggregate(t(Xn), list(g), env$triMean)
tAgg <- as.numeric(difftime(Sys.time(), tC, units="secs"))
cat("aggregate x1       :", round(tAgg,3), "s  -> x nboot=", NB, ":", round(tAgg*NB/60,2), "min\n", sep="")
cat("maxmem             :", round(sum(gc()[,2]),1), "MB\n")