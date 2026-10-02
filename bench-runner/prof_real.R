.libPaths("/scratch/mdra00001/cellchat-rust/.rlib")
suppressMessages({library(collapse); library(Matrix); library(dplyr)})
CC <- "/scratch/mdra00001/tmp/opencode/CellChat"
methods::setClass("CellChatProbe", representation(data="ANY", data.signaling="ANY", data.smooth="ANY",
  idents="ANY", meta="ANY", images="ANY", DB="ANY", LR="ANY", net="ANY", netP="ANY", options="ANY"))
env <- new.env(); source(file.path(CC,"R","modeling.R"), local=env)
e <- new.env(); load(file.path(CC, "data", paste0("CellChatDB.", Sys.getenv("SPECIES","human"), ".rda")), envir=e); DB <- get(ls(e)[1], e)

f <- Sys.getenv("FIX", "/scratch/mdra00001/tmp/opencode/data/humanSkin.rda")
NB <- as.integer(Sys.getenv("NB","10")); MAXLR <- as.integer(Sys.getenv("MAXLR","0"))
raw <- new.env(); load(f, envir=raw); o <- raw[[ls(raw)[1]]]
counts <- o$data
grp <- if (!is.null(o$meta)) factor(o$meta$labels) else factor(o$labels)
nC <- ncol(counts)
cat(sprintf("FIXTURE %s : %d genes x %d cells, K=%d, nnz=%.1f%%\n", basename(f), nrow(counts), nC,
            nlevels(grp), 100*Matrix::nnzero(counts)/length(counts)))

extractSig <- function(DB) {
  ii <- DB$interaction
  cs <- as.matrix(DB$complex[, grepl("^subunit", colnames(DB$complex)), drop=FALSE])
  cof <- unique(c(ii$agonist, ii$antagonist, ii$co_A_receptor, ii$co_I_receptor))
  cof <- cof[!is.na(cof) & cof != ""]
  cfm <- DB$cofactor[match(cof, rownames(DB$cofactor), nomatch=0), , drop=FALSE]
  cm <- as.matrix(cfm[, grepl("^cofactor", colnames(cfm)), drop=FALSE])
  unique(c(ii$ligand, ii$receptor,
           unique(as.character(cs)[cs != ""]), unique(as.character(cm)[cm != ""])))
}
genes.sig <- sort(intersect(extractSig(DB), rownames(counts)))
cat("data.signaling genes:", length(genes.sig), "\n")
# emulate CellChat's normalised input: log1p CPM-style scaling (sparse -> dense)
lib <- colSums(counts)
X <- as.matrix(counts[genes.sig, , drop=FALSE])
X <- X / rep(lib, each=nrow(X)) * 1e4
X <- log1p(X)
present <- c(rownames(counts), rownames(DB$complex))
subOK <- function(g) {
  if (g %in% rownames(counts)) return(TRUE)
  if (!(g %in% rownames(DB$complex))) return(FALSE)
  s <- as.character(DB$complex[g, ]); s <- s[!is.na(s) & s != ""]
  length(s) > 0 && all(s %in% rownames(counts))
}
LR <- DB$interaction[DB$interaction$annotation %in% c("Secreted Signaling","ECM-Receptor","Cell-Cell Contact"), ]
keep <- vapply(LR$ligand, subOK, TRUE) & vapply(LR$receptor, subOK, TRUE)
LR <- LR[keep, ]
cat("LR pairs with all subunits present:", nrow(LR), "\n")
if (MAXLR > 0) LR <- LR[seq_len(min(MAXLR, nrow(LR))), ]
cat("nLR:", nrow(LR), " density(data.signaling):", round(100*Matrix::nnzero(as.matrix(counts[genes.sig,]))/length(counts[genes.sig,]),1), "%\n")
M <- Matrix(X, sparse=FALSE); M <- as(M, "dgCMatrix")
obj <- methods::new("CellChatProbe", data=M, data.signaling=M, data.smooth=NULL, idents=grp,
  meta=data.frame(samples=factor(rep("s1", nC))), images=list(),
  DB=DB, LR=list(LRsig=LR), net=list(), netP=list(), options=list(datatype="RNA"))
gc(reset=TRUE, full=FALSE)
t0 <- Sys.time(); o2 <- env$computeCommunProb(obj, nboot=NB, seed.use=1)
el <- as.numeric(difftime(Sys.time(), t0, units="secs"))
cat(sprintf("RESULT nC=%d K=%d nGenes=%d nLR=%d nboot=%d total=%.2fs maxmem=%.0fMB\n",
            nC, nlevels(grp), length(genes.sig), nrow(LR), NB, el, sum(gc()[,2])))
pr <- o2@net$prob
cat("prob: dim=", paste(dim(pr), collapse="x"), " nonzero=", sum(pr>0), " frac=",
    round(mean(pr>0),4), " max=", signif(max(pr),5), "\n")
cat("pval: unique=", length(unique(o2@net$pval)), " min=", min(o2@net$pval), " frac<0.05=",
    round(mean(o2@net$pval<0.05),4), "\n")