#!/usr/bin/env Rscript
## The object family for the configuration matrix: the gene universe, the L-R table, and the
## builders that turn one configuration into a `computeCommunProb`-shaped S4 object.
##
## Split out of `check_matrix.R` because a divergence report has to be reproducible from a single
## configuration, and re-deriving the object in a scratch probe is exactly where "it worked when I
## ran it by hand" comes from. A probe sources this file and calls `build_object()`; the matrix
## runner does the same, so both sides of a comparison are built by the same code.

## ------------------------------------------------------------------ the pinned upstream
##
## `computeRegionDistance` is replaced by the exact-neighbour oracle so the spatial half of the
## matrix can run at all; `BiocNeighbors` is not installable here and Annoy is approximate. The
## oracle's arithmetic is upstream's.
up <- get("cellchatrs_upstream_cached", envir = asNamespace("cellchatrs"))()
source(file.path(ROOT, "tests", "parity", "exact_neighbour_oracle.R"), local = TRUE)
assign("computeRegionDistance", compute_region_distance_exact, envir = up)

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)

setClass("MatrixChat", representation(
  data.signaling = "ANY", data.smooth = "ANY", LR = "list", LRsig = "data.frame",
  DB = "ANY", idents = "factor", meta = "data.frame", images = "list",
  options = "list", net = "ANY"))

## ------------------------------------------------------------------ the gene universe
##
## Deliberately the same L-R table and gene universe as `gen_prob_golden.R`, so the matrix overlaps
## the hand fixture and a regression shows up in both. The pre-flight `stopifnot` at the bottom of
## that file is reproduced here: a gene the kernel resolves but the matrix lacks makes *upstream*
## die with "subscript out of bounds", which reads as a port bug and is a fixture bug.
LR_TABLE <- data.frame(
  ligand       = c("G1", "Activin AB", "IL12AB", "G3", "G4", "INHBA", "IL23A", "G5"),
  receptor     = c("G2", "IL12B", "G6", "Activin AB", "G7", "FST", "G8", "IL12A"),
  agonist      = c("", "", "TGFb agonist", "", "", "TGFb antagonist", "", ""),
  antagonist   = c("", "ACTIVIN antagonist", "", "NODAL agonist", "", "", "", ""),
  co_A_receptor = c("", "", "", "TGFb agonist", "", "", "", ""),
  co_I_receptor = c("", "", "TGFb inhibition receptor", "", "", "", "", ""),
  annotation   = c("Secreted Signaling", "Secreted Signaling", "Secreted Signaling",
                   "ECM-Receptor", "ECM-Receptor", "Non-protein Signaling",
                   "Secreted Signaling", "Secreted Signaling"),
  stringsAsFactors = FALSE)
LR_TABLE$ligand_antagonist <- ifelse(LR_TABLE$antagonist == "", LR_TABLE$ligand,
                                     paste0(LR_TABLE$ligand, " - ", LR_TABLE$antagonist))
LR_TABLE$interaction_name <- paste(LR_TABLE$ligand, LR_TABLE$receptor, sep = "_")
rownames(LR_TABLE) <- LR_TABLE$interaction_name

subunits_of <- function(n, which) {
  tab <- if (which == "complex") DB$complex else DB$cofactor
  cols <- if (which == "complex") grep("^subunit", colnames(tab)) else grep("cofactor", colnames(tab))
  if (n %in% rownames(tab)) as.character(tab[n, cols]) else NA_character_
}
all_subunits <- unique(c(
  unlist(lapply(c(LR_TABLE$ligand, LR_TABLE$receptor), subunits_of, which = "complex")),
  unlist(lapply(c(LR_TABLE$agonist, LR_TABLE$antagonist,
                  LR_TABLE$co_A_receptor, LR_TABLE$co_I_receptor),
                 subunits_of, which = "cofactor"))))
all_subunits <- all_subunits[!is.na(all_subunits) & all_subunits != ""]
plain_lr <- setdiff(c(LR_TABLE$ligand, LR_TABLE$receptor), rownames(DB$complex))
GENES <- unique(c(paste0("G", 1:8), plain_lr, all_subunits))

## ------------------------------------------------------------------ object builders
##
## Each returns the object **and** the effective axes it actually ran with. `note` records any
## adaptation, so the coverage table is computed from reality rather than intent.

## Which L-R rows a structure keeps, as a 1-based index into `LR_TABLE`.
lr_rows <- function(struct, nlr) {
  plain <- which(LR_TABLE$ligand %in% plain_lr & LR_TABLE$receptor %in% plain_lr)
  cx <- which(apply(LR_TABLE, 1, function(r) r[["ligand"]] %in% rownames(DB$complex) ||
                                          r[["receptor"]] %in% rownames(DB$complex)))
  cof <- which(LR_TABLE$co_A_receptor != "" | LR_TABLE$co_I_receptor != "")
  ago <- which(LR_TABLE$agonist != "")
  ant <- which(LR_TABLE$antagonist != "")
  pick <- switch(struct,
    plain = plain,
    complex = cx,
    cofactor = cof,
    agonist = ago,
    antagonist = ant,
    mixed = sort(unique(c(plain, cx, cof, ago, ant))),
    ## A complex whose subunit is removed from the matrix -- upstream's "subscript out of bounds".
    missing_subunit = cx,
    plain)
  if (!length(pick)) pick <- plain
  ## Cycle if the structure has fewer rows than `nLR` asks for, so `nLR` stays a real axis.
  idx <- rep(pick, length.out = nlr)
  sort(unique(idx))[seq_len(min(nlr, length(unique(idx))))]
}

## `%||%`, so the builders default an absent axis rather than failing on a config built by hand.
`%||%` <- function(a, b) if (is.null(a)) b else a

## The expression matrix, with the pathological class applied.
build_data <- function(cfg) {
  nc <- as.integer(cfg$nC); k <- as.integer(cfg$K)
  levs <- paste0("g", seq_len(k))
  ## Balanced-ish groups, every one non-empty (upstream's own `nlevels` check).
  grp <- factor(rep(levs, length.out = nc), levels = levs)
  set.seed(1000L + nc * 31L + k)
  m <- matrix(stats::runif(length(GENES) * nc, 0.01, 1), nrow = length(GENES), ncol = nc,
              dimnames = list(GENES, paste0("c", seq_len(nc))))
  ## Structure that makes the branches reachable rather than drowned in noise, exactly as the
  ## hand fixture does: an all-zero gene, a single-group gene, a constant gene, and a zero pair.
  m["G1", ] <- 0
  if (k >= 2) m["G2", grp != "g2"] <- 0
  m["G3", ] <- 1
  m["G4", ] <- 0.5
  switch(as.character(cfg$data),
    all_zero = m[] <- 0,
    with_na = { m[1, 1] <- NA_real_; m[2, 1] <- NA },
    with_nan = { m[1, 1] <- NaN },
    with_inf = { m[1, 1] <- Inf; m[2, 1] <- -Inf },
    ## `dup_rownames` is applied in `build_object`, *after* the gene subsetting: R's
    ## `d$expr[genes, , drop = FALSE]` indexes by name, and a duplicate row name makes that raise
    ## "subscript out of bounds" inside the fixture builder -- an error naming neither the
    ## duplicate nor the builder.
    NULL)
  list(expr = m, group = grp, ncell = nc)
}

build_object <- function(cfg) {
  eff <- cfg
  notes <- character(0)
  d <- build_data(cfg)
  ## `single_cell` has to collapse the group axis too: a one-cell factor cannot have four levels
  ## without an empty one, and upstream's `nlevels != length(unique)` check rejects that outright.
  ## Recorded as an adaptation so the coverage table does not claim K's other levels ran here.
  if (as.integer(cfg$nC) == 1L && as.integer(cfg$K) > 1L) {
    eff$K <- 1L
    notes <- c(notes, "K forced to 1: a one-cell object cannot have more than one non-empty level")
    d <- build_data(eff)
  }
  rows <- lr_rows(as.character(cfg$lrstruct), as.integer(cfg$nLR))
  LR <- LR_TABLE[rows, , drop = FALSE]
  rownames(LR) <- LR$interaction_name

  genes <- GENES
  if (as.character(cfg$lrstruct) == "missing_subunit") {
    ## Drop one subunit of the first complex the kept rows name. Upstream's `computeExpr_complex`
    ## indexes the matrix with it and dies with "subscript out of bounds"; the port must match.
    need <- subunits_of(LR$ligand[LR$ligand %in% rownames(DB$complex)][1], "complex")
    need <- need[!is.na(need) & need %in% genes]
    if (length(need)) {
      genes <- setdiff(genes, need[1])
      notes <- c(notes, paste0("dropped subunit ", need[1], ": upstream's subscript-out-of-bounds"))
    }
  }
  expr <- d$expr[genes, , drop = FALSE]
  ## The `scale` axis: a constant factor on the matrix. Upstream divides by `max(data)` before
  ## aggregating, so this must not move `Prob` at all -- and making it move is exactly how the
  ## missing division was found. `eff$scale` is recorded so the coverage table counts what ran.
  sc <- as.numeric(cfg$scale %||% 1)
  if (is.finite(sc) && sc != 1 && sc != 0) expr <- expr * sc
  if (identical(as.numeric(cfg$scale), 0)) {
    ## `max == 0` makes `data/max(data)` all `NaN`, which upstream turns into
    ## "missing value where TRUE/FALSE needed". Recorded as its own level for that reason.
    expr[] <- 0
  }

  ## Duplicate row names, last: the object the kernel sees is the only thing that should carry them.
  ## Renaming row 2 to row 1's name can *remove* a subunit from the universe under its own name,
  ## which is how `dup_rownames` reaches upstream's "subscript out of bounds" rather than being a
  ## mere label collision.
  if (identical(as.character(cfg$data), "dup_rownames") && nrow(expr) > 1L) {
    rownames(expr)[2] <- rownames(expr)[1]
  }

  ## `raw.use = FALSE` reads `data.smooth`, which upstream creates by `normalizeData`. A smooth
  ## matrix that is *not* the raw one is the point of the axis, so it is perturbed rather than
  ## copied -- an identical copy would make the two axes indistinguishable.
  smooth <- expr * 0.9 + 0.01

  coord <- NULL
  if (identical(as.character(cfg$datatype), "spatial")) {
    ## Groups in well-separated blocks, so the exact-neighbour substitution in the oracle is sound
    ## (a cell's nearest neighbour in another group is decided by a wide margin).
    nc <- ncol(expr)
    coord <- matrix(0, nrow = nc, ncol = 2)
    for (i in seq_len(nc)) {
      g <- match(as.character(d$group[i]), levels(d$group))
      j <- (i - 1L) %% max(1L, nc %/% max(1L, length(levels(d$group))))
      coord[i, ] <- c((g - 1L) * 25 + (j %% 8L) * 0.3, (j %/% 8L) * 0.25)
    }
    rownames(coord) <- colnames(expr)
  }

  meta <- data.frame(samples = factor(rep("s1", ncol(expr)), levels = "s1"),
                     row.names = colnames(expr))
  o <- new("MatrixChat",
           data.signaling = as(expr, "dgCMatrix"), data.smooth = as(smooth, "dgCMatrix"),
           LR = list(LRsig = LR), LRsig = LR,
           DB = list(complex = DB$complex, cofactor = DB$cofactor),
           idents = d$group, meta = meta,
           images = if (is.null(coord)) list() else
             list(coordinates = coord,
                  spatial.factors = list(ratio = 1, tol = 0.5)),
           options = list(datatype = as.character(cfg$datatype), mode = "single",
                          db = normalizePath(DBDIR), population.size = FALSE),
           net = list())
  list(object = o, effective = eff, note = notes)
}

