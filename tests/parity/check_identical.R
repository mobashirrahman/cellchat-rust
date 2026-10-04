# Compare serialized bytes, including signed zero and NaN payloads.
identical <- function(x, y, ...) base::identical(serialize(x, NULL, version=3L),
                                               serialize(y, NULL, version=3L))
# Differential gate: run the Rust-backed `computeCommunProb` and the pinned upstream one on
# the same objects and require `identical()`.
#
# This is the acceptance criterion from the objective: not "close", but `identical()` on
# `net$prob`, `net$pval` and `options$parameter`. Wall-clock `options$run.time` is
# necessarily excluded and is reported separately.
#
# Usage:  R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R
suppressWarnings(suppressMessages({library(Matrix); library(collapse); library(dplyr)}))

CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
DBDIR <- Sys.getenv("CELLCHATRS_DB", "tests/fixtures/db_human")

suppressWarnings(suppressMessages(library(CellChat)))

## ---------------------------------------------------------------- minimal S4 object
setClass("MiniCellChat", representation(
  data.signaling = "Matrix",
  LR = "list", LRsig = "data.frame", DB = "ANY",
  idents = "factor", options = "list", net = "ANY"
))

E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1], E)

## ---------------------------------------------------------------- the fixture
## Same shape as `gen_prob_golden.R`, so a failure here is a real regression and not a
## difference in the test data. Sourced from the generator's preamble rather than
## duplicated.
src <- readLines("tests/parity/gen_prob_golden.R")
stop_at <- grep("^q <- file", src)[1]
eval(parse(text = paste(src[seq_len(stop_at - 1)], collapse = "\n")))

## Snapshot the RNA fixture's inputs under private names. Later blocks in this file reassign the
## globals `LRsig` and friends -- `computeCommunProbPathway` and `subsetCommunication` each build
## their own small `LRsig` in the *global* environment -- so anything added later that reaches for
## `LRsig` silently gets a 1-row, ligand-less frame from whichever block ran last. That is how the
## spatial block came to pass `nrow(pairLR.use) = 1, ligand = 0` and fail with a kernel length
## complaint that named neither the fixture nor the clobber.
RNA_LRsig <- LRsig
RNA_data <- data.signaling
RNA_group <- cell_group

mk <- function() {
  new("MiniCellChat",
      data.signaling = data.signaling,
      LR = list(LRsig = LRsig), LRsig = LRsig,
      DB = list(complex = DB$complex, cofactor = DB$cofactor),
      idents = cell_group,
      ## `mode` is read by `subsetCommunication`, which is not in the pinned version of
      ## `computeCommunProb`'s own preconditions but is set by `createCellChat`.
      options = list(datatype = "RNA", mode = "single", db = normalizePath(DBDIR)))
}

## ---------------------------------------------------------------- configurations
configs <- list(
  list(name = "tri_ps0",    type = "triMean",         trim = 0.1, pop = FALSE, nboot = 5, seed = 1L, Kh = 0.5, n = 1),
  list(name = "tri_ps1",    type = "triMean",         trim = 0.1, pop = TRUE,  nboot = 5, seed = 1L, Kh = 0.5, n = 1),
  list(name = "trim_ps0",   type = "truncatedMean",   trim = 0.2, pop = FALSE, nboot = 7, seed = 2L, Kh = 0.5, n = 1),
  list(name = "thresh_ps1", type = "thresholdedMean", trim = 0.3, pop = TRUE,  nboot = 4, seed = 3L, Kh = 0.5, n = 1),
  list(name = "median_ps0", type = "median",          trim = 0.1, pop = FALSE, nboot = 6, seed = 4L, Kh = 0.5, n = 1),
  list(name = "hill2",      type = "triMean",         trim = 0.1, pop = FALSE, nboot = 5, seed = 1L, Kh = 0.5, n = 2),
  list(name = "kh1e3",      type = "triMean",         trim = 0.1, pop = TRUE,  nboot = 8, seed = 7L, Kh = 1e3, n = 1),
  ## `match.arg` partial matching (R-ism 9): "tri" is a unique prefix of "triMean", and R
  ## records the *matched* name in options$parameter$type.mean, not the argument.
  list(name = "matcharg",   type = "tri",             trim = 0.1, pop = FALSE, nboot = 3, seed = 9L, Kh = 0.5, n = 1),
  ## Ambiguous prefix: "t" matches three of the four choices, so `match.arg` must error --
  ## and must error the same way from Rust.
  list(name = "ambiguous",  type = "t",               trim = 0.1, pop = FALSE, nboot = 3, seed = 9L, Kh = 0.5, n = 1)
)

run_up <- function(obj, cfg) {
  suppressWarnings(suppressMessages(cellchatrs_upstream_computeCommunProb(
    obj, type = cfg$type, trim = cfg$trim, population.size = cfg$pop,
    nboot = cfg$nboot, seed.use = cfg$seed, Kh = cfg$Kh, n = cfg$n)))
}
run_rs <- function(obj, cfg) {
  suppressWarnings(suppressMessages(computeCommunProb(
    obj, type = cfg$type, trim = cfg$trim, population.size = cfg$pop,
    nboot = cfg$nboot, seed.use = cfg$seed, Kh = cfg$Kh, n = cfg$n)))
}

## 1-based linear index -> array subscripts, as `which(arr.ind = TRUE)` reports them.
arrayInd <- function(ind, dims) {
  subs <- vector("list", length(dims))
  for (i in seq_along(dims)) {
    subs[[i]] <- ((ind - 1) %% dims[i]) + 1L
    ind <- (ind - 1L) %/% dims[i] + 1L
  }
  do.call(cbind, subs)
}

`%||%` <- function(a, b) if (is.null(a)) b else a

fails <- 0L
## ---------------------------------------------------------------------------------------
## Precondition: the shim has to actually be the thing under test.
##
## This gate compares the shim against pinned upstream, so it is trivially green if every call
## routes to upstream -- and that is not hypothetical. `computeCommunProb` gated its Rust path on
## `nzchar(Sys.getenv("CELLCHATRS_FALLBACK", "0"))`, and the unset default is the string `"0"`,
## which has four characters, so the predicate was **always true**: the escape hatch was permanently
## engaged, every `computeCommunProb` call went to upstream, and the whole run compared upstream
## with itself while reporting `IDENTICAL`. Nothing in the output distinguished the two cases.
##
## So the escape hatch is now an explicit affirmative (`1`/`true`/`yes`/`on`) and this assertion
## runs before any comparison. `tests/parity/check_rust_path.R` is the complementary check: it
## makes the two paths distinguishable from outside the package, by pointing `options$db` at a
## directory that does not exist and requiring the Rust-only db error.
if (tolower(Sys.getenv("CELLCHATRS_FALLBACK", "0")) %in% c("1", "true", "yes", "on")) {
  stop("CELLCHATRS_FALLBACK is set: every shim call would be delegated to pinned upstream and\n",
       "this gate would compare upstream against itself. Unset it to run the real comparison.",
       call. = FALSE)
}

n_cmp <- 0L

## `mk2` clones a computed object, so the upstream and Rust calls see identical inputs.
mk2 <- function(net_obj) {
  o <- mk()
  o@net <- net_obj@net
  o
}


## An error is part of the contract: `match.arg`'s message must match byte for byte, and so
## must upstream's own aborts.
try_both <- function(cfg) {
  up_err <- tryCatch({ run_up(mk(), cfg); NULL }, error = conditionMessage)
  rs_err <- tryCatch({ run_rs(mk(), cfg); NULL }, error = conditionMessage)
  if (!is.null(up_err) || !is.null(rs_err)) {
    n_cmp <<- n_cmp + 1L
    ok <- identical(up_err, rs_err)
    cat(sprintf("%-10s error=%-5s  up: %s\n", cfg$name, ok,
                if (is.null(up_err)) "<no error>" else up_err))
    if (!ok) {
      fails <<- fails + 1L
      cat(sprintf("                        rs: %s\n",
                  if (is.null(rs_err)) "<no error>" else rs_err))
    }
    return(invisible(NULL))
  }
  set.seed(314159L)
  up <- run_up(mk(), cfg)
  up_rng <- .Random.seed
  up_next <- runif(4L)
  set.seed(314159L)
  rs <- run_rs(mk(), cfg)
  rs_rng <- .Random.seed
  rs_next <- runif(4L)
  ok_prob <- identical(up@net$prob, rs@net$prob)
  ok_pval <- identical(up@net$pval, rs@net$pval)
  ok_par  <- identical(up@options$parameter, rs@options$parameter)
  ok_dim  <- identical(dimnames(up@net$prob), dimnames(rs@net$prob))
  ok_rng <- identical(up_rng, rs_rng) && identical(up_next, rs_next)
  n_cmp <<- n_cmp + 5L
  cat(sprintf("%-10s prob=%-5s pval=%-5s parameter=%-5s dimnames=%-5s rng=%-5s  run.time up=%.4fs rs=%.4fs\n",
              cfg$name, ok_prob, ok_pval, ok_par, ok_dim, ok_rng,
              up@options$run.time, rs@options$run.time))
  if (!ok_prob) {
    fails <<- fails + 1L
    i <- which(abs(up@net$prob - rs@net$prob) > 0)
    cat(sprintf("   %d of %d entries differ; first: up=%.17g rs=%.17g\n",
                length(i), length(up@net$prob), up@net$prob[i[1]], rs@net$prob[i[1]]))
  }
  if (!ok_pval) fails <<- fails + 1L
  if (!ok_dim)   fails <<- fails + 1L
  if (!ok_rng) {
    fails <<- fails + 1L
    cat("   R's post-call RNG state or next draws differ\n")
  }
  if (!ok_par) {
    fails <<- fails + 1L
    cat("   parameter objects:\n")
    d <- setdiff(names(up@options$parameter), names(rs@options$parameter))
    if (length(d)) cat("   parameter fields only in upstream:", paste(d, collapse = ", "), "\n")
    for (nm in intersect(names(up@options$parameter), names(rs@options$parameter))) {
      if (!identical(up@options$parameter[[nm]], rs@options$parameter[[nm]])) {
        cat(sprintf("   %s: up=%s rs=%s\n", nm,
                    paste(format(up@options$parameter[[nm]]), collapse = ","),
                    paste(format(rs@options$parameter[[nm]]), collapse = ",")))
      }
    }
  }
  invisible(NULL)
}

invisible(lapply(configs, try_both))

## Alternate sample generators are delegated wholesale; the port's permutation implementation is
## intentionally pinned to R's default Mersenne-Twister + Rejection stream.
check_alternate_rng <- function() {
  old_kind <- base::RNGkind()
  on.exit(do.call(base::RNGkind, as.list(old_kind)), add = TRUE)
  base::RNGkind(kind = "L'Ecuyer-CMRG", normal.kind = old_kind[[2L]],
                sample.kind = old_kind[[3L]])
  cfg <- configs[[1L]]
  set.seed(2718L)
  up <- run_up(mk(), cfg)
  up_rng <- .Random.seed
  up_next <- runif(4L)
  set.seed(2718L)
  rs <- run_rs(mk(), cfg)
  rs_rng <- .Random.seed
  rs_next <- runif(4L)
  identical(up@net, rs@net) && identical(up_rng, rs_rng) && identical(up_next, rs_next)
}
alternate_rng_ok <- check_alternate_rng()
n_cmp <<- n_cmp + 1L
cat(sprintf("alternate sample.kind delegated with matching output and RNG state: %s\n",
            alternate_rng_ok))
if (!alternate_rng_ok) fails <<- fails + 1L

## ---------------------------------------------------------------- computeAveExpr
## Fed the real `data` slot, with the same three `type` choices plus the `match.arg` error.
## `features` uses R's `intersect` order (caller's order, duplicates collapsed, absent names
## dropped), which is the part a naive port gets wrong.
setClass("MiniAve2", representation(
  data = "matrix", data.signaling = "ANY", DB = "ANY", idents = "factor",
  meta = "list", var.features = "list", options = "list"
))
mk_ave <- function() {
  new("MiniAve2", data = expr, DB = list(complex = DB$complex, cofactor = DB$cofactor),
      idents = cell_group, meta = list(), var.features = list(),
      options = list(datatype = "RNA"))
}
ave_specs <- list(
  list(name = "tri",    type = "triMean",       trim = NULL, features = NULL),
  list(name = "trim01", type = "truncatedMean", trim = 0.1,  features = NULL),
  list(name = "trim05", type = "truncatedMean", trim = 0.5,  features = NULL),
  list(name = "median", type = "median",        trim = NULL, features = NULL),
  list(name = "featsel", type = "triMean",      trim = NULL,
       features = c(genes[10], genes[1], genes[10], "NOPE", genes[2])),
  list(name = "matcharg", type = "tri",         trim = NULL, features = NULL)
)
for (spec in ave_specs) {
  a <- suppressWarnings(suppressMessages(cellchatrs_upstream_computeAveExpr(
    mk_ave(), features = spec$features, type = spec$type, trim = spec$trim,
    slot.name = "data")))
  b <- suppressWarnings(suppressMessages(computeAveExpr(
    mk_ave(), features = spec$features, type = spec$type, trim = spec$trim,
    slot.name = "data")))
  n_cmp <<- n_cmp + 4L
  ok_dim <- identical(dim(a), dim(b))
  ok_rn <- identical(rownames(a), rownames(b))
  ok_cn <- identical(colnames(a), colnames(b))
  ok_v <- isTRUE(all.equal(as.numeric(a), as.numeric(b), tolerance = 0))
  cat(sprintf("%-14s dim=%-5s rownames=%-5s colnames=%-5s values=%-5s\n",
              paste0("ave_", spec$name), ok_dim, ok_rn, ok_cn, ok_v))
  if (!(ok_dim && ok_rn && ok_cn && ok_v)) {
    fails <<- fails + 1L
    if (!ok_rn) cat("   upstream rownames:", paste(head(rownames(a), 5), collapse = ","),
                    "\n   rust rownames:    ", paste(head(rownames(b), 5), collapse = ","), "\n")
    if (!ok_v) cat("   max abs diff:", max(abs(as.numeric(a) - as.numeric(b)), na.rm = TRUE), "\n")
  }
}
## `match.arg` must reject an ambiguous prefix with upstream's exact message.
ave_err <- tryCatch({ computeAveExpr(mk_ave(), type = "t"); NULL },
                    error = function(e) conditionMessage(e))
ave_err_up <- tryCatch({ cellchatrs_upstream_computeAveExpr(mk_ave(), type = "t"); NULL },
                       error = function(e) conditionMessage(e))
n_cmp <<- n_cmp + 1L
cat(sprintf("ave_argerr      equal=%-5s  %s\n", identical(ave_err, ave_err_up),
            if (is.null(ave_err)) "<no error>" else ave_err))
if (!identical(ave_err, ave_err_up)) fails <<- fails + 1L

## ---------------------------------------------------------------- subsetDB
for (spec in list(
  list(name = "default", search = NULL, non_protein = FALSE),
  list(name = "nonprotein", search = NULL, non_protein = TRUE),
  list(name = "contact_only", search = "Cell-Cell Contact", non_protein = FALSE),
  list(name = "empty", search = character(0), non_protein = FALSE)
)) {
  a <- suppressMessages(cellchatrs_upstream_subsetDB(
    list(interaction = DB$interaction), search = spec$search, key = "annotation",
    non_protein = spec$non_protein))
  b <- suppressMessages(subsetDB(
    list(interaction = DB$interaction), search = spec$search, key = "annotation",
    non_protein = spec$non_protein))
  n_cmp <<- n_cmp + 2L
  ok_n <- identical(nrow(a$interaction), nrow(b$interaction))
  ok_i <- identical(a$interaction$interaction_name, b$interaction$interaction_name)
  cat(sprintf("%-14s nrow=%-5d names=%-5s\n", paste0("subdb_", spec$name),
              nrow(a$interaction), ok_i))
  if (!ok_n || !ok_i) fails <<- fails + 1L
}
## An unknown key must give upstream's message, byte for byte.
k_up <- tryCatch(cellchatrs_upstream_subsetDB(list(interaction = DB$interaction),
                                              search = "Secreted Signaling", key = "nope"),
                 error = function(e) conditionMessage(e))
k_rs <- tryCatch(subsetDB(list(interaction = DB$interaction),
                          search = "Secreted Signaling", key = "nope"),
                 error = function(e) conditionMessage(e))
n_cmp <<- n_cmp + 1L
cat(sprintf("subdb_keyerr    equal=%-5s\n", identical(k_up, k_rs)))
if (!identical(k_up, k_rs)) { fails <<- fails + 1L; cat("   up:", k_up, "\n   rs:", k_rs, "\n") }

## ---------------------------------------------------------------- downstream functions
## `aggregateNet` and `subsetCommunication` are fed a *real* `net` produced by the upstream
## `computeCommunProb`, so this is a genuine differential test of the downstream functions
## on data the upstream kernel itself produced -- not a synthetic array.
## ---------------------------------------------------------------------------------
downstream_specs <- list(
  list(name = "agg_thresh005", thresh = 0.05),
  list(name = "agg_thresh01",  thresh = 0.01),
  list(name = "agg_thresh1",   thresh = 1)
)
base_cfg <- configs[[1]]
up0 <- run_up(mk(), base_cfg)

for (spec in downstream_specs) {
  a <- suppressWarnings(suppressMessages(cellchatrs_upstream_aggregateNet(mk2(up0), thresh = spec$thresh)))
  b <- suppressWarnings(suppressMessages(aggregateNet(mk2(up0), thresh = spec$thresh)))
  n_cmp <<- n_cmp + 3L
  ok_c <- identical(a@net$count, b@net$count)
  ok_w <- identical(a@net$weight, b@net$weight)
  ok_d <- identical(dimnames(a@net$count), dimnames(b@net$count))
  cat(sprintf("%-14s count=%-5s weight=%-5s dimnames=%-5s\n", spec$name, ok_c, ok_w, ok_d))
  if (!ok_c) { fails <<- fails + 1L
    cat("   max rel diff:", max(abs(a@net$count - b@net$count)), "\n") }
  if (!ok_w) { fails <<- fails + 1L
    d <- abs(a@net$weight - b@net$weight); d <- d[d > 0]
    cat("   weight max rel diff:", if (length(d)) max(d) else 0, "\n")
    w <- which(a@net$weight != b@net$weight, arr.ind = TRUE)
    if (nrow(w)) cat(sprintf("   first cell [%d,%d]: up=%.17g rs=%.17g\n", w[1,1], w[1,2],
                              a@net$weight[w[1,1], w[1,2]], b@net$weight[w[1,1], w[1,2]])) }
  if (!ok_d) fails <<- fails + 1L
}

## A Net with at least one surviving row, so `subsetCommunication` returns a table.
sel <- up0
sel@net$prob[1, 1, ] <- ifelse(seq_len(dim(sel@net$prob)[3]) %% 2 == 1, 0.9, 0.0)
sel@net$pval[1, 1, ] <- ifelse(seq_len(dim(sel@net$pval)[3]) %% 2 == 1, 0.01, 0.9)
sel@net$prob[2, 3, ] <- 0.7
sel@net$pval[2, 3, ] <- 0.02

for (case in list(
  list(name = "sub_all",      s = NULL, t = NULL),
  list(name = "sub_src_g1",   s = "g1", t = NULL),
  list(name = "sub_src_tgt",  s = c("g1", "g2"), t = "g3")
)) {
  da <- suppressWarnings(suppressMessages(cellchatrs_upstream_subsetCommunication(
    object = mk2(sel), sources.use = case$s, targets.use = case$t, thresh = 0.05)))
  db <- suppressWarnings(suppressMessages(subsetCommunication(
    object = mk2(sel), sources.use = case$s, targets.use = case$t, thresh = 0.05)))
  n_cmp <<- n_cmp + 2L
  ok_cols <- identical(colnames(da), colnames(db))
  ok_lv <- identical(levels(da$source), levels(db$source)) &&
          identical(levels(da$target), levels(db$target))
  same <- function(col) isTRUE(all.equal(da[[col]], db[[col]], tolerance = 0))
  ok_num <- same("prob") && same("pval")
  ok_chr <- all(vapply(c("source", "target", "ligand", "receptor", "interaction_name",
                         "annotation", "pathway_name"),
                       function(col) identical(as.character(da[[col]]), as.character(db[[col]])), TRUE))
  ok_rn <- identical(rownames(da), rownames(db))
  cat(sprintf("%-14s nrow=%-3d cols=%-5s levels=%-5s num=%-5s chr=%-5s rownames=%-5s\n",
              case$name, nrow(da), ok_cols, ok_lv, ok_num, ok_chr, ok_rn))
  if (!ok_cols || !ok_lv || !ok_num || !ok_chr || !ok_rn) {
    fails <<- fails + 1L
    if (!ok_cols) cat("   upstream cols:", paste(colnames(da), collapse = ","),
                      "\n   rust cols:    ", paste(colnames(db), collapse = ","), "\n")
    if (!ok_lv) cat("   levels differ\n")
    if (!ok_num) cat("   prob/pval differ\n")
    if (!ok_chr) {
      cand <- c("source", "target", "ligand", "receptor", "interaction_name",
                "annotation", "pathway_name", "interaction_name_2", "evidence")
      for (col in cand) {
        ua <- as.character(da[[col]]); ub <- as.character(db[[col]])
        if (identical(ua, ub)) next
        ## NA-aware: `NA != NA` is NA, so a naive mask misses exactly the cells that
        ## differ while skipping over the NA==NA ones.
        bad <- which(!(is.na(ua) & is.na(ub)) & (is.na(ua) != is.na(ub) | (!is.na(ua) & ua != ub)))
        cat(sprintf("   %s: %d of %d differ; first row %s\n", col, length(bad), length(ua),
                    if (length(bad)) paste0(bad[1], ": up=", dQuote(ua[bad[1]]),
                                            " rs=", dQuote(ub[bad[1]])) else "-"))
      }
    }
    if (!ok_rn) cat("   rownames differ\n")
  }
}


## ------------------------------------------------------- identifyOverExpressedGenes
## `do.fast = FALSE` only: the presto branch is a different algorithm (and an optional
## C++ dependency), so it is a fallback rather than a parity claim.
##
## The gate is `identical()` on the *whole* `var.features` list, not on the marker table
## alone, because that is what the S4 object stores -- and because it is the only way the
## `0x1` collapsed frame and the `dplyr::filter`/named-`rowSums` row-name rules get
## checked at all.
de_specs <- list(
  list(name = "oeg_default",  pc = 0,   fc = 0,   p = 0.05,  pos = TRUE,  DE = TRUE,  mc = 10),
  list(name = "oeg_fc05",     pc = 0,   fc = 0.5, p = 0.1,   pos = TRUE,  DE = TRUE,  mc = 10),
  list(name = "oeg_pc01",     pc = 0.1, fc = 0,   p = 0.05,  pos = TRUE,  DE = TRUE,  mc = 10),
  list(name = "oeg_ptight",   pc = 0,   fc = 0,   p = 0.001, pos = TRUE,  DE = TRUE,  mc = 10),
  list(name = "oeg_twosided", pc = 0,   fc = 0,   p = 0.05,  pos = FALSE, DE = TRUE,  mc = 10),
  list(name = "oeg_node",     pc = 0.2, fc = 0,   p = 1,     pos = TRUE,  DE = FALSE, mc = 10),
  list(name = "oeg_mincells", pc = 0,   fc = 0,   p = 0.05,  pos = TRUE,  DE = TRUE,  mc = 1000)
)

## A 3-group expression matrix with the shapes the function has to survive: an all-zero
## row (a gene that can never be a marker, and whose Wilcoxon p-value is `NaN` because
## `SIGMA == 0`), a constant row, a half-zero row, and a row enriched in group 1.
## A *separate* class: `MiniCellChat` is (re)defined by the `gen_prob_golden.R` preamble
## evaluated above, so adding slots to it here is silently undone. And the slots are
## required: `identifyOverExpressedGenes`'s first statement is
## `if (!is.list(object@var.features)) stop("Please update your CellChat object via
## \`updateCellChat()\`")`, so without them *both* implementations stop with the same
## message -- a "passing" comparison that tests nothing.
setClass("MiniDE", representation(
  data.signaling = "ANY", LR = "list", LRsig = "data.frame", DB = "ANY",
  idents = "factor", options = "list", net = "ANY", var.features = "list", meta = "list"
))

set.seed(20240411)
NG <- 30L; NC <- 80L
de_genes <- sprintf("G%d", seq_len(NG))
de_group <- factor(c(rep("g1", 30), rep("g2", 30), rep("g3", 20))[sample.int(NC)])
de_X <- matrix(rpois(NG * NC, 3), nrow = NG, dimnames = list(de_genes, NULL))
de_X[1, ] <- 0
de_X[2, seq_len(NC / 2)] <- 0
de_X[3, ] <- 1
de_X[4, de_group == "g1"] <- 6
mkde <- function() {
  new("MiniDE", data.signaling = as(de_X, "dgCMatrix"), LR = list(LRsig = LRsig),
      LRsig = LRsig, DB = list(complex = DB$complex, cofactor = DB$cofactor),
      idents = de_group, options = list(datatype = "RNA", mode = "single"),
      var.features = list())
}

for (sp in de_specs) {
  args <- list(group.by = NULL, idents.use = NULL, invert = FALSE,
               features.name = "features", only.pos = sp$pos, features = NULL,
               return.object = TRUE, thresh.pc = sp$pc, thresh.fc = sp$fc, thresh.p = sp$p,
               do.DE = sp$DE, do.fast = FALSE, min.cells = sp$mc)
  ea <- tryCatch(suppressWarnings(suppressMessages(do.call(
        cellchatrs_upstream_identifyOverExpressedGenes, c(list(mkde()), args)))),
    error = function(e) conditionMessage(e))
  eb <- tryCatch(suppressWarnings(suppressMessages(do.call(
        identifyOverExpressedGenes, c(list(mkde()), args)))),
    error = function(e) conditionMessage(e))
  n_cmp <<- n_cmp + 2L
  if (is.character(ea) || is.character(eb)) {
    ok <- identical(ea, eb)
    ## `ea` may be a returned S4 object rather than a message, and `paste` on an S4 object
    ## is "no method for coercing this S4 class to a vector" -- which would replace the
    ## real diagnostic with a confusing one.
    brief <- function(x) if (is.character(x)) substr(x, 1, 70) else paste0("<", class(x)[1], ">")
    cat(sprintf("%-14s error=%-5s  up: %s | rs: %s\n", sp$name, ok, brief(ea), brief(eb)))
    if (!ok) fails <<- fails + 1L
    next
  }
  ok_vf <- identical(ea@var.features, eb@var.features)
  ok_ret <- identical(ea, eb)
  ia <- ea@var.features[["features.info"]]
  cat(sprintf("%-14s var.features=%-5s object=%-5s  info=%dx%d cols=%d\n", sp$name,
              ok_vf, ok_ret, nrow(ia), ncol(ia), ncol(ia)))
  if (!ok_vf) {
    fails <<- fails + 1L
    for (nm in union(names(ea@var.features), names(eb@var.features))) {
      if (identical(ea@var.features[[nm]], eb@var.features[[nm]])) next
      cat("   differs:", nm, "\n")
      print(all.equal(ea@var.features[[nm]], eb@var.features[[nm]]))
    }
  }
  if (!ok_ret) fails <<- fails + 1L
}

## ------------------------------------------------- identifyOverExpressedGenes, group.dataset
##
## `group.dataset` was filed as blocked on presto and is not: upstream honours it in **both**
## halves of `if (do.fast)` -- the presto half at `utilities.R:429-484` and the Wilcoxon half at
## `:512-520`. Only the first needs presto, so `do.fast = FALSE` with `group.dataset` set is
## testable, and it is now ported.
##
## What `group.dataset` changes is *only* how the two cell sets are chosen per group:
##
## ```r
## cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
## cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
## ```
##
## or, with `group.DE.combined = TRUE`, both pooled across groups with no `labels` term at all.
## The percentage filter, `mean.fxn`, `wilcox.test` and the Bonferroni multiplier are untouched,
## which is why the kernel takes the selection as a parameter rather than growing a flag.
##
## The cases below cover the three things that are easy to get wrong and that a "does it work"
## test would not see:
##
##   * **the `toString` collapse.** `labels.dataset[labels.dataset != pos.dataset] <-
##     toString(setdiff(unique(labels.dataset), pos.dataset))` puts *every* non-positive dataset
##     into one level, so a four-dataset comparison still has two levels, one of them named
##     `"D2, D3, D4"`. A port that keeps the datasets separate produces a different
##     `datasets` factor and a different row order.
##   * **the empty `cell.use2`.** When a group's cells are all in the positive dataset,
##     `cell.use2` is empty; `pct.2` is then `NaN`, every feature is dropped by
##     `which(alpha.min > thresh.pc)`, and `markers.all` is still the bare `data.frame()`.
##     Upstream then fails on `-markers.all$logFC` -- *after* adding a `datasets` column to a
##     0 x 0 frame. Substituting the collapsed 1-column frame returns a zero-row marker table
##     where upstream returns nothing at all.
##   * **the row names after the reorder.** `markers.all[order(datasets, pvalues, -logFC), ]`
##     keeps them, and after `rbind` they are R's duplicate-uniquified feature names
##     (`gene7`, `gene7.1`, `gene7.2`), not the feature column's repeated `gene7`.
de_ds_specs <- list(
  list(name = "oegds_2ds",        conds = c("D1", "D2"),         pos = "D1", comb = FALSE, seed = 101),
  list(name = "oegds_2ds_comb",   conds = c("D1", "D2"),         pos = "D1", comb = TRUE, seed = 102),
  list(name = "oegds_2ds_posD2",  conds = c("D1", "D2"),         pos = "D2", comb = FALSE, seed = 103),
  list(name = "oegds_4ds",        conds = c("D1", "D2", "D3", "D4"), pos = "D1", comb = FALSE, seed = 104),
  list(name = "oegds_4ds_comb",   conds = c("D1", "D2", "D3", "D4"), pos = "D1", comb = TRUE, seed = 105),
  list(name = "oegds_4ds_posD3",  conds = c("D1", "D2", "D3", "D4"), pos = "D3", comb = FALSE, seed = 106),
  ## Every cell of every group is in the positive dataset, so `cell.use2` is empty and
  ## upstream raises. Aligned labels are what produce that, and it is the case a shuffled
  ## fixture never reaches.
  list(name = "oegds_aligned_empty", conds = c("ctrl", "treatA", "treatB"), pos = "ctrl",
       comb = FALSE, align = TRUE, seed = 107),
  list(name = "oegds_aligned_comb",  conds = c("ctrl", "treatA", "treatB"), pos = "ctrl",
       comb = TRUE, align = TRUE, seed = 108)
)
## Seeded per spec. The first version of this fixture drew from the ambient RNG, so the object
## depended on how many random numbers earlier blocks had consumed -- a failure that cannot be
## reproduced, and a *pass* that depends on the block's position in the file.
mkds <- function(sp) {
  set.seed(sp$seed)
  nc <- 80L; ng <- 3L
  cond <- rep(sp$conds, length.out = nc)
  grp <- factor(if (isTRUE(sp$align)) rep(paste0("g", seq_len(ng)), length.out = nc)
                else sample(rep(paste0("g", seq_len(ng)), length.out = nc)))
  nf <- 40L
  X <- matrix(rbinom(nf * nc, 1, 0.12), nf, nc,
              dimnames = list(sprintf("gene%d", seq_len(nf)), sprintf("c%d", seq_len(nc))))
  for (b in seq_len(ng)) {
    blk <- ((b - 1L) * 12L + 1L):(b * 12L)
    X[blk, grp == paste0("g", b) & cond == sp$conds[1]] <- 9
  }
  new("MiniDE", data.signaling = as(X, "dgCMatrix"), LR = list(LRsig = LRsig),
      LRsig = LRsig, DB = list(complex = DB$complex, cofactor = DB$cofactor),
      idents = grp, options = list(datatype = "RNA", mode = "single"),
      var.features = list(features = rownames(X)),
      meta = list(cond = cond))
}
for (sp in de_ds_specs) {
  base <- list(group.dataset = "cond", pos.dataset = sp$pos,
               group.DE.combined = sp$comb)
  for (variant in list(list(), list(only.pos = FALSE), list(thresh.fc = 99),
                       list(thresh.p = 0.99), list(idents.use = "g1"))) {
    ## `modifyList`, not `c`. `c(list(thresh.fc = 0), list(thresh.fc = 99))` is a list of
    ## **length 2** with two `thresh.fc` names, and `do.call` then answers `formal argument
    ## "thresh.fc" matched by multiple actual arguments` -- on *both* sides, which the
    ## `error=TRUE` check scores as a pass. Nine "matching errors" that test nothing, and the
    ## real `thresh.fc = 99` case never ran.
    args <- modifyList(
      list(group.by = NULL, invert = FALSE, features.name = "features", features = NULL,
           return.object = TRUE, thresh.pc = 0, thresh.fc = 0, thresh.p = 0.05,
           do.DE = TRUE, do.fast = FALSE, min.cells = 10),
      c(base, variant))
    ea <- tryCatch(suppressWarnings(suppressMessages(do.call(
          cellchatrs_upstream_identifyOverExpressedGenes, c(list(mkds(sp)), args)))),
      error = function(e) conditionMessage(e))
    eb <- tryCatch(suppressWarnings(suppressMessages(do.call(
          identifyOverExpressedGenes, c(list(mkds(sp)), args)))),
      error = function(e) conditionMessage(e))
    n_cmp <<- n_cmp + 1L
    lbl <- sprintf("%s%s", sp$name,
                   if (length(variant)) paste0("/", names(variant)[1L], variant[[1L]]) else "")
    if (is.character(ea) || is.character(eb)) {
      ok <- identical(ea, eb)
      cat(sprintf("%-30s error=%-5s  %s\n", lbl, ok,
                  if (is.character(ea)) ea else paste0("<", class(ea)[1L], ">")))
      if (!ok) {
        fails <<- fails + 1L
        cat(sprintf("%-30s   port: %s\n", "",
                    if (is.character(eb)) eb else paste0("<", class(eb)[1L], ">")))
      }
    } else {
      ok <- identical(ea@var.features, eb@var.features)
      cat(sprintf("%-30s var.features=%-5s  rows=%d\n", lbl, ok,
                  nrow(ea@var.features[["features.info"]])))
      if (!ok) {
        fails <<- fails + 1L
        for (nm in union(names(ea@var.features), names(eb@var.features))) {
          if (identical(ea@var.features[[nm]], eb@var.features[[nm]])) next
          cat("   differs:", nm, "\n")
          print(all.equal(ea@var.features[[nm]], eb@var.features[[nm]]))
        }
      }
    }
  }
}
## `pos.dataset` that names nothing in `group.dataset`: upstream `cat()`s the available names
## and then calls a **bare** `stop()`, so the error message is empty. A shim that writes a
## helpful message is a divergence, and a shim that writes a helpful *warning* is worse.
ds_bad <- tryCatch(suppressWarnings(suppressMessages(CellChat::cellchatrs_upstream_identifyOverExpressedGenes(
  mkds(de_ds_specs[[1]]), group.dataset = "cond", pos.dataset = "nope", do.fast = FALSE))),
  error = function(e) conditionMessage(e))
ds_bad_p <- tryCatch(suppressWarnings(suppressMessages(CellChat::identifyOverExpressedGenes(
  mkds(de_ds_specs[[1]]), group.dataset = "cond", pos.dataset = "nope", do.fast = FALSE))),
  error = function(e) conditionMessage(e))
n_cmp <<- n_cmp + 1L
ok_bad <- identical(ds_bad, ds_bad_p) && identical(ds_bad, "")
cat(sprintf("%-30s bare stop()=%-5s  (message must be empty)\n", "oegds_bad_pos", ok_bad))
if (!ok_bad) {
  fails <<- fails + 1L
  cat("   up:", dQuote(ds_bad), " port:", dQuote(ds_bad_p), "\n")
}

## `do.fast = TRUE` must *fall back* to upstream -- which then requires presto and stops.
## The point of the check is that the Rust kernel is not silently used for a different
## algorithm, so the expected outcome is upstream's own presto `stop()`, not a result.
fb_up <- tryCatch(suppressWarnings(suppressMessages(
  cellchatrs_upstream_identifyOverExpressedGenes(mkde(), do.fast = TRUE, thresh.p = 0.05))),
  error = function(e) conditionMessage(e))
fb_rs <- tryCatch(suppressWarnings(suppressMessages(
  identifyOverExpressedGenes(mkde(), do.fast = TRUE, thresh.p = 0.05))),
  error = function(e) conditionMessage(e))
n_cmp <<- n_cmp + 1L
## The shim must produce the *same* outcome as upstream for `do.fast = TRUE`, which without
## presto installed is upstream's own `stop()`. Comparing the two is stronger than asserting
## a substring: it also catches a shim that runs the Wilcoxon kernel and returns a result.
ok_fb <- identical(fb_up, fb_rs)
cat(sprintf("%-14s fell back=%-5s  %s\n", "oeg_fast", ok_fb,
            if (is.character(fb_rs)) substr(gsub("[[:space:]]+", " ", fb_rs), 1, 60)
            else paste0("<", class(fb_rs)[1], ">")))
if (!ok_fb) {
  fails <<- fails + 1L
  cat("   up:", if (is.character(fb_up)) fb_up else paste0("<", class(fb_up)[1], ">"), "\n")
  cat("   rs:", if (is.character(fb_rs)) fb_rs else paste0("<", class(fb_rs)[1], ">"), "\n")
}

## ------------------------------------------------------- computeCommunProbPathway
## The corpus is read from `tests/fixtures/pathway_inputs.tsv` -- the generator's own dump,
## and the same file `src/rust/crates/r-core/tests/pathway_parity.rs` reads, so a Rust-side and an
## R-side failure cannot be about different data.
##
## It is deliberately **not** obtained by `eval`-ing the generator's preamble: that preamble
## contains `q <- file("tests/fixtures/pathway_golden.txt", "wt")`, so evaluating it without
## the run loop truncates the golden corpus to zero bytes. It failed exactly that way once.
## Loading the dump is also this repo's standing rule -- a fixture that needs an RNG is
## written out, never re-derived.
pw_lines <- readLines("tests/fixtures/pathway_inputs.tsv")
pw_at <- function(tag, name) {
  hit <- grep(sprintf("^%s\t%s\t", tag, name), pw_lines)
  if (!length(hit)) return(NA_character_)
  sub(sprintf("^%s\t%s\t", tag, name), "", pw_lines[hit[1]])
}
setClass("MiniNet2", representation(net = "list", netP = "list", LR = "list", LRsig = "data.frame", idents = "factor", options = "list"))
mk_net <- function(prob, pval, levels, lr) {
  dim(prob) <- c(length(levels), length(levels), length(lr))
  dimnames(prob) <- list(levels, levels, lr)
  dim(pval) <- c(length(levels), length(levels), length(lr))
  new("MiniNet2", net = list(
    prob = prob, pval = pval, prob.dim = dim(prob),
    dimnames = list(list(source = levels, target = levels),
                    list(source = levels, target = levels),
                    list(interaction_name = lr))
  ))
}
specs <- list()
pw_i <- 1L
while (pw_i <= length(pw_lines)) {
  l <- pw_lines[pw_i]
  if (!startsWith(l, "spec\t")) { pw_i <- pw_i + 1L; next }
  f <- strsplit(l, "\t")[[1]]
  nm <- f[2]
  k <- as.integer(f[3]); nlr <- as.integer(f[4])
  lev <- strsplit(f[5], ",")[[1]]; lr <- strsplit(f[6], ",")[[1]]
  pathv <- strsplit(pw_at("pathways", nm), ",")[[1]]
  thr <- as.numeric(pw_at("thresh", nm))
  stopifnot(identical(pw_lines[pw_i + 3L], "prob"))
  pv <- as.numeric(strsplit(pw_lines[pw_i + 4L], " ")[[1]])
  stopifnot(identical(pw_lines[pw_i + 5L], "pval"))
  qv <- as.numeric(strsplit(pw_lines[pw_i + 6L], " ")[[1]])
  stopifnot(length(pv) == k * k * nlr, length(qv) == k * k * nlr)
  specs[[length(specs) + 1L]] <- list(name = nm, prob = pv, pval = qv, lev = lev, lr = lr,
                                      thresh = thr)
  pw_i <- pw_i + 7L
}

for (sp in specs) {
  o <- mk_net(sp$prob, sp$pval, sp$lev, sp$lr)
  nlr <- length(sp$lr)
  pwv <- strsplit(pw_at("pathways", sp$name), ",")[[1]]
  LRsig <- data.frame(
    interaction_name = sp$lr,
    pathway_name = vapply(seq_len(nlr), function(i) pwv[[((i - 1L) %/% 2L) + 1L]], character(1)),
    stringsAsFactors = FALSE
  )
  ea <- tryCatch(suppressWarnings(suppressMessages(
        cellchatrs_upstream_computeCommunProbPathway(object = o, pairLR.use = LRsig,
                                                    thresh = sp$thresh))),
    error = function(e) conditionMessage(e))
  eb <- tryCatch(suppressWarnings(suppressMessages(
        computeCommunProbPathway(object = o, pairLR.use = LRsig, thresh = sp$thresh))),
    error = function(e) conditionMessage(e))
  n_cmp <<- n_cmp + 4L
  if (is.character(ea) || is.character(eb)) {
    ok <- identical(ea, eb)
    cat(sprintf("pw_%-10s error=%-5s  %s\n", sp$name, ok,
                if (is.character(eb)) substr(gsub("[[:space:]]+", " ", eb), 1, 50) else "<result>"))
    if (!ok) { fails <<- fails + 1L
      cat("   up:", if (is.character(ea)) ea else "<result>",
          "\n   rs:", if (is.character(eb)) eb else "<result>", "\n") }
    next
  }
  ## Both return shapes: the S4 object (which writes `net$LRs` and `netP`) and the bare
  ## `list(pathways, prob)` for `object = NULL`. `identical()` on the whole object also
  ## covers the `dimnames` and the `0 x 0` empty-network case.
  bare_a <- suppressWarnings(suppressMessages(
    cellchatrs_upstream_computeCommunProbPathway(object = NULL, net = o@net,
                                                 pairLR.use = LRsig, thresh = sp$thresh)))
  bare_b <- suppressWarnings(suppressMessages(
    computeCommunProbPathway(object = NULL, net = o@net, pairLR.use = LRsig, thresh = sp$thresh)))
  ok_obj <- identical(ea, eb)
  ok_bare <- identical(bare_a, bare_b)
  ok_vp <- identical(ea@netP$pathways, eb@netP$pathways)
  ok_lr <- identical(ea@net$LRs, eb@net$LRs)
  cat(sprintf("pw_%-10s object=%-5s bare=%-5s pathways=%-5s LRs=%-5s  dim=%s\n", sp$name,
              ok_obj, ok_bare, ok_vp, ok_lr, paste(dim(eb@netP$prob), collapse = "x")))
  if (!(ok_obj && ok_bare && ok_vp && ok_lr)) fails <<- fails + 1L
}

## --------------------------------------------------------------- aggregateNet, filtered
## The filtered branch is broken upstream -- `stringr::str_split(key, "|")` takes a regex and
## `|` is alternation, so `a[, 1]` is `""` and `a[, 2]` is the key's first character, and the
## result is `0 x 0` or a `k x k` matrix of zeros for every filter. Reproduced, and the shape
## is checked against upstream on every case.
##
## The metadata is read from the generator's own dump, one `|`-separated line per case, so
## this gate and `src/rust/crates/r-core/tests/netfiltered_parity.rs` cannot disagree about the data.
## The generator is deliberately *not* `eval`-ed here: its preamble opens
## `tests/fixtures/netfiltered_golden.txt` with "wt", which would truncate the corpus.
nf <- readLines("tests/fixtures/netfiltered_inputs.tsv")
dec_filter <- function(kind, val) {
  switch(kind, none = NULL,
         vec = strsplit(val, ",", fixed = TRUE)[[1]],
         df = data.frame(interaction_name = strsplit(val, ",", fixed = TRUE)[[1]],
                         stringsAsFactors = FALSE),
         stop("unknown filter encoding"))
}
nf_cases <- list()
nf_i <- 1L
while (nf_i <= length(nf)) {
  if (!startsWith(nf[nf_i], "case\t")) { nf_i <- nf_i + 1L; next }
  f <- strsplit(nf[nf_i], "\t")[[1]]
  ## The trailing "|" is deliberate: the last filter value is empty for most cases, and R's
  ## \ drops a trailing empty field, which would silently shift every later field.
  m <- strsplit(paste0(f[3], "|"), "|", fixed = TRUE)[[1]][seq_len(14L)]
  if (length(m) != 14L) {
    cat("DIAG i=", nf_i, " nchar(line)=", nchar(nf[nf_i]), " len(m)=", length(m),
        " line=[", nf[nf_i], "]
", sep = "")
  }
  ## NB plain `stop(...)`, not `stopifnot(cond, "message", extra)`: `stopifnot` checks
  ## *every* argument and a string argument is not a length-1 logical, so
  ## `stopifnot(TRUE, "oops")` fails. Putting the text where a message belongs makes a
  ## guard that always trips, and the error names the string rather than the condition.
  if (length(m) != 14L) {
    stop("malformed netfiltered_inputs case line: ", nf[nf_i], call. = FALSE)
  }
  if (!identical(nf[nf_i + 1L], "prob")) {
    stop("expected a prob marker after ", nf[nf_i], call. = FALSE)
  }
  nlr <- as.integer(m[2])
  k <- as.integer(m[1])
  lev <- strsplit(m[3], ",", fixed = TRUE)[[1]]
  nf_cases[[length(nf_cases) + 1L]] <- list(
    name = f[2], k = k, levels = lev, lr = strsplit(m[4], ",", fixed = TRUE)[[1]],
    thresh = as.numeric(m[5]), remove = identical(m[6], "TRUE"),
    sources = dec_filter(m[7], m[8]), targets = dec_filter(m[9], m[10]),
    signaling = dec_filter(m[11], m[12]), pair_lr = dec_filter(m[13], m[14]),
    prob = as.numeric(strsplit(nf[nf_i + 2L], " ")[[1]]),
    pval = as.numeric(strsplit(nf[nf_i + 4L], " ")[[1]])
  )
  if (!identical(nf[nf_i + 3L], "pval")) {
    stop("expected a pval marker after the prob row of ", f[2], call. = FALSE)
  }
  nf_i <- nf_i + 5L
}

for (cs in nf_cases) {
  o <- mk_net(array(cs$prob, dim = c(cs$k, cs$k, length(cs$lr)),
                   dimnames = list(cs$levels, cs$levels, cs$lr)),
              array(cs$pval, dim = c(cs$k, cs$k, length(cs$lr))),
              cs$levels, cs$lr)
  ## The `LRsig` the generator used: two pathways, so a `signaling` filter has both a keeper
  ## and a discarder, and the full optional column set (a two-column table makes
  ## `subsetCommunication`'s `signaling` branch fail on `dplyr::select`).
  nlr <- length(cs$lr)
  o@idents <- factor(cs$levels, levels = cs$levels)
  o@options <- list(datatype = "RNA", mode = "single")
  o@LR$LRsig <- data.frame(
    interaction_name = cs$lr, interaction_name_2 = cs$lr,
    pathway_name = rep(c("PWA", "PWB"), length.out = nlr),
    ligand = paste0("L", seq_len(nlr)), receptor = paste0("R", seq_len(nlr)),
    annotation = rep("Secreted Signaling", nlr), evidence = NA_character_,
    stringsAsFactors = FALSE)
  o@LRsig <- o@LR$LRsig
  run <- function(f) suppressMessages(suppressWarnings(
    f(o, sources.use = cs$sources, targets.use = cs$targets, signaling = cs$signaling,
      pairLR.use = cs$pair_lr, remove.isolate = cs$remove, thresh = cs$thresh,
      return.object = TRUE)))
  ea <- tryCatch(run(cellchatrs_upstream_aggregateNet), error = function(e) conditionMessage(e))
  eb <- tryCatch(run(aggregateNet), error = function(e) conditionMessage(e))
  n_cmp <<- n_cmp + 1L
  if (is.character(ea) || is.character(eb)) {
    ok <- identical(ea, eb)
    cat(sprintf("nf_%-13s error=%-5s  %s\n", cs$name, ok,
                if (is.character(eb)) substr(gsub("[[:space:]]+", " ", eb), 1, 44) else "<result>"))
    if (!ok) fails <<- fails + 1L
    next
  }
  ok <- identical(ea, eb)
  cat(sprintf("nf_%-13s identical=%-5s dim=%s\n", cs$name, ok,
              paste(dim(eb@net$count), collapse = "x")))
  if (!ok) {
    fails <<- fails + 1L
    cat("   up dim:", paste(dim(ea@net$count), collapse = "x"),
        " rs dim:", paste(dim(eb@net$count), collapse = "x"), "\n")
    cat("   up rn:", paste(rownames(ea@net$count), collapse = ","),
        " rs rn:", paste(rownames(eb@net$count), collapse = ","), "\n")
  }
}

## ------------------------------------------------- subsetCommunication, the DEG branch
## The corpus is the generator's own dump, so a Rust-side and an R-side failure cannot be
## about different data. The internal function is called through the exported
## `subsetCommunication`, with `net` handed in as the data frame the corpus describes -- which
## is the only form in which the DEG columns can exist at all (the melt attaches only `LR`'s six
## metadata columns, so with an array `net` every threshold raises).
sb_lines <- readLines("tests/fixtures/subset_golden.txt")
sb_rec <- function(tag, name) {
  hit <- grep(sprintf("^%s\t%s\t", tag, name), sb_lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  sub(sprintf("^%s\t%s\t", tag, name), "", hit[length(hit)])
}
## One corpus cell. Both markers come back as **text**, never as R's logical `NA`/`NaN`: the
## numeric columns go through `as.numeric()` afterwards, which turns "NA" into `NA_real_` and
## "NaN" into `NaN` and so preserves the distinction, whereas a logical would collapse the two
## and make the column a logical vector.
sb_cell <- function(s) {
  s <- gsub("<NL>", "\n", s, fixed = TRUE)
  if (identical(s, "<NA>")) NA_character_ else s
}
## Row `i` of a `<tag>_row` record set. The index is part of the record key, so it has to be
## matched: taking the *last* match silently builds a table whose every row is the last row of
## the corpus, which still "works" and produces a plausible-looking reference.
sb_row <- function(tag, name, i) {
  hit <- grep(sprintf("^%s\t%s\t%d\t", tag, name, i), sb_lines, value = TRUE)
  if (!length(hit)) return(NULL)
  strsplit(hit[1], "\t")[[1]][-seq_len(3L)]
}
## Only the *comparisons* are reconstructed from the corpus here; the `slot`, the two group
## filters and the eight thresholds are re-declared, because the point of the installed-shim
## gate is to prove the binding marshals them faithfully, and re-deriving them from the corpus
## would make the gate agree with the corpus by construction.
sb_args <- list(
  no_thresholds = list(), no_thresholds_nodatasets = list(), all_na_row_no_threshold = list(),
  ligand_pvalues = list(ligand.pvalues = 0.05), ligand_pvalues_none = list(ligand.pvalues = 0),
  ligand_logfc_pos = list(ligand.logFC = 0.5), ligand_logfc_neg = list(ligand.logFC = -0.5),
  ligand_logfc_zero = list(ligand.logFC = 0), ligand_pct1 = list(ligand.pct.1 = 0.5),
  ligand_pct2 = list(ligand.pct.2 = 0.5), receptor_pvalues = list(receptor.pvalues = 0.02),
  receptor_logfc_pos = list(receptor.logFC = 0.5), receptor_logfc_neg = list(receptor.logFC = -0.5),
  receptor_pct1 = list(receptor.pct.1 = 0.5), receptor_pct2 = list(receptor.pct.2 = 0.5),
  all_eight_permissive = list(ligand.pvalues = 0.5, ligand.logFC = 0, ligand.pct.1 = 0.1,
                              ligand.pct.2 = 0.1, receptor.pvalues = 0.1, receptor.logFC = 0,
                              receptor.pct.1 = 0.1, receptor.pct.2 = 0.1),
  single_threshold_blanking = list(ligand.pvalues = 0.5),
  all_eight_tight = list(ligand.pvalues = 0.05, ligand.logFC = 0.25, ligand.pct.1 = 0.85,
                         ligand.pct.2 = 0.05, receptor.pvalues = 0.05, receptor.logFC = 0.4,
                         receptor.pct.1 = 0.75, receptor.pct.2 = 0.15),
  datasets_d1 = list(datasets = "d1"), datasets_d2_nodeg = list(datasets = "d2"),
  sources_targets = list(sources.use = "g1", targets.use = c("g2", "g10")),
  sources_missing = list(sources.use = "gNA"),
  err_ligand_pvalues_missing = list(ligand.pvalues = 0.05),
  err_ligand_logfc_missing = list(ligand.logFC = 0.5),
  err_receptor_pct2_missing = list(receptor.pct.2 = 0.5),
  err_datasets_missing = list(datasets = "d1"),
  err_two_absent_order = list(datasets = "d1", ligand.pvalues = 0.05, receptor.logFC = 0.5),
  err_empty_result = list(ligand.pct.1 = 2),
  netp_plain = list(), netp_after_deg = list(ligand.pvalues = 0.05, receptor.logFC = 0.4),
  netp_sources_targets = list(sources.use = c("g1", "g10"), targets.use = "g2"),
  netp_allna_row = list(), unknown_slot = list(), unknown_slot_deg = list(ligand.pvalues = 0.05)
)
sb_slot <- c(netp_plain = "netP", netp_after_deg = "netP", netp_sources_targets = "netP",
             netp_allna_row = "netP", unknown_slot = "foo", unknown_slot_deg = "foo")
mk_net_df <- function(name) {
  in_cols <- strsplit(sb_rec("in_cols", name), ",")[[1]]
  in_n <- as.integer(sb_rec("in_nrow", name))
  num_cols <- c("prob", "pval", "ligand.pvalues", "ligand.logFC", "ligand.pct.1",
                "ligand.pct.2", "receptor.pvalues", "receptor.logFC", "receptor.pct.1",
                "receptor.pct.2")
  cols <- lapply(in_cols, function(nm) {
    v <- vapply(seq_len(in_n), function(i) {
      f <- sb_row("in_row", name, i)
      sb_cell(f[match(nm, in_cols)])
    }, "")
    ## 17-digit text back to a double, or kept as character. `as.numeric("NA")` is `NA_real_`
    ## and `as.numeric("NaN")` is `NaN`, which is exactly the pair the function distinguishes.
    if (nm %in% num_cols) suppressWarnings(as.numeric(v)) else v
  })
  names(cols) <- in_cols
  as.data.frame(cols, stringsAsFactors = FALSE, check.names = FALSE)
}
setClass("MiniSub", representation(net = "list", idents = "factor", LR = "list", options = "list"))
LEVELS3 <- c("g1", "g2", "g10")
for (nm in names(sb_args)) {
  slot <- if (nm %in% names(sb_slot)) sb_slot[[nm]] else "net"
  df <- mk_net_df(nm)
  mkobj <- function(d) new("MiniSub", net = list(prob = array(0, c(1, 1, 0)),
    pval = array(0, c(1, 1, 0))), idents = factor(LEVELS3, levels = LEVELS3),
    LR = list(LRsig = data.frame()), options = list(mode = "single"))
  args <- c(list(object = mkobj(df), net = df, slot.name = slot), sb_args[[nm]])
  ## `do.call` for both sides: `f(**args)` needs a newer R than the one this runs on, and
  ## using it for only one side would make the comparison asymmetric.
  ## Upstream's warnings are *captured*, not suppressed. Suppressing them and then comparing
  ## the port's warning list against an always-empty one makes every warning a false failure --
  ## and, worse, a real divergence invisible.
  ea <- tryCatch({
    w <- character(0)
    v <- withCallingHandlers(
      suppressMessages(do.call(cellchatrs_upstream_subsetCommunication, args)),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    list(w = w, v = v)
  }, error = function(e) list(err = conditionMessage(e)))
  eb <- tryCatch({
    w <- character(0)
    v <- withCallingHandlers(suppressMessages(do.call(subsetCommunication, args)),
      warning = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleWarning") })
    list(w = w, v = v)
  }, error = function(e) list(err = conditionMessage(e)))
  n_cmp <<- n_cmp + 2L
  want_err <- !is.na(sb_rec("error", nm))
  if (want_err) {
    ok <- identical(ea$err, eb$err)
    cat(sprintf("sb_%-28s error=%-5s  rs=%s | up=%s\n", nm, ok,
                if (is.null(eb$err)) "<result>" else substr(gsub("[[:space:]]+", " ", eb$err), 1, 30),
                if (is.null(ea$err)) "<result>" else substr(gsub("[[:space:]]+", " ", ea$err), 1, 30)))
    if (!ok) { fails <<- fails + 1L
      cat("   up:", if (is.null(ea$err)) "<result>" else ea$err,
          "\n   rs:", if (is.null(eb$err)) "<result>" else eb$err, "\n") }
    next
  }
  ok_v <- identical(ea$v, eb$v)
  ok_rn <- identical(rownames(ea$v), rownames(eb$v))
  ## Upstream's own dplyr/tidyselect deprecations are dependency chatter, not CellChat's; the
  ## port does not use dplyr, so it cannot emit another package's version-specific text. The
  ## one CellChat warning must match.
  cellchat_warn <- "No significant signaling interactions are inferred!"
  ok_w <- identical(any(eb$w == cellchat_warn), any(ea$w == cellchat_warn))
  if (!ok_w) cat("   warnings up:", length(ea$w), " port:", length(eb$w),
                 " port has the empty-table warning:", any(eb$w == cellchat_warn), "\n")
  cat(sprintf("sb_%-28s value=%-5s rownames=%-5s warn=%-5s (%d x %d)\n", nm, ok_v, ok_rn, ok_w,
              nrow(eb$v), ncol(eb$v)))
  if (!ok_v) {
    fails <<- fails + 1L
    d <- setdiff(names(ea$v), names(eb$v)); if (length(d)) cat("   only upstream:", d, "\n")
    d <- setdiff(names(eb$v), names(ea$v)); if (length(d)) cat("   only port   :", d, "\n")
    ## Report the *first* difference of any kind -- a column value, a column class, or an
    ## attribute -- because `identical()` failing with no differing column is the confusing
    ## case, and it means the difference is in the column *types*.
    for (cl in intersect(names(ea$v), names(eb$v))) {
      if (!identical(ea$v[[cl]], eb$v[[cl]])) {
        cat("   column", cl, "differs (class", class(ea$v[[cl]])[1], "vs",
            class(eb$v[[cl]])[1], ")\n")
        cat("   upstream:"); print(utils::head(ea$v[[cl]], 4))
        cat("   port    :"); print(utils::head(eb$v[[cl]], 4))
        break
      }
    }
    if (!identical(attributes(ea$v), attributes(eb$v))) {
      cat("   attributes differ:", paste(names(attributes(ea$v)), collapse = ","), "vs",
          paste(names(attributes(eb$v)), collapse = ","), "\n")
      for (at in union(names(attributes(ea$v)), names(attributes(eb$v)))) {
        if (!identical(attr(ea$v, at), attr(eb$v, at)))
          cat("     ", at, ":", paste(utils::head(attr(ea$v, at), 4), collapse = ","), "vs",
              paste(utils::head(attr(eb$v, at), 4), collapse = ","), "\n")
      }
    }
  }
  if (!ok_rn || !ok_w) fails <<- fails + 1L
}

## --------------------------------------------- computeCommunProb, the spatial branch
##
## Compare against the unmodified pinned spatial function, including its Annoy queries.
sp_env <- get("cellchatrs_upstream_cached", envir = asNamespace("CellChat"))()

setClass("SpChat", representation(
  data.signaling = "Matrix", LRsig = "data.frame", DB = "ANY", idents = "factor",
  meta = "data.frame", images = "list", LR = "list", net = "list", options = "list"))

## Both paths resolve the supplied object database.
SP_DB <- DB
## The RNA fixture's own cells, laid out spatially: `nlev` blocks `sep` apart with the cells inside
## a block spread over `per` columns. The number of cells and the group sizes come from
## `cell_group`, so the expression side is exactly what the passing RNA configurations already
## compare.
sp_obj <- function(per = 8L, sep = 25, samples = 2L) {
  ds <- RNA_data
  ncell <- ncol(ds)
  ## Every level of `cell_group` gets a block. Subsetting to the first `nlev` levels while keeping
  ## all the cells puts `NA` coordinates in the matrix for the dropped groups, and upstream then
  ## dies with "subscript out of bounds" from somewhere unrelated -- a fixture that looks like a
  ## crash in the spatial branch and is a crash in the layout.
  levs <- levels(RNA_group)
  nlev <- length(levs)
  group <- droplevels(RNA_group)
  co <- matrix(0, nrow = ncell, ncol = 2)
  for (i in seq_len(ncell)) {
    g <- match(as.character(group[i]), levs)
    k <- (i - 1) %% (ncell / nlev)
    co[i, ] <- c((g - 1) * sep + (k %% per) * 0.3, (k %/% per) * 0.25)
  }
  rownames(co) <- colnames(ds)
  new("SpChat", data.signaling = ds, LRsig = RNA_LRsig,
      DB = list(complex = SP_DB$complex, cofactor = SP_DB$cofactor),
      idents = group,
      meta = data.frame(samples = factor(rep(paste0("s", seq_len(samples)),
                                             each = ncell / samples),
                                        levels = paste0("s", seq_len(samples))),
                        row.names = colnames(ds)),
      images = list(coordinates = co,
                    spatial.factors = list(ratio = rep(1, samples), tol = rep(0.5, samples))),
      LR = list(LRsig = RNA_LRsig), net = list(),
      options = list(datatype = "spatial", mode = "single", db = normalizePath(DBDIR),
                     population.size = FALSE))
}

sp_cases <- list(
  list(name = "distance_use", args = list(distance.use = TRUE, scale.distance = 10,
                                          k.min = 1, contact.range = 60,
                                          interaction.range = 100)),
  list(name = "no_distance_use", args = list(distance.use = FALSE, scale.distance = 10,
                                              k.min = 1, contact.range = 60,
                                              interaction.range = 100)),
  ## `d.min < 1` after scaling: upstream stops, and the message embeds `1/d.min` at 2 digits.
  list(name = "scale_too_small", args = list(distance.use = TRUE, scale.distance = 0.001,
                                             k.min = 1, contact.range = 60,
                                             interaction.range = 100)),
  ## `contact.dependent.forced` takes the first branch and multiplies by `adj.contact` for all.
  list(name = "forced_contact", args = list(distance.use = TRUE, scale.distance = 10,
                                            k.min = 1, contact.range = 60,
                                            interaction.range = 100,
                                            contact.dependent.forced = TRUE)),
  ## A contact-only LR set: the second `nLR1` branch, and a different `cat`.
  list(name = "all_contact", lr = "Cell-Cell Contact",
       args = list(distance.use = TRUE, scale.distance = 10, k.min = 1,
                   contact.range = 60, interaction.range = 100)),
  ## Both kinds: the fourth branch, `nLR1 <- max(which(annotation %in% diffusible))`.
  list(name = "mixed", lr = c("Cell-Cell Contact", "Secreted Signaling"),
       args = list(distance.use = TRUE, scale.distance = 10, k.min = 1,
                   contact.range = 60, interaction.range = 100)),
  ## `contact.dependent = FALSE` reaches the final `else`, whose `cat` is a fourth distinct text.
  list(name = "contact_false", args = list(distance.use = TRUE, scale.distance = 10,
                                           k.min = 1, contact.range = 60,
                                           interaction.range = 100, contact.dependent = FALSE)),
  ## `k.min` above the number of distinct reachable cells: everything non-adjacent, `d.spatial`
  ## all `NaN`, so `min(na.rm = TRUE)` is `Inf` and the `d.min < 1` stop does *not* fire.
  list(name = "k_min_too_large", args = list(distance.use = TRUE, scale.distance = 10,
                                             k.min = 99, contact.range = 60,
                                             interaction.range = 100)),
  ## A single sample, so `ratio`/`tol` are length 1 and the per-sample axis collapses.
  list(name = "one_sample", samples = 1L,
       args = list(distance.use = TRUE, scale.distance = 10, k.min = 1,
                   contact.range = 60, interaction.range = 100)),
  ## `contact.knn.k` instead of `contact.range`, which zeroes `adj.contact` before the swap.
  list(name = "knn_k", args = list(distance.use = TRUE, scale.distance = 10, k.min = 1,
                                   contact.knn.k = 2L, interaction.range = 100))
)
for (cs in sp_cases) {
  o <- sp_obj(samples = cs$samples %||% 2L)
  if (!is.null(cs$lr)) {
    ## Vary the `annotation` column only: it is what selects the four `nLR1` branches, and the
    ## ligand/receptor columns have to stay as the fixture defines them or the two sides resolve
    ## different gene pairs (see `SP_DB` above).
    ## `LRsig` has one row per interaction, so a single annotation name has to be recycled --
    ## assigning a length-1 vector to an 8-row column *works* and assigning a length-2 vector
    ## does not, which makes the `mixed` case fail with a recycling error from `$<-.data.frame`.
    o@LR$LRsig$annotation <- rep(cs$lr, length.out = nrow(o@LR$LRsig))
  }
  call_args <- c(list(object = o, nboot = 3L, seed.use = 1L, type = "triMean"), cs$args)
  run_up <- function() {
    txt <- character(0)
    v <- withCallingHandlers(
      do.call(get("computeCommunProb", envir = sp_env), call_args),
      warning = function(w) invokeRestart("muffleWarning"),
      message = function(m) { txt <<- c(txt, conditionMessage(m)); invokeRestart("muffleMessage") })
    list(txt = txt, v = v)
  }
  run_rs <- function() {
    txt <- character(0)
    v <- withCallingHandlers(do.call(computeCommunProb, call_args),
      warning = function(w) invokeRestart("muffleWarning"),
      message = function(m) { txt <<- c(txt, conditionMessage(m)); invokeRestart("muffleMessage") })
    list(txt = txt, v = v)
  }
  ea <- tryCatch(run_up(), error = function(e) list(err = conditionMessage(e)))
  eb <- tryCatch(run_rs(), error = function(e) list(err = conditionMessage(e)))
  n_cmp <<- n_cmp + 3L
  if (!is.null(ea$err) || !is.null(eb$err)) {
    ok <- identical(ea$err, eb$err)
    cat(sprintf("sp_%-18s error=%-5s  up=%s | rs=%s\n", cs$name, ok,
                if (is.null(ea$err)) "<result>" else substr(gsub("[[:space:]]+", " ", ea$err), 1, 70),
                if (is.null(eb$err)) "<result>" else substr(gsub("[[:space:]]+", " ", eb$err), 1, 70)))
    if (!ok) fails <<- fails + 1L
    next
  }
  ok_p <- identical(ea$v@net$prob, eb$v@net$prob)
  ok_v <- identical(ea$v@net$pval, eb$v@net$pval)
  ## `run.time` is wall-clock and is deliberately excluded, as everywhere else.
  ok_par <- identical(ea$v@options$parameter, eb$v@options$parameter)
  ok_txt <- identical(paste(ea$txt, collapse = ""), paste(eb$txt, collapse = ""))
  ## The printed `Sys.time()` lines are compared with the timestamps removed, since they are the
  ## clock rather than the behaviour.
  strip_time <- function(x) gsub("\\[[0-9-]+ [0-9:.]+\\]", "[TIME]", paste(x, collapse = ""))
  ok_print <- identical(strip_time(ea$txt), strip_time(eb$txt))
  cat(sprintf("sp_%-18s prob=%-5s pval=%-5s parameter=%-5s msgs=%-5s (%d positive)\n",
              cs$name, ok_p, ok_v, ok_par, ok_print, sum(eb$v@net$prob > 0)))
  if (!ok_p) { fails <<- fails + 1L
    cat("   prob dim up:", paste(dim(ea$v@net$prob), collapse = "x"),
        " rs:", paste(dim(eb$v@net$prob), collapse = "x"), "\n")
    cat("   sum up:", sum(ea$v@net$prob), " rs:", sum(eb$v@net$prob), "\n")
    ## `arr.ind = TRUE` indexes as `[i, , j]` for an array, so a linear difference index has to be
    ## converted rather than dropped. `head(d)` on an empty selection also prints noise, so the
    ## count gates the report.
    pu <- as.vector(ea$v@net$prob); pr <- as.vector(eb$v@net$prob)
    bad <- which(!(isTRUE(all.equal(pu, pr)) | (pu == pr | (is.na(pu) & is.na(pr)))))
    cat("   differing entries:", length(bad), "of", length(pu), "\n")
    dn <- dimnames(ea$v@net$prob)
    for (idx in head(bad, 8)) {
      ij <- arrayInd(idx, dim(ea$v@net$prob))
      cat(sprintf("     [%s,%s | %s]  up=%s  rs=%s\n",
                  dn[[1]][ij[1]], dn[[2]][ij[2]], dn[[3]][ij[3]],
                  format(pu[idx], digits = 10), format(pr[idx], digits = 10)))
    }
    cat("   ratio rs/up (nonzero):",
        paste(format(pr[pr != 0 & pu != 0] / pu[pr != 0 & pu != 0], digits = 6), collapse = " "), "\n")
  }
  if (!ok_v) { fails <<- fails + 1L; cat("   pval differs\n") }
  if (!ok_par) {
    fails <<- fails + 1L
    cat("   parameter objects:\n")
    for (f in union(names(ea$v@options$parameter), names(eb$v@options$parameter)))
      if (!identical(ea$v@options$parameter[[f]], eb$v@options$parameter[[f]]))
        cat("   parameter", f, ": up =", paste(format(ea$v@options$parameter[[f]]), collapse = ","),
            " rs =", paste(format(eb$v@options$parameter[[f]]), collapse = ","), "\n")
  }
  if (!ok_print) {
    fails <<- fails + 1L
    cat("   up:", strip_time(ea$txt), "\n   rs:", strip_time(eb$txt), "\n")
  }
}

## One class for both rankNet fixtures. Declared before the comparison block because the
## single-mode block below already had it further down, and a class used before its
## `setClass` is a bare "getClass: object not found" that names neither block.
setClass("MiniRank", representation(netP = "list", net = "list", LR = "list",
                                    options = "list", DB = "list"))

## ----------------------------------------------------- rankNet, comparison mode
##
## `ggplot2` is attached for this block, and not incidentally. `cellchatrs_upstream_env()` sources
## upstream with `parent = globalenv()`, so upstream's body resolves its *unqualified* names --
## `ggplot`, `theme`, `element_text`, and `CellChat_theme_opts()`'s internals -- from the caller's
## search path. A real CellChat session has `ggplot2` attached because CellChat imports it. This
## gate does not, and without the attach the *reference* fails with "could not find function
## \"ggplot\"" while the port, which calls `ggplot2::` explicitly, succeeds: the failure looks
## like a port bug and is a fixture-environment bug.
if (!requireNamespace("ggplot2", quietly = TRUE)) {
  stop("this block needs ggplot2 attached: upstream's body resolves it unqualified", call. = FALSE)
}
suppressPackageStartupMessages(library(ggplot2))
##
## Driven from the `cmp_*` records of `ranknet_golden.txt`, so the fixture is the same object the
## Rust suite checks and a divergence here is a divergence there too.
##
## Two comparisons of what is compared:
##
##   * `signaling.contribution` with `identical()`. This is the numeric contract: the union of the
##     pathway vocabularies, the row order from the multi-key `order()`, the factor levels, the
##     `rbind`-uniquified row names, and the one-significant-digit relative contributions.
##   * `ggplot_build(gg)$data` rather than `identical(gg.obj)`. `identical()` on the ggplot object
##     **cannot** hold: `geom_bar`'s `position` carries a quosure -- an expression plus the
##     environment it was captured in -- and upstream's plot code is sourced into
##     `cellchatrs_upstream_env()` while the shim's runs in the package namespace. What can hold,
##     and does, is the rendered geometry. (PLAN.md 14.1 keeps visualization in R precisely because
##     it is not the numeric surface; this is the concrete reason a plot object is compared by its
##     build rather than by identity.)
rn_cmp_root <- Sys.getenv("CELLCHATRS_ROOT", unset = getwd())
rn_cmp_lines <- readLines(file.path(rn_cmp_root, "tests", "fixtures", "ranknet_golden.txt"),
                          warn = FALSE)
## A local field extractor. The region block below defines `rg_field` for the same job, but it
## comes *later* in this file, and a helper used before its definition is a "object not found"
## that names the helper rather than the block.
rn_field <- function(lines, i) vapply(strsplit(lines, "\t", fixed = TRUE), function(f) f[i], "")
rn_cmp_recs <- function(tag) rn_cmp_lines[rn_field(rn_cmp_lines, 1) == tag]
## Same "filter before rbind" rule as the region corpus: `rn_field` gives the second field, and
## records for different cases have different widths.
rn_cmp_rows <- function(tag, nm) {
  f <- rn_cmp_recs(tag)
  if (!length(f)) return(character(0))
  sp <- strsplit(f, "\t", fixed = TRUE)
  sp <- sp[vapply(sp, function(g) length(g) > 2L && identical(g[2], nm), logical(1))]
  if (!length(sp)) return(character(0))
  do.call(rbind, sp)
}
## The payload only: field 1 is the record tag and field 2 the case name. Returning `r[1, -1]`
## would hand back the case name as the first value, which for a three-field record is half the
## row and shifts every subsequent index by one -- an error that shows up as a `dimnames` length
## complaint several lines later.
rn_cmp_one <- function(tag, nm) {
  r <- rn_cmp_rows(tag, nm)
  if (!length(r)) return(NULL)
  r[1, -(1:2)]
}
rn_cmp_vec <- function(tag, nm) {
  r <- rn_cmp_rows(tag, nm)
  if (!length(r)) return(numeric(0))
  as.numeric(r[1, -(1:2)])
}
rn_cmp_names_of <- function(f) {
  if (is.na(f)) return(NULL)
  f
}
rn_cmp_case_names <- rn_field(rn_cmp_recs("cmp_case"), 2)
for (nm in rn_cmp_case_names) {
  meta <- rn_cmp_one("cmp_meta", nm)
  stopifnot(length(meta) >= 4L)
  ncomp <- as.integer(meta[3])
  levs <- strsplit(rn_cmp_one("cmp_levels", nm), ",", fixed = TRUE)[[1]]
  names.list <- lapply(strsplit(meta[4], ";", fixed = TRUE)[[1]],
                       function(g) strsplit(g, ",", fixed = TRUE)[[1]])
  prob.list <- pval.list <- vector("list", ncomp)
  ## The `cmp_prob` records are flat, column-major, one per comparison: `k * k * n` values for
  ## that comparison's own pathway count. Rebuilding the array is what puts the *names* back --
  ## `apply(prob, 3, sum)` returns a named vector and every downstream name-based step depends on
  ## it, so a fixture whose dimnames are absent would compare a different question.
  for (i in seq_len(ncomp)) {
    tag <- sprintf("%s#%d", nm, i)
    np <- length(names.list[[i]])
    for (what in c("cmp_prob", "cmp_pval")) {
      r <- rn_cmp_rows(what, tag)
      stopifnot(nrow(r) == 1L, ncol(r) == 2L + 16L * np)
    }
    r <- rn_cmp_rows("cmp_prob", tag)
    v <- as.numeric(r[1, -(1:3)])
    prob.list[[i]] <- array(v, dim = c(4, 4, np), dimnames = list(levs, levs, names.list[[i]]))
    r <- rn_cmp_rows("cmp_pval", tag)
    v <- as.numeric(r[1, -(1:3)])
    pval.list[[i]] <- array(v, dim = c(4, 4, np), dimnames = list(levs, levs, names.list[[i]]))
  }
  o <- new("MiniRank")
  o@netP <- lapply(seq_len(ncomp), function(i)
    list(group = levs[1], prob = prob.list[[i]], pval = pval.list[[i]]))
  names(o@netP) <- paste0("d", seq_len(ncomp))
  o@net <- list(prob = prob.list[[1]], pval = pval.list[[1]])
  o@LR <- list(LRsig = data.frame())
  o@options <- list(mode = "merged")

  ## The corpus writes the explicit marker `-` for an absent argument, not an empty field --
  ## `sprintf` recycles its format against its arguments, so an omitted value silently shifts the
  ## rest of the record. `!is.na()` is not the test: `"-"` is a string, and passing it as
  ## `sources.use` is a *character* argument that fails the group-name check.
  absent <- function(x) is.null(x) || length(x) == 0L || identical(x, "-")
  srcs <- rn_cmp_one("cmp_sources", nm)
  tgts <- rn_cmp_one("cmp_targets", nm)
  args <- list(object = o, slot.name = "netP", measure = meta[1], mode = "comparison",
               comparison = seq_len(ncomp), thresh = as.numeric(meta[2]),
               color.use = c("#4C72B0", "#DD8452", "#55A868", "#C44E52"),
               return.data = TRUE)
  if (!absent(srcs)) args$sources.use <- srcs
  if (!absent(tgts)) args$targets.use <- tgts

  ## Both sides go through a wrapper that puts the return value under `v`, so the error case and
  ## the value case share one shape. Calling `rankNet` directly and then reading `ea$v` leaves
  ## `ea$v` NULL, and `sprintf("%d", NULL)` returns `character(0)` -- so `cat(sprintf(...))` prints
  ## **nothing at all** and six passing cases are simply absent from the report with the summary
  ## still reading "IDENTICAL". A silent hole in a gate is worse than a failing one.
  wrap <- function(f) tryCatch(
    suppressWarnings(suppressMessages(list(v = do.call(f, args)))),
    error = function(e) list(err = conditionMessage(e)))
  ea <- wrap(cellchatrs_upstream_rankNet)
  eb <- wrap(rankNet)
  ## Two oracles, and this block is the one that decides between them.
  ##
  ## The corpus records the error the **harness** -- upstream's `mode = "comparison"` body lifted
  ## out of the function -- produces. The Rust suite asserts that record verbatim, and it is the
  ## right oracle for the *kernel*. This block's oracle is **upstream's exported `rankNet`**, and
  ## that is the contract the shim has to meet.
  ##
  ## For two fixtures the two disagree: the harness stops and `rankNet` does not. `rankNet` wraps
  ## its body in `options(warn = -1)` and is reached through `mergeCellChat`'s object, so an
  ## intermediate it computes differently no longer propagates. Asserting the harness record here
  ## would therefore report two failures that are not port bugs, and *skipping* the case would
  ## hide a real regression. So: when either side errors, the two must produce the **same** error;
  ## when neither does, the case is compared normally and the unreproduced harness record is
  ## printed rather than silently dropped.
  want_err <- rn_cmp_one("cmp_error", nm)
  if (!is.null(ea$err) || !is.null(eb$err)) {
    ok <- identical(ea$err, eb$err) &&
      (is.null(want_err) || identical(ea$err, want_err))
    cat(sprintf("rc_%-22s error=%-5s  up=%s | rs=%s\n", nm, ok,
                if (is.null(ea$err)) "<result>" else substr(gsub("[[:space:]]+", " ", ea$err), 1, 44),
                if (is.null(eb$err)) "<result>" else substr(gsub("[[:space:]]+", " ", eb$err), 1, 44)))
    if (!ok) fails <<- fails + 1L
    next
  }
  if (!is.null(want_err)) {
    cat(sprintf("rc_%-22s harness-error-not-reproduced  upstream returns a value; \
the corpus records %s for the lifted body\n", nm, substr(want_err, 1, 44)))
  }
  da <- ea$v$signaling.contribution
  db <- eb$v$signaling.contribution
  ok_v <- identical(da, db)
  ## The rendered geometry, which is the strongest claim available for the plot.
  ## `ggplot_build()` is an S3 generic in ggplot2 4.x, so it dispatches on the *package* rather
  ## than being called directly -- `ggplot2::ggplot_build(p)` works, but `UseMethod("ggplot_build")`
  ## failing means the object is not a ggplot at all, so the call is guarded and reported rather
  ## than allowed to abort the whole gate.
  built <- function(p) tryCatch(suppressWarnings(ggplot2::ggplot_build(p)$data),
                                error = function(e) paste("build failed:", conditionMessage(e)))
  ok_gg <- identical(built(ea$v$gg.obj), built(eb$v$gg.obj))
  cat(sprintf("rc_%-22s value=%-5s built=%-5s (%d x %d, rows %s)\n", nm, ok_v, ok_gg,
              nrow(db), ncol(db), paste(head(rownames(db), 4), collapse = ",")))
  if (!ok_v) {
    fails <<- fails + 1L
    cat("   upstream:\n"); print(da)
    cat("   port    :\n"); print(db)
    for (cl in intersect(names(da), names(db))) {
      if (!identical(da[[cl]], db[[cl]]))
        cat("   column", cl, "up=", paste(format(head(da[[cl]], 6)), collapse = ","),
            " rs=", paste(format(head(db[[cl]], 6)), collapse = ","), "\n")
    }
  }
  if (!ok_gg) {
    fails <<- fails + 1L
    ba <- list(data = built(ea$v$gg.obj))
    bb <- list(data = built(eb$v$gg.obj))
    for (i in seq_along(ba$data)) {
      if (!identical(ba$data[[i]], bb$data[[i]]))
        cat("   built layer", i, "differs:",
            paste(utils::head(colnames(ba$data[[i]]), 5), collapse = ","), "\n")
    }
  }
}

## ------------------------------------------------------------------- rankNet, single mode
## The corpus is the generator's own dump. `rankNet`'s plot is a ggplot, which is not the
## numeric surface, so the gate compares `signaling.contribution` with `identical()` and the
## plot's `$data` -- the frame ggplot actually built the bars from -- rather than the ggplot
## object itself, whose environment pointer is not a meaningful thing to compare.
rn_lines <- readLines("tests/fixtures/ranknet_golden.txt")
rn_rec <- function(tag, name) {
  hit <- grep(sprintf("^%s\t%s\t", tag, name), rn_lines, value = TRUE)
  if (!length(hit)) return(NA_character_)
  sub(sprintf("^%s\t%s\t", tag, name), "", hit[length(hit)])
}
rn_row <- function(tag, name) { strsplit(rn_rec(tag, name), "\t")[[1]] }
rn_cell <- function(s) if (identical(s, "<NA>")) NA else if (identical(s, "<NaN>")) NaN else s
rn_vec <- function(tag, name) vapply(rn_row(tag, name), rn_cell, "")

rn_args <- list(
  weight_degenerate = list(), weight_three_degenerate = list(), weight_all_degenerate = list(),
  thresh_cuts = list(), count_measure = list(), sources_filter = list(sources.use = "g1"),
  targets_filter = list(targets.use = "g1"), both_filters = list(sources.use = "g1", targets.use = "g2"),
  all_zero = list(), thresh_empties = list(), order_ties = list()
)
for (nm in names(rn_args)) {
  meta <- rn_row("meta", nm)
  measure <- meta[1]; thresh <- as.numeric(meta[2])
  ## `meta[3]` is the *pathway* axis (`dimnames(prob)[[3]]`); the two group axes are `k` long and
  ## live in their own record. Reading the pathway names as the group levels works for a square
  ## k == n fixture and fails for every one where they differ, which is the worst kind of
  ## confusion: it passes on the first case and misreports the rest.
  pnames <- strsplit(meta[3], ",")[[1]]
  levs <- strsplit(rn_rec("levels", nm), ",")[[1]]
  ## `rn_rec` strips `meta<TAB><name><TAB>`, so the remaining fields are
  ## measure, thresh, names, k, n -- five of them, not six.
  k <- as.integer(meta[4]); n <- as.integer(meta[5])
  ## `as.numeric`, not the raw corpus text: a character array makes *upstream* fail with
  ## "invalid 'type' (character) of argument" on `prob[pval > thresh] <- 0`, which reads like a
  ## port bug and is really a fixture-construction bug. `as.numeric("NA")` is `NA_real_` and
  ## `as.numeric("NaN")` is `NaN`, so the two markers stay apart.
  pv <- array(as.numeric(rn_vec("prob", nm)), dim = c(k, k, n),
              dimnames = list(levs, levs, pnames))
  qv <- array(as.numeric(rn_vec("pval", nm)), dim = c(k, k, n))
  net <- list(prob = pv, pval = qv)
  slotname <- "netP"
  o <- new("MiniRank"); o@netP <- list(group = "G1", prob = pv, pval = qv); o@net <- net
  a <- rn_args[[nm]]
  run_up <- function() cellchatrs_upstream_rankNet(
    object = o, slot.name = slotname, measure = measure, mode = "single",
    comparison = 1, thresh = thresh, return.data = TRUE)
  run_rs <- function() rankNet(
    object = o, slot.name = slotname, measure = measure, mode = "single",
    comparison = 1, thresh = thresh, return.data = TRUE)
  ea <- tryCatch(suppressWarnings(suppressMessages(run_up())),
                 error = function(e) list(err = conditionMessage(e)))
  eb <- tryCatch(suppressWarnings(suppressMessages(run_rs())),
                 error = function(e) list(err = conditionMessage(e)))
  n_cmp <<- n_cmp + 2L
  want_err <- !is.na(rn_rec("error", nm))
  if (want_err) {
    ok <- identical(ea$err, eb$err)
    cat(sprintf("rn_%-24s error=%-5s  %s\n", nm, ok,
                if (is.null(eb$err)) "<result>" else substr(gsub("[[:space:]]+", " ", eb$err), 1, 40)))
    if (!ok) { fails <<- fails + 1L
      cat("   up:", if (is.null(ea$err)) "<result>" else ea$err,
          "\n   rs:", if (is.null(eb$err)) "<result>" else eb$err, "\n") }
    next
  }
  ok_df <- identical(ea$v$signaling.contribution, eb$v$signaling.contribution)
  ## The plot's data, not the plot: a ggplot carries an environment whose pointer is not a
  ## reproducible thing to compare, but the frame it was built from is the numerics.
  ok_gg <- identical(ea$v$gg.obj$data, eb$v$gg.obj$data)
  cat(sprintf("rn_%-24s data=%-5s plot.data=%-5s (%d x %d)\n", nm, ok_df, ok_gg,
              nrow(eb$v$signaling.contribution), ncol(eb$v$signaling.contribution)))
  if (!ok_df) {
    fails <<- fails + 1L
    a1 <- ea$v$signaling.contribution; b1 <- eb$v$signaling.contribution
    cat("   dim:", paste(dim(a1), collapse = "x"), "vs", paste(dim(b1), collapse = "x"), "\n")
    cat("   rn :", paste(rownames(a1), collapse = ","), "|",
        paste(rownames(b1), collapse = ","), "\n")
    for (cl in intersect(names(a1), names(b1))) {
      if (!identical(a1[[cl]], b1[[cl]])) {
        cat("   column", cl, "differs\n   up:"); print(utils::head(a1[[cl]], 5))
        cat("   rs:"); print(utils::head(b1[[cl]], 5)); break
      }
    }
    d <- setdiff(names(a1), names(b1)); if (length(d)) cat("   only upstream:", d, "\n")
    d <- setdiff(names(b1), names(a1)); if (length(d)) cat("   only port   :", d, "\n")
  }
  if (!ok_gg) fails <<- fails + 1L
}

## -------------------------------------------------------- rankNetPairwise, against upstream
##
## `rankNetPairwise` is one `order(pval, -prob)` per `(i, j)` group pair plus a great deal of
## data-frame construction. The ordering is Rust's; the `data.frame()`, its `row.names`, the nested
## `list()` and the `names(temp) <- colnames(prob)` are R's, and the gate compares the whole
## returned object so both halves are covered.
##
## The cases are chosen for the parts that are easy to get wrong and invisible in a "does it work"
## test:
##
##   * **ties.** `order` is a stable radix sort, so equal `(pval, -prob)` keeps input order. A
##     non-stable sort produces a different *row order* and an otherwise perfectly plausible result.
##   * **a missing value in the second key.** `order(c(1,1,2), c(NA,5,3))` is `2 1 3` in R: the
##     `NA` sorts last *within* its group of equal first keys. The first version of the Rust
##     ordering compared a missing value as equal to a finite one, which is `1 2 3` -- a different
##     answer that only shows up when a `pval` repeats and a `prob` is `NaN`.
##   * **non-square `prob`.** `numCluster <- dim(prob)[1]` and upstream loops `1:numCluster` in both
##     dimensions, so it reads the leading `k x k` block; then `names(temp) <- colnames(prob)` has
##     the wrong length and raises. The port must reach the same error from its own code, not by
##     noticing and handing over.
##   * **`LR.use` shorter than `dim(prob)[3]`.** Upstream does not reconcile the two: its
##     `data.frame(row.names = rownames(pairLR.use))` raises, because `probij` has `dim(prob)[3]`
##     values. Slicing `prob` to `LR.use` -- the obvious reading -- turns that into a success.
##
## The shim warns when it falls back to upstream, and the block counts those warnings: a
## pass-through and a port are both "identical" to upstream, so a gate that cannot tell them apart
## reports success for a function that is not ported at all. That is not hypothetical -- it is what
## hid the kernel reading R's `dim()` as a list of slice lengths.
pw_fell <- 0L
pw_cmp <- function(lbl, o, lr.use = NULL) {
  up <- tryCatch(cellchatrs_upstream_rankNetPairwise(o, lr.use),
                 error = function(e) structure(conditionMessage(e), class = "pw.err"))
  po <- withCallingHandlers(
    tryCatch(rankNetPairwise(o, lr.use), error = function(e) structure(conditionMessage(e), class = "pw.err")),
    warning = function(w) {
      if (grepl("^cellchatrs: rankNetPairwise answered from upstream", conditionMessage(w))) {
        pw_fell <<- pw_fell + 1L
      }
      invokeRestart("muffleWarning")
    })
  lbl0 <- sprintf("rankNetPairwise/%s", lbl)
  if (inherits(up, "pw.err") && inherits(po, "pw.err")) {
    ## Both raise: the *message* is the contract, and it is upstream's own error text.
    n_cmp <<- n_cmp + 1L
    ok <- identical(as.character(up), as.character(po))
    cat(sprintf("%-46s error=%-5s  up: %s\n", lbl0, ok, up))
    if (!ok) {
      fails <<- fails + 1L
      cat(sprintf("%-46s       port: %s\n", "", po))
    }
    return(invisible(NULL))
  } else if (inherits(up, "pw.err") || inherits(po, "pw.err")) {
    n_cmp <<- n_cmp + 1L
    fails <<- fails + 1L
    cat(sprintf("%-46s ONE-SIDED-ERROR  up=%s | port=%s\n", lbl0,
                if (inherits(up, "pw.err")) up else "an object",
                if (inherits(po, "pw.err")) po else "an object"))
  } else if (identical(up, po)) {
    n_cmp <<- n_cmp + 1L
    cat(sprintf("%-46s pairwiseRank=%s\n", lbl0, TRUE))
  } else {
    d <- character(0)
    for (i in seq_along(up@net$pairwiseRank)) for (j in seq_along(up@net$pairwiseRank[[i]])) {
      a <- up@net$pairwiseRank[[i]][[j]]; b <- po@net$pairwiseRank[[i]][[j]]
      if (!identical(a, b)) {
        d <- c(d, sprintf("[%d][%d] upstream=%s port=%s", i, j,
                          paste(rownames(a), collapse = ","), paste(rownames(b), collapse = ",")))
      }
    }
    n_cmp <<- n_cmp + 1L
    fails <<- fails + 1L
    cat(sprintf("%-46s DIFFERS\n%s\n", lbl0,
                paste0("      ", d, collapse = "\n")))
  }
}
pw_obj <- function(d, prob, pval, lr.rows = d[3], dimnames = TRUE) {
  ## `array()` is a no-op on an array of the same shape, so this accepts both the arrays the
  ## cases mutate with `[i, j, l]` and the plain vectors they pass for constant cases. Without it
  ## those cases fail with `'dimnames' applied to non-array` on *both* sides, which the block
  ## scores as a pass -- three "errors matched" that were really fixture bugs.
  prob <- array(prob, d); pval <- array(pval, d)
  if (dimnames) dimnames(prob) <- list(paste0("g", seq_len(d[1])), paste0("g", seq_len(d[2])),
                                       paste0("P", seq_len(d[3])))
  o <- new("MiniRank")
  o@net <- list(prob = prob, pval = pval)
  o@LR <- list(LRsig = data.frame(
    interaction_name = paste0("I", seq_len(lr.rows)), interaction_name_2 = paste0("J", seq_len(lr.rows)),
    pathway_name = paste0("P", seq_len(lr.rows)), ligand = paste0("L", seq_len(lr.rows)),
    receptor = paste0("R", seq_len(lr.rows)), row.names = paste0("row", seq_len(lr.rows))))
  o
}
## `pw_arr` returns a **filled array**, not a vector: the cases below mutate it with `[i, j, l]`
## indexing, and mutating the flat vector is a subscript error that stops the gate dead rather
## than failing a case.
pw_arr <- function(d, f = NULL) {
  a <- array(if (is.null(f)) round(runif(prod(d)), 2) else f, d)
  a
}
pw_r <- function(d) array(round(runif(prod(d)), 2), d)
set.seed(1L)
pw_cmp("random",                pw_obj(c(3, 3, 6), pw_r(c(3, 3, 6)), pw_r(c(3, 3, 6))))
pw_cmp("all ties",              pw_obj(c(3, 3, 4), rep(0.5, 36), rep(0.5, 36)))
p <- pw_r(c(3, 3, 6)); p[c(1, 1, 1, 2, 2, 3)] <- 0.3
pw_cmp("prob ties, varying pval", pw_obj(c(3, 3, 6), p, pw_r(c(3, 3, 6))))
v <- pw_r(c(3, 3, 6)); v[1, 1, 1] <- NA_real_
pw_cmp("one NA pval",           pw_obj(c(3, 3, 6), pw_r(c(3, 3, 6)), v))
b1 <- pw_r(c(3, 3, 6)); b1[1, 1, 2] <- NaN
pw_cmp("one NaN prob",          pw_obj(c(3, 3, 6), b1, pw_r(c(3, 3, 6))))
v2 <- pw_r(c(3, 3, 6)); v2[1, 1, 1] <- NaN; b2 <- pw_r(c(3, 3, 6)); b2[2, 2, 3] <- Inf
pw_cmp("NaN pval + Inf prob",   pw_obj(c(3, 3, 6), b2, v2))
b3 <- pw_r(c(3, 3, 6)); b3[1, 1, 1] <- NaN
## The exact `order` case that caught the second-key `NA` bug: equal pval, one `NaN` prob.
pw_cmp("tied pval + NaN prob",  pw_obj(c(3, 3, 6), b3, rep(0.49, 54)))
pw_cmp("k=1 n=3",               pw_obj(c(1, 1, 3), c(0.2, 0.5, 0.1), c(0.3, 0.1, 0.2)))
o7 <- pw_obj(c(3, 3, 6), pw_r(c(3, 3, 6)), pw_r(c(3, 3, 6)))
pw_cmp("LR.use shorter than n", o7, o7@LR$LRsig[1:2, ])
o8 <- pw_obj(c(3, 3, 6), pw_r(c(3, 3, 6)), pw_r(c(3, 3, 6)))
pw_cmp("LR.use reversed",       o8, o8@LR$LRsig[6:1, ])
pw_cmp("no dimnames",           pw_obj(c(3, 3, 4), pw_r(c(3, 3, 4)), pw_r(c(3, 3, 4)), dimnames = FALSE))
o10 <- pw_obj(c(3, 3, 6), pw_r(c(3, 3, 6)), pw_r(c(3, 3, 6))); o10@LR$LRsig <- o10@LR$LRsig[0, ]
pw_cmp("empty LRsig",           o10)
pw_cmp("k=4 n=2",               pw_obj(c(4, 4, 2), pw_r(c(4, 4, 2)), pw_r(c(4, 4, 2))))
pw_cmp("k=2 n=3",               pw_obj(c(2, 2, 3), pw_r(c(2, 2, 3)), pw_r(c(2, 2, 3))))
pw_cmp("non-square 2x5x3",      pw_obj(c(2, 5, 3), pw_r(c(2, 5, 3)), pw_r(c(2, 5, 3))))
pw_cmp("k=1 n=1",               pw_obj(c(1, 1, 1), 0.5, 0.5))
pw_cmp("5x3x4",                 pw_obj(c(5, 3, 4), pw_r(c(5, 3, 4)), pw_r(c(5, 3, 4))))
pw_cmp("3x1x2",                 pw_obj(c(3, 1, 2), pw_r(c(3, 1, 2)), pw_r(c(3, 1, 2))))
if (pw_fell > 0L) {
  fails <<- fails + pw_fell
  cat(sprintf("MISMATCH: rankNetPairwise fell back to upstream on %d of 18 cases\n", pw_fell))
}

## ------------------------------------------------------------- filterCommunication
## Both halves, both `type.mean` choices that change the per-sample means, `rare.keep` both
## ways, and both upstream failures. The `cat()` output is compared too: it contains
## `scales::percent()` and a literal tab, and it is the only place the *counts* surface.
##
## The fixture definitions come from the generator's source with its output connection
## removed. `eval`-ing the preamble unmodified truncates the committed corpus, because the
## preamble opens the golden file with "wt" -- that has happened here twice, so the generator
## now writes to a `.part` file and renames.
fc_src <- readLines("tests/parity/gen_filter_golden.R")
fc_src <- fc_src[!grepl("^\\s*(q|qi) <- file", fc_src)]
fc_src <- fc_src[!grepl("^\\s*(cat|writeLines|close|file.rename)\\(", fc_src)]
fc_stop <- grep("^## -{10,} run and record", fc_src)[1]
eval(parse(text = paste(fc_src[seq_len(fc_stop - 1)], collapse = "\n")))

setClass("MiniFilter", representation(
  net = "list", idents = "factor", meta = "list", options = "list",
  data.signaling = "ANY", DB = "list", data.smooth = "ANY"
))

for (cs in cases) {
  fc_data <- mk_data(cs$genes, cs$ncell, seed = 20240440 + length(cases))
  fc_tm <- if (is.null(cs$type_mean)) "triMean" else cs$type_mean
  ## `raw.use` and `data.smooth`, read from the corpus rather than assumed. The generator resolves
  ## the smooth matrices once, up front, for *both* of its dump loops -- a value built inside one
  ## loop is invisible to the other, and the first version of this block read a `data.smooth`
  ## that only two of the cases have.
  fc_ru <- isTRUE(cs$raw_use)
  fc_smooth <- cs$smooth
  fc_drop <- isTRUE(cs$drop_smooth)
  fcmk <- function() {
    o <- mk(cs$prob, cs$pval, cs$lev, cs$lr, cs$ident, cs$samples, fc_data, DB,
            type_mean = fc_tm, raw_use = fc_ru, with_smooth = !fc_ru && !fc_drop,
            smooth = fc_smooth)
    o
  }
  a <- fcmk()
  b <- fcmk()
  ## Keep the *returned* object. S4 objects are copied, so discarding the return value would
  ## compare the two unfiltered inputs and every case would pass trivially -- which is exactly
  ## what the first version of this block did.
  ea <- tryCatch({
    txt <- capture.output(oa <- cellchatrs_upstream_filterCommunication(
      a, min.cells = cs$min_cells, min.samples = cs$min_samples,
      rare.keep = cs$rare_keep, nonFilter.keep = cs$nonfilter_keep))
    list(msg = txt, obj = oa)
  }, error = function(e) list(err = conditionMessage(e)))
  eb <- tryCatch({
    txt <- capture.output(ob <- filterCommunication(
      b, min.cells = cs$min_cells, min.samples = cs$min_samples,
      rare.keep = cs$rare_keep, nonFilter.keep = cs$nonfilter_keep))
    list(msg = txt, obj = ob)
  }, error = function(e) list(err = conditionMessage(e)))
  n_cmp <<- n_cmp + 2L
  if (!is.null(ea$err) || !is.null(eb$err)) {
    ok <- identical(ea$err, eb$err)
    cat(sprintf("fc_%-18s error=%-5s  up=%s | rs=%s\n", cs$name, ok,
                if (is.null(ea$err)) "<result>" else substr(ea$err, 1, 44),
                if (is.null(eb$err)) "<result>" else substr(eb$err, 1, 44)))
    if (!ok) fails <<- fails + 1L
    next
  }
  ## `capture.output` keeps a `cat`'s tab-and-newline pair as one line per call, so compare
  ## the concatenation: the tab is inside the first `cat` and the percentage inside the second.
  ok_obj <- identical(ea$obj, eb$obj) && identical(ea$obj@net, eb$obj@net)
  ok_msg <- identical(paste(ea$msg, collapse = ""), paste(eb$msg, collapse = ""))
  cat(sprintf("fc_%-18s object=%-5s msgs=%-5s (%d lines)\n", cs$name, ok_obj, ok_msg,
              length(eb$msg)))
  if (!ok_obj) {
    fails <<- fails + 1L
    cat("   upstream positives:", sum(ea$obj@net$prob > 0),
        " rust:", sum(eb$obj@net$prob > 0), "\n")
  }
  if (!ok_msg) {
    fails <<- fails + 1L
    cat("   upstream:", paste(ea$msg, collapse = "<NL>"), "\n")
    cat("   rust    :", paste(eb$msg, collapse = "<NL>"), "\n")
  }
}

## --------------------------------------------------- computeCellDistance, against upstream
##
## The only spatial function whose upstream form *can* run here: it is `collapse::fdist` plus a
## threshold, with no neighbour query, so upstream is a valid oracle and this is an ordinary
## differential test. `computeRegionDistance` gets the corpus-based gate below instead, because
## upstream's own cannot run without `BiocNeighbors`.
##
## What is under test is the *shim*, not `fdist`: the Rust `fdist` already has its own parity suite.
## The two-column error, the `colnames` assignment upstream makes before anything can fail, `ratio`
## being applied only when non-`NULL`, the message's leading newline, and `tol` being an **addend**
## rather than a tolerance are all R-side and none of them are visible from Rust.
cd_cases <- list(
  list(name = "plain",        coords = cbind(c(0, 3, 1), c(0, 4, 0))),
  list(name = "ratio",        coords = cbind(c(0, 3, 1), c(0, 4, 0)), ratio = 2),
  list(name = "zero_ratio",   coords = cbind(c(0, 3, 1), c(0, 4, 0)), ratio = 0),
  ## `ratio = numeric(0)` is *not* `NULL`: `if (!is.null(ratio))` is TRUE and `dm * numeric(0)` is
  ## an empty matrix, so the result keeps its dim and gains no values.
  list(name = "empty_ratio",  coords = cbind(c(0, 3, 1), c(0, 4, 0)), ratio = numeric(0)),
  ## `tol` alone does not fire the threshold: the guard is `!is.null(range) & !is.null(tol)`.
  list(name = "tol_only",     coords = cbind(c(0, 3, 1), c(0, 4, 0)), tol = 1),
  ## Strictly greater than, so a distance exactly on the threshold survives.
  list(name = "boundary",     coords = cbind(c(0, 5, 0), c(0, 0, 0)),
                             interaction.range = 4, tol = 1),
  list(name = "range_and_tol", coords = cbind(c(0, 3, 1), c(0, 4, 0)),
                             interaction.range = 4, tol = 1),
  ## The scaling happens *before* the threshold, so `ratio` moves what is on the boundary.
  list(name = "ratio_then_threshold", coords = cbind(c(0, 5, 0), c(0, 0, 0)),
                             interaction.range = 4, tol = 1, ratio = 2),
  list(name = "rownames",     coords = { m <- cbind(c(0, 3, 1), c(0, 4, 0))
                                        rownames(m) <- c("a", "b", "c"); m }),
  list(name = "duplicates",   coords = cbind(c(0, 0, 0, 1), c(0, 0, 0, 1))),
  list(name = "single_cell",  coords = cbind(0, 0)),
  list(name = "negatives",    coords = cbind(c(-5, 5, 0), c(-5, 5, 12))),
  ## The error: three columns, and upstream's own message, raised before the `colnames` line.
  list(name = "three_cols",   coords = cbind(1:3, 1:3, 1:3)),
  list(name = "one_col",      coords = matrix(1:3, ncol = 1)),
  list(name = "no_col",       coords = matrix(numeric(0), nrow = 0, ncol = 0))
)
for (cs in cd_cases) {
  run_cd <- function(f) {
    w <- character(0)
    v <- withCallingHandlers(
      suppressMessages(f(cs$coords, interaction.range = cs$interaction.range,
                        ratio = cs$ratio, tol = cs$tol)),
      message = function(x) { w <<- c(w, conditionMessage(x)); invokeRestart("muffleMessage") })
    list(w = w, v = v)
  }
  ea <- tryCatch(run_cd(cellchatrs_upstream_computeCellDistance), error = function(e) list(err = conditionMessage(e)))
  eb <- tryCatch(run_cd(computeCellDistance), error = function(e) list(err = conditionMessage(e)))
  n_cmp <<- n_cmp + 2L
  if (!is.null(ea$err) || !is.null(eb$err)) {
    ok <- identical(ea$err, eb$err)
    cat(sprintf("cd_%-20s error=%-5s  up=%s | rs=%s\n", cs$name, ok,
                if (is.null(ea$err)) "<result>" else substr(ea$err, 1, 44),
                if (is.null(eb$err)) "<result>" else substr(eb$err, 1, 44)))
    if (!ok) fails <<- fails + 1L
    next
  }
  ## The message carries upstream's leading newline, so the comparison has to keep it -- a
  ## `trimws` here would hide a missing or extra `\n`, which is exactly what a `message()` move
  ## between the shim and the kernel would introduce.
  ok_obj <- identical(ea$v, eb$v)
  ok_msg <- identical(ea$w, eb$w)
  cat(sprintf("cd_%-20s value=%-5s msgs=%-5s (%d x %d)\n", cs$name, ok_obj, ok_msg,
              nrow(eb$v), ncol(eb$v)))
  if (!ok_obj) {
    fails <<- fails + 1L
    cat("   upstream:\n"); print(ea$v)
    cat("   port    :\n"); print(eb$v)
  }
  if (!ok_msg) {
    fails <<- fails + 1L
    cat("   upstream msg:", paste(sprintf("[%s]", ea$w), collapse = ""), "\n")
    cat("   port     msg:", paste(sprintf("[%s]", eb$w), collapse = ""), "\n")
  }
}

## ----------------------------------------------- computeRegionDistance, against the corpus
##
## This corpus tests the optional exact-neighbour algorithm; it is not an upstream compatibility gate.
## Historical corpus provenance: upstream was not called here -- it evaluates `BiocNeighbors::queryKNN(...,
## AnnoyParam())` unconditionally and the package is not installable without network. So this is
## not a differential test; it replays `tests/fixtures/region_golden.txt` through the *shim*, whose
## job is the R-side marshalling the Rust suite cannot see: `levels(factor)` order, the
## declared-but-absent level surviving into the dimensions, `match()` on the factor, the dimnames,
## and recycling a scalar `ratio`/`tol`.
##
## The arithmetic and the neighbour query are covered by `src/rust/crates/r-core/tests/region_parity.rs`.
## A fixture passing in both places is the whole path.
rg_path <- file.path(Sys.getenv("CELLCHATRS_ROOT", unset = getwd()),
                     "tests", "fixtures", "region_golden.txt")
## Records are tab-separated and variable-length, so field access has to go through `strsplit`.
## `` `[1]` `` is *not* valid R for "element 1 of each" -- the parser reads the backtick as opening
## a name and `[1]` as its argument, giving a call to an undefined function -- so the extraction
## goes through an explicit closure. Getting this wrong fails at run time with an error that names
## nothing useful.
rg_field <- function(lines, i) {
  vapply(strsplit(lines, "\t", fixed = TRUE), function(f) f[i], "")
}
rg_recs <- function(lines, tag) lines[rg_field(lines, 1) == tag]
## Every record for one case, as a character matrix split into fields.
rg_rows <- function(lines, tag, nm) {
  f <- rg_recs(lines, tag)
  if (!length(f)) return(character(0))
  sp <- strsplit(f, "\t", fixed = TRUE)
  ## Filter to this case *before* `rbind`. Records for different cases of the same tag have
  ## different widths -- `group` has one field per cell, so a 15-cell case pads every other case's
  ## row out to 15 with `NA` -- and `rbind` on the whole tag then hands back a row whose width is
  ## the maximum across *all* cases. The values past this case's own length are `NA`, and the
  ## failure surfaces much later as a length mismatch in `data.frame`, with nothing pointing at the
  ## corpus reader.
  sp <- sp[vapply(sp, function(g) length(g) > 2L && identical(g[2], nm), logical(1))]
  if (!length(sp)) return(character(0))
  do.call(rbind, sp)
}
rg_num <- function(x) if (x %in% c("-", "", NA)) NA_real_ else as.numeric(x)

if (!file.exists(rg_path)) {
  cat("crd_ SKIPPED: region_golden.txt not found at", rg_path, "\n")
} else {
  rg_lines <- readLines(rg_path, warn = FALSE)
  rg_case_lines <- rg_recs(rg_lines, "case")
  for (nm in rg_field(rg_case_lines, 2)) {
    meta <- rg_rows(rg_lines, "meta", nm)
    opt  <- rg_rows(rg_lines, "opt", nm)
    flg  <- rg_rows(rg_lines, "flags", nm)
    ra   <- rg_rows(rg_lines, "ratio", nm)
    to   <- rg_rows(rg_lines, "tol", nm)
    grp  <- rg_rows(rg_lines, "group", nm)
    smp  <- rg_rows(rg_lines, "samples", nm)
    craw <- rg_rows(rg_lines, "coord_row", nm)
    err  <- rg_rows(rg_lines, "error", nm)
    stopifnot(nrow(meta) == 1L, nrow(opt) == 1L, nrow(flg) == 1L, nrow(ra) == 1L, nrow(to) == 1L)

    n <- as.integer(meta[1, 3]); d <- as.integer(meta[1, 4])
    levs <- strsplit(meta[1, 5], ",", fixed = TRUE)[[1]]
    nsamp <- as.integer(meta[1, 7])
    ## Coordinates arrive one row per record, so `nrow(craw)` is `n` and each record's fields 3.. are
    ## that row's `d` values. `vapply` over the rows rather than `do.call(rbind, ...)` so a zero-row
    ## fixture gives a `0 x d` matrix instead of `NULL`.
    coords <- matrix(0, nrow = nrow(craw), ncol = d)
    if (nrow(craw)) {
      for (r in seq_len(nrow(craw))) {
        coords[r, ] <- as.numeric(craw[r, -(1:3)])
      }
    }
    dimnames(coords) <- list(NULL, paste0("v", seq_len(d)))

    ## The corpus stores 0-based indices into its own level vectors, so the labels have to be put
    ## back before the factor is built. Deriving `levels()` from `unique()` instead would drop the
    ## declared-but-absent level and shrink the result -- the exact failure the fixture exists for.
    grp_lab <- levs[as.integer(grp[1, -(1:2)]) + 1L]
    smp_lev <- sprintf("s%d", seq_len(nsamp))
    smp_lab <- smp_lev[as.integer(smp[1, -(1:2)]) + 1L]
    meta_df <- data.frame(group = factor(grp_lab, levels = levs),
                          samples = factor(smp_lab, levels = smp_lev))

    args <- list(coordinates = coords, meta = meta_df,
                 interaction.range = rg_num(opt[1, 3]),
                 ratio = as.numeric(strsplit(ra[1, 3], ",", fixed = TRUE)[[1]]),
                 tol = as.numeric(strsplit(to[1, 3], ",", fixed = TRUE)[[1]]),
                 k.min = as.integer(meta[1, 8]),
                 contact.dependent = identical(meta[1, 9], "TRUE"),
                 contact.range = rg_num(opt[1, 4]),
                 contact.knn.k = if (is.na(rg_num(opt[1, 5]))) NULL else as.integer(opt[1, 5]),
                 do.symmetric = identical(flg[1, 3], "TRUE"))
    args$interaction.range <- if (is.na(args$interaction.range)) NULL else args$interaction.range
    args$contact.range <- if (is.na(args$contact.range)) NULL else args$contact.range

    eb <- tryCatch(do.call(cellchatrs_computeRegionDistance_exact, args), error = function(e) list(err = conditionMessage(e)))
    n_cmp <<- n_cmp + 2L
    if (length(err) >= 1L) {
      want <- err[1, 3]
      ok <- identical(if (is.null(eb$err)) "<result>" else eb$err, want)
      cat(sprintf("crd_%-28s error=%-5s  %s\n", nm, ok,
                  if (is.null(eb$err)) "<result>" else substr(eb$err, 1, 44)))
      if (!ok) { fails <<- fails + 1L
        cat("   want:", want, "\n   got :", if (is.null(eb$err)) "<result>" else eb$err, "\n") }
      next
    }
    if (!is.null(eb$err)) {
      fails <<- fails + 1L
      cat(sprintf("crd_%-28s value=%-5s  unexpected error: %s\n", nm, FALSE, eb$err))
      next
    }
    ## The corpus writes matrices one row per record, with `NaN` for the non-adjacent entries.
    read_mat <- function(tag) {
      rr <- rg_rows(rg_lines, tag, nm)
      if (!nrow(rr)) return(NULL)
      k <- length(levs)
      out <- matrix(0, nrow = k, ncol = k)
      for (r in seq_len(nrow(rr))) {
        out[as.integer(rr[r, 3]), ] <- as.numeric(rr[r, -(1:3)])
      }
      out
    }
    want_d <- read_mat("ds_row")
    want_a <- read_mat("ac_row")
    ## `all.equal(tolerance = 0)` on the un-named values: `NaN` compares equal to `NaN` there, and
    ## the dimnames are checked separately because the corpus does not record them.
    same <- function(x, y) {
      identical(dim(x), dim(y)) && isTRUE(all.equal(unname(x), unname(y), tolerance = 0))
    }
    ok_d <- same(eb$d.spatial, want_d)
    ok_a <- same(eb$adj.contact, want_a)
    ok_nm <- identical(rownames(eb$d.spatial), levs) && identical(colnames(eb$d.spatial), levs)
    cat(sprintf("crd_%-28s d.spatial=%-5s adj.contact=%-5s dimnames=%-5s (%d x %d)\n",
                nm, ok_d, ok_a, ok_nm, nrow(eb$d.spatial), ncol(eb$d.spatial)))
    if (!ok_d) { fails <<- fails + 1L; cat("   want:\n"); print(want_d); cat("   got:\n"); print(eb$d.spatial) }
    if (!ok_a) { fails <<- fails + 1L; cat("   want:\n"); print(want_a); cat("   got:\n"); print(eb$adj.contact) }
    if (!ok_nm) { fails <<- fails + 1L
      cat("   want dimnames:", paste(levs, collapse = ","), "\n   got  dimnames:",
          paste(rownames(eb$d.spatial), collapse = ","), "\n") }
  }
}

cat(sprintf("\n%s: %d failing comparisons out of %d, over %d configurations\n",
            if (fails == 0L) "IDENTICAL" else "MISMATCH", fails, n_cmp, length(configs)))
cat("note: options$run.time is wall-clock and is deliberately excluded from the gate.\n")
quit(status = if (fails == 0L) 0L else 1L)
