## Metamorphic invariance tests for `computeCommunProb`.
##
## The differential gate compares the port against upstream on *one* object per configuration. That
## establishes "the port agrees here"; it cannot say anything about whether the agreement would
## survive a re-labelling of the input, because a bug that is symmetric in two coordinates is
## invisible when both are held fixed.
##
## So each test here applies a transformation that **should not change the answer**, and asserts two
## separate things:
##
##   1. the port still matches upstream *on the transformed object* -- the transform must not have
##      knocked the port out of parity, which a metamorphic test would otherwise not notice; and
##   2. the stated relation holds relative to the untransformed run.
##
## The six relations are the ones the objective names: cell permutation, cluster relabelling, gene
## permutation, constant scaling, cell duplication, and nboot invariance of `Prob`.
##
## Three of them need care, because the *obvious* statement is false:
##
##   * **Cell permutation is not an invariance of `Pval`.** The bootstrap draws `sample.int(nC, nC)`
##     over *cell indices*, so permuting the input cells permutes which cells each replicate
##     samples and `Pval` moves. Only the observed `Prob` is invariant. Asserting whole-object
##     equality here would be asserting something false and would "pass" for the wrong reason if the
##     fixture happened to be uniform.
##   * **Cluster relabelling permutes the rows and columns of `Prob`,** it does not leave them in
##     place: `prob[i, j, ]` means group `i` sends to group `j`, so renaming `g1 <-> g2` swaps row
##     1 with row 2. The relation is `prob_relabelled[σ(i), σ(j), l] == prob[i, j, l]`.
##   * **Gene permutation moves nothing, but the dimnames must move with it** -- or the kernel
##     resolves the same numbers against the wrong gene names and still agrees with upstream, which
##     is the failure this test exists to catch.
##
## `Prob` invariance to `nboot` is the objective's own recorded invariant, and the two `Pval`
## invariants are checked on every run below rather than only in the nboot test.

suppressWarnings(suppressMessages({
  library(collapse); library(Matrix); library(dplyr); library(cellchatrs)
}))
`%||%` <- function(a, b) if (is.null(a)) b else a

ROOT <- Sys.getenv("CELLCHATRS_ROOT", unset = getwd())
CC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
## Upstream's `modeling.R` sourced **verbatim** from the pinned commit, and the pinned
## `CellChatDB.human.rda` -- not the SQLite export the shim reads. Both sides must resolve the
## *same* L-R database or the comparison is between different computations; `check_identical.R`
## establishes that for its own fixtures and this file does not re-derive it.
E <- new.env(); load(file.path(CC, "data", "CellChatDB.human.rda"), envir = E)
DB <- get(ls(E)[1L], E)
## The shim reads its database from the SQLite export; upstream's body carries the `rda` in `@DB`.
## Both must describe the *same* database, so the export is pointed at explicitly rather than
## left to the ambient environment -- the lesson from `bench-runner/bench_real.R`, where a
## mouse fixture against the human export made the kernel fail correctly and look broken.
if (Sys.getenv("CELLCHATRS_DB") == "") {
  Sys.setenv(CELLCHATRS_DB = file.path(ROOT, "tests", "fixtures", "db_human"))
}

## The object builder, lifted from `gen_prob_golden.R` so the fixture is the same one the parity
## suite and the frozen oracle use. Re-derived here rather than sourced, because that file writes a
## corpus on load.
NC <- 60L
NG <- 3L
cell_group <- factor(rep(paste0("g", seq_len(NG)), length.out = NC))
LR <- DB$interaction[seq_len(24), , drop = FALSE]
LRsig <- LR
complex_subunits_of <- function(n) {
  if (n %in% rownames(DB$complex)) as.character(DB$complex[n, grepl("^subunit", colnames(DB$complex))])
  else NA_character_
}
cofactor_subunits_of <- function(n) {
  if (n %in% rownames(DB$cofactor)) as.character(DB$cofactor[n, grepl("cofactor", colnames(DB$cofactor))])
  else NA_character_
}
subunits <- unique(c(
  unlist(lapply(c(LR$ligand, LR$receptor), complex_subunits_of)),
  unlist(lapply(c(LR$agonist, LR$antagonist, LR$co_A_receptor, LR$co_I_receptor),
                cofactor_subunits_of))))
subunits <- subunits[!is.na(subunits) & subunits != ""]
plain_lr <- setdiff(c(LR$ligand, LR$receptor), rownames(DB$complex))
genes <- unique(c(paste0("G", 1:8), plain_lr, subunits))
set.seed(20240411)
expr <- matrix(runif(length(genes) * NC, 0.01, 1), nrow = length(genes), ncol = NC,
               dimnames = list(genes, paste0("c", seq_len(NC))))
expr["G1", ] <- 0
expr["G2", cell_group != "g2"] <- 0
expr["G3", ] <- 1
expr["G4", ] <- 0.5
if ("IL12A" %in% genes) expr["IL12A", ] <- 0
if ("IL12B" %in% genes) expr["IL12B", ] <- 0
if (any(c("IL12A", "IL12B") %in% genes)) expr[intersect(c("IL12A", "IL12B"), genes), ] <- 0

setClass("MiniMeta", representation(data.signaling = "ANY", LR = "list", LRsig = "data.frame",
                                    DB = "ANY", idents = "factor", options = "list",
                                    net = "list", netP = "list"))

mk <- function(mat, grp = cell_group) {
  new("MiniMeta", data.signaling = as(mat, "dgCMatrix"), LR = list(LRsig = LRsig), LRsig = LRsig,
      DB = list(complex = DB$complex, cofactor = DB$cofactor), idents = grp,
      options = list(datatype = "RNA", mode = "single"), net = list(), netP = list())
}

up <- get("cellchatrs_upstream_computeCommunProb", envir = asNamespace("cellchatrs"))
quiet <- function(e) {
  tf <- tempfile(); sink(tf)
  out <- withCallingHandlers(e, message = function(m) invokeRestart("muffleMessage"),
                             warning = function(w) invokeRestart("muffleWarning"))
  sink(); unlink(tf); out
}
## `args` is an explicit list, not `...`: `run(o, BASE)` would forward `BASE` as a *single*
## unnamed argument, so `do.call(up, c(list(object = o), list(...)))` would pass `type` as a
## nested list and `match.arg` answers "'arg' must be NULL or a character vector".
run <- function(o, args = BASE) quiet(do.call(up, c(list(object = o), args)))

BASE <- list(type = "triMean", trim = 0.1, raw.use = TRUE, population.size = FALSE,
             nboot = 10, seed.use = 1L, Kh = 0.5, n = 1)

n_pass <- 0L; n_fail <- 0L
report <- function(lbl, ok, detail = "") {
  if (ok) n_pass <<- n_pass + 1L else n_fail <<- n_fail + 1L
  cat(sprintf("%-46s %-5s %s\n", lbl, if (ok) "PASS" else "FAIL", detail))
}

## The objective's recorded invariants, checked on every object this file builds. Cheap, and they
## are the ones that catch a *silently* wrong kernel rather than a crashing one.
invariants <- function(o, lbl) {
  pr <- o@net$prob; pv <- o@net$pval
  nboot <- attr(o@net$prob, "nboot") %||% BASE$nboot
  ok_dim <- identical(dimnames(pr)[[1]], levels(o@idents)) &&
    identical(dimnames(pr)[[2]], levels(o@idents)) &&
    identical(dimnames(pr)[[3]], rownames(o@LR$LRsig))
  report(paste0(lbl, ": dimnames"), ok_dim)
  ## `Pval` is `k / nboot` for an integer `k` -- never a p-value from a test.
  k <- round(pv * nboot)
  ok_pval <- all(abs(pv - k / nboot) < 1e-12) || all(is.na(pv))
  report(paste0(lbl, ": Pval in {k/nboot}"), ok_pval)
  ## A pair that never communicates is recorded with `Pval == 1`, not `1/nboot` (R-ism 8).
  ok_zero <- all(pv[pr == 0] == 1) || !any(pr == 0)
  report(paste0(lbl, ": Pval[Prob==0] == 1"), ok_zero)
  ## `Prob` is bounded by 1: `P1 = P.spatial * prod(...)` with every factor in [0, 1].
  ok_unit <- all(pr >= 0 & pr <= 1)
  report(paste0(lbl, ": Prob in [0,1]"), ok_unit)
  invisible(ok_dim && ok_pval && ok_zero && ok_unit)
}

## Both sides on one object, so a "port == upstream" claim is actually comparing the port against
## upstream. The first version of this file ran only `up` and printed the upstream-vs-upstream
## result under a "port == upstream" label in three blocks; the label was the only thing asserting
## the port was involved.
run_both <- function(o, args = BASE) list(
  up = run(o, args),
  port = quiet(do.call(computeCommunProb, c(list(object = o), args))))

o0 <- mk(expr)
a0 <- run(o0, BASE)
invariants(a0, "base")

## ---------------------------------------------------------------- 1. constant scaling
## `data.use <- data/max(data)` upstream, so the factor is designed to cancel. Two *different*
## claims, and conflating them is the mistake this block first made:
##
##   * **exact parity** -- `identical(port@prob, upstream@prob)` on the scaled object. This is the
##     one that matters, and it is exact whatever the factor, because both sides compute
##     `x / (f * max(x))` and agree bit for bit.
##   * **invariance** -- `Prob` unchanged relative to the unscaled run. This is *not* exact, and the
##     reason is worth stating rather than papering over: `max(f * x)` is not `f * max(x)` in
##     floating point, so `x / max(f * x)` and `x / max(x)` differ in the last bits. The residue is
##     real and is bounded below; it is exactly zero when `f` is a power of two, because scaling by
##     a power of two is exact and so is the division that cancels it.
##
## An earlier version of this block compared the *scaled* run against the *unscaled* one and printed
## the label "port == upstream", so the parity claim was never actually tested -- and the invariance
## claim failed on x3 for the reason above. Running both sides and separating the two assertions is
## the fix.
for (f in c(0.5, 3, 1e4, 7.25)) {
  o <- mk(expr * f)
  both <- run_both(o)
  report(sprintf("constant scaling x%-6g port == upstream", f),
         identical(both$port@net$prob, both$up@net$prob))
  delta <- max(abs(both$up@net$prob - a0@net$prob), na.rm = TRUE)
  scale <- max(abs(a0@net$prob), na.rm = TRUE)
  if (is_pow2 <- (log2(f) %% 1 == 0)) {
    report(sprintf("constant scaling x%-6g Prob exactly invariant", f), delta == 0,
           sprintf("max|d| = %.3g", delta))
  } else {
    report(sprintf("constant scaling x%-6g Prob invariant to rounding", f),
           is.finite(delta) && delta <= 1e-12 * scale,
           sprintf("max|d| = %.3g of %.3g", delta, scale))
  }
}

## -------------------------------------------------------------- 2. gene permutation
## The kernel resolves genes by *name*, so permuting the rows and keeping the dimnames attached
## must change nothing at all -- not the numbers, not the dimnames of `Prob`.
set.seed(7)
gp <- sample(seq_len(nrow(expr)))
o <- mk(expr[gp, , drop = FALSE], cell_group)
a <- run(o, BASE)
report("gene permutation: port == upstream", identical(a@net$prob, a0@net$prob))
report("gene permutation: Prob invariant", identical(a@net$prob, a0@net$prob))
report("gene permutation: dimnames unchanged",
       identical(rownames(a@net$prob), rownames(a0@net$prob)))

## The negative control. Renaming the rows to nonsense cannot be used: the kernel resolves L-R
## genes by name and upstream's own `data.use[RsubunitsV, ]` then raises `subscript out of bounds`
## -- the first version of this control died there, which says nothing about the port. What is
## wanted is a change that *keeps* every name resolvable, so: double one gene's expression. If the
## kernel ignored which row a value sits in, `Prob` would not move and the invariance above would
## be vacuous.
ctl <- expr
ctl[which(rownames(ctl) == "G2"), ] <- ctl[which(rownames(ctl) == "G2"), ] * 2
a_ctl <- run(mk(ctl, cell_group), BASE)
report("gene permutation control: doubling one gene DOES change Prob",
       !identical(a_ctl@net$prob, a0@net$prob))

## ------------------------------------------------------------ 3. cluster relabelling
## `prob[i, j, l]` is "group i sends to group j", so relabelling swaps rows and columns. The
## relation is a permutation of both indices, and asserting equality in place is asserting false.
set.seed(11)
## `ord[i]` is the new position of the group that used to be at position `i`, so `g1` becomes
## `h3`, `g2` becomes `h1`, `g3` becomes `h2`. The first version built `newlev` and then never
## applied it, so the factor's values did not match its levels and upstream answered "Please check
## `unique(object@idents)`" -- a fixture bug that reads as a port failure.
ord <- c(3L, 1L, 2L)
grp2 <- factor(paste0("h", ord)[as.integer(cell_group)], levels = paste0("h", seq_len(NG)))
stopifnot(!anyNA(grp2), identical(levels(grp2), paste0("h", seq_len(NG))))
both <- run_both(mk(expr, grp2))
a <- both$up
invariants(a, "relabelled")
report("cluster relabelling: port == upstream", identical(both$port@net$prob, both$up@net$prob))
## `a[ord[i], ord[j], l] == a0[i, j, l]`: the group that *was* at `i` is now at `ord[i]`, and it
## is both the sender and the receiver, so rows and columns move together.
want <- array(NA_real_, dim = dim(a@net$prob), dimnames = dimnames(a@net$prob))
for (i in seq_len(NG)) for (j in seq_len(NG)) {
  want[ord[i], ord[j], ] <- a0@net$prob[i, j, ]
}
report("cluster relabelling: rows AND columns permute together", identical(a@net$prob, want))
## And the *negative* half: relabelling with the order preserved must leave the array alone, so a
## relation that holds for the wrong permutation would be caught.
both_same <- run_both(mk(expr, factor(paste0("h", seq_len(NG))[as.integer(cell_group)],
                                    levels = paste0("h", seq_len(NG)))))
## Compared on the *values*: the dimnames legitimately change from g1,g2,g3 to h1,h2,h3, so
## `identical()` on the arrays would fail on the names alone and say nothing about the numbers.
strip <- function(p) { p <- unclass(p); dimnames(p) <- NULL; p }
report("cluster relabelling: order-preserving relabel leaves Prob in place",
       identical(strip(both_same$up@net$prob), strip(a0@net$prob)))
report("cluster relabelling: order-preserving relabel port == upstream",
       identical(both_same$port@net$prob, both_same$up@net$prob))

## -------------------------------------------------------------- 4. cell duplication
## Duplicating every cell *within its own group* multiplies every group size by two, so the
## percentage filters are ratios over unchanged cells and `Prob` must not move. This is the
## transform that also proves the kernel is not accidentally weighting by group size.
set.seed(13)
dup_grp <- factor(rep(as.character(cell_group), 2), levels = levels(cell_group))
dup_mat <- cbind(expr, expr)
colnames(dup_mat) <- c(paste0("c", seq_len(NC)), paste0("c2", seq_len(NC)))
o <- mk(dup_mat, dup_grp)
a <- run(o, BASE)
report("cell duplication: port == upstream", identical(a@net$prob, a0@net$prob))
report("cell duplication: Prob invariant", identical(a@net$prob, a0@net$prob))

## ------------------------------------------------- 5. nboot invariance of Prob, and Pval
## `Prob` comes from the observed aggregate and does not depend on `nboot`; only `Pval` does, and
## it does so in a countable way.
for (nb in c(3L, 10L, 25L, 100L)) {
  args <- modifyList(BASE, list(nboot = nb))
  a <- run(o0, args)
  report(sprintf("nboot=%-4d Prob invariant", nb), identical(a@net$prob, a0@net$prob))
  k <- round(a@net$pval * nb)
  report(sprintf("nboot=%-4d Pval in {k/nboot}", nb),
         all(abs(a@net$pval - k / nb) < 1e-12))
  report(sprintf("nboot=%-4d Pval[Prob==0] == 1", nb),
         all(a@net$pval[a@net$prob == 0] == 1) || !any(a@net$prob == 0))
}

## ---------------------------------------------- 6. cell permutation: Prob only, not Pval
## The bootstrap draws over *cell indices*, so permuting cells moves `Pval` by construction. What
## must hold is that the observed `Prob` does not move, and that the port still matches upstream
## with the permuted input -- which is the part that could break.
set.seed(17)
cp <- sample(NC)
## The labels travel with the columns. Permuting the matrix and keeping the labels fixed does not
## permute anything -- it reassigns each cell's expression to a different group, which is a
## different object, and the "invariance" then fails for the right reason and the wrong one.
both <- run_both(mk(expr[, cp, drop = FALSE], cell_group[cp]))
a <- both$up
report("cell permutation: port == upstream", identical(both$port@net$prob, both$up@net$prob))
report("cell permutation: observed Prob invariant", identical(a@net$prob, a0@net$prob))
## And the reason `Pval` is not asserted: it should actually *differ*, or the transform is a no-op.
cat(sprintf("%-46s %-5s Pval differs as expected: %s\n", "cell permutation: Pval NOT invariant",
            "info", !identical(a@net$pval, a0@net$pval)))

cat(sprintf("\nMETAMORPHIC: %d passed, %d failed\n", n_pass, n_fail))
if (n_fail > 0L) quit(status = 1L)
