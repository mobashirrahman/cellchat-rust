#!/usr/bin/env Rscript
## Golden corpus for `subsetCommunication_internal` -- the post-melt half: the eight DEG
## thresholds, the all-`NA` row drop, the `netP` aggregation and the final column selection.
##
## The internal function is called directly rather than through the exported
## `subsetCommunication`, because the exported one needs a `CellChat` S4 object and this
## branch never looks at one: given `net` as a data.frame it is pure row filtering, and the
## point of the corpus is the filtering.
##
## The inputs are **hand-built, not generated**. There is no RNG here -- every value is
## literal -- which is the whole reason this corpus is reproducible: the failures it pins
## (the `NA`-blanking rule, the `paste`-turns-`NA`-into-`"NA"` rule, the byte-wise group
## order, the `pval`-in-last-column order) are all deterministic given the table.
##
## Output is tab-separated, one record per line, NA escaped as \x01 and NaN as \x02 so that
## an empty string and a missing cell cannot be confused.
suppressWarnings(suppressMessages({
  library(Matrix); library(collapse); library(dplyr)
}))

## The generators are run from the repository root (`scripts/gen_fixtures.sh` does it, and
## so does every invocation in the docs). Deriving the root from `commandArgs()` is fragile:
## under `Rscript` the `--file=` argument is present but its spelling has changed between
## R versions, and a generator that cannot find its own output is worse than one that
## requires a known working directory.
ROOT <- Sys.getenv("CELLCHATRS_ROOT", unset = getwd())
if (!dir.exists(file.path(ROOT, "R"))) {
  stop("run from the repository root, or set CELLCHATRS_ROOT")
}
SRC <- Sys.getenv("CELLCHAT_SRC", "../CellChat")
env <- new.env(parent = globalenv())
for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R")) {
  sys.source(file.path(SRC, "R", f), envir = env, keep.source = FALSE)
}
internal <- get("subsetCommunication_internal", envir = env)

# Numeric columns come back as 17-digit text; re-parse them so the corpus round-trips
# bit-exactly rather than through R's 7-significant-digit default.
cell_txt <- function(df, i, nm) {
  v <- df[[nm]][i]
  if (is.na(v) && !is.nan(v)) return("<NA>")
  if (is.nan(v)) return("<NaN>")
  if (is.character(v)) return(gsub("\t", "<TAB>", v, fixed = TRUE))
  format(as.numeric(v), digits = 17, scientific = TRUE, trim = TRUE)
}

## ---------------------------------------------------------------- the fixture table
## Deliberately awkward in every way the function has an opinion about:
##   * group names g1, g2, g10 -- `group_by` orders these byte-wise, so "g10" < "g1"
##   * `NA` and `NaN` in the threshold columns, which are the same to R's comparisons
##   * one row that is entirely `NA`, dropped even with no threshold supplied
##   * `datasets` present, and a second table where it is present but `ligand.logFC` is not
mk_net <- function(extra = TRUE) {
  source <- c("g1", "g1", "g2", "g2", "g10", "g10", "g1", "g2", "g10", "g1", "g1", NA)
  target <- c("g1", "g2", "g1", "g2", "g1", "g2", "g10", "g10", "g10", "g2", "g1", "g1")
  n <- length(source)
  df <- data.frame(
    source = source, target = target,
    ligand = rep(c("L1", "L2"), length.out = n),
    receptor = rep(c("R1", "R2"), length.out = n),
    prob = c(0.5, 0.25, 0.125, 1, 0.0625, 0.03125, 2, 0.5, 0.75, 1/3, 0.2, 0.4),
    pval = c(0.01, 0.02, 0.5, 0.001, 0.2, 0.04, NA, 0.03, 0.9, 0.005, NaN, 0.06),
    interaction_name = paste0("LR", seq_len(n)),
    interaction_name_2 = rep(c("a", "b"), length.out = n),
    pathway_name = rep(c("pA", "pB", "pC"), length.out = n),
    annotation = rep(c("anno1", "anno2"), length.out = n),
    evidence = rep(c("e1", NA), length.out = n),
    stringsAsFactors = FALSE
  )
  if (extra) {
    df$ligand.pvalues <- c(0.01, 0.2, 0.3, 0.04, NA, 0.5, 0.02, 0.02, NaN, 0.1, 0.1, 0.1)
    df$ligand.logFC  <- c(0.5, -0.5, 1.5, 0, -2, 0.25, NA, 0.75, 0.1, 0, -1, 0)
    df$ligand.pct.1  <- c(0.9, 0.1, 0.5, 0.6, 0.2, 0.8, 0.4, 0.4, 0.4, 0.7, 0.3, 0.5)
    df$ligand.pct.2  <- c(0.1, 0.9, 0.2, 0.3, 0.8, 0.1, 0.5, 0.5, 0.5, 0.2, 0.6, 0.4)
    df$receptor.pvalues <- c(0.05, 0.06, 0.01, 0.02, 0.3, 0.04, 0.02, 0.02, 0.02, 0.01, 0.1, 0.1)
    df$receptor.logFC  <- c(0.5, 0.5, -1, 0.5, 0.5, 0.5, 0, 0.5, 0.5, 0.5, NA, 0.5)
    df$receptor.pct.1  <- c(0.8, 0.8, 0.4, 0.9, 0.1, 0.9, 0.5, 0.5, 0.5, 0.6, 0.2, 0.6)
    df$receptor.pct.2  <- c(0.2, 0.2, 0.6, 0.1, 0.9, 0.1, 0.5, 0.5, 0.5, 0.4, 0.8, 0.4)
    df$datasets <- c(rep("d1", 6), rep("d2", 6))
  }
  df
}
## A second table: `datasets` present, no `ligand.logFC`. Upstream's final `if/else` only
## reaches the `datasets` branch when `ligand.logFC` is *also* present, so this one drops
## `datasets` from the output. It looks like a bug and is not.
mk_net_nodeg <- function() {
  df <- mk_net(extra = FALSE)
  df$datasets <- c(rep("d1", 6), rep("d2", 6))
  df
}
## One row entirely `NA`, which `rowSums(is.na(x)) != ncol(x)` drops with no threshold set.
mk_net_allna <- function() {
  df <- mk_net()
  df[nrow(df), ] <- NA
  df
}
LEVELS <- c("g1", "g2", "g10")

cases <- list(
  list(name = "no_thresholds",          net = mk_net(),            slot = "net",  args = list()),
  list(name = "no_thresholds_nodatasets", net = mk_net_nodeg(),   slot = "net",  args = list()),
  list(name = "all_na_row_no_threshold", net = mk_net_allna(),    slot = "net",  args = list()),
  list(name = "ligand_pvalues",         net = mk_net(),            slot = "net",
       args = list(ligand.pvalues = 0.05)),
  list(name = "ligand_pvalues_none",    net = mk_net(),            slot = "net",
       args = list(ligand.pvalues = 0)),
  list(name = "ligand_logfc_pos",       net = mk_net(),            slot = "net",
       args = list(ligand.logFC = 0.5)),
  list(name = "ligand_logfc_neg",       net = mk_net(),            slot = "net",
       args = list(ligand.logFC = -0.5)),
  list(name = "ligand_logfc_zero",      net = mk_net(),            slot = "net",
       args = list(ligand.logFC = 0)),
  list(name = "ligand_pct1",            net = mk_net(),            slot = "net",
       args = list(ligand.pct.1 = 0.5)),
  list(name = "ligand_pct2",            net = mk_net(),            slot = "net",
       args = list(ligand.pct.2 = 0.5)),
  list(name = "receptor_pvalues",       net = mk_net(),            slot = "net",
       args = list(receptor.pvalues = 0.02)),
  list(name = "receptor_logfc_pos",     net = mk_net(),            slot = "net",
       args = list(receptor.logFC = 0.5)),
  list(name = "receptor_logfc_neg",     net = mk_net(),            slot = "net",
       args = list(receptor.logFC = -0.5)),
  list(name = "receptor_pct1",          net = mk_net(),            slot = "net",
       args = list(receptor.pct.1 = 0.5)),
  list(name = "receptor_pct2",          net = mk_net(),            slot = "net",
       args = list(receptor.pct.2 = 0.5)),
  ## All eight at once, twice. `permissive` sets every p-value threshold to the column's
  ## maximum and every logFC/pct threshold to its minimum, so the comparisons exclude
  ## nothing and the *only* rows lost are the ones blanked by an `NA`/`NaN` -- which is the
  ## case that isolates the blanking rule from the comparisons. `tight` keeps exactly one
  ## row, and which one is a direct read of all eight comparison directions at once.
  ## The logFC thresholds are `0`, not the column minimum. A *negative* argument takes the
  ## `<=` branch, so the minimum (-2) would keep exactly one row -- and that row is one of the
  ## poisoned ones, which empties the table and turns the case into an error. `0` takes the
  ## `>=` branch and keeps every reachable value, which is the most permissive setting that
  ## exists for this column.
  list(name = "all_eight_permissive",  net = mk_net(),            slot = "net",
       args = list(ligand.pvalues = 0.5, ligand.logFC = 0, ligand.pct.1 = 0.1,
                   ligand.pct.2 = 0.1, receptor.pvalues = 0.1, receptor.logFC = 0,
                   receptor.pct.1 = 0.1, receptor.pct.2 = 0.1)),
  ## One threshold, set to the column maximum, so no comparison excludes anything. The rows
  ## that go are exactly the two with a missing `ligand.pvalues` -- one `NA`, one `NaN`. This
  ## is the case that isolates the blanking rule with nothing else in play.
  list(name = "single_threshold_blanking", net = mk_net(),         slot = "net",
       args = list(ligand.pvalues = 0.5)),
  list(name = "all_eight_tight",       net = mk_net(),            slot = "net",
       args = list(ligand.pvalues = 0.05, ligand.logFC = 0.25, ligand.pct.1 = 0.85,
                   ligand.pct.2 = 0.05, receptor.pvalues = 0.05, receptor.logFC = 0.4,
                   receptor.pct.1 = 0.75, receptor.pct.2 = 0.15)),
  list(name = "datasets_d1",            net = mk_net(),            slot = "net",
       args = list(datasets = "d1")),
  list(name = "datasets_d2_nodeg",      net = mk_net_nodeg(),      slot = "net",
       args = list(datasets = "d2")),
  list(name = "sources_targets",        net = mk_net(),            slot = "net",
       args = list(sources.use = "g1", targets.use = c("g2", "g10"))),
  list(name = "sources_missing",        net = mk_net(),            slot = "net",
       args = list(sources.use = "gNA")),
  ## Errors. The two `stop()` wordings differ and both are pinned.
  list(name = "err_ligand_pvalues_missing", net = mk_net(extra = FALSE), slot = "net",
       args = list(ligand.pvalues = 0.05)),
  list(name = "err_ligand_logfc_missing",   net = mk_net(extra = FALSE), slot = "net",
       args = list(ligand.logFC = 0.5)),
  list(name = "err_receptor_pct2_missing",  net = mk_net(extra = FALSE), slot = "net",
       args = list(receptor.pct.2 = 0.5)),
  list(name = "err_datasets_missing",       net = mk_net(extra = FALSE), slot = "net",
       args = list(datasets = "d1")),
  ## Two absent columns at once: the *first* test's message is the one raised.
  list(name = "err_two_absent_order",       net = mk_net(extra = FALSE), slot = "net",
       args = list(datasets = "d1", ligand.pvalues = 0.05, receptor.logFC = 0.5)),
  list(name = "err_empty_result",           net = mk_net(),            slot = "net",
       args = list(ligand.pct.1 = 2)),
  ## `netP`: the aggregation, its column order, the byte-wise group order, the LONG_DOUBLE
  ## `mean`/`sum`, and `paste` turning a missing cell into the literal key "NA".
  list(name = "netp_plain",             net = mk_net(),             slot = "netP", args = list()),
  list(name = "netp_after_deg",         net = mk_net(),             slot = "netP",
       args = list(ligand.pvalues = 0.05, receptor.logFC = 0.4)),
  list(name = "netp_sources_targets",   net = mk_net(),             slot = "netP",
       args = list(sources.use = c("g1", "g10"), targets.use = "g2")),
  list(name = "netp_allna_row",         net = mk_net_allna(),      slot = "netP", args = list()),
  ## `slot.name` that is neither: upstream's `==` comparison falls through both branches
  ## of the final `if/else`, so every column survives.
  list(name = "unknown_slot",           net = mk_net(),             slot = "foo",  args = list()),
  list(name = "unknown_slot_deg",       net = mk_net(),             slot = "foo",
       args = list(ligand.pvalues = 0.05))
)

## ---------------------------------------------------------------- run and record
out <- file.path(ROOT, "tests/fixtures/subset_golden.txt")
q <- file(paste0(out, ".part"), "wt")
writeLines(sprintf("n_cases\t%d", length(cases)), q)
for (cs in cases) {
  df <- cs$net
  writeLines(sprintf("case\t%s", cs$name), q)
  writeLines(sprintf("in_cols\t%s\t%s", cs$name, paste(colnames(df), collapse = ",")), q)
  writeLines(sprintf("in_nrow\t%s\t%d", cs$name, nrow(df)), q)
  for (i in seq_len(nrow(df))) {
    writeLines(sprintf("in_row\t%s\t%d\t%s", cs$name, i,
      paste(vapply(colnames(df), function(nm) cell_txt(df, i, nm), ""), collapse = "\t")), q)
  }
  writeLines(sprintf("slot\t%s\t%s", cs$name, cs$slot), q)
  ## Thresholds are named after the corpus so the Rust side can read them back without a
  ## parallel hand-written list that could drift.
  a <- cs$args
  writeLines(sprintf("arg\t%s\tsources.use\t%s", cs$name,
    if (is.null(a$sources.use)) "-" else paste(a$sources.use, collapse = ",")), q)
  writeLines(sprintf("arg\t%s\ttargets.use\t%s", cs$name,
    if (is.null(a$targets.use)) "-" else paste(a$targets.use, collapse = ",")), q)
  for (nm in c("datasets")) {
    writeLines(sprintf("arg\t%s\t%s\t%s", cs$name, nm,
      if (is.null(a[[nm]])) "-" else paste(a[[nm]], collapse = ",")), q)
  }
  for (nm in c("ligand.pvalues", "ligand.logFC", "ligand.pct.1", "ligand.pct.2",
               "receptor.pvalues", "receptor.logFC", "receptor.pct.1", "receptor.pct.2")) {
    writeLines(sprintf("arg\t%s\t%s\t%s", cs$name, nm,
      if (is.null(a[[nm]])) "-" else format(as.numeric(a[[nm]]), digits = 17,
                                            scientific = TRUE, trim = TRUE)), q)
  }

  ## Warnings are captured, not suppressed: the empty-after-`sources.use` path *warns* where
  ## the empty-after-threshold path *errors*, and "identical cat/print progress messages"
  ## covers both. `withCallingHandlers` keeps the value while still seeing the condition.
  warns <- character(0)
  res <- tryCatch(
    withCallingHandlers(
      do.call(internal, c(list(net = df, LR = data.frame(), cells.level = LEVELS,
                               slot.name = cs$slot), a)),
      warning = function(w) {
        warns <<- c(warns, conditionMessage(w))
        invokeRestart("muffleWarning")
      }),
    error = function(e) structure(list(msg = conditionMessage(e)), class = "subset_error"))
  writeLines(sprintf("warn	%s\t%s", cs$name,
    ## Newlines escaped: tidyselect's deprecation message is two lines, and a raw newline
    ## would split the record in half, leaving a continuation line that no longer parses as a
    ## record at all -- which is a confusing failure three files away from the cause.
    if (length(warns)) {
      paste(gsub("\n", "<NL>", gsub("\t", "<TAB>", warns, fixed = TRUE)), collapse = " || ")
    }
    else "-"), q)

  if (inherits(res, "subset_error")) {
    writeLines(sprintf("error\t%s\t%s", cs$name,
      gsub("\t", "<TAB>", res$msg, fixed = TRUE)), q)
  } else {
    writeLines(sprintf("out_cols\t%s\t%s", cs$name, paste(colnames(res), collapse = ",")), q)
    writeLines(sprintf("out_nrow\t%s\t%d", cs$name, nrow(res)), q)
    for (i in seq_len(nrow(res))) {
      writeLines(sprintf("out_row\t%s\t%d\t%s", cs$name, i,
        paste(vapply(colnames(res), function(nm) cell_txt(res, i, nm), ""),
              collapse = "\t")), q)
    }
    ## `rownames` are part of the return value and are assigned *after* the empty check, so
    ## a zero-row result has none.
    writeLines(sprintf("out_rn\t%s\t%s", cs$name,
      paste(rownames(res), collapse = ",")), q)
  }
}
close(q)
if (!file.rename(paste0(out, ".part"), out)) stop("could not install ", out)
cat("wrote", length(cases), "cases to", out, "\n")
