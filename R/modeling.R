## Drop-in replacement for CellChat's `computeCommunProb`, backed by the Rust kernel.
##
## The contract is *bit-identical*, not "close": `net$prob` and `net$pval` must satisfy
## `identical()` against the pinned upstream implementation, and `options$parameter` and
## `options$run.time` are part of the object too. That is checked by
## `tests/parity/check_identical.R`, which runs both and compares.
##
## ## Why the R shim does so little
##
## It marshals and nothing else. `data.use <- data/max(data)`, the aggregation, the
## bootstrap permutation, the Hill functions and the p-value counting are all in `r-core`,
## tested bit-exactly against upstream in `src/rust/crates/r-core/tests/prob_parity.rs` over eight
## parameter configurations. Anything that stays here is something the kernel has not been
## given yet, and each such case is listed in `NOT_PORTED` below so it is visible rather
## than silently slow.
##
## ## Fallback
##
## `CELLCHATRS_FALLBACK=1` routes to `cellchatrs_upstream_computeCommunProb`, which is
## the pinned upstream body, verbatim. That is how the differential test obtains its
## reference, and it is the escape hatch if the Rust kernel is ever wrong on a machine we
## have not tested. The exact upstream sources are bundled with the package; see
## `cellchatrs_upstream_env()`.

#' @useDynLib CellChat, .registration = TRUE
NULL

## The pinned upstream body, kept verbatim for the differential test and as a fallback.
## Sourced from the pinned files shipped in inst/upstream; never modified.
cellchatrs_upstream_env <- function() {
  src <- system.file("upstream/CellChat-75253cd0", package = "CellChat")
  f <- file.path(src, "R", "modeling.R")
  if (!nzchar(src) || !file.exists(f)) {
    stop("the pinned CellChat source bundled with cellchatrs is missing", call. = FALSE)
  }
  ## Resolve imports exactly as the package namespace does, without attaching dependencies
  ## or inheriting bindings from the caller's global environment.
  env <- new.env(parent = parent.env(environment(cellchatrs_upstream_env)))
  ## `utilities.R` for `identifyOverExpressedGenes`/`extractGene`/`subsetData`/`subsetDB`,
  ## `database.R` for `CellChatDB`. Sourcing them all into `env` means upstream's own
  ## `subsetData` shadows nothing here but *is* what `cellchatrs_upstream_subsetData` and
  ## the `identifyOverExpressedGenes` fallback return, which is what the differential test
  ## needs. Anything upstream's copies call that this package also defines resolves to
  ## `env` first -- upstream, not the shim -- so the reference is genuinely upstream.
  ##
  ## `analysis.R` as well as `modeling.R`: `subsetCommunication` lives there, and omitting
  ## it makes `get("subsetCommunication", envir = env)` fall through to this package's
  ## namespace -- where the shim's own `subsetCommunication` lives. The result is silent,
  ## infinite recursion, and a C stack overflow rather than an error.
  ## `visualization.R` is here because the *reference* needs it, not because the port uses it:
  ## `rankNet(mode = "comparison")` builds a grouped bar chart and calls `CellChat_theme_opts()`,
  ## which is defined there, so the comparison branch cannot even be evaluated without it. The
  ## 14.1 scope decision keeps visualization in R; that is a statement about what the Rust core
  ## reimplements, not about which files the pinned reference is allowed to load. Sourcing it
  ## defines 31 functions and runs no top-level code.
  ## `CellChat_class.R` for `mergeCellChat` and `createCellChat`. It defines the `CellChat` S4
  ## class and both of those live in it, not in `analysis.R` despite what the objective's file list
  ## suggests. Sourcing it runs two `setClassUnion` calls and defines the class; no other
  ## top-level code. Without it, `mergeCellChat` is unreachable from the reference environment and
  ## the comparison-analysis path cannot be tested at all.
  ##
  ## `CellChat` and its unions are registered by this package at load time. Evaluate only the
  ## upstream functions here; rerunning `setClass` in .GlobalEnv would leak implementation
  ## bindings into the caller and emit Matrix class-registration warnings.
  for (f in c("modeling.R", "analysis.R", "utilities.R", "database.R", "visualization.R", "app.R")) {
    p <- file.path(src, "R", f)
    if (file.exists(p)) sys.source(p, envir = env, keep.source = FALSE)
  }
  ## Evaluate the functions from CellChat_class.R into the private reference environment. The
  ## class definitions and show method already live in the package namespace from load time.
  cf <- file.path(src, "R", "CellChat_class.R")
  class_expr <- parse(cf)
  for (i in 5L:length(class_expr)) eval(class_expr[[i]], envir = env)
  ## The upstream R body calls this helper by name. Its native implementation is supplied
  ## by this package, so buildSNN and network analysis need no upstream DLL.
  env$ComputeSNN <- ComputeSNN
  data_dir <- file.path(src, "data")
  for (data_file in list.files(data_dir, pattern = "\\.rda$", full.names = TRUE)) {
    load(data_file, envir = env)
  }
  env
}

#' The pinned upstream `computeCommunProb`, verbatim.
#'
#' Used as the reference in the differential test and as the fallback when
#' `CELLCHATRS_FALLBACK=1`. Sourced into a fresh environment on first use.
#' @export
cellchatrs_upstream_computeCommunProb <- function(object, ...) {
  env <- cellchatrs_upstream_cached()
  get("computeCommunProb", envir = env)(object, ...)
}

## Per-session cache for the sourced upstream environment. Sourcing `modeling.R` is cheap
## but not free, and the differential test calls this once per configuration.
## Uses an environment rather than a variable because a package namespace is locked.
.cellchatrs_cache <- new.env(parent = emptyenv())

cellchatrs_upstream_cached <- function() {
  if (is.null(.cellchatrs_cache$env)) .cellchatrs_cache$env <- cellchatrs_upstream_env()
  .cellchatrs_cache$env
}

#' The pinned upstream `mergeCellChat`, verbatim.
#'
#' `mergeCellChat` is **not** reimplemented in Rust, and the gate below says so at a named rung
#' rather than leaving it untested. It is pure S4 slot assembly -- `rbind` of `meta`, `cbind` of
#' `data`, a `union` of factor levels, and two `stop()` calls -- with no arithmetic anywhere, so
#' there is nothing to move to `r-core` and the 14.1 scope decision covers it.
#'
#' What is worth testing is what it *produces*: a merged object whose `data`, `idents` and `meta`
#' become the input to `computeCommunProb` in the comparison-analysis workflow, so
#' `tests/parity/check_merge.R` runs the Rust kernel on a merged object and requires
#' `identical()` to upstream on that same object. That is a real test, and unlike comparing this
#' function with itself it is not vacuous.
#' @export
cellchatrs_upstream_mergeCellChat <- function(...) {
  get("mergeCellChat", envir = cellchatrs_upstream_cached())(...)
}

#' The pinned upstream `createCellChat`, verbatim. Used only to build fixtures for
#' `check_merge.R`; the port has no counterpart because the shim never creates the object.
#' @export
cellchatrs_upstream_createCellChat <- function(...) {
  get("createCellChat", envir = cellchatrs_upstream_cached())(...)
}

#' The pinned upstream `computeAveExpr`, verbatim.
#' @export
cellchatrs_upstream_computeAveExpr <- function(...) {
  get("computeAveExpr", envir = cellchatrs_upstream_cached())(...)
}

#' The pinned upstream `filterCommunication`, verbatim.
#'
#' The reference for the differential test and the fallback under
#' `CELLCHATRS_FALLBACK=1`. Wrapped rather than reached for directly because
#' `cellchatrs_upstream_cached` is internal, and `get("filterCommunication", envir = env)`
#' from a script would otherwise find *this package's* shim -- the two functions share a
#' name, and the whole point of the comparison is that they differ.
#' @export
cellchatrs_upstream_filterCommunication <- function(object, ...) {
  get("filterCommunication", envir = cellchatrs_upstream_cached())(object, ...)
}

#' The pinned upstream `subsetData`, verbatim.
#' @export
cellchatrs_upstream_subsetData <- function(...) {
  get("subsetData", envir = cellchatrs_upstream_cached())(...)
}

#' The pinned upstream `subsetDB`, verbatim.
#' @export
cellchatrs_upstream_subsetDB <- function(...) {
  get("subsetDB", envir = cellchatrs_upstream_cached())(...)
}

#' The pinned upstream `aggregateNet`, verbatim.
#' @export
cellchatrs_upstream_aggregateNet <- function(...) {
  get("aggregateNet", envir = cellchatrs_upstream_cached())(...)
}

#' The pinned upstream `subsetCommunication`, verbatim.
#' @export
cellchatrs_upstream_subsetCommunication <- function(...) {
  get("subsetCommunication", envir = cellchatrs_upstream_cached())(...)
}

#' The pinned upstream `identifyOverExpressedGenes`, verbatim.
#' @export
cellchatrs_upstream_identifyOverExpressedGenes <- function(...) {
  get("identifyOverExpressedGenes", envir = cellchatrs_upstream_cached())(...)
}

#' `aggregateNet`, backed by the Rust kernel.
#'
#' Both branches are ported. The unfiltered one aggregates the probability array directly.
#'
#' The filtered one (`sources.use` / `targets.use` / `signaling` / `pairLR.use`) is ported
#' because it turned out to be *broken upstream* at this commit, and reproducing a bug is
#' only possible if you know it is one. `stringr::str_split(key, "|", simplify = TRUE)`
#' takes a regular expression, `|` is the alternation metacharacter, so `a[, 1]` is `""` and
#' `a[, 2]` is the key's first character; the level restriction then finds no cell group and
#' `tapply` yields a `0 x 0` matrix under `remove.isolate = TRUE` and a `k x k` matrix of
#' zeros under `remove.isolate = FALSE`. `net$count` and `net$weight` therefore carry no
#' information for **any** filter. See `docs/SEMANTICS.md` and
#' `r_core::net::aggregate_net_filtered`.
#'
#' The `message()` about isolate cell groups is emitted for `remove.isolate = TRUE` in the
#' filtered branch only, as upstream does.
#' @export
aggregateNet <- function(object, sources.use = NULL, targets.use = NULL, signaling = NULL,
                         pairLR.use = NULL, remove.isolate = TRUE, thresh = 0.05,
                         return.object = TRUE) {
  if (!is.null(sources.use) || !is.null(targets.use) || !is.null(signaling) ||
      !is.null(pairLR.use)) {
    ## The filtered branch. It is ported rather than deferred because it is on the analysis
    ## path -- and because what it computes is worth knowing: `stringr::str_split(key, "|")`
    ## takes a *regex*, `|` is alternation, so `a[, 1]` is `""` and `a[, 2]` is the key's
    ## first character. `net$count` and `net$weight` therefore come back as `0 x 0` under
    ## `remove.isolate = TRUE` and as a `k x k` matrix of zeros under `remove.isolate = FALSE`,
    ## for every filter. See `r_core::net::aggregate_net_filtered` and
    ## `docs/SEMANTICS.md`. Reproduced exactly, because a drop-in replacement that returned
    ## real numbers here would differ from upstream.
    df <- subsetCommunication(object, slot.name = "net", sources.use = sources.use,
                              targets.use = targets.use, signaling = signaling,
                              pairLR.use = pairLR.use, thresh = thresh)
    if (remove.isolate) {
      message("Isolate cell groups without any interactions are removed. To block it, set `remove.isolate = FALSE`")
    }
    cells.level <- levels(object@idents)
    res <- aggregate_net_filtered(
      table_source = as.character(df$source),
      table_target = as.character(df$target),
      table_prob = as.numeric(df$prob),
      cells_level = cells.level,
      remove_isolate = remove.isolate
    )
    net <- object@net
    net$count <- matrix(res$count, nrow = res$dim[1], ncol = res$dim[2], dimnames = res$dimnames)
    net$weight <- matrix(res$weight, nrow = res$dim[1], ncol = res$dim[2], dimnames = res$dimnames)
    object@net <- net
    return(if (return.object) object else net)
  }
  net <- object@net
  k <- dim(net$prob)[1]
  res <- aggregate_net(
    prob = as.numeric(net$prob),
    pval = as.numeric(net$pval),
    group_levels = dimnames(net$prob)[[1]],
    interaction_names = dimnames(net$prob)[[3]],
    thresh = thresh
  )
  count <- array(res$count, dim = res$dim)
  weight <- array(res$weight, dim = res$dim)
  dimnames(count) <- dimnames(weight) <- res$dimnames
  net$count <- count
  net$weight <- weight
  if (return.object) {
    object@net <- net
    object
  } else {
    net
  }
}

## `NaN`, **not** `NA_real_`, as the "argument absent" marker for the
## `subset_communication_deg` binding.
##
## `extendr` 0.9 rejects `NA_real_` for an `f64` parameter: the generated wrapper raises
## `Error::MustNotBeNA` during *argument conversion*, before the Rust body runs, and its
## message is the fixed string "Must not be NA." -- which names neither the argument nor the
## function. A `NA_real_` sentinel therefore aborts the call instead of reaching the kernel,
## and the resulting error looks like a bug in the port rather than in the sentinel. `NaN`
## converts cleanly, and every threshold upstream accepts is finite.
.cellchatrs_opt <- function(x) if (is.null(x)) NaN else as.numeric(x)

## Rebuild the data frame from the kernel's dynamically named list of parallel vectors.
##
## Two things have to be right here. The **column type** comes from the column name, not from
## the vector: `prob`, `pval` and the eight DEG columns are numeric and everything else is
## character, because a `format()`-ed number coming back as a string would make
## `identical()` fail on the column type even when every value matches. And the **empty** case
## *warns* rather than errors, because upstream's final `nrow(net) == 0` check runs after the
## `sources.use` filter, where the earlier threshold-stage `stop()` no longer applies.
.cellchatrs_subset_df <- function(res, net_levels, slot_name = "net") {
  ## The kernel reports its own errors as a value (see the binding): extendr 0.9 discards the
  ## text of an `Err` and raises "Must not be NA." instead, and the text is the contract.
  if (!is.null(res[["__error"]])) stop(res[["__error"]], call. = FALSE)
  n <- res[["__nrow"]]
  cols <- names(res)[!names(res) %in% c("__nrow", "__levels", "__error")]
  numeric_cols <- c("prob", "pval", "ligand.pvalues", "ligand.logFC", "ligand.pct.1",
                    "ligand.pct.2", "receptor.pvalues", "receptor.logFC", "receptor.pct.1",
                    "receptor.pct.2")
  if (n == 0L) {
    warning("No significant signaling interactions are inferred!")
    out <- as.data.frame(
      lapply(cols, function(nm) {
        v <- res[[nm]]
        if (nm %in% numeric_cols) numeric(0) else character(0)
      }), stringsAsFactors = FALSE)
    names(out) <- cols
    ## `BiocGenerics::as.data.frame` on the zero-row table, then no row names are assigned --
    ## the `if (nrow(net) == 0) warning(...) else rownames(net) <- 1:nrow(net)` split.
    return(out)
  }
  out <- as.data.frame(
    lapply(cols, function(nm) {
      v <- res[[nm]]
      if (nm %in% numeric_cols) as.numeric(v) else as.character(v)
    }), stringsAsFactors = FALSE, optional = TRUE)
  names(out) <- cols
  rownames(out) <- 1:nrow(out)
  if (!slot_name %in% c("net", "netP")) {
    attributes(out) <- attributes(out)[c("names", "row.names", "class")]
  }
  out
}

## `rankNet`, with the numerics backed by the Rust kernel.
#'
#' Only `mode = "single"` is Rust-backed, and only for the part that decides *which value goes in
#' which row*: the `thresh` cut, the `measure = "count"` binarisation, the group filters, the
#' per-pathway `apply(prob, 3, sum)`, the `-1/log` rescaling with its degenerate-entry
#' reassignment, and the `order()`. The ggplot stays here (`PLAN.md` 14.1) because the plot is
#' the deliverable in R and the numerics are what `identical()` can check.
#'
#' `mode = "comparison"`, `signaling.type`, and the `signaling` / `pairLR` filters still route to
#' upstream. The comparison mode needs the Wilcoxon kernel for `do.stat`, which is ported but not
#' yet wired in here; routing it out is honest and the gap is in `NOT_PORTED`.
#' @export
rankNet <- function(object, slot.name = "netP", measure = c("weight", "count"),
                    mode = c("comparison", "single"), comparison = c(1, 2), color.use = NULL,
                    stacked = FALSE, sources.use = NULL, targets.use = NULL, signaling = NULL,
                    pairLR = NULL, signaling.type = NULL, do.stat = FALSE, paired.test = TRUE,
                    cutoff.pvalue = 0.05, tol = 0.05, thresh = 0.05, show.raw = FALSE,
                    return.data = FALSE, x.rotation = 90, title = NULL, bar.w = 0.75,
                    font.size = 8, do.flip = TRUE, x.angle = NULL, y.angle = 0,
                    x.hjust = 1, y.hjust = 1, axis.gap = FALSE, ylim = NULL,
                    segments = NULL, tick_width = NULL, rel_heights = c(0.9, 0, 0.1)) {
  ## `match.arg` is R's, and it has to be: partial matching, the recorded matched name, and the
  ## ambiguous-prefix error are all part of the contract. Doing it here rather than in Rust keeps
  ## `options$parameter$measure` semantics intact for free.
  measure <- match.arg(measure)
  mode <- match.arg(mode)
  ## `comparison` is only meaningful in comparison mode. Routing on it here would send every
  ## `mode = "comparison"` call upstream -- `comparison` *defaults* to `c(1, 2)`, so the check
  ## `any(comparison != comparison[1])` is TRUE by default and upstream won unconditionally.
  needs_upstream <- !is.null(signaling) || !is.null(pairLR) || !is.null(signaling.type) ||
    (identical(mode, "single") && !is.null(comparison) && any(comparison != comparison[1]))
  if (needs_upstream) {
    return(cellchatrs_upstream_rankNet(
      object, slot.name, measure, mode, comparison, color.use, stacked, sources.use, targets.use,
      signaling, pairLR, signaling.type, do.stat, paired.test, cutoff.pvalue, tol, thresh,
      show.raw, return.data, x.rotation, title, bar.w, font.size))
  }

  if (identical(mode, "comparison")) {
    ## `thresh` is forwarded explicitly rather than read from the enclosing frame: the helper is
    ## a separate function, so reaching for `thresh` inside it is an unbound-variable error that
    ## only appears on the comparison path.
    return(.cellchatrs_ranknet_comparison(
      object, slot.name, measure, comparison, color.use, sources.use, targets.use,
      do.stat, paired.test, cutoff.pvalue, return.data, bar.w, font.size, title, thresh,
      stacked, show.raw, tol, do.flip, x.angle, y.angle, x.hjust, y.hjust, axis.gap, ylim,
      segments, tick_width, rel_heights))
  }

  object1 <- methods::slot(object, slot.name)
  object.names <- names(methods::slot(object, slot.name))
  prob <- object1$prob
  pval <- object1$pval
  dn <- dimnames(prob)
  k <- dim(prob)[1]
  n <- dim(prob)[3]

  ## `sources.use` and `targets.use` are validated against `dimnames(prob)[[1]]` -- the *first*
  ## axis -- for both. A `targets.use` naming a group that exists only on the second axis is
  ## rejected, and one naming a first-axis group that is not a target is accepted and then
  ## matches nothing. Both reproduced, because "fixing" it would change which calls work.
  resolve_groups <- function(x, which) {
    if (is.null(x)) return(NULL)
    if (is.character(x)) {
      if (!all(x %in% dn[[1]])) {
        stop(paste0("The input `", which, "` should be cell group names or a numerical vector!"),
             call. = FALSE)
      }
      return(match(x, dn[[1]]) - 1L)
    }
    if (is.numeric(x)) {
      ## `x - 1L` unchecked, as upstream: an out-of-range index is not an error there either,
      ## it just contributes no cells.
      return(as.integer(x) - 1L)
    }
    stop(paste0("The input `", which, "` should be cell group names or a numerical vector!"),
         call. = FALSE)
  }
  src <- resolve_groups(sources.use, "sources.use")
  tgt <- resolve_groups(targets.use, "targets.use")

  res <- ranknet_information_flow(
    prob = as.numeric(prob), pval = as.numeric(pval),
    group_levels = dn[[1]], names = dn[[3]],
    thresh = thresh, measure = measure,
    sources_use = if (is.null(src)) integer(0) else as.integer(src),
    targets_use = if (is.null(tgt)) integer(0) else as.integer(tgt))
  if (!is.null(res[["__error"]])) stop(res[["__error"]], call. = FALSE)

  pSum.original <- res[["pSum_original"]]
  pSum <- res[["pSum"]]
  pair.name <- dn[[3]]

  df <- data.frame(name = pair.name, contribution = pSum.original,
                   contribution.scaled = pSum, group = object.names[comparison[1]],
                   stringsAsFactors = FALSE)
  ## `idx <- with(df, order(df$contribution))`, on the **original** totals. `order()` on a
  ## double vector is radix and therefore stable, which is what makes the factor levels below
  ## reproducible when two pathways tie.
  df <- df[order(df$contribution), , drop = FALSE]
  df$name <- factor(df$name, levels = as.character(df$name))
  ## The zero-dropping loop. `df[-which(...), ]` keeps the *original* integer row names, so the
  ## gaps are observable in the result and `identical()` sees them.
  for (i in seq_along(pair.name)) {
    df.t <- df[df$name == pair.name[i], "contribution"]
    if (sum(df.t) == 0) {
      df <- df[-which(df$name == pair.name[i]), , drop = FALSE]
    }
  }

  ylabel <- if (measure == "weight") "Information flow" else "Number of interactions"
  gg <- ggplot2::ggplot(df, ggplot2::aes(x = .data$name, y = .data$contribution.scaled)) +
    ggplot2::geom_bar(stat = "identity", width = bar.w) +
    ggplot2::theme_classic() +
    ggplot2::theme(axis.text = ggplot2::element_text(size = font.size),
                   axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   axis.title.y = ggplot2::element_text(size = 10)) +
    ggplot2::xlab("") + ggplot2::ylab(ylabel) + ggplot2::coord_flip()
  if (!is.null(title)) {
    gg <- gg + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    ## `slot.name = "netP"` makes `signaling.contribution` a **tibble**, not a data frame.
    ##
    ## It is the one branch that passes through `dplyr::group_by |> summarize`, and that is what
    ## produces a `tbl_df`; every other branch ends at `BiocGenerics::as.data.frame` and stays a
    ## plain `data.frame`. The class is part of the return value, so `identical()` sees it and a
    ## port that returns a `data.frame` fails on the class alone with every value matching.
    ## Reproducing it needs `tibble::as_tibble` rather than a hand-built class vector, because a
    ## tibble's `row.names` are a compact `c(NA, -n)` integer vector rather than `1:n`.
    if (identical(slot.name, "netP")) {
      df <- tibble::as_tibble(df, .name_repair = "minimal")
    }
    return(list(signaling.contribution = df, gg.obj = gg))
  }
  gg
}

## `rankNet(mode = "comparison")`, backed by Rust.
##
## The numbers and the row order come from `r-core`; the data-frame assembly stays in R, which is
## where upstream does it and why:
##
## * `df[[i]] <- data.frame(name = pair.name.all, ..., row.names = pair.name.all)` followed by
##   `df[[i]][pair.name[[i]], 2] <- pSum.original[[i]]` is a **row-name-indexed** assignment, so a
##   pathway absent from a comparison keeps its 0. That is what makes a union of vocabularies
##   across datasets mean anything.
## * `df[[i]] <- df[[i]][idx, ]` keeps the *permuted* row names, and
##   `df[[1]]$contribution.data2 <- NULL` drops a column only the first comparison carried.
## * `do.call(rbind, df)` makes the row names unique across the comparisons -- observed to append
##   a bare `1`, `2`, ... with no separator (`X`, `B`, ..., `X1`, `B1`, ...), which is
##   `make.unique` behaviour and not what a hand-written `paste0(name, i)` would produce for the
##   first block.
## * `df$group <- factor(df$group, levels = rev(levels(df$group)))`: the group levels come back
##   **reversed**, because the grouped bar chart is built bottom-up.
##
## `do.stat` is the one thing still routed upstream: it runs `wilcox.test` / `kruskal.test` per
## pathway, which is a statistical test rather than a kernel and has no bit-identical requirement
## of its own.
.cellchatrs_ranknet_comparison <- function(object, slot.name, measure, comparison, color.use,
                                           sources.use, targets.use, do.stat, paired.test,
                                           cutoff.pvalue, return.data, bar.w, font.size, title,
                                           thresh, stacked, show.raw, tol, do.flip, x.angle,
                                           y.angle, x.hjust, y.hjust, axis.gap, ylim,
                                           segments, tick_width, rel_heights) {
  if (isTRUE(do.stat)) {
    ## Not a fallback so much as an unimplemented branch, and it is announced rather than
    ## silently skipped: a caller asking for p-values must not get a frame without them.
    message("cellchatrs: rankNet(do.stat = TRUE) -> pinned upstream (per-pathway tests)")
    return(cellchatrs_upstream_rankNet(
      object, slot.name, measure, "comparison", comparison, color.use, FALSE, sources.use,
      targets.use, NULL, NULL, NULL, do.stat, paired.test, cutoff.pvalue, 0.05, 0.05, FALSE,
      return.data, 90, title, bar.w, font.size))
  }
  oopts <- options(warn = -1)
  on.exit(options(oopts))
  slot_list <- methods::slot(object, slot.name)
  object.names <- names(slot_list)
  comparison <- as.integer(comparison)
  if (any(comparison < 1L) || any(comparison > length(slot_list))) {
    stop("`comparison` should be within the range of the compared objects!", call. = FALSE)
  }
  k <- dim(slot_list[[comparison[1]]]$prob)[1]

  ## `sources.use` / `targets.use` validated against `dimnames(prob)[[1]]` for **both**, as
  ## upstream does -- so a target-only name is accepted and then matches nothing.
  resolve <- function(x, which) {
    if (is.null(x)) return(integer(0))
    lev <- dimnames(slot_list[[comparison[1]]]$prob)[[1]]
    if (is.character(x)) {
      if (!all(x %in% lev)) {
        stop(paste0("The input `", which, "` should be cell group names or a numerical vector!"),
             call. = FALSE)
      }
      return(as.integer(match(x, lev)) - 1L)
    }
    if (is.numeric(x)) return(as.integer(x) - 1L)
    stop(paste0("The input `", which, "` should be cell group names or a numerical vector!"),
         call. = FALSE)
  }
  src <- resolve(sources.use, "sources.use")
  tgt <- resolve(targets.use, "targets.use")

  prob_list <- lapply(comparison, function(i) as.numeric(slot_list[[i]]$prob))
  pval_list <- lapply(comparison, function(i) as.numeric(slot_list[[i]]$pval))
  pair_names <- lapply(comparison, function(i) dimnames(slot_list[[i]]$prob)[[3]])

  res <- ranknet_comparison_flow(
    prob = unlist(prob_list, use.names = FALSE), prob_sizes = lengths(prob_list),
    pval = unlist(pval_list, use.names = FALSE), pval_sizes = lengths(pval_list),
    k = as.integer(k),
    pair_names = unlist(pair_names, use.names = FALSE),
    pair_name_sizes = lengths(pair_names),
    measure = measure, thresh = thresh,
    sources_use = as.integer(src), targets_use = as.integer(tgt))
  if (!is.null(res[["__error"]])) stop(res[["__error"]], call. = FALSE)

  ncomp <- length(comparison)
  pair.name.all <- res[["pair_names_all"]]

  ## `df[[i]] <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0,
  ## group = ..., row.names = pair.name.all)`, then the row-name-indexed assignment
  ## `df[[i]][pair.name[[i]], 2] <- pSum.original[[i]]`. Pathways the comparison does not have keep
  ## their 0 -- which is what makes a union of vocabularies across datasets mean anything.
  df <- lapply(seq_len(ncomp), function(i) {
    d <- data.frame(name = pair.name.all, contribution = 0, contribution.scaled = 0,
                    group = object.names[comparison[i]], row.names = pair.name.all,
                    stringsAsFactors = FALSE)
    nm <- res[[paste0("c", i - 1L, "_names")]]
    if (!is.null(nm) && length(nm)) {
      d[nm, 3] <- res[[paste0("c", i - 1L, "_scaled")]]
      d[nm, 2] <- res[[paste0("c", i - 1L, "_original")]]
    }
    d
  })

  ## `contribution.relative[[i]] <- as.numeric(format(df[[ncomp-i+1]]$contribution /
  ## df[[1]]$contribution, digits = 1))`, then `is.na -> 0`.
  ##
  ## `format` is applied **here**, in R, rather than in the kernel, and the reason is worth
  ## recording: `digits` is not per element. `format(1.2, digits = 1)` is `"1"`, but
  ## `format(c(0.04, 1.2), digits = 1)` is `c("0.04", "1.20")` -- the pair shares a two-decimal
  ## scale, so `as.numeric` returns `1.2` for the second element. No element-wise rule, including
  ## `signif`, reproduces that, and the `order()` below consumes the *formatted* values, so a
  ## kernel-side approximation would change which row comes first, not merely a digit.
  contribution.relative <- vector("list", max(0L, ncomp - 1L))
  for (j in seq_len(res[["n_ratio"]])) {
    r <- res[[paste0("ratio", j - 1L)]]
    fmt <- suppressWarnings(as.numeric(format(r, digits = 1)))
    fmt[is.na(fmt)] <- 0
    contribution.relative[[j]] <- fmt
    for (i in seq_len(ncomp)) df[[i]][[paste0("contribution.relative.", j)]] <- fmt
  }
  ## `idx <- with(df[[1]], order(...))` with the key count chosen by `length(comparison)`. R's
  ## `order`, because it is what upstream calls and it is radix-and-therefore-stable, which is
  ## what decides the row order when every sort key ties.
  keys <- list(-contribution.relative[[1]])
  if (ncomp == 2L) {
    keys <- list(-contribution.relative[[1]], df[[1]]$contribution, -df[[ncomp]]$contribution)
  } else if (ncomp == 3L) {
    keys <- list(-contribution.relative[[1]], -contribution.relative[[2]],
                 df[[1]]$contribution, -df[[ncomp]]$contribution)
  } else if (ncomp >= 4L) {
    keys <- list(-contribution.relative[[1]], -contribution.relative[[2]], -contribution.relative[[3]],
                 df[[1]]$contribution, -df[[ncomp]]$contribution)
  } else {
    stop("`rankNet` needs at least two objects in `comparison`", call. = FALSE)
  }
  idx <- do.call(order, keys)

  ## `df[[1]]$contribution.data2 <- df[[length(comparison)]]$contribution` -- only the first
  ## comparison carries it, and it is dropped again after the sort.
  df[[1]]$contribution.data2 <- df[[ncomp]]$contribution
  for (i in seq_len(ncomp)) {
    df[[i]] <- df[[i]][idx, , drop = FALSE]
    df[[i]]$name <- factor(as.character(df[[i]]$name), levels = as.character(df[[i]]$name))
  }
  df[[1]]$contribution.data2 <- NULL
  df <- do.call(rbind, df)
  df$group <- factor(as.character(df$group),
                     levels = object.names[comparison])
  df$group <- factor(as.character(df$group), levels = rev(levels(df$group)))

  ## The zero-dropping loop, shared by both modes and easy to miss because it sits in the tail
  ## *after* the `do.stat` block. `df` here is the rbind'd frame, so a pathway is dropped only if
  ## it contributes zero in **every** comparison -- and it is dropped from all of them at once,
  ## which changes `nrow` and leaves the row names with gaps. Omitting it is invisible on any
  ## fixture where every pathway has a positive total, and produces extra rows otherwise.
  for (i in seq_along(pair.name.all)) {
    df.t <- df[df$name == pair.name.all[i], "contribution"]
    if (sum(df.t) == 0) {
      df <- df[-which(df$name == pair.name.all[i]), , drop = FALSE]
    }
  }

  ## Upstream reverses `color.use` *before* it computes `colors.text`, so the first entry of the
  ## reversed vector is what a `contribution.relative > 1 + tol` row is painted. Reversing after
  ## would swap the two colours on every row that is coloured at all.
  if (is.null(color.use)) color.use <- .cellchatrs_gg_palette(ncomp)
  color.use <- rev(color.use)

  if (length(comparison) > 2L) {
    message("The text on the y-axis will not be colored for the number of compared datasets larger than 2!")
    colors.text <- NULL
  } else if (nrow(df) > 0L) {
    ## `ifelse(df$contribution.relative < 1 - tol, color.use[2], ifelse(> 1 + tol, color.use[1],
    ## "black"))`. Only with `do.stat` do the p-values gate the colour.
    rel <- df$contribution.relative.1
    colors.text <- ifelse(rel < 1 - tol, color.use[2],
                          ifelse(rel > 1 + tol, color.use[1], "black"))
  } else {
    colors.text <- character(0)
  }

  ## Upstream's plot construction, lifted. Visualization stays in R (PLAN.md 14.1), so this is a
  ## copy of upstream's ggplot2 calls rather than a reimplementation -- the point is that a
  ## drop-in returns the *same* `gg.obj` under `return.data = TRUE`, and a plot that merely looks
  ## right is a different object.
  ylabel <- if (measure == "weight") "Information flow" else "Number of interactions"
  ## The `y` aesthetic and the bar position are branched rather than computed, because upstream
  ## writes two separate `geom_bar()` calls and the deparsed call is stored in the layer's
  ## `constructor` field. A single call with an `if` inside it renders identically and still
  ## compares unequal.
  bar <- function(y_aes, pos) {
    ## `.data$name` rather than a bare `name`: inside `aes()` both are resolved in the data mask,
    ## but a bare symbol is also a *global* binding as far as `R CMD check` is concerned, and the
    ## result is a NOTE on every build. The mapping's quosure expression changes, which is why the
    ## gate compares `ggplot_build()$data` and not `identical(gg.obj)`.
    ggplot2::ggplot(df, ggplot2::aes(x = .data$name, y = .data[[y_aes]], fill = .data$group)) +
      ggplot2::geom_bar(stat = "identity", width = bar.w, position = pos)
  }
  if (stacked) {
    gg <- bar("contribution", "fill")
  } else if (isTRUE(show.raw)) {
    gg <- bar("contribution", ggplot2::position_dodge(0.8))
  } else {
    gg <- bar("contribution.scaled", ggplot2::position_dodge(0.8))
  }
  if (stacked) {
    gg <- gg + ggplot2::xlab("") +
      ggplot2::ylab(if (measure == "weight") "Relative information flow"
                    else "Relative number of interactions") +
      ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", color = "grey50", size = 0.5)
  } else {
    gg <- gg + ggplot2::xlab("") + ggplot2::ylab(ylabel)
  }
  if (isTRUE(axis.gap)) {
    ## `gg.gap` is an optional package upstream only reaches when `axis.gap = TRUE`. Announced
    ## rather than ignored: a caller who asked for the gap would otherwise get a plot without it
    ## and no indication why.
    message("cellchatrs: rankNet(axis.gap = TRUE) needs the optional 'gg.gap' package; the gap \
is not drawn")
  }
  ## `CellChat_theme_opts()` is visualization, so it is read from the pinned source rather than
  ## reimplemented -- which also guarantees the theme cannot drift from upstream's.
  gg <- gg + get("CellChat_theme_opts", envir = cellchatrs_upstream_cached())() +
    ggplot2::theme_classic()
  if (isTRUE(do.flip)) {
    if (is.null(x.angle)) x.angle <- 0
    gg <- gg + ggplot2::coord_flip() +
      ggplot2::theme(axis.text.y = ggplot2::element_text(colour = colors.text))
  } else {
    if (is.null(x.angle)) x.angle <- 45
    gg <- gg + ggplot2::scale_x_discrete(limits = rev) +
      ggplot2::theme(axis.text.x = ggplot2::element_text(colour = rev(colors.text)))
  }
  gg <- gg +
    ggplot2::theme(axis.text = ggplot2::element_text(size = font.size),
                   axis.title.y = ggplot2::element_text(size = font.size)) +
    ggplot2::scale_fill_manual(name = "", values = color.use) +
    ggplot2::guides(fill = ggplot2::guide_legend(reverse = TRUE)) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = x.angle, hjust = x.hjust),
                   axis.text.y = ggplot2::element_text(angle = y.angle, hjust = y.hjust))
  if (!is.null(title)) {
    gg <- gg + ggplot2::ggtitle(title) +
      ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5))
  }
  if (return.data) {
    df$contribution <- abs(df$contribution)
    df$contribution.scaled <- abs(df$contribution.scaled)
    return(list(signaling.contribution = df, gg.obj = gg))
  }
  gg
}

## Upstream's `ggPalette`, which lives in `visualization.R` and is therefore not exported here.
## The shim needs it only to give `ggplot` the same default fills the upstream plot would, and
## `ggPalette`'s first entries are fixed, so the common case is identical by construction.
.cellchatrs_gg_palette <- function(n) {
  grDevices::hcl.colors(max(1L, n), palette = "YlOrRd", rev = TRUE)
}

#' The pinned upstream `rankNet`, verbatim.
#' @export
cellchatrs_upstream_rankNet <- function(...) {
  get("rankNet", envir = cellchatrs_upstream_cached())(...)
}

## `subsetCommunication` for `mode = "single"`, backed by Rust.
#'
#' Two paths, split at the melt exactly as upstream splits it:
#'
#' * `net` is a 3-d array, no DEG threshold and `slot.name = "net"`: the melt, the `thresh`
#'   cut, the L-R join and the final column selection.
#' * otherwise -- `net` is already a data frame, or any of the eight DEG thresholds is given,
#'   or `slot.name = "netP"` -- the same melt, then the thresholds, the `netP` aggregation and
#'   the final selection. The melt runs *unselected* in this path, because upstream's
#'   `rowSums(is.na(net)) != ncol(net)` counts over the full column set.
#'
#' `mode = "merged"` and `pairLR.use` / `signaling` still route to upstream: `signaling` is a
#' database query (`searchPair`) rather than a numeric kernel, and `pairLR.use` needs the
#' tryCatch-over-two-columns fallback.
#' @export
subsetCommunication <- function(object = NULL, net = NULL, slot.name = "net",
                                sources.use = NULL, targets.use = NULL, signaling = NULL,
                                pairLR.use = NULL, thresh = 0.05, datasets = NULL,
                                ligand.pvalues = NULL, ligand.logFC = NULL,
                                ligand.pct.1 = NULL, ligand.pct.2 = NULL,
                                receptor.pvalues = NULL, receptor.logFC = NULL,
                                receptor.pct.1 = NULL, receptor.pct.2 = NULL) {
  if (!identical(object@options$mode, "single")) {
    return(cellchatrs_upstream_subsetCommunication(object, net, slot.name, sources.use,
                                                   targets.use, signaling, pairLR.use, thresh,
                                                   datasets, ligand.pvalues, ligand.logFC,
                                                   ligand.pct.1, ligand.pct.2,
                                                   receptor.pvalues, receptor.logFC,
                                                   receptor.pct.1, receptor.pct.2))
  }
  any_deg <- !is.null(datasets) || !is.null(ligand.pvalues) || !is.null(ligand.logFC) ||
    !is.null(ligand.pct.1) || !is.null(ligand.pct.2) || !is.null(receptor.pvalues) ||
    !is.null(receptor.logFC) || !is.null(receptor.pct.1) || !is.null(receptor.pct.2)
  ## `signaling` and `pairLR.use` still route to upstream: `signaling` goes through
  ## `searchPair`, which is a database query rather than a numeric kernel, and `pairLR.use`
  ## needs the tryCatch-over-two-columns fallback that upstream relies on.
  if (!is.null(pairLR.use) || !is.null(signaling)) {
    return(cellchatrs_upstream_subsetCommunication(object, net, slot.name, sources.use,
                                                   targets.use, signaling, pairLR.use, thresh,
                                                   datasets, ligand.pvalues, ligand.logFC,
                                                   ligand.pct.1, ligand.pct.2,
                                                   receptor.pvalues, receptor.logFC,
                                                   receptor.pct.1, receptor.pct.2))
  }
  if (is.null(net)) net <- object@net

  ## ------------------------------------------------------------------ the DEG / netP path
  ##
  ## Reached when `net` is already a data frame, or when any DEG threshold is supplied, or when
  ## `slot.name = "netP"`. The threshold stage runs on the *melted* table, so with an array
  ## `net` the kernel melts first and hands the result on unselected -- upstream's
  ## `rowSums(is.na(net)) != ncol(net)` counts over the full column set, and selecting first
  ## would change which rows count as all-`NA`.
  if (is.data.frame(net) || any_deg || !identical(slot.name, "net")) {
    LR <- object@LR$LRsig
    cells.level <- levels(object@idents)
    if (is.numeric(sources.use)) sources.use <- cells.level[sources.use]
    if (is.numeric(targets.use)) targets.use <- cells.level[targets.use]
    if (is.data.frame(net)) {
      cols <- colnames(net)
      ## `as.matrix()` on a data frame with any character column stringifies the *whole*
      ## frame, so numeric columns would go through `as.character` too. That is deliberate and
      ## lossless -- 17 significant digits round-trip every double -- and it is the only way to
      ## keep `NA` and `NaN` apart in one flat vector. `format(..., digits = 17)` is what
      ## guarantees the width; the default 7 would not.
      txt <- matrix(
        vapply(cols, function(nm) {
          v <- net[[nm]]
          if (is.numeric(v)) {
            out <- format(as.numeric(v), digits = 17, scientific = TRUE, trim = TRUE)
            out[is.na(v) & !is.nan(v)] <- "NA"
            out
          } else {
            out <- as.character(v)
            out[is.na(v)] <- "NA"
            out
          }
        }, character(nrow(net))),
        nrow = nrow(net), ncol = length(cols), byrow = FALSE)
      res <- subset_communication_deg(
        colnames = cols, cell_text = as.vector(txt), n_rows = nrow(net),
        prob = numeric(0), pval = numeric(0),
        group_levels = character(0), interaction_names = character(0),
        lr_columns_present = character(0),
        lr_interaction_name = character(0), lr_interaction_name_2 = character(0),
        lr_pathway_name = character(0), lr_ligand = character(0), lr_receptor = character(0),
        lr_annotation = character(0), lr_evidence = character(0),
        thresh = thresh,
        datasets = if (is.null(datasets)) character(0) else as.character(datasets),
        ligand_pvalues = .cellchatrs_opt(ligand.pvalues),
        ligand_logfc = .cellchatrs_opt(ligand.logFC),
        ligand_pct1 = .cellchatrs_opt(ligand.pct.1),
        ligand_pct2 = .cellchatrs_opt(ligand.pct.2),
        receptor_pvalues = .cellchatrs_opt(receptor.pvalues),
        receptor_logfc = .cellchatrs_opt(receptor.logFC),
        receptor_pct1 = .cellchatrs_opt(receptor.pct.1),
        receptor_pct2 = .cellchatrs_opt(receptor.pct.2),
        sources_use = if (is.null(sources.use)) character(0) else as.character(sources.use),
        targets_use = if (is.null(targets.use)) character(0) else as.character(targets.use),
        slot_name = slot.name)
    } else {
      present <- intersect(c("interaction_name_2", "pathway_name", "ligand", "receptor",
                             "annotation", "evidence"), colnames(LR))
      col_or_na <- function(nm) {
        if (!nm %in% present) return(rep("", nrow(LR)))
        v <- LR[[nm]]; out <- as.character(v); out[is.na(out)] <- ""; out
      }
      res <- subset_communication_deg(
        colnames = character(0), cell_text = character(0), n_rows = 0L,
        prob = as.numeric(net$prob), pval = as.numeric(net$pval),
        group_levels = dimnames(net$prob)[[1]],
        interaction_names = dimnames(net$prob)[[3]],
        lr_columns_present = present,
        lr_interaction_name = as.character(LR$interaction_name),
        lr_interaction_name_2 = col_or_na("interaction_name_2"),
        lr_pathway_name = col_or_na("pathway_name"),
        lr_ligand = col_or_na("ligand"), lr_receptor = col_or_na("receptor"),
        lr_annotation = col_or_na("annotation"), lr_evidence = col_or_na("evidence"),
        thresh = thresh,
        datasets = if (is.null(datasets)) character(0) else as.character(datasets),
        ligand_pvalues = .cellchatrs_opt(ligand.pvalues),
        ligand_logfc = .cellchatrs_opt(ligand.logFC),
        ligand_pct1 = .cellchatrs_opt(ligand.pct.1),
        ligand_pct2 = .cellchatrs_opt(ligand.pct.2),
        receptor_pvalues = .cellchatrs_opt(receptor.pvalues),
        receptor_logfc = .cellchatrs_opt(receptor.logFC),
        receptor_pct1 = .cellchatrs_opt(receptor.pct.1),
        receptor_pct2 = .cellchatrs_opt(receptor.pct.2),
        sources_use = if (is.null(sources.use)) character(0) else as.character(sources.use),
        targets_use = if (is.null(targets.use)) character(0) else as.character(targets.use),
        slot_name = slot.name)
    }
    ## Note the deliberate *absence* of a `tibble::as_tibble` wrap here, which is the mirror
    ## image of the one in `rankNet`. Measured against pinned upstream:
    ## `subsetCommunication(slot.name = "netP")` returns a plain `data.frame`, and
    ## `rankNet(slot.name = "netP", return.data = TRUE)` returns a `tbl_df`. Two netP paths,
    ## two different classes, and both are part of the return value, so `identical()` sees them.
    ## Adding the wrap here to "match" the `rankNet` one is exactly backwards: it made four
    ## fixtures fail on the class alone with every value, every rowname and every column class
    ## already matching. Concluding this from the `rankNet` evidence alone is the trap -- the two
    ## branches share a name and nothing else.
    out <- .cellchatrs_subset_df(res, net_levels = cells.level, slot_name = slot.name)
    return(out)
  }

  LR <- object@LR$LRsig
  cells.level <- levels(object@idents)
  ## Upstream resolves numeric indices against `cells.level` before comparing.
  if (is.numeric(sources.use)) sources.use <- cells.level[sources.use]
  if (is.numeric(targets.use)) targets.use <- cells.level[targets.use]

  ## `intersect(c(...), colnames(net))` sees which columns the L-R table *has*, not which
  ## cells are non-NA, so a table without an `evidence` column must drop it from the output
  ## while a table with an all-NA `evidence` column must keep it. `lr_columns_present`
  ## carries that distinction; `""` stands for NA within a present column.
  present <- intersect(c("interaction_name_2", "pathway_name", "ligand", "receptor",
                         "annotation", "evidence"), colnames(LR))
  col_or_na <- function(nm) {
    if (!nm %in% present) return(rep("", nrow(LR)))
    v <- LR[[nm]]
    out <- as.character(v)
    out[is.na(out)] <- ""
    out
  }
  res <- subset_communication(
    prob = as.numeric(net$prob),
    pval = as.numeric(net$pval),
    group_levels = dimnames(net$prob)[[1]],
    interaction_names = dimnames(net$prob)[[3]],
    lr_columns_present = present,
    ## `LR$interaction_name`, **not** `rownames(LR)`. The melt's last dimension is labelled
    ## from the column; `rownames` of a default data.frame is `"1"`, `"2"`, ... , so using it
    ## silently mislabels every interaction. It went unnoticed because the fixtures used for
    ## `subsetCommunication` parity all had a row-name column that happened to agree, and
    ## because `identical()` on the *values* passed -- the L-R name only appears in the
    ## `interaction_name` output column, which the earlier gate did compare, but with a
    ## fixture whose `rownames` equalled `interaction_name`.
    lr_interaction_name = as.character(LR$interaction_name),
    lr_interaction_name_2 = col_or_na("interaction_name_2"),
    lr_pathway_name = col_or_na("pathway_name"),
    lr_ligand = col_or_na("ligand"),
    lr_receptor = col_or_na("receptor"),
    lr_annotation = col_or_na("annotation"),
    lr_evidence = col_or_na("evidence"),
    thresh = thresh,
    sources_use = if (is.null(sources.use)) character(0) else as.character(sources.use),
    targets_use = if (is.null(targets.use)) character(0) else as.character(targets.use)
  )
  if (res$nrow == 0L) {
    stop("No significant signaling interactions are inferred based on the input!", call. = FALSE)
  }
  ## Three columns are factors upstream, not two, and getting this wrong is invisible in the row
  ## count. Upstream melts the `K x K x N` array and calls `var.convert` on the result, which factors
  ## every column that came from a *dimname*: `source` and `target` from the first two dimensions,
  ## and `interaction_name` from the third. Everything else -- `ligand`, `receptor`, `annotation`,
  ## `evidence` -- arrives from the join with `LR`, is not a dimname, and stays character.
  ##
  ## This was found by `tests/parity/check_merge.R`, not by inspection: the port returned
  ## `interaction_name` as character and the gate reported a class mismatch on a table whose values,
  ## columns, row count and row names all agreed. An earlier version of this comment claimed only
  ## `source` and `target` were factors and that the rest "stays character" -- correct about the join
  ## columns, wrong about `interaction_name`, and the kind of claim that survives because nothing
  ## compares `class()`.
  ##
  ## The factor levels are the **whole** vector, not the subset that survived the threshold:
  ## `interaction_levels` is the array's third dimnames in order, so an interaction filtered out of
  ## the table still contributes an unused level. Reproducing that matters -- `droplevels` on the
  ## result would not be `identical()`.
  df <- as.data.frame(lapply(res[res$colnames], identity), stringsAsFactors = FALSE)
  names(df) <- res$colnames
  df$source <- factor(df$source, levels = res$source_levels)
  df$target <- factor(df$target, levels = res$target_levels)
  if (!is.null(res$interaction_levels) && "interaction_name" %in% names(df)) {
    df$interaction_name <- factor(df$interaction_name, levels = res$interaction_levels)
  }
  rownames(df) <- seq_len(nrow(df))
  ## Upstream ends with `BiocGenerics::as.data.frame(net, stringsAsFactors = FALSE)`, whose
  ## attribute pairlist is ordered `names, class, row.names`. Each `df$col <- ...` above is a
  ## `$<-.data.frame`, which rebuilds the frame and leaves it `names, row.names, class`.
  ## `identical()` treats attributes as a set so both orders compare equal in principle, but the
  ## gate at the end of `tutorial_repro.R` observed `identical(a, b)` FALSE on frames whose every
  ## attribute and column matched by name -- and `serialize()` disagrees outright, since it writes
  ## the pairlist in stored order. Matching upstream's order makes both comparisons agree, and
  ## costs one reorder.
  at <- attributes(df)
  attributes(df) <- at[c("names", "class", "row.names")]
  df
}

#' `computeAveExpr`, backed by the Rust kernel.
#'
#' Handles `type` in `c("triMean", "truncatedMean", "median")` with R's `match.arg`
#' partial matching, and `features` in R's `intersect` order.
#' @export
computeAveExpr <- function(object, features = NULL, group.by = NULL,
                           type = c("triMean", "truncatedMean", "median"), trim = NULL,
                           slot.name = c("data.signaling", "data"), data.use = NULL) {
  type <- match.arg(type)
  slot.name <- match.arg(slot.name)
  if (is.null(data.use)) data.use <- slot(object, slot.name)
  if (is.null(features)) {
    features.use <- rownames(data.use)
  } else {
    features.use <- intersect(features, rownames(data.use))
  }
  data.use <- as.matrix(data.use)
  if (is.null(group.by)) {
    labels <- object@idents
    if (!is.factor(labels)) {
      message("Use the joint cell labels from the merged CellChat object")
      labels <- object@idents$joint
    }
  } else {
    labels <- object@meta[[group.by]]
  }
  if (!is.factor(labels)) labels <- factor(labels)
  ## `trim = NULL` is never substituted, so `mean(x, trim = NULL, na.rm = TRUE)` is R's
  ## `trim = 0` -- i.e. no trimming. The kernel's own default is 0.1 and is a different
  ## function; do not unify them.
  res <- compute_ave_expr(
    data = as.numeric(t(data.use)),
    genes = rownames(data.use),
    group = as.integer(labels) - 1L,
    group_levels = levels(labels),
    features = if (is.null(features)) character(0) else as.character(features),
    type_ = type,
    trim = if (is.null(trim)) 0 else trim
  )
  out <- matrix(res$values, nrow = res$dim[1], ncol = res$dim[2],
                dimnames = list(res$features, res$groups))
  out
}

#' `subsetDB`, backed by the Rust kernel.
#'
#' Only `key = "annotation"` is routed here; any other key falls back, because the message
#' for an unknown key is part of the contract.
#' @export
subsetDB <- function(CellChatDB, search = c(), key = "annotation", non_protein = FALSE) {
  if (!identical(key, "annotation")) {
    return(cellchatrs_upstream_subsetDB(CellChatDB, search, key, non_protein))
  }
  ## The default substitution is `is.null(search)`, **not** `length(search) == 0`:
  ##   if (is.null(search) & non_protein == FALSE & any(key == "annotation")) {
  ##     search <- c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
  ##   } else if (is.null(search) & non_protein == TRUE & any(key == "annotation")) { ... }
  ## So `search = character(0)` -- the *declared* default of the Rust-facing helper, and
  ## what a caller who writes `search = c()` gets -- keeps its empty value and selects
  ## **zero** rows, while `search = NULL` selects 2239. Collapsing the two is silent and
  ## wrong; it is done here, in R, because only R can tell `NULL` from `character(0)`.
  if (is.null(search)) {
    search <- if (isTRUE(non_protein)) {
      c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact", "Non-protein Signaling")
    } else {
      c("Secreted Signaling", "ECM-Receptor", "Cell-Cell Contact")
    }
  }
  ann <- as.character(CellChatDB$interaction$annotation)
  res <- subset_db_by_annotation(
    annotations = ann,
    search = as.character(search),
    key_is_annotation = TRUE,
    non_protein = non_protein
  )
  if (isTRUE(res$non_protein) && !non_protein) {
    message("The non-protein signaling is now included for CellChat analysis, which is usually used for neuron-neuron and metabolic communication!")
  }
  CellChatDB$interaction <- CellChatDB$interaction[res$keep, , drop = FALSE]
  CellChatDB
}

#' `subsetData`, backed by the Rust kernel for its gene list and annotation reordering.
#' @export
subsetData <- function(object, features = NULL) {
  ## The annotation reordering is a no-op on a `CellChatDB` that is already in factor order,
  ## which is what `cellchatrs`'s export is. It is still done here so the object is correct
  ## for a user-supplied DB, and so the behaviour matches for a DB that is not sorted.
  if (object@options$datatype != "RNA") {
    if (!("annotation" %in% colnames(object@DB$interaction))) {
      warning("A column named `annotation` is required in `object@DB$interaction` when running CellChat on spatial transcriptomics! The `annotation` column is now automatically added and all L-R pairs are assigned as `Secreted Signaling`, which means that these L-R pairs are assumed to mediate diffusion-based cellular communication.")
      object@DB$interaction$annotation <- "Secreted Signaling"
    }
  }
  if ("annotation" %in% colnames(object@DB$interaction)) {
    if (length(unique(object@DB$interaction$annotation)) > 1) {
      a <- object@DB$interaction$annotation
      lev <- c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling", "Cell-Cell Contact")
      object@DB$interaction$annotation <- factor(a, levels = lev)
      object@DB$interaction <- object@DB$interaction[order(object@DB$interaction$annotation), , drop = FALSE]
      object@DB$interaction$annotation <- as.character(object@DB$interaction$annotation)
    }
  }
  gene.use.input <- .cellchatrs_gene_list(object@DB)
  gene.use <- subset_data_gene_use(
    gene_use_input = gene.use.input,
    data_rownames = rownames(object@data),
    features = if (is.null(features)) character(0) else as.character(features)
  )
  object@data.signaling <- object@data[rownames(object@data) %in% gene.use, ]
  object
}

## `extractGene` for whatever `DB` the object carries. Uses the pinned upstream when the
## object holds a real `CellChatDB`, and falls back to it rather than reimplementing.
.cellchatrs_gene_list <- function(DB) {
  if (!is.null(get0("extractGene", envir = cellchatrs_upstream_cached()))) {
    return(get("extractGene", envir = cellchatrs_upstream_cached())(DB))
  }
  stop("the pinned upstream is needed for extractGene; set CELLCHAT_SRC", call. = FALSE)
}

## What the Rust kernel does not cover yet, so the shim still has to do it in R.
## Each entry is a real behavioural difference, not a stylistic one, and each is on the
## critical path for the drop-in claim.
NOT_PORTED <- c(
  spatial        = "computeRegionDistance preserves the original Annoy algorithm; cellchatrs_computeRegionDistance_exact is an explicit alternative, outside the compatibility contract",
  rankNet        = "rankNet(mode = 'comparison') is Rust-backed for the information flow and stays in R for `format`/row assembly, as does rankNetPairwise's data-frame construction; `rankNetPairwise`'s ordering is Rust. mergeCellChat is still R throughout",
  mergeCellChat  = "mergeCellChat: still R throughout -- a list of S4 objects, slot combination and net reindexing, with no arithmetic to move",
  oeg_fast       = "identifyOverExpressedGenes(do.fast = TRUE): presto, a different algorithm, falls back",
  oeg_dataset    = "identifyOverExpressedGenes(group.dataset != NULL) is ported for do.fast = FALSE -- the Wilcoxon branch honours group.dataset too (utilities.R:512-520), so it was never blocked on presto. The do.fast = TRUE dataset branch is presto's and still falls back",
  progress_bar   = "txtProgressBar output: not reproduced (not part of identical())",
  run_time       = "options$run.time is wall-clock and therefore necessarily not identical"
)

#' `computeCommunProb`, backed by the Rust kernel.
#'
#' The whole numeric surface is in `r-core`: the group means, the `i x b` bootstrap
#' permutation test, the Hill functions, the p-value counting and the `P.spatial`
#' loop-carried mutation. See the file header for why the shim does so little, and
#' `docs/SEMANTICS.md` for the R-isms that are reproduced verbatim.
#'
#' `datatype = "spatial"` routes to the pinned upstream body; see `NOT_PORTED`.
#'
#' `raw.use = FALSE` is Rust-backed. It reads `object@data.smooth` rather than
#' `object@data.signaling` and then takes exactly the same path, because upstream's own
#' `if (raw.use) ... else ...` is two lines that choose a matrix and there is no separate kernel
#' afterwards. What this package does *not* implement is `projectData`, which is what produces
#' `data.smooth` in the first place -- so the caller has to supply the slot, and the corpus
#' exercises the path with an independently constructed one.
#' @export
computeCommunProb <- function(object,
                              type = c("triMean", "truncatedMean", "thresholdedMean", "median"),
                              trim = 0.1, LR.use = NULL, raw.use = TRUE, population.size = FALSE,
                              distance.use = TRUE, interaction.range = 250,
                              scale.distance = 0.01, k.min = 10, contact.dependent = TRUE,
                              contact.range = NULL, contact.knn.k = NULL,
                              contact.dependent.forced = FALSE, do.symmetric = TRUE,
                              nboot = 100, seed.use = 1L, Kh = 0.5, n = 1) {
  ## `match.arg` first, so the error message and the recorded `type.mean` are upstream's.
  type <- match.arg(type)
  ## `CELLCHATRS_FALLBACK=1` sends every call to pinned upstream. The spatial branch used to route
  ## here unconditionally because `P.spatial` had no port; it now has one, so the escape hatch is
  ## the explicit opt-in it should have been rather than the only way to run spatial data.
  ## `nzchar(Sys.getenv("CELLCHATRS_FALLBACK", "0"))` looks right and is **always TRUE**: the
  ## unset default is the string `"0"`, which has four characters, so the escape hatch was
  ## permanently engaged and every call went to pinned upstream. The differential gate compares
  ## the shim against upstream, so it stayed green for the whole run -- comparing upstream with
  ## itself -- and every `prob.*` and `computeCommunProb` quantity was vacuous. Silent, and
  ## invisible from a green gate: the only symptom was that the Rust path was never executed.
  ## The test is an explicit affirmative, not a presence test.
  fallback <- tolower(Sys.getenv("CELLCHATRS_FALLBACK", "0")) %in% c("1", "true", "yes", "on")
  rng_kind <- base::RNGkind()
  ## The Rust bootstrap reproduces R's Mersenne-Twister + rejection-sampling stream. Preserve
  ## upstream behavior for alternate RNG/sample kinds by delegating the whole call.
  rng_supported <- identical(rng_kind[[1L]], "Mersenne-Twister") &&
    identical(rng_kind[[3L]], "Rejection")
  if (fallback || !rng_supported) {
    return(cellchatrs_upstream_computeCommunProb(object, type = type, trim = trim,
                                                 LR.use = LR.use, raw.use = raw.use,
                                                 population.size = population.size,
                                                 distance.use = distance.use,
                                                 interaction.range = interaction.range,
                                                 scale.distance = scale.distance,
                                                 k.min = k.min,
                                                 contact.dependent = contact.dependent,
                                                 contact.range = contact.range,
                                                 contact.knn.k = contact.knn.k,
                                                 contact.dependent.forced = contact.dependent.forced,
                                                 do.symmetric = do.symmetric,
                                                 nboot = nboot, seed.use = seed.use,
                                                 Kh = Kh, n = n))
  }

  ptm <- Sys.time()
  cat(type, "is used for calculating the average gene expression per cell group.", "\n")

  data <- if (raw.use) as.matrix(object@data.signaling) else as.matrix(object@data.smooth)
  pairLR.use <- if (is.null(LR.use)) object@LR$LRsig else .order_lr_use(LR.use)
  group <- object@idents
  nLR <- nrow(pairLR.use)
  numCluster <- nlevels(group)
  if (numCluster != length(unique(group))) {
    stop("Please check `unique(object@idents)` and ensure that the factor levels are correct!\n",
         "You may need to drop unused levels using 'droplevels' function. e.g.,\n",
         "`meta$labels = droplevels(meta$labels, exclude = setdiff(levels(meta$labels),unique(meta$labels)))`",
         call. = FALSE)
  }
  nC <- ncol(data)

  complex_input <- object@DB$complex
  cofactor_input <- object@DB$cofactor
  complex_cols <- grep("^subunit", names(complex_input), ignore.case = TRUE)
  cofactor_cols <- grep("^cofactor", names(cofactor_input), ignore.case = TRUE)

  ## The spatial constraint, upstream's `if (object@options$datatype != "RNA")` block in full.
  ##
  ## Preserve upstream Annoy neighbours before applying the Rust probability kernel.
  ##
  ## The R-side transformation is kept in R rather than pushed into the kernel because it is mostly
  ## *messages*: three `cat`s, two `print`s whose text embeds `Sys.time()`, and a `stop()` whose
  ## message embeds a computed `format(1/d.min, digits = 2)`. `options$parameter` and the saved
  ## `d.spatial` are part of the drop-in object, so the branch also has to leave both exactly as
  ## upstream does -- including NULLED spatial parameters on the RNA side.
  P.spatial <- NULL
  adj.contact <- NULL
  d.spatial <- NULL
  nLR1 <- nLR
  if (object@options$datatype != "RNA") {
    data.spatial <- object@images$coordinates
    if ("spatial.factors" %in% names(object@images)) {
      ratio <- object@images$spatial.factors$ratio
      tol <- object@images$spatial.factors$tol
    } else {
      stop("`object@images$spatial.factors` is missing. Please update the object via `updateCellChat`! \n")
    }
    meta.t <- data.frame(group = group, samples = object@meta$samples,
                         row.names = rownames(object@meta))
    res.sp <- computeRegionDistance(
      coordinates = data.spatial, meta = meta.t, interaction.range = interaction.range,
      ratio = ratio, tol = tol, k.min = k.min, contact.dependent = contact.dependent,
      contact.range = contact.range, contact.knn.k = contact.knn.k)
    d.spatial <- res.sp$d.spatial   # NaN if no nearby cell pairs
    adj.contact <- res.sp$adj.contact # zeros if no nearby cell pairs
    if (distance.use) {
      print(paste0(">>> Run CellChat on spatial transcriptomics data using distances as constraints of the computed communication probability <<< [", Sys.time(), "]"))
      d.spatial <- d.spatial * scale.distance
      diag(d.spatial) <- NaN
      d.min <- min(d.spatial, na.rm = TRUE)
      if (d.min < 1) {
        cat("The suggested minimum value of scaled distances is in [1,2], and the calculated value here is ", d.min, "\n")
        stop("Please increase the value of `scale.distance` and use a value that is slighly smaller than ", format(1/d.min, digits = 2), "\n")
      }
      P.spatial <- 1/d.spatial
      P.spatial[is.na(d.spatial)] <- 0
      diag(P.spatial) <- max(P.spatial)
      d.spatial <- d.spatial/scale.distance # This is only for saving the data
    } else {
      print(paste0(">>> Run CellChat on spatial transcriptomics data without distance values as constraints of the computed communication probability <<< [", Sys.time(), "]"))
      P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
      P.spatial[is.na(d.spatial)] <- 0 # diagonal is 1
    }
    ## Upstream's `if (object@options$datatype == "RNA") { nLR1 <- nLR } else { ... }`, which is
    ## what decides whether `P.spatial` is multiplied by `adj.contact` for *all* L-R pairs or only
    ## the contact-dependent ones past `nLR1`. The `cat` texts are part of the observable output.
    if (isTRUE(contact.dependent.forced)) {
      cat("Force to run CellChat in a `contact-dependent` manner for all L-R pairs including secreted signaling.\n")
      P.spatial <- P.spatial * adj.contact
      nLR1 <- nLR
    } else if (isTRUE(contact.dependent) && length(unique(pairLR.use$annotation)) > 0) {
      ann <- unique(pairLR.use$annotation)
      if (all(ann %in% c("Cell-Cell Contact"))) {
        cat("All the input L-R pairs are `Cell-Cell Contact` signaling. Run CellChat in a contact-dependent manner. \n")
        P.spatial <- P.spatial * adj.contact
        nLR1 <- nLR
      } else if (all(ann %in% c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling"))) {
        cat("Molecules of the input L-R pairs are diffusible. Run CellChat in a diffusion manner based on the `interaction.range`.\n")
        nLR1 <- nLR
      } else {
        cat("The input L-R pairs have both secreted signaling and contact-dependent signaling. Run CellChat in a contact-dependent manner for `Cell-Cell Contact` signaling, and in a diffusion manner based on the `interaction.range` for other L-R pairs. \n")
        nLR1 <- max(which(pairLR.use$annotation %in%
                            c("Secreted Signaling", "ECM-Receptor", "Non-protein Signaling")))
      }
    } else {
      cat("Run CellChat in a diffusion manner based on the `interaction.range` for all L-R pairs. Setting `contact.dependent = TRUE` if preferring a contact-dependent manner for `Cell-Cell Contact` signaling. \n")
      nLR1 <- nLR
    }
  } else {
    print(paste0(">>> Run CellChat on sc/snRNA-seq data <<< [", Sys.time(), "]"))
    d.spatial <- matrix(NaN, nrow = numCluster, ncol = numCluster)
    P.spatial <- matrix(1, nrow = numCluster, ncol = numCluster)
    adj.contact <- matrix(1, nrow = numCluster, ncol = numCluster)
    contact.dependent <- FALSE; contact.dependent.forced <- FALSE
    contact.range <- NULL; contact.knn.k <- NULL
    ## Upstream NULLs these five on the RNA branch, and they are recorded in
    ## `options$parameter`, so the NULLED values are the observable outcome rather than the
    ## caller's arguments.
    distance.use <- NULL; interaction.range <- NULL; ratio <- NULL; tol <- NULL; k.min <- NULL
  }

  ## `t(data.use)` -- the kernel consumes cells x genes, column-major, which is what the
  ## upstream `aggregate` call sees. `data.use <- data/max(data)` happens inside the
  ## kernel, because `max` must span the whole matrix including zeros (R-ism 5).
  ## Bound to locals before the call, not inline as arguments.
  ##
  ## R evaluates call arguments lazily, so `lr_agonist = lr_col(pairLR.use, "agonist", nLR)` is a
  ## promise that the callee forces at an arbitrary later point -- after the spatial branch has
  ## rebound `nLR1`, after `nLR` has been used for other things, and with no way to see from the
  ## kernel's error message ("lr_agonist has 1 entries but lr_ligand has 8") that the *caller*
  ## computed it. Binding them first makes the values eager, named, and inspectable, and it is how
  ## a kernel that validates its inputs deserves to be called.
  lr_ligand <- as.character(pairLR.use$ligand)
  lr_receptor <- as.character(pairLR.use$receptor)
  lr_agonist <- lr_col(pairLR.use, "agonist", nLR)
  lr_antagonist <- lr_col(pairLR.use, "antagonist", nLR)
  lr_co_a <- lr_col(pairLR.use, "co_A_receptor", nLR)
  lr_co_i <- lr_col(pairLR.use, "co_I_receptor", nLR)
  lr_label <- rownames(pairLR.use)
  ## A loud, R-level invariant. The kernel checks the same thing and refuses with
  ## "lr_agonist has N entries but lr_ligand has M", which names neither the fixture nor the
  ## caller; this says what was actually handed over, and it fires *before* the spatial messages
  ## have been printed so the failure is not mistaken for a kernel crash mid-run.
  lr_lens <- c(nrow(pairLR.use), length(lr_ligand), length(lr_receptor), length(lr_agonist),
               length(lr_antagonist), length(lr_co_a), length(lr_co_i), length(lr_label))
  if (any(lr_lens != lr_lens[1])) {
    stop(sprintf(paste("LR vector lengths must all be nrow(pairLR.use) = %d, got",
                       "nrow=%d ligand=%d receptor=%d agonist=%d antagonist=%d co_A=%d co_I=%d",
                       "label=%d"), lr_lens[1], lr_lens), call. = FALSE)
  }

  ## Upstream calls set.seed(seed.use) and draws all bootstrap permutations before its LR loop.
  ## Rust uses its own exact MT19937 implementation to compute those same permutations, so replay
  ## the R draws here to preserve the caller-visible .Random.seed without allocating a second
  ## nC-by-nboot permutation matrix.
  res <- compute_commun_prob(
    data = as.numeric(t(data)),
    dim = c(nrow(data), ncol(data)),
    genes = rownames(data),
    group = as.integer(group) - 1L,
    group_levels = levels(group),
    lr_ligand = lr_ligand,
    lr_receptor = lr_receptor,
    ## Every LR vector must be `nLR` long: the kernel indexes them positionally against
    ## `lr_label`, and a zero-length one is a panic, not a missing-optional. `LR.use` is a
    ## user-supplied subset of `LRsig` and need not carry all six metadata columns, so an absent
    ## column becomes `NA` rather than a shorter vector. `as.character(NULL)` is `character(0)`,
    ## which is how this used to end in "lr_agonist has 0 entries but lr_ligand has 1" -- a Rust
    ## panic surfacing as an R error, from a path the permanently-engaged fallback had been
    ## hiding.
    lr_agonist = lr_agonist,
    lr_antagonist = lr_antagonist,
    lr_co_a = lr_co_a,
    lr_co_i = lr_co_i,
    lr_label = lr_label,
    complex_names = rownames(complex_input) %||% character(0),
    complex_subunits = as.character(as.matrix(complex_input[, complex_cols, drop = FALSE])),
    complex_n_cols = length(complex_cols),
    cofactor_names = rownames(cofactor_input) %||% character(0),
    cofactor_subunits = as.character(as.matrix(cofactor_input[, cofactor_cols, drop = FALSE])),
    cofactor_n_cols = length(cofactor_cols),
    nboot = as.integer(nboot),
    seed = as.integer(seed.use),
    kh = Kh,
    n = n,
    mean_type = type,
    trim = trim,
    population_size = population.size,
    p_spatial = as.numeric(P.spatial),
    adj_contact = as.numeric(adj.contact),
    n_lr1 = as.integer(nLR1),
    advance_rng = function() cellchatrs_advance_rng(seed.use, nboot, nC)
  )

  Prob <- array(res$prob, dim = res$dim)
  Pval <- array(res$pval, dim = res$dim)
  dimnames(Prob) <- list(levels(group), levels(group), rownames(pairLR.use))
  dimnames(Pval) <- dimnames(Prob)
  object@net <- list("prob" = Prob, "pval" = Pval)

  ## `options$parameter` is part of the drop-in object and must match upstream field for
  ## field, including the RNA branch's NULLED spatial parameters. Upstream's RNA branch
  ## sets `distance.use <- interaction.range <- ratio <- tol <- k.min <- NULL`.
  object@options$run.time <- as.numeric(Sys.time() - ptm, units = "secs")
  ## Recorded from the values *in force*, not from constants. Upstream's RNA branch NULLs
  ## `distance.use`, `interaction.range`, `ratio`, `tol` and `k.min` and forces
  ## `contact.dependent`/`contact.dependent.forced` to `FALSE` before this list is built, so the
  ## RNA object records NULLs; the spatial branch leaves them alone, so it records what the caller
  ## passed. Writing the RNA constants unconditionally makes every spatial run report
  ## `distance.use = NULL` while having just used `TRUE`, and `parity.json` would call it a match.
  object@options$parameter <- list(
    type.mean = type, trim = trim, raw.use = raw.use,
    population.size = population.size, nboot = nboot, seed.use = seed.use,
    Kh = Kh, n = n,
    distance.use = distance.use, interaction.range = interaction.range,
    ratio = ratio, tol = tol, k.min = k.min,
    contact.dependent = contact.dependent, contact.range = contact.range,
    contact.knn.k = contact.knn.k,
    contact.dependent.forced = contact.dependent.forced
  )

  print(paste0(">>> CellChat inference is done. Parameter values are stored in `object@options$parameter` <<< [",
               Sys.time(), "]"))
  object
}

cellchatrs_advance_rng <- function(seed, nboot, n_cells) {
  base::set.seed(seed)
  for (i in seq_len(nboot)) sample.int(n_cells, size = n_cells)
  invisible(NULL)
}

## One LR metadata column as a character vector of exactly `n` entries.
##
## Upstream's `pairLRsig` always carries all six because `subsetDB` produces them, so this only
## matters for a hand-built or column-subset `LR.use`. Returning `NA` rather than a zero-length
## vector keeps every LR vector the same length, which is what the kernel requires; the `NA`s are
## then handled by the same `which(!is.na(...) & ... != "")` tests upstream uses.
lr_col <- function(x, nm, n) {
  v <- if (nm %in% names(x)) as.character(x[[nm]]) else rep("", n)
  ## `NA` becomes `""`, which is exactly equivalent: every use upstream makes of these columns is
  ## `which(!is.na(v) & v != "")` or `v == ""`, and both reject the two spellings identically.
  ## `NA_character_` also cannot cross into a Rust `Vec<String>` -- `extendr` turns it into a
  ## `Must not be NA` failure -- so passing the `NA` through is not even available as a choice.
  v[is.na(v)] <- ""
  ## The kernel indexes every LR vector positionally against `lr_label` and checks the lengths
  ## first, so a column that is present but the wrong length is a Rust panic. Upstream would index
  ## out of bounds instead; neither is a useful answer, and normalising to `n` keeps the failure a
  ## plain R one.
  if (length(v) != n) v <- rep("", n)
  v
}


## `LR.use` handling from upstream: reorder by the annotation factor, then back to
## character. R-ism 9 territory; kept here because it is pure data-frame bookkeeping and
## not on the hot path.
.order_lr_use <- function(LR.use) {
  if (length(unique(LR.use$annotation)) > 1) {
    LR.use$annotation <- factor(LR.use$annotation,
                                levels = c("Secreted Signaling", "ECM-Receptor",
                                           "Non-protein Signaling", "Cell-Cell Contact"))
    LR.use <- LR.use[order(LR.use$annotation), , drop = FALSE]
    LR.use$annotation <- as.character(LR.use$annotation)
  }
  LR.use
}

## ----------------------------------------------------------------------------
## identifyOverExpressedGenes
## ----------------------------------------------------------------------------

#' `identifyOverExpressedGenes`, backed by the Rust kernel.
#'
#' Routes to Rust only where the numerics are a *pure function of the input matrix*:
#' the `do.DE = FALSE` branch (percentage and `min.cells` only) and the
#' `do.DE = TRUE, do.fast = FALSE` branch (the Wilcoxon path), with or without
#' `group.dataset` / `group.DE.combined`. Everything else -- `do.fast = TRUE` (presto, an
#' optional C++ dependency that may not be installed, and a *different algorithm* rather than a
#' faster one) -- goes to the pinned upstream body.
#'
#' `group.dataset` was long filed as blocked on presto and is not: upstream honours it in **both**
#' halves of `if (do.fast)`, at `utilities.R:429-484` for presto and `:512-520` for the Wilcoxon
#' test. Only the first needs presto. What it changes is the *cell selection* --
#'
#' ```r
#' cell.use1 <- which((labels == level.use[i]) & (labels.dataset == pos.dataset))
#' cell.use2 <- which((labels == level.use[i]) & (labels.dataset != pos.dataset))
#' ```
#'
#' or, with `group.DE.combined = TRUE`, both pooled across groups with no `labels` term at all --
#' so the kernel takes the selection as a parameter instead of growing a flag through the path.
#'
#' The argument resolution (`group.by` / `idents.use` / `invert` / `features.use`) and the
#' two `stop()` calls stay in R, so the error messages and the `message()` for a non-factor
#' `group.by` are upstream's. `features.use <- intersect(features, row.names(X))` is
#' R's `intersect`, which keeps the *first* argument's order -- so a caller-supplied
#' `features` vector reorders the output relative to the rownames, and that has to happen
#' in R, not in Rust.
#'
#' ## The schema collapses when nothing passes
#'
#' Upstream builds `markers.all <- data.frame()` and only `rbind`s a group's rows once
#' they survive `pvalues < thresh.p`. If no group does, the frame stays zero-column, and
#' the next line, `markers.all$features <- as.character(markers.all$features)`, *adds* a
#' single zero-length column: the "marker table" is `0x1` with only `features`, not `0x7`.
#' The kernel therefore returns `n_before_only_pos` so the shim knows which of the two
#' shapes to build; see `r_core::wilcox::GroupMarkers`.
#' @export
identifyOverExpressedGenes <- function(object, data.use = NULL, group.by = NULL, idents.use = NULL,
                                      invert = FALSE, group.dataset = NULL,
                                      pos.dataset = NULL, group.DE.combined = FALSE,
                                      features.name = "features",
                                      only.pos = TRUE, features = NULL,
                                      return.object = TRUE, thresh.pc = 0, thresh.fc = 0,
                                      thresh.p = 0.05, do.DE = TRUE, do.fast = TRUE,
                                      min.cells = 10) {
  if (!is.list(object@var.features)) {
    stop("Please update your CellChat object via `updateCellChat()`")
  }
  input_data_use <- data.use
  if (is.null(data.use)) {
    X <- object@data.signaling
    if (nrow(X) < 3) {
      stop("Please check `object@data.signaling` and ensure that you have run `subsetData` and that the data matrix `object@data.signaling` looks OK.")
    }
  } else {
    X <- data.use
  }
  if (is.null(features)) {
    features.use <- row.names(X)
  } else {
    features.use <- intersect(features, row.names(X))
  }
  data.use <- X[features.use, ]

  if (do.fast || !do.DE || is.null(dim(data.use))) {
    ## `do.fast = TRUE` needs presto; upstream checks with `rlang::is_installed` and
    ## emits a multi-line `message()` before stopping, so let it run.
    if (!do.DE && !do.fast && !is.null(dim(data.use))) {
      ## The `else` branch, which is reachable without presto and without the Wilcoxon:
      ## ```r
      ## markers.all <- data.frame(features = as.character(rownames(data.use)),
      ##                           nCells = rowSums(data.use > 0))
      ## markers.all <- dplyr::filter(markers.all, nCells >= min.cells)
      ## ```
      ## `rowSums(data.use > 0)` is `dplyr::filter`'s predicate plus a count, so it goes to
      ## the kernel -- which also sidesteps a real hazard: `data.use` is a `dgCMatrix` here,
      ## so the expression is `rowSums` of an `lgCMatrix`, an S4 method in `Matrix` that is
      ## not reachable from this package's namespace in a bare session. Upstream gets it
      ## because its `env` chains to the global search path.
      res <- expressed_in_cells(
        data = as.numeric(data.use),
        features = as.character(row.names(data.use)),
        min_cells = as.integer(min.cells)
      )
      ## Matrix::rowSums returns a named double vector. Let data.frame adopt its row names
      ## during construction to preserve both storage type and serialized attribute order.
      ##
      ## Storage type is input-dependent, and `identical()` sees it: `base::rowSums` on a
      ## dense logical matrix returns double, while `Matrix::rowSums` on an `lgCMatrix`
      ## returns integer. The kernel counts in integers (`res$n_cells` is integer), so for
      ## dense input coerce to double and for sparse input keep integer -- matching what
      ## `rowSums(data.use > 0)` itself returns in each case. Getting this backwards fails
      ## `identical()` while `all.equal()` stays silent, which is how the `oeg_node`
      ## configuration (sparse `data.use`, 29 retained genes) caught it after the dense
      ## `test-contract-inputs.R` case had already passed.
      counts <- stats::setNames(res$n_cells, res$features)
      if (!inherits(data.use, "sparseMatrix")) {
        counts <- stats::setNames(as.numeric(counts), res$features)
      }
      markers.all <- data.frame(features = res$features, nCells = counts,
                                stringsAsFactors = FALSE)
      ## Upstream's next statement is `dplyr::filter(markers.all, nCells >= min.cells)`, and
      ## the predicate is the only part of it that reached Rust -- so the frame itself is
      ## built here instead. `dplyr::filter` rebuilds the frame through vctrs and leaves its
      ## attribute pairlist in the order `row.names, names, class`; `data.frame()` produces
      ## `names, class, row.names`.
      ##
      ## That difference is invisible to `identical()`, which looks attributes up by name and
      ## so calls the two frames equal, and visible to `serialize()`, which writes the pairlist
      ## in stored order. It is therefore observable in any `saveRDS` of a CellChat object, and
      ## it is what `tests/test-contract-inputs.R` compares. Reorder rather than rebuild:
      ## `c("row.names", setdiff(...))` keeps whatever attributes a future R version adds,
      ## where naming all three would silently drop them.
      at <- attributes(markers.all)
      attributes(markers.all) <- at[c("row.names", setdiff(names(at), "row.names"))]
      object@var.features[[features.name]] <- markers.all$features
      object@var.features[[paste0(features.name, ".info")]] <- markers.all
      return(if (return.object) object else markers.all)
    }
    ## Named, not positional. Upstream's signature is
    ## `(object, group.by, idents.use, invert, features.name, only.pos, features,
    ##   return.object, thresh.pc, thresh.fc, thresh.p, do.DE, do.fast, min.cells)`
    ## and this function's is the same, so a positional call happens to line up -- but a
    ## one-argument reordering in either file would then silently feed `invert` a string
    ## and surface as `argument is not interpretable as logical` deep inside upstream's
    ## `do.DE` branch, which is exactly what happened while this was being written.
    return(cellchatrs_upstream_identifyOverExpressedGenes(
      object = object, data.use = input_data_use, group.by = group.by, idents.use = idents.use, invert = invert,
      group.dataset = group.dataset, pos.dataset = pos.dataset,
      group.DE.combined = group.DE.combined,
      features.name = features.name, only.pos = only.pos, features = features,
      return.object = return.object, thresh.pc = thresh.pc, thresh.fc = thresh.fc,
      thresh.p = thresh.p, do.DE = do.DE, do.fast = do.fast, min.cells = min.cells))
  }

  ## `do.DE = TRUE, do.fast = FALSE`.
  data.use <- as.matrix(data.use)
  if (is.null(group.by)) {
    labels <- object@idents
    if (!is.factor(labels)) {
      message("Use the joint cell labels from the merged CellChat object")
      labels <- object@idents$joint
    }
  } else {
    labels <- object@meta[[group.by]]
  }
  if (!is.factor(labels)) labels <- factor(labels)
  level.use <- levels(labels)[levels(labels) %in% unique(labels)]
  if (!is.null(idents.use)) {
    level.use <- if (invert) level.use[!(level.use %in% idents.use)] else
      level.use[level.use %in% idents.use]
  }
  if (length(level.use) < 1L) {
    markers.all <- data.frame()
    object@var.features[[features.name]] <- character(0)
    object@var.features[[paste0(features.name, ".info")]] <-
      data.frame(features = character(0), stringsAsFactors = FALSE)
    return(if (return.object) object else object@var.features[[paste0(features.name, ".info")]])
  }
  all_levels <- levels(labels)

  ## Upstream's `labels.dataset`, built before the group loop (utilities.R:420-426):
  ## ```r
  ## labels.dataset <- as.character(object@meta[[group.dataset]])
  ## if (!(pos.dataset %in% unique(labels.dataset))) {
  ##   cat("Please set pos.dataset to be one of the following dataset names: ", unique(as.character(labels.dataset)))
  ##   stop()
  ## }
  ## labels.dataset[labels.dataset != pos.dataset] <- toString(setdiff(unique(labels.dataset), pos.dataset))
  ## labels.dataset <- factor(labels.dataset, levels = c(pos.dataset, setdiff(unique(labels.dataset), pos.dataset)))
  ## ```
  ## Two details are load-bearing. The bare `stop()` has **no message**, so R's error text is empty
  ## and any message here would be a divergence; and the collapse means a three-dataset
  ## `group.dataset` still has exactly **two** levels, with the other two joined by `toString`'s
  ## ", " into one. `toString` is what puts the comma-space in the level name.
  dataset_levels <- NULL
  if (!is.null(group.dataset)) {
    labels.dataset <- as.character(object@meta[[group.dataset]])
    if (!(pos.dataset %in% unique(labels.dataset))) {
      cat("Please set pos.dataset to be one of the following dataset names: ",
          unique(as.character(labels.dataset)))
      stop()
    }
    labels.dataset[labels.dataset != pos.dataset] <-
      toString(setdiff(unique(labels.dataset), pos.dataset))
    dataset_levels <- c(pos.dataset, setdiff(unique(labels.dataset), pos.dataset))
    labels.dataset <- factor(labels.dataset, levels = dataset_levels)
  }

  ## NB the kernel's argument names are the Rust identifiers (`thresh_pc`), which R does
  ## **not** partial-match to this function's `thresh.pc`. Passing `thresh_pc = thresh.pc`
  ## would silently bind the wrong thing -- or fail to resolve at all -- so the values are
  ## bound to locals first.
  .pc <- thresh.pc
  .fc <- thresh.fc
  .p <- thresh.p
  res <- if (is.null(group.dataset)) {
    identify_over_expressed_genes(
      data = as.numeric(data.use),
      features = as.character(row.names(data.use)),
      ## 0-based, like `compute_ave_expr`'s `group`.
      group_index = as.integer(labels) - 1L,
      group_levels = all_levels,
      thresh_pc = .pc,
      thresh_fc = .fc,
      thresh_p = .p,
      only_pos = only.pos,
      ## `nrow(X)`, not `nrow(data.use)`: upstream writes `n = nrow(X)`.
      n_adjust = nrow(X)
    )
  } else {
    ## The only difference `group.dataset` makes to the Wilcoxon branch is how
    ## `cell.use1` / `cell.use2` are chosen per group; the percentage filter, `mean.fxn`,
    ## `wilcox.test` and the Bonferroni multiplier are untouched. The kernel takes the
    ## selection as a parameter rather than growing a `group.dataset` flag through the path.
    identify_over_expressed_genes_dataset(
      data = as.numeric(data.use),
      features = as.character(row.names(data.use)),
      group_index = as.integer(labels) - 1L,
      group_levels = all_levels,
      ## 0 == the positive dataset. Upstream has already collapsed every other dataset into
      ## one level, so a single bit is all the selection needs.
      dataset_index = as.integer(labels.dataset) - 1L,
      thresh_pc = .pc,
      thresh_fc = .fc,
      thresh_p = .p,
      only_pos = only.pos,
      n_adjust = nrow(X),
      group_de_combined = isTRUE(group.DE.combined)
    )
  }
  ## groups it excluded. (The kernel iterates every level of `labels`; filtering here keeps
  ## `n_before_only_pos` an over-estimate in exactly one case -- every excluded group would
  ## have contributed a marker and no included group did -- which only affects the schema,
  ## and only when `n_before_only_pos` is already positive.)
  if (length(level.use) < length(all_levels)) {
    sel <- res$clusters %in% level.use
    res$clusters <- res$clusters[sel]
    for (cc in c("features", "pvalues", "logFC", "pct1", "pct2", "pvalues_adj")) {
      res[[cc]] <- res[[cc]][sel]
    }
    res$nrow <- length(res$features)
  }
  ## `length(res$features)`, not `res$n_before_only_pos`. The counter is an over-estimate by
  ## design once `level.use` has filtered rows away, and building frames from an *empty* result
  ## is worse than a wrong count: `starts <- which(c(TRUE, <logical(0)>))` is `1` and
  ## `ends <- c(numeric(0), 0)` is `0`, so the frame loop indexes `res$features[c(1, 0)]` and
  ## `data.frame` answers `row names contain missing values`. Upstream's guard is per group --
  ## `if (nrow(gde) > 0) markers.all <- rbind(markers.all, gde)` -- so an empty result leaves
  ## `markers.all` as the bare `data.frame()` it started as.
  if (length(res$features) > 0L) {
    ## Row names come from upstream's
    ## `data.frame(clusters, features, pvalues, logFC, pct.1, pct.2, pvalues.adj,
    ##             data.alpha[features, , drop = FALSE], pvalues.adj)`: the unnamed
    ## `data.alpha[features, ]` **matrix** argument carries `rownames == features`, and
    ## `data.frame()` adopts them. So the marker table's row names are the feature names --
    ## which is why upstream's `features.info` prints with `features` down the side.
    ## Rebuilding per group and `rbind`-ing reproduces R's duplicate-name uniquification
    ## ("G14", "G14.1") for a feature that is a marker in more than one group.
    frames <- list()
    g <- res$clusters
    starts <- which(c(TRUE, g[-1L] != g[-length(g)]))
    ends <- c(starts[-1L] - 1L, length(g))
    for (k in seq_along(starts)) {
      ii <- starts[k]:ends[k]
      frames[[k]] <- data.frame(
        clusters = g[ii], features = res$features[ii], pvalues = res$pvalues[ii],
        logFC = res$logFC[ii], pct.1 = res$pct1[ii], pct.2 = res$pct2[ii],
        pvalues.adj = res$pvalues_adj[ii], stringsAsFactors = FALSE,
        row.names = res$features[ii]
      )
    }
    ## Upstream's loop, not `do.call(rbind, frames)`: `markers.all` starts as
    ## `data.frame()` and is rbind-ed **one group at a time**, and the difference is observable.
    ## R's row-name uniquification renames a collision to `x.1` only while `x.1` is free, and
    ## falls back to appending a bare digit (`x` -> `x1`) once it is taken -- so the *number* of
    ## earlier rbind passes changes the names. `do.call(rbind, ...)` uniquifies the whole
    ## concatenation once and produces a different set.
    markers.all <- data.frame()
    for (k in seq_along(frames)) {
      markers.all <- rbind(markers.all, frames[[k]])
    }
  } else if (is.null(group.dataset)) {
    ## Upstream's collapsed frame: `data.frame()` plus the one `features` column.
    markers.all <- data.frame(features = character(0), stringsAsFactors = FALSE)
  } else {
    ## ... but only that. With `group.dataset` set, upstream is still holding a **0 x 0**
    ## `data.frame()` at this point, because `markers.all <- data.frame()` and the `rbind`
    ## is conditional on a group having rows. The next statement adds a `datasets` column to
    ## it, and the one after orders by `-markers.all$logFC`, which does not exist:
    ##
    ## ```text
    ## Error in -markers.all$logFC : invalid argument to unary operator
    ## ```
    ##
    ## Substituting the 1-column frame here would return a zero-row marker table where
    ## upstream refuses to return anything, so the shape is part of the contract.
    markers.all <- data.frame()
  }
  if (!is.null(group.dataset)) {
    ## utilities.R:594-598, in order and including the parts that look redundant:
    ## ```r
    ## markers.all$datasets[markers.all$logFC > 0] <- pos.dataset
    ## markers.all$datasets[markers.all$logFC < 0] <- setdiff(unique(labels.dataset), pos.dataset)
    ## markers.all$datasets <- factor(markers.all$datasets, levels = levels(labels.dataset))
    ## markers.all <- markers.all[order(markers.all$datasets, markers.all$pvalues, -markers.all$logFC), ]
    ## ```
    ## The first assignment *creates* the column, as `NA` everywhere it does not match. With
    ## `only.pos = TRUE` the `logFC > 0` filter has already removed the zero rows, so nothing
    ## is left as `NA`; with `only.pos = FALSE` a `logFC == 0` row keeps `NA`, and the
    ## subsequent `factor()` turns that into a real level `NA` rather than dropping the row.
    ## `setdiff` on a factor returns a plain character here -- confirmed, not assumed, because
    ## assigning a factor into a character vector is a silent-`NA` trap.
    markers.all$datasets[markers.all$logFC > 0] <- pos.dataset
    markers.all$datasets[markers.all$logFC < 0] <- setdiff(unique(labels.dataset), pos.dataset)
    markers.all$datasets <- factor(markers.all$datasets, levels = levels(labels.dataset))
    ## The reorder keeps the row names, and they are *not* the feature names after `rbind`:
    ## a feature that is a marker in three groups gets `gene7`, `gene7.1`, `gene7.2` from R's
    ## duplicate-name uniquification, and the feature column still says `gene7` three times.
    ## Resetting them to `NULL` would drop exactly the information that distinguishes the rows.
    markers.all <- markers.all[order(markers.all$datasets, markers.all$pvalues,
                                      -markers.all$logFC), ]
  }
  markers.all$features <- as.character(markers.all$features)
  object@var.features[[features.name]] <- markers.all$features
  object@var.features[[paste0(features.name, ".info")]] <- markers.all
  if (return.object) object else markers.all
}

## ----------------------------------------------------------------------------
## computeCommunProbPathway
## ----------------------------------------------------------------------------

#' The pinned upstream `computeCommunProbPathway`, verbatim.
#' @export
cellchatrs_upstream_computeCommunProbPathway <- function(...) {
  get("computeCommunProbPathway", envir = cellchatrs_upstream_cached())(...)
}

#' `computeCommunProbPathway`, backed by the Rust kernel.
#'
#' The aggregation, the LONG_DOUBLE sums, the `!= 0` selection and the decreasing-total
#' reordering are all in `r-core`; only the two `return.object` shapes and the S4 slot
#' writes are here.
#'
#' ## One upstream failure is preserved, not fixed
#'
#' ```r
#' prob.pathways <- aperm(apply(prob, c(1, 2), by, group, sum), c(2, 3, 1))
#' ```
#'
#' `apply(prob, c(1, 2), by, group, sum)` returns a **2-d matrix** when `group` has a single
#' level and a **3-d array** otherwise, so `aperm(x, c(2, 3, 1))` raises
#' `'perm' is of wrong length 3 (!= 2)` for a one-pathway `LRsig`. That is reachable from
#' ordinary use -- an `LRsig` where every interaction belongs to the same pathway -- and a
#' drop-in replacement has to fail identically, so this case routes to upstream. The Rust
#' core computes a value for it; it is the shim that must not claim it. Pinned by
#' `src/rust/crates/r-core/tests/pathway_parity.rs::single_pathway_errors_like_upstream`.
#' @export
computeCommunProbPathway <- function(object = NULL, net = NULL, pairLR.use = NULL, thresh = 0.05) {
  if (is.null(net)) net <- object@net
  if (is.null(pairLR.use)) pairLR.use <- object@LR$LRsig
  pathway_name <- as.character(pairLR.use$pathway_name)
  if (length(unique(pathway_name)) < 2L) {
    return(cellchatrs_upstream_computeCommunProbPathway(object, net, pairLR.use, thresh))
  }
  k <- dim(net$prob)[1]
  lr <- dimnames(net$prob)[[3]]
  res <- compute_commun_prob_pathway(
    prob = as.numeric(net$prob),
    pval = as.numeric(net$pval),
    group_levels = dimnames(net$prob)[[1]],
    interaction_names = if (is.null(lr)) character(0) else lr,
    pathway_name = pathway_name,
    thresh = thresh
  )
  prob.pwp <- array(res$prob, dim = res$dim,
                    dimnames = list(dimnames(net$prob)[[1]], dimnames(net$prob)[[2]],
                                    res$pathways))
  if (is.null(object)) {
    list(pathways = res$pathways, prob = prob.pwp)
  } else {
    object@net$LRs <- res$lrs
    object@netP$pathways <- res$pathways
    object@netP$prob <- prob.pwp
    object
  }
}

## ----------------------------------------------------------------------------
## filterCommunication
## ----------------------------------------------------------------------------

#' `scales::percent(x, accuracy = .1)`, inlined.
#'
#' `scales` is a heavy dependency for one formatter, and the exact string is part of the
#' contract -- `filterCommunication` prints it. For finite `x` this is
#' `paste0(format(round(100 * x, 1), nsmall = 1, scientific = FALSE), "%")`.
#' For missing input (`NA` or `NaN`, e.g. `(0 - 0) / 0` on an all-zero network) it returns
#' the bare string `"NA"`, because `scales::percent(NA)` is `NA` and the caller pastes it --
#' `paste0(NA, " interactions...")` renders the two characters "NA".
#' @noRd
cellchatrs_percent <- function(x, accuracy = 0.1) {
  v <- 100 * x
  digits <- as.integer(log10(1 / accuracy))
  ## `scales::percent(NA)` is `NA` (not "NA%", not "NaN%"): the caller pastes it with
  ## `paste0(pct, " interactions are removed!")`, and `paste0` renders `NA` as the two characters
  ## "NA". So this returns the bare string "NA" and lets the caller supply the rest -- returning
  ## "NA%" here would print "NA%% ...", and the old code returned "NaN%", which matched neither
  ## `NA` (from `NA_real_` input) nor `NA` (from `NaN` input; `scales::percent(NaN)` is also `NA`).
  ## Both missingness flavors therefore map to "NA", exactly as upstream's call chain does.
  if (length(v) == 0L) {
    return(character(0))
  }
  if (is.na(v)) {
    return("NA")
  }
  paste0(format(round(v, digits), nsmall = digits, scientific = FALSE, trim = TRUE), "%")
}

#' `filterCommunication`, backed by the Rust kernel.
#'
#' The numerics -- the `min.cells` zeroing, the per-sample group means, the outer products,
#' the binarisation and the cross-sample consistency mask -- are in `r-core`. The `cat()`
#' output stays here because it is `scales::percent` and `toString(levels(...))`, i.e. R
#' formatting, and because *where* the messages appear is part of the observable behaviour.
#'
#' ## Two upstream behaviours that are reproduced rather than repaired
#'
#' * **`nonFilter.keep = TRUE` is a no-op.** `net <- object@net` is bound *before* the two
#'   `object@net$prob.nonFilter <- ...` assignments, and the function ends with
#'   `object@net <- net`, which throws them away. The `cat()` still fires, so the flag has an
#'   observable effect and no other one. See `src/rust/crates/r-core/tests/filter_parity.rs`.
#' * **An all-zero `net$prob` with `min.samples >= 2` is `"subscript out of bounds"`**, from
#'   `for (jj in 1:length(LR.nonzero))` with `LR.nonzero` empty being `1:0` and then
#'   `score.LR[, , 0, i] <- ...`. The kernel returns that as an error carrying upstream's
#'   text, so it is re-raised unchanged.
#'
#' The `raw.use` setting is read back out of `options$parameter`, not taken as an argument, so a
#' `filterCommunication` on a projected object filters the *projected* data. The corpus exercises
#' both settings, because the two take visibly different paths through the kernel: the per-sample
#' means come from a different matrix, so the zero-fill and the consistency mask can differ even
#' when the network structure is identical.
#' @export
filterCommunication <- function(object, min.cells = 10, min.samples = NULL, rare.keep = FALSE,
                                nonFilter.keep = FALSE) {
  net <- object@net
  raw_use <- isTRUE(object@options$parameter$raw.use)
  if (!raw_use) {
    ## Upstream writes `"data.smooth" %in% methods::slotNames(object) == FALSE` -- with **no**
    ## `!`. `%in%` binds tighter than `==`, so this is `(in_set) == FALSE`, i.e. "the slot is
    ## *absent*", and the stop fires exactly when it should. An earlier version of this shim
    ## added a `!` on the reading that the condition was doubly negated, which inverted it: every
    ## `raw.use = FALSE` call on a *valid* object then failed with the "is missing" message, and
    ## the corpus could not tell the difference because it only exercised the broken direction.
    if ("data.smooth" %in% methods::slotNames(object) == FALSE) {
      stop("`object@data.smooth` is missing. Please update the CellChat object via `updateCellChat`! \n")
    }
    return(cellchatrs_upstream_filterCommunication(object, min.cells, min.samples, rare.keep,
                                                   nonFilter.keep))
  }
  idents <- object@idents
  ## A missing `meta$samples` means one sample, not zero: upstream reads
  ## `sample.id <- levels(sample.info)` (NULL here) and its `length(sample.id) >= 2` guard then
  ## skips the whole multi-sample branch. The kernel needs a per-cell index vector it can actually
  ## receive -- `Vec<i32>` rejects NULL -- so a single `sample1` level stands in, and with
  ## `n_samples == 1` the kernel's own `n_samples >= 2` guard skips identically. No warning:
  ## upstream's filter prints none for this (the "assumed to belong to `sample1`" warning belongs
  ## to `createCellChat`, which already printed it if it applied).
  sample.info <- object@meta$samples
  if (is.null(sample.info)) {
    sample.info <- factor(rep("sample1", length(idents)))
  }
  if (is.null(min.samples)) min.samples <- 1L
  if (nonFilter.keep) {
    ## Printed, then discarded -- see the header. Reproduced because the message is real
    ## output and the slots genuinely never appear.
    cat("The non-filtered cell-cell communication is stored in `object@net$prob.nonFilter` and `object@net$pval.nonFilter`. \n")
  }
  lr <- dimnames(net$prob)[[3]]
  DB <- object@DB
  cx <- DB$complex
  subunit_cols <- grep("^subunit", colnames(cx), value = TRUE)
  ## An empty complex table (0 columns, as in a minimal or synthetic DB) gives NULL rownames
  ## and an empty subunit selection, both of which the `Vec<String>` conversion rejects with
  ## "Expected Strings got Null". Upstream never notices: its per-sample path looks complexes
  ## up by name and finds none. Coerce to empty vectors here so the kernel sees "no complexes",
  ## which is what an empty table means, rather than failing at the boundary.
  cx_names <- rownames(cx) %||% character(0)
  cx_cells <- if (length(subunit_cols)) as.character(t(cx[, subunit_cols, drop = FALSE])) else
    character(0)
  res <- filter_communication(
    prob = as.numeric(net$prob),
    group_levels = levels(idents),
    interaction_names = lr,
    group_index = as.integer(idents) - 1L,
    sample_levels = levels(sample.info),
    sample_index = as.integer(sample.info) - 1L,
    data = as.numeric(as.matrix(object@data.signaling)),
    gene_names = rownames(as.matrix(object@data.signaling)),
    ligand = as.character(DB$interaction$ligand[match(lr, DB$interaction$interaction_name)]),
    receptor = as.character(DB$interaction$receptor[match(lr, DB$interaction$interaction_name)]),
    complex_names = cx_names,
    ## `t()` first: the binding reads `complex_subunits[r * n_subunit_cols + c]`, i.e.
    ## **row-major** over the complexes, and `as.character()` on a matrix is column-major.
    ## Passing it straight through transposes the subunit table, so every complex gets the
    ## next one's subunits -- and `computeExpr_complex`'s geometric mean is then over the
    ## wrong genes. The consistency mask then keeps a different set of pairs and the net
    ## comes out with 21 significant entries where upstream has 5.
    complex_subunits = cx_cells,
    complex_n_subunit_cols = length(subunit_cols),
    symbols = as.character(DB$geneInfo$Symbol) %||% character(0),
    min_cells = as.integer(min.cells),
    min_samples = as.integer(min.samples),
    rare_keep = rare.keep,
    mean_type = object@options$parameter$type.mean,
    trim = if (is.null(object@options$parameter$trim)) 0 else object@options$parameter$trim
  )

  ## The message order is upstream's and is not the order of the computation:
  ##   1. `nonFilter.keep`          (before anything)
  ##   2. the `min.cells` line      (a `cat` with a `'\t'` separator, so the percentage lands
  ##                                 on the *same* line -- the corpus records the tab)
  ##   3. the sample-count `stop()`  <- here, not before (2)
  ##   4. one line per sample with too few cells, in sample order
  ##   5. the cross-sample percentage
  if (length(res$cell_excludes) > 0L) {
    cat("The cell-cell communication related with the following cell groups are excluded due to the few number of cells: ",
        toString(levels(idents)[res$cell_excludes + 1L]), "!", "\t")
    cat(paste0(cellchatrs_percent((res$n_interaction0 - res$n_interaction1) / res$n_interaction0,
                                  accuracy = 0.1), " interactions are removed!", "\n"))
  }
  if (!isTRUE(res$min_samples_ok)) {
    stop(paste0("There are only ", res$n_samples,
                " samples in the data. Please change the value of `min.samples`! "))
  }
  ## Index arithmetic rather than `split()`. `split(x, f)` factorises `f`, so the groups
  ## come out in *sorted* factor order rather than the order the cuts appear, and with a
  ## flat integer vector of group indices that silently reorders the group names in the
  ## message: upstream prints "g1, g3" and a `split()`-based shim printed "g3, g1".
  at <- 0L
  for (i in seq_along(levels(sample.info))) {
    len <- res$sample_excluded_len[i]
    g <- if (len > 0L) res$sample_excluded[(at + 1L):(at + len)] else integer(0)
    at <- at + len
    if (length(g) > 0L) {
      cat(paste0("The number of cells of the following cell groups in ", levels(sample.info)[i],
                 " sample are less than ", min.cells, " cells: ",
                 toString(levels(idents)[g + 1L]), "!", "\n"))
    }
  }
  if (isTRUE(res$ran_sample_filter)) {
    cat(paste0(cellchatrs_percent((res$n_interaction1 - res$n_interaction2) / res$n_interaction1,
                                  accuracy = 0.1),
               " interactions are removed due to their inconsistence across ", min.samples,
               " samples!", "\n"))
  }
  k <- length(levels(idents))
  nlr <- length(lr)
  ## `array()`, not `matrix()`: `net$prob` is a 3-d array and a `k x (k*nLR)` matrix has the
  ## wrong *number of dimnames*, so `dimnames=` rejects it. The flat order is R's
  ## column-major with the L-R axis fastest, which is exactly what `as.numeric()` delivered.
  net$prob <- array(res$prob, dim = c(k, k, nlr), dimnames = dimnames(net$prob))
  object@net <- net
  object
}

## ---------------------------------------------------------------------------------------
## `computeCellDistance`, backed by Rust.
## ---------------------------------------------------------------------------------------

#' The pinned upstream `computeCellDistance`, verbatim.
#' @export
cellchatrs_upstream_computeCellDistance <- function(...) {
  get("computeCellDistance", envir = cellchatrs_upstream_cached())(...)
}

#' @rdname cellchatrs_upstream_computeCellDistance
#'
#' `computeCellDistance` from pinned upstream, with `collapse::fdist` replaced by the Rust kernel.
#'
#' Returns a `dist` object, as upstream does -- `collapse::fdist` does not return a matrix, and the
#' rewrap reproduces the three places `stats::as.dist` disagrees with it (no `call` attribute,
#' `method = "euclidean"` present, `Labels` `NULL` rather than synthesised). Also reproduces the
#' two-column check, `fdist`'s own "at least 2 rows" error, `ratio` being applied whenever it is
#' non-`NULL`, and the threshold that needs **both** `interaction.range` and `tol` with `tol` as an
#' addend. See `docs/SEMANTICS.md`.
#' @export
computeCellDistance <- function(coordinates, interaction.range = NULL, ratio = NULL,
                                tol = NULL) {
  ## The two-column check comes first and is upstream's own message, so it is written here rather
  ## than produced by the kernel. It runs *before* `colnames(coordinates) <- c("x_cent","y_cent")`,
  ## which upstream also does first and which is why the assignment can never be the thing that
  ## fails.
  if (ncol(coordinates) != 2) {
    stop("Please check the input 'coordinates' and make sure it is a two column matrix.",
         call. = FALSE)
  }
  colnames(coordinates) <- c("x_cent", "y_cent")

  ## `collapse::fdist` returns a full `n x n` matrix, *not* a `dist`. The core does the same
  ## arithmetic; the dimnames come from the input's row names exactly as `fdist` sets them, and
  ## `fdist` copies them to the columns too.
  n <- nrow(coordinates)
  ## `fdist` needs at least two rows, and its own message is what a one-cell input produces. The
  ## error is raised here rather than in the kernel because it is a property of `collapse`, not of
  ## the arithmetic -- and because the kernel has no way to know which of upstream's shapes the
  ## caller meant. It fires *after* the two-column check, which is why a one-column matrix of one
  ## row reports the column error instead.
  if (n < 2L) {
    stop("If v is left empty, x needs to be a matrix with at least 2 rows", call. = FALSE)
  }
  ## `byrow = TRUE`: the kernel returns the `n x n` values row-major and `matrix()` fills
  ## column-major. The matrix is symmetric so the transpose is invisible here, but relying on that
  ## would hide a real ordering bug the moment anything asymmetric shares the path.
  ## `t()` is not a no-op detour: `as.numeric(matrix)` is **column-major** and both kernels index
  ## the flat vector row-major. Passing `as.numeric(coordinates)` unchanged hands the kernel the
  ## transpose. `fdist`'s output is symmetric, so that mistake is invisible here and visible only
  ## in `computeRegionDistance` -- which is the worse place to find it.
  full <- matrix(spatial_fdist(as.numeric(t(coordinates)), n, ncol(coordinates)),
                 nrow = n, ncol = n, byrow = TRUE)

  ## `collapse::fdist` returns a **`dist` object**, not a matrix, and `computeCellDistance` returns
  ## it unchanged. The rewrap has to reproduce three things `stats::as.dist` gets wrong, so it is
  ## written out rather than delegated:
  ##
  ##   * `as.dist` records a `call` attribute; `fdist` has none. `identical()` sees attributes.
  ##   * `as.dist` leaves `method` unset; `fdist` sets it to `"euclidean"`.
  ##   * `as.dist` synthesises `Labels = as.character(seq_len(n))` when the input has no row names;
  ##     `fdist` leaves `Labels` `NULL`. `as.character(1:3)` and `NULL` are different objects.
  ##
  ## Everything else -- the lower triangle, `Diag`, `Upper`, `Size` -- matches.
  rewrap <- function(m, rn) {
    structure(m[lower.tri(m, diag = FALSE)],
              Size = nrow(m), Labels = rn, Diag = FALSE, Upper = FALSE,
              method = "euclidean", class = "dist")
  }
  d.spatial <- rewrap(full, rownames(coordinates))

  ## `ratio` is a scalar multiplier here, not the per-sample vector `computeRegionDistance` takes.
  ## `if (!is.null(ratio))`, so a zero-length `ratio` *is* multiplied in -- unlike
  ## `computeRegionDistance`, where `ratio = NULL` empties the vector.
  if (!is.null(ratio)) {
    d.spatial <- d.spatial * ratio
  }

  ## The threshold only fires when **both** `interaction.range` and `tol` are given, and `tol` is
  ## an addend here rather than a tolerance: `> interaction.range + tol`. The message is upstream's,
  ## with its leading newline, and it is emitted before the assignment.
  if (!is.null(interaction.range) && !is.null(tol)) {
    message("\n Apply a predefined spatial distance threshold based on the interaction length...")
    d.spatial[d.spatial > (interaction.range + tol)] <- NaN
  }
  d.spatial
}

## ---------------------------------------------------------------------------------------
## `computeRegionDistance`, backed by Rust with the exact k-d tree in place of `AnnoyParam`.
## ---------------------------------------------------------------------------------------

#' The pinned upstream `computeRegionDistance`, verbatim.
#' @export
cellchatrs_upstream_computeRegionDistance <- function(...) {
  get("computeRegionDistance", envir = cellchatrs_upstream_cached())(...)
}

#' Spatial distances with the original Annoy nearest-neighbour algorithm.
#' @export
computeRegionDistance <- function(coordinates, meta,
                                  interaction.range = NULL, ratio = NULL, tol = NULL,
                                  k.min = 10, contact.dependent = TRUE,
                                  contact.range = NULL, contact.knn.k = NULL,
                                  do.symmetric = TRUE) {
  cellchatrs_upstream_computeRegionDistance(coordinates, meta, interaction.range,
    ratio, tol, k.min, contact.dependent, contact.range, contact.knn.k, do.symmetric)
}

#' Exact spatial distances, an explicit alternative to CellChat compatibility.
#' @export
cellchatrs_computeRegionDistance_exact <- function(coordinates, meta,
                                  interaction.range = NULL, ratio = NULL, tol = NULL,
                                  k.min = 10, contact.dependent = TRUE,
                                  contact.range = NULL, contact.knn.k = NULL,
                                  do.symmetric = TRUE) {
  group <- meta$group
  samples <- meta$samples
  ## `levels(factor)`, not `unique(as.character(x))`. Upstream sizes its arrays with
  ## `nlevels(group)` and then filters `level.use` against `unique(group)`, so a *declared* level
  ## with no cells survives into the output dimensions. Rebuilding the levels from the observed
  ## values would silently shrink the result.
  group_levels <- levels(group)
  samples_levels <- levels(samples)

  ## 0-based indices into those level vectors. `match(..., levels(x))` rather than `as.integer(f)`:
  ## both agree for a factor, but `match` also gives a usable answer for the `NA` group upstream
  ## would drop, instead of a silent 1-based position.
  gi <- match(as.character(group), group_levels) - 1L
  si <- match(as.character(samples), samples_levels) - 1L

  ## `as.numeric()` on a matrix is column-major; the core indexes row-major. Same trap as
  ## `computeCellDistance`, and unlike it here the transpose is *not* invisible -- it moves every
  ## cell to a different place and the k-d tree answers confidently about the wrong layout.
  res <- spatial_region_distance(
    as.numeric(t(coordinates)), nrow(coordinates), ncol(coordinates),
    group_levels, gi, samples_levels, si,
    ## `extendr` derives the R wrapper's parameter names from the Rust ones, so they are
    ## snake_case. Passing `interaction.range` here is an "unused argument" error, not a silent
    ## no-op -- and it is raised from *inside* the shim, so the traceback points at the kernel call
    ## rather than at the mismatch.
    interaction_range = interaction.range,
    ratio = if (is.null(ratio)) numeric(0) else as.numeric(ratio),
    tol = if (is.null(tol)) numeric(0) else as.numeric(tol),
    k_min = as.integer(k.min),
    contact_dependent = isTRUE(contact.dependent),
    contact_range = contact.range,
    contact_knn_k = if (is.null(contact.knn.k)) NULL else as.integer(contact.knn.k),
    do_symmetric = isTRUE(do.symmetric))
  if (!is.null(res[["__error"]])) stop(res[["__error"]], call. = FALSE)

  k <- length(group_levels)
  ## `rownames(d.spatial) <- levels(group); colnames(d.spatial) <- levels(group)`, unconditionally
  ## -- so with a declared-but-absent level the names do **not** line up with the data, because
  ## upstream's loop writes A's row into slot 1 while labelling it with the first level. That
  ## mislabelling is reproduced rather than corrected; see `region_distance` in the Rust core.
  d.spatial <- matrix(res[["d_spatial"]], nrow = k, ncol = k,
                      dimnames = list(group_levels, group_levels))
  adj.contact <- matrix(res[["adj_contact"]], nrow = k, ncol = k)
  list(d.spatial = d.spatial, adj.contact = adj.contact)
}

#' The pinned upstream `rankNetPairwise`, verbatim.
#' @export
cellchatrs_upstream_rankNetPairwise <- function(...) {
  get("rankNetPairwise", envir = cellchatrs_upstream_cached())(...)
}

## The shim's own pass-through. A pass-through is a legitimate outcome -- an `LR.use` shape the
## port does not accept, or a kernel that declined -- but it must never be *silent*. The first
## version of this shim fell back on a shape mismatch and every differential test still passed,
## because upstream's answer is by definition identical to upstream's: a gate that cannot tell a
## port from a pass-through reports "identical" for a function that is not ported at all. That is
## exactly what happened -- the kernel was reading R's `dim()` as a list of slice lengths,
## returning `__error`, and the shim quietly answering from upstream.
.cellchatrs_ranknet_pairwise_fallback <- function(object, LR.use, why) {
  warning(sprintf("cellchatrs: rankNetPairwise answered from upstream (%s)", why), call. = FALSE)
  cellchatrs_upstream_rankNetPairwise(object, LR.use)
}

#' `rankNetPairwise`, backed by the Rust kernel's ordering.
#'
#' Upstream's body is one `order` and a great deal of R plumbing:
#'
#' ```r
#' for (i in 1:numCluster) for (j in 1:numCluster) {
#'   data <- data.frame(pathway_index = index, interaction_name = ..., ..., row.names = rownames(pairLR.use))
#'   temp[[j]] <- data[with(data, order(pval, -prob)), ]
#' }
#' object@net$pairwiseRank <- pairwiseLR
#' ```
#'
#' The **ordering** is the only arithmetic, and it is done in Rust: `order(pval, -prob)` with two
#' numeric keys is a stable radix sort, so the permutation is a total function of
#' `(pval, -prob, index)` and `order_f64_multi` is a faithful transcription of it. `na.last` is R's
#' default `TRUE`, which is why a missing `pval` sorts last rather than first.
#'
#' Everything else stays in R, for the same reason `format` and `rankNet`'s row assembly do: the
#' `data.frame()` column order, `row.names = rownames(pairLR.use)`, the nested `list()` and the
#' `names(temp) <- colnames(prob)` are language behaviour rather than numerics, and reproducing
#' them in Rust would buy nothing but new ways to be wrong.
#'
#' `pairLR.use` supplies only the data-frame's **columns**; the ordering is over the whole third
#' dimension of `prob`/`pval`. Upstream does not reconcile the two, and that is load-bearing: with
#' `LR.use` shorter than `dim(prob)[3]`, its `data.frame(..., row.names = rownames(pairLR.use))`
#' raises `row names supplied are of the wrong length`. Slicing `prob` down to `LR.use` first --
#' which reads like a faithful interpretation of "use this subset of interactions" -- silently
#' turns that error into a successful shorter result.
#'
#' Upstream's loops are `for (i in 1:numCluster) for (j in 1:numCluster)` with
#' `numCluster <- dim(prob)[1]`, so it reads the leading `k x k` block and never touches the rest.
#' `rownames(prob)` and `colnames(prob)` are still used for the names, which is where a non-square
#' `prob` raises.
#'
#' @export
rankNetPairwise <- function(object, LR.use = NULL) {
  if (is.null(LR.use)) {
    pairLR.use <- object@LR$LRsig
  } else {
    pairLR.use <- LR.use
  }
  net <- object@net
  prob <- net$prob
  pval <- net$pval
  numCluster <- dim(prob)[1]
  ## The ordering is over the **whole** third dimension, and `pairLR.use` supplies only the
  ## data-frame's columns. Upstream does not reconcile the two, and that is load-bearing: with
  ## `LR.use` shorter than `dim(prob)[3]`, its
  ## `data.frame(..., row.names = rownames(pairLR.use))` raises
  ## `row names supplied are of the wrong length`, because `probij` has `dim(prob)[3]` values and
  ## the row names do not. Slicing `prob` to `LR.use` first -- which looked like a faithful
  ## reading of "use this subset of interactions" -- silently turned that hard error into a
  ## successful shorter result.
  n <- dim(prob)[3]
  ## Upstream's loops are `for (i in 1:numCluster) for (j in 1:numCluster)` with
  ## `numCluster <- dim(prob)[1]`, so it reads the **leading `k x k` block** and never touches the
  ## rest. Taking that block here keeps the kernel's `k * k` orderings matching the number the shim
  ## consumes when the second dimension is larger than the first. The block is taken on `prob`/
  ## `pval` only: `rownames(prob)` and `colnames(prob)` below stay the originals, because
  ## `names(temp) <- colnames(prob)` on a wider array is precisely where upstream raises
  ## `'names' attribute [5] must be the same length as the vector [2]`.
  if (dim(prob)[2] != numCluster || dim(pval)[2] != numCluster) {
    blk <- seq_len(numCluster)
    prob.n <- prob[blk, blk, , drop = FALSE]
    pval.n <- pval[blk, blk, , drop = FALSE]
  } else {
    prob.n <- prob
    pval.n <- pval
  }
  if (length(dim(prob.n)) != 3L || !identical(dim(prob.n), dim(pval.n))) {
    return(.cellchatrs_ranknet_pairwise_fallback(object, LR.use,
      "prob and pval have different dimensions"))
  }

  ## Flattened here rather than in Rust, following every other binding in this package:
  ## `extendr` has no `TryFrom<Robj>` for a 3-d array, and `as.vector` is R's own column-major
  ## flattening, so the two sides agree on the element order by construction rather than by
  ## convention.
  dp <- as.integer(dim(prob.n))
  ord <- ranknet_pairwise_orders(as.vector(prob.n), dp, as.vector(pval.n), dp)
  ## `[["__error"]]`, not `$__error`: R's parser rejects a `$` immediately followed by an
  ## underscore, so the latter is a *parse* error rather than a runtime one.
  if (is.list(ord) && !is.null(ord[["__error"]])) {
    return(.cellchatrs_ranknet_pairwise_fallback(object, LR.use,
      paste0("the kernel declined: ", ord[["__error"]])))
  }
  ord <- unlist(ord, use.names = FALSE)
  if (n == 0L || length(ord) != numCluster * numCluster * n) {
    return(.cellchatrs_ranknet_pairwise_fallback(object, LR.use,
      sprintf("the kernel returned %d orderings, expected %d",
              length(ord), numCluster * numCluster * n)))
  }

  row.names <- rownames(pairLR.use)
  pairwiseLR <- list()
  for (i in 1:numCluster) {
    temp <- list()
    for (j in 1:numCluster) {
      ## The kernel returns one permutation per `(i, j)` in **column-major** order -- `j` outer,
      ## `i` inner -- because that is the order the slices appear in the flattened array. A shared
      ## counter advancing in the `i`-outer loop order below looks right and is not: it transposes
      ## the `(i, j)` indexing, which shows up only when a group pair's ordering is not the identity,
      ## and never as an error.
      at <- ((j - 1L) * numCluster + (i - 1L)) * n
      pvalij <- as.vector(pval.n[i, j, ])
      probij <- as.vector(prob.n[i, j, ])
      index <- 1:n
      data <- data.frame(pathway_index = index,
                         interaction_name = pairLR.use$interaction_name,
                         interaction_name_2 = pairLR.use$interaction_name_2,
                         pathway_name = pairLR.use$pathway_name,
                         ligand = pairLR.use$ligand,
                         receptor = pairLR.use$receptor,
                         prob = probij, pval = pvalij, row.names = row.names)
      temp[[j]] <- data[ord[(at + 1L):(at + n)], , drop = FALSE]
    }
    names(temp) <- colnames(prob)
    pairwiseLR[[i]] <- temp
  }
  names(pairwiseLR) <- rownames(prob)
  object@net$pairwiseRank <- pairwiseLR
  return(object)
}

## ---------------------------------------------------------------------------------------
## `netAnalysis_computeCentrality`, split between Rust and igraph.
## ---------------------------------------------------------------------------------------

#' The pinned upstream `netAnalysis_computeCentrality`, verbatim.
#' @export
cellchatrs_upstream_netAnalysis_computeCentrality <- function(...) {
  get("netAnalysis_computeCentrality", envir = cellchatrs_upstream_cached())(...)
}

#' `netAnalysis_computeCentrality`, with the deterministic measures in Rust.
#'
#' Upstream's `computeCentralityLocal` computes eleven measures per pathway network. Four are
#' iterative solvers whose output is not stable across runs -- `hub_score` (deprecated in igraph
#' 2.0.3), `authority_score`, `eigen_centrality` (ARPACK) and `page_rank` (PRPACK) all disagree
#' with themselves on identical input, measured with no seed set. Bit-parity with a
#' nondeterministic oracle is meaningless, so those four stay in R and call igraph directly --
#' the same package, the same function, the same input, which is identical by construction rather
#' than by reimplementation.
#'
#' The other seven are pure functions of the input and live in `r-core::centrality`, verified
#' bit-for-bit against installed igraph 2.3.4 on a 77-case corpus (`centrality_parity`):
#' unweighted in/out degrees, weighted strengths (plain sequential `f64` accumulation in edge-ID
#' order -- *not* R's long-double `rowSums`, which agrees only 3802/4669), and betweenness
#' (Dijkstra with dist-plus-one encoding, exact 2-way-heap tie rules, epsilon comparisons at
#' 1e-10, Brandes accumulation). `flowbet`/`infocent` are `sna` calls in upstream's own
#' `tryCatch`, copied verbatim including the zeros fallback: `sna` is installed here, so both
#' sides run the real functions and agree on nonzero values; a host without `sna` gets zeros
#' from both sides through the same `tryCatch`.
#'
#' The preamble (`pval[prob == 0] <- 1`, `prob[pval >= thresh] <- 0`, `net.name`, the
#' `pbapply`/`future.apply` dispatch, the `names(centr.all)` assignment and the slot writeback)
#' is upstream verbatim, because none of it is arithmetic the kernel owns. `igraph`'s
#' `hub_score` deprecation warning fires here exactly as often as upstream -- same function,
#' same lifecycle state -- so warning behaviour matches by construction rather than by silencing.
#' @export
netAnalysis_computeCentrality <- function(object = NULL, slot.name = "netP", net = NULL,
                                          net.name = NULL, thresh = 0.05) {
  if (is.null(net)) {
    prob <- methods::slot(object, slot.name)$prob
    pval <- methods::slot(object, slot.name)$pval
    pval[prob == 0] <- 1
    prob[pval >= thresh] <- 0
    net = prob
  }
  if (is.null(net.name)) {
    net.name <- dimnames(net)[[3]]
  }
  if (length(dim(net)) == 3) {
    nrun <- dim(net)[3]
    my.sapply <- ifelse(
      test = future::nbrOfWorkers() == 1,
      yes = pbapply::pbsapply,
      no = future.apply::future_sapply
    )
    centr.all = my.sapply(
      X = 1:nrun,
      FUN = function(x) {
        net0 <- net[ , , x]
        return(.cellchatrs_centrality_one(net0))
      },
      simplify = FALSE
    )
  } else {
    centr.all <- as.list(.cellchatrs_centrality_one(net))
  }
  names(centr.all) <- net.name
  if (is.null(object)) {
    return(centr.all)
  } else {
    slot(object, slot.name)$centr <- centr.all
    return(object)
  }
}

## One pathway network: the five deterministic measures from Rust, the four solver measures from
## igraph, the two `sna` measures from upstream's own `tryCatch`, assembled in upstream's order
## under upstream's names.
.cellchatrs_centrality_one <- function(net0) {
  ## `net[,,x]` drops dimensions, so a 1x1xN pathway arrives here as a length-1 *vector*, not a
  ## matrix. Upstream fails on it in `rowSums(net > 0)` -- the first line of its
  ## `computeCentralityLocal` -- with "'x' must be an array of at least two dimensions". Evaluate
  ## the same expression rather than inventing text: the message is then upstream's by construction,
  ## and any future change to the phrasing propagates automatically. (What the Rust side would say
  ## instead -- an extendr conversion error -- names neither the function nor the condition.)
  if (is.null(dim(net0))) {
    rowSums(net0 > 0)
  }
  k <- nrow(net0)
  ## Row-major flatten: `as.numeric()` on a matrix is column-major and the core indexes
  ## `m[row * k + col]`, so the transpose goes over the boundary, as elsewhere in this package.
  res <- centrality_deterministic(as.numeric(t(net0)), k)
  if (!is.null(res[["__error"]])) stop(res[["__error"]], call. = FALSE)
  ## igraph's warning text verbatim (measured, including the lack of a trailing call note being
  ## irrelevant: gates compare `conditionMessage()`). Raised here rather than in Rust because a
  ## warning's call belongs to R.
  if (isTRUE(res[["tiny_weights"]])) {
    warning(paste0("Some weights are smaller than epsilon, calculations may suffer from ",
                   "numerical precision issues.\nSource: centrality/betweenness.c:441"),
            call. = FALSE)
  }
  G <- igraph::graph_from_adjacency_matrix(net0, mode = "directed", weighted = TRUE)
  nm <- rownames(net0)
  named <- function(v) { names(v) <- nm; v }
  centr <- vector("list")
  centr$outdeg_unweighted <- named(res$outdeg_unweighted)
  centr$indeg_unweighted <- named(res$indeg_unweighted)
  centr$outdeg <- named(res$outdeg)
  centr$indeg <- named(res$indeg)
  centr$hub <- igraph::hub_score(G)$vector
  centr$authority <- igraph::authority_score(G)$vector
  centr$eigen <- igraph::eigen_centrality(G)$vector
  centr$page_rank <- igraph::page_rank(G)$vector
  igraph::E(G)$weight <- 1 / igraph::E(G)$weight
  centr$betweenness <- named(res$betweenness)
  centr$flowbet <- tryCatch({
    sna::flowbet(net0)
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net0), ncol = 1))
  })
  centr$info <- tryCatch({
    sna::infocent(net0, diag = TRUE, rescale = TRUE, cmode = "lower")
  }, error = function(e) {
    as.vector(matrix(0, nrow = nrow(net0), ncol = 1))
  })
  centr
}

## Compatible replacement for the original Eigen helper (not part of the public R API).
ComputeSNN <- function(nn_ranked, prune) {
  graph <- compute_snn(as.integer(nn_ranked), nrow(nn_ranked), ncol(nn_ranked), prune)
  methods::new("dgCMatrix", i = graph$i, p = graph$p, x = graph$x,
               Dim = as.integer(rep(nrow(nn_ranked), 2L)))
}
