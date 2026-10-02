#!/usr/bin/env python3
"""Emit R/delegate.R: verbatim-upstream delegating wrappers for drop-in completeness.

The shim ports the numeric kernel to Rust and overrides the functions it accelerates. Every other
*computational* name upstream exports must still resolve when a user writes
`library(cellchatrs)` instead of `library(CellChat)` -- otherwise the package is not a drop-in
replacement, it is a namespace that fails on the tutorial's third line. Each wrapper below calls
the pinned upstream body verbatim through `cellchatrs_upstream_cached()`.

What this is and is not, stated once so no wrapper needs to restate it:

* It is NOT the rankNetPairwise-fallback failure mode. That fallback *claimed* to be a port and
  silently answered from upstream on shape mismatch, so the gate compared upstream with itself.
  These wrappers make no port claim: their docs say verbatim-upstream, the parity ledger records
  them at the pass-through rung (never `exact`), and where a delegated function feeds a ported
  one the gate compares the ported one's output, not the delegation.
* Heavy-dependency functions (NMF/Seurat/IRLBA-backed manifold code, presto-gated branches) fail
  exactly as upstream fails when the dependency is absent -- same function body, same error --
  because delegation is lazy: nothing is evaluated until called.
* Visualization.R and app.R contents are deliberately absent (locked decision 14.1), as are the
  `netAnalysis_*` and `compareInteractions` functions that return ggplot objects. The full
  exclusion list with reasons is in the header this script emits.

Usage: python3 scripts/gen_delegate.py  # rewrites R/delegate.R and patches NAMESPACE
Then:  python3 scripts/gen_man.py      # .Rd files for the new exports
"""

import re
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# name -> (upstream file, note). Note is appended to the roxygen block when the reason is specific;
# the default prose covers the rest.
DELEGATE = {
    # Constructors and object surgery: no numerics, must resolve for any pipeline to start.
    "createCellChat": ("CellChat_class.R", "Upstream's constructor, verbatim. Takes no numeric "
                       "arguments, so there is nothing to accelerate; the reference constructor is "
                       "the constructor."),
    "setIdent": ("utilities.R", ""),
    "updateCellChat": ("CellChat_class.R", ""),
    "subsetCellChat": ("CellChat_class.R", ""),
    "updateClusterLabels": ("utilities.R", ""),
    "updateCCC_score": ("utilities.R", ""),
    "updateCellChatDB": ("database.R", ""),
    "showDatabaseCategory": ("database.R", ""),
    "checkGeneSymbol": ("database.R", ""),
    "addMeta": ("utilities.R", ""),
    "addReduction": ("utilities.R", ""),
    "mergeCellChat": ("CellChat_class.R", "Bare name delegating to the documented "
                      "cellchatrs_upstream_mergeCellChat wrapper, so there is exactly one code path. "
                      "Pass-through rung; see tests/parity/check_merge.R."),
    # Per-gene expression kernels and means: reimplemented in r-core and used by the kernel, but the
    # bare R-level names are upstream's R semantics verbatim (match.arg, NA handling, dimnames).
    "computeExpr_LR": ("modeling.R", ""),
    "computeExpr_complex": ("modeling.R", ""),
    "computeExpr_coreceptor": ("modeling.R", ""),
    "computeExpr_agonist": ("modeling.R", ""),
    "computeExpr_antagonist": ("modeling.R", ""),
    "computeExprGroup_agonist": ("modeling.R", ""),
    "computeExprGroup_antagonist": ("modeling.R", ""),
    "triMean": ("modeling.R", ""),
    "geometricMean": ("modeling.R", ""),
    # Utilities numerics: deterministic table surgery, no floating point worth porting.
    "extractGene": ("database.R", ""),
    "extractGeneSubset": ("database.R", ""),
    "extractGeneSubsetFromPair": ("database.R", ""),
    "extractLRfromGenes": ("database.R", ""),
    "extractEnrichedLR": ("analysis.R", ""),
    "findEnrichedSignaling": ("analysis.R", ""),
    "identifyEnrichedInteractions": ("modeling.R", "Downstream of the ported filterCommunication; the enrichment statistic itself is a rank test over the interaction table."),
    "identifyOverExpressedInteractions": ("utilities.R", "Pure set-membership filtering over the DB "
                          "tables; no floating point, no randomness. Covered end-to-end by the "
                          "tutorial reproduction, which routes the tutorial's own call through here."),
    "identifyOverExpressedLigandReceptor": ("utilities.R", ""),
    "mergeInteractions": ("analysis.R", ""),
    "searchPair": ("database.R", ""),
    "getMaxWeight": ("analysis.R", ""),
    "normalizeData": ("utilities.R", ""),
    "scaleData": ("utilities.R", ""),
    "scaleMat": ("utilities.R", ""),
    "smoothData": ("utilities.R", ""),
    "sketchData": ("utilities.R", "Needs Seurat at call time, exactly as upstream does."),
    "selectK": ("analysis.R", ""),
    "computeEnrichmentScore": ("analysis.R", ""),
    "computeNetD_structure": ("analysis.R", ""),
    "computeNetSimilarity": ("analysis.R", ""),
    "computeNetSimilarityPairwise": ("analysis.R", ""),
    "rankSimilarity": ("analysis.R", ""),
    "liftCellChat": ("CellChat_class.R", ""),
    # Manifold / NMF-backed analysis: lazy delegation fails exactly as upstream fails.
    "identifyCommunicationPatterns": ("analysis.R", "Needs NMF at call time, exactly as upstream does."),
    "netClustering": ("analysis.R", "Needs NMF/IRLBA at call time, exactly as upstream does."),
    "netEmbedding": ("analysis.R", "Needs Seurat/IRLBA at call time, exactly as upstream does."),
    "netMappingDEG": ("analysis.R", ""),
    "computeLaplacian": ("analysis.R", ""),
    "computeEigengap": ("analysis.R", ""),
    "buildSNN": ("analysis.R", "Needs Seurat at call time, exactly as upstream does."),
    "runPCA": ("utilities.R", "Needs Seurat at call time, exactly as upstream does."),
    "runUMAP": ("utilities.R", "Needs Seurat at call time, exactly as upstream does."),
    "preProcMultiomics": ("utilities.R", "Needs Seurat/Signac at call time, exactly as upstream does."),
}

# Deliberately absent, with the reason. Visualization.R/app.R contents per locked decision 14.1;
# the plot-returning analysis functions with them, even though "comparison" is in requirement 1's
# list: compareInteractions and netAnalysis_contribution both RETURN ggplot objects, so porting
# their numerics would still leave a plot function, and the tutorial reproduction covers the numeric
# slots feeding the plots instead.
EXCLUDED = {
    "compareInteractions": "returns a ggplot object; the numeric comparison workflow it fronts "
                           "(mergeCellChat + rankNet comparison mode + rankNetPairwise) is ported "
                           "and gated, and the tutorial covers the feeding slots.",
    "netAnalysis_contribution": "builds and returns a ggplot object.",
    "netAnalysis_signalingRole_network": "returns a ggplot object.",
    "netVisual* / plotGeneExpression / dotPlot / barPlot / pieChart / spatial*Plot / StackedVlnPlot":
        "visualization.R contents, locked decision 14.1.",
    "netAnalysis_signalingRole_scatter/heatmap, river, dot, diffInteraction/scatter plots":
        "visualization.R contents, locked decision 14.1.",
    "scPalette/ggPalette/colorRamp3/CellChat_theme_opts": "plot theming helpers, 14.1.",
    "runCellChatApp": "app.R Shiny application, 14.1.",
}

HEADER = '''# Verbatim-upstream delegating wrappers: drop-in completeness for the computational surface.
#
# The shim ports the numeric kernel to Rust and overrides the functions it accelerates. Every other
# *computational* name upstream exports must still resolve, or `library(cellchatrs)` fails on the
# tutorial's third line. Each wrapper below calls the pinned upstream body verbatim (commit
# 75253cd0c9e68410e6e721a6d3a0419a1d7e358f) through `cellchatrs_upstream_cached()`, which sources the
# pinned files once per session.
#
# This is not the rankNetPairwise-fallback failure mode: nothing here claims to be a port. The docs
# say verbatim-upstream, the parity ledger records delegated names at the pass-through rung (never
# `exact`), and where a delegated function feeds a ported one the gates compare the ported output.
# Heavy-dependency callees (NMF/Seurat/IRLBA/presto-gated branches) fail exactly as upstream fails
# when the dependency is absent -- same body, same error -- because delegation evaluates nothing
# until called.
#
# Deliberately absent (locked decision PLAN.md 14.1 keeps visualization.R and app.R in R):
#
%s
#
# Generated by scripts/gen_delegate.py -- do not hand-edit; edit the table there and re-run.
# (Hand edits would be overwritten, which is the point: the export list must stay in sync with
# NAMESPACE and man/ by construction rather than by discipline.)
''' % "\n".join(f"# - {k}: {v}" for k, v in EXCLUDED.items())

BODY = '''
#' The pinned upstream `%s`, verbatim.
#'
#' Drop-in delegation: resolves the name for `library(cellchatrs)` users by calling the pinned
#' upstream body (see the header of this file for what delegation does and does not claim).
#' Upstream file: `R/%s`.%s
#' @export
%s <- function(...) {
  get("%s", envir = cellchatrs_upstream_cached())(...)
}
'''

MERGE_SPECIAL = '''
#' `mergeCellChat`, resolving to the documented upstream wrapper.
#'
#' Bare-name alias for `cellchatrs_upstream_mergeCellChat`, so there is exactly one code path and
#' one documented rung (pass-through; see tests/parity/check_merge.R). Upstream file:
#' `R/CellChat_class.R`.
#' @export
mergeCellChat <- function(...) {
  cellchatrs_upstream_mergeCellChat(...)
}
'''

CREATE_SPECIAL = '''
#' `createCellChat`, resolving to the documented upstream wrapper.
#'
#' Bare-name alias for `cellchatrs_upstream_createCellChat`: the constructor takes no numeric
#' arguments, so there is nothing to accelerate and the reference constructor is the constructor.
#' Upstream file: `R/CellChat_class.R`.
#' @export
createCellChat <- function(...) {
  cellchatrs_upstream_createCellChat(...)
}
'''


def main() -> None:
    # Verify every delegated name actually exists in the pinned source before emitting anything, so a
    # typo or an upstream rename fails here with the name attached rather than as a lazy
    # "object not found" at first call.
    missing = []
    for name, (f, _) in DELEGATE.items():
        if name in ("mergeCellChat", "createCellChat"):
            continue
        src = (ROOT / ".." / "tmp" / "opencode" / "CellChat" / "R" / f).resolve()
        if not src.is_file():
            missing.append(f"{name}: {f} not found")
            continue
        text = src.read_text()
        if not re.search(rf"^{re.escape(name)} <- function", text, re.M):
            missing.append(f"{name}: no definition in R/{f}")
    if missing:
        raise SystemExit("delegation target check failed:\n  " + "\n  ".join(missing))

    parts = [HEADER]
    for name, (f, note) in DELEGATE.items():
        if name == "mergeCellChat":
            parts.append(MERGE_SPECIAL)
            continue
        if name == "createCellChat":
            parts.append(CREATE_SPECIAL)
            continue
        extra = f" {note}" if note else ""
        parts.append(BODY % (name, f, extra, name, name))
    (ROOT / "R" / "delegate.R").write_text("\n".join(parts) + "\n")

    ns = (ROOT / "NAMESPACE").read_text().rstrip() + "\n"
    added = []
    for name in DELEGATE:
        if not re.search(rf"^export\({re.escape(name)}\)$", ns, re.M):
            ns += f"export({name})\n"
            added.append(name)
    (ROOT / "NAMESPACE").write_text(ns)
    print(f"R/delegate.R written ({len(DELEGATE)} wrappers); NAMESPACE +{len(added)} exports")


if __name__ == "__main__":
    main()
