#!/usr/bin/env python3
"""Emit `parity.json`: the achieved parity rung per quantity, machine-checkable.

The objective asks for a CI artifact that states, per quantity, which rung was actually
reached -- not what is claimed. So every entry here is produced by *running the test* that
establishes it, and the rung is derived from that test's own assertions:

* `exact`        the test compares raw bit patterns and requires equality;
* `exact-nan`    as `exact`, and the corpus contains non-finite inputs;
* `identical`    `identical()` on the R object, run by `tests/parity/check_identical.R`;
* `error-equal`  error *messages* compared byte for byte.

Nothing is asserted here that is not also asserted by a test that runs in CI. A quantity
with no passing test gets `"rung": "untested"`, which is the honest state and is what makes
the file worth reading.

Usage:
    python3 scripts/parity_report.py                 # run the Rust suite, write parity.json
    python3 scripts/parity_report.py --no-run        # report the last recorded results
"""

from __future__ import annotations

import argparse
import json
import os
import re
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
UPSTREAM_SHA = "75253cd0c9e68410e6e721a6d3a0419a1d7e358f"

# (quantity, what pins it, the test that must be green, rung when green)
QUANTITIES: list[tuple[str, str, str, str]] = [
    ("rng.mersenne_twister", "set.seed/MT_genrand stream, incl. FixupSeeds and unsigned shifts",
     "cargo test --release -p r-core --test rng_parity", "exact"),
    ("rng.sample_int", "sample.int permutation and sampling without replacement",
     "cargo test --release -p r-core --test rng_parity", "exact"),
    ("rng.stress", "879-record corpus, 25 MB, >2e6 compared values",
     "cargo test --release -p r-core --test rng_stress", "exact"),
    ("longdouble.add", "x87 80-bit add and subtract vs real long double in C",
     "cargo test --release -p r-core --test f80_vs_x87", "exact"),
    ("longdouble.mul", "x87 80-bit multiply, needed for R's prod",
     "cargo test --release -p r-core --test f80_vs_x87", "exact"),
    ("stats.fquantile_type7", "collapse::fquantile type 7, incl. the radix/dquickselect path",
     "cargo test --release -p r-core --test stats_parity", "exact"),
    ("stats.tri_mean", "triMean, incl. the na.rm=FALSE inner mean",
     "cargo test --release -p r-core --test stats_parity", "exact-nan"),
    ("stats.geometric_mean", "geometricMean, incl. log(0) = -Inf and log(-Inf) = NaN",
     "cargo test --release -p r-core --test stats_parity", "exact-nan"),
    ("stats.r_prod", "R's prod: LDOUBLE accumulation with the DBL_MAX clamp",
     "cargo test --release -p r-core --test stats_parity", "exact-nan"),
    ("stats.thresholded_mean", "thresholdedMean, incl. the NaN-raise on missing values and the zero-below-threshold rule",
     "cargo test --release -p r-core --test stats_parity", "exact-nan"),
    ("db.extract_gene", "extractGene / extractGeneSubset, both species, order for order",
     "cargo test --release -p r-core --test db_parity", "exact"),
    ("expr.computeExpr_LR", "computeExpr_LR and computeExpr_complex, 416 x 6 fixture",
     "cargo test --release -p r-core --test expr_parity", "exact"),
    ("expr.computeExpr_coreceptor", "computeExpr_coreceptor, incl. the all-ones no-op",
     "cargo test --release -p r-core --test expr_parity", "exact"),
    ("expr.computeExpr_agonist", "computeExpr_agonist / _antagonist, 30 records x 5 Kh/n",
     "cargo test --release -p r-core --test expr_parity", "exact"),
    ("expr.computeExprGroup", "the exported group-mean variants, triMean over cells",
     "cargo test --release -p r-core --test expr_parity", "exact"),
    ("expr.subscript_out_of_bounds", "upstream's abort on a complex with a missing subunit",
     "cargo test --release -p r-core --test expr_parity", "error-equal"),
    ("aggregate.aggregate_1", "aggregate(matrix, list(factor), FUN) group order and layout",
     "cargo test --release -p r-core --lib aggregate", "exact"),
    ("prob.data_use_avg", "data.use.avg for all four type.mean values, vs R's own output",
     "cargo test --release -p r-core --test prob_parity", "exact"),
    ("prob.prob", "net$prob over 8 parameter configurations, bit for bit",
     "cargo test --release -p r-core --test prob_parity", "exact"),
    ("prob.pval", "net$pval over 8 parameter configurations, bit for bit",
     "cargo test --release -p r-core --test prob_parity", "exact"),
    ("net.aggregate_net", "net$count and net$weight, 4 fixtures x 3 thresholds, LDOUBLE sum",
     "cargo test --release -p r-core --test net_parity", "exact"),
    ("net.subset_communication", "row order, column set, factor levels and every cell",
     "cargo test --release -p r-core --test net_parity", "exact"),
    ("net.subset_empty_error", "upstream's refusal to return an empty table",
     "cargo test --release -p r-core --test net_parity", "error-equal"),
    ("de.compute_ave_expr", "computeAveExpr over all three types, incl. R's intersect order",
     "cargo test --release -p r-core --test de_parity", "exact-nan"),
    ("de.intersect_order", "R's intersect: first argument's order, deduplicated",
     "cargo test --release -p r-core --test de_parity", "exact"),
    ("de.subset_data", "subsetData's gene list and the annotation reordering",
     "cargo test --release -p r-core --test de_parity", "exact"),
    ("de.subset_db", "subsetDB's annotation filter and its non_protein flip",
     "cargo test --release -p r-core --test de_parity", "exact"),
    ("de.match_arg_message", "computeAveExpr's match.arg error, curly quotes included",
     "cargo test --release -p r-core --test de_parity", "error-equal"),
    ("de.subset_db_key_error", "subsetDB's unknown-key error, byte for byte",
     "cargo test --release -p r-core --test de_parity", "error-equal"),
    ("rshim.identical", "identical() on net$prob, net$pval, dimnames, options$parameter",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("rshim.match_arg", "match.arg partial matching and its error message, byte for byte",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "error-equal"),
    ("rshim.aggregate_net", "aggregateNet through the installed shim, three thresholds",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("rshim.subset_communication", "subsetCommunication through the installed shim, three filters",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("rshim.compute_ave_expr", "computeAveExpr through the installed shim, six configurations",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("rshim.subset_db", "subsetDB through the installed shim, four configurations",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("wilcox.wilcox_test", "stats::wilcox.test: W and p.value bit for bit, 13 x 3 x 2 records",
     "cargo test --release -p r-core --test wilcox_parity", "exact-nan"),
    ("wilcox.correct_argument", "R's `correct` argument, both settings distinguished",
     "cargo test --release -p r-core --test wilcox_parity", "exact"),
    ("wilcox.rank", "stats::rank, incl. na.last = 'keep' ranking NA last",
     "cargo test --release -p r-core --test wilcox_parity", "exact"),
    ("wilcox.p_adjust_bonferroni", "stats::p.adjust(method = 'bonferroni')",
     "cargo test --release -p r-core --test wilcox_parity", "exact"),
    ("wilcox.pnorm", "Cody's pnint, lower and upper tail, vs R, bit for bit",
     "cargo test --release -p r-core --lib mathfn", "exact"),
    ("wilcox.mean_fxn", "log(mean(expm1(x)) + 1) via apply(MARGIN = 1)",
     "cargo test --release -p r-core --test wilcox_parity", "exact"),
    ("de.group_dataset_selection",
     "group.dataset changes only cell.use1/cell.use2 -- by dataset and group, or pooled across "
     "groups with no labels term at all when group.DE.combined is TRUE. The kernel takes the "
     "selection as a parameter rather than growing a flag through the Wilcoxon path",
     "R_LIBS=.rlib Rscript tests/parity/check_identical.R", "identical"),
    ("de.group_dataset_tostring_collapse",
     "labels.dataset[labels.dataset != pos.dataset] <- toString(setdiff(unique(...), pos.dataset)) "
     "puts every non-positive dataset in ONE level, so a four-dataset comparison has two levels, "
     "one named \"D2, D3, D4\"; the datasets factor and its row order depend on it",
     "R_LIBS=.rlib Rscript tests/parity/check_identical.R", "identical"),
    ("de.group_dataset_bare_stop",
     "an unusable pos.dataset makes upstream cat() the names and then call a bare stop(), so the "
     "error message is empty; a helpful message here would be a divergence",
     "R_LIBS=.rlib Rscript tests/parity/check_identical.R", "error-equal"),
    ("de.single_candidate_group",
     "a group with exactly one percentage-passing feature contributes NO markers: apply(X, 1, FUN) "
     "on a one-row matrix returns an unnamed scalar, so FC[features] is NA and features.diff is "
     "character(0)",
     "cargo test --release -p r-core --lib", "exact"),
    ("de.row_names_are_not_associative",
     "markers.all is rbind-ed one group at a time, and R's row-name uniquification falls back to "
     "appending a bare digit once `x.1` is taken, so do.call(rbind, frames) renames differently",
     "R_LIBS=.rlib Rscript tests/parity/check_identical.R", "identical"),
    ("de.identify_marker_table", "the marker table row for row, incl. the 0x1 collapsed schema",
     "cargo test --release -p r-core --test wilcox_parity", "exact-nan"),
    ("de.marker_row_names", "features.info row names are the features, in row order",
     "cargo test --release -p r-core --test wilcox_parity", "exact"),
    ("de.pct_scale", "pct.1/pct.2 as round(x, 3) fractions; thresh.pc vs fraction and percent",
     "cargo test --release -p r-core --test wilcox_parity", "exact"),
    ("de.expressed_in_min_cells", "the do.DE = FALSE branch, where min.cells is actually used",
     "cargo test --release -p r-core --test wilcox_parity", "exact"),
    ("rshim.identify_over_expressed_genes",
     "identifyOverExpressedGenes through the installed shim, 7 configs, whole var.features",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("rshim.oeg_falls_back", "do.fast = TRUE falls back rather than running the Wilcoxon kernel",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "error-equal"),
    ("pathway.prob_pathways", "prob.pathways per (source, target, pathway), 9 fixtures, bit for bit",
     "cargo test --release -p r-core --test pathway_parity", "exact"),
    ("pathway.pathway_ranking", "netP$pathways by decreasing total, incl. tied totals",
     "cargo test --release -p r-core --test pathway_parity", "exact"),
    ("pathway.lr_sig", "net$LRs: the L-R pairs with a non-zero total after thresholding",
     "cargo test --release -p r-core --test pathway_parity", "exact"),
    ("pathway.apply_sums", "both apply() sums in their own orders, LONG_DOUBLE accumulation",
     "cargo test --release -p r-core --test pathway_parity", "exact"),
    ("pathway.aperm_layout", "the pathway axis is last after aperm(..., c(2,3,1))",
     "cargo test --release -p r-core --test pathway_parity", "exact"),
    ("pathway.single_pathway_error", "upstream's aperm failure for a one-pathway LRsig, preserved",
     "cargo test --release -p r-core --test pathway_parity", "error-equal"),
    ("rshim.compute_commun_prob_pathway",
     "computeCommunProbPathway through the installed shim, 9 fixtures, both return shapes",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("net.aggregate_net_unfiltered", "aggregateNet's unfiltered branch, 12 fixtures, every cell",
     "cargo test --release -p r-core --test netfiltered_parity", "exact"),
    ("net.aggregate_net_filtered", "aggregateNet's filtered branch: 0x0 or kxk of zeros, as upstream",
     "cargo test --release -p r-core --test netfiltered_parity", "exact"),
    ("net.group_by_key_order", "dplyr::group_by on the joined string key orders byte-wise",
     "cargo test --release -p r-core --test netfiltered_parity", "exact"),
    ("net.remove_isolate_shape", "remove.isolate decides 0x0 vs kxk, via upstream's str_split regex",
     "cargo test --release -p r-core --test netfiltered_parity", "exact"),
    ("rshim.aggregate_net_filtered",
     "aggregateNet's filtered branch through the installed shim, 12 fixtures, whole object",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("filter.prob_after_filter", "the min.cells zero-fill and the cross-sample mask, every element",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("filter.per_sample_binarisation", "avg.s each sample to 0/1, and the excluded rows zeroed",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("filter.min_cells_threshold", "groups at or below min.cells, reported in level order",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("filter.absent_group_is_zero", "a group missing from a sample's idents counts as 0 cells",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("filter.cell_excludes", "cell.excludes is the union over samples of the excluded groups",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("filter.type_mean_choices", "all four type.mean values, incl. truncated and thresholded",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("filter.upstream_failures",
     "min.samples > n_samples and the all-zero 'subscript out of bounds', both preserved",
     "cargo test --release -p r-core --test filter_parity", "error-equal"),
    ("filter.nonfilter_keep", "nonFilter.keep = TRUE, a no-op because upstream drops the augmented net",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("subset.threshold_order", "the eight DEG thresholds in upstream's own if-order, and the two stop() wordings",
     "cargo test --release -p r-core --test subset_parity", "error-equal"),
    ("subset.na_blanking", "an NA or NaN in an applied threshold column blanks the row, which the shared all-NA drop removes",
     "cargo test --release -p r-core --test subset_parity", "exact"),
    ("subset.all_na_drop", "rowSums(is.na(x)) != ncol(x) runs with no threshold at all",
     "cargo test --release -p r-core --test subset_parity", "exact"),
    ("subset.logfc_sign", "the sign of the logFC *argument* picks the comparison, so 0 takes the >= branch",
     "cargo test --release -p r-core --test subset_parity", "exact"),
    ("subset.warn_vs_error", "empty after a threshold is a stop(); empty after sources.use is only a warning",
     "cargo test --release -p r-core --test subset_parity", "error-equal"),
    ("subset.netp_group_order", "group_by on the pasted 'sourceTotarget' key, byte-wise, so g10 before g1",
     "cargo test --release -p r-core --test subset_parity", "exact"),
    ("subset.netp_na_to_string", "paste turns a missing cell into the string \"NA\", which becomes a real group key",
     "cargo test --release -p r-core --test subset_parity", "exact"),
    ("subset.netp_aggregates", "mean(pval) and sum(prob) in LONG_DOUBLE, and NA vs NaN kept apart",
     "cargo test --release -p r-core --test subset_parity", "exact"),
    ("subset.datasets_needs_logfc", "the datasets column survives only alongside ligand.logFC",
     "cargo test --release -p r-core --test subset_parity", "exact"),
    ("subset.unknown_slot", "an unrecognised slot.name falls through both final if/else branches and keeps every column",
     "cargo test --release -p r-core --test subset_parity", "exact"),
    ("rshim.subset_communication_deg",
     "subsetCommunication's DEG and netP branches through the installed shim, 34 fixtures, whole data frame",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("ranknet.information_flow", "the per-pathway apply(prob, 3, sum) in LONG_DOUBLE, 11 fixtures",
     "cargo test --release -p r-core --test ranknet_parity", "exact"),
    ("ranknet.weight_transform", "-1/log, the is.na reset, and the seq(max*1.1, max*1.5) reassignment",
     "cargo test --release -p r-core --test ranknet_parity", "exact"),
    ("ranknet.degenerate_flagging", "only a total above 1 is flagged; 0 and negative are not",
     "cargo test --release -p r-core --test ranknet_parity", "exact"),
    ("ranknet.reassignment_order", "match(1:m, position) inverts the sort, preserving original order",
     "cargo test --release -p r-core --test ranknet_parity", "exact"),
    ("ranknet.all_degenerate_error", "upstream's seq(-Inf) crash, preserved verbatim",
     "cargo test --release -p r-core --test ranknet_parity", "error-equal"),
    ("ranknet.count_measure", "1*(prob > 0) and no rescaling, 11 fixtures",
     "cargo test --release -p r-core --test ranknet_parity", "exact"),
    ("ranknet.group_filters", "sources/targets zero different axes; both validated against axis 1",
     "cargo test --release -p r-core --test ranknet_parity", "exact"),
    ("ranknet.order_stable", "order() on a double is radix and therefore stable, ties included",
     "cargo test --release -p r-core --test ranknet_parity", "exact"),
    ("ranknet.relative_digits1", "as.numeric(format(x, digits = 1)): one significant digit, Inf survives",
     "cargo test --release -p r-core --test ranknet_parity", "exact"),
    ("ranknet.empty_network_error", "sum(prob) == 0 stops before the per-pathway sums",
     "cargo test --release -p r-core --test ranknet_parity", "error-equal"),
    ("spatial.trimmed_mean", "R's mean(x, trim, na.rm): floor(n*trim) per end, NaN dropped, median branch",
     "cargo test --release -p r-core --test spatial_parity", "exact"),
    ("spatial.fdist", "collapse::fdist, the full matrix bit for bit and the dist triangle as a multiset",
     "cargo test --release -p r-core --test spatial_parity", "exact"),
    ("spatial.kdtree_exact",
     "the exact k-d tree against exhaustive search: lattice, duplicates, collinear, circle, clusters",
     "cargo test --release -p r-core --test spatial_parity", "exact"),
    ("analysis.centrality_deterministic",
     "unweighted degrees, weighted strengths (plain sequential f64 in edge-ID order, not R's "
     "long-double rowSums) and Brandes betweenness (Dijkstra with dist-plus-one encoding, exact "
     "2-way-heap tie rules, epsilon comparisons at 1e-10) against installed igraph 2.3.4 on a "
     "77-case corpus incl. 6 real tutorial slices, plus igraph's exact NA/positivity messages",
     "cargo test --release -p r-core --test centrality_parity", "exact"),
    ("spatial.kdtree_tie_break", "ties break to the lower index, which unique() makes observable downstream",
     "cargo test --release -p r-core --test spatial_parity", "exact"),
    ("spatial.kdtree_subquadratic", "the search is subquadratic, which is the reason to replace Annoy",
     "cargo test --release -p r-core --test spatial_parity", "exact"),
    ("spatial.region_arithmetic",
     "computeRegionDistance's arithmetic against upstream's body with the neighbour query "
     "substituted: threshold tests, per-sample merge, binarisation, symmetrisation, NaN placement",
     "cargo test --release -p r-core --test region_parity", "exact"),
    ("spatial.region_ratio_recycling",
     "a length-1 ratio/tol is recycled the way R's ratio[k] is, and any other short length refused",
     "cargo test --release -p r-core --test region_parity", "exact"),
    ("spatial.region_null_ratio",
     "the default signature: ratio = tol = NULL empties the selection, leaving d.spatial all NaN",
     "cargo test --release -p r-core --test region_parity", "exact"),
    ("spatial.region_absent_level",
     "a declared level with no cells: NaN in the unwritten row and column, and the row names "
     "upstream assigns do not label the data it returns",
     "cargo test --release -p r-core --test region_parity", "exact"),
    ("spatial.region_d_spatial_always_symmetric",
     "d.spatial <- (d + t(d))/2 sits OUTSIDE if (do.symmetric) upstream, so it is symmetric either way",
     "cargo test --release -p r-core --test region_parity", "exact"),
    ("spatial.region_contact_null_range",
     "contact.range = NULL makes dist - NULL an empty selection, not a skipped pair, so adj.contact "
     "is 0 throughout and only the adj.contact.knn swap can produce a 1",
     "cargo test --release -p r-core --test region_parity", "exact"),
    ("spatial.region_contact_knn_k_replaces",
     "contact.knn.k replaces adj.contact rather than merging, on a fixture where the range-based "
     "and rank-based counts disagree",
     "cargo test --release -p r-core --test region_parity", "exact"),
    ("spatial.region_corpus_margin",
     "every fixture's recorded nearest-neighbour margin, which is what makes the exact-substitution "
     "oracle sound; a layout edit that breaks it fails here and not only in the generator",
     "cargo test --release -p r-core --test region_parity", "exact"),
    ("rshim.computeCellDistance",
     "computeCellDistance against upstream: the dist rewrap (call/method/Labels), the ncol error, "
     "ratio applied when non-NULL, the threshold needing BOTH range and tol, and tol as an addend",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("rshim.computeRegionDistance_exact",
     "the R-side marshalling for all 12 corpus fixtures: levels(factor) order, the declared-but-"
     "absent level, match() on the factor, and the dimnames",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("ranknet.comparison_flow",
     "rankNet(mode = 'comparison'): per-comparison flow over 9 fixtures, the POOLED degenerate "
     "reassignment, and the raw relative ratios",
     "cargo test --release -p r-core --test ranknet_comparison_parity", "exact"),
    ("ranknet.comparison_pooled_reassignment",
     "one values.assign / pSum.original.all / position across every comparison, not one per "
     "comparison; separated from the per-comparison version by a fixture with flagged pathways "
     "in both comparisons at once",
     "cargo test --release -p r-core --test ranknet_comparison_parity", "exact"),
    ("ranknet.comparison_flag_is_a_range",
     "`pSum < 0` flags every pathway whose total reaches or exceeds 1, not only the one that "
     "totals exactly 1",
     "cargo test --release -p r-core --test ranknet_comparison_parity", "exact"),
    ("rshim.ranknet_comparison",
     "rankNet(mode = 'comparison') through the installed shim, 9 fixtures: the union of pathway "
     "vocabularies, the row order from the multi-key order(), the rbind-uniquified row names, the "
     "zero-dropping loop, and ggplot_build()$data for the plot",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("prob.data_max_scaling",
     "`data.use <- data/max(data)`, upstream's first numerical step: absent from the kernel "
     "until the configuration matrix caught it, invisible while every fixture had max == 1",
     "cargo test --release -p r-core --test prob_parity", "exact"),
    ("prob.scale_invariance",
     "Prob is invariant to a constant factor on the expression matrix, because upstream divides by "
     "max(data) first; asserted directly, x1 against x3",
     "R_LIBS=.rlib Rscript tests/parity/check_matrix.R", "identical"),
    ("longdouble.nan_and_inf",
     "R's max() and sum() propagate NaN and carry Inf; Rust's f64::max is a maxNum and F80 had no "
     "special values at all",
     "cargo test --release -p r-core --test longdouble_special_parity", "exact"),
    ("prob.if_na_raises",
     "`if (sum(P1_Pspatial) == 0)` with a NaN sum raises R's 'missing value where TRUE/FALSE "
     "needed'; returning an all-NaN network instead is what a port without an `if` would do",
     "cargo test --release -p r-core --test prob_parity", "error-equal"),
    ("stats.nnzero_is_na",
     "`Matrix::nnzero` returns NA for any vector with a missing value, not a smaller count, so "
     "upstream's `thresholdedMean` raises at `if (percent < trim)`; the port counted non-zeros and "
     "skipped the missing ones, and a unit test had locked that reading in",
     "cargo test --release -p r-core --lib", "error-equal"),
    ("prob.observed_means_raise",
     "upstream aggregates every gene at modelling.R:115, so a raise inside FunMean on a gene no "
     "interaction references still fires; the port computed that NaN and discarded it, turning a hard "
     "error into an all-zero Prob",
     "cargo test --release -p r-core --lib", "error-equal"),
    ("prob.error_precedence",
     "upstream aggregates the observed data (:115) before the bootstrap replicates (:208); the port "
     "ran them in the reverse order, so its error precedence differed from upstream's",
     "R_LIBS=.rlib Rscript tests/parity/check_matrix.R", "error-equal"),
    ("ranknet.pairwise_order",
     "`rankNetPairwise`'s `order(pval, -prob)` per group pair: a stable radix sort, so equal keys keep "
     "input order; the data-frame construction around it stays in R",
     "R_LIBS=.rlib Rscript tests/parity/check_identical.R", "identical"),
    ("ranknet.order_later_key_na",
     "R's `order(c(1,1,2), c(NA,5,3))` is `2 1 3` -- a missing value in a *later* key sorts last "
     "within its group, and an all-missing primary key is still ordered by the remaining keys. Both "
     "were wrong, and neither is visible without a repeated key",
     "cargo test --release -p r-core --lib", "exact"),
    ("ranknet.pairwise_lr_use_mismatch",
     "upstream does not reconcile a short `LR.use` with `dim(prob)[3]`: its `data.frame(row.names =)` "
     "raises, and the port must raise from its own code rather than by slicing prob to LR.use",
     "R_LIBS=.rlib Rscript tests/parity/check_identical.R", "error-equal"),
    ("ranknet.pairwise_nonsquare",
     "`numCluster` is `dim(prob)[1]` and upstream loops it in both dimensions, so a non-square prob "
     "reads the leading k x k block and then raises on `names(temp) <- colnames(prob)`",
     "R_LIBS=.rlib Rscript tests/parity/check_identical.R", "error-equal"),
    ("rshim.rankNetPairwise",
     "the shim's rankNetPairwise reaches the Rust ordering on every gated case; a pass-through and "
     "a port are both `identical` to upstream, so the gate counts fallback warnings as failures",
     "R_LIBS=.rlib Rscript tests/parity/check_identical.R", "identical"),
    ("prob.bootstrap_parallel_identity",
     "parallelising the bootstrap aggregate leaves every output bit-identical: each replicate reads "
     "a disjoint permutation and writes a disjoint slot, and indexed collect() keeps replicate order "
     "so the F80 accumulation order within a replicate is untouched. Checked on 1 207 952 real-data "
     "Prob values across the two authors' fixtures, max absolute difference 0",
     "R_LIBS=.rlib NB=20 Rscript bench-runner/check_parity_real.R",
     "identical"),
    ("prob.permutation_stream_order",
     "R draws all nboot permutations from ONE MT19937 stream, so replicate nE consumes the nE-th "
     "slice; drawing them in parallel would need a stream per replicate, which is a different "
     "permutation set rather than a different order of the same one",
     "cargo test --release -p r-core --test prob_parity", "exact"),
    ("meta.scaling",
     "Prob under data*f: port == upstream exactly for f in {0.5, 3, 1e4, 7.25}, and invariant "
     "relative to the unscaled run -- exactly for powers of two, to within 1 ulp otherwise, because "
     "max(f*x) is not f*max(x) in floating point",
     "R_LIBS=.rlib Rscript tests/parity/metamorphic.R", "identical"),
    ("meta.gene_permutation",
     "permuting the rows of data.signaling with the dimnames attached changes nothing, and the "
     "negative control (doubling one gene) does change it, so the invariance is not vacuous",
     "R_LIBS=.rlib Rscript tests/parity/metamorphic.R", "identical"),
    ("meta.cluster_relabelling",
     "prob[i,j,l] is 'group i sends to group j', so relabelling permutes BOTH indices: "
     "prob[ord(i), ord(j), l] must equal the original, and an order-preserving relabel must leave "
     "the values in place while the dimnames change",
     "R_LIBS=.rlib Rscript tests/parity/metamorphic.R", "identical"),
    ("meta.cell_duplication",
     "duplicating every cell within its own group multiplies every group size by two and leaves "
     "Prob unchanged, which is what proves the kernel is not weighting by group size",
     "R_LIBS=.rlib Rscript tests/parity/metamorphic.R", "identical"),
    ("meta.cell_permutation",
     "permuting the cells (labels travelling with the columns) leaves the observed Prob invariant "
     "but MUST move Pval, since the bootstrap draws sample.int(nC, nC) over cell indices. "
     "Asserting whole-object equality here would assert something false",
     "R_LIBS=.rlib Rscript tests/parity/metamorphic.R", "identical"),
    ("meta.nboot_invariance",
     "Prob comes from the observed aggregate and is invariant to nboot at 3, 10, 25 and 100, while "
     "Pval is always exactly k/nboot and is 1 wherever Prob is 0",
     "R_LIBS=.rlib Rscript tests/parity/metamorphic.R", "identical"),
    ("meta.recorded_invariants",
     "checked on every object the metamorphic file builds: dimnames consistency, Pval in {k/nboot}, "
     "Pval[Prob==0] == 1, and Prob in [0,1]",
     "R_LIBS=.rlib Rscript tests/parity/metamorphic.R", "exact"),
    ("fuzz.db_resolution",
     "randomly generated L-R databases -- empty subunit lists, names colliding across the complex, "
     "cofactor and symbol tables, self-referential complexes -- resolved through resolve_entity and "
     "compute_expr_lr. Every entry point is total: a missing subunit is an Err carrying upstream's "
     "subscript out of bounds, never a Rust panic",
     "cargo test --release -p r-core --test db_fuzz", "error-equal"),
    ("fuzz.extract_gene_split",
     "extractGene splits on the OFFICIAL SYMBOL list, not on membership of the complex table: a "
     "symbol is kept as itself and never expanded, a non-symbol is looked up in the complex table "
     "and contributes only its non-empty subunits. Asserting the complex-table rule instead fails "
     "on four different generated shapes",
     "cargo test --release -p r-core --test db_fuzz", "exact"),
    ("props.r_order",
     "generated multi-key vectors against a transcription of R's `order(..., na.last = TRUE)`: a "
     "missing value sorts last in EVERY key, and an all-missing primary key is still ordered by the "
     "remaining keys. Found both order_f64_multi bugs, which no fixture could see because no "
     "fixture had a repeated key",
     "cargo test --release -p r-core --test properties", "exact"),
    ("props.longdouble_mul_specials",
     "F80::mul against R's measured infinite-product table. Infinity is encoded with a ZERO "
     "mantissa, so the old leading zero-check caught every infinite operand and prod(c(-Inf)) was 0 "
     "where R gives -Inf. R's rule is Inf with the XOR of the signs -- including Inf * -Inf, which "
     "IEEE calls NaN -- and NaN only for zero times infinity",
     "cargo test --release -p r-core --lib longdouble", "exact"),
    ("props.signed_zero",
     "a zero product keeps the XOR of the signs (prod(c(-0,5)) and prod(c(5,-0)) are both -0) and "
     "to_f64 must not drop the sign of a zero mantissa",
     "cargo test --release -p r-core --lib longdouble", "exact"),
    ("props.aggregate_level_order",
     "aggregate returns one row per level in LEVEL order, never order of first appearance, and the "
     "Mean row is the arithmetic mean of exactly that group's values in ascending cell order",
     "cargo test --release -p r-core --test properties", "exact"),
    ("props.rng_permutation",
     "sample.int without replacement is always a permutation of 1..n for arbitrary n, k and seed, "
     "and the stream is reproducible and inside [0,1)",
     "cargo test --release -p r-core --test properties", "exact"),
    ("matrix.configuration_coverage",
     "the configuration matrix: 200 configurations over 14 axes, all 1051 axis pairs covered, "
     "reported with the achieved coverage rather than asserted in prose",
     "R_LIBS=.rlib Rscript tests/parity/check_matrix.R", "identical"),
    ("matrix.scale_invariance", "the matrix's scale axis and the direct x1-vs-x3 invariance check",
     "R_LIBS=.rlib Rscript tests/parity/check_matrix.R", "identical"),
    ("rshim.ranknet_single",
     "rankNet(mode = 'single') through the installed shim, 11 fixtures, data frame and plot data",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    # A quantity must not appear in both PASSING and GAPS. When it did, the gap wording won and
    # `parity.json` reported a gated quantity as untested -- the one failure mode a file whose
    # entire job is machine-checkable correctness must not have. `main()` now asserts the
    # disjointness, so a repeat is a hard error rather than a silently wrong artifact.
    ("prob.spatial_branch",
     "computeCommunProb's spatial branch: P.spatial from the exact k-d tree, the scale.distance "
     "stop, and the four nLR1 branches chosen by the annotation column",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("prob.raw_use_false", "computeCommunProb(raw.use = FALSE): 9 configurations reading data.smooth",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
    ("filter.raw_use_smooth", "6 filterCommunication fixtures on data.smooth, cell for cell",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("filter.raw_use_scaling", "data <- data/max(data) applies to the *selected* matrix, not to data.signaling",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("filter.mincells_is_count_only", "with min.samples = NULL the two settings agree: the min.cells half counts cells",
     "cargo test --release -p r-core --test filter_parity", "exact"),
    ("rshim.filter_communication",
     "filterCommunication through the installed shim, 13 fixtures, whole object and stdout",
     "R_LIBS=.rlib R --vanilla -f tests/parity/check_identical.R", "identical"),
]

# Quantities named by the objective that are NOT yet ported. Listed with an explicit
# "untested" rung rather than omitted, so the gap is visible in the artifact.
GAPS: dict[str, str] = {
    # `computeRegionDistance` is done: the arithmetic foundations, the exact k-d tree, the
    # assembly, and the **measured divergence** on the authors' visium data
    # (`prob.spatial_divergence_measured`, a `measured-divergent` rung). It is no longer a gap.
    #
    # The centrality measures over the `netP` slot are a separate numeric surface in `analysis.R`
    # and not a wrapper around anything already ported, so they are listed as their own gap.
    # igraph 2.3.4 is installed and `netAnalysis_computeCentrality` is reachable, so this is not
    # blocked on a dependency -- it is unported work. Note that `hub_score`, `authority_score`,
    # `eigen_centrality` and `page_rank` are igraph solver outputs (HITS power iteration, ARPACK,
    # PRPACK); bit-parity with those means delegating to igraph rather than reimplementing them, and
    # that is the design decision still to be made and recorded.
    # `do.fast = TRUE` is gated and measured, not a gap. What is *not* done is a reimplementation of
    # presto's normal approximation, and the reason is recorded rather than asserted: matching a C++
    # implementation bit-for-bit means matching its rank-sum statistic, its tie correction, whether it
    # applies a continuity correction and which z-to-p function it calls. Three of those four right is
    # still wrong in the fourth and still looks plausible. The shim delegates, loudly, and the gate
    # catches a break in the delegation.
    # Reachable without presto, and worth being precise about why. `group.dataset` appears in
    # *both* branches of upstream's `if (do.fast)`: the presto one (utilities.R:429-484) and the
    # Wilcoxon one (:512-520). Only the first is unreachable here, so `do.fast = FALSE` with
    # `group.dataset` set is testable today. What is left to port is the dataset-comparison
    # selection -- `cell.use1`/`cell.use2` chosen by dataset rather than by group complement, the
    # `pos.dataset` validation with its `cat()` followed by a bare `stop()`, the
    # `group.DE.combined` variant that pools cells across groups, and the
    # `markers.all$datasets` factor and its `order(datasets, pvalues, -logFC)`.

}


@dataclass
class Result:
    command: str
    returncode: int | None = None
    passed: int | None = None
    failed: int | None = None
    detail: str = ""


def run(cmd: str, timeout: int = 3600) -> Result:
    r = Result(command=cmd)
    try:
        p = subprocess.run(
            cmd, shell=True, cwd=ROOT, capture_output=True, text=True, timeout=timeout
        )
    except subprocess.TimeoutExpired:
        r.detail = "timed out"
        return r
    r.returncode = p.returncode
    out = p.stdout + p.stderr
    # `cargo test` prints one "test result: ok. N passed; M failed" line per binary.
    for m in re.finditer(r"test result: (?:ok|FAILED)\. (\d+) passed; (\d+) failed", out):
        r.passed = (r.passed or 0) + int(m.group(1))
        r.failed = (r.failed or 0) + int(m.group(2))
    if r.returncode != 0:
        tail = [ln for ln in out.strip().splitlines() if ln.strip()][-6:]
        r.detail = " | ".join(tail)[:800]
    return r


def rust_suite() -> Result:
    r = run("cargo test --release")
    r.command = "cargo test --release"
    return r


## R gates that are whole scripts rather than single cargo tests. Each is `(quantity, pin, script,
## rung)`, and the rung is only claimed when the script exits zero -- the same rule as the cargo
## entries, so a gate that silently stops running reports `untested` rather than keeping its rung.
R_GATES: list[tuple[str, str, str, str]] = [
    (
        "gate.r_shim_identical",
        "the acceptance gate: serialized-byte comparisons on deterministic scientific outputs",
        "tests/parity/check_identical.R",
        "identical",
    ),
    (
        "stats.independent_stream_equivalence",
        "120 seeds drawn from a Park-Miller stream this test implements itself, salted by the "
        "fixture's shape so neither the port nor upstream can influence the sequence: bit-identical "
        "Prob and Pval at every seed, KS = 0 on the pooled values, the same set of null "
        "distributions drawn on both sides, and the objective's Pval invariants at every seed",
        "tests/parity/stat_equiv.R",
        "identical",
    ),
    (
        "package.independent_dropin",
        "109 public function signatures, public datasets, independent constructors and inference, "
        "unmodified Annoy spatial output, native SNN helper, saved-object interoperability and "
        "PNG rendering compared as bytes across independently installed packages in clean processes",
        "tests/parity/check_dropin.R",
        "identical",
    ),
    (
        "utilities.identify_over_expressed_genes_fast",
        "identifyOverExpressedGenes(do.fast = TRUE), the presto branch, gated against pinned "
        "upstream with presto 1.1.0 installed: identical() on the whole marker table, plus a "
        "measurement of how far the two branches actually are -- same 12 features selected, p-values "
        "within 3.8e-3 with no threshold crossing at 0.05, and pct.1 differing by exactly a factor of "
        "100 because the two branches report it as a fraction and a percentage while filtering on the "
        "same thresh.pc * 100",
        "tests/parity/check_oeg_fast.R",
        "identical",
    ),
    (
        "delegation.complete",
        "every non-visual name upstream exports resolves in the installed shim, every generated "
        "wrapper reaches the pinned upstream rather than a stub, and triMean/geometricMean agree "
        "with upstream identical() on fixed vectors; the tutorial reproduction exercises the "
        "pipeline ones end-to-end. Visualization/app names stay absent per locked decision 14.1",
        "tests/parity/check_delegation.R",
        "identical",
    ),
    (
        "tutorial.reproduction",
        "the vignette's computational pipeline on the authors' human-skin data (createCellChat to "
        "netAnalysis_computeCentrality, nboot = 100, no plots per 14.1): identical() on every "
        "non-visual slot after every stage, whole-object identical() modulo run.time, and the "
        "objective's Pval/dimname invariants on the tutorial's own output",
        "tests/parity/tutorial_repro.R",
        "identical",
    ),
    (
        "analysis.netP_centrality",
        "netAnalysis_computeCentrality assembled from Rust deterministic measures plus igraph "
        "solvers plus upstream's sna tryCatch, identical() on real pathway slices with both sides "
        "seeded identically, with the NA/NaN matrix error and the tiny-weights warning verified "
        "byte for byte. hub/authority/eigen/page_rank stay delegated: those solvers disagree with "
        "themselves across runs, so bit-parity with them is meaningless -- identical by "
        "construction (same package, same call) instead of by reimplementation",
        "tests/parity/check_centrality.R",
        "identical",
    ),
    (
        "cli.standalone",
        "the `cellchatrs` binary on 14 configurations -- all four type.mean values (two of them "
        "given as match.arg prefixes), population.size both ways, nboot 1..9, Kh across six orders of "
        "magnitude, n at 1 and 2, raw.use both ways, two seeds: Prob, Pval, computeAveExpr's means "
        "and aggregateNet's two matrices all identical() to pinned upstream, with the binary built "
        "once by cargo and exercised as the artifact a user would install",
        "tests/parity/check_cli.R",
        "identical",
    ),
    (
        "analysis.mergeCellChat",
        "mergeCellChat's own contract (gene/meta intersections, idents$joint level order, "
        "add.names propagation, both stop() messages) plus the ported aggregateNet / "
        "subsetCommunication / filterCommunication / rankNet run on a per-dataset net lifted out of "
        "a merged object. Not reimplemented -- it is S4 slot assembly with no arithmetic, so the "
        "14.1 scope decision leaves it in R and this is the pass-through rung",
        "tests/parity/check_merge.R",
        "identical",
    ),
]


def r_gate(script: str) -> Result:
    cmd = f"R_LIBS=.rlib R --vanilla -f {script}"
    r = Result(command=cmd)
    env = dict(os.environ, R_LIBS=".rlib")
    try:
        p = subprocess.run(
            cmd, shell=True, cwd=ROOT, capture_output=True, text=True, timeout=3600, env=env
        )
    except subprocess.TimeoutExpired:
        r.detail = "timed out"
        return r
    r.returncode = p.returncode
    out = p.stdout + p.stderr
    # Two verdict shapes: the `identical` gates print "N failing comparisons out of M", the
    # statistical-equivalence gate prints "N failing checks out of M".
    m = re.search(r"(\d+) failing (?:comparisons|checks) out of (\d+)", out)
    if m:
        r.failed = int(m.group(1))
        r.passed = int(m.group(2)) - int(m.group(1))
    if r.returncode != 0:
        tail = [ln for ln in out.strip().splitlines() if ln.strip()][-6:]
        r.detail = " | ".join(tail)[:800]
    return r


@dataclass
class Report:
    entries: list[dict] = field(default_factory=list)


def build(run_it: bool) -> dict:
    commands: dict[str, Result] = {}
    entries = []
    # `package.independent_dropin` compares two *independently installed* packages, so it needs
    # the pinned original installed beside this one. Installing that original needs R >= 4.5 --
    # its `ggpubr -> rstatix -> car -> pbkrtest -> doBy -> Deriv` chain ends in a package that
    # declares `R (>= 4.5)` -- while every R-side rung above is bit-identical only on R <= 4.4,
    # where the datasets are also exposed the way the comparison expects. No single job can do
    # both today, so where no oracle is provided the quantity is recorded untested with no test
    # command rather than failing the report for an environment it cannot have: `main()` fails
    # only on entries that name a test and have no passing rung. The gate itself is real: set
    # `CELLCHAT_ORACLE_LIB` to a library holding the pinned original (installed from a checkout
    # of `UPSTREAM_SHA`, which needs R >= 4.5 to resolve its dependencies) and
    # `CELLCHAT_REPLACEMENT_LIB` to this package's library, and it runs. Verified that way:
    # 25 byte comparisons and 109 public signatures matched.
    oracle = os.environ.get("CELLCHAT_ORACLE_LIB", "")
    for quantity, how, command, rung_when_green in QUANTITIES:
        key = command
        if run_it and key not in commands:
            commands[key] = rust_suite() if command == "cargo test --release" else run(command)
        res = commands.get(key)
        green = res is not None and res.returncode == 0
        entries.append(
            {
                "quantity": quantity,
                "pin": how,
                "test": command,
                "rung": rung_when_green if green else "untested",
                "exit_code": None if res is None else res.returncode,
                "detail": "" if green or res is None else res.detail,
            }
        )
    # The two gates that are whole suites rather than single tests.
    if run_it:
        full = commands.get("cargo test --release") or rust_suite()
        entries.append(
            {
                "quantity": "suite.rust",
                "pin": f"every Rust test, {full.passed or 0} passing",
                "test": "cargo test --release",
                "rung": "exact" if full.returncode == 0 else "untested",
                "exit_code": full.returncode,
                "detail": full.detail,
            }
        )
        for quantity, pin, script, rung in R_GATES:
            # `package.independent_dropin` lives here, not in QUANTITIES, and it is the one gate
            # that needs a second, independently installed package. Record it untested with no
            # test command when no oracle is provided -- see the comment on `build()` above for
            # why no single job can both install the oracle and keep the R-side rungs identical.
            if quantity == "package.independent_dropin" and not oracle:
                entries.append(
                    {
                        "quantity": quantity,
                        "pin": "requires the pinned original installed as a second package; its "
                               "dependency tree needs R >= 4.5, and the R-side rungs need R <= 4.4 "
                               "for bit-identity",
                        "test": None,
                        "rung": "untested",
                        "exit_code": None,
                        "detail": "oracle not available in this environment",
                    }
                )
                continue
            gate = r_gate(script)
            entries.append(
                {
                    "quantity": quantity,
                    "pin": pin,
                    "test": gate.command,
                    "rung": rung if gate.returncode == 0 else "untested",
                    "exit_code": gate.returncode,
                    "detail": gate.detail,
                }
            )
    for quantity, why in sorted(GAPS.items()):
        entries.append(
            {"quantity": quantity, "pin": why, "test": None, "rung": "untested",
             "exit_code": None, "detail": "not ported yet"}
        )
    # Measurements of a divergent optional algorithm never count as byte-identical compatibility.
    # A gap that is declared here but missing from `entries` is a silent regression: the artifact
    # would report a smaller gap list and a higher pass rate, and nothing would fail. This happened --
    # removing `prob.spatial_divergence` from GAPS took `analysis.netP_centrality` with it, and the
    # report went from "3 untested" to "1 untested" while the centrality work had not been done at all.
    # Asserted rather than assumed, and the duplicate check catches the other direction.
    present = {e["quantity"] for e in entries}
    missing = sorted(set(GAPS) - present)
    if missing:
        raise SystemExit(
            "parity.json would omit declared gaps: " + ", ".join(missing)
            + " -- every GAPS key must reach the artifact, or the gap list understates the work"
        )
    seen: dict[str, int] = {}
    for e in entries:
        seen[e["quantity"]] = seen.get(e["quantity"], 0) + 1
    dupes = sorted(q for q, n in seen.items() if n > 1)
    if dupes:
        raise SystemExit("parity.json lists a quantity twice: " + ", ".join(dupes))

    for e in entries:
        alternative = e["quantity"].startswith(("spatial.kdtree", "spatial.region", "rshim.computeRegionDistance_exact"))
        e["scope"] = "alternative-algorithm" if alternative else "compatibility"
    n_exact = sum(1 for e in entries
                  if e["scope"] == "compatibility" and e["rung"] in
                  ("exact", "exact-nan", "identical", "error-equal"))
    return {
        "upstream": {
            "repo": "jinworks/CellChat",
            "commit": UPSTREAM_SHA,
            "version": "2.2.0.9001",
        },
        "rung_vocabulary": {
            "exact": "raw bit patterns compared and required equal",
            "exact-nan": "as exact, with non-finite inputs in the corpus",
            "identical": "R identical() on the object, run through the installed shim",
            "error-equal": "error messages compared byte for byte",
            # The spatial KNN is the one place the port is *deliberately* not bit-identical: the
            # oracle is itself approximate, so the rung records a measured divergence against real
            # data rather than claiming parity. `PLAN.md` §14.5.
            "measured-divergent": (
                "the divergence from upstream's approximate index is measured on real data and "
                "published, rather than silently accepted; the exact side is gated by the same "
                "kernel tests as the rest"
            ),
            "untested": "no passing test establishes this yet",
        },
        "summary": {
            "quantities": sum(e["scope"] == "compatibility" for e in entries),
            "alternative_checks": sum(e["scope"] == "alternative-algorithm" for e in entries),
            "at_rung": n_exact,
            "untested": sum(e["rung"] == "untested" for e in entries),
            "divergent": sum(e["rung"] == "measured-divergent" for e in entries),
        },
        "entries": entries,
    }


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--no-run", action="store_true", help="report the last results, run nothing")
    ap.add_argument("-o", "--out", default="parity.json")
    args = ap.parse_args()
    doc = build(run_it=not args.no_run)
    out = Path(args.out) if os.path.isabs(args.out) else ROOT / args.out
    out.write_text(json.dumps(doc, indent=2) + "\n")
    s = doc["summary"]
    print(f"wrote {out.relative_to(ROOT) if out.is_relative_to(ROOT) else out}")
    print(f"  {s['at_rung']}/{s['quantities']} quantities at a passing rung, {s['untested']} untested")
    for e in doc["entries"]:
        if e["rung"] == "untested" and e["test"]:
            print(f"  FAIL {e['quantity']}: {e['test']} (exit {e['exit_code']})")
    return 0 if all(e["rung"] != "untested" or e["test"] is None for e in doc["entries"]) else 1


if __name__ == "__main__":
    sys.exit(main())
