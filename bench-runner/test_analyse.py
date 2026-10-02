#!/usr/bin/env python3
"""Self-tests for `bench-runner/analyse.py`.

The benchmark analysis step had no tests at all, and that is how it went broken in a way nothing
noticed: `analyse.py` assumed the fixture metadata lived under a `fixture:` key, the way
`bench_real.R` writes it. `bench_more_datasets.R` puts those fields at the top level and
`bench_visium.R` calls cells `n_spots`, so the moment either of those runners wrote a file into
`bench-runner/results/` the documented `python3 bench-runner/analyse.py` step raised `KeyError`.
`tables.md` stayed a day out of date while every runner reported green, because nothing read the
step's output.

Stdlib `unittest` only, for the same reason `analyse.py` hand-writes its SVG and `bench_real.R`
hand-writes its JSON: the pinned host has no `pytest` and no `jsonlite`, and a benchmark artefact
that cannot be regenerated without installing something is an artefact nobody regenerates.

Run:  python3 bench-runner/test_analyse.py
"""
from __future__ import annotations

import importlib.util
import io
import json
import contextlib
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

HERE = Path(__file__).resolve().parent

_spec = importlib.util.spec_from_file_location("analyse", HERE / "analyse.py")
analyse = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(analyse)


def side(**over):
    """A `bench_real.R`-shaped timing block."""
    d = {"raw_secs": [1.0, 1.1, 0.9, 1.05, 0.95],
         "median": 1.0, "ci_lo": 0.9, "ci_hi": 1.1, "min": 0.9, "max": 1.1}
    d.update(over)
    return d


class TestNormFixture(unittest.TestCase):
    """Three runners, three schemas, one reader."""

    def test_bench_real(self):
        got = analyse.norm_fixture({"fixture": {"fixture": "humanSkin.rda", "n_cells": 7563,
                                                "n_genes_signaling": 1094, "n_groups": 12,
                                                "n_lr": 1583}})
        self.assertEqual(got["name"], "humanSkin.rda")
        self.assertEqual(got["n_cells"], 7563)
        self.assertEqual(got["kind"], "scrnaseq")
        self.assertEqual(got["n_lr"], 1583)

    def test_bench_more_datasets_is_flat(self):
        """This is the schema that raised `KeyError: 'fixture'`."""
        got = analyse.norm_fixture({"dataset": "nl", "species": "human", "n_cells": 2552,
                                    "n_genes_signaling": 593, "n_groups": 12, "n_lr": 492})
        self.assertEqual(got["name"], "nl")
        self.assertEqual(got["n_cells"], 2552)
        self.assertEqual(got["kind"], "scrnaseq")
        self.assertEqual(got["n_genes_signaling"], 593)

    def test_bench_visium_spatial(self):
        got = analyse.norm_fixture({"fixture": {"dataset": "visium", "n_spots": 1073,
                                                "n_genes_signaling": 433, "n_groups": 8,
                                                "n_lr": 134}})
        self.assertEqual(got["name"], "visium")
        self.assertEqual(got["n_cells"], 1073)
        self.assertEqual(got["kind"], "spatial")

    def test_missing_counts_do_not_raise(self):
        """A runner that grows a field must not be able to break the reader."""
        got = analyse.norm_fixture({"dataset": "x", "n_cells": 10})
        self.assertIsNone(got["n_lr"])
        self.assertEqual(got["n_cells"], 10)


class TestFitAffine(unittest.TestCase):
    def test_exact_line(self):
        a, b, r2 = analyse.fit_affine([1.0, 2.0, 3.0], [3.0, 5.0, 7.0])
        self.assertAlmostEqual(a, 1.0)
        self.assertAlmostEqual(b, 2.0)
        self.assertAlmostEqual(r2, 1.0)

    def test_flat(self):
        a, b, r2 = analyse.fit_affine([1.0, 2.0, 3.0], [4.0, 4.0, 4.0])
        self.assertAlmostEqual(b, 0.0)
        self.assertAlmostEqual(a, 4.0)

    def test_too_few_is_nan_not_a_number(self):
        self.assertTrue(all(v != v for v in analyse.fit_affine([1.0], [2.0])))


class TestEndToEnd(unittest.TestCase):
    """Run the real script over a directory holding all three schemas at once."""

    def setUp(self):
        self.dir = Path(tempfile.mkdtemp())
        res = self.dir
        # one per schema, plus a thread sweep and the files the reader must skip
        (res / "human_skin.json").write_text(json.dumps({
            "schema": 1, "nboot": 100, "repeats": 5, "seed": 1,
            "fixture": {"fixture": "humanSkin.rda", "n_genes_total": 17328, "n_cells": 7563,
                        "n_groups": 12, "n_genes_signaling": 1094, "n_lr": 1583},
            "r_upstream": side(median=93.0), "rust": side(median=2.2),
            "speedup": 42.4, "speedup_ci_lo": 41.4, "speedup_ci_hi": 45.1,
            "peak_rss_mb": 1250.0, "parity_net": True, "parity_prob": True, "max_abs_diff": 0.0,
            "scatter": {"log10_prob_r": [0.5, 1.0, 1.5],
                        "log10_prob_rs": [0.5, 1.0, 1.5]}}))
        for stem, n, med in (("dataset_nl", 2552, 0.48), ("dataset_ls", 5011, 0.66)):
            (res / f"{stem}.json").write_text(json.dumps({
                "schema": 1, "nboot": 100, "repeats": 5, "seed": 1,
                "dataset": stem.split("_")[1], "species": "human", "n_cells": n,
                "n_genes_total": 1000, "n_genes_signaling": 593, "n_groups": 12,
                "n_lr": 500, "r_upstream": side(median=30.0), "rust": side(median=med),
                "speedup": 60.0, "speedup_ci_lo": 55.0, "speedup_ci_hi": 65.0,
                "peak_rss_mb": 550.0, "parity_net": True, "parity_prob": True,
                "max_abs_diff": 0.0}))
        (res / "visium.json").write_text(json.dumps({
            "schema": 1, "nboot": 100, "repeats": 5, "seed": 1,
            "fixture": {"dataset": "visium", "n_spots": 1073, "n_genes_total": 648,
                        "n_genes_signaling": 433, "n_groups": 8, "n_lr": 134,
                        "n_lr_total": 437},
            "r_upstream": side(median=17.1), "rust": side(median=0.30),
            "speedup": 57.2, "speedup_ci_lo": 54.3, "speedup_ci_hi": 62.0,
            "peak_rss_mb": 553.0, "parity_net": False, "parity_prob": False,
            "max_abs_diff": 0.00177,
            "scatter": {"log10_prob_r": [-1.0, -1.1, -2.0],
                        "log10_prob_rs": [-1.02, -1.09, -2.0]}}))
        for t in (1, 2, 4, 8):
            (res / f"sweep_t{t}.json").write_text(json.dumps({
                "threads": t, "nboot": 100, "repeats": 5,
                "rust": side(median=10.0 / t), "rust_only": True}))
        (res / "scaling.json").write_text("{}")
        (res / "spatial_divergence.json").write_text("{}")

    def run_script(self, *extra):
        out = Path(tempfile.mkdtemp())
        p = subprocess.run([sys.executable, str(HERE / "analyse.py"),
                            "--results", str(self.dir), "--outdir", str(out), *extra],
                           capture_output=True, text=True)
        return p, out

    def test_runs_over_every_schema(self):
        p, out = self.run_script()
        self.assertEqual(p.returncode, 0, p.stderr)
        md = (out / "tables.md").read_text()
        for row in ("humanSkin.rda", "nl", "ls", "visium"):
            self.assertIn(row, md, f"{row} missing from the table")
        self.assertIn("57.20", md)  # nothing was reformatted into a different number

    def test_spatial_divergence_is_not_labelled_identical(self):
        """The one row that must never read `identical` is the exact-KNN one."""
        p, out = self.run_script()
        md = (out / "tables.md").read_text()
        visium = [ln for ln in md.splitlines() if "| visium |" in ln][0]
        self.assertNotIn("identical", visium)
        self.assertIn("divergence by design", visium)

    def test_fit_over_three_or_more_points_only(self):
        p, out = self.run_script()
        fit = json.loads((out / "scaling.json").read_text())["scaling"]["fixed_nboot_100_fit"]
        self.assertEqual(fit["n_points"], 3)
        self.assertNotIn("status", fit)
        self.assertIn("r2", fit)
        # visium is spatial and must not be in a scRNA-seq cost fit
        self.assertIn("visium", fit["excluded_spatial"])
        self.assertNotIn("visium", [q["fixture"] for q in fit["points"]])

    def test_two_points_withhold_the_fit(self):
        """A line through two points has R^2 == 1 by construction; it must not be reported."""
        for stem in ("dataset_nl.json", "dataset_ls.json"):
            (self.dir / stem).unlink()
        (self.dir / "visium.json").unlink()
        p, out = self.run_script()
        self.assertEqual(p.returncode, 0, p.stderr)
        fit = json.loads((out / "scaling.json").read_text())["scaling"]["fixed_nboot_100_fit"]
        self.assertEqual(fit["n_points"], 1)
        self.assertIn("withheld", fit["status"])
        self.assertNotIn("r2", fit)

    def test_protocol_violation_is_reported(self):
        d = json.loads((self.dir / "dataset_nl.json").read_text())
        d["repeats"] = 3
        (self.dir / "dataset_nl.json").write_text(json.dumps(d))
        p, out = self.run_script()
        self.assertIn("PROTOCOL", p.stderr)
        self.assertIn("dataset_nl", p.stderr)
        prot = json.loads((out / "scaling.json").read_text())["scaling"]["protocol"]
        self.assertEqual(prot["fixtures_below_min_repeats"], {"dataset_nl": 3})
        self.assertIn("| 3 |", (out / "tables.md").read_text())

    def test_trim_scatter_keeps_the_summary_and_the_svg(self):
        p, out = self.run_script("--trim-scatter")
        self.assertEqual(p.returncode, 0, p.stderr)
        d = json.loads((self.dir / "human_skin.json").read_text())
        self.assertNotIn("scatter", d)
        s = d["scatter_summary"]
        self.assertEqual(s["n_values"], 3)
        self.assertEqual(s["max_abs_log10_difference"], 0.0)
        self.assertTrue((out / "parity_human_skin.svg").exists())
        self.assertTrue((out / "parity_visium.svg").exists())
        # and the trimmed file still round-trips as a results file
        self.assertEqual(analyse.norm_fixture(d)["n_cells"], 7563)

    def test_committed_artifacts_are_trimmed_and_current(self):
        """The committed JSONs must not carry raw vectors, and must record the summary.

        A hand-trimmed artifact is how `parity_human_skin.svg` ended up unreproducible: the
        vectors were deleted after plotting by hand, so re-running `analyse.py` on the committed
        file did nothing and the pipeline's output could not be regenerated from the repository.
        """
        res = HERE / "results"
        if not res.is_dir():
            self.skipTest("no committed results directory")
        # The runners that emit a raw `scatter`. `bench_more_datasets.R` does not, so its files
        # have neither the vectors nor a summary and asserting a summary there would be wrong.
        scatter_emitters = {"human_skin", "mouse_wound", "visium"}
        seen = set()
        for p in res.glob("*.json"):
            if p.name in ("scaling.json",) or p.name.startswith("sweep"):
                continue
            d = json.loads(p.read_text())
            if "rust" not in d or "r_upstream" not in d:
                continue
            self.assertNotIn("scatter", d, f"{p.name} still carries the raw scatter vectors")
            if p.stem in scatter_emitters:
                seen.add(p.stem)
                self.assertIn("scatter_summary", d,
                              f"{p.name} is trimmed but has no summary to regenerate from")
                s = d["scatter_summary"]
                self.assertIn("max_abs_log10_difference", s)
                self.assertTrue((res.parent / s["svg"]).exists()
                                or (HERE.parent / s["svg"]).exists(),
                                f"{p.name} names an SVG that is not in the repository")
        self.assertEqual(seen, scatter_emitters,
                         "a scatter-emitting result file has gone missing")

    def test_committed_tables_are_not_stale(self):
        """`tables.md` must cover every committed results file that has timings."""
        res = HERE / "results"
        out = Path(tempfile.mkdtemp())
        p = subprocess.run([sys.executable, str(HERE / "analyse.py"),
                            "--results", str(res), "--outdir", str(out)],
                           capture_output=True, text=True)
        self.assertEqual(p.returncode, 0, p.stderr)
        regenerated = (out / "tables.md").read_text()
        committed = (res / "tables.md").read_text()
        self.assertEqual(regenerated, committed,
                         "tables.md is stale: re-run analyse.py --trim-scatter")


if __name__ == "__main__":
    unittest.main(verbosity=2)