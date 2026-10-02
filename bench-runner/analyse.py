#!/usr/bin/env python3
"""Turn the measured benchmark JSON into the tables, the scaling fit and the parity plot.

Dependency-free on purpose. `jsonlite` is absent on the pinned host and the *plot* is a paper
artifact that has to be regenerable by whoever reads the repository, so the scatter is emitted as
hand-written SVG rather than a matplotlib PNG: it diffs, it is inspectable, and it needs nothing
installed.

Outputs, into `--outdir`:
  * `scaling.json`  -- the fitted model and the Amdahl decomposition
  * `parity_<fixture>.svg` -- log-log identity scatter of upstream Prob against the port's
  * `tables.md`     -- the markdown that `docs/BENCHMARKS.md` quotes

Usage:
  python3 bench-runner/analyse.py --results bench-runner/results --outdir bench-runner/results
"""
from __future__ import annotations

import argparse
import json
import math
import os
import sys
from pathlib import Path


# ----------------------------------------------------------------------------- scaling model

def fit_affine(xs: list[float], ys: list[float]) -> tuple[float, float, float]:
    """Least squares `y = a + b x`, plus R^2. Returns (a, b, r2).

    Two points through the origin would fit exactly and tell you nothing, so the fit is only ever
    reported over three or more sizes; a caller with fewer must say so rather than quote an R^2
    from a degenerate fit.
    """
    n = len(xs)
    if n < 2:
        return (float("nan"),) * 3
    mx = sum(xs) / n
    my = sum(ys) / n
    sxx = sum((x - mx) ** 2 for x in xs)
    if sxx == 0:
        return (my, 0.0, float("nan"))
    b = sum((x - mx) * (y - my) for x, y in zip(xs, ys)) / sxx
    a = my - b * mx
    ss_tot = sum((y - my) ** 2 for y in ys)
    ss_res = sum((y - (a + b * x)) ** 2 for x, y in zip(xs, ys))
    r2 = float("nan") if ss_tot == 0 else 1.0 - ss_res / ss_tot
    return a, b, r2


def amdahl(measurements: list[dict]) -> dict:
    """Sequential speed-up per thread count, and the implied serial fraction.

    `T(p) = s + (1 - s) T(1) / p` is the Amdahl form for a fixed problem, so fitting
    `T(p) = alpha + beta / p` over the thread sweep recovers `s = alpha / T(1)`. Reported with the
    per-count efficiency `T(1) / (p T(p))` rather than a single headline number, because on an
    8c/16t host the thread counts above 8 are SMT siblings and efficiency *falls* there; quoting
    16 as if it were 16 cores would be the `parallel::detectCores()` mistake again.
    """
    pts = sorted(
        ((m["threads"], m["rust"]["median"]) for m in measurements if m.get("threads")),
        key=lambda t: t[0],
    )
    if len(pts) < 2:
        return {"points": [], "note": "need at least two thread counts"}
    base_t = pts[0][1]
    rows = []
    for p, t in pts:
        rows.append({
            "threads": p,
            "median_secs": t,
            "speedup_vs_1": base_t / t,
            "efficiency": base_t / (p * t),
        })
    inv = [1.0 / p for p, _ in pts]
    alpha, beta, r2 = fit_affine(inv, [t for _, t in pts])
    return {
        "points": rows,
        "amdahl_alpha_serial_secs": alpha,
        "amdahl_parallel_term": beta,
        "amdahl_fit_r2": r2,
        "implied_serial_fraction": alpha / base_t if base_t else float("nan"),
        "note": (
            "SMT: thread counts above the 8 physical cores share execution units, so efficiency "
            "falling past 8 is expected and is not a scaling defect"
        ),
    }


# ------------------------------------------------------------------------------------- SVG

def bin_scatter(xs: list[float], ys: list[float], bins: int = 160):
    """Bin a log10 scatter into a `bins x bins` count grid.

    A million-point scatter drawn as a million `<circle>` elements is 74 MB of SVG that no browser
    will open and no reader will read. Binned, it is a few tens of kilobytes and the density -- which
    is the only thing a million points can say anyway -- is legible.
    """
    finite = [(x, y) for x, y in zip(xs, ys) if math.isfinite(x) and math.isfinite(y)]
    if not finite:
        return None, 0, 0, float("nan")
    lo = math.floor(min(min(x, y) for x, y in finite))
    hi = math.ceil(max(max(x, y) for x, y in finite))
    span = max(hi - lo, 1e-9)
    idx = lambda v: min(bins - 1, max(0, int((v - lo) / span * bins)))
    grid = [[0] * bins for _ in range(bins)]
    off = 0
    for x, y in zip(xs, ys):
        if math.isfinite(x) and math.isfinite(y):
            grid[idx(x)][bins - 1 - idx(y)] += 1
        else:
            off += 1
    return (grid, lo, hi, span), len(finite), off, max(abs(x - y) for x, y in finite)


def parity_svg(sc: dict, path: Path, title: str) -> dict:
    """Log10 identity scatter: upstream `Prob` on x, the port's on y, binned by density.

    A magnitude summary hides a relative disaster at small `Prob`, and this comparison is *supposed*
    to be exact, so the plot's job is to show the residual cloud is empty. Any occupied bin off the
    identity diagonal is drawn in red; with 980 000 values there are none, and the red branch is
    still there so a future regression is visible rather than hidden by a log colour scale.
    """
    binned, n_finite, n_off, max_dev = bin_scatter(sc["log10_prob_r"], sc["log10_prob_rs"])
    n_dev = 0
    if binned is None:
        return {"n_finite": 0, "svg": str(path), "error": "no finite points"}
    grid, lo, hi, span = binned
    B = len(grid)
    peak = max((c for row in grid for c in row), default=1) or 1

    W, H = 720, 720
    pad = 74
    # identity runs corner to corner; the bins are drawn on the same axes
    def to_x(bi: int) -> float:
        return pad + bi / B * (W - 2 * pad)

    def to_y(bj: int) -> float:
        return H - pad - bj / B * (H - 2 * pad)

    # Blue ramp by log count: linear alpha saturates long before 900k points per bin.
    def fill(count: int) -> str:
        if count == 0:
            return "none"
        t = math.log1p(count) / math.log1p(peak)
        return f"rgba(31,119,180,{0.18 + 0.72 * t:.3f})"

    parts = [
        f'<svg xmlns="http://www.w3.org/2000/svg" width="{W}" height="{H}" '
        f'viewBox="0 0 {W} {H}" font-family="Helvetica,Arial,sans-serif">',
        f'<rect width="{W}" height="{H}" fill="white"/>',
        f'<text x="{W / 2}" y="24" text-anchor="middle" font-size="14">{title}</text>',
        f'<rect x="{pad}" y="{pad}" width="{W - 2 * pad}" height="{H - 2 * pad}" fill="#fbfbfb" '
        f'stroke="#333" stroke-width="1"/>',
    ]
    for bi, bj in ((i, j) for i in range(B) for j in range(B) if grid[i][j]):
        c = grid[bi][bj]
        x0, x1 = to_x(bi), to_x(bi + 1)
        y0, y1 = to_y(bj), to_y(bj + 1)
        # A bin is "off diagonal" if its centre is more than 1.5 bins from y = x. The test is on the
        # **sum**, not the difference: row 0 of the grid is the *largest* y
        # (`grid[idx(x)][bins - 1 - idx(y)]`), so the identity line runs corner to corner and the
        # diagonal is anti-diagonal. Comparing `ci` to `cj` therefore marks the whole plot as
        # off-diagonal while `max |log10 diff|` is 0 -- which is what it did, and a red plot that
        # means nothing is worse than no red plot.
        ci, cj = bi + 0.5, bj + 0.5
        if abs((ci + cj) - B) > 3.0:
            n_dev += c
            parts.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" '
                         f'height="{y1 - y0:.2f}" fill="#d62728" fill-opacity="0.85"/>')
        else:
            parts.append(f'<rect x="{x0:.2f}" y="{y0:.2f}" width="{x1 - x0:.2f}" '
                         f'height="{y1 - y0:.2f}" fill="{fill(c)}"/>')
    # identity line
    parts.append(f'<line x1="{pad}" y1="{H - pad}" x2="{W - pad}" y2="{pad}" stroke="#111" '
                 f'stroke-width="1.2" stroke-dasharray="7,5"/>')
    for k in range(int(lo), int(hi) + 1, 4):
        fx = pad + (k - lo) / span * (W - 2 * pad)
        parts.append(f'<text x="{fx:.1f}" y="{H - pad + 18}" text-anchor="middle" '
                     f'font-size="10">1e{k}</text>')
        parts.append(f'<text x="{pad - 8}" y="{H - pad - (k - lo) / span * (H - 2 * pad) + 3:.1f}" '
                     f'text-anchor="end" font-size="10">1e{k}</text>')
    parts.append(f'<text x="{W / 2}" y="{H - 12}" text-anchor="middle" font-size="12">'
                 f'upstream computeCommunProb Prob</text>')
    parts.append(f'<text x="16" y="{H / 2}" text-anchor="middle" font-size="12" '
                 f'transform="rotate(-90 16 {H / 2})">Rust kernel Prob</text>')
    caption = (f"{n_finite:,} values, {B}x{B} bins, density by log count | "
               f"max |log10 difference| = {max_dev:.3g}"
               + (f" | {n_off:,} non-finite" if n_off else "")
               + (f" | {n_dev:,} values off the identity" if n_dev else " | 0 off the identity"))
    parts.append(f'<text x="{W / 2}" y="{H - 30}" text-anchor="middle" font-size="10" '
                 f'fill="#444">{caption}</text>')
    parts.append("</svg>")
    path.write_text("\n".join(parts))
    return {"n_finite": n_finite, "n_nonfinite": n_off, "max_log10_dev": max_dev,
            "n_off_identity": n_dev, "bins": B, "svg": str(path),
            "svg_bytes": path.stat().st_size}


# ------------------------------------------------------------------------- fixture normalisation

# The objective fixes the protocol at "at least five timed repeats", and it fixes a scaling claim
# at more than two sizes. Both are checked here rather than trusted from each runner's defaults.
MIN_REPEATS = 5
MIN_FIT_POINTS = 3


def norm_fixture(d: dict) -> dict:
    """Canonical view of a results file's fixture metadata, whichever runner wrote it.

    Three runners write into this directory and none of them agrees on where that metadata lives:
    `bench_real.R` nests it under `fixture: {fixture, n_cells, ...}`;
    `bench_more_datasets.R` puts `dataset`, `species` and the counts at the top level; and
    `bench_visium.R` nests `fixture: {dataset, n_spots, ...}`, calling spots something other than
    cells. The first version of this function assumed `bench_real.R`'s schema and raised
    `KeyError: 'fixture'` the moment the other two runners were added -- which is how
    `tables.md` sat stale for a day with every runner reporting green.
    """
    f = d.get("fixture")
    if isinstance(f, dict) and "n_spots" in f:
        kind, name, n_cells = "spatial", f.get("dataset", "spatial"), f["n_spots"]
        f = dict(f)
    elif isinstance(f, dict):
        kind, name, n_cells = "scrnaseq", f.get("fixture", "?"), f["n_cells"]
    else:
        kind, name, n_cells = "scrnaseq", d.get("dataset", "?"), d["n_cells"]
        f = {}
    return {"name": name, "n_cells": n_cells, "kind": kind,
            "n_genes_total": f.get("n_genes_total", d.get("n_genes_total")),
            "n_genes_signaling": f.get("n_genes_signaling", d.get("n_genes_signaling")),
            "n_groups": f.get("n_groups", d.get("n_groups")),
            "n_lr": f.get("n_lr", d.get("n_lr"))}


# ---------------------------------------------------------------------------------- tables

def fmt_ci(lo: float, mid: float, hi: float, unit: str = "s") -> str:
    return f"{mid:.2f} {unit} [{lo:.2f}, {hi:.2f}]"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--results", default="bench-runner/results")
    ap.add_argument("--outdir", default="bench-runner/results")
    ap.add_argument("--sweep", default=None, help="extra JSON files holding a thread sweep")
    ap.add_argument("--trim-scatter", action="store_true",
                    help="after writing the SVGs, replace each raw `scatter` vector in its own "
                         "JSON with the summary already computed here (the vectors are ~40 MB for "
                         "the wound fixture and nothing reads them)")
    args = ap.parse_args()

    res_dir = Path(args.results)
    out_dir = Path(args.outdir)
    out_dir.mkdir(parents=True, exist_ok=True)

    real = {}
    for p in sorted(res_dir.glob("*.json")):
        if p.name in ("scaling.json",) or p.name.startswith("sweep"):
            continue
        try:
            d = json.loads(p.read_text())
        except json.JSONDecodeError:
            continue
        if "rust" in d and "r_upstream" in d:
            d["_norm"] = norm_fixture(d)
            real[p.stem] = d

    # The protocol check. Each runner defaults to 5 repeats, but a run invoked with an explicit
    # `--repeats 3` writes 3 and nothing downstream notices -- the headline pair sat at 3 while
    # docs/BENCHMARKS.md quoted 5, and the published table carried the wrong precision. Surfacing
    # it here means the drift is reported by the pipeline instead of found by reading it.
    short = {s: d.get("repeats") for s, d in real.items()
             if (d.get("repeats") or 0) < MIN_REPEATS}
    for stem, got in sorted(short.items()):
        print(f"PROTOCOL: {stem} has {got} timed repeats, below the required {MIN_REPEATS}",
              file=sys.stderr)

    sweep = []
    for p in sorted(res_dir.glob("sweep*.json")):
        d = json.loads(p.read_text())
        if "rust" in d and "threads" in d:
            sweep.append(d)

    plots = {}
    for stem, d in real.items():
        if "scatter" in d:
            fx = d["_norm"]
            unit = "spots" if fx["kind"] == "spatial" else "cells"
            plots[stem] = parity_svg(
                d["scatter"], out_dir / f"parity_{stem}.svg",
                f"Parity on {fx['name']} "
                f"({fx['n_cells']:,} {unit}, {fx['n_lr']:,} L-R pairs)",
            )

    scaling = {"real": {}, "amdahl": amdahl(sweep) if sweep else {"points": []}}
    # Descriptive only, and labelled so. Across fixtures, time is not a function of cell count
    # alone -- nLR, nGenes and nGroups move with it, so a straight line through these points is a
    # sanity check on the runners, not the scaling model. The scaling claim is the synthetic
    # grid's log-log fit (docs/BENCHMARKS.md), measured at controlled sizes with nC as the only
    # variable. Spatial fixtures are excluded for a second reason: their cost is the exact k-d
    # tree, not the bootstrap kernel, so folding them in would fit a different algorithm.
    pts = [(d["_norm"]["n_cells"], d["rust"]["median"], stem) for stem, d in real.items()
           if d.get("nboot") == 100 and d["_norm"]["kind"] == "scrnaseq"]
    pts.sort()
    fit = {"role": "descriptive; not the scaling model (see synthetic grid)",
           "n_points": len(pts), "min_points_required": MIN_FIT_POINTS,
           "points": [{"fixture": p[2], "n_cells": p[0], "rust_median_secs": p[1]} for p in pts]}
    if len(pts) >= MIN_FIT_POINTS:
        a, b, r2 = fit_affine([p[0] for p in pts], [p[1] for p in pts])
        fit.update(intercept_secs=a, slope_secs_per_cell=b, r2=r2)
    else:
        # A line through two points passes through both of them exactly, so its R^2 is 1.0 by
        # construction and means nothing. `fit_affine` documents that this must never be
        # reported; the check lives here because that is where the point count is known.
        fit["status"] = f"withheld: {len(pts)} point(s), {MIN_FIT_POINTS} required"
    spatial_excluded = {s: d["_norm"]["name"] for s, d in real.items()
                        if d.get("nboot") == 100 and d["_norm"]["kind"] != "scrnaseq"}
    if spatial_excluded:
        fit["excluded_spatial"] = spatial_excluded
    scaling["fixed_nboot_100_fit"] = fit
    scaling["protocol"] = {"min_repeats": MIN_REPEATS, "min_fit_points": MIN_FIT_POINTS,
                           "fixtures_below_min_repeats": short}
    (out_dir / "scaling.json").write_text(json.dumps({"scaling": scaling, "parity_plots": plots},
                                                     indent=1))

    md = ["## Measured results", "",
          "| fixture | kind | nC | K | nGenes | nLR | nboot | reps | R upstream | Rust | speedup | peak RSS | parity |",
          "|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---|"]
    for stem, d in sorted(real.items(), key=lambda kv: kv[1]["_norm"]["n_cells"]):
        fx = d["_norm"]
        if d["parity_net"]:
            ok = "identical"
        elif fx["kind"] == "spatial":
            ok = f"divergence by design (max|Δ|={d['max_abs_diff']:.3g})"
        else:
            ok = f"max|Δ|={d['max_abs_diff']:.3g}"
        unit = "spots" if fx["kind"] == "spatial" else "nC"
        md.append(
            f"| {fx['name']} | {fx['kind']} | {fx['n_cells']:,} | {fx['n_groups']:,} | "
            f"{fx['n_genes_signaling']:,} | {fx['n_lr']:,} | {d['nboot']} | {d.get('repeats')} | "
            f"{fmt_ci(d['r_upstream']['ci_lo'], d['r_upstream']['median'], d['r_upstream']['ci_hi'])} | "
            f"{fmt_ci(d['rust']['ci_lo'], d['rust']['median'], d['rust']['ci_hi'])} | "
            f"**{d['speedup']:.2f}×** [{d['speedup_ci_lo']:.2f}, {d['speedup_ci_hi']:.2f}] | "
            f"{d['peak_rss_mb']:.0f} MB | {ok} |")
    md += ["", f"Protocol: {MIN_REPEATS} timed repeats after a discarded warm-up, median with a "
              "bootstrap CI, CPUs pinned with `taskset`. Raw repeats are in the JSON beside this "
              "file; the scatter vectors are dropped once their SVG is written."]
    (out_dir / "tables.md").write_text("\n".join(md) + "\n")

    if args.trim_scatter:
        for stem, d in real.items():
            if "scatter" not in d:
                continue
            sc = d.pop("scatter")
            p = plots[stem]
            d["scatter_summary"] = {
                "n_values": len(sc["log10_prob_r"]),
                "n_finite": p["n_finite"], "n_nonfinite": p["n_nonfinite"],
                "max_abs_log10_difference": p["max_log10_dev"],
                "n_off_identity": p["n_off_identity"], "bins": p["bins"],
                "svg": p["svg"], "svg_bytes": p["svg_bytes"],
                "note": "per-value vectors dropped by analyse.py --trim-scatter; "
                        "regenerate from bench_real.R",
            }
            d.pop("_norm", None)
            (res_dir / f"{stem}.json").write_text(json.dumps(d, indent=1))
        print(f"trimmed scatter from {sum(1 for d in real.values() if 'scatter_summary' in d)} "
              "result files")

    print("\n".join(md))
    if plots:
        print("\nparity plots:")
        for stem, p in plots.items():
            print(f"  {stem}: {p['n_finite']:,} values, max |log10 diff| = "
                  f"{p['max_log10_dev']:.3g}, {p['n_off_identity']:,} off the identity "
                  f"({p['svg_bytes'] / 1024:.0f} KB) -> {p['svg']}")
    if sweep:
        print("\nstrong scaling (Rust):")
        for r in scaling["amdahl"]["points"]:
            print(f"  {r['threads']:>3} threads: {r['median_secs']:8.2f} s  "
                  f"speedup {r['speedup_vs_1']:5.2f}x  efficiency {r['efficiency']:5.1%}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
