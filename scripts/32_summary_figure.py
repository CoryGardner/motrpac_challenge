#!/usr/bin/env python
"""figures/summary_figure.png — a three-panel composite for slides (300 dpi): the transfer ladder, the
ten-gene stable core, and the identifiability heatmap. Reads site/data/*.json (the provenance-checked
exports), never results/ directly, and uses the site's design tokens."""
from __future__ import annotations

import json
import textwrap
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np

ROOT = Path(__file__).resolve().parents[1]
D = ROOT / "site" / "data"
OUT = ROOT / "figures" / "summary_figure.png"
T = {"page": "#f9f9f7", "surface": "#fcfcfb", "ink": "#0b0b0b", "ink2": "#52514e", "muted": "#898781", "grid": "#e1e0d9", "axis": "#c3c2b7",
     "c1": "#2a78d6", "c2": "#eb6834", "seq": ["#cde2fb", "#9ec5f4", "#6da7ec", "#3987e5", "#256abf", "#184f95", "#0d366b"]}
plt.rcParams.update({"font.family": ["DejaVu Sans", "sans-serif"], "font.size": 9, "axes.edgecolor": T["axis"], "axes.labelcolor": T["ink2"],
                     "xtick.color": T["ink2"], "ytick.color": T["ink2"], "axes.titlecolor": T["ink"], "axes.titleweight": "semibold", "axes.titlesize": 10})


def load(name):
    return json.loads((D / name).read_text())


SUB_LINES = 3   # every panel's subtitle is padded to this many lines so the three bold titles sit on one row


def title(ax, head, sub):
    """A one-line panel title (≤ 60 characters) above a smaller wrapped subtitle; the title pad is the same on every panel."""
    assert len(head) <= 60, (len(head), head)
    lines = textwrap.wrap(sub, 64)
    assert len(lines) <= SUB_LINES, (len(lines), sub)
    lines += [""] * (SUB_LINES - len(lines))
    ax.set_title(head, loc="left", fontsize=10, fontweight="bold", pad=10 + 11.5 * SUB_LINES)
    ax.text(0, 1.015, "\n".join(lines), transform=ax.transAxes, fontsize=8, color=T["ink2"], va="bottom", ha="left")


def main():
    H, SC, N, G = load("headline.json"), load("stable_core.json"), load("nesting.json"), load("genes.json")
    fig, axes = plt.subplots(1, 3, figsize=(16, 6.2), facecolor=T["page"], gridspec_kw={"width_ratios": [1.45, 1, 1.15], "wspace": 0.62})
    for ax in axes:
        ax.set_facecolor(T["surface"])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="x", color=T["grid"], linewidth=0.8)
        ax.set_axisbelow(True)

    # (a) transfer ladder, k20 marginal pooled
    rungs = [("in_distribution", "held-out\nanimals"), ("train_control_test_trained", "trained\n(fit on\ncontrols)"), ("train_male_test_female", "held-out\nsex\n(M → F)"),
             ("different_lab", "other lab\n(rat\nBodyMap)"), ("different_species", "other\nspecies\n(GTEx)")]
    rows = [next(r for r in H["ladder"] if r["rung_id"] == rid and r["model"] == "k20" and r["variant"] == "marginal" and r.get("calibration", "pooled") == "pooled") for rid, _ in rungs]
    ax = axes[0]
    x = np.arange(len(rungs))
    w = 0.36
    acc = [r["accuracy"] for r in rows]
    cov = [r["coverage"] for r in rows]
    ax.bar(x - w / 2, acc, w, color=T["c1"], label="accuracy", edgecolor=T["surface"], linewidth=1.5)
    ax.bar(x + w / 2, cov, w, color=T["c2"], label="coverage of the 90 % set", edgecolor=T["surface"], linewidth=1.5)
    for i, (a, c) in enumerate(zip(acc, cov)):
        ax.text(i - w / 2, a + 0.015, f"{a:.3f}", ha="center", va="bottom", fontsize=8, color=T["ink2"])
        ax.text(i + w / 2, c + 0.015, f"{c:.3f}", ha="center", va="bottom", fontsize=8, color=T["ink2"])
    ax.axhline(0.9, color=T["ink2"], linewidth=1, linestyle="--")
    ax.text(len(rungs) - 0.55, 0.905, "1 − α = 0.90", ha="right", va="bottom", fontsize=8, color=T["ink2"])
    ax.set_xticks(x, [lab for _, lab in rungs], fontsize=7)
    ax.set_ylim(0, 1.28)
    ax.set_ylabel("fraction")
    ax.grid(axis="x", visible=False)
    ax.grid(axis="y", color=T["grid"], linewidth=0.8)
    ax.legend(frameon=False, loc="upper left", fontsize=8, ncol=2)
    title(ax, "a  Accuracy holds at every rung of the shift ladder", "The 90 % guarantee holds within the study and is restored by three animals from a new laboratory (20-gene panel, α = 0.10, source-calibrated)")

    # (b) stable core: effect size per gene
    core = sorted(SC["core"], key=lambda r: r["effect_size"])
    ax = axes[1]
    y = np.arange(len(core))
    ax.barh(y, [r["effect_size"] for r in core], color=T["c1"], height=0.62, edgecolor=T["surface"], linewidth=1.5)
    ax.set_yticks(y, [f"{r['gene_symbol']} · {r['marker_tissue']}" for r in core], fontsize=8)
    for i, r in enumerate(core):
        ax.text(r["effect_size"] + 0.15, i, f"{r['effect_size']:.1f}  (freq {r['selection_frequency']:.2f})", va="center", fontsize=7.5, color=T["ink2"])
    ax.set_xlim(0, max(r["effect_size"] for r in core) * 1.55)
    ax.set_xlabel("log2 CPM above the next-highest tissue")
    title(ax, "b  The ten-gene stable core", "Chosen in ≥ 80 % of 50 animal bootstraps: single markers, mostly textbook; bar = log2 CPM above the next-highest tissue")

    # (c) batch measured directly on the bridging reference pools: Σ V_batch / Σ V_tissue per pool, all genes
    ax = axes[2]
    pools = sorted((N.get("bridge") or {}).get("pools") or [], key=lambda r: r["sum_ratio_batch_over_tissue"])
    y = np.arange(len(pools))
    x = [100 * r["sum_ratio_batch_over_tissue"] for r in pools]
    both = [";" in r["sites"] for r in pools]
    ax.hlines(y, 0, x, color=T["grid"], linewidth=2)
    ax.scatter(x, y, s=70, color=[T["c1"] if b else T["c2"] for b in both], edgecolor=T["surface"], linewidth=1.5, zorder=3)
    for i, r in enumerate(pools):
        ax.text(x[i] + 0.25, i, f"{x[i]:.1f} %  ·  {r['n_plates']} plates{', both sites' if both[i] else ''}", va="center", fontsize=7.5, color=T["ink2"])
    ax.set_yticks(y, [f"{r['pool_tissue'].replace(' Powder', '').lower()} pool {r['pool_type']}" for r in pools], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlim(0, 10)
    ax.set_xlabel("batch variance as % of the tissue-separating variance (all genes; tissue signal = 100 %)")
    ax.scatter([], [], color=T["c1"], label="both sequencing sites")
    ax.scatter([], [], color=T["c2"], label="one site")
    ax.legend(frameon=False, loc="center right", fontsize=8)
    title(ax, "c  Batch is a few percent of the tissue signal",
          "Between-plate variance of a reference RNA pool over the variance of the 19 tissue means, all genes; the two gastrocnemius pools ran on both sequencing sites")
    fig.text(0.01, 0.01, "Sources: results/05_panels, 06_conformal, 08_shift, 12_bodymap, 13_gtex, 16_identifiability via site/data/*.json (provenance in site/data/provenance.json)", fontsize=7, color=T["muted"])
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor=T["page"])
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
