#!/usr/bin/env python
"""figures/summary_figure.png — a three-panel composite for slides (300 dpi): the transfer ladder, the
ten-gene stable core, and the identifiability heatmap. Reads site/data/*.json (the provenance-checked
exports), never results/ directly, and uses the site's design tokens."""
from __future__ import annotations

import json
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


def main():
    H, SC, N = load("headline.json"), load("stable_core.json"), load("nesting.json")
    fig, axes = plt.subplots(1, 3, figsize=(15, 5.2), facecolor=T["page"], gridspec_kw={"width_ratios": [1.25, 1, 1.1], "wspace": 0.55})
    for ax in axes:
        ax.set_facecolor(T["surface"])
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
        ax.grid(axis="x", color=T["grid"], linewidth=0.8)
        ax.set_axisbelow(True)

    # (a) transfer ladder, k20 marginal pooled
    rungs = [("in_distribution", "in-distribution\n(held-out animals)"), ("train_male_test_female", "held-out sex\n(M → F)"),
             ("different_lab", "other laboratory\n(rat BodyMap)"), ("different_species", "other species\n(human GTEx)")]
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
    ax.set_xticks(x, [lab for _, lab in rungs], fontsize=8)
    ax.set_ylim(0, 1.1)
    ax.set_ylabel("fraction")
    ax.grid(axis="x", visible=False)
    ax.grid(axis="y", color=T["grid"], linewidth=0.8)
    ax.legend(frameon=False, loc="upper right", fontsize=8, ncol=1, bbox_to_anchor=(1.0, 1.02))
    ax.set_title("a  Accuracy survives every shift; the guarantee does not\n    (20-gene panel, α = 0.10, source calibration)", loc="left", fontsize=9.5)

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
    ax.set_title("b  The ten-gene stable core (≥ 80 % of 50\n    animal bootstraps): single markers, mostly textbook", loc="left", fontsize=9.5)

    # (c) identifiability heatmap: Cramér's V and fraction of pairs sharing a level, RNA-seq variables + one row per other layer
    keep = [("TRNSCRPT", "RNA_extr_plate_ID", "RNA-seq · extraction plate"), ("TRNSCRPT", "Lib_batch_ID", "RNA-seq · library batch"), ("TRNSCRPT", "Seq_flowcell_ID", "RNA-seq · flowcell"),
            ("TRNSCRPT", "Seq_flowcell_lane", "RNA-seq · lane"), ("METHYL", "DNA_extr_plate_ID", "RRBS · extraction plate"), ("METHYL", "Seq_flowcell_ID", "RRBS · flowcell"),
            ("ATAC", "Seq_flowcell_ID", "ATAC · flowcell"), ("PROT", "plex_id", "proteomics · TMT plex"), ("PROT", "tmt11_channel", "proteomics · TMT channel"), ("IMMUNO", "plate_id", "immunoassay · plate")]
    rows_c = []
    for assay, var, label in keep:
        r = next((x for x in N["nesting"].get(assay, []) if x["variable"] == var), None)
        if r:
            rows_c.append((label, r["cramers_v"], r["n_pairs_sharing_level"] / r["n_pairs_total"], f"{r['n_pairs_sharing_level']}/{r['n_pairs_total']}"))
    ax = axes[2]
    z = np.array([[v, f] for _, v, f, _ in rows_c])
    cmap = LinearSegmentedColormap.from_list("seq", [T["surface"]] + T["seq"])
    ax.imshow(z, cmap=cmap, vmin=0, vmax=1, aspect="auto")
    ax.set_xticks([0, 1], ["Cramér's V\nwith tissue", "tissue pairs\nsharing a level"], fontsize=8)
    ax.set_yticks(range(len(rows_c)), [r[0] for r in rows_c], fontsize=8)
    for i, (_, v, f, txt) in enumerate(rows_c):
        ax.text(0, i, f"{v:.2f}", ha="center", va="center", fontsize=8, color="white" if v > 0.6 else T["ink"])
        ax.text(1, i, txt, ha="center", va="center", fontsize=8, color="white" if f > 0.6 else T["ink"])
    ax.grid(False)
    for s in ("top", "right", "left", "bottom"):
        ax.spines[s].set_visible(False)
    est = next(r for r in N["estimable_pairs"] if r["assay"] == "TRNSCRPT")
    ax.set_title(f"c  Batch is nested in tissue: {est['n_pairs_estimable']} of {est['n_pairs_total']} RNA-seq tissue pairs\n    share a plate, library batch and flowcell", loc="left", fontsize=9.5)
    fig.text(0.01, 0.01, "Sources: results/05_panels, 06_conformal, 08_shift, 12_bodymap, 13_gtex, 16_identifiability via site/data/*.json (provenance in site/data/provenance.json)", fontsize=7, color=T["muted"])
    OUT.parent.mkdir(exist_ok=True)
    fig.savefig(OUT, dpi=300, bbox_inches="tight", facecolor=T["page"])
    print(f"wrote {OUT} ({OUT.stat().st_size / 1e6:.2f} MB)")


if __name__ == "__main__":
    main()
