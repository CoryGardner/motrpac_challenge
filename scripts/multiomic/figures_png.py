#!/usr/bin/env python
"""Render the multiomic page's figures as PNG for slides, from site/data/multiomic.json (the same numbers the page reads).
Output: results_multiomic/figures/*.png. Colours follow the site's light theme tokens (site/assets/theme.css)."""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from tfp import config as C

ROOT = C.ROOT
OUT = ROOT / "results_multiomic" / "figures"
PAL = ["#2a78d6", "#eb6834", "#1baf7a", "#eda100", "#e87ba4", "#008300", "#4a3aa7", "#e34948"]
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e1e0d9"
plt.rcParams.update({"font.family": "sans-serif", "font.size": 10, "axes.edgecolor": GRID, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2, "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.8,
                     "axes.spines.top": False, "axes.spines.right": False, "axes.spines.left": False, "axes.spines.bottom": False, "figure.dpi": 150})


def save(fig, name):
    fig.tight_layout()
    fig.savefig(OUT / name, dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", OUT / name)


TISSUE_NAMES = {"CORTEX": "cerebral cortex", "HEART": "heart", "KIDNEY": "kidney", "LIVER": "liver", "LUNG": "lung", "SKM-GN": "gastrocnemius", "WAT-SC": "white adipose"}


def fig0(M):
    """The proteomics problem in one picture: per-vial PC1 vs PC2 on both scales, tissue colours (the page's order), 1600 px wide."""
    P = M["pca_two_scales"]
    tissues = sorted(set(P["rii"]["tissue"]) | set(P["ratio"]["tissue"]))
    col = {t: PAL[k % len(PAL)] for k, t in enumerate(tissues)}
    fig, axes = plt.subplots(1, 2, figsize=(8, 4.1))
    for ax, tag, title in ((axes[0], "ratio", "As distributed: each sample ÷ a reference\npool of the same tissue"),
                           (axes[1], "rii", "Rebuilt from the reporter-ion intensities:\neach sample ÷ its own total signal")):
        d = P[tag]
        for t in tissues:
            ii = [i for i, x in enumerate(d["tissue"]) if x == t]
            ax.scatter([d["x"][i] for i in ii], [d["y"][i] for i in ii], s=9, color=col[t], label=TISSUE_NAMES.get(t, t), linewidths=0)
        ax.set_xlabel(f"PC1 ({100 * d['explained'][0]:.0f} % of variance)"); ax.set_ylabel(f"PC2 ({100 * d['explained'][1]:.0f} %)")
        ax.set_title(title, loc="left", color=INK, fontsize=9.5)
        r2 = d["r2_tissue_pc1"]
        ax.text(0, -0.25, f"tissue explains R² {r2:.4f} of PC1" if r2 < 0.01 else f"tissue explains R² {r2:.3f} of PC1", transform=ax.transAxes, fontsize=8.5, color=INK2)
    axes[1].legend(frameon=False, fontsize=7, loc="upper left", bbox_to_anchor=(1.0, 1.0), markerscale=1.6)
    fig.text(0.01, 0.01, "Source: results_multiomic/01_rii/pca_scores_ratio.csv, pca_scores_rii.csv (scripts/multiomic/01b_pca_scores.py); "
             f"{P['ratio']['n_vials']} and {P['rii']['n_vials']} vials", fontsize=6.5, color=INK2)
    fig.tight_layout(rect=(0, 0.05, 1, 1))
    fig.savefig(OUT / "fig0_pca_two_scales.png", dpi=200, facecolor="white")
    plt.close(fig)
    print("wrote", OUT / "fig0_pca_two_scales.png")


def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", default=None, help="render one figure (fig0)")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    M = json.loads((ROOT / "site" / "data" / "multiomic.json").read_text())
    fig0(M)
    if args.only == "fig0":
        return
    # 1 same proteins two scales
    S = M["scales"]["rows"]
    fig, ax = plt.subplots(figsize=(7, 3.6))
    x = np.arange(len(S)); w = 0.25
    ax.bar(x - w, [r["r2_ratio"] for r in S], w, color=PAL[7], label="distributed ratios")
    ax.bar(x, [r["r2_rii"] for r in S], w, color=PAL[0], label="reporter-ion log2 ppm")
    ax.bar(x + w, [r["r2_complete"] for r in S], w, color=PAL[2], label="reporter-ion, complete proteins")
    ax.scatter(x, [r["r2_rii_null95"] for r in S], marker="_", s=300, color=INK2, label="permutation null 95th pct", zorder=3)
    ax.set_xticks(x, [r["PC"] for r in S]); ax.set_ylim(0, 1.08); ax.set_ylabel("tissue R² of the component"); ax.legend(frameon=False, fontsize=8, ncol=2, loc="upper right")
    ax.set_title("Same proteins, two scales: tissue R² of PC1–3", loc="left", color=INK, fontsize=11)
    save(fig, "fig1_same_proteins_two_scales.png")
    # 2 ladder protein vs RNA (k20)
    W, X = M["ladder_within"], M["ladder_species"]
    rungs = [("in_distribution", "held-out\nanimals"), ("train_control_test_trained", "trained animals\n(fit on controls)"), ("train_male_test_female", "held-out sex\nM→F"), ("train_female_test_male", "held-out sex\nF→M"), ("species", "other species\nJiang / GTEx"), ("same42", "same 42 human\nsamples, 7-class")]
    def row(rid, layer):
        if rid in ("species", "same42"):
            if rid == "species":
                return next((r for r in X if r["model"] == "k20" and (r["target"].startswith("protein → Jiang 2020 (cleaned") if layer == "protein" else r["target"].startswith("RNA → GTEx"))), None)
            return next((r for r in X if r["model"] == "k20" and r["layer"] == layer and r["target"].startswith(f"{layer}, 7-class")), None)
        return next((r for r in W if r["rung_id"] == rid and r["layer"] == layer and r["model"] == "k20"), None)
    fig, ax = plt.subplots(figsize=(9, 4))
    x = np.arange(len(rungs)); w = 0.2
    for i, (layer, key, col) in enumerate([("protein", "accuracy", PAL[0]), ("protein", "coverage", PAL[2]), ("RNA", "accuracy", PAL[1]), ("RNA", "coverage", PAL[3])]):
        vals = [(row(r, layer) or {}).get(key) for r, _ in rungs]
        ax.bar(x + (i - 1.5) * w, [v if v is not None else 0 for v in vals], w, color=col, label=f"{layer} {key}")
    ax.axhline(0.9, color=INK2, lw=1, ls="--"); ax.text(len(rungs) - 0.5, 0.905, "0.90", color=INK2, fontsize=8, ha="right", va="bottom")
    ax.set_xticks(x, [l for _, l in rungs], fontsize=8); ax.set_ylim(0, 1.1); ax.set_ylabel("fraction"); ax.legend(frameon=False, fontsize=8, ncol=4, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    ax.set_title("The protein ladder beside the RNA ladder (20-feature panels, α = 0.10, source calibration)", loc="left", color=INK, fontsize=11, pad=22)
    save(fig, "fig2_protein_ladder_beside_rna.png")
    # 3 recalibration coverage and set size
    RC, R0 = M["recalibration"], M["recalibration_n0"]
    SRC = [("Phase 3 protein → Jiang 2020 (cleaned relative)", "protein → Jiang (cleaned), 7 cls", PAL[0]), ("Phase 3 protein → Jiang 2020 (raw ppm)", "protein → Jiang (raw ppm), 7 cls", PAL[2]), ("frozen phase 13 RNA → GTEx", "RNA → GTEx, 19 cls", PAL[1]),
           ("Phase 4 deep_mw", "metabolites → aging atlas, 9 cls", PAL[3]), ("Phase 4 deep_sato", "metabolites → Sato (deep), 9 cls", PAL[4]), ("Phase 4 hilic_sato", "metabolites → Sato (HILIC+), 19 cls", PAL[6])]
    fig, axes = plt.subplots(1, 2, figsize=(10, 3.8))
    for src, name, col in SRC:
        n0 = next((r for r in R0 if r["source"] == src and r["model"] == "k20"), None)
        pts = [(0, n0["coverage"], n0["set_size"])] if n0 else []
        for n in (3, 5):
            r = next((r for r in RC if r["source"].startswith(src) and r["model"] == "k20" and r["n_recal"] == n), None)
            if r:
                pts.append((n, r["coverage_recalibrated"], r["set_size_recalibrated"]))
        axes[0].plot([p[0] for p in pts], [p[1] for p in pts], marker="o", color=col, label=name, lw=2)
        axes[1].plot([p[0] for p in pts], [p[2] for p in pts], marker="o", color=col, label=name, lw=2)
    axes[0].axhline(0.9, color=INK2, lw=1, ls="--"); axes[0].set_ylim(0, 1.05); axes[0].set_ylabel("coverage"); axes[0].set_xlabel("target individuals used to recalibrate"); axes[0].set_xticks([0, 3, 5])
    axes[1].set_ylabel("mean set size (classes per set)"); axes[1].set_xlabel("target individuals used to recalibrate"); axes[1].set_xticks([0, 3, 5]); axes[1].set_ylim(bottom=0)
    axes[0].set_title("Observed coverage after recalibration…", loc="left", color=INK, fontsize=11); axes[1].set_title("…at these set sizes", loc="left", color=INK, fontsize=11)
    axes[1].legend(frameon=False, fontsize=7, loc="upper right")
    save(fig, "fig3_recalibration_coverage_and_set_size.png")
    # 4 fusion
    FU = [r for r in M["fusion"] if r["model"] == "k20"]
    order = ["RNA", "protein", "late_mean", "stacked_LR"]; lab = {"RNA": "RNA alone", "protein": "protein alone", "late_mean": "late fusion", "stacked_LR": "stacked fusion"}
    rows = [next(r for r in FU if r["layer"] == l) for l in order]
    fig, ax = plt.subplots(figsize=(7, 3.6)); x = np.arange(4); w = 0.25
    ax.bar(x - w, [r["accuracy"] for r in rows], w, color=PAL[0], label="accuracy", yerr=[[r["accuracy"] - r["acc_ci95_low_donor_boot"] for r in rows], [r["acc_ci95_high_donor_boot"] - r["accuracy"] for r in rows]], capsize=2, ecolor=INK2)
    ax.bar(x, [r["coverage_motrpac_cal"] for r in rows], w, color=PAL[1], label="coverage (MoTrPAC calibration)")
    ax.bar(x + w, [r["frac_empty"] for r in rows], w, color=PAL[7], label="empty sets")
    ax.axhline(0.9, color=INK2, lw=1, ls="--"); ax.set_xticks(x, [lab[l] for l in order]); ax.set_ylim(0, 1.1); ax.legend(frameon=False, fontsize=8, ncol=3, loc="upper center", bbox_to_anchor=(0.5, 1.02))
    ax.set_title(f"Fusion under the species shift: same {rows[0]['n_mapped']} Jiang samples, {rows[0]['n_donors']} donors (k20)", loc="left", color=INK, fontsize=11, pad=22)
    save(fig, "fig4_fusion_under_shift.png")
    # 5 design comparison
    D = sorted(M["design"], key=lambda r: (r["cramers_v"] if r["cramers_v"] is not None else 0))
    labels = [f"{r['dataset']} · {r['batch_variable']}" for r in D]
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.6), sharey=True)
    axes[0].barh(labels, [r["cramers_v"] or 0 for r in D], color=PAL[0]); axes[0].set_xlim(0, 1.05); axes[0].set_xlabel("Cramér's V (tissue × batch)"); axes[0].set_title("Cramér's V: 1 = one tissue per batch", loc="left", color=INK, fontsize=11)
    axes[1].barh(labels, [r["frac_pairs_estimable"] or 0 for r in D], color=PAL[2]); axes[1].set_xlim(0, 1.05); axes[1].set_xlabel("estimable tissue pairs / all pairs"); axes[1].set_title("Tissue pairs sharing a batch", loc="left", color=INK, fontsize=11)
    for i, r in enumerate(D):
        axes[1].text((r["frac_pairs_estimable"] or 0) + 0.02, i, f"{r['n_pairs_estimable']}/{r['n_pairs_total']}", va="center", fontsize=7, color=INK2)
    axes[0].tick_params(axis="y", labelsize=8)
    save(fig, "fig5_design_comparison.png")
    # 6 RNA markers at the protein level
    RP = M["rna_protein"]
    fig, ax = plt.subplots(figsize=(7, 3.6))
    edges = RP["hist_edges"]; centers = [(edges[i] + edges[i + 1]) / 2 for i in range(len(edges) - 1)]
    ax.bar(centers, RP["hist_counts"], width=0.095, color=PAL[0])
    for v, t, ls in ((RP["spearman_median"], f"median {RP['spearman_median']:.2f}", "-"), (RP["null_median"], f"null median {RP['null_median']:.2f}", ":"), (RP["null_q95"], f"null 95 % {RP['null_q95']:.2f}", "--")):
        ax.axvline(v, color=INK2, ls=ls, lw=1.5); ax.text(v, ax.get_ylim()[1] * 0.98, t, rotation=90, va="top", ha="right", fontsize=7, color=INK2)
    ax.set_xlabel("Spearman correlation of RNA and protein across the 7 tissues"); ax.set_ylabel("genes"); ax.set_xlim(-1, 1)
    ax.set_title(f"RNA vs protein across tissues, {RP['n_genes']} genes; {RP['n_same_marker']} of {RP['n_testable']} RNA panel markers keep their tissue as proteins", loc="left", color=INK, fontsize=10)
    save(fig, "fig6_rna_markers_at_protein_level.png")
    # keep any description already written after a figure's name in the README
    old = (OUT / "README.md").read_text() if (OUT / "README.md").exists() else ""
    desc = {m.group(1): m.group(2) for m in __import__("re").finditer(r"^- `([^`]+)`(.*)$", old, __import__("re").M)}
    (OUT / "README.md").write_text("# Figures for slides\n\nRendered by `scripts/multiomic/figures_png.py` from `site/data/multiomic.json` (the same numbers as `site/multiomic.html`; provenance in `site/data/provenance.json`, ids `mo_*`).\n\n"
                                   + "\n".join(f"- `{p.name}`{desc.get(p.name, '')}" for p in sorted(OUT.glob("*.png"))) + "\n")


if __name__ == "__main__":
    main()
