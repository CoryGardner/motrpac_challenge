#!/usr/bin/env python
"""Figure for parts 5/6 (reads the CSVs written by tc56_fingerprint_shift.py and tc6_fingerprint_seeds.py)."""
from __future__ import annotations

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from motrpac import config as C

OUT = C.RESULTS_DIR / "15_time_course" / "5_6_fingerprint"
COL = {"full": "#2a78d6", "panel_k20": "#eb6834"}
LAB = {"full": "full model", "panel_k20": "k = 20 panel"}


def main() -> None:
    p5 = pd.read_csv(OUT / "part5_by_duration.csv")
    ms = pd.read_csv(OUT / "part5_matched_summary.csv")
    p6 = pd.read_csv(OUT / "part6_extrapolation.csv")
    s6 = pd.read_csv(OUT / "part6_seeds.csv")
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
    order = ["control", "1w", "2w", "4w", "8w"]
    x = np.arange(len(order))
    for ax, metric, title in ((axes[0], "coverage", "A  Conformal coverage (α = 0.10), control-trained model"),
                              (axes[1], "accuracy", "B  Tissue accuracy, control-trained model")):
        for j, arm in enumerate(("full", "panel_k20")):
            off = (j - 0.5) * 0.22
            c = COL[arm]
            # held-out controls (matched folds, pooled) as the in-distribution reference
            hc = ms[(ms["arm"] == arm) & ms["test_group"].str.startswith("control_heldout")].iloc[0]
            ax.errorbar(0 + off, hc[metric], yerr=[[hc[metric] - hc[f"{metric}_lo"]], [hc[f"{metric}_hi"] - hc[metric]]],
                        fmt="s", mfc="white", color=c, ms=8, lw=2, capsize=3)
            d = p5[(p5["design"] == "control_only_fit7_cal3") & (p5["arm"] == arm) & p5["test_group"].isin(order[1:])]
            d = d.set_index("test_group").loc[order[1:]]
            ax.errorbar(x[1:] + off, d[metric], yerr=[d[metric] - d[f"{metric}_lo"], d[f"{metric}_hi"] - d[metric]],
                        fmt="o-", color=c, ms=8, lw=2, capsize=3, label=f"{LAB[arm]} (fit 7 / cal 3 controls)")
            dx = p5[(p5["design"].str.contains("excl")) & (p5["arm"] == arm) & p5["test_group"].isin(["1w", "2w"])]
            dx = dx.set_index("test_group").loc[["1w", "2w"]]
            ax.plot(x[1:3] + off, dx[metric], "x", color=c, ms=9, mew=2,
                    label=f"{LAB[arm]}, excl. consortium-flagged VENACV vials" if metric == "accuracy" else None)
        if metric == "coverage":
            ax.axhline(0.9, color="#888", lw=1, ls="--")
            ax.text(4.35, 0.901, "1 − α", va="bottom", ha="right", fontsize=9, color="#555")
        ax.set_xticks(x, ["held-out\ncontrols*", "1w†", "2w†", "4w†", "8w"])
        ax.set_xlabel("test group (10 animals each, ~180 vials)")
        ax.set_ylabel(metric + " (vials)")
        ax.set_title(title, fontsize=10, loc="left")
        ax.grid(axis="y", color="#e5e5e5"); ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylim(0.78, 1.0); axes[1].set_ylim(0.84, 1.005)
    axes[0].legend(fontsize=8, loc="lower left", frameon=False)
    axes[1].legend(fontsize=8, loc="lower right", frameon=False)
    ax = axes[2]
    des = ["d_control_only", "a_control+1w", "c_control+4w", "b_control+1w+2w+4w"]
    lbl = ["control\n(10)", "control+1w\n(20)", "control+4w\n(20)", "control+1w\n+2w+4w (40)"]
    for j, arm in enumerate(("full", "panel_k20")):
        off = (j - 0.5) * 0.22
        d = s6[s6["arm"] == arm]
        for i, dn in enumerate(des):
            v = d[d["design"] == dn]["coverage"].to_numpy()
            ax.scatter(np.full(len(v), i + off) + np.linspace(-0.05, 0.05, len(v)), v, s=14, color=COL[arm], alpha=0.35, lw=0)
        m = d.groupby("design")["coverage"].agg(["mean", "std"]).loc[des]
        ax.errorbar(np.arange(4) + off, m["mean"], yerr=m["std"], fmt="D", color=COL[arm], ms=7, lw=2, capsize=3,
                    label=f"{LAB[arm]}: mean ± sd over 10 fit/cal splits")
    ax.axhline(0.9, color="#888", lw=1, ls="--")
    ax.set_xticks(np.arange(4), lbl, fontsize=8)
    ax.set_xlabel("fit+cal groups (n animals); test: 10 held-out 8w animals")
    ax.set_ylabel("coverage on 8w vials")
    ax.set_ylim(0.84, 0.98)
    ax.set_title("C  Held-out 8w: calibration source vs coverage", fontsize=10, loc="left")
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    ax.grid(axis="y", color="#e5e5e5"); ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.text(0.01, 0.005, "Bars: 95% cluster-bootstrap over test animals (conditional on one fit/cal split). "
             "* held-out controls: 5 matched folds (fit 6 / cal 2), pooled. "
             "† 1w/2w/4w differ from controls in arrival cohort, sacrifice season and date as well as training; only 8w is cohort-matched.",
             fontsize=7.5, color="#444")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    fig.savefig(OUT / "fingerprint_time_course.png", dpi=150)
    print(OUT / "fingerprint_time_course.png")


if __name__ == "__main__":
    main()
