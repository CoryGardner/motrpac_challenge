"""Matplotlib-only figures. One idea per figure; PNG at 150 dpi."""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402

from . import config as C  # noqa: E402

DPI = 150


def _save(fig, path: Path):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)
    return path


def pca_scatter(scores: np.ndarray, meta: pd.DataFrame, color: str, shape: str | None,
                title: str, path: Path, explained: tuple[float, float] | None = None,
                pcs: tuple[int, int] = (0, 1)) -> Path:
    fig, ax = plt.subplots(figsize=(7, 6))
    cats = sorted(meta[color].astype(str).unique())
    cmap = plt.get_cmap("tab20" if len(cats) > 10 else "tab10")
    markers = {}
    if shape:
        for i, s in enumerate(sorted(meta[shape].astype(str).unique())):
            markers[s] = ["o", "^", "s", "D", "v", "P", "X"][i % 7]
    for i, c in enumerate(cats):
        m = (meta[color].astype(str) == c).to_numpy()
        if shape:
            for s, mk in markers.items():
                mm = m & (meta[shape].astype(str) == s).to_numpy()
                if mm.any():
                    ax.scatter(scores[mm, pcs[0]], scores[mm, pcs[1]], s=22, color=cmap(i % 20), marker=mk,
                               edgecolor="k", linewidth=0.3, label=c if s == list(markers)[0] else None)
        else:
            ax.scatter(scores[m, pcs[0]], scores[m, pcs[1]], s=22, color=cmap(i % 20), edgecolor="k", linewidth=0.3, label=c)
    xl = f"PC{pcs[0]+1}" + (f" ({explained[0]:.0%})" if explained else "")
    yl = f"PC{pcs[1]+1}" + (f" ({explained[1]:.0%})" if explained else "")
    ax.set_xlabel(xl)
    ax.set_ylabel(yl)
    ax.set_title(title + (f"  [marker = {shape}]" if shape else ""))
    ax.legend(fontsize=7, ncol=2, frameon=False, bbox_to_anchor=(1.02, 1), loc="upper left")
    return _save(fig, path)


def confusion_heatmap(conf: pd.DataFrame, title: str, path: Path, normalize: bool = True) -> Path:
    M = conf.to_numpy(dtype=float)
    if normalize:
        M = M / np.maximum(M.sum(axis=1, keepdims=True), 1)
    fig, ax = plt.subplots(figsize=(0.45 * len(conf) + 3, 0.45 * len(conf) + 2.5))
    im = ax.imshow(M, cmap="Blues", vmin=0, vmax=1 if normalize else None)
    ax.set_xticks(range(len(conf.columns)))
    ax.set_xticklabels(conf.columns, rotation=90, fontsize=7)
    ax.set_yticks(range(len(conf.index)))
    ax.set_yticklabels(conf.index, fontsize=7)
    ax.set_xlabel("predicted")
    ax.set_ylabel("true")
    ax.set_title(title)
    for i in range(M.shape[0]):
        for j in range(M.shape[1]):
            if conf.iat[i, j] > 0 and i != j:
                ax.text(j, i, int(conf.iat[i, j]), ha="center", va="center", fontsize=6, color="red")
    fig.colorbar(im, ax=ax, fraction=0.03)
    return _save(fig, path)


def panel_curve_plot(per_fold: pd.DataFrame, metric: str, title: str, path: Path,
                     certified_k: int | None = None) -> Path:
    agg = per_fold.groupby("k")[metric].agg(["mean", "std"]).reset_index()
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.errorbar(agg["k"], agg["mean"], yerr=agg["std"].fillna(0), marker="o", capsize=3)
    ax.set_xscale("log")
    ax.set_xlabel("panel size k (features)")
    ax.set_ylabel(metric.replace("_", " "))
    ax.set_title(title)
    if certified_k:
        ax.axvline(certified_k, color="red", ls="--", label=f"certified k = {certified_k}")
        ax.legend(frameon=False)
    ax.grid(alpha=0.3)
    return _save(fig, path)


def certificate_plot(table: pd.DataFrame, alpha: float, title: str, path: Path) -> Path:
    fig, ax = plt.subplots(figsize=(6.5, 4.2))
    ax.plot(table["k"], table["err_hat"], marker="o", label="calibration error (point est.)")
    ax.plot(table["k"], table["ucb"], marker="s", label="upper confidence bound")
    ax.axhline(alpha, color="red", ls="--", label=f"α = {alpha}")
    cert = table[table["certified"]]
    if len(cert):
        ax.axvline(cert["k"].min(), color="green", ls=":", label=f"smallest certified k = {int(cert['k'].min())}")
    ax.set_xscale("log")
    ax.set_xlabel("panel size k")
    ax.set_ylabel("error rate")
    ax.set_title(title)
    ax.legend(frameon=False, fontsize=8)
    ax.grid(alpha=0.3)
    return _save(fig, path)


def coverage_bars(df: pd.DataFrame, x: str, title: str, path: Path, alpha: float | None = None) -> Path:
    fig, ax = plt.subplots(figsize=(max(5, 0.4 * len(df) + 2), 4))
    ax.bar(df[x].astype(str), df["coverage"], color="steelblue")
    if alpha is not None:
        ax.axhline(1 - alpha, color="red", ls="--", label=f"target 1−α = {1-alpha:.2f}")
        ax.legend(frameon=False)
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("empirical coverage")
    ax.set_title(title)
    ax.tick_params(axis="x", rotation=90, labelsize=7)
    return _save(fig, path)


def histogram(values: np.ndarray, title: str, xlabel: str, path: Path, bins: int = 40) -> Path:
    fig, ax = plt.subplots(figsize=(6, 4))
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    ax.hist(v, bins=bins, color="gray", edgecolor="white")
    ax.axvline(np.median(v), color="red", ls="--", label=f"median = {np.median(v):.2f}")
    ax.set_xlabel(xlabel)
    ax.set_ylabel("count")
    ax.set_title(title)
    ax.legend(frameon=False)
    return _save(fig, path)


def missingness_heatmap(present: pd.DataFrame, title: str, path: Path) -> Path:
    """present: tissues × assays (values = n samples)."""
    fig, ax = plt.subplots(figsize=(1.1 * present.shape[1] + 3, 0.4 * present.shape[0] + 2))
    im = ax.imshow(present.to_numpy(dtype=float), cmap="viridis")
    ax.set_xticks(range(present.shape[1]))
    ax.set_xticklabels(present.columns, rotation=45, ha="right")
    ax.set_yticks(range(present.shape[0]))
    ax.set_yticklabels(present.index, fontsize=8)
    for i in range(present.shape[0]):
        for j in range(present.shape[1]):
            v = present.iat[i, j]
            if v > 0:
                ax.text(j, i, int(v), ha="center", va="center", fontsize=7, color="white")
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.03, label="n samples")
    return _save(fig, path)
