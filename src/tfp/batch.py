"""Batch nesting: how processing variables (extraction plates, library batches, flowcells, TMT plexes,
platforms) relate to tissue, and which tissue pairs stay *estimable*.

A pair of tissues is estimable when the two tissues share a level of EVERY processing variable given,
so that a contrast between them exists inside one batch. When each tissue was processed as its own
batch, no such contrast exists and the tissue effect cannot be separated from the batch effect
without an external dataset.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd


def cramers_v(a, b) -> float:
    """Cramér's V between two categorical series (0 = independent, 1 = one determines the other).
    A variable with a single level has no association: 0."""
    ct = pd.crosstab(pd.Series(a).astype(str).to_numpy(), pd.Series(b).astype(str).to_numpy())
    if ct.shape[0] < 2 or ct.shape[1] < 2:
        return 0.0
    obs = ct.to_numpy(dtype=float)
    n = obs.sum()
    exp = np.outer(obs.sum(axis=1), obs.sum(axis=0)) / n
    chi2 = float(((obs - exp) ** 2 / exp).sum())
    return float(np.sqrt(chi2 / (n * (min(ct.shape) - 1))))


def _shared_pairs(meta: pd.DataFrame, tissue_col: str, variable: str) -> list[tuple[str, str]]:
    ct = pd.crosstab(meta[tissue_col].astype(str), meta[variable].astype(str)) > 0
    tissues = sorted(ct.index)
    return [(a, b) for a, b in itertools.combinations(tissues, 2) if bool((ct.loc[a] & ct.loc[b]).any())]


def estimable_pairs(meta: pd.DataFrame, tissue_col: str, variables: list[str]) -> tuple[list[tuple[str, str]], int]:
    """Tissue pairs that share a level of every variable in `variables`, and the number of pairs in total."""
    tissues = sorted(meta[tissue_col].astype(str).unique())
    total = len(tissues) * (len(tissues) - 1) // 2
    pairs: set[tuple[str, str]] | None = None
    for v in variables:
        m = meta[[tissue_col, v]].dropna()
        s = set(_shared_pairs(m, tissue_col, v))
        pairs = s if pairs is None else pairs & s
    return sorted(pairs or []), total


def nesting_table(meta: pd.DataFrame, tissue_col: str, variables: list[str]) -> pd.DataFrame:
    """One row per processing variable: how many levels it has, how many levels a tissue spans, how many
    tissues a level holds, how many levels are shared by more than one tissue, Cramér's V with tissue, and
    how many tissue pairs share at least one level (rows with a missing value are dropped per variable)."""
    rows = []
    for v in variables:
        m = meta[[tissue_col, v]].dropna()
        ct = pd.crosstab(m[tissue_col].astype(str), m[v].astype(str)) > 0
        levels_per_tissue = ct.sum(axis=1)
        tissues_per_level = ct.sum(axis=0)
        n_t = int(ct.shape[0])
        rows.append({"variable": v, "n_levels": int(ct.shape[1]), "n_tissues": n_t,
                     "median_levels_per_tissue": float(levels_per_tissue.median()),
                     "max_levels_per_tissue": int(levels_per_tissue.max()),
                     "max_tissues_per_level": int(tissues_per_level.max()),
                     "n_levels_shared": int((tissues_per_level > 1).sum()),
                     "tissues_in_one_level": int((levels_per_tissue == 1).sum()),
                     "cramers_v": cramers_v(m[tissue_col], m[v]),
                     "n_pairs_sharing_level": len(_shared_pairs(m, tissue_col, v)),
                     "n_pairs_total": n_t * (n_t - 1) // 2})
    return pd.DataFrame(rows)
