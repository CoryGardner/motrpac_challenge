"""Individual-level bootstrap intervals for transfer metrics (added for the overnight multiomic run).

Transfer accuracies and coverages are means over target samples that are clustered in individuals (donors, mice). The
spread reported beside them is a percentile bootstrap over INDIVIDUALS: individuals are resampled with replacement, all
their samples come along, the metric is recomputed. Nothing here is imported by the frozen pipeline.
"""
from __future__ import annotations

import numpy as np


def individual_bootstrap_ci(values, groups, n_boot: int = 1000, seed: int = 0, level: float = 0.95) -> tuple[float, float, int]:
    """Percentile CI of the mean of `values` under resampling of `groups` (individuals) with replacement.
    Returns (low, high, n_individuals). With one individual the interval is (mean, mean)."""
    values = np.asarray(values, dtype=float)
    groups = np.asarray(groups)
    ok = ~np.isnan(values)
    values, groups = values[ok], groups[ok]
    uniq = np.unique(groups)
    if len(values) == 0:
        return float("nan"), float("nan"), 0
    if len(uniq) < 2:
        m = float(values.mean())
        return m, m, int(len(uniq))
    idx_by_g = {u: np.flatnonzero(groups == u) for u in uniq}
    rng = np.random.default_rng(seed)
    means = np.empty(n_boot)
    for b in range(n_boot):
        chosen = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([idx_by_g[c] for c in chosen])
        means[b] = values[idx].mean()
    a = (1 - level) / 2
    return float(np.quantile(means, a)), float(np.quantile(means, 1 - a)), int(len(uniq))
