"""Split-conformal prediction sets for classification, and an LTT-style certified panel size.

References: Angelopoulos & Bates (2021) for split conformal / LAC / APS; Angelopoulos et al.
(2021) "Learn then Test" and Bates et al. (2021) RCPS for the risk-control view.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy import stats

from . import config as C


# ---------------------------------------------------------------------------------------
# Nonconformity scores
# ---------------------------------------------------------------------------------------
def lac_scores(proba: np.ndarray, y_idx: np.ndarray) -> np.ndarray:
    """1 - p̂(true class). Smallest average set size among simple scores."""
    return 1.0 - proba[np.arange(len(y_idx)), y_idx]


def aps_scores(proba: np.ndarray, y_idx: np.ndarray, rng: np.random.Generator | None = None) -> np.ndarray:
    """Adaptive prediction sets: cumulative mass of classes ranked above the true class
    (+ randomized share of the true class's mass). Better conditional coverage, larger sets."""
    rng = rng or np.random.default_rng(C.SEED)
    order = np.argsort(-proba, axis=1)
    sorted_p = np.take_along_axis(proba, order, axis=1)
    cum = np.cumsum(sorted_p, axis=1)
    ranks = np.argmax(order == y_idx[:, None], axis=1)
    cum_before = cum[np.arange(len(y_idx)), ranks] - sorted_p[np.arange(len(y_idx)), ranks]
    u = rng.uniform(size=len(y_idx))
    return cum_before + u * sorted_p[np.arange(len(y_idx)), ranks]


def conformal_quantile(scores: np.ndarray, alpha: float) -> float:
    """Textbook split-conformal threshold: the ceil((n+1)(1-α))-th smallest of the n calibration
    scores. When that rank exceeds n (too few calibration points for level α) there is no finite
    threshold: return +inf, so every class enters every set (LAC and APS alike).

    Fixed 2026-09-25. The previous rule, np.quantile(scores, min(ceil((n+1)(1-α))/n, 1), "higher"),
    returned the largest score instead of +inf when the rank exceeded n, and one rank too high
    otherwise (FINDINGS_REPORT §14, item 13(a))."""
    n = len(scores)
    if n == 0:
        return float("inf")
    k = int(np.ceil((n + 1) * (1 - alpha) - 1e-9))   # tolerance: (n+1)(1-α) can be an integer + ε
    if k > n:
        return float("inf")
    return float(np.sort(np.asarray(scores, dtype=float))[max(k, 1) - 1])


def predict_sets(proba: np.ndarray, qhat: float, method: str = "lac") -> np.ndarray:
    """Boolean matrix (n × classes): class j is in sample i's set."""
    if method == "lac":
        return proba >= (1.0 - qhat)
    # APS: include classes in descending order until cumulative mass exceeds qhat
    order = np.argsort(-proba, axis=1)
    sorted_p = np.take_along_axis(proba, order, axis=1)
    cum = np.cumsum(sorted_p, axis=1)
    include_sorted = (cum - sorted_p) <= qhat  # classes whose "cum before" is within budget
    sets = np.zeros_like(proba, dtype=bool)
    np.put_along_axis(sets, order, include_sorted, axis=1)
    return sets


def class_indices(y: np.ndarray, classes: list[str]) -> tuple[np.ndarray, np.ndarray]:
    """Map labels to column positions of a fitted estimator's `classes_`.

    Returns (idx, seen): idx is -1 where the label is not among `classes`, and `seen` marks the
    rows that have a probability column at all. Never index a probability matrix with idx
    without masking on `seen`: -1 would silently pick the last class.
    """
    pos = {c: i for i, c in enumerate(classes)}
    idx = np.array([pos.get(v, -1) for v in np.asarray(y)], dtype=int)
    return idx, idx >= 0


def calibration_scores(proba: np.ndarray, y: np.ndarray, classes: list[str], method: str = "lac",
                       rng: np.random.Generator | None = None) -> tuple[np.ndarray, int]:
    """Nonconformity scores for the calibration samples whose label the estimator has seen.

    Samples whose label is absent from `classes` (a class missing from the fit animals) are
    DROPPED and counted, never scored: they have no probability column, and the conformal
    quantile must come from scores that are exchangeable with those of coverable test points.
    Returns (scores over seen samples, number of dropped samples).
    """
    proba = np.asarray(proba, dtype=float)
    if proba.ndim != 2 or proba.shape[1] != len(classes):
        raise ValueError(f"proba has {proba.shape[1] if proba.ndim == 2 else '?'} columns, "
                         f"but {len(classes)} classes were given")
    idx, seen = class_indices(y, classes)
    p, yi = proba[seen], idx[seen]
    scores = lac_scores(p, yi) if method == "lac" else aps_scores(p, yi, rng)
    return scores, int((~seen).sum())


def coverage_seen(sets: np.ndarray, y: np.ndarray, classes: list[str]) -> tuple[float, int]:
    """Empirical coverage over samples whose label is among `classes`; the others are dropped
    and counted (they cannot be covered by any set). Returns (coverage, number dropped)."""
    idx, seen = class_indices(y, classes)
    if not seen.any():
        return float("nan"), int(len(idx))
    covered = sets[np.flatnonzero(seen), idx[seen]]
    return float(covered.mean()), int((~seen).sum())


def conformal_quantile_per_class(scores: np.ndarray, y_idx: np.ndarray, alpha: float, n_classes: int,
                                 fallback: float | None = None, floor: float | None = None) -> np.ndarray:
    """Mondrian (class-conditional) conformal: one finite-sample quantile per true class, from the
    calibration samples of that class. A class with no calibration sample gets `fallback` (use the
    marginal quantile) or NaN. Coverage is then guaranteed within each class, at the price of
    larger sets for the classes the model confuses.

    `floor` (marginal-floor Mondrian): every class quantile is raised to at least `floor` (pass the
    marginal quantile). The resulting sets contain both the marginal and the plain Mondrian sets,
    so the marginal 1 − α guarantee is kept alongside the per-class one; the price is paid in set
    size on the classes whose own quantile was below the marginal one (the easy classes, whose
    Mondrian threshold is tighter than the marginal because a few well-classified calibration
    samples make a small per-class quantile). With few calibration samples per class the plain
    Mondrian quantile is noisy in both directions; the floor removes the downward noise."""
    q = np.full(n_classes, np.nan)
    for c in range(n_classes):
        sc = scores[y_idx == c]
        if len(sc):
            q[c] = conformal_quantile(sc, alpha)
    if fallback is not None:
        q[np.isnan(q)] = fallback
    if floor is not None:
        q = np.fmax(q, float(floor))
    return q


def predict_sets_conditional(proba: np.ndarray, qhat_per_class: np.ndarray, method: str = "lac") -> np.ndarray:
    """Boolean matrix (n × classes) with a class-specific threshold: class j is in sample i's set
    when its class-j nonconformity is within class j's quantile."""
    q = np.asarray(qhat_per_class, dtype=float)[None, :]
    if method == "lac":
        return proba >= (1.0 - q)
    order = np.argsort(-proba, axis=1)
    sorted_p = np.take_along_axis(proba, order, axis=1)
    cum_before_sorted = np.cumsum(sorted_p, axis=1) - sorted_p
    cum_before = np.empty_like(proba)
    np.put_along_axis(cum_before, order, cum_before_sorted, axis=1)
    return cum_before <= q


def evaluate_sets(sets: np.ndarray, y_idx: np.ndarray, classes: list[str],
                  extra: pd.DataFrame | None = None) -> tuple[dict, pd.DataFrame]:
    """Coverage / size overall and per true class (conditional coverage)."""
    covered = sets[np.arange(len(y_idx)), y_idx]
    size = sets.sum(axis=1)
    overall = {"coverage": float(covered.mean()), "avg_set_size": float(size.mean()),
               "frac_singleton": float((size == 1).mean()), "frac_empty": float((size == 0).mean()),
               "n": int(len(y_idx))}
    df = pd.DataFrame({"y_true": [classes[i] for i in y_idx], "covered": covered, "size": size})
    if extra is not None:
        df = pd.concat([df.reset_index(drop=True), extra.reset_index(drop=True)], axis=1)
    per_class = df.groupby("y_true").agg(coverage=("covered", "mean"), avg_set_size=("size", "mean"),
                                         n=("covered", "size")).reset_index()
    return overall, per_class


# ---------------------------------------------------------------------------------------
# One-sample-per-animal calibration (exchangeability hygiene)
# ---------------------------------------------------------------------------------------
def one_per_group(idx: np.ndarray, groups: np.ndarray, seed: int = C.SEED) -> np.ndarray:
    """Subsample idx so each group (animal) contributes one calibration point."""
    rng = np.random.default_rng(seed)
    g = groups[idx]
    out = []
    for u in np.unique(g):
        members = idx[g == u]
        out.append(rng.choice(members))
    return np.array(sorted(out))


# ---------------------------------------------------------------------------------------
# Split-conformal wrapper
# ---------------------------------------------------------------------------------------
@dataclass
class ConformalResult:
    alpha: float
    method: str
    qhat: float
    n_cal: int
    overall: dict
    per_class: pd.DataFrame
    qhat_per_class: np.ndarray | None = None


CONDITIONAL_VARIANTS = ("marginal", "mondrian", "floored")


def calibrate_and_evaluate(proba_cal: np.ndarray, y_cal_idx: np.ndarray, proba_test: np.ndarray,
                           y_test_idx: np.ndarray, classes: list[str], alpha: float,
                           method: str = "lac", extra_test: pd.DataFrame | None = None,
                           conditional: bool | str = False) -> ConformalResult:
    """Split conformal. conditional=True or "mondrian" uses Mondrian (per-class) quantiles with the
    marginal quantile as fallback for classes absent from the calibration set; "floored" is the
    marginal-floor Mondrian variant (per-class quantiles never below the marginal one)."""
    variant = {False: "marginal", True: "mondrian"}.get(conditional, conditional)
    if variant not in CONDITIONAL_VARIANTS:
        raise ValueError(f"conditional must be one of {CONDITIONAL_VARIANTS} (or a bool), got {conditional!r}")
    scores = lac_scores(proba_cal, y_cal_idx) if method == "lac" else aps_scores(proba_cal, y_cal_idx)
    qhat = conformal_quantile(scores, alpha)
    if variant != "marginal":
        qs = conformal_quantile_per_class(scores, y_cal_idx, alpha, len(classes), fallback=qhat,
                                          floor=qhat if variant == "floored" else None)
        sets = predict_sets_conditional(proba_test, qs, method)
        overall, per_class = evaluate_sets(sets, y_test_idx, classes, extra_test)
        return ConformalResult(alpha, method, float(np.nanmean(qs)), len(scores), overall, per_class, qs)
    sets = predict_sets(proba_test, qhat, method)
    overall, per_class = evaluate_sets(sets, y_test_idx, classes, extra_test)
    return ConformalResult(alpha, method, qhat, len(scores), overall, per_class)


# ---------------------------------------------------------------------------------------
# Certified panel size (fixed-sequence testing, LTT flavour)
# ---------------------------------------------------------------------------------------
def clopper_pearson_upper(n_err: int, n: int, delta: float) -> float:
    """Exact one-sided (1-δ) upper confidence bound on a binomial error rate."""
    if n == 0:
        return 1.0
    if n_err >= n:
        return 1.0
    return float(stats.beta.ppf(1 - delta, n_err + 1, n - n_err))


def hoeffding_upper(err_hat: float, n: int, delta: float) -> float:
    if n == 0:
        return 1.0
    return float(min(1.0, err_hat + np.sqrt(np.log(1 / delta) / (2 * n))))


def certify_panel_size(errors: dict[int, tuple[int, int]], alpha: float, delta: float,
                       bound: str = "clopper_pearson") -> tuple[int | None, pd.DataFrame]:
    """errors: {k: (n_errors_on_calibration, n_calibration)}.

    Fixed-sequence testing from the largest k downward: certify k while UCB(err_k) <= alpha,
    stop at the first failure. Returns (smallest certified k or None, table).
    The order must be fixed before looking at the data — it is (descending k).
    """
    ks = sorted(errors, reverse=True)
    rows, certified, stopped = [], None, False
    for k in ks:
        n_err, n = errors[k]
        err_hat = n_err / n if n else np.nan
        ucb = clopper_pearson_upper(n_err, n, delta) if bound == "clopper_pearson" else hoeffding_upper(err_hat, n, delta)
        passes = (ucb <= alpha) and not stopped
        if passes:
            certified = k
        else:
            stopped = True
        rows.append({"k": k, "n_cal": n, "n_err": n_err, "err_hat": err_hat, "ucb": ucb,
                     "certified": passes})
    table = pd.DataFrame(rows).sort_values("k").reset_index(drop=True)
    return certified, table
