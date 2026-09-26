"""Baselines, tuned pipelines, grouped-CV evaluation, compact panels, and fusion.

Everything data-dependent (imputation, scaling, variance filter, feature selection, tuning)
lives inside an sklearn Pipeline so it is fit on the training fold only.
"""
from __future__ import annotations

import time
import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin, clone
from sklearn.ensemble import RandomForestClassifier
from sklearn.feature_selection import SelectKBest, f_classif
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import balanced_accuracy_score, f1_score, log_loss, roc_auc_score
from sklearn.model_selection import GridSearchCV, GroupKFold
from sklearn.neighbors import NearestCentroid
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

from . import config as C
from .io import OmicsMatrix
from .splits import assert_no_group_leak, grouped_kfold

BASELINES = ["centroid", "logreg_l1", "logreg_l2", "rf"]


# ---------------------------------------------------------------------------------------
# Transformers
# ---------------------------------------------------------------------------------------
class VarianceTopK(BaseEstimator, TransformerMixin):
    """Keep the k highest-variance features (fit on the training fold, no labels used)."""

    def __init__(self, k: int | None = 5000):
        self.k = k

    def fit(self, X, y=None):
        X = np.asarray(X, dtype=float)
        v = np.nanvar(X, axis=0)
        v = np.where(np.isnan(v), -np.inf, v)
        if self.k is None or self.k >= X.shape[1]:
            self.idx_ = np.arange(X.shape[1])
        else:
            self.idx_ = np.sort(np.argsort(v)[::-1][: self.k])
        return self

    def transform(self, X):
        return np.asarray(X, dtype=float)[:, self.idx_]

    def get_support_indices(self):
        return self.idx_


class RoundRobinSelector(BaseEstimator, TransformerMixin):
    """Class-aware panel selection: rank features per class by one-vs-rest effect size
    (|mean_in − mean_out| / pooled sd), then take the best remaining feature of each class in
    turn until k unique features are chosen.

    Why not SelectKBest(f_classif)? For multiclass problems the top-k F-scores tend to be
    k markers of the *same* easy class, so a k=20 "panel" can leave most classes
    unrepresented. The synthetic smoke test reproduces that failure (accuracy ~0.3 at k=20).
    """

    def __init__(self, k: int = 10, min_sd: float = 1e-6, exclude: np.ndarray | None = None):
        self.k = k
        self.min_sd = min_sd
        self.exclude = exclude  # (classes × features) bool, classes in np.unique(y) order

    def fit(self, X, y):
        X = np.asarray(X, dtype=float)
        y = np.asarray(y)
        self.classes_ = np.unique(y)
        scores = np.zeros((len(self.classes_), X.shape[1]))
        for i, c in enumerate(self.classes_):
            m_in, m_out = y == c, y != c
            if m_in.sum() < 2 or m_out.sum() < 2:
                continue
            mu_in, mu_out = np.nanmean(X[m_in], axis=0), np.nanmean(X[m_out], axis=0)
            sd = np.sqrt((np.nanvar(X[m_in], axis=0) + np.nanvar(X[m_out], axis=0)) / 2) + self.min_sd
            scores[i] = np.abs(mu_in - mu_out) / sd
        scores = np.nan_to_num(scores, nan=0.0)
        self.scores_ = scores
        self.order_ = roundrobin_order(scores, self.k, self.exclude)  # pick order; idx_ for any k' <= k is sorted(order_[:k'])
        self.idx_ = np.array(sorted(self.order_))
        return self

    def transform(self, X):
        return np.asarray(X, dtype=float)[:, self.idx_]

    def get_support(self, indices: bool = False):
        if indices:
            return self.idx_
        mask = np.zeros(self.scores_.shape[1], dtype=bool)
        mask[self.idx_] = True
        return mask


def roundrobin_order(scores: np.ndarray, k: int, exclude: np.ndarray | None = None) -> np.ndarray:
    """The order in which the round-robin selector picks features from a (classes × features) score
    matrix: the best remaining feature of each class in turn. The choice for a smaller k is the
    prefix of the choice for a larger k, so one call with k_max serves a whole panel-size grid.
    `exclude` (classes × features, bool) marks features a class may not pick — e.g. genes that are
    training-regulated or QC-correlated in that tissue (transferability-aware selection)."""
    if exclude is not None:
        exclude = np.asarray(exclude, dtype=bool)
        ranked = [[j for j in np.argsort(-scores[i]) if not exclude[i, j]] for i in range(scores.shape[0])]
    else:
        ranked = [list(np.argsort(-scores[i])) for i in range(scores.shape[0])]
    chosen: list[int] = []
    k = min(k, scores.shape[1])
    while len(chosen) < k:
        progressed = False
        for r in ranked:
            while r and r[0] in chosen:
                r.pop(0)
            if r and len(chosen) < k:
                chosen.append(r.pop(0))
                progressed = True
        if not progressed:
            break
    return np.array(chosen, dtype=int)


def make_selector(selector: str, k: int):
    if selector == "roundrobin":
        return RoundRobinSelector(k=k)
    if selector == "fclassif":
        return SelectKBest(f_classif, k=k)
    raise ValueError(f"unknown selector {selector}")


class BlockScaler(BaseEstimator, TransformerMixin):
    """Z-score features, then divide each block by sqrt(block size) so blocks contribute
    equally to distances/penalties regardless of feature count (for early fusion)."""

    def __init__(self, block_sizes: list[int] | None = None):
        self.block_sizes = block_sizes

    def fit(self, X, y=None):
        X = np.asarray(X, dtype=float)
        self.mean_ = np.nanmean(X, axis=0)
        self.std_ = np.nanstd(X, axis=0)
        self.std_[self.std_ == 0] = 1.0
        sizes = self.block_sizes or [X.shape[1]]
        w = np.concatenate([np.full(s, 1.0 / np.sqrt(s)) for s in sizes])
        self.w_ = w[: X.shape[1]]
        return self

    def transform(self, X):
        X = np.asarray(X, dtype=float)
        return (X - self.mean_) / self.std_ * self.w_


# ---------------------------------------------------------------------------------------
# Pipelines and grids
# ---------------------------------------------------------------------------------------
def make_pipeline(kind: str, k: int | None = None, prefilter: int | None = 5000,
                  seed: int = C.SEED, n_jobs: int = -1, selector: str = "roundrobin") -> Pipeline:
    """impute → variance prefilter → scale → panel selector(k) → classifier."""
    if kind == "centroid":
        clf = NearestCentroid()
    elif kind == "logreg_l1":
        clf = LogisticRegression(penalty="l1", solver="saga", C=0.1, max_iter=3000, tol=1e-3)
    elif kind == "logreg_l2":
        clf = LogisticRegression(penalty="l2", solver="lbfgs", C=0.1, max_iter=3000)
    elif kind == "rf":
        clf = RandomForestClassifier(n_estimators=500, max_features="sqrt", random_state=seed,
                                     n_jobs=n_jobs, class_weight="balanced")
    else:
        raise ValueError(f"unknown model kind {kind}")
    steps = [
        ("impute", SimpleImputer(strategy="median")),
        ("prefilter", VarianceTopK(prefilter)),
        ("scale", StandardScaler()),
        ("select", make_selector(selector, k) if k else "passthrough"),
        ("clf", clf),
    ]
    return Pipeline(steps)


def param_grid(kind: str, quick: bool = False) -> dict:
    if kind == "centroid":
        return {"clf__shrink_threshold": [None, 0.5]} if not quick else {}
    if kind in ("logreg_l1", "logreg_l2"):
        return {"clf__C": [0.01, 0.1, 1.0]} if not quick else {"clf__C": [0.1]}
    if kind == "rf":
        return {"clf__max_features": ["sqrt", 0.05]} if not quick else {}
    return {}


def fit_tuned(kind: str, X: np.ndarray, y: np.ndarray, groups: np.ndarray, k: int | None = None,
              prefilter: int | None = 5000, inner_splits: int = 3, quick: bool = False,
              seed: int = C.SEED, selector: str = "roundrobin"):
    """Fit a pipeline with an inner GroupKFold grid search (on the training fold only)."""
    pipe = make_pipeline(kind, k=k, prefilter=prefilter, seed=seed, selector=selector)
    grid = param_grid(kind, quick)
    n_groups = len(np.unique(groups))
    if not grid or n_groups < inner_splits + 1:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            return pipe.fit(X, y)
    gs = GridSearchCV(pipe, grid, cv=GroupKFold(n_splits=inner_splits), scoring="balanced_accuracy",
                      n_jobs=1 if kind == "rf" else -1, refit=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        gs.fit(X, y, groups=groups)
    return gs.best_estimator_


def has_proba(est) -> bool:
    return hasattr(est, "predict_proba")


# ---------------------------------------------------------------------------------------
# Metrics
# ---------------------------------------------------------------------------------------
def classification_metrics(y_true, y_pred, proba: np.ndarray | None, classes) -> dict:
    out = {
        "balanced_accuracy": balanced_accuracy_score(y_true, y_pred),
        "macro_f1": f1_score(y_true, y_pred, average="macro"),
        "accuracy": float(np.mean(np.asarray(y_true) == np.asarray(y_pred))),
    }
    if proba is not None:
        try:
            out["log_loss"] = log_loss(y_true, proba, labels=list(classes))
        except ValueError:
            out["log_loss"] = np.nan
        if len(classes) == 2:
            pos = list(classes)[1]
            out["auroc"] = roc_auc_score((np.asarray(y_true) == pos).astype(int), proba[:, 1])
    return out


# ---------------------------------------------------------------------------------------
# Grouped-CV evaluation
# ---------------------------------------------------------------------------------------
@dataclass
class CVResult:
    per_fold: pd.DataFrame          # one row per (model, fold)
    predictions: pd.DataFrame       # one row per (model, sample): y_true, y_pred, fold, proba cols
    confusion: dict[str, pd.DataFrame]

    def summary(self) -> pd.DataFrame:
        num = [c for c in self.per_fold.columns if c not in ("model", "fold", "fit_seconds")]
        agg = self.per_fold.groupby("model")[num].agg(["mean", "std"])
        agg.columns = [f"{a}_{b}" for a, b in agg.columns]
        agg["n_folds"] = self.per_fold.groupby("model")["fold"].nunique()
        return agg.reset_index()


def evaluate_cv(om: OmicsMatrix, label: str = "tissue", kinds: list[str] | None = None,
                n_splits: int = 5, k: int | None = None, prefilter: int | None = 5000,
                quick: bool = False, seed: int = C.SEED, verbose: bool = True,
                max_folds: int | None = None) -> CVResult:
    """Animal-grouped outer CV for each model kind; tuning is inside each training fold.
    max_folds limits the number of outer folds actually run (for timing a projection)."""
    kinds = kinds or BASELINES
    X = om.X.to_numpy(dtype=float)
    y = om.meta[label].astype(str).to_numpy()
    g = om.groups()
    classes = sorted(np.unique(y))
    rows, preds, conf = [], [], {}
    for fold, (tr, te) in enumerate(grouped_kfold(om.meta, label, n_splits, seed)):
        if max_folds is not None and fold >= max_folds:
            break
        assert_no_group_leak(om.meta, tr, te)
        for kind in kinds:
            t0 = time.time()
            est = fit_tuned(kind, X[tr], y[tr], g[tr], k=k, prefilter=prefilter, quick=quick, seed=seed)
            y_pred = est.predict(X[te])
            proba = est.predict_proba(X[te]) if has_proba(est) else None
            m = classification_metrics(y[te], y_pred, proba, list(est.classes_) if proba is not None else classes)
            m.update({"model": kind, "fold": fold, "n_train_animals": int(om.meta.iloc[tr]["pid"].nunique()),
                      "n_test_animals": int(om.meta.iloc[te]["pid"].nunique()),
                      "n_test": len(te), "fit_seconds": round(time.time() - t0, 1)})
            rows.append(m)
            p = pd.DataFrame({"model": kind, "fold": fold, "viallabel": om.meta.index[te],
                              "pid": g[te], "y_true": y[te], "y_pred": y_pred})
            if proba is not None:
                for j, c in enumerate(est.classes_):
                    p[f"p_{c}"] = proba[:, j]
            preds.append(p)
            if verbose:
                print(f"  fold {fold} {kind:10s} bal.acc={m['balanced_accuracy']:.3f} "
                      f"macroF1={m['macro_f1']:.3f} ({m['fit_seconds']}s, {m['n_test_animals']} test animals)")
    per_fold = pd.DataFrame(rows)
    predictions = pd.concat(preds, ignore_index=True)
    for kind in kinds:
        pk = predictions[predictions["model"] == kind]
        conf[kind] = pd.crosstab(pk["y_true"], pk["y_pred"]).reindex(index=classes, columns=classes, fill_value=0)
    return CVResult(per_fold, predictions, conf)


# ---------------------------------------------------------------------------------------
# Compact panels
# ---------------------------------------------------------------------------------------
def panel_curve(om: OmicsMatrix, label: str = "tissue", grid: list[int] | None = None,
                kind: str = "logreg_l2", n_splits: int = 5, prefilter: int | None = 5000,
                quick: bool = False, seed: int = C.SEED, verbose: bool = True,
                selector: str = "roundrobin", return_predictions: bool = False):
    """Accuracy vs panel size k. Returns (per-fold table, selected features per fold/k), plus the
    per-sample predictions table (fold, k, viallabel, y_true, y_pred) when return_predictions=True."""
    grid = grid or C.PANEL_GRID
    X = om.X.to_numpy(dtype=float)
    y = om.meta[label].astype(str).to_numpy()
    g = om.groups()
    feat_ids = np.asarray(om.X.columns)
    rows, sel, preds = [], [], []
    for fold, (tr, te) in enumerate(grouped_kfold(om.meta, label, n_splits, seed)):
        for k in grid:
            if k > X.shape[1]:
                continue
            est = fit_tuned(kind, X[tr], y[tr], g[tr], k=k, prefilter=prefilter, quick=quick, seed=seed,
                            selector=selector)
            y_pred = est.predict(X[te])
            proba = est.predict_proba(X[te]) if has_proba(est) else None
            m = classification_metrics(y[te], y_pred, proba, list(est.classes_) if proba is not None else None)
            m.update({"k": k, "fold": fold, "model": kind, "selector": selector,
                      "n_train_animals": int(om.meta.iloc[tr]["pid"].nunique()),
                      "n_test_animals": int(om.meta.iloc[te]["pid"].nunique())})
            rows.append(m)
            chosen = selected_feature_ids(est, feat_ids)
            sel.extend({"fold": fold, "k": k, "feature_ID": f} for f in chosen)
            preds.append(pd.DataFrame({"fold": fold, "k": k, "viallabel": om.meta.index[te], "pid": g[te],
                                       "y_true": y[te], "y_pred": y_pred}))
            if verbose:
                print(f"  fold {fold} k={k:4d} bal.acc={m['balanced_accuracy']:.3f}")
    if return_predictions:
        return pd.DataFrame(rows), pd.DataFrame(sel), pd.concat(preds, ignore_index=True)
    return pd.DataFrame(rows), pd.DataFrame(sel)


def selected_feature_ids(est: Pipeline, feat_ids: np.ndarray) -> list[str]:
    """Trace a fitted pipeline's prefilter + SelectKBest back to feature IDs."""
    idx = np.arange(len(feat_ids))
    pre = est.named_steps.get("prefilter")
    if pre is not None and hasattr(pre, "idx_"):
        idx = idx[pre.idx_]
    sel = est.named_steps.get("select")
    if sel is not None and sel != "passthrough" and hasattr(sel, "get_support"):
        idx = idx[sel.get_support()]
    return list(feat_ids[idx])


def stability_selection(om: OmicsMatrix, label: str = "tissue", k: int = 10, n_boot: int = 50,
                        prefilter: int | None = 5000, seed: int = C.SEED,
                        selector: str = "roundrobin") -> pd.DataFrame:
    """Selection frequency of each feature over bootstrap resamples of *animals*."""
    rng = np.random.default_rng(seed)
    X = om.X.to_numpy(dtype=float)
    y = om.meta[label].astype(str).to_numpy()
    pids = om.meta["pid"].astype(str).to_numpy()
    uniq = np.unique(pids)
    counts = pd.Series(0, index=om.X.columns, dtype=float)
    for _ in range(n_boot):
        chosen = rng.choice(uniq, size=len(uniq), replace=True)
        idx = np.concatenate([np.flatnonzero(pids == p) for p in chosen])
        if len(np.unique(y[idx])) < 2:
            continue
        pipe = Pipeline([("impute", SimpleImputer(strategy="median")), ("prefilter", VarianceTopK(prefilter)),
                         ("scale", StandardScaler()), ("select", make_selector(selector, k))])
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            pipe.fit(X[idx], y[idx])
        for f in selected_feature_ids(pipe, np.asarray(om.X.columns)):
            counts[f] += 1
    freq = (counts / n_boot).sort_values(ascending=False)
    out = freq[freq > 0].rename("selection_frequency").reset_index().rename(columns={"index": "feature_ID"})
    return out


# ---------------------------------------------------------------------------------------
# Fusion
# ---------------------------------------------------------------------------------------
def early_fusion_matrix(blocks: dict[str, pd.DataFrame]) -> tuple[pd.DataFrame, list[int]]:
    """Concatenate aligned blocks (same row order). Returns (X, block_sizes)."""
    names = list(blocks)
    idx = blocks[names[0]].index
    for n in names[1:]:
        if not blocks[n].index.equals(idx):
            raise ValueError(f"block {n} is not aligned with {names[0]} — use align_by_animal first")
    X = pd.concat([blocks[n].add_prefix(f"{n}::") for n in names], axis=1)
    return X, [blocks[n].shape[1] for n in names]


def make_fusion_pipeline(kind: str, block_sizes: list[int], seed: int = C.SEED) -> Pipeline:
    base = make_pipeline(kind, k=None, prefilter=None, seed=seed)
    clf = base.named_steps["clf"]
    return Pipeline([("impute", SimpleImputer(strategy="median")),
                     ("scale", BlockScaler(block_sizes)), ("clf", clone(clf))])


def late_fusion_proba(estimators: list, Xs: list[np.ndarray]) -> tuple[np.ndarray, np.ndarray]:
    """Average predict_proba across per-block estimators (all must share classes_)."""
    classes = estimators[0].classes_
    P = np.zeros((Xs[0].shape[0], len(classes)))
    for est, X in zip(estimators, Xs):
        p = est.predict_proba(X)
        order = [list(est.classes_).index(c) for c in classes]
        P += p[:, order]
    P /= len(estimators)
    return P, classes
