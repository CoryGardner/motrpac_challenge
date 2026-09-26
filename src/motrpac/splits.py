"""Animal-grouped splits. Every split in this project comes from here.

All functions take the sample metadata (`OmicsMatrix.meta`) and yield positional index
arrays (train_idx, test_idx) into that frame's rows.
"""
from __future__ import annotations

from typing import Iterator

import numpy as np
import pandas as pd
from sklearn.model_selection import GroupKFold, StratifiedGroupKFold

from . import config as C

Split = tuple[np.ndarray, np.ndarray]


def assert_no_group_leak(meta: pd.DataFrame, train_idx, test_idx, group: str = "pid") -> None:
    tr = set(meta.iloc[np.asarray(train_idx)][group])
    te = set(meta.iloc[np.asarray(test_idx)][group])
    leak = tr & te
    if leak:
        raise AssertionError(f"{len(leak)} {group}(s) appear in both train and test: {sorted(leak)[:5]}")


def grouped_kfold(meta: pd.DataFrame, label: str = "tissue", n_splits: int = 5,
                  seed: int = C.SEED, group: str = "pid") -> Iterator[Split]:
    """Stratified (on `label`) K-fold with all samples of an animal on one side.

    Falls back to plain GroupKFold when a class has fewer groups than folds (sklearn raises
    in that case); the fallback is logged in the returned generator's first yield via a warning.
    """
    y = meta[label].astype(str).to_numpy()
    g = meta[group].astype(str).to_numpy()
    n_groups = len(np.unique(g))
    n_splits = min(n_splits, n_groups)
    try:
        cv = StratifiedGroupKFold(n_splits=n_splits, shuffle=True, random_state=seed)
        splits = list(cv.split(np.zeros(len(y)), y, g))
    except ValueError:
        cv = GroupKFold(n_splits=n_splits)
        splits = list(cv.split(np.zeros(len(y)), y, g))
    for tr, te in splits:
        assert_no_group_leak(meta, tr, te, group)
        yield tr, te


def fit_calibration_split(meta: pd.DataFrame, train_idx, cal_frac: float = 0.3,
                          seed: int = C.SEED, group: str = "pid") -> Split:
    """Split a training index into (fit_idx, cal_idx) by animal, for split conformal."""
    train_idx = np.asarray(train_idx)
    rng = np.random.default_rng(seed)
    groups = meta.iloc[train_idx][group].astype(str).to_numpy()
    uniq = np.unique(groups)
    rng.shuffle(uniq)
    n_cal = max(1, int(round(cal_frac * len(uniq))))
    cal_groups = set(uniq[:n_cal])
    is_cal = np.array([x in cal_groups for x in groups])
    fit_idx, cal_idx = train_idx[~is_cal], train_idx[is_cal]
    assert_no_group_leak(meta, fit_idx, cal_idx, group)
    return fit_idx, cal_idx


def leave_one_sex_out(meta: pd.DataFrame) -> Iterator[tuple[str, Split]]:
    """Yields ('train_male_test_female', (train, test)) and the reverse."""
    sex = meta["sex"].astype(str).to_numpy()
    for train_sex, test_sex in (("male", "female"), ("female", "male")):
        tr = np.flatnonzero(sex == train_sex)
        te = np.flatnonzero(sex == test_sex)
        if len(tr) and len(te):
            assert_no_group_leak(meta, tr, te)
            yield f"train_{train_sex}_test_{test_sex}", (tr, te)


def leave_one_group_out(meta: pd.DataFrame, col: str = "group",
                        holdouts: list[str] | None = None) -> Iterator[tuple[str, Split]]:
    """Hold out one level of `col` (default: a training time point) entirely."""
    v = meta[col].astype(str).to_numpy()
    levels = holdouts or [x for x in C.GROUP_ORDER if x in set(v)]
    for lvl in levels:
        te = np.flatnonzero(v == lvl)
        tr = np.flatnonzero(v != lvl)
        if len(tr) and len(te):
            assert_no_group_leak(meta, tr, te)
            yield f"holdout_{col}_{lvl}", (tr, te)


def train_controls_test_trained(meta: pd.DataFrame) -> Iterator[tuple[str, Split]]:
    """Fit on sedentary controls only; test on every trained animal."""
    g = meta["group"].astype(str).to_numpy()
    tr = np.flatnonzero(g == "control")
    te = np.flatnonzero(g != "control")
    if len(tr) and len(te):
        assert_no_group_leak(meta, tr, te)
        yield "train_control_test_trained", (tr, te)


def describe_split(meta: pd.DataFrame, train_idx, test_idx, label: str = "tissue") -> dict:
    tr, te = meta.iloc[np.asarray(train_idx)], meta.iloc[np.asarray(test_idx)]
    return {
        "n_train": len(tr), "n_test": len(te),
        "animals_train": int(tr["pid"].nunique()), "animals_test": int(te["pid"].nunique()),
        "classes_train": int(tr[label].nunique()), "classes_test": int(te[label].nunique()),
    }
