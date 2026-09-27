"""Leakage tests: every split keeps all samples of an animal on one side."""
import numpy as np
import pandas as pd
import pytest

from tfp import splits


def _meta(n_animals=30, tissues=("A", "B", "C", "D"), seed=0):
    rng = np.random.default_rng(seed)
    rows = []
    for i in range(n_animals):
        sex = "male" if i % 2 else "female"
        group = ["control", "1w", "2w", "4w", "8w"][i % 5]
        for t in tissues:
            rows.append({"viallabel": f"{i:05d}{rng.integers(100000, 999999)}", "pid": f"P{i:03d}", "bid": f"{i:05d}",
                         "sex": sex, "group": group, "tissue": t, "assay": "X"})
    return pd.DataFrame(rows).set_index("viallabel")


def test_grouped_kfold_no_leak_and_full_coverage():
    meta = _meta()
    seen = np.zeros(len(meta), dtype=int)
    for tr, te in splits.grouped_kfold(meta, "tissue", n_splits=5):
        splits.assert_no_group_leak(meta, tr, te)
        assert len(set(tr) & set(te)) == 0
        seen[te] += 1
    assert (seen == 1).all(), "every sample is tested exactly once"


def test_assert_detects_leak():
    meta = _meta()
    with pytest.raises(AssertionError):
        splits.assert_no_group_leak(meta, np.arange(0, 10), np.arange(5, 15))


def test_calibration_split_by_animal():
    meta = _meta()
    tr = np.arange(len(meta))
    fit, cal = splits.fit_calibration_split(meta, tr, cal_frac=0.3)
    splits.assert_no_group_leak(meta, fit, cal)
    assert len(fit) + len(cal) == len(tr)
    n_cal_animals = meta.iloc[cal]["pid"].nunique()
    assert 0.2 <= n_cal_animals / meta["pid"].nunique() <= 0.4


def test_shift_splits_are_clean():
    meta = _meta()
    for name, (tr, te) in list(splits.leave_one_sex_out(meta)) + list(splits.leave_one_group_out(meta)) \
            + list(splits.train_controls_test_trained(meta)):
        splits.assert_no_group_leak(meta, tr, te)
        assert len(tr) and len(te), name
    for name, (tr, te) in splits.leave_one_sex_out(meta):
        assert meta.iloc[tr]["sex"].nunique() == 1 and meta.iloc[te]["sex"].nunique() == 1
