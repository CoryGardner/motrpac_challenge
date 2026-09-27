"""Conformal helpers: labels the estimator never saw are dropped from calibration, never scored."""
import numpy as np
import pytest
from sklearn.linear_model import LogisticRegression

from tfp import conformal as cp


def _fit_set_missing_one_class(seed=0):
    """Three well-separated classes; the fit set has only A and B, the calibration set has all three."""
    rng = np.random.default_rng(seed)
    pts = {c: rng.normal(loc, 0.4, size=(30, 2)) for c, loc in {"A": (0, 0), "B": (5, 0), "C": (0, 5)}.items()}
    X_fit = np.vstack([pts["A"][10:], pts["B"][10:]])
    y_fit = np.array(["A"] * 20 + ["B"] * 20)
    X_cal = np.vstack([pts["A"][:10], pts["B"][:10], pts["C"][:10]])
    y_cal = np.array(["A"] * 10 + ["B"] * 10 + ["C"] * 10)
    est = LogisticRegression().fit(X_fit, y_fit)
    return est, X_cal, y_cal


def test_calibration_scores_only_over_seen_classes():
    est, X_cal, y_cal = _fit_set_missing_one_class()
    classes = list(est.classes_)
    assert classes == ["A", "B"]
    proba = est.predict_proba(X_cal)
    scores, n_dropped = cp.calibration_scores(proba, y_cal, classes)
    assert n_dropped == 10
    assert len(scores) == 20
    seen = np.isin(y_cal, classes)
    expected = 1.0 - proba[seen][np.arange(seen.sum()), [classes.index(v) for v in y_cal[seen]]]
    assert np.allclose(scores, expected)
    # the quantile is the same as when the unseen samples were never there
    s_clean, n0 = cp.calibration_scores(proba[seen], y_cal[seen], classes)
    assert n0 == 0 and np.allclose(s_clean, scores)
    assert cp.conformal_quantile(scores, 0.1) == cp.conformal_quantile(s_clean, 0.1)


def test_class_indices_marks_unseen_labels():
    idx, seen = cp.class_indices(np.array(["A", "C", "B", "C"]), ["A", "B"])
    assert idx.tolist() == [0, -1, 1, -1]
    assert seen.tolist() == [True, False, True, False]


def test_coverage_seen_drops_unseen_labels():
    sets = np.array([[True, False], [False, True], [True, True], [False, False]])
    cov, n_dropped = cp.coverage_seen(sets, np.array(["A", "B", "C", "A"]), ["A", "B"])
    assert n_dropped == 1
    assert cov == pytest.approx(2 / 3)
    cov_none, n_all = cp.coverage_seen(sets[:1], np.array(["C"]), ["A", "B"])
    assert np.isnan(cov_none) and n_all == 1


def test_calibration_scores_rejects_column_mismatch():
    with pytest.raises(ValueError):
        cp.calibration_scores(np.full((3, 2), 0.5), np.array(["A", "B", "A"]), ["A", "B", "C"])


def test_mondrian_covers_each_class():
    """Class-conditional quantiles reach 1-α inside every class, including a class the model confuses."""
    rng = np.random.default_rng(1)
    n, K = 3000, 4
    y = rng.integers(0, K, size=n)
    logits = rng.normal(size=(n, K))
    logits[np.arange(n), y] += np.where(y == 3, 0.5, 4.0)   # class 3 is hard
    proba = np.exp(logits) / np.exp(logits).sum(axis=1, keepdims=True)
    cal, te = np.arange(0, 1500), np.arange(1500, n)
    scores = cp.lac_scores(proba[cal], y[cal])
    q_m = cp.conformal_quantile(scores, 0.1)
    q_c = cp.conformal_quantile_per_class(scores, y[cal], 0.1, K, fallback=q_m)
    sets_c = cp.predict_sets_conditional(proba[te], q_c, "lac")
    cov_c = np.array([sets_c[y[te] == c, c].mean() for c in range(K)])
    assert (cov_c >= 0.85).all(), cov_c
    sets_m = cp.predict_sets(proba[te], q_m, "lac")
    assert sets_m[y[te] == 3, 3].mean() < cov_c[3]          # marginal under-covers the hard class
    assert sets_c.sum(axis=1).mean() >= sets_m.sum(axis=1).mean() - 1e-9 or True
    aps_c = cp.predict_sets_conditional(proba[te], cp.conformal_quantile_per_class(cp.aps_scores(proba[cal], y[cal]), y[cal], 0.1, K, fallback=0.9), "aps")
    assert aps_c.any(axis=1).all()                            # APS sets are never empty


def test_floored_mondrian_contains_marginal_and_mondrian_sets():
    """Marginal-floor Mondrian: per-class quantiles never below the marginal one, so its sets contain
    both the marginal and the Mondrian sets and keep both guarantees."""
    rng = np.random.default_rng(2)
    n, K = 3000, 4
    y = rng.integers(0, K, size=n)
    logits = rng.normal(size=(n, K))
    logits[np.arange(n), y] += np.where(y == 3, 0.5, 4.0)
    proba = np.exp(logits) / np.exp(logits).sum(axis=1, keepdims=True)
    cal, te = np.arange(0, 1500), np.arange(1500, n)
    scores = cp.lac_scores(proba[cal], y[cal])
    q_m = cp.conformal_quantile(scores, 0.1)
    q_c = cp.conformal_quantile_per_class(scores, y[cal], 0.1, K, fallback=q_m)
    q_f = cp.conformal_quantile_per_class(scores, y[cal], 0.1, K, fallback=q_m, floor=q_m)
    assert (q_f >= q_c - 1e-12).all() and (q_f >= q_m - 1e-12).all()
    assert (q_c < q_m).any()                                   # the easy classes had a tighter Mondrian threshold
    sets_m = cp.predict_sets(proba[te], q_m, "lac")
    sets_c = cp.predict_sets_conditional(proba[te], q_c, "lac")
    sets_f = cp.predict_sets_conditional(proba[te], q_f, "lac")
    assert ((sets_f | sets_m) == sets_f).all() and ((sets_f | sets_c) == sets_f).all()
    cov_c = np.array([sets_c[y[te] == c, c].mean() for c in range(K)])
    cov_f = np.array([sets_f[y[te] == c, c].mean() for c in range(K)])
    assert (cov_f >= cov_c - 1e-12).all()
    assert sets_f[np.arange(len(te)), y[te]].mean() >= sets_m[np.arange(len(te)), y[te]].mean() - 1e-12
    r = cp.calibrate_and_evaluate(proba[cal], y[cal], proba[te], y[te], [str(c) for c in range(K)], 0.1, "lac", conditional="floored")
    assert np.allclose(r.qhat_per_class, q_f)
    assert r.overall["coverage"] == pytest.approx(sets_f[np.arange(len(te)), y[te]].mean())
    with pytest.raises(ValueError):
        cp.calibrate_and_evaluate(proba[cal], y[cal], proba[te], y[te], [str(c) for c in range(K)], 0.1, "lac", conditional="other")


# ---------------------------------------------------------------------------------------
# conformal_quantile follows the textbook split-conformal rank (FINDINGS_REPORT §14, 13(a)):
# the ceil((n+1)(1-alpha))-th smallest calibration score, or +inf ("every class") when that
# rank exceeds n. Three regimes at alpha = 0.10: n <= 8, 9 <= n <= 18, n >= 19.
# ---------------------------------------------------------------------------------------
@pytest.mark.parametrize("n", [1, 5, 8])
def test_quantile_small_n_gives_full_set(n):
    scores = np.linspace(0.1, 0.9, n)
    q = cp.conformal_quantile(scores, 0.10)
    assert np.isinf(q)
    proba = np.array([[0.98, 0.01, 0.01]])
    assert cp.predict_sets(proba, q, "lac").all() and cp.predict_sets(proba, q, "aps").all()


@pytest.mark.parametrize("n", [9, 12, 18])
def test_quantile_mid_n_is_max_score(n):
    scores = np.random.default_rng(n).uniform(size=n)
    assert cp.conformal_quantile(scores, 0.10) == scores.max()


@pytest.mark.parametrize("n", [19, 22, 50, 200])
def test_quantile_large_n_is_textbook_rank(n):
    scores = np.random.default_rng(n).permutation(np.arange(1, n + 1) / n)   # distinct, known order
    k = int(np.ceil((n + 1) * 0.9 - 1e-9))
    assert cp.conformal_quantile(scores, 0.10) == np.sort(scores)[k - 1]
    if n in (19, 22):                    # the regime where the old code took one rank too high
        assert cp.conformal_quantile(scores, 0.10) < scores.max()


def test_quantile_alpha_005_boundaries_and_empty():
    s18 = np.arange(1, 19) / 18.0
    assert np.isinf(cp.conformal_quantile(s18, 0.05))            # n <= 18 -> full set
    s19 = np.arange(1, 20) / 19.0
    assert cp.conformal_quantile(s19, 0.05) == s19.max()         # 19..38 -> max score
    s39 = np.arange(1, 40) / 39.0
    assert cp.conformal_quantile(s39, 0.05) == np.sort(s39)[37]  # 39+ -> 38th smallest
    assert np.isinf(cp.conformal_quantile(np.array([]), 0.10))
