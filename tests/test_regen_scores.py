"""The --save-scores regeneration runs (results/31_site_regen/) must (a) reproduce the aggregates of the
published results/ phases and (b) contain per-sample scores from which those aggregates can be recomputed
with the library's conformal functions. Skipped when results/ is absent (CI without data)."""
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from tfp import config as C, conformal as cp

RES = C.results_root()          # results/ when a complete run is present, else the committed results_frozen/
REGEN = RES / "31_site_regen"

pytestmark = pytest.mark.skipif(not (RES / "06_conformal" / "TRNSCRPT" / "coverage.csv").exists(),
                                reason="results/ not present (needs the data and a full pipeline run)")


def _same_csv(a: Path, b: Path, tol: float = 1e-6):
    da, db = pd.read_csv(a), pd.read_csv(b)
    assert list(da.columns) == list(db.columns), f"columns differ: {a.name}"
    assert len(da) == len(db), f"row count differs: {a.name}"
    for c in da.columns:
        if pd.api.types.is_numeric_dtype(da[c]) and pd.api.types.is_numeric_dtype(db[c]):
            x, y = da[c].to_numpy(dtype=float), db[c].to_numpy(dtype=float)
            both_inf = np.isinf(x) & np.isinf(y)
            both_nan = np.isnan(x) & np.isnan(y)
            m = ~(both_inf | both_nan)
            assert np.allclose(x[m], y[m], atol=tol, rtol=0), f"{a.name}: column {c} differs"
        else:
            assert da[c].astype(str).tolist() == db[c].astype(str).tolist(), f"{a.name}: column {c} differs"


def test_06_regen_reproduces_published_coverage():
    d = REGEN / "06_conformal" / "TRNSCRPT"
    assert (d / "coverage.csv").exists(), "regeneration of phase 06 with --save-scores has not been run"
    # tables computed before AND after the extra --save-scores panel fit: none may move
    for f in ("coverage.csv", "per_tissue_marginal_vs_mondrian_alpha0.1.csv", "coverage_marginal_vs_mondrian.csv", "per_class_coverage.csv",
              "certificate_validity.csv", "certificate_by_fold.csv", "certificate_distribution.csv", "certificate_final.csv",
              "certificate_final_repeats.csv", "certificate_alpha_delta_grid.csv", "sizing_table.csv"):
        _same_csv(d / f, RES / "06_conformal" / "TRNSCRPT" / f)


def test_06_scores_recompute_full_model_marginal_coverage():
    d = REGEN / "06_conformal" / "TRNSCRPT"
    probs = pd.read_csv(d / "scores_test_probs.csv", dtype={"viallabel": str, "pid": str})
    cal = pd.read_csv(d / "scores_calibration.csv", dtype={"viallabel": str, "pid": str})
    classes = json.loads((d / "classes.json").read_text())
    cov = pd.read_csv(d / "coverage.csv")
    pcols = [f"p_{c}" for c in classes]
    for fold in sorted(probs["fold"].unique()):
        pt = probs[(probs["fold"] == fold) & (probs["model"] == "full")]
        P = pt[pcols].to_numpy(dtype=float)
        y_idx = np.array([classes.index(t) for t in pt["tissue"]])
        for mode in ("pooled", "one_per_animal"):
            sc = cal[(cal["fold"] == fold) & (cal["model"] == "full") & (cal["calibration"] == mode)]["score_lac"].to_numpy()
            for alpha in C.ALPHAS:
                q = cp.conformal_quantile(sc, alpha)
                sets = cp.predict_sets(P, q, "lac")
                got = float(sets[np.arange(len(y_idx)), y_idx].mean())
                row = cov[(cov["fold"] == fold) & (cov["calibration"] == mode) & (cov["conformal"] == "marginal")
                          & (cov["method"] == "lac") & (np.isclose(cov["alpha"], alpha))]
                assert len(row) == 1
                assert got == pytest.approx(float(row["coverage"].iloc[0]), abs=1e-9), (fold, mode, alpha)


@pytest.mark.parametrize("phase,files", [
    ("12_bodymap", ["conformal_transfer.csv", "recalibration.csv", "age_shift_accuracy.csv", "age_shift.csv", "coverage_by_organ.csv", "ood_sets.csv",
                    "accuracy_by_organ.csv", "panel_gene_check.csv", "panel_survival.csv", "native_panel.csv", "gene_overlap.csv",
                    "confusion_k20_adult.csv", "confusion_k50_adult.csv", "confusion_full_adult.csv"]),
    ("13_gtex", ["conformal_transfer.csv", "recalibration.csv", "accuracy_overall.csv", "accuracy_by_tissue.csv", "coverage_by_tissue.csv",
                 "panel_gene_check.csv", "panel_survival.csv", "native_panel.csv", "gene_overlap.csv", "confusion_k20.csv", "confusion_k50.csv", "confusion_full.csv"]),
])
def test_transfer_regen_reproduces_published_tables(phase, files):
    d = REGEN / phase
    assert (d / files[0]).exists(), f"regeneration of {phase} with --save-scores has not been run"
    for f in files:
        _same_csv(d / f, RES / phase / f)


@pytest.mark.parametrize("phase,stage_col,primary", [("12_bodymap", "stage_weeks", 21), ("13_gtex", "stage", "adult")])
def test_transfer_scores_recompute_marginal_coverage(phase, stage_col, primary):
    d = REGEN / phase
    probs = pd.read_csv(d / "scores_target_probs.csv", dtype=str)
    cal = pd.read_csv(d / "scores_calibration.csv", dtype={"viallabel": str, "pid": str})
    classes = json.loads((d / "classes.json").read_text())
    organ_map = json.loads((d / "organ_map.json").read_text())
    ct = pd.read_csv(d / "conformal_transfer.csv")
    pcols = [f"p_{c}" for c in classes]
    prim = probs[probs[stage_col].astype(str) == str(primary)]
    for model in ("k20", "k50", "full"):
        pm = prim[prim["model"] == model]
        P = pm[pcols].to_numpy(dtype=float)
        sc = cal[cal["model"] == model]["score_lac"].to_numpy(dtype=float)
        q = cp.conformal_quantile(sc, 0.1)
        sets = cp.predict_sets(P, q, "lac")
        mapped = np.array([bool(organ_map.get(o)) for o in pm["organ"]])
        covered = np.array([bool(organ_map.get(o)) and any(sets[i, classes.index(t)] for t in organ_map[o] if t in classes)
                            for i, o in enumerate(pm["organ"])])
        got = float(covered[mapped].mean())
        row = ct[(ct[stage_col].astype(str) == str(primary)) & (ct["model"] == model) & (ct["conformal"] == "marginal")]
        assert len(row) == 1
        assert got == pytest.approx(float(row["coverage_mapped"].iloc[0]), abs=1e-9), (phase, model)


def test_recal_thresholds_exported_for_each_model_and_n():
    for phase in ("12_bodymap", "13_gtex"):
        rt = pd.read_csv(REGEN / phase / "recal_thresholds.csv")
        assert set(rt["model"]) == {"k20", "k50", "full"}
        assert set(rt["n_recal"]) == {3, 5}
        assert (rt.groupby(["model", "n_recal"])["draw"].min() == 0).all()


def test_08_regen_reproduces_shift_table_and_per_vial_sets_agree():
    d = REGEN / "08_shift_k20"
    assert (d / "shift_table.csv").exists(), "regeneration of phase 08 (k = 20) with --save-scores has not been run"
    _same_csv(d / "shift_table.csv", RES / "08_shift" / "TRNSCRPT" / "shift_table.csv")
    _same_csv(d / "shift_recalibration.csv", RES / "08_shift" / "TRNSCRPT" / "shift_recalibration.csv")
    v = pd.read_csv(d / "scores_target_vials.csv", dtype={"viallabel": str, "pid": str})
    st = pd.read_csv(d / "shift_table.csv")
    for _, r in st.iterrows():
        seen = v[(v["split"] == r["split"]) & (v["arm"] == r["arm"]) & v["seen"]]
        assert abs(seen["covered_marginal"].mean() - r["coverage_target_seen"]) < 1e-9, (r["split"], r["arm"])
    k50 = REGEN / "08_shift_k50" / "shift_table.csv"
    assert k50.exists(), "the k = 50 shift run is missing"
    assert set(pd.read_csv(k50)["arm"]) == {"full", "panel_k50"}


def test_15_fingerprint_by_duration_matches_regen():
    """Phase 15b's primary design is the phase-08 controls→trained split scored per duration: its per-group accuracy,
    coverage and empty-set rate must equal the numbers recomputed from the regenerated per-vial sets."""
    p5_path = RES / "15_time_course" / "fingerprint_by_duration" / "part5_by_duration.csv"
    if not p5_path.exists():
        pytest.skip("phase 15b has not been run (results/15_time_course/fingerprint_by_duration/part5_by_duration.csv)")
    p5 = pd.read_csv(p5_path)
    p5 = p5[(p5["design"] == "control_only_fit7_cal3") & p5["test_group"].isin(["1w", "2w", "4w", "8w"])]
    v = pd.read_csv(REGEN / "08_shift_k20" / "scores_target_vials.csv", dtype={"viallabel": str, "pid": str})
    v = v[v["split"] == "train_control_test_trained"]
    assert set(p5["arm"]) == {"full", "panel_k20"} and len(p5) == 8
    for _, r in p5.iterrows():
        d = v[(v["arm"] == r["arm"]) & (v["group"] == r["test_group"])]
        assert len(d) == r["n_test_vials"], (r["arm"], r["test_group"])
        assert float((d["y_pred"] == d["tissue"]).mean()) == pytest.approx(r["accuracy"], abs=1e-9), (r["arm"], r["test_group"], "accuracy")
        assert float(d.loc[d["seen"], "covered_marginal"].mean()) == pytest.approx(r["coverage"], abs=1e-9), (r["arm"], r["test_group"], "coverage")
        assert float((d["size_marginal"] == 0).mean()) == pytest.approx(r["empty_rate"], abs=1e-9), (r["arm"], r["test_group"], "empty")
