"""The "Check samples" page data (site/data/product.json): every pv_* provenance entry reproduces from the
results_product/ or results_frozen/ file it names, copied tables match their CSVs, the validation outputs are internally
consistent, and the reference map re-derives from site/data/expr_motrpac.json and the model."""
import json
import math
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site" / "data"
PV = ROOT / "results_product" / "40_product"
pytestmark = pytest.mark.skipif(not (SITE / "product.json").exists(), reason="site/data/product.json absent")


def _load(n):
    return json.loads((SITE / n).read_text())


def test_every_pv_entry_reproduces():
    prov = _load("provenance.json")
    entries = [e for e in prov["entries"] if e["id"].startswith("pv_")]
    assert len(entries) >= 80
    for e in entries:
        df = pd.read_csv(ROOT / e["file"])
        for k, v in e["where"].items():
            df = df[df[k].astype(str) == str(v)]
        assert len(df) == 1, e["id"]
        got = df[e["column"]].iloc[0]
        if isinstance(e["value"], bool):
            assert bool(got) == e["value"], e["id"]
        elif isinstance(e["value"], (int, float)):
            assert math.isclose(float(got), float(e["value"]), abs_tol=e["tol"]), (e["id"], got, e["value"])
        elif e["value"] is None:
            assert pd.isna(got), e["id"]
        else:
            assert str(got) == str(e["value"]), e["id"]
    for t in [t for t in prov["tables"] if t["id"].startswith("pv_")]:
        n = sum(1 for _ in open(ROOT / t["file"])) - 1
        assert n == t["n_rows"], t


def test_page_values_equal_the_ledger():
    p = _load("product.json")
    byid = {e["id"]: e["value"] for e in _load("provenance.json")["entries"] if e["id"].startswith("pv_")}
    assert p["scaling_key"]["adult_21wk.reference.coverage"] == byid["pv_scaling_adult_21wk_reference_coverage"]
    assert p["flag_key"]["motrpac_heldout.0.1.false_mismatch_rate"] == byid["pv_flag_motrpac_heldout_0.1_false_mismatch_rate"]
    assert p["ood"]["Thymus"]["frac_empty"] == byid["pv_ood_Thymus_frac_empty"]
    assert len(p["flag"]) == len(pd.read_csv(PV / "flag_rates.csv")) and len(p["venacv"]) == len(pd.read_csv(PV / "venacv_cases.csv"))


def test_validation_outputs_are_consistent():
    f = pd.read_csv(PV / "flag_rates.csv")
    assert set(f["alpha"]) == {0.05, 0.1, 0.2} and (f["n_swaps"] == 1000).all()
    s = f[["consistent_rate", "false_mismatch_rate", "cant_confirm_rate"]].sum(axis=1)
    assert np.allclose(s, 1.0)
    # a larger α gives smaller sets: fewer empty-set abstentions never rise as α falls
    for _, g in f.groupby("setting"):
        g = g.sort_values("alpha")
        assert g["cant_confirm_rate"].is_monotonic_increasing
    sc = pd.read_csv(PV / "scaling.csv")
    assert set(sc["mode"]) == {"within_all", "reference", "within_organ_alone"}
    v = pd.read_csv(PV / "venacv_cases.csv", dtype={"viallabel": str})
    h = _load("headline.json")["extras"]
    assert len(v) == h["venacv_bat_calls"] and int(v["consortium_flagged"].sum()) == h["venacv_bat_calls_flagged"]


def test_reference_map_rederives():
    model, em, p = _load("panel_model.json"), _load("expr_motrpac.json"), _load("product.json")
    gi = {g: i for i, g in enumerate(em["genes"])}
    ids = [g["id"] for g in model["genes"]]
    X = np.array([em["values"][gi[g]] for g in ids], dtype=float).T
    Z = (X - np.array(model["source_z"]["mean"])) / np.array(model["source_z"]["scale"])
    C = (Z - np.array(p["pca"]["centre"])) @ np.array(p["pca"]["loadings"]).T
    assert np.abs(C[:, 0] - np.array(p["pca"]["reference"]["x"])).max() < 1e-5
    assert np.abs(C[:, 1] - np.array(p["pca"]["reference"]["y"])).max() < 1e-5
    L = np.array(p["pca"]["loadings"])
    assert np.allclose(L @ L.T, np.eye(2), atol=1e-9)


def test_recalibration_draws_reproduce_the_published_means():
    """Per-draw recalibrated coverage (replayed from results_frozen) averages to recalibration.csv; the draws that have no
    finite threshold are counted; every value the pages read has a pv_ entry."""
    d = pd.read_csv(PV / "recal_draws.csv")
    s = pd.read_csv(PV / "recal_draws_summary.csv")
    for ds, sub in (("bodymap", "12_bodymap"), ("gtex", "13_gtex")):
        rc = pd.read_csv(ROOT / "results_frozen" / sub / "recalibration.csv")
        for n in (3, 5):
            g = d[(d["dataset"] == ds) & (d["model"] == "k20") & (d["n_recal"] == n)]
            pub = rc[(rc["model"] == "k20") & (rc["n_recal"] == n)]["coverage_recalibrated"].iloc[0]
            assert len(g) == 20 and math.isclose(g["coverage"].mean(), pub, abs_tol=1e-9), (ds, n)
            r = s[(s["dataset"] == ds) & (s["n_recal"] == n)].iloc[0]
            assert r["n_infinite"] == int((~g["finite"]).sum()) and r["n_draws_below_0.90"] == int((g["coverage"] < 0.9).sum())
            assert (g.loc[~g["finite"], "avg_set_size"] == g["n_classes"].iloc[0]).all()      # an infinite threshold puts every class in the set
    p = _load("product.json")["recal_draws"]
    byid = {e["id"]: e["value"] for e in _load("provenance.json")["entries"]}
    assert p["bodymap.3"]["min_coverage"] == byid["pv_recal_bodymap_3_min_coverage"] and p["gtex.3"]["n_infinite"] == byid["pv_recal_gtex_3_n_infinite"]


def test_composition_sensitivity():
    c = pd.read_csv(PV / "composition.csv")
    s = pd.read_csv(PV / "composition_summary.csv")
    assert len(c) == 200 and c["n"].between(8, 80).all() and (c["n_organs"] >= 3).all() and (c["n_mapped"] >= 1).all()
    for _, r in s.iterrows():
        v = c[f"{r['mode']}_{r['metric']}"]
        assert math.isclose(r["p05"], v.quantile(0.05), abs_tol=1e-12) and math.isclose(r["p95"], v.quantile(0.95), abs_tol=1e-12)
    comp = _load("product.json")["composition"]
    assert comp["within.coverage.p05"] == float(s[(s["mode"] == "within") & (s["metric"] == "coverage")]["p05"].iloc[0])


def test_venacv_both_directions():
    """All held-out vena cava vials, split by the consortium's brown-fat flag, reproduce from venacv_all.csv, and the
    page's copy equals the ledger."""
    a = pd.read_csv(PV / "venacv_all.csv", dtype={"viallabel": str})
    s = pd.read_csv(PV / "venacv_summary.csv")
    h = _load("headline.json")["extras"]
    flagged = set(pd.read_csv(ROOT / "results_frozen" / "15_time_course" / "design" / "flagged_vials.csv", dtype=str)["viallabel"])
    assert len(a) == h["venacv_vials"] and int(a["consortium_flagged"].sum()) == len(set(a["viallabel"]) & flagged)
    for fl, g in a.groupby("consortium_flagged"):
        r = s[s["consortium_flagged"] == fl].iloc[0]
        assert r["n_vials"] == len(g) and r["n_called_bat"] == int((g["call"] == "BAT").sum())
        assert r["n_consistent"] + r["n_mismatch"] + r["n_cant_confirm"] == len(g)
    v = pd.read_csv(PV / "venacv_cases.csv", dtype={"viallabel": str})
    assert set(v["viallabel"]) == set(a[(a["call"] == "BAT")]["viallabel"])      # the 7 cases are the BAT calls among all vena cava vials
    p = _load("product.json")["venacv_all"]
    byid = {e["id"]: e["value"] for e in _load("provenance.json")["entries"]}
    assert p["flagged"]["n_mismatch"] == byid["pv_venacv_flagged_n_mismatch"] and p["unflagged"]["n_consistent"] == byid["pv_venacv_unflagged_n_consistent"]
