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
