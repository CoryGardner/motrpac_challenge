"""Provenance of the site data: every number in site/data/*.json traces to a results/ file (or is marked
pending with a reason), the home-page tiles have provenance entries, files stay small, and the per-sample
exports have the expected shapes. The value checks against results/ are skipped when results/ is absent."""
import json
import math
from pathlib import Path

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site" / "data"
RES = ROOT / "results"

# Without results/ (CI) there is nothing to check against; with results/ present, a missing export is a failure.
pytestmark = pytest.mark.skipif(not (RES / "06_conformal").exists() and not (SITE / "manifest.json").exists(),
                                reason="neither results/ nor site/data present")


def _load(name):
    return json.loads((SITE / name).read_text())


def test_every_json_is_small_and_carries_meta():
    files = sorted(SITE.glob("*.json"))
    assert files, "no site data"
    for f in files:
        assert f.stat().st_size < 3_000_000, f"{f.name} exceeds 3 MB"
        d = json.loads(f.read_text())
        if f.name not in ("manifest.json", "conformal_fixtures.json"):
            assert "_meta" in d and "sources" in d["_meta"], f"{f.name} lacks _meta.sources"


def test_manifest_lists_phases_and_hash():
    m = _load("manifest.json")
    assert m["git_hash"] and m["generated"]
    assert "06_conformal" in m["phases"] and "12_bodymap" in m["phases"] and "13_gtex" in m["phases"]
    for name, ph in m["phases"].items():
        for f in ph["files"]:
            assert set(f) >= {"path", "bytes", "sha256"}


def test_home_tiles_have_provenance():
    h = _load("headline.json")
    prov = {p["id"]: p for p in _load("provenance.json")["entries"]}
    assert len(h["tiles"]) == 4
    for t in h["tiles"]:
        assert t["id"] in prov, f"tile {t['id']} has no provenance entry"
        assert t["value"] == prov[t["id"]]["value"]


def test_ladder_rungs_carry_n_and_spread_or_pending():
    h = _load("headline.json")
    assert h["ladder"], "ladder is empty"
    for r in h["ladder"]:
        if r.get("pending"):
            assert r.get("reason"), r
            continue
        assert r["n_samples"] > 0 and r["n_individuals"] > 0, r
        assert 0 <= r["accuracy"] <= 1 and 0 <= r["coverage"] <= 1, r


@pytest.mark.skipif(not (RES / "06_conformal").exists(), reason="results/ absent")
def test_provenance_values_match_results_files():
    prov = _load("provenance.json")
    entries = prov["entries"]
    checked = 0
    for e in entries:
        if e.get("pending"):
            assert e.get("reason"), e
            continue
        if e.get("agg") == "recomputed":
            assert e.get("note") and e.get("files"), e
            continue
        if not e.get("file"):
            continue
        p = ROOT / e["file"]
        assert p.exists(), f"{e['id']}: {e['file']} missing"
        if p.suffix == ".csv":
            df = pd.read_csv(p)
            sel = df
            for k, v in (e.get("where") or {}).items():
                sel = sel[sel[k].astype(str) == str(v)]
            agg = e.get("agg", "value")
            if agg == "value":
                assert len(sel) == 1, f"{e['id']}: selector {e.get('where')} matched {len(sel)} rows in {e['file']}"
                got = sel[e["column"]].iloc[0]
            elif agg == "mean":
                got = sel[e["column"]].astype(float).mean()
            elif agg == "std":
                got = sel[e["column"]].astype(float).std()
            elif agg == "sum":
                got = sel[e["column"]].astype(float).sum()
            elif agg == "count":
                got = len(sel)
            else:
                raise AssertionError(f"{e['id']}: unknown agg {agg}")
            if isinstance(e["value"], (int, float)) and not isinstance(e["value"], bool):
                assert math.isclose(float(got), float(e["value"]), abs_tol=e.get("tol", 1e-6)), (e["id"], got, e["value"])
            else:
                assert str(got) == str(e["value"]), (e["id"], got, e["value"])
            checked += 1
        elif p.suffix == ".json":
            d = json.loads(p.read_text())
            for k in e["column"].split("."):
                d = d[k]
            assert d == e["value"], (e["id"], d, e["value"])
            checked += 1
    assert checked >= 20, f"only {checked} provenance entries were checkable"
    # copied tables: the JSON copy has the same number of rows as the CSV it cites (sample tables are quantised copies)
    for t in prov["tables"]:
        p = ROOT / t["file"]
        assert p.exists(), t
        n = sum(1 for _ in open(p)) - 1
        assert n == t["n_rows"] or t.get("matrix") or t.get("quantised"), (t["id"], n, t["n_rows"])


def test_sample_exports_have_expected_shapes():
    m = _load("samples_motrpac.json")
    assert len(m["classes"]) == 19 and len(m["samples"]) == 899
    assert set(m["samples"][0]["p"]) == {"full", "k20", "k50"}
    assert all(len(s["p"]["k20"]) == 19 for s in m["samples"][:5])
    assert len(m["calibration"]) == 5
    b = _load("samples_bodymap.json")
    assert len(b["samples"]) == 316 and {"Thymus", "Uterus"} <= {s["organ"] for s in b["samples"]}
    g = _load("samples_gtex.json")
    assert len(g["samples"]) == 2485
    for d in (b, g):
        for model in ("full", "k20", "k50"):
            assert len(d["calibration"][model]) > 0
            assert d["recal_thresholds"][model].keys() >= {"3", "5"}


def test_gene_export_covers_the_required_sets():
    g = _load("genes.json")
    ids = {x["id"] for x in g["genes"]}
    syms = {x["symbol"] for x in g["genes"]}
    assert 20 <= len(ids) <= 110
    assert set(g["sets"]["k20"]) <= ids and set(g["sets"]["core"]) <= ids and set(g["sets"]["k50"]) <= ids
    assert {"Pgk2", "Gnb3", "Mybph", "Hbq1b", "Fcrl5"} <= syms
    for name in ("expr_motrpac.json", "expr_bodymap.json", "expr_gtex.json"):
        e = _load(name)
        assert set(e["genes"]) <= ids
        assert len(e["values"]) == len(e["genes"]) and len(e["values"][0]) == len(e["samples"])
