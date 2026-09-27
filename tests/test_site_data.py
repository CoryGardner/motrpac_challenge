"""Provenance of the site data: every number in site/data/*.json traces to a results/ file (or is marked
pending with a reason), the home-page tiles have provenance entries, files stay small, and the per-sample
exports have the expected shapes. The value checks against results/ are skipped when results/ is absent."""
import json
import math
from pathlib import Path

import pandas as pd
import pytest

from tfp import config as C

ROOT = Path(__file__).resolve().parents[1]
SITE = ROOT / "site" / "data"
RES = C.results_root()          # results/ when a complete run is present, else the committed results_frozen/


def _res(file: str) -> Path:
    """A canonical provenance path (results/...) resolved against the results root in use; paths under
    results_multiomic/ or results_frozen/ (the multiomic follow-up page) are repository-relative."""
    if file.startswith("results_multiomic/") or file.startswith("results_frozen/"):
        return ROOT / file
    return RES / Path(file).relative_to("results")

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
    ids = [t["id"] for t in h["tiles"]]
    assert h["home_tiles"] == ["tile_acc_k20", "tile_bodymap_k20", "tile_bodymap_recal_k20", "tile_bridge"]
    assert set(h["home_tiles"]) <= set(ids) and "tile_bodymap_cov_k20" in ids   # the coverage tile stays for the Identifiability page
    for t in h["tiles"]:
        assert t["id"] in prov, f"tile {t['id']} has no provenance entry"
        assert t["value"] == prov[t["id"]]["value"]


def test_recalibration_tile_matches_the_table():
    """The home tile's 0.943 / 1.00 tissue per set are the k20, n_recal = 3 row of results/12_bodymap/recalibration.csv."""
    import pandas as pd
    h = _load("headline.json")
    t = next(t for t in h["tiles"] if t["id"] == "tile_bodymap_recal_k20")
    r = pd.read_csv(_res("results/12_bodymap/recalibration.csv"))
    r = r[(r["model"] == "k20") & (r["n_recal"] == 3)].iloc[0]
    assert abs(t["value"] - r["coverage_recalibrated"]) < 1e-9
    assert t["line"] == f"mean of {int(r['draws'])} draws · {r['set_size_recalibrated']:.2f} tissue per set"
    assert f"{int(r['n_test_individuals'])} test animals" in t["sub"]
    prov = {p["id"]: p for p in _load("provenance.json")["entries"]}
    assert prov["tile_bodymap_recal_size_k20"]["value"] == r["set_size_recalibrated"]


def test_ladder_rungs_are_complete_with_n_and_uncertainty():
    """Every rung × model × variant is present (no pending rows), carries n, and carries a fold sd or a bootstrap interval."""
    h = _load("headline.json")
    assert h["ladder"], "ladder is empty"
    pending = [r for r in h["ladder"] if r.get("pending")]
    assert not pending, pending
    for r in h["ladder"]:
        assert r["n_samples"] > 0 and r["n_individuals"] > 0, r
        assert 0 <= r["accuracy"] <= 1 and 0 <= r["coverage"] <= 1, r
        if r["variant"] == "marginal":
            assert (r.get("accuracy_sd") is not None) or (r.get("accuracy_ci") is not None), r
            assert (r.get("coverage_sd") is not None) or (r.get("coverage_ci") is not None), r
            assert r.get("wrong_non_empty") is not None and r["coverage"] + r["empty"] + r["wrong_non_empty"] <= 1 + 1e-9, r
    rungs = {(r["rung_id"], r["model"], r["variant"]) for r in h["ladder"]}
    assert {("train_male_test_female", "k50", "marginal"), ("different_species", "k50", "floored"),
            ("train_control_test_trained", "k20", "marginal"), ("train_control_test_trained", "k50", "floored"), ("train_control_test_trained", "full", "mondrian")} <= rungs
    # the training-state rung sits right after in-distribution: it is the mildest shift (same study, same laboratory)
    assert h["rung_order"].index("train_control_test_trained") == 1, h["rung_order"]
    # classes the source model knows: 18 on the held-out sex (one sex-specific tissue absent), all 19 on the training-state split
    seen = {r["rung_id"]: r["n_classes_seen"] for r in h["ladder"] if r["rung_id"].startswith("train_")}
    assert seen == {"train_male_test_female": 18, "train_female_test_male": 18, "train_control_test_trained": 19}, seen
    tr = next(r for r in h["ladder"] if r["rung_id"] == "train_control_test_trained" and r["model"] == "k20" and r["variant"] == "marginal")
    assert tr["n_individuals"] == 40 and tr["n_source_animals"] == 10 and tr["n_calibration_animals"] == 3, tr
    assert tr.get("unseen") is None and tr["n_samples_coverage"] == tr["n_samples"], tr


def test_qc_only_baseline_uses_balanced_accuracy():
    """The QC-only baseline is set beside the gene model's balanced accuracy, so it must be balanced accuracy too."""
    prov = _load("provenance.json")
    byid = {e["id"]: e for e in prov["entries"]}
    for k in ("qc_technical", "qc_composition", "qc_all"):
        assert byid[k]["column"] == "bal_acc_mean", (k, byid[k]["column"])
        assert byid[k + "_sd"]["column"] == "bal_acc_sd", (k, byid[k + "_sd"]["column"])


def test_bodymap_age_curve_carries_n_and_interval():
    """Every point of the BodyMap age curve carries its n (samples, animals) and a bootstrap interval that brackets the estimate."""
    b = _load("bodymap.json")
    assert len(b["age_accuracy"]) == 4
    for r in b["age_accuracy"]:
        assert r["n_samples"] > 0 and r["n_animals"] > 0, r
        for m in ("k20", "k50", "full"):
            lo, hi = r[f"{m}_ci"]
            assert 0 <= lo <= r[m] <= hi <= 1, (r["stage_weeks"], m, lo, r[m], hi)


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
        p = _res(e["file"])
        assert p.exists(), f"{e['id']}: {e['file']} missing under {RES}"
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
        p = _res(t["file"])
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


def test_every_extras_key_a_page_reads_is_exported():
    """A page that reads headline.extras.<key> must find it: a missing key would render as 'pending'."""
    import re
    ex = _load("headline.json")["extras"]
    missing = []
    for js in (ROOT / "site" / "assets" / "pages").glob("*.js"):
        src = js.read_text()
        for key in set(re.findall(r"\bex\.([A-Za-z0-9_]+)", src)) | set(re.findall(r"extras\[[\"']([^\"']+)[\"']\]", src)):
            if key not in ex:
                missing.append(f"{js.name}: {key}")
    assert not missing, missing


def test_manifest_lists_the_site_data_files():
    m = _load("manifest.json")
    names = {f["name"] for f in m.get("site_data_files", [])}
    assert "headline.json" in names and "provenance.json" in names and len(names) >= 20


def test_held_out_sex_rungs_carry_both_denominators():
    """accuracy_all counts the unseen-class vials (never right); coverage is over seen-class vials: both n's are shown."""
    h = _load("headline.json")
    rows = [r for r in h["ladder"] if r["rung_id"] in ("train_male_test_female", "train_female_test_male") and not r.get("pending")]
    assert rows
    for r in rows:
        assert r["n_samples_coverage"] < r["n_samples"], r
        assert 0 <= r["accuracy_seen"] <= 1, r


def test_provenance_records_the_reconciliation_counts():
    meta = _load("provenance.json")["_meta"]
    rec = meta.get("reconciliation")
    assert rec and rec["comparable"] > 0 and 0 <= rec["changed"] <= rec["comparable"]


def test_manifest_lists_only_phases_the_site_reads():
    """The manifest describes the results the export read, not every directory under results/."""
    m = _load("manifest.json")
    allowed = {"03_eda", "04_baselines", "05_panels", "06_conformal", "07_fusion", "08_shift", "12_bodymap", "13_gtex",
               "14_transfer", "14_transfer_cpm", "15_time_course", "16_identifiability", "34_panel_model"}
    assert set(m["phases"]) <= allowed, sorted(set(m["phases"]) - allowed)
    assert "09_discordance" not in m["phases"] and "absent_phases" not in m
    assert m["phases_used"] and set(m["phases_used"]) >= {"05_panels", "06_conformal", "31_site_regen"}
    assert m["results_dir"] in ("results", "results_frozen")
    for k in ("motrpac", "bodymap", "gtex"):
        assert m["data_access"][k]["date"], k
    assert all("regenerated" in f for f in m["site_data_files"])


def test_no_float_counts_in_tile_text():
    """Counts in the tile subtitles are integers ("899 vials", never "899.0 vials")."""
    import re
    h = _load("headline.json")
    for t in h["tiles"]:
        assert not re.search(r"\d+\.0 (vials|animals|samples|organs)", t.get("sub", "")), t["sub"]
    for r in h["ladder"]:
        assert isinstance(r["n_samples"], int), r["n_samples"]


def test_headline_carries_the_design_constants():
    d = _load("headline.json")["design"]
    assert d["n_outer_folds"] == 5 and d["n_train_animals"] + d["n_test_animals"] == 50
    assert d["n_fit_animals"] == 18 and d["n_cal_animals"] == 22 and d["alpha"] == 0.1 and d["k_panel"] == 20
    assert d["variance_prefilter"] == 5000 and d["C_grid"] == [0.01, 0.1, 1.0] and d["inner_splits"] == 3


def test_no_orphan_site_data_files():
    """Every JSON in site/data is read by some page script (or is the manifest, provenance or fixtures)."""
    import re
    keep = {"manifest.json", "provenance.json", "conformal_fixtures.json"}
    refs = set()
    for js in (ROOT / "site" / "assets").rglob("*.js"):
        refs |= set(re.findall(r'"data/([\w.-]+\.json)"', js.read_text()))
    orphans = sorted(p.name for p in SITE.glob("*.json") if p.name not in keep and p.name not in refs)
    assert not orphans, orphans


# ---- the Exercise page (exercise.json) --------------------------------------------------------------------
def test_exercise_json_shape():
    x = _load("exercise.json")
    for k in ("separability", "covariates", "within_tissue", "fingerprint_by_duration", "summary"):
        assert k in x, k
    rows = [r for r in x["fingerprint_by_duration"]["rows"] if not r.get("pending")]
    assert len(rows) == 12, len(rows)                      # 3 models × 4 training durations
    for r in rows:
        assert r["n_individuals"] == 10 and r["n_samples"] in (179, 180), r
        assert r["accuracy_ci"][0] <= r["accuracy"] <= r["accuracy_ci"][1], r
        assert r["coverage_ci"][0] <= r["coverage"] <= r["coverage_ci"][1], r
        assert r["coverage"] + r["empty"] + r["wrong_non_empty"] <= 1 + 1e-9, r
    assert len(x["separability"]["duration"]) == 28 and len(x["covariates"]["auroc"]) == 84
    assert len(x["within_tissue"]["TRNSCRPT"]) == 19 and len(x["within_tissue"]["PROT"]) == 7


def test_exercise_summary_has_provenance_and_no_pending():
    """Every scalar the Exercise page reads carries a provenance entry with the same value; pending is allowed only
    for the parts whose result files are absent."""
    x = _load("exercise.json")
    prov = {e["id"]: e for e in _load("provenance.json")["entries"]}
    optional = ("fbd_ref_", "fbd_acc_excl_", "fbd_cov_excl_", "design_", "physio_", "ptr_")
    for k, v in x["summary"].items():
        assert k in prov, f"summary.{k} has no provenance entry"
        if prov[k].get("pending"):
            assert k.startswith(optional), f"{k} is pending"
            continue
        assert v == prov[k]["value"], (k, v, prov[k]["value"])
    have15 = (RES / "15_time_course" / "fingerprint_by_duration" / "part5_by_duration.csv").exists()
    have05 = (RES / "05_panels" / "TRNSCRPT" / "panel_training_summary.csv").exists()
    for k, e in prov.items():
        if e.get("pending") and k.startswith(("fbd_ref_", "design_", "physio_")):
            assert not have15, f"{k} pending although phase 15 results are present"
        if e.get("pending") and k.startswith("ptr_"):
            assert not have05, f"{k} pending although the panel training response is present"


def test_exercise_page_keys_are_exported():
    import re
    x = _load("exercise.json")
    keys = set()
    ex_js = ROOT / "site" / "assets" / "pages" / "exercise.js"
    if ex_js.exists():                       # the page aliases X.summary as S
        keys |= set(re.findall(r"\bS\.([A-Za-z0-9_]+)", ex_js.read_text()))
    for js in ((ROOT / "site" / "assets" / "pages" / "home.js"),):
        if js.exists():
            keys |= set(re.findall(r"\bX\.summary\.([A-Za-z0-9_]+)", js.read_text()))
    missing = sorted(k for k in keys if k not in x["summary"])
    assert not missing, missing


def test_fbd_rows_reproduce_the_shift_table():
    """The per-duration rows recomputed from the per-vial regeneration pool back to the published accuracy."""
    x = _load("exercise.json")
    rows = [r for r in x["fingerprint_by_duration"]["rows"] if r["model"] == "k20" and not r.get("pending")]
    pooled = sum(r["accuracy"] * r["n_samples"] for r in rows) / sum(r["n_samples"] for r in rows)
    st = pd.read_csv(_res("results/08_shift/TRNSCRPT/shift_table.csv"))
    ref = float(st[(st["split"] == "train_control_test_trained") & (st["arm"] == "panel_k20")]["accuracy_all"].iloc[0])
    assert abs(pooled - ref) < 1e-9, (pooled, ref)


def test_exercise_covariate_scalars_use_the_logistic_rows():
    prov = {e["id"]: e for e in _load("provenance.json")["entries"]}
    for k, e in prov.items():
        if k.startswith("cov_") and k.endswith(("_logreg", "_null95", "_p")):
            assert e["where"].get("model") == "logreg", (k, e["where"])


# ---- the multiomic follow-up page (branch multiomic-overnight) --------------------------------------------------------
MO = SITE / "multiomic.json"


@pytest.mark.skipif(not MO.exists(), reason="site/data/multiomic.json absent")
def test_multiomic_numbers_match_their_csvs():
    """Every mo_* provenance entry reproduces from the results_multiomic/ or results_frozen/ CSV it names, and the
    values the page reads from multiomic.json equal the provenance values."""
    prov = _load("provenance.json")
    entries = [e for e in prov["entries"] if str(e.get("id", "")).startswith("mo_")]
    assert len(entries) >= 150, f"only {len(entries)} multiomic provenance entries"
    checked = 0
    for e in entries:
        p = _res(e["file"])
        assert p.exists(), f"{e['id']}: {e['file']} missing"
        df = pd.read_csv(p)
        sel = df
        for k, v in (e.get("where") or {}).items():
            sel = sel[sel[k].astype(str) == str(v)]
        assert len(sel) == 1, f"{e['id']}: selector {e.get('where')} matched {len(sel)} rows in {e['file']}"
        got = sel[e["column"]].iloc[0]
        if e.get("pending"):
            assert e.get("reason") and pd.isna(got), (e["id"], got)
            checked += 1
            continue
        if isinstance(e["value"], bool):
            assert bool(got) == e["value"], (e["id"], got, e["value"])
        elif isinstance(e["value"], (int, float)):
            assert math.isclose(float(got), float(e["value"]), abs_tol=e.get("tol", 1e-6)), (e["id"], got, e["value"])
        elif e["value"] is None:
            assert pd.isna(got), (e["id"], got)
        else:
            assert str(got) == str(e["value"]), (e["id"], got, e["value"])
        checked += 1
    assert checked == len(entries)
    # the JSON the page reads carries the same values as the ledger, for the numbers the page quotes in its text and tiles
    m = _load("multiomic.json")
    byid = {e["id"]: e["value"] for e in entries}
    pc1 = next(r for r in m["scales"]["rows"] if r["PC"] == "PC1")
    assert pc1["r2_rii"] == byid["mo_r2_rii_PC1"] and pc1["r2_ratio"] == byid["mo_r2_ratio_PC1"]
    jr = next(r for r in m["ladder_species"] if r["target"].startswith("protein → Jiang 2020 (cleaned") and r["model"] == "k20")
    assert jr["accuracy"] == byid["mo_jiang_relative_k20_accuracy"] and jr["coverage"] == byid["mo_jiang_relative_k20_coverage"]
    assert jr["accuracy_ci"] == [byid["mo_jiang_relative_k20_acc_lo"], byid["mo_jiang_relative_k20_acc_hi"]]
    assert m["design_jiang"]["n_pairs_estimable"] == byid["mo_jiang_pairs_est"] and m["design_jiang"]["n_pairs_total"] == byid["mo_jiang_pairs_total"]
    mw = next(x for x in m["metabolites"] if x["leg"] == "deep_mw")
    assert mw["acc_k20"] == byid["mo_deep_mw_acc_k20"] and mw["n_mapped"] == byid["mo_deep_mw_n_mapped"]
    assert m["sato"]["n_exercised_mice_union"] == byid["mo_sato_n_exercised_mice_union"]
    # every recalibrated coverage on the page carries its set size and label-space size
    for r in m["recalibration"]:
        assert r["set_size_recalibrated"] is not None and r["n_classes_label_space"] in (7, 9, 19), r
    # copied tables: row counts match the CSVs (subsets are flagged as quantised)
    for t in prov["tables"]:
        if not str(t["id"]).startswith("mo_"):
            continue
        p = _res(t["file"])
        assert p.exists(), t
        n = sum(1 for _ in open(p)) - 1
        assert n == t["n_rows"] or t.get("quantised"), (t["id"], n, t["n_rows"])


@pytest.mark.skipif(not MO.exists(), reason="site/data/multiomic.json absent")
def test_multiomic_page_is_wired():
    html = (ROOT / "site" / "multiomic.html").read_text()
    assert 'src="assets/pages/multiomic.js"' in html and "vendor/plotly-cartesian-2.35.2.min.js" in html
    js = (ROOT / "site" / "assets" / "site.js").read_text()
    assert '"multiomic.html"' in js, "the page is not in the nav"
    page_js = (ROOT / "site" / "assets" / "pages" / "multiomic.js").read_text()
    assert "data/multiomic.json" in page_js
    # the page never types a number: the only numeric literals in its script are formatting digits, chance denominators and layout
    import re
    literals = {float(x) for x in re.findall(r"(?<![\w.])(0\.\d+)(?![\w.])", page_js)}
    assert not (literals - {0.9, 0.095, 0.08, 0.1, 0.3, 0.35}), f"numeric literals in multiomic.js that are not layout constants: {sorted(literals)}"
