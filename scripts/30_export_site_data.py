#!/usr/bin/env python
"""Phase 30 — export everything the site shows into site/data/*.json, with provenance.

Rules (docs/superpowers/specs/2026-09-26-submission-site-spec.md):
- No number on the site is typed by hand: every value comes from a results/ file (this script records
  file, row selector, column and aggregation in site/data/provenance.json) or is recomputed here from
  per-sample outputs of the pipeline (recorded as agg = "recomputed" with a note).
- A number that cannot be sourced is exported as pending with a reason.
- results/ is read only. Per-sample scores come from the --save-scores regeneration runs under
  results/31_site_regen/ (phases 06, 12, 13) and the identifiability recompute results/16_identifiability/.

Usage: python scripts/30_export_site_data.py [--skip-expr] [--check-anchors] [--reconciliation]
       [--readme-table] [--abstract]
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import platform
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from motrpac import config as C, conformal as cp, io, transfer

ROOT = C.ROOT
RES = C.RESULTS_DIR
REGEN = RES / "31_site_regen"
SITE = ROOT / "site" / "data"
WORKSPACE = ROOT.parents[1]
PREFIX_BACKUP = WORKSPACE / "backup" / "pipeline_history" / "results_pre_quantile_fix_2026-09-25"
PLOTLY_VERSION = "2.35.2"

TISSUE_NAMES = {
    "ADRNL": "adrenal gland", "BAT": "brown adipose tissue", "BLOOD": "whole blood", "COLON": "colon",
    "CORTEX": "cerebral cortex", "HEART": "heart", "HIPPOC": "hippocampus", "HYPOTH": "hypothalamus",
    "KIDNEY": "kidney", "LIVER": "liver", "LUNG": "lung", "OVARY": "ovary", "PLASMA": "plasma",
    "SKM-GN": "gastrocnemius (skeletal muscle)", "SKM-VL": "vastus lateralis (skeletal muscle)",
    "SMLINT": "small intestine", "SPLEEN": "spleen", "TESTES": "testes", "VENACV": "vena cava", "WAT-SC": "white adipose (subcutaneous)",
}
# organ systems, in the palette order fixed by the design spec (index 1..8)
ORGAN_SYSTEMS = [
    ("brain", ["CORTEX", "HIPPOC", "HYPOTH"]), ("muscle", ["SKM-GN", "SKM-VL", "HEART"]),
    ("adipose", ["WAT-SC", "BAT"]), ("gut", ["COLON", "SMLINT"]), ("gonad", ["OVARY", "TESTES"]),
    ("circulation & immune", ["BLOOD", "SPLEEN", "VENACV"]), ("visceral", ["LIVER", "KIDNEY", "LUNG"]), ("endocrine", ["ADRNL"]),
]
DEVELOPMENTAL = ["Pgk2", "Prm1", "Tnp1", "Tnp2", "Acrv1", "Mybph", "Klf1", "Alas2", "Hbq1b", "Hbb", "Fcrl5"]
MODELS = ["k20", "k50", "full"]
VARIANTS = ["marginal", "mondrian", "floored"]


# ---------------------------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------------------------
def git_hash() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True, stderr=subprocess.DEVNULL).strip()
    except Exception:
        return "n/a"


def sha256(p: Path) -> str:
    h = hashlib.sha256()
    with p.open("rb") as f:
        for chunk in iter(lambda: f.read(1 << 20), b""):
            h.update(chunk)
    return h.hexdigest()


def jsonable(x):
    """numpy / pandas scalars → JSON: NaN → None, ±inf → "inf" / "-inf"."""
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer, int)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        v = float(x)
        if math.isnan(v):
            return None
        if math.isinf(v):
            return "inf" if v > 0 else "-inf"
        return v
    if isinstance(x, (np.ndarray, list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if x is None or isinstance(x, str):
        return x
    if pd.isna(x):
        return None
    return str(x)


def q6(x):
    return None if x is None or (isinstance(x, float) and math.isnan(x)) else round(float(x), 6)


def records(df: pd.DataFrame) -> list[dict]:
    return [jsonable(r) for r in df.to_dict("records")]


def rel(p: Path) -> str:
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


class Prov:
    """Provenance ledger: one entry per headline number, one per copied table."""

    def __init__(self):
        self.entries: list[dict] = []
        self.tables: list[dict] = []
        self.sources: set[str] = set()
        self._cache: dict[str, pd.DataFrame] = {}

    def read(self, file: str | Path) -> pd.DataFrame:
        p = RES / file if not str(file).startswith("/") else Path(file)
        key = str(p)
        if key not in self._cache:
            self._cache[key] = pd.read_csv(p)
        self.sources.add(rel(p))
        return self._cache[key]

    def val(self, id: str, file: str, column: str, where: dict | None = None, agg: str = "value",
            note: str | None = None, tol: float = 1e-6):
        """Read one value from a results CSV (filter rows by `where`, aggregate with `agg`) and record it."""
        df = self.read(file)
        sel = df
        for k, v in (where or {}).items():
            sel = sel[sel[k].astype(str) == str(v)]
        if agg == "value":
            if len(sel) != 1:
                raise ValueError(f"{id}: selector {where} matched {len(sel)} rows in {file}")
            v = sel[column].iloc[0]
        elif agg == "mean":
            v = sel[column].astype(float).mean()
        elif agg == "std":
            v = sel[column].astype(float).std()
        elif agg == "sum":
            v = sel[column].astype(float).sum()
        elif agg == "count":
            v = len(sel)
        else:
            raise ValueError(agg)
        v = jsonable(v)
        self.entries.append({"id": id, "value": v, "file": rel(RES / file), "where": where or {}, "column": column,
                             "agg": agg, "tol": tol, **({"note": note} if note else {})})
        return v

    def recomputed(self, id: str, value, files: list[str], note: str):
        self.entries.append({"id": id, "value": jsonable(value), "file": None, "files": [rel(RES / f) for f in files],
                             "agg": "recomputed", "note": note})
        for f in files:
            self.sources.add(rel(RES / f))
        return jsonable(value)

    def pending(self, id: str, reason: str):
        self.entries.append({"id": id, "value": None, "pending": True, "reason": reason})
        return {"pending": True, "reason": reason}

    def table(self, id: str, file: str, json_file: str, json_path: str) -> list[dict]:
        df = self.read(file)
        self.tables.append({"id": id, "file": rel(RES / file), "json_file": json_file, "json_path": json_path, "n_rows": int(len(df))})
        return records(df)


class Writer:
    def __init__(self, prov: Prov, generated: str, ghash: str):
        self.prov, self.generated, self.ghash = prov, generated, ghash
        self.written: dict[str, int] = {}

    def write(self, name: str, obj: dict, sources: list[str]):
        obj = {"_meta": {"generated": self.generated, "git_hash": self.ghash, "sources": sorted(set(sources))}, **obj}
        p = SITE / name
        p.write_text(json.dumps(jsonable(obj), separators=(",", ":"), ensure_ascii=False))
        self.written[name] = p.stat().st_size
        if p.stat().st_size > 3_000_000:
            raise RuntimeError(f"{name} is {p.stat().st_size / 1e6:.1f} MB (limit 3 MB)")
        print(f"  wrote {name} ({p.stat().st_size / 1e3:.0f} kB)")


# ---------------------------------------------------------------------------------------------
# A1 manifest
# ---------------------------------------------------------------------------------------------
def export_manifest(w: Writer):
    phases = {}
    for d in sorted(p for p in RES.iterdir() if p.is_dir() and p.name != "31_site_regen"):
        files = [{"path": rel(f), "bytes": f.stat().st_size, "sha256": sha256(f)} for f in sorted(d.rglob("*")) if f.is_file()]
        phases[d.name] = {"files": files, "n_files": len(files), "bytes": int(sum(f["bytes"] for f in files))}
    regen = {}
    if REGEN.exists():
        for d in sorted(p for p in REGEN.iterdir() if p.is_dir() and p.name != "logs"):
            regen[d.name] = [{"path": rel(f), "bytes": f.stat().st_size, "sha256": sha256(f)} for f in sorted(d.rglob("*")) if f.is_file()]
    expected_absent = {"17_training_transfer": "17_*", "21_identifiability_audit": "21_*", "22_26_decomposition": "2[2-6]_*"}
    absent = {k: v for k, v in expected_absent.items() if not list(RES.glob(v))}
    import sklearn
    versions = {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
                "scikit-learn": sklearn.__version__, "plotly.js": PLOTLY_VERSION}
    try:
        versions["node"] = subprocess.check_output(["node", "--version"], text=True).strip()
    except Exception:
        pass
    root_files = [{"path": rel(f), "bytes": f.stat().st_size, "sha256": sha256(f)} for f in sorted(RES.iterdir()) if f.is_file()]
    m = {"generated": w.generated, "git_hash": w.ghash, "results_dir": rel(RES), "phases": phases, "root_files": root_files,
         "regeneration_runs": regen, "absent_phases": absent, "versions": versions,
         "data": {"motrpac": "MotrpacRatTraining6moData 2.0.0 (= portal release c1.0, rn6)",
                  "bodymap": "rat BodyMap GSE53960 via bodymapRat (316 samples)", "gtex": "GTEx v8 (2017-06-05) open-access expression"}}
    (SITE / "manifest.json").write_text(json.dumps(m, indent=1))
    print(f"  wrote manifest.json: {len(phases)} phases, absent {list(absent)}")
    return m


# ---------------------------------------------------------------------------------------------
# recomputation from the --save-scores runs
# ---------------------------------------------------------------------------------------------
def coverage_from_scores(regen06: Path):
    """In-distribution coverage of the k20 / k50 / full models per fold and conformal variant, from the phase-06
    regeneration (fit 18 / calibrate 22 / test 10 animals per fold), pooled-vial and one-vial-per-animal calibration,
    LAC, α ∈ ALPHAS. The full-model marginal rows must equal coverage.csv (asserted)."""
    probs = pd.read_csv(regen06 / "scores_test_probs.csv", dtype={"viallabel": str, "pid": str})
    cal = pd.read_csv(regen06 / "scores_calibration.csv", dtype={"viallabel": str, "pid": str})
    classes = json.loads((regen06 / "classes.json").read_text())
    cov_ref = pd.read_csv(regen06 / "coverage.csv")
    pcols = [f"p_{c}" for c in classes]
    rows = []
    for fold in sorted(probs["fold"].unique()):
        for model in MODELS:
            pt = probs[(probs["fold"] == fold) & (probs["model"] == model)]
            P = pt[pcols].to_numpy(dtype=float)
            y_idx = np.array([classes.index(t) for t in pt["tissue"]])
            for mode in ("pooled", "one_per_animal"):
                cc = cal[(cal["fold"] == fold) & (cal["model"] == model) & (cal["calibration"] == mode)]
                sc = cc["score_lac"].to_numpy(dtype=float)
                yc = np.array([classes.index(t) for t in cc["tissue"]])
                for alpha in C.ALPHAS:
                    q = cp.conformal_quantile(sc, alpha)
                    for variant in VARIANTS:
                        if variant == "marginal":
                            sets = cp.predict_sets(P, q, "lac")
                        else:
                            qs = cp.conformal_quantile_per_class(sc, yc, alpha, len(classes), fallback=q,
                                                                 floor=q if variant == "floored" else None)
                            sets = cp.predict_sets_conditional(P, qs, "lac")
                        covered = sets[np.arange(len(y_idx)), y_idx]
                        size = sets.sum(axis=1)
                        rows.append({"fold": int(fold), "model": model, "calibration": mode, "conformal": variant, "alpha": alpha,
                                     "coverage": float(covered.mean()), "avg_set_size": float(size.mean()),
                                     "frac_empty": float((size == 0).mean()), "frac_singleton": float((size == 1).mean()),
                                     "n": int(len(y_idx)), "n_cal": int(len(sc)), "n_test_animals": int(pt["pid"].nunique()),
                                     "qhat": float(q)})
    rec = pd.DataFrame(rows)
    # sanity: the full-model marginal rows reproduce the phase-06 coverage.csv exactly
    ref = cov_ref[(cov_ref["method"] == "lac") & (cov_ref["conformal"] == "marginal")]
    for _, r in ref.iterrows():
        got = rec[(rec["fold"] == r["fold"]) & (rec["model"] == "full") & (rec["calibration"] == r["calibration"])
                  & (rec["conformal"] == "marginal") & np.isclose(rec["alpha"], r["alpha"])]
        assert len(got) == 1 and abs(got["coverage"].iloc[0] - r["coverage"]) < 1e-9, ("recompute mismatch", dict(r))
    return rec, probs, cal, classes


# ---------------------------------------------------------------------------------------------
# aggregates
# ---------------------------------------------------------------------------------------------
def export_aggregates(w: Writer, prov: Prov, rec: pd.DataFrame | None):
    P = prov
    # ---- panel curve ------------------------------------------------------------------------
    rr = P.read("05_panels/TRNSCRPT/panel_curve.csv")
    fc = P.read("05_panels/TRNSCRPT/panel_curve_fclassif.csv")
    curve = []
    for k in sorted(rr["k"].unique()):
        a, b = rr[rr["k"] == k], fc[fc["k"] == k]
        curve.append({"k": int(k), "roundrobin_mean": a["balanced_accuracy"].mean(), "roundrobin_sd": a["balanced_accuracy"].std(),
                      "roundrobin_folds": a.sort_values("fold")["balanced_accuracy"].tolist(),
                      "fclassif_mean": b["balanced_accuracy"].mean() if len(b) else None, "fclassif_sd": b["balanced_accuracy"].std() if len(b) else None,
                      "fclassif_folds": b.sort_values("fold")["balanced_accuracy"].tolist(),
                      "n_folds": int(a["fold"].nunique()), "n_train_animals": int(a["n_train_animals"].iloc[0]), "n_test_animals": int(a["n_test_animals"].iloc[0])})
    sel = P.read("05_panels/TRNSCRPT/selected_by_fold.csv")
    selected = {str(k): {str(f): d[d["fold"] == f]["feature_ID"].tolist() for f in sorted(d["fold"].unique())}
                for k, d in sel.groupby("k")}
    sym_sel = io.map_to_gene_symbols(sorted(sel["feature_ID"].unique()))
    stab = P.read("05_panels/TRNSCRPT/stability_k20_annotated.csv").set_index("feature_ID")
    # the round-robin selector picks the best remaining gene of each class in turn, classes in sorted order
    # (models.roundrobin_order): the k-th pick belongs to class (k − 1) mod n_classes. Verified below on the annotated genes.
    classes_sorted = sorted(C.TISSUES)
    classes_sorted = [t for t in classes_sorted if t != "PLASMA"]
    marker_known = {g: stab.loc[g, "marker_tissue"] for g in sel["feature_ID"].unique() if g in stab.index}
    verify = []
    for k in (5, 10, 15, 20):
        for f in (0, 1, 2, 3, 4):
            genes_k = selected[str(k)][str(f)]
            covered = {marker_known[g] for g in genes_k if g in marker_known}
            expected = set(classes_sorted[:min(k, len(classes_sorted))])
            verify.append({"k": k, "fold": f, "n_annotated": len([g for g in genes_k if g in marker_known]),
                           "annotated_markers_within_first_k_classes": bool(covered <= expected)})
    w.write("panel_curve.json", {
        "design": "5-fold animal-grouped CV on 899 TRNSCRPT vials (50 animals); selection inside the fold; logreg_l2 (C = 0.1) on the selected genes",
        "curve": curve, "selected_by_fold": selected,
        "gene_symbols": {g: (s if isinstance(s, str) else g) for g, s in sym_sel.items()},
        "marker_tissue_known": marker_known,
        "roundrobin_class_order": classes_sorted,
        "roundrobin_note": "the selector takes the best remaining gene of each class in turn, classes in sorted order, so the k-th pick is a marker of class (k − 1) mod 19; tissues covered at k are the first min(k, 19) classes in this order",
        "roundrobin_verification": verify,
        "consistency": P.table("selection_consistency", "05_panels/TRNSCRPT/selection_consistency.csv", "panel_curve.json", "consistency"),
        "baselines": P.table("baselines_summary", "04_baselines/TRNSCRPT/summary.csv", "panel_curve.json", "baselines"),
        "baselines_per_fold": P.table("baselines_per_fold", "04_baselines/TRNSCRPT/per_fold.csv", "panel_curve.json", "baselines_per_fold"),
    }, ["results/05_panels/TRNSCRPT/panel_curve.csv", "results/05_panels/TRNSCRPT/panel_curve_fclassif.csv",
        "results/05_panels/TRNSCRPT/selected_by_fold.csv", "results/04_baselines/TRNSCRPT/summary.csv"])

    # ---- stable core and confusions ---------------------------------------------------------
    conf = {}
    for name, f in (("k20", "05_panels/TRNSCRPT/confusion_k20.csv"), ("k30", "05_panels/TRNSCRPT/confusion_k30.csv"),
                    ("full", "04_baselines/TRNSCRPT/confusion_logreg_l2.csv")):
        df = P.read(f)
        first = df.columns[0]
        mat = df.set_index(first)
        conf[name] = {"rows": list(mat.index), "cols": list(mat.columns), "counts": mat.to_numpy().tolist(),
                      "per_tissue_accuracy": {t: float(mat.loc[t, t] / mat.loc[t].sum()) for t in mat.index if t in mat.columns},
                      "source": rel(RES / f)}
        P.tables.append({"id": f"confusion_{name}", "file": rel(RES / f), "json_file": "confusion_motrpac.json", "json_path": f"confusion.{name}.counts", "n_rows": int(len(mat)),
                         "matrix": True})
    w.write("confusion_motrpac.json", {"confusion": conf,
                                       "confusable_k20": P.table("confusable_k20", "05_panels/TRNSCRPT/confusable_pairs_k20.csv", "confusion_motrpac.json", "confusable_k20"),
                                       "confusable_full": P.table("confusable_full", "04_baselines/TRNSCRPT/confusable_pairs_logreg_l2.csv", "confusion_motrpac.json", "confusable_full")},
            [c["source"] for c in conf.values()])
    w.write("stable_core.json", {
        "core": P.table("stable_core", "05_panels/TRNSCRPT/candidate_panel_annotated.csv", "stable_core.json", "core"),
        "stability_k20": P.table("stability_k20_annotated", "05_panels/TRNSCRPT/stability_k20_annotated.csv", "stable_core.json", "stability_k20"),
        "rule": "genes selected in ≥ 80 % of 50 animal-bootstrap resamples at k = 20 (round-robin selector)"},
        ["results/05_panels/TRNSCRPT/candidate_panel_annotated.csv", "results/05_panels/TRNSCRPT/stability_k20_annotated.csv"])

    # ---- certificate and in-distribution conformal ------------------------------------------
    cert = {
        "distribution": P.table("certificate_distribution", "06_conformal/TRNSCRPT/certificate_distribution.csv", "certificate.json", "distribution"),
        "alpha_delta_grid": P.table("certificate_alpha_delta_grid", "06_conformal/TRNSCRPT/certificate_alpha_delta_grid.csv", "certificate.json", "alpha_delta_grid"),
        "sizing": P.table("sizing_table", "06_conformal/TRNSCRPT/sizing_table.csv", "certificate.json", "sizing"),
        "validity": P.table("certificate_validity", "06_conformal/TRNSCRPT/certificate_validity.csv", "certificate.json", "validity"),
        "final_first_repeat": P.table("certificate_final", "06_conformal/TRNSCRPT/certificate_final.csv", "certificate.json", "final_first_repeat"),
        "readout": P.table("readout_certificates", "08_shift/readout_certificates.csv", "certificate.json", "readout"),
        "coverage_summary": P.table("coverage_marginal_vs_mondrian", "06_conformal/TRNSCRPT/coverage_marginal_vs_mondrian.csv", "certificate.json", "coverage_summary"),
        "per_tissue_alpha01": P.table("per_tissue_marginal_vs_mondrian", "06_conformal/TRNSCRPT/per_tissue_marginal_vs_mondrian_alpha0.1.csv", "certificate.json", "per_tissue_alpha01"),
        "coverage_per_fold": P.table("coverage_per_fold", "06_conformal/TRNSCRPT/coverage.csv", "certificate.json", "coverage_per_fold"),
        "design": "per outer fold: 18 fit / 22 calibration / 10 test animals; LAC and APS; pooled vials and one vial per animal; α ∈ {0.05, 0.10, 0.20}",
    }
    if rec is not None:
        cert["coverage_recomputed_by_model"] = records(rec)
        cert["coverage_recomputed_note"] = ("recomputed from results/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv with the library's conformal "
                                            "functions: the same fold design, extended to the k = 20 / 50 panels of the same fit animals")
    w.write("certificate.json", cert, [rel(RES / f) for f in ("06_conformal/TRNSCRPT/coverage.csv", "06_conformal/TRNSCRPT/coverage_marginal_vs_mondrian.csv",
                                                             "06_conformal/TRNSCRPT/certificate_distribution.csv", "06_conformal/TRNSCRPT/sizing_table.csv",
                                                             "08_shift/readout_certificates.csv")])

    # ---- shift ----------------------------------------------------------------------------------
    w.write("shift.json", {"table": P.table("shift_table", "08_shift/TRNSCRPT/shift_table.csv", "shift.json", "table"),
                           "recalibration": P.table("shift_recalibration", "08_shift/TRNSCRPT/shift_recalibration.csv", "shift.json", "recalibration"),
                           "per_class": P.table("shift_per_class", "08_shift/TRNSCRPT/shift_per_class.csv", "shift.json", "per_class"),
                           "metab_core9": P.table("shift_table_metab", "08_shift/METAB_core9/shift_table.csv", "shift.json", "metab_core9"),
                           "design": "train + calibrate on the source split (pooled vials), test on the target; α = 0.10, LAC; k = 20 panel and full model"},
            ["results/08_shift/TRNSCRPT/shift_table.csv", "results/08_shift/TRNSCRPT/shift_recalibration.csv"])

    # ---- BodyMap ------------------------------------------------------------------------------
    bm = {k: P.table(f"bodymap_{k}", f"12_bodymap/{f}", "bodymap.json", k) for k, f in (
        ("age_accuracy", "age_shift_accuracy.csv"), ("age_shift", "age_shift.csv"), ("conformal", "conformal_transfer.csv"),
        ("recalibration", "recalibration.csv"), ("coverage_by_organ", "coverage_by_organ.csv"), ("ood_sets", "ood_sets.csv"),
        ("panel_survival", "panel_survival.csv"), ("gene_check", "panel_gene_check.csv"), ("accuracy_by_organ", "accuracy_by_organ.csv"),
        ("native", "native_panel.csv"), ("gene_overlap", "gene_overlap.csv"), ("juvenile_markers", "juvenile_marker_check.csv"))}
    for m in MODELS:
        df = P.read(f"12_bodymap/confusion_{m}_adult.csv")
        mat = df.set_index(df.columns[0])
        bm[f"confusion_{m}_adult"] = {"rows": list(mat.index), "cols": list(mat.columns), "counts": mat.to_numpy().tolist()}
    bm["organ_map"] = json.loads((REGEN / "12_bodymap" / "organ_map.json").read_text()) if (REGEN / "12_bodymap" / "organ_map.json").exists() else None
    bm["design"] = ("panels fit on all 50 MoTrPAC animals for accuracy; for conformal transfer the models are refit on 35 animals and calibrated on the other 15 "
                    "(pooled vials); BodyMap genes z-scored within BodyMap; muscle and brain scored as super-classes; thymus and uterus out-of-distribution; "
                    "the animal id is assumed from the replicate index (GEO carries none)")
    w.write("bodymap.json", bm, [rel(RES / "12_bodymap" / f) for f in ("age_shift_accuracy.csv", "conformal_transfer.csv", "recalibration.csv", "coverage_by_organ.csv",
                                                                         "ood_sets.csv", "panel_gene_check.csv", "accuracy_by_organ.csv", "native_panel.csv", "juvenile_marker_check.csv")])

    # ---- GTEx -----------------------------------------------------------------------------------
    gt = {k: P.table(f"gtex_{k}", f"13_gtex/{f}", "gtex.json", k) for k, f in (
        ("accuracy_by_tissue", "accuracy_by_tissue.csv"), ("accuracy_overall", "accuracy_overall.csv"), ("conformal", "conformal_transfer.csv"),
        ("recalibration", "recalibration.csv"), ("coverage_by_tissue", "coverage_by_tissue.csv"), ("gene_check", "panel_gene_check.csv"),
        ("panel_survival", "panel_survival.csv"), ("native", "native_panel.csv"), ("gene_overlap", "gene_overlap.csv"))}
    for m in MODELS:
        df = P.read(f"13_gtex/confusion_{m}.csv")
        mat = df.set_index(df.columns[0])
        gt[f"confusion_{m}"] = {"rows": list(mat.index), "cols": list(mat.columns), "counts": mat.to_numpy().tolist()}
    gt["organ_map"] = json.loads((REGEN / "13_gtex" / "organ_map.json").read_text()) if (REGEN / "13_gtex" / "organ_map.json").exists() else None
    rt = REGEN / "13_gtex" / "recal_thresholds.csv"
    if rt.exists():
        r = pd.read_csv(rt)
        frac = {}
        for (m, n), d in r.groupby(["model", "n_recal"]):
            frac[f"{m}_n{int(n)}"] = {"frac_infinite": float(np.isinf(d["q_t"]).mean()), "draws": int(len(d)),
                                      "n_cal_scores_median": float(d["n_cal_scores"].median())}
        gt["recal_infinite_draws"] = P.recomputed("gtex_recal_infinite_draws", frac, ["31_site_regen/13_gtex/recal_thresholds.csv"],
                                                  "fraction of recalibration draws whose threshold is +∞ (too few mapped samples for a finite rank), from the per-draw thresholds")
    else:
        gt["recal_infinite_draws"] = P.pending("gtex_recal_infinite_draws", "per-draw thresholds not regenerated (phase 13 --save-scores)")
    gt["design"] = ("1:1 orthologs only; genes z-scored within species; panels re-selected in ortholog space on all 50 MoTrPAC animals for accuracy; conformal transfer with "
                    "models refit on 35 animals and rat calibration on 15; Artery - Aorta → VENACV is an imperfect mapping; every GTEx split is grouped on the donor")
    w.write("gtex.json", gt, [rel(RES / "13_gtex" / f) for f in ("accuracy_by_tissue.csv", "accuracy_overall.csv", "conformal_transfer.csv", "recalibration.csv",
                                                                   "coverage_by_tissue.csv", "panel_gene_check.csv", "native_panel.csv", "gene_overlap.csv")])

    # ---- representations (14) --------------------------------------------------------------------
    w.write("representations.json", {
        "target_summary": P.table("rep_target_summary", "14_transfer/target_summary.csv", "representations.json", "target_summary"),
        "target_accuracy": P.table("rep_target_accuracy", "14_transfer/target_accuracy.csv", "representations.json", "target_accuracy"),
        "recalibration": P.table("rep_recalibration", "14_transfer/recalibration.csv", "representations.json", "recalibration"),
        "motrpac_cv_cost": P.table("rep_cv_cost", "14_transfer/motrpac_cv_cost.csv", "representations.json", "motrpac_cv_cost"),
        "panels": P.table("rep_panels", "14_transfer/panels.csv", "representations.json", "panels"),
        "cpm_target_summary": P.table("rep_cpm_target_summary", "14_transfer_cpm/target_summary.csv", "representations.json", "cpm_target_summary"),
        "cpm_target_accuracy": P.table("rep_cpm_target_accuracy", "14_transfer_cpm/target_accuracy.csv", "representations.json", "cpm_target_accuracy"),
        "design": "representations: per-gene z-score within dataset, within-sample rank of the panel genes, top-scoring pairs; standard vs transfer-aware selector; GTEx on log2 TPM and, in the CPM run, on log2 CPM from read counts"},
        [rel(RES / f) for f in ("14_transfer/target_summary.csv", "14_transfer/recalibration.csv", "14_transfer/motrpac_cv_cost.csv", "14_transfer_cpm/target_summary.csv")])

    # ---- identifiability (16) and batch verdict (07) ---------------------------------------------
    d16 = RES / "16_identifiability"
    nesting = {}
    for f in sorted(d16.glob("nesting_*.csv")):
        assay = f.stem.replace("nesting_", "")
        nesting[assay] = P.table(f"nesting_{assay}", f"16_identifiability/{f.name}", "nesting.json", f"nesting.{assay}")
    layers = json.loads((d16 / "layers.json").read_text())
    w.write("nesting.json", {"nesting": nesting,
                             "estimable_pairs": P.table("estimable_pairs", "16_identifiability/estimable_pairs.csv", "nesting.json", "estimable_pairs"),
                             "batch_counts": P.table("batch_counts", "16_identifiability/batch_counts.csv", "nesting.json", "batch_counts"),
                             "layers": layers["layers"],
                             "definition": "a tissue pair is estimable when the two tissues share a level of every processing variable listed in variables_used, so a within-batch contrast exists",
                             "notes": (d16 / "NOTES.md").read_text()},
            [rel(d16 / "estimable_pairs.csv"), rel(d16 / "batch_counts.csv")] + [rel(f) for f in d16.glob("nesting_*.csv")])
    w.write("qc_baseline.json", {"summary": P.table("qc_only_summary", "16_identifiability/qc_only_summary.csv", "qc_baseline.json", "summary"),
                                 "per_fold": P.table("qc_only_per_fold", "16_identifiability/qc_only_per_fold.csv", "qc_baseline.json", "per_fold"),
                                 "info": layers["qc_only"],
                                 "feature_sets": {"technical": ["RIN", "r_260_280", "r_260_230", "pct_adapter_detected", "pct_trimmed", "pct_GC", "pct_dup_sequence",
                                                                "pct_umi_dup", "pct_multimapped", "median_5_3_bias", "reads_log10", "avg_input_read_length"],
                                                  "composition": ["pct_rRNA", "pct_globin", "pct_chrM", "pct_chrX", "pct_chrY", "pct_mrna", "pct_coding", "pct_utr",
                                                                  "pct_intronic", "pct_intergenic"]},
                                 "caveat": "composition fractions (mitochondrial, globin, rRNA, intronic reads, chrX/chrY) are read biologically by MoTrPAC itself; the technical set is closer to pure processing, but RIN and duplication also depend on the tissue's RNA"},
            [rel(d16 / "qc_only_summary.csv"), rel(d16 / "qc_only_per_fold.csv")])
    tc = WORKSPACE / "docs" / "TIME_COURSE_INVESTIGATION.md"
    w.write("batch_verdict.json", {
        "conclusion": P.table("batch_conclusion", "07_fusion/batch_conclusion.csv", "batch_verdict.json", "conclusion"),
        "covariate_auroc": P.table("batch_covariate_auroc", "07_fusion/batch_covariate_auroc.csv", "batch_verdict.json", "covariate_auroc"),
        "time_course_revision": {
            "status": "interpreted" if tc.exists() else "absent",
            "source": "../../docs/TIME_COURSE_INVESTIGATION.md §1 (numbers in results/15_time_course/)" if tc.exists() else None,
            "note": ("The phase-07 verdict rests on a QC covariate set that includes biologically-read composition fractions. The phase-15 "
                     "time-course investigation (2026-09-26) revises the reading per tissue from the duration profiles of RNA, protein and "
                     "metabolite separation; its verdicts are interpretations and are quoted, not measured here."),
            "verdicts": ({"SKM-GN": "training", "HEART": "training in protein; RNA 8w contrast partly unexplained",
                          "KIDNEY": "RNA and metabolomics: cohort/batch contributes; protein: training-consistent",
                          "LIVER": "training at 8w", "WAT-SC": "training in protein and metabolites; none in RNA",
                          "LUNG": "cohort/batch", "CORTEX": "no training signal in RNA"} if tc.exists() else None)},
        "part5_by_duration": (P.table("tc_part5_by_duration", "15_time_course/5_6_fingerprint/part5_by_duration.csv", "batch_verdict.json", "part5_by_duration")
                              if (RES / "15_time_course/5_6_fingerprint/part5_by_duration.csv").exists() else None)},
        [rel(RES / "07_fusion/batch_conclusion.csv"), rel(RES / "07_fusion/batch_covariate_auroc.csv")])

    # ---- fusion and discordance ---------------------------------------------------------------------
    w.write("fusion.json", {"taskA": P.table("fusion_taskA", "07_fusion/taskA_summary.csv", "fusion.json", "taskA"),
                            "taskB": P.table("fusion_taskB", "07_fusion/taskB_best_vs_null.csv", "fusion.json", "taskB"),
                            "taskB_duration": P.table("fusion_taskB_duration", "07_fusion/taskB_duration_summary.csv", "fusion.json", "taskB_duration"),
                            "n_fusion_beats_single": P.val("fusion_n_beats_single", "07_fusion/taskB_best_vs_null.csv", "fusion_beats_single_by_gt_sd", agg="sum",
                                                           note="count of tissues where a fusion arm beats the best single omic by more than its fold sd"),
                            "n_beats_null": P.val("fusion_n_beats_null", "07_fusion/taskB_best_vs_null.csv", "best_beats_null_p95", agg="sum"),
                            "n_tissues": P.val("fusion_n_tissues", "07_fusion/taskB_best_vs_null.csv", "tissue", agg="count")},
            [rel(RES / "07_fusion/taskA_summary.csv"), rel(RES / "07_fusion/taskB_best_vs_null.csv")])
    dr = P.read("09_discordance/discordance_rates.csv")
    pa = P.read("09_discordance/prediction_auroc.csv")
    pa = pa[pa["auroc"].notna()]
    auroc = pa.groupby(["target", "model"])["auroc"].agg(["mean", "std", "count"]).reset_index()
    w.write("discordance.json", {"summary": P.table("discordance_summary", "09_discordance/summary.csv", "discordance.json", "summary"),
                                 "rates": P.table("discordance_rates", "09_discordance/discordance_rates.csv", "discordance.json", "rates"),
                                 "n_sig_either": P.val("disc_n_sig_either", "09_discordance/discordance_rates.csv", "n_sig_either", agg="sum"),
                                 "n_one_layer_only": P.val("disc_n_one_layer", "09_discordance/discordance_rates.csv", "n_one_layer_only", agg="sum"),
                                 "n_sign_discordant": P.val("disc_n_sign", "09_discordance/discordance_rates.csv", "n_sign_discordant", agg="sum"),
                                 "n_concordant": P.val("disc_n_concordant", "09_discordance/discordance_rates.csv", "n_concordant", agg="sum"),
                                 "auroc": records(auroc),
                                 "auroc_with_flag_logreg": P.val("disc_auroc_with_flag", "09_discordance/prediction_auroc.csv", "auroc", where={"target": "any_one_layer", "model": "logreg"}, agg="mean"),
                                 "auroc_without_flag_logreg": P.val("disc_auroc_without_flag", "09_discordance/prediction_auroc.csv", "auroc", where={"target": "any_one_layer (no regulated flag)", "model": "logreg"}, agg="mean")},
            [rel(RES / "09_discordance/summary.csv"), rel(RES / "09_discordance/discordance_rates.csv"), rel(RES / "09_discordance/prediction_auroc.csv")])

    # ---- inventory ----------------------------------------------------------------------------------
    inv = json.loads((RES / "02_inventory/summary.json").read_text())
    w.write("inventory.json", {"summary": inv, "tissue_by_assay": P.table("tissue_by_assay", "02_inventory/tissue_by_assay.csv", "inventory.json", "tissue_by_assay"),
                               "inventory": P.table("inventory", "02_inventory/inventory.csv", "inventory.json", "inventory"),
                               "tissues": [{"code": t, "name": TISSUE_NAMES.get(t, t), "system": next((i + 1 for i, (_, ts) in enumerate(ORGAN_SYSTEMS) if t in ts), None),
                                            "system_name": next((n for n, ts in ORGAN_SYSTEMS if t in ts), None)} for t in C.TISSUES],
                               "organ_systems": [{"index": i + 1, "name": n, "tissues": ts} for i, (n, ts) in enumerate(ORGAN_SYSTEMS)]},
            [rel(RES / "02_inventory/summary.json"), rel(RES / "02_inventory/tissue_by_assay.csv")])

    # ---- EDA: variance and batch partition of the leading PCs (03) ------------------------------------
    eda = {}
    for assay in ("TRNSCRPT", "PROT", "METAB"):
        eda[f"variance_{assay}"] = P.table(f"variance_partition_{assay}", f"03_eda/variance_partition_{assay}.csv", "eda.json", f"variance_{assay}")
        eda[f"batch_{assay}"] = P.table(f"batch_partition_{assay}", f"03_eda/batch_partition_{assay}.csv", "eda.json", f"batch_{assay}")
    eda["readout_variance"] = P.table("readout_variance_pc", "03_eda/readout_variance_pc1-3.csv", "eda.json", "readout_variance")
    eda["prot_diagnostic"] = P.table("prot_diagnostic_accuracy", "04_baselines/PROT/diagnostic_accuracy.csv", "eda.json", "prot_diagnostic")
    w.write("eda.json", eda, [rel(RES / f"03_eda/variance_partition_{a}.csv") for a in ("TRNSCRPT", "PROT", "METAB")] + [rel(RES / "04_baselines/PROT/diagnostic_accuracy.csv")])

    # ---- beyond: status of the parallel phases ------------------------------------------------------
    beyond = {"training_transfer": {"phase": "17", "present": bool(list(RES.glob("17_*"))), "dirs": [rel(p) for p in RES.glob("17_*")]},
              "decomposition": {"phases": "22–26", "present": bool(list(RES.glob("2[2-6]_*"))), "dirs": [rel(p) for p in RES.glob("2[2-6]_*")]},
              "identifiability_audit": {"phase": "21", "present": bool(list(RES.glob("21_*"))), "dirs": [rel(p) for p in RES.glob("21_*")]},
              "preregistration": {"present": False, "note": "no pre-registration file exists in this copy of the workspace (searched for 'prereg' and 'pre-registration')"},
              "time_course": {"phase": "15", "present": (RES / "15_time_course").exists(), "note": "training-axis work: is the fingerprint invariant to training duration? (parts 5–6)"}}
    w.write("beyond.json", beyond, [])


# ---------------------------------------------------------------------------------------------
# headline: tiles and the transfer ladder
# ---------------------------------------------------------------------------------------------
def export_headline(w: Writer, prov: Prov, rec: pd.DataFrame | None):
    P = prov
    tiles = [
        {"id": "tile_acc_k20", "value": P.val("tile_acc_k20", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 20}, agg="mean",
                                             note="mean over 5 animal-grouped folds; round-robin selector + logreg_l2"),
         "sd": P.val("tile_acc_k20_sd", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 20}, agg="std"),
         "label": "balanced accuracy of a 20-gene panel", "sub": "19 rat tissues, 899 vials, 50 animals, 5 animal-grouped folds", "format": "3",
         "source": "results/05_panels/TRNSCRPT/panel_curve.csv (k = 20, mean ± sd over folds)"},
        {"id": "tile_bodymap_k20", "value": P.val("tile_bodymap_k20", "12_bodymap/age_shift_accuracy.csv", "k20", where={"stage_weeks": 21}),
         "label": "adult organs named correctly in another lab's rats", "sub": "rat BodyMap, 21-week adults, 68 mapped samples, 8 animals, 20-gene panel", "format": "3",
         "source": "results/12_bodymap/age_shift_accuracy.csv (stage 21, k20)"},
        {"id": "tile_bodymap_cov_k20", "value": P.val("tile_bodymap_cov_k20", "12_bodymap/conformal_transfer.csv", "coverage_mapped",
                                                      where={"stage_weeks": 21, "model": "k20", "conformal": "marginal"}),
         "empty": P.val("tile_bodymap_empty_k20", "12_bodymap/conformal_transfer.csv", "frac_empty_mapped", where={"stage_weeks": 21, "model": "k20", "conformal": "marginal"}),
         "label": "coverage of the 90 % guarantee there", "sub": "MoTrPAC-calibrated prediction sets on the same adults, α = 0.10; the shortfall is empty sets", "format": "3",
         "source": "results/12_bodymap/conformal_transfer.csv (stage 21, k20, marginal)"},
        {"id": "tile_estimable", "value": P.val("tile_estimable", "16_identifiability/estimable_pairs.csv", "n_pairs_estimable", where={"assay": "TRNSCRPT"}),
         "total": P.val("tile_estimable_total", "16_identifiability/estimable_pairs.csv", "n_pairs_total", where={"assay": "TRNSCRPT"}),
         "pairs": P.val("tile_estimable_pairs", "16_identifiability/estimable_pairs.csv", "estimable_pairs", where={"assay": "TRNSCRPT"}),
         "label": "tissue pairs whose contrast exists inside one processing batch", "sub": "RNA-seq: extraction plate, library batch and flowcell; the one pair is ovary vs testes, the sex contrast", "format": "of",
         "source": "results/16_identifiability/estimable_pairs.csv (TRNSCRPT)"},
    ]
    # in-distribution accuracy per model
    acc = {"k20": (P.val("acc_k20", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 20}, agg="mean"),
                   P.val("acc_k20_sd", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 20}, agg="std")),
           "k50": (P.val("acc_k50", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 50}, agg="mean"),
                   P.val("acc_k50_sd", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 50}, agg="std")),
           "full": (P.val("acc_full", "04_baselines/TRNSCRPT/summary.csv", "balanced_accuracy_mean", where={"model": "logreg_l2"}),
                    P.val("acc_full_sd", "04_baselines/TRNSCRPT/summary.csv", "balanced_accuracy_std", where={"model": "logreg_l2"}))}
    P.val("acc_fclassif_k20", "05_panels/TRNSCRPT/panel_curve_fclassif.csv", "balanced_accuracy", where={"k": 20}, agg="mean", note="F-test selector at the same k")
    ladder = []

    def rung(**kw):
        ladder.append(kw)

    # 1. in-distribution: coverage of the full model from phase 06; k20/k50 recomputed from the regeneration
    cov06 = P.read("06_conformal/TRNSCRPT/coverage.csv")
    for model in MODELS:
        for variant in VARIANTS:
            for calib in ("pooled", "one_per_animal"):
                if model == "full" and (calib == "pooled" or variant == "marginal"):
                    where = {"calibration": calib, "conformal": variant, "method": "lac", "alpha": 0.1}
                    cov = P.val(f"cov_id_{model}_{variant}_{calib}", "06_conformal/TRNSCRPT/coverage.csv", "coverage", where=where, agg="mean")
                    sd = P.val(f"cov_id_{model}_{variant}_{calib}_sd", "06_conformal/TRNSCRPT/coverage.csv", "coverage", where=where, agg="std")
                    empty = P.val(f"empty_id_{model}_{variant}_{calib}", "06_conformal/TRNSCRPT/coverage.csv", "frac_empty", where=where, agg="mean")
                    size = P.val(f"size_id_{model}_{variant}_{calib}", "06_conformal/TRNSCRPT/coverage.csv", "avg_set_size", where=where, agg="mean")
                    src = "results/06_conformal/TRNSCRPT/coverage.csv"
                elif rec is not None:
                    d = rec[(rec["model"] == model) & (rec["conformal"] == variant) & (rec["calibration"] == calib) & np.isclose(rec["alpha"], 0.1)]
                    cov = P.recomputed(f"cov_id_{model}_{variant}_{calib}", d["coverage"].mean(), ["31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv", "31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv"],
                                       f"mean over 5 folds of the coverage recomputed from the phase-06 design scores ({model}, {variant}, {calib})")
                    sd, empty, size = float(d["coverage"].std()), float(d["frac_empty"].mean()), float(d["avg_set_size"].mean())
                    src = "results/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv (recomputed)"
                else:
                    rung(rung_id="in_distribution", label="in-distribution (held-out animals)", model=model, variant=variant, calibration=calib,
                         **P.pending(f"cov_id_{model}_{variant}_{calib}", "phase-06 regeneration with --save-scores not available"))
                    continue
                rung(rung_id="in_distribution", label="in-distribution (held-out animals)", model=model, variant=variant, calibration=calib,
                     accuracy=acc[model][0], accuracy_sd=acc[model][1], coverage=cov, coverage_sd=sd, empty=empty, set_size=size,
                     n_samples=899, n_individuals=50, n_calibration_animals=22, source=[src, "results/05_panels/TRNSCRPT/panel_curve.csv" if model != "full" else "results/04_baselines/TRNSCRPT/summary.csv"])
    # 2. held-out sex (phase 08): k20 and full only
    st = P.read("08_shift/TRNSCRPT/shift_table.csv")
    for split, label in (("train_male_test_female", "held-out sex: trained on males, tested on females"), ("train_female_test_male", "held-out sex: trained on females, tested on males")):
        for model in MODELS:
            arm = {"k20": "panel_k20", "full": "full"}.get(model)
            for variant in VARIANTS:
                if arm is None:
                    rung(rung_id=split, label=label, model=model, variant=variant, calibration="pooled",
                         **P.pending(f"cov_{split}_{model}_{variant}", "phase 08 ran the k = 20 panel and the full model only"))
                    continue
                col = {"marginal": "coverage_target_seen", "mondrian": "coverage_target_seen_mondrian", "floored": "coverage_target_seen_floored"}[variant]
                where = {"split": split, "arm": arm}
                cov = P.val(f"cov_{split}_{model}_{variant}", "08_shift/TRNSCRPT/shift_table.csv", col, where=where)
                size_col = {"marginal": "avg_set_size_target", "mondrian": "avg_set_size_target_mondrian", "floored": "avg_set_size_target_floored"}[variant]
                rung(rung_id=split, label=label, model=model, variant=variant, calibration="pooled",
                     accuracy=P.val(f"acc_{split}_{model}", "08_shift/TRNSCRPT/shift_table.csv", "accuracy_all", where=where) if variant == "marginal" else acc_cache[(split, model)],
                     accuracy_sd=None, coverage=cov, coverage_sd=None,
                     empty=P.val(f"empty_{split}_{model}", "08_shift/TRNSCRPT/shift_table.csv", "lac_frac_empty_target", where=where) if variant == "marginal" else None,
                     set_size=P.val(f"size_{split}_{model}_{variant}", "08_shift/TRNSCRPT/shift_table.csv", size_col, where=where),
                     coverage_source=P.val(f"covsrc_{split}_{model}", "08_shift/TRNSCRPT/shift_table.csv", "coverage_source_id", where=where) if variant == "marginal" else None,
                     recal_n3=P.val(f"recal3_{split}_{model}", "08_shift/TRNSCRPT/shift_table.csv", "cov_target_recal_N3", where=where) if variant == "marginal" else None,
                     unseen=P.val(f"unseen_{split}_{model}", "08_shift/TRNSCRPT/shift_table.csv", "unseen_classes", where=where) if variant == "marginal" else None,
                     n_samples=int(st[(st["split"] == split) & (st["arm"] == arm)]["n_test"].iloc[0]),
                     n_individuals=int(st[(st["split"] == split) & (st["arm"] == arm)]["n_test_animals"].iloc[0]),
                     n_calibration_animals=int(st[(st["split"] == split) & (st["arm"] == arm)]["n_cal_animals"].iloc[0]),
                     source=["results/08_shift/TRNSCRPT/shift_table.csv"])
                if variant == "marginal":
                    acc_cache[(split, model)] = ladder[-1]["accuracy"]
    # 3. different lab: BodyMap adults
    rc12 = P.read("12_bodymap/recalibration.csv")
    n_adult_animals = int(rc12["n_recal"].iloc[0] + rc12["n_test_individuals"].iloc[0])
    for model in MODELS:
        a = P.val(f"acc_bodymap_{model}", "12_bodymap/age_shift_accuracy.csv", model, where={"stage_weeks": 21})
        for variant in VARIANTS:
            where = {"stage_weeks": 21, "model": model, "conformal": variant}
            rung(rung_id="different_lab", label="different laboratory: rat BodyMap adults (same strain)", model=model, variant=variant, calibration="pooled",
                 accuracy=a, accuracy_sd=None,
                 coverage=P.val(f"cov_bodymap_{model}_{variant}", "12_bodymap/conformal_transfer.csv", "coverage_mapped", where=where), coverage_sd=None,
                 empty=P.val(f"empty_bodymap_{model}_{variant}", "12_bodymap/conformal_transfer.csv", "frac_empty_mapped", where=where),
                 set_size=P.val(f"size_bodymap_{model}_{variant}", "12_bodymap/conformal_transfer.csv", "avg_set_size_mapped", where=where),
                 recal_n3=P.val(f"recal3_bodymap_{model}", "12_bodymap/recalibration.csv", "coverage_recalibrated", where={"model": model, "n_recal": 3}) if variant == "marginal" else None,
                 recal_n3_size=P.val(f"recal3size_bodymap_{model}", "12_bodymap/recalibration.csv", "set_size_recalibrated", where={"model": model, "n_recal": 3}) if variant == "marginal" else None,
                 n_samples=P.val(f"n_bodymap_{model}_{variant}", "12_bodymap/conformal_transfer.csv", "n_mapped", where=where), n_individuals=n_adult_animals, n_calibration_animals=15,
                 source=["results/12_bodymap/age_shift_accuracy.csv", "results/12_bodymap/conformal_transfer.csv", "results/12_bodymap/recalibration.csv"])
    # 4. different species: GTEx
    go = P.read("13_gtex/gene_overlap.csv")
    for model in MODELS:
        a = P.val(f"acc_gtex_{model}", "13_gtex/accuracy_overall.csv", "accuracy_sample_weighted", where={"model": model})
        for variant in VARIANTS:
            where = {"stage": "adult", "model": model, "conformal": variant}
            rung(rung_id="different_species", label="different species: human GTEx v8 (1:1 orthologs)", model=model, variant=variant, calibration="pooled",
                 accuracy=a, accuracy_sd=None,
                 coverage=P.val(f"cov_gtex_{model}_{variant}", "13_gtex/conformal_transfer.csv", "coverage_mapped", where=where), coverage_sd=None,
                 empty=P.val(f"empty_gtex_{model}_{variant}", "13_gtex/conformal_transfer.csv", "frac_empty_mapped", where=where),
                 set_size=P.val(f"size_gtex_{model}_{variant}", "13_gtex/conformal_transfer.csv", "avg_set_size_mapped", where=where),
                 recal_n3=P.val(f"recal3_gtex_{model}", "13_gtex/recalibration.csv", "coverage_recalibrated", where={"model": model, "n_recal": 3}) if variant == "marginal" else None,
                 recal_n3_size=P.val(f"recal3size_gtex_{model}", "13_gtex/recalibration.csv", "set_size_recalibrated", where={"model": model, "n_recal": 3}) if variant == "marginal" else None,
                 n_samples=P.val(f"n_gtex_{model}_{variant}", "13_gtex/conformal_transfer.csv", "n_mapped", where=where), n_individuals=int(go["gtex_donors"].iloc[0]), n_calibration_animals=15,
                 source=["results/13_gtex/accuracy_overall.csv", "results/13_gtex/conformal_transfer.csv", "results/13_gtex/recalibration.csv"])
    # a few more headline numbers used in prose
    extras = {
        "gtex_heart_k20": P.val("gtex_heart_k20", "13_gtex/accuracy_by_tissue.csv", "accuracy", where={"model": "k20", "gtex_tissue": "Heart - Left Ventricle"}),
        "gtex_heart_k20_top": P.val("gtex_heart_k20_top", "13_gtex/accuracy_by_tissue.csv", "top_prediction", where={"model": "k20", "gtex_tissue": "Heart - Left Ventricle"}),
        "gtex_heart_k20_top_frac": P.val("gtex_heart_k20_top_frac", "13_gtex/accuracy_by_tissue.csv", "top_prediction_frac", where={"model": "k20", "gtex_tissue": "Heart - Left Ventricle"}),
        "gtex_ovary_k20": P.val("gtex_ovary_k20", "13_gtex/accuracy_by_tissue.csv", "accuracy", where={"model": "k20", "gtex_tissue": "Ovary"}),
        "gtex_ovary_k20_top": P.val("gtex_ovary_k20_top", "13_gtex/accuracy_by_tissue.csv", "top_prediction", where={"model": "k20", "gtex_tissue": "Ovary"}),
        "gtex_ovary_k20_top_frac": P.val("gtex_ovary_k20_top_frac", "13_gtex/accuracy_by_tissue.csv", "top_prediction_frac", where={"model": "k20", "gtex_tissue": "Ovary"}),
        "gtex_native_k20": P.val("gtex_native_k20", "13_gtex/native_panel.csv", "native_gtex_panel_k20_accuracy"),
        "bodymap_native_k20": P.val("bodymap_native_k20", "12_bodymap/native_panel.csv", "native_panel_k20_accuracy_9_mapped_organs"),
        "bodymap_imported_k20_all_ages": P.val("bodymap_imported_k20_all_ages", "12_bodymap/native_panel.csv", "imported_k20_accuracy_9_mapped_organs_all_ages"),
        "bodymap_acc_2wk_k20": P.val("bodymap_acc_2wk_k20", "12_bodymap/age_shift_accuracy.csv", "k20", where={"stage_weeks": 2}),
        "bodymap_acc_6wk_k20": P.val("bodymap_acc_6wk_k20", "12_bodymap/age_shift_accuracy.csv", "k20", where={"stage_weeks": 6}),
        "bodymap_acc_104wk_k20": P.val("bodymap_acc_104wk_k20", "12_bodymap/age_shift_accuracy.csv", "k20", where={"stage_weeks": 104}),
        "bodymap_floored_k20": P.val("bodymap_floored_k20", "12_bodymap/conformal_transfer.csv", "coverage_mapped", where={"stage_weeks": 21, "model": "k20", "conformal": "floored"}),
        "shared_genes_bodymap": P.val("shared_genes_bodymap", "12_bodymap/gene_overlap.csv", "shared_ensembl"),
        "motrpac_genes": P.val("motrpac_genes", "12_bodymap/gene_overlap.csv", "motrpac_genes"),
        "orthologs_1to1": P.val("orthologs_1to1", "13_gtex/gene_overlap.csv", "rat_genes_with_1to1_human_ortholog"),
        "orthologs_in_gtex": P.val("orthologs_in_gtex", "13_gtex/gene_overlap.csv", "orthologs_present_in_gtex"),
        "gtex_samples": P.val("gtex_samples", "13_gtex/gene_overlap.csv", "gtex_samples"),
        "gtex_donors": P.val("gtex_donors", "13_gtex/gene_overlap.csv", "gtex_donors"),
        "qc_technical": P.val("qc_technical", "16_identifiability/qc_only_summary.csv", "acc_mean", where={"features": "technical"}),
        "qc_technical_sd": P.val("qc_technical_sd", "16_identifiability/qc_only_summary.csv", "acc_sd", where={"features": "technical"}),
        "qc_composition": P.val("qc_composition", "16_identifiability/qc_only_summary.csv", "acc_mean", where={"features": "composition"}),
        "qc_composition_sd": P.val("qc_composition_sd", "16_identifiability/qc_only_summary.csv", "acc_sd", where={"features": "composition"}),
        "qc_all": P.val("qc_all", "16_identifiability/qc_only_summary.csv", "acc_mean", where={"features": "all"}),
        "qc_all_sd": P.val("qc_all_sd", "16_identifiability/qc_only_summary.csv", "acc_sd", where={"features": "all"}),
        "n_plates": P.val("n_plates", "16_identifiability/batch_counts.csv", "n_levels", where={"assay": "TRNSCRPT", "variable": "RNA_extr_plate_ID"}),
        "n_lib_batches": P.val("n_lib_batches", "16_identifiability/batch_counts.csv", "n_levels", where={"assay": "TRNSCRPT", "variable": "Lib_batch_ID"}),
        "n_flowcells": P.val("n_flowcells", "16_identifiability/batch_counts.csv", "n_levels", where={"assay": "TRNSCRPT", "variable": "Seq_flowcell_ID"}),
        "cert_k20_frac": P.val("cert_k20_frac", "06_conformal/TRNSCRPT/certificate_distribution.csv", "frac_certified", where={"k": 20}),
        "cert_none_frac": 1 - P.val("cert_any_frac", "06_conformal/TRNSCRPT/certificate_alpha_delta_grid.csv", "frac_any_certified", where={"alpha": 0.1, "delta": 0.1}),
        "n_zero_error_0505": P.val("n_zero_error_0505", "06_conformal/TRNSCRPT/sizing_table.csv", "n_zero_error_needed", where={"alpha": 0.05, "delta": 0.05}),
        "n_zero_error_1005": P.val("n_zero_error_1005", "06_conformal/TRNSCRPT/sizing_table.csv", "n_zero_error_needed", where={"alpha": 0.1, "delta": 0.05}),
        "n_zero_error_1010": P.val("n_zero_error_1010", "06_conformal/TRNSCRPT/sizing_table.csv", "n_zero_error_needed", where={"alpha": 0.1, "delta": 0.1}),
        "n_zero_error_2010": P.val("n_zero_error_2010", "06_conformal/TRNSCRPT/sizing_table.csv", "n_zero_error_needed", where={"alpha": 0.2, "delta": 0.1}),
        "skmgn_cov_marginal": P.val("skmgn_cov_marginal", "06_conformal/TRNSCRPT/per_tissue_marginal_vs_mondrian_alpha0.1.csv", "coverage_marginal", where={"y_true": "SKM-GN"}),
        "skmgn_cov_mondrian": P.val("skmgn_cov_mondrian", "06_conformal/TRNSCRPT/per_tissue_marginal_vs_mondrian_alpha0.1.csv", "coverage_mondrian", where={"y_true": "SKM-GN"}),
        "pgk2_testes_2wk": P.val("pgk2_testes_2wk", "12_bodymap/juvenile_marker_check.csv", "mean_log2_cpm", where={"gene_symbol": "Pgk2", "organ": "Testes", "stage_weeks": 2}),
        "pgk2_testes_21wk": P.val("pgk2_testes_21wk", "12_bodymap/juvenile_marker_check.csv", "mean_log2_cpm", where={"gene_symbol": "Pgk2", "organ": "Testes", "stage_weeks": 21}),
    }
    # super-class fractions from the confusion matrix (the spec's "heart → skeletal muscle 0.907" is SKM-GN + SKM-VL)
    ck = P.read("13_gtex/confusion_k20.csv").set_index("organ")
    heart = ck.loc["Heart - Left Ventricle"]
    extras["gtex_heart_k20_to_skm_frac"] = P.recomputed("gtex_heart_k20_to_skm_frac", float((heart["SKM-GN"] + heart["SKM-VL"]) / heart.sum()),
                                                        ["13_gtex/confusion_k20.csv"], "fraction of GTEx heart samples called either skeletal muscle class (SKM-GN + SKM-VL) by the k20 panel")
    rt = REGEN / "13_gtex" / "recal_thresholds.csv"
    if rt.exists():
        r = pd.read_csv(rt)
        d = r[(r["model"] == "k20") & (r["n_recal"] == 3)]
        extras["gtex_recal_k20_n3_frac_inf"] = P.recomputed("gtex_recal_k20_n3_frac_inf", float(np.isinf(d["q_t"]).mean()), ["31_site_regen/13_gtex/recal_thresholds.csv"],
                                                            "fraction of the 20 three-donor recalibration draws (k20) whose threshold is +∞")
    else:
        extras["gtex_recal_k20_n3_frac_inf"] = P.pending("gtex_recal_k20_n3_frac_inf", "per-draw thresholds not regenerated")
    w.write("headline.json", {"question": "Can a molecular signature identify a tissue reliably?",
                              "tiles": tiles, "accuracy": {m: {"mean": a[0], "sd": a[1]} for m, a in acc.items()}, "ladder": ladder, "extras": extras,
                              "rung_order": ["in_distribution", "train_male_test_female", "train_female_test_male", "different_lab", "different_species"]},
            sorted(prov.sources))


acc_cache: dict = {}


# ---------------------------------------------------------------------------------------------
# per-sample exports
# ---------------------------------------------------------------------------------------------
def export_samples(w: Writer, prov: Prov, probs: pd.DataFrame, cal: pd.DataFrame, classes: list[str]):
    pcols = [f"p_{c}" for c in classes]
    by = {m: probs[probs["model"] == m].set_index("viallabel") for m in MODELS}
    samples = []
    for vial, r in by["full"].iterrows():
        samples.append({"id": vial, "pid": r["pid"], "tissue": r["tissue"], "sex": r["sex"], "group": r["group"], "fold": int(r["fold"]),
                        "p": {m: [q6(v) for v in by[m].loc[vial, pcols].to_numpy(dtype=float)] for m in MODELS}})
    calibration = {}
    for fold in sorted(cal["fold"].unique()):
        calibration[str(int(fold))] = {}
        for mode in ("pooled", "one_per_animal"):
            d = cal[(cal["fold"] == fold) & (cal["calibration"] == mode)]
            ref = d[d["model"] == "full"]
            calibration[str(int(fold))][mode] = {"vials": ref["viallabel"].tolist(), "y_idx": [classes.index(t) for t in ref["tissue"]],
                                                 "scores": {m: [q6(v) for v in d[d["model"] == m]["score_lac"].to_numpy(dtype=float)] for m in MODELS}}
            for m in MODELS:
                assert d[d["model"] == m]["viallabel"].tolist() == ref["viallabel"].tolist()
    w.write("samples_motrpac.json", {"classes": classes, "samples": samples, "calibration": calibration,
                                     "design": "phase-06 design: per fold, models fit on 18 animals, calibrated on 22 (pooled vials or one vial per animal), tested on 10 held-out animals; every vial appears once as a test vial; probabilities quantised to 6 decimals"},
            ["results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv", "results/31_site_regen/06_conformal/TRNSCRPT/scores_calibration.csv"])
    prov.tables.append({"id": "samples_motrpac", "file": "results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv", "json_file": "samples_motrpac.json", "json_path": "samples", "n_rows": len(samples), "quantised": 6})
    return samples, calibration


def export_transfer_samples(w: Writer, prov: Prov, phase: str, name: str, keep: list[str], id_col_out: str):
    d = REGEN / phase
    probs = pd.read_csv(d / "scores_target_probs.csv", dtype=str)
    cal = pd.read_csv(d / "scores_calibration.csv", dtype={"viallabel": str, "pid": str})
    classes = json.loads((d / "classes.json").read_text())
    organ_map = json.loads((d / "organ_map.json").read_text())
    rt = pd.read_csv(d / "recal_thresholds.csv")
    pcols = [f"p_{c}" for c in classes]
    by = {m: probs[probs["model"] == m].set_index("sample") for m in MODELS}
    samples = []
    for sid, r in by["full"].iterrows():
        s = {"id": sid, **{k: (int(r[k]) if k == "stage_weeks" else r[k]) for k in keep},
             "p": {m: [q6(float(v)) for v in by[m].loc[sid, pcols].to_numpy()] for m in MODELS},
             "pred_all_animals": {m: by[m].loc[sid, "pred_all_animals"] for m in MODELS}}
        samples.append(s)
    calibration = {}
    ref = cal[cal["model"] == "full"]
    for m in MODELS:
        cm = cal[cal["model"] == m]
        assert cm["viallabel"].tolist() == ref["viallabel"].tolist()
        calibration[m] = {"y_idx": [classes.index(t) for t in cm["tissue"]], "scores": [q6(v) for v in cm["score_lac"].to_numpy(dtype=float)],
                          "n_animals": int(cm["pid"].nunique()), "n_vials": int(len(cm))}
    recal = {}
    for m in MODELS:
        recal[m] = {}
        for n in sorted(rt["n_recal"].unique()):
            dd = rt[(rt["model"] == m) & (rt["n_recal"] == n)].sort_values("draw")
            first = dd.iloc[0]
            recal[m][str(int(n))] = {"draw0": {"q": jsonable(float(first["q_t"])), "chosen": str(first["chosen"]).split(";"), "n_cal_scores": int(first["n_cal_scores"])},
                                    "draws_q": [jsonable(float(v)) for v in dd["q_t"]], "frac_infinite": float(np.isinf(dd["q_t"]).mean()),
                                    "mean_coverage": None, "mean_set_size": None}
    rc = pd.read_csv(RES / phase / "recalibration.csv")
    for _, r in rc.iterrows():
        recal[r["model"]][str(int(r["n_recal"]))]["mean_coverage"] = float(r["coverage_recalibrated"])
        recal[r["model"]][str(int(r["n_recal"]))]["mean_set_size"] = float(r["set_size_recalibrated"])
    w.write(name, {"classes": classes, "organ_map": organ_map, "samples": samples, "calibration": calibration, "recal_thresholds": recal,
                   "design": ("models refit on 35 MoTrPAC animals and calibrated on the other 15 (pooled vials, LAC); probabilities quantised to 6 decimals; "
                              "pred_all_animals = the call of the all-animal model the accuracy tables use; recal_thresholds: per-draw thresholds of the "
                              "recalibration on n target individuals (draw 0 is the one the explorer offers)")},
            [f"results/31_site_regen/{phase}/scores_target_probs.csv", f"results/31_site_regen/{phase}/scores_calibration.csv", f"results/31_site_regen/{phase}/recal_thresholds.csv", f"results/{phase}/recalibration.csv"])
    prov.tables.append({"id": f"samples_{phase}", "file": f"results/31_site_regen/{phase}/scores_target_probs.csv", "json_file": name, "json_path": "samples", "n_rows": len(samples), "quantised": 6})
    return samples, calibration, classes, organ_map


# ---------------------------------------------------------------------------------------------
# genes and expression
# ---------------------------------------------------------------------------------------------
def gene_set(prov: Prov, sym: pd.Series, f2g: pd.DataFrame):
    """Union of the required gene sets, in priority order, capped at ~100."""
    P = prov
    core = P.read("05_panels/TRNSCRPT/candidate_panel_annotated.csv")["feature_ID"].tolist()
    stab = P.read("05_panels/TRNSCRPT/stability_k20_annotated.csv")
    sel = P.read("05_panels/TRNSCRPT/selected_by_fold.csv")
    k20cv = sel[sel["k"] == 20]["feature_ID"].value_counts()
    k50cv = sel[sel["k"] == 50]["feature_ID"].value_counts()
    gc12 = P.read("12_bodymap/panel_gene_check.csv")
    gc13 = P.read("13_gtex/panel_gene_check.csv")
    surv = P.read("12_bodymap/panel_survival.csv").set_index("panel")
    k20_all = surv.loc["k=20 (all-animal fit)", "genes"].split(";")
    k50_all = surv.loc["k=50 (all-animal fit)", "genes"].split(";")
    # symbol → transcript gene id (ENSRNOG only: feature_to_gene also lists protein and metabolite feature ids);
    # ids already named in the results tables win, then the stacked-matrix symbol map, then feature_to_gene
    sym_to_id = {}
    for fid, s in sym.items():                       # genes of the stacked MoTrPAC matrix
        if isinstance(s, str) and s not in sym_to_id:
            sym_to_id[s] = fid
    m = f2g.dropna(subset=["feature_ID"])
    m = m[m["feature_ID"].astype(str).str.startswith("ENSRNOG")]
    symcol = "gene_symbol" if "gene_symbol" in m.columns else "symbol"
    for fid, s in zip(m["feature_ID"], m[symcol]):
        if isinstance(s, str) and s not in sym_to_id:
            sym_to_id[s] = fid
    for fid, s in zip(gc12["feature_ID"], gc12["gene_symbol"]):
        sym_to_id[s] = fid
    for fid, s in zip(stab["feature_ID"], stab["gene_symbol"]):
        if isinstance(s, str):
            sym_to_id[s] = fid
    ids_k20_all = [sym_to_id.get(s) for s in k20_all]
    ids_k50_all = [sym_to_id.get(s) for s in k50_all]
    assert all(ids_k20_all) and len(set(ids_k20_all)) == 20, ("k20 panel symbols did not all map", k20_all)
    dev_ids = [sym_to_id.get(s) for s in DEVELOPMENTAL]
    ordered = []
    for group in (core, list(k20cv.index), gc12["feature_ID"].tolist(), gc13["feature_ID"].tolist(), [i for i in ids_k20_all if i],
                  [i for i in dev_ids if i], stab["feature_ID"].tolist(), [i for i in ids_k50_all if i], list(k50cv.index)):
        for g in group:
            if g and g not in ordered:
                ordered.append(g)
    cap = 100
    ordered = ordered[:cap]
    sets = {"core": core, "k20": [i for i in ids_k20_all if i], "k50": [i for i in ids_k50_all if i and i in ordered],
            "k20_cv_any_fold": list(k20cv.index), "developmental": [i for i in dev_ids if i], "stability_k20": stab["feature_ID"].tolist()}
    return ordered, sets, {"k20cv": k20cv, "k50cv": k50cv, "stab": stab.set_index("feature_ID"), "gc12": gc12.set_index("feature_ID"),
                           "gc13": gc13.set_index("feature_ID"), "core": set(core), "dev": {sym_to_id.get(s): s for s in DEVELOPMENTAL}}


def export_genes(w: Writer, prov: Prov, skip_expr: bool):
    P = prov
    f2g = io.load_feature_to_gene()
    pheno = io.load_pheno()
    om = io.stack_tissues("TRNSCRPT", source="counts", pheno=pheno, verbose=False)
    sym = io.map_to_gene_symbols(om.X.columns, f2g)
    genes, sets, aux = gene_set(P, sym, f2g)
    sym_all = io.map_to_gene_symbols(genes, f2g)
    rows = []
    for g in genes:
        s = aux["stab"].loc[g] if g in aux["stab"].index else None
        c12 = aux["gc12"].loc[g] if g in aux["gc12"].index else None
        c13 = aux["gc13"].loc[g] if g in aux["gc13"].index else None
        rows.append({"id": g, "symbol": sym_all.get(g) if isinstance(sym_all.get(g), str) else aux["dev"].get(g, g),
                     "marker_tissue": (s["marker_tissue"] if s is not None else (c12["marker_tissue"] if c12 is not None else (c13["marker_tissue"] if c13 is not None else None))),
                     "freq": float(s["selection_frequency"]) if s is not None else None,
                     "effect_size": float(s["effect_size"]) if s is not None else None,
                     "ovr_score": float(s["ovr_score"]) if s is not None else None,
                     "next_highest_tissue": s["next_highest_tissue"] if s is not None else None,
                     "r_pct_mrna": float(s["r_pct_mrna_in_marker_tissue"]) if s is not None else None,
                     "regulated": bool(s["risk_T7_regulated"]) if s is not None else (bool(c12["risk_T7_regulated"]) if c12 is not None and pd.notna(c12["risk_T7_regulated"]) else None),
                     "qc_flag": bool(s["risk_qc_correlated"]) if s is not None else (bool(c12["risk_qc_correlated"]) if c12 is not None and pd.notna(c12["risk_qc_correlated"]) else None),
                     "fails_bodymap": (bool(c12["fails_in_bodymap"]) if c12 is not None and pd.notna(c12["fails_in_bodymap"]) else None),
                     "weakened_bodymap": (bool(c12["weakened"]) if c12 is not None and pd.notna(c12["weakened"]) else None),
                     "bodymap_top_organ": c12["bodymap_top_organ"] if c12 is not None and isinstance(c12["bodymap_top_organ"], str) else None,
                     "fails_gtex": (bool(c13["fails_in_target"]) if c13 is not None and pd.notna(c13["fails_in_target"]) else None),
                     "weakened_gtex": (bool(c13["weakened"]) if c13 is not None and pd.notna(c13["weakened"]) else None),
                     "gtex_top_tissue": c13["target_top_organ"] if c13 is not None and isinstance(c13["target_top_organ"], str) else None,
                     "human_gene": c13["human_gene"] if c13 is not None else None,
                     "in_k20": g in sets["k20"], "in_k50": g in sets["k50"], "in_core": g in aux["core"], "in_developmental": g in sets["developmental"],
                     "n_folds_selected_k20": int(aux["k20cv"].get(g, 0)), "n_folds_selected_k50": int(aux["k50cv"].get(g, 0)),
                     "in_motrpac_matrix": g in set(om.X.columns)})
    w.write("genes.json", {"genes": rows, "sets": sets,
                           "annotation_sources": ["results/05_panels/TRNSCRPT/stability_k20_annotated.csv", "results/05_panels/TRNSCRPT/candidate_panel_annotated.csv",
                                                  "results/12_bodymap/panel_gene_check.csv", "results/13_gtex/panel_gene_check.csv", "results/12_bodymap/panel_survival.csv",
                                                  "results/05_panels/TRNSCRPT/selected_by_fold.csv"]},
            ["results/05_panels/TRNSCRPT/stability_k20_annotated.csv", "results/12_bodymap/panel_gene_check.csv", "results/13_gtex/panel_gene_check.csv", "results/12_bodymap/panel_survival.csv"])
    if skip_expr:
        return
    # MoTrPAC log2 CPM (the stacked counts matrix of the pipeline; genes absent from it are reported as absent)
    present = [g for g in genes if g in om.X.columns]
    Xm = om.X[present]
    w.write("expr_motrpac.json", {"unit": "log2 CPM (total library size), stacked-filter matrix of the pipeline", "genes": present,
                                  "symbols": [rows[genes.index(g)]["symbol"] for g in present],
                                  "samples": [{"id": v, "tissue": t, "sex": s, "group": gr} for v, t, s, gr in zip(om.X.index, om.meta["tissue"], om.meta["sex"], om.meta["group"])],
                                  "values": [[round(float(x), 2) for x in Xm[g].to_numpy()] for g in present],
                                  "absent": [g for g in genes if g not in om.X.columns]},
            ["data/raw/counts/TRNSCRPT__*.csv via io.stack_tissues"])
    # BodyMap log2 CPM
    Xb = pd.read_csv(C.EXTERNAL_DIR / "bodymap_counts.csv", index_col=0)
    Xb.index = Xb.index.astype(str)
    mb = pd.read_csv(C.EXTERNAL_DIR / "bodymap_meta.csv", dtype=str).set_index("sample").loc[Xb.columns]
    lb = io.log_cpm(Xb.T, log=True)
    present_b = [g for g in genes if g in lb.columns]
    mb["replicate_index"] = mb["replicate"].astype(str).str.rsplit("_", n=1).str[-1]
    sexmap = {"F": "female", "M": "male"}
    w.write("expr_bodymap.json", {"unit": "log2 CPM (total library size), technical runs summed", "genes": present_b,
                                  "symbols": [rows[genes.index(g)]["symbol"] for g in present_b],
                                  "samples": [{"id": s, "organ": o, "age_weeks": int(a), "sex": sexmap.get(x, x), "animal": f"{sexmap.get(x, x)[0]}_{a}_{ri}"}
                                              for s, o, a, x, ri in zip(mb.index, mb["organ"], mb["stage_weeks"], mb["sex"], mb["replicate_index"])],
                                  "values": [[round(float(x), 2) for x in lb[g].to_numpy()] for g in present_b],
                                  "absent": [g for g in genes if g not in lb.columns]},
            ["data/external/bodymap_counts.csv", "data/external/bodymap_meta.csv"])
    # GTEx log2 TPM through 1:1 orthologs
    orth = transfer.one_to_one_orthologs()
    orth = orth[orth["RAT_ENSEMBL_ID"].isin(genes)]
    header = pd.read_csv(C.EXTERNAL_DIR / "gtex_tpm_subset.csv", nrows=0)
    gsym = transfer.gtex_symbols()
    pairs, _, _ = transfer.match_gtex_orthologs(orth, header.columns[1:], gsym)
    rat_to_hum = dict(zip(pairs["RAT_ENSEMBL_ID"], pairs["HUMAN_ORTHOLOG_ENSEMBL_ID"]))
    hum_cols = [rat_to_hum[g] for g in genes if g in rat_to_hum]
    Xg = pd.read_csv(C.EXTERNAL_DIR / "gtex_tpm_subset.csv", index_col=0, usecols=[header.columns[0]] + hum_cols)
    mg = pd.read_csv(C.EXTERNAL_DIR / "gtex_meta.csv", dtype=str).set_index("SAMPID").loc[Xg.index]
    present_g = [g for g in genes if g in rat_to_hum]
    w.write("expr_gtex.json", {"unit": "log2(TPM + 1), GTEx v8, 1:1 orthologs", "genes": present_g, "human_genes": [rat_to_hum[g] for g in present_g],
                               "symbols": [rows[genes.index(g)]["symbol"] for g in present_g],
                               "samples": [{"id": s, "tissue": t, "donor": d} for s, t, d in zip(mg.index, mg["SMTSD"], mg["donor"])],
                               "values": [[round(float(x), 2) for x in Xg[rat_to_hum[g]].to_numpy()] for g in present_g],
                               "absent": [g for g in genes if g not in rat_to_hum]},
            ["data/external/gtex_tpm_subset.csv", "data/external/gtex_meta.csv", "data/raw/rat_to_human_gene.csv"])


# ---------------------------------------------------------------------------------------------
# fixtures for the JS conformal port
# ---------------------------------------------------------------------------------------------
def export_fixtures(w: Writer, motrpac: tuple, bodymap: tuple, gtex: tuple):
    samples_m, cal_m = motrpac
    samples_b, cal_b, classes, organ_map_b = bodymap
    samples_g, cal_g, _, organ_map_g = gtex
    calibrations, cases = {}, []
    n_cls = len(classes)
    for mode in ("pooled", "one_per_animal"):
        for m in MODELS:
            calibrations[f"motrpac_f0_{mode}_{m}"] = {"scores": cal_m["0"][mode]["scores"][m], "y_idx": cal_m["0"][mode]["y_idx"], "n_classes": n_cls}
    for m in MODELS:
        calibrations[f"bodymap_{m}"] = {"scores": cal_b[m]["scores"], "y_idx": cal_b[m]["y_idx"], "n_classes": n_cls}
        calibrations[f"gtex_{m}"] = {"scores": cal_g[m]["scores"], "y_idx": cal_g[m]["y_idx"], "n_classes": n_cls}
    # test-only truncated calibrations: n = 8 (+∞ at α = 0.05 and 0.10, the largest score at 0.20) and n = 18 (+∞ at α = 0.05)
    for n_small in (8, 18):
        for m in MODELS:
            src = cal_m["0"]["one_per_animal"]
            calibrations[f"motrpac_f0_opa_first{n_small}_{m}"] = {"scores": src["scores"][m][:n_small], "y_idx": src["y_idx"][:n_small], "n_classes": n_cls,
                                                                  "note": f"test-only: the first {n_small} one-vial-per-animal calibration scores of fold 0"}
    f0 = [s for s in samples_m if s["fold"] == 0]
    chosen_m, seen = [], set()
    for s in f0:
        if s["tissue"] not in seen:
            chosen_m.append(s)
            seen.add(s["tissue"])
    chosen_m = chosen_m[:19] + [s for s in f0 if s not in chosen_m][:1]
    adults = [s for s in samples_b if s["age_weeks"] == 21]
    chosen_b = [next(s for s in adults if s["organ"] == "Thymus"), next(s for s in adults if s["organ"] == "Uterus")]
    for o in ("Adrenal", "Brain", "Heart", "Kidney", "Liver", "Lung", "Muscle", "Spleen"):
        chosen_b.append(next(s for s in adults if s["organ"] == o))
    chosen_g, seen = [], set()
    for s in samples_g:
        if s["tissue"] not in seen and len(chosen_g) < 10:
            chosen_g.append(s)
            seen.add(s["tissue"])

    def add(sid, p, key, alpha, variant, true_idx):
        sc = np.array(calibrations[key]["scores"], dtype=float)
        yc = np.array(calibrations[key]["y_idx"], dtype=int)
        q = cp.conformal_quantile(sc, alpha)
        if variant == "marginal":
            sets = cp.predict_sets(np.array([p]), q, "lac")[0]
        else:
            qs = cp.conformal_quantile_per_class(sc, yc, alpha, n_cls, fallback=q, floor=q if variant == "floored" else None)
            sets = cp.predict_sets_conditional(np.array([p]), qs, "lac")[0]
        cases.append({"id": f"{sid}|{key}|{variant}|{alpha}", "sample": sid, "calibration": key, "variant": variant, "alpha": alpha,
                      "p": p, "expected_set": [int(b) for b in sets], "expected_q": jsonable(float(q)),
                      "expected_true_idx": true_idx})
    for s in chosen_m:
        for m in MODELS:
            for mode in ("pooled", "one_per_animal", "opa_first8", "opa_first18"):
                for a in C.ALPHAS:
                    for v in VARIANTS:
                        add(s["id"], s["p"][m], f"motrpac_f0_{mode}_{m}", a, v, classes.index(s["tissue"]))
    for s in chosen_b:
        for m in MODELS:
            for a in C.ALPHAS:
                for v in VARIANTS:
                    add(s["id"], s["p"][m], f"bodymap_{m}", a, v, None if not organ_map_b.get(s["organ"]) else [classes.index(t) for t in organ_map_b[s["organ"]]])
    for s in chosen_g:
        for m in MODELS:
            for a in C.ALPHAS:
                for v in VARIANTS:
                    add(s["id"], s["p"][m], f"gtex_{m}", a, v, [classes.index(t) for t in organ_map_g[s["tissue"]]])
    (SITE / "conformal_fixtures.json").write_text(json.dumps(jsonable({"classes": classes, "calibrations": calibrations, "cases": cases}), separators=(",", ":")))
    print(f"  wrote conformal_fixtures.json: {len(cases)} cases, {sum(1 for c in cases if c['expected_q'] == 'inf')} with an infinite threshold")


# ---------------------------------------------------------------------------------------------
# anchors, reconciliation, readme table, abstract
# ---------------------------------------------------------------------------------------------
ANCHORS = [  # (id in provenance, spec value, tolerance)
    ("acc_full", 0.995, 0.0005), ("acc_k20", 0.976, 0.0005), ("acc_fclassif_k20", 0.399, 0.0005), ("acc_bodymap_k20", 1.000, 1e-9),
    ("bodymap_imported_k20_all_ages", 0.9104, 0.00005), ("bodymap_acc_2wk_k20", 0.765, 0.0005), ("bodymap_acc_6wk_k20", 0.970, 0.0005),
    ("bodymap_acc_104wk_k20", 0.908, 0.0005), ("bodymap_native_k20", 0.9925, 0.00005), ("cov_bodymap_k20_marginal", 0.618, 0.0005),
    ("empty_bodymap_k20_marginal", 0.382, 0.0005), ("bodymap_floored_k20", 0.691, 0.0005), ("recal3_bodymap_k20", 0.943, 0.0005),
    ("recal3size_bodymap_k20", 0.998, 0.0005), ("acc_gtex_k20", 0.654, 0.0005), ("acc_gtex_k50", 0.781, 0.0005), ("acc_gtex_full", 0.855, 0.0005),
    ("gtex_native_k20", 0.979, 0.0005), ("gtex_heart_k20", 0.033, 0.0005), ("gtex_heart_k20_to_skm_frac", 0.907, 0.0005), ("gtex_ovary_k20", 0.0, 1e-9),
    ("gtex_ovary_k20_top_frac", 0.953, 0.0005), ("cov_gtex_k20_marginal", 0.364, 0.0005), ("empty_gtex_k20_marginal", 0.616, 0.0005),
    ("cov_gtex_full_marginal", 0.062, 0.0005), ("empty_gtex_full_marginal", 0.938, 0.0005), ("recal3_gtex_k20", 0.954, 0.0005),
    ("recal3size_gtex_k20", 11.70, 0.005), ("gtex_recal_k20_n3_frac_inf", 0.45, 0.005), ("qc_technical", 0.873, 0.0005), ("qc_technical_sd", 0.027, 0.0005),
    ("qc_composition", 0.949, 0.0005), ("qc_composition_sd", 0.011, 0.0005), ("qc_all", 0.975, 0.0005), ("qc_all_sd", 0.020, 0.0005),
    ("n_plates", 17, 0), ("n_lib_batches", 17, 0), ("n_flowcells", 4, 0), ("shared_genes_bodymap", 21040, 0), ("motrpac_genes", 21193, 0),
    ("orthologs_1to1", 14609, 0), ("orthologs_in_gtex", 14569, 0), ("disc_n_sig_either", 1948, 0), ("disc_n_one_layer", 1899, 0), ("disc_n_sign", 4, 0),
    ("disc_auroc_with_flag", 0.780, 0.0005), ("disc_auroc_without_flag", 0.742, 0.0005), ("fusion_n_beats_single", 0, 0), ("fusion_n_beats_null", 7, 0),
    ("tile_estimable", 1, 0), ("tile_estimable_total", 171, 0), ("cov_id_full_marginal_one_per_animal", 0.916, 0.0005), ("cov_id_full_marginal_pooled", 0.908, 0.0005),
]


def check_anchors(prov_entries: list[dict]) -> bool:
    byid = {e["id"]: e for e in prov_entries}
    ok_all = True
    print(f"{'anchor':42s} {'spec':>10s} {'exported':>12s}  status")
    for id_, spec, tol in ANCHORS:
        e = byid.get(id_)
        if e is None or e.get("pending") or e["value"] is None:
            print(f"{id_:42s} {spec:>10} {'—':>12}  MISSING/PENDING")
            ok_all = False
            continue
        v = e["value"]
        ok = abs(float(v) - float(spec)) <= tol + 1e-12
        ok_all &= ok
        print(f"{id_:42s} {spec:>10} {float(v):>12.4f}  {'ok' if ok else 'MISMATCH'}")
    print("all anchors ok" if ok_all else "ANCHOR MISMATCH — investigate before proceeding")
    return ok_all


def reconciliation(prov_entries: list[dict], out_md: Path):
    """docs/NUMBERS_RECONCILIATION.md: every headline number, its source, post-fix value and, where a pre-fix
    snapshot of the same file exists, the pre-fix value; plus where the stale documents still show pre-fix values."""
    rows = []
    for e in prov_entries:
        if e.get("pending") or e.get("agg") == "recomputed" or not e.get("file"):
            rows.append({"id": e["id"], "value": e.get("value"), "source": e.get("file") or ", ".join(e.get("files", [])),
                         "where": json.dumps(e.get("where", {})), "column": e.get("column", ""), "agg": e.get("agg", "pending"), "pre_fix": "", "changed": ""})
            continue
        pre = ""
        relf = e["file"].replace("results/", "", 1)
        bp = PREFIX_BACKUP / relf
        if bp.exists():
            df = pd.read_csv(bp)
            sel = df
            for k, v in e["where"].items():
                sel = sel[sel[k].astype(str) == str(v)]
            if e["agg"] == "value" and len(sel) == 1:
                pre = jsonable(sel[e["column"]].iloc[0])
            elif e["agg"] == "mean":
                pre = jsonable(sel[e["column"]].astype(float).mean())
            elif e["agg"] == "std":
                pre = jsonable(sel[e["column"]].astype(float).std())
            elif e["agg"] == "sum":
                pre = jsonable(sel[e["column"]].astype(float).sum())
            elif e["agg"] == "count":
                pre = len(sel)
        changed = ""
        if pre != "" and pre is not None and isinstance(pre, (int, float)) and isinstance(e["value"], (int, float)):
            changed = "yes" if abs(float(pre) - float(e["value"])) > 1e-9 else "no"
        elif pre != "" and pre is not None:
            changed = "yes" if str(pre) != str(e["value"]) else "no"
        rows.append({"id": e["id"], "value": e["value"], "source": e["file"], "where": json.dumps(e["where"]), "column": e["column"], "agg": e["agg"],
                     "pre_fix": pre, "changed": changed})
    df = pd.DataFrame(rows)

    def fmt(v):
        if isinstance(v, float):
            return f"{v:.4f}"
        return "" if v is None else str(v)

    lines = ["# Numbers reconciliation — what the site shows, where it comes from, and what moved with the conformal-quantile fix",
             "", f"Generated {datetime.now(timezone.utc):%Y-%m-%d %H:%M UTC} by `scripts/30_export_site_data.py --reconciliation` from `site/data/provenance.json`.",
             "`results/` (post-fix, 2026-09-25) is the truth; the pre-fix snapshot is "
             "`../../backup/pipeline_history/results_pre_quantile_fix_2026-09-25/`. Values are shown to 4 decimals; the JSON holds them unrounded.", "",
             "## 1. Headline numbers (the home-page tiles and the transfer ladder)", ""]
    head_ids = ["tile_acc_k20", "tile_acc_k20_sd", "acc_k50", "acc_full", "acc_fclassif_k20", "tile_bodymap_k20", "tile_bodymap_cov_k20", "tile_bodymap_empty_k20",
                "bodymap_floored_k20", "recal3_bodymap_k20", "recal3size_bodymap_k20", "tile_estimable", "tile_estimable_total",
                "cov_id_full_marginal_pooled", "cov_id_full_marginal_one_per_animal", "cov_id_full_mondrian_pooled", "cov_id_full_floored_pooled",
                "cov_train_male_test_female_k20_marginal", "cov_train_male_test_female_full_marginal", "cov_train_female_test_male_k20_marginal",
                "acc_gtex_k20", "acc_gtex_k50", "acc_gtex_full", "cov_gtex_k20_marginal", "empty_gtex_k20_marginal", "cov_gtex_full_marginal",
                "recal3_gtex_k20", "recal3size_gtex_k20", "gtex_recal_k20_n3_frac_inf", "qc_technical", "qc_composition", "qc_all"]
    lines += ["| id | post-fix value | pre-fix value | changed | source | selector | column | agg |", "|---|---|---|---|---|---|---|---|"]
    for i in head_ids:
        r = df[df["id"] == i]
        if len(r):
            r = r.iloc[0]
            lines.append(f"| `{i}` | {fmt(r['value'])} | {fmt(r['pre_fix'])} | {r['changed']} | `{r['source']}` | `{r['where']}` | `{r['column']}` | {r['agg']} |")
    ch = df[df["changed"] == "yes"]
    lines += ["", f"## 2. Every exported number that moved with the fix ({len(ch)} of {int((df['changed'] != '').sum())} comparable entries)", "",
              "| id | post-fix | pre-fix | source | selector | column |", "|---|---|---|---|---|---|"]
    for _, r in ch.iterrows():
        lines.append(f"| `{r['id']}` | {fmt(r['value'])} | {fmt(r['pre_fix'])} | `{r['source']}` | `{r['where']}` | `{r['column']}` |")
    # stale documents: which pre-fix headline strings still appear
    stale = {"in-distribution coverage, one vial per animal (0.962 → 0.916)": ("0.962", "0.916"),
             "Mondrian pooled α = 0.10 (0.960 → 0.919)": ("0.960", "0.919"), "floored pooled α = 0.10 (0.986 → 0.970)": ("0.986", "0.970"),
             "BodyMap k20 Mondrian (0.279 → 0.338)": ("0.279", "0.338"), "BodyMap recalibrated n = 3, k20 (0.970 → 0.943)": ("0.970", "0.943"),
             "GTEx marginal k20 (0.368 → 0.364)": ("0.368", "0.364"), "GTEx floored k20 (0.404 → 0.522)": ("0.404", "0.522"),
             "GTEx recalibrated n = 3, k20 coverage (0.880 → 0.954)": ("0.880", "0.954"), "GTEx recalibrated n = 3, k20 set size (4.15 → 11.70)": ("4.15", "11.70"),
             "male → female full coverage (0.819 → 0.807)": ("0.819", "0.807"), "METAB male → female full coverage (0.662 → 0.569)": ("0.662", "0.569")}
    docs = {"results/SUMMARY.md": RES / "SUMMARY.md", "results/ABSTRACT.md": RES / "ABSTRACT.md",
            "docs/findings/FINDINGS_REPORT.md (workspace)": WORKSPACE / "docs" / "findings" / "FINDINGS_REPORT.md"}
    lines += ["", "## 3. The three stale documents", "",
              "These documents are not edited (a banner at the top of each points here). For each headline that moved, the count of lines in each "
              "document containing the pre-fix string and the post-fix string. A pre-fix count above zero marks a passage that is stale "
              "(a string can also occur by coincidence; the counts are a locator, not a verdict).", "",
              "| number | " + " | ".join(f"{d} pre / post" for d in docs) + " |", "|---|" + "---|" * len(docs)]
    for label, (pre_s, post_s) in stale.items():
        cells = []
        for name, p in docs.items():
            if p.exists():
                txt = p.read_text(errors="ignore").splitlines()
                cells.append(f"{sum(pre_s in l for l in txt)} / {sum(post_s in l for l in txt)}")
            else:
                cells.append("absent")
        lines.append(f"| {label} | " + " | ".join(cells) + " |")
    lines += ["", "## 4. Recomputed and pending entries", "", "| id | value | files / reason |", "|---|---|---|"]
    for e in prov_entries:
        if e.get("pending"):
            lines.append(f"| `{e['id']}` | pending | {e['reason']} |")
        elif e.get("agg") == "recomputed":
            lines.append(f"| `{e['id']}` | {fmt(e['value']) if not isinstance(e['value'], dict) else 'table'} | {', '.join('`' + f + '`' for f in e['files'])}: {e['note']} |")
    out_md.write_text("\n".join(lines) + "\n")
    print(f"  wrote {out_md} ({len(ch)} changed entries)")


def readme_table(prov_entries: list[dict]) -> str:
    byid = {e["id"]: e.get("value") for e in prov_entries}

    def f(i, d=3):
        v = byid.get(i)
        return "pending" if v is None else (f"{v:.{d}f}" if isinstance(v, float) else str(v))
    rows = [
        ("20-gene panel, balanced accuracy (19 tissues, 5 animal-grouped folds)", f"{f('acc_k20')} ± {f('acc_k20_sd')}", "results/05_panels/TRNSCRPT/panel_curve.csv"),
        ("50-gene panel / all genes", f"{f('acc_k50')} / {f('acc_full')}", "results/05_panels/TRNSCRPT/panel_curve.csv, results/04_baselines/TRNSCRPT/summary.csv"),
        ("F-test selector at k = 20 (why the selector matters)", f('acc_fclassif_k20'), "results/05_panels/TRNSCRPT/panel_curve_fclassif.csv"),
        ("Coverage of 90 % sets in-distribution (pooled / one vial per animal)", f"{f('cov_id_full_marginal_pooled')} / {f('cov_id_full_marginal_one_per_animal')}", "results/06_conformal/TRNSCRPT/coverage.csv"),
        ("BodyMap adults (another lab): accuracy k20 / coverage / empty sets", f"{f('acc_bodymap_k20')} / {f('cov_bodymap_k20_marginal')} / {f('empty_bodymap_k20_marginal')}", "results/12_bodymap/"),
        ("BodyMap recalibrated on 3 animals: coverage at set size", f"{f('recal3_bodymap_k20')} at {f('recal3size_bodymap_k20', 2)}", "results/12_bodymap/recalibration.csv"),
        ("GTEx (human): accuracy k20 / k50 / full", f"{f('acc_gtex_k20')} / {f('acc_gtex_k50')} / {f('acc_gtex_full')}", "results/13_gtex/accuracy_overall.csv"),
        ("GTEx coverage k20 / empty; recalibrated on 3 donors: coverage at set size", f"{f('cov_gtex_k20_marginal')} / {f('empty_gtex_k20_marginal')}; {f('recal3_gtex_k20')} at {f('recal3size_gtex_k20', 2)}", "results/13_gtex/"),
        ("Estimable tissue pairs within study (RNA-seq)", f"{f('tile_estimable')} of {f('tile_estimable_total')}", "results/16_identifiability/estimable_pairs.csv"),
        ("QC covariates alone: technical / composition / all", f"{f('qc_technical')} / {f('qc_composition')} / {f('qc_all')}", "results/16_identifiability/qc_only_summary.csv"),
    ]
    return "\n".join(["| result | value | source |", "|---|---|---|"] + [f"| {a} | {b} | `{c}` |" for a, b, c in rows])


def abstract(prov_entries: list[dict]) -> str:
    byid = {e["id"]: e.get("value") for e in prov_entries}
    g = lambda i, d=3: ("pending" if byid.get(i) is None else (f"{byid[i]:.{d}f}" if isinstance(byid[i], float) else str(byid[i])))
    return (
        f"Can a compact molecular signature identify a tissue reliably? In the MoTrPAC 6-month rat endurance-training transcriptomes "
        f"(19 tissues, 899 vials, 50 animals) a 20-gene panel chosen by a class-aware round-robin selector reaches {g('acc_k20')} balanced "
        f"accuracy under animal-grouped cross-validation (k = 50: {g('acc_k50')}; all genes: {g('acc_full')}); the same k with an F-test "
        f"selector reaches {g('acc_fclassif_k20')}. Reliability has two further parts. A conformal prediction set that promises the true "
        f"tissue 90 % of the time covers {g('cov_id_full_marginal_pooled')} of held-out vials in-distribution, but with the same calibration "
        f"only {g('cov_bodymap_k20_marginal')} of adult rat BodyMap organs from another laboratory, where the panel still names "
        f"{g('acc_bodymap_k20')} of organs correctly, and {g('cov_gtex_k20_marginal')} of human GTEx samples (accuracy {g('acc_gtex_k20')}); "
        f"the loss is abstention (empty sets), not confident error. Three target animals restore coverage within species "
        f"({g('recal3_bodymap_k20')} at {g('recal3size_bodymap_k20', 2)} tissues per set); across species they do not "
        f"({g('recal3_gtex_k20')} at {g('recal3size_gtex_k20', 1)} tissues per set). Finally, within one multi-tissue study the tissue axis "
        f"is confounded with processing: only {g('tile_estimable')} of {g('tile_estimable_total')} tissue pairs shares an RNA extraction plate, "
        f"library batch and flowcell, and library QC numbers alone classify tissue at {g('qc_all')}. Within-study accuracy is therefore not "
        f"evidence of biology; transfer to an independently processed cohort is. Code and site: MIT."
    )


# ---------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="export site data with provenance")
    ap.add_argument("--skip-expr", action="store_true", help="skip the per-sample expression exports (slow GTEx read)")
    ap.add_argument("--check-anchors", action="store_true")
    ap.add_argument("--reconciliation", action="store_true", help="also write docs/NUMBERS_RECONCILIATION.md")
    ap.add_argument("--readme-table", action="store_true", help="print the README key-results table and exit")
    ap.add_argument("--abstract", action="store_true", help="print the abstract and exit")
    args = ap.parse_args()
    SITE.mkdir(parents=True, exist_ok=True)
    if args.readme_table or args.abstract:
        entries = json.loads((SITE / "provenance.json").read_text())["entries"]
        print(readme_table(entries) if args.readme_table else abstract(entries))
        return
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    ghash = git_hash()
    prov = Prov()
    w = Writer(prov, generated, ghash)
    print("== manifest")
    export_manifest(w)
    regen06 = REGEN / "06_conformal" / "TRNSCRPT"
    rec = probs = cal = classes = None
    if (regen06 / "scores_test_probs.csv").exists():
        print("== recompute in-distribution coverage by model from the phase-06 scores")
        rec, probs, cal, classes = coverage_from_scores(regen06)
    print("== aggregates")
    export_aggregates(w, prov, rec)
    print("== headline")
    export_headline(w, prov, rec)
    if probs is not None:
        print("== samples")
        motrpac = export_samples(w, prov, probs, cal, classes)
        bodymap = export_transfer_samples(w, prov, "12_bodymap", "samples_bodymap.json", ["organ", "stage_weeks", "sex", "animal_id"], "id")
        gtex = export_transfer_samples(w, prov, "13_gtex", "samples_gtex.json", ["organ", "donor", "stage"], "id")
        for s in bodymap[0]:
            s["age_weeks"] = s.pop("stage_weeks")
            s["animal"] = s.pop("animal_id")
        for s in gtex[0]:
            s["tissue"] = s.pop("organ")
            s.pop("stage", None)
        # rewrite with the renamed keys
        for name, data in (("samples_bodymap.json", bodymap), ("samples_gtex.json", gtex)):
            d = json.loads((SITE / name).read_text())
            d["samples"] = jsonable(data[0])
            (SITE / name).write_text(json.dumps(d, separators=(",", ":")))
        print("== fixtures")
        export_fixtures(w, motrpac, bodymap, gtex)
    else:
        print("!! phase-06 regeneration missing: per-sample exports skipped")
    print("== genes")
    export_genes(w, prov, args.skip_expr)
    (SITE / "provenance.json").write_text(json.dumps(jsonable({"_meta": {"generated": generated, "git_hash": ghash, "sources": sorted(prov.sources)},
                                                              "entries": prov.entries, "tables": prov.tables}), indent=0))
    total = sum(p.stat().st_size for p in SITE.glob("*.json"))
    print(f"== {len(prov.entries)} provenance entries, {len(prov.tables)} tables; site/data total {total / 1e6:.2f} MB; largest "
          f"{max(SITE.glob('*.json'), key=lambda p: p.stat().st_size).name}")
    if args.reconciliation:
        reconciliation(prov.entries, ROOT / "docs" / "NUMBERS_RECONCILIATION.md")
    if args.check_anchors:
        ok = check_anchors(prov.entries)
        if not ok:
            sys.exit(2)


if __name__ == "__main__":
    main()
