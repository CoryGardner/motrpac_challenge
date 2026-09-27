#!/usr/bin/env python
"""Phase 30 — export everything the site shows into site/data/*.json, with provenance.

Rules:
- No number on the site is typed by hand: every value comes from a results/ file (this script records
  file, row selector, column and aggregation in site/data/provenance.json) or is recomputed here from
  per-sample outputs of the pipeline (recorded as agg = "recomputed" with a note).
- A number that cannot be sourced is exported as pending with a reason.
- results/ is read only. Per-sample scores come from the --save-scores regeneration runs under
  results/31_site_regen/ (phases 06, 12, 13) and the identifiability recompute results/16_identifiability/.

Usage: python scripts/30_export_site_data.py [--results DIR] [--out DIR] [--skip-expr] [--check-anchors]
       [--reconciliation] [--readme-table] [--abstract] [--write-prefix-table --prefix-backup DIR]
Results are read from --results, else $TFP_RESULTS, else results/ when complete, else results_frozen/ (tfp.config.results_root).
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

from tfp import config as C, conformal as cp, io, transfer

ROOT = C.ROOT
RES = C.results_root()
REGEN = RES / "31_site_regen"
SITE = ROOT / "site" / "data"
# pre-fix values of every provenance entry, generated once from the pre-quantile-fix results of 2026-09-25 (not in the repository)
PREFIX_TABLE = ROOT / "docs" / "reconciliation" / "pre_quantile_fix_values.csv"
PLOTLY_VERSION = "2.35.2"
DATA_ACCESS = {
    "motrpac": {"what": "MotrpacRatTraining6moData 2.0.0 (portal release c1.0, rn6); GitHub commit f831a4fe421ec11687640452484a8247137aa74a",
                "how": "R/export_motrpac.R (make export)", "date": "2026-09-17", "documented_in": "docs/DATA_GUIDE.md"},
    "bodymap": {"what": "rat BodyMap GSE53960, Bioconductor bodymapRat 1.28.0 (ExperimentHub)", "how": "R/export_bodymap.R (make bodymap)",
                "date": "2026-09-17", "documented_in": "docs/EXTERNAL_VALIDATION.md"},
    "gtex": {"what": "GTEx v8 open-access (release 2017-06-05, RNASeQCv1.1.9): gene TPM GCT, sample attributes, gene read-count GCT",
             "how": "portal download + scripts/11_gtex_prepare.py (make gtex)", "date": "2026-09-17 (TPM, attributes); 2026-09-18 (read counts)",
             "documented_in": "docs/GTEX_TRANSFER.md"},
}


def configure(results: Path | None = None, out: Path | None = None):
    """Set the results root and the output directory (module globals used by every exporter)."""
    global RES, REGEN, SITE
    RES = C.results_root(results)
    REGEN = RES / "31_site_regen"
    if out:
        SITE = Path(out)

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
    """Canonical provenance path: files under the results root are always written as results/<path>, whichever
    root (results/ or results_frozen/) the export read them from; other paths are relative to the repository."""
    p = Path(p)
    if p == RES:
        return "results"
    try:
        return "results/" + str(p.relative_to(RES))
    except ValueError:
        pass
    try:
        return str(p.relative_to(ROOT))
    except ValueError:
        return str(p)


def gene_symbols(ids) -> dict:
    """feature_ID → gene symbol from the package annotation (data/raw/feature_to_gene.csv); without the data
    (a re-export from the snapshot) the symbols already committed in site/data/panel_curve.json are reused."""
    if (C.RAW_DIR / "feature_to_gene.csv").exists():
        return io.map_to_gene_symbols(ids)
    committed = ROOT / "site" / "data" / "panel_curve.json"
    known = json.loads(committed.read_text()).get("gene_symbols", {}) if committed.exists() else {}
    return {g: known.get(g, g) for g in ids}


def cluster_boot(ok, groups, n_boot: int = 2000, seed: int = C.SEED) -> list[float]:
    """95 % interval of a mean over individuals (animals / donors): a cluster percentile bootstrap, or, when every
    individual is at the same boundary (all right or all wrong) and the bootstrap cannot move, the exact
    Clopper–Pearson interval over the number of individuals whose samples are all right."""
    from scipy import stats
    ok = np.asarray(ok, dtype=float)
    groups = np.asarray(groups)
    u = np.unique(groups)
    idx = {g: np.flatnonzero(groups == g) for g in u}
    per_ind = np.array([ok[idx[g]].mean() for g in u])
    if np.all(per_ind == per_ind[0]) and per_ind[0] in (0.0, 1.0):
        k, n = int((per_ind == 1.0).sum()), len(u)
        lo = float(stats.beta.ppf(0.025, k, n - k + 1)) if k > 0 else 0.0
        hi = float(stats.beta.ppf(0.975, k + 1, n - k)) if k < n else 1.0
        return [lo, hi]
    rng = np.random.default_rng(seed)
    vals = []
    for _ in range(n_boot):
        pick = rng.choice(u, size=len(u), replace=True)
        sel = np.concatenate([idx[g] for g in pick])
        vals.append(ok[sel].mean())
    return [float(np.percentile(vals, 2.5)), float(np.percentile(vals, 97.5))]


def load_transfer_scores(phase: str):
    """The --save-scores outputs of phase 12 or 13: target probabilities (refit models), calibration scores, classes, organ map."""
    d = REGEN / phase
    probs = pd.read_csv(d / "scores_target_probs.csv", dtype=str)
    cal = pd.read_csv(d / "scores_calibration.csv", dtype={"viallabel": str, "pid": str})
    classes = json.loads((d / "classes.json").read_text())
    organ_map = json.loads((d / "organ_map.json").read_text())
    return probs, cal, classes, organ_map


def transfer_rung_stats(phase: str, model: str, stage_col: str, primary, indiv_col: str, alpha: float = 0.10):
    """Per target sample (primary stage, mapped organs): correctness of the all-animal model's call, correctness of the
    refit model's call, and the α-marginal set (refit model, source calibration): covered / empty / wrong; cluster
    bootstrap intervals over individuals."""
    probs, cal, classes, organ_map = load_transfer_scores(phase)
    pcols = [f"p_{c}" for c in classes]
    pm = probs[(probs["model"] == model) & (probs[stage_col].astype(str) == str(primary))].copy()
    mapped = pm["organ"].map(lambda o: bool(organ_map.get(o))).to_numpy()
    pm = pm[mapped]
    P = pm[pcols].to_numpy(dtype=float)
    sc = cal[cal["model"] == model]["score_lac"].to_numpy(dtype=float)
    q = cp.conformal_quantile(sc, alpha)
    sets = cp.predict_sets(P, q, "lac")
    ok_all = np.array([pa in organ_map[o] for pa, o in zip(pm["pred_all_animals"], pm["organ"])])
    pred_refit = np.array(classes)[P.argmax(axis=1)]
    ok_refit = np.array([pr in organ_map[o] for pr, o in zip(pred_refit, pm["organ"])])
    covered = np.array([any(sets[i, classes.index(t)] for t in organ_map[o] if t in classes) for i, o in enumerate(pm["organ"])])
    size = sets.sum(axis=1)
    groups = pm[indiv_col].to_numpy()
    return {"n_samples": int(len(pm)), "n_individuals": int(pm[indiv_col].nunique()),
            "accuracy": float(ok_all.mean()), "accuracy_ci": cluster_boot(ok_all, groups), "accuracy_refit": float(ok_refit.mean()),
            "coverage": float(covered.mean()), "coverage_ci": cluster_boot(covered, groups), "empty": float((size == 0).mean()),
            "wrong_non_empty": float(((~covered) & (size > 0)).mean()), "set_size": float(size.mean())}


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
def export_manifest(w: Writer, prov: Prov):
    """site/data/manifest.json: the result files this export actually read (from the provenance ledger), grouped by
    phase, each with size and sha256; the regeneration runs; software versions; data releases and access dates;
    and the site data files with a flag saying whether this run wrote them."""
    def entry(src: str) -> dict:
        f = RES / src.split("/", 1)[1]
        return {"path": src, "bytes": f.stat().st_size, "sha256": sha256(f)}
    sources = sorted(s for s in prov.sources if s.startswith("results/"))
    phases: dict[str, dict] = {}
    regen: dict[str, list] = {}
    for src in sources:
        parts = src.split("/")
        if parts[1] == "31_site_regen":
            regen.setdefault(parts[2], []).append(entry(src))
        else:
            phases.setdefault(parts[1], []).append(entry(src))
    phases = {k: {"files": v, "n_files": len(v), "bytes": int(sum(f["bytes"] for f in v))} for k, v in phases.items()}
    import sklearn
    versions = {"python": platform.python_version(), "pandas": pd.__version__, "numpy": np.__version__,
                "scikit-learn": sklearn.__version__, "plotly.js": PLOTLY_VERSION}
    try:
        versions["node"] = subprocess.check_output(["node", "--version"], text=True).strip()
    except Exception:
        pass
    results_dir = "results" if RES == C.RESULTS_DIR else ("results_frozen" if RES == C.FROZEN_DIR else str(RES))
    snapshot = None
    if (RES / "MANIFEST.json").exists():
        fm = json.loads((RES / "MANIFEST.json").read_text())
        snapshot = {k: fm.get(k) for k in ("snapshot_date", "frozen_from_commit", "n_files")}
    m = {"generated": w.generated, "git_hash": w.ghash, "results_dir": results_dir, "results_snapshot": snapshot,
         "phases": phases, "regeneration_runs": regen, "phases_used": sorted({s.split("/")[1] for s in sources}),
         "versions": versions,
         "data": {"motrpac": "MotrpacRatTraining6moData 2.0.0 (= portal release c1.0, rn6)",
                  "bodymap": "rat BodyMap GSE53960 via bodymapRat 1.28.0 (316 samples)", "gtex": "GTEx v8 (2017-06-05) open-access expression"},
         "data_access": DATA_ACCESS,
         "site_data_files": [{"name": p.name, "bytes": p.stat().st_size, "regenerated": p.name in w.written or p.name in ("provenance.json", "conformal_fixtures.json")}
                             for p in sorted(SITE.glob("*.json")) if p.name != "manifest.json"]}
    stale = [f["name"] for f in m["site_data_files"] if not f["regenerated"]]
    if stale:
        print(f"  WARNING: site/data files not written by this export: {', '.join(stale)}")
    (SITE / "manifest.json").write_text(json.dumps(m, indent=1))
    print(f"  wrote manifest.json: {len(phases)} phases, {len(regen)} regeneration runs, {len(m['site_data_files'])} site data files")
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
    # balanced accuracy of the same models (argmax of the exported test probabilities), per fold
    from sklearn.metrics import balanced_accuracy_score
    acc_rows = []
    for fold in sorted(probs["fold"].unique()):
        for model in MODELS:
            pt = probs[(probs["fold"] == fold) & (probs["model"] == model)]
            pred = np.array(classes)[pt[pcols].to_numpy(dtype=float).argmax(axis=1)]
            acc_rows.append({"fold": int(fold), "model": model, "balanced_accuracy": float(balanced_accuracy_score(pt["tissue"], pred)),
                             "n_test_animals": int(pt["pid"].nunique()), "n_test": int(len(pt))})
    rec.attrs["accuracy"] = pd.DataFrame(acc_rows)
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
    sym_sel = gene_symbols(sorted(sel["feature_ID"].unique()))
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
    # the age curve: n (samples, animals) and a 95 % cluster-bootstrap interval over animals per point, from the exported per-sample calls
    if (REGEN / "12_bodymap" / "scores_target_probs.csv").exists():
        for row in bm["age_accuracy"]:
            for m in MODELS:
                st = transfer_rung_stats("12_bodymap", m, "stage_weeks", row["stage_weeks"], "animal_id")
                assert abs(st["accuracy"] - row[m]) < 1e-9, ("accuracy from the exported calls disagrees with age_shift_accuracy.csv", m, row["stage_weeks"])
                row[f"{m}_ci"] = P.recomputed(f"bodymap_acc_{row['stage_weeks']}wk_{m}_ci", st["accuracy_ci"], ["31_site_regen/12_bodymap/scores_target_probs.csv"],
                                              f"95 % cluster bootstrap over the {st['n_individuals']} BodyMap animals at {row['stage_weeks']} weeks ({m})")
                row["n_samples"], row["n_animals"] = st["n_samples"], st["n_individuals"]
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

    # ---- identifiability (16) ---------------------------------------------
    d16 = RES / "16_identifiability"
    nesting = {}
    for f in sorted(d16.glob("nesting_*.csv")):
        assay = f.stem.replace("nesting_", "")
        nesting[assay] = P.table(f"nesting_{assay}", f"16_identifiability/{f.name}", "nesting.json", f"nesting.{assay}")
    layers = json.loads((d16 / "layers.json").read_text())
    bridge = {"status": "pending", "reason": "phase 16 was not run with --bridge (needs the portal RNA-seq count files with the reference-standard vials)"}
    if (d16 / "bridge_variance.csv").exists() and layers.get("bridge") and layers["bridge"].get("status") == "recomputed":
        bridge = {"status": "recomputed", "info": layers["bridge"],
                  "summary": P.table("bridge_variance", "16_identifiability/bridge_variance.csv", "nesting.json", "bridge.summary"),
                  "per_gene": P.table("bridge_variance_per_gene", "16_identifiability/bridge_variance_per_gene.csv", "nesting.json", "bridge.per_gene")}
        for bid, ptype in ((80001, 99), (80000, 88)):
            for gs in ("all_genes", "all_genes_expressed_in_pool", "panel_k20", "panel_k20_expressed_in_pool"):
                where = {"pool_bid": bid, "gene_set": gs}
                P.val(f"bridge_sum_ratio_{gs}_pool{ptype}", "16_identifiability/bridge_variance.csv", "sum_ratio_batch_over_tissue", where=where,
                      note=f"gastrocnemius-derived reference pool {bid} (type {ptype}) on 6 plates at both sites: Σ V_batch / Σ V_tissue over {gs}")
                P.val(f"bridge_median_ratio_{gs}_pool{ptype}", "16_identifiability/bridge_variance.csv", "median_ratio_batch_over_tissue", where=where)
                P.val(f"bridge_n_genes_{gs}_pool{ptype}", "16_identifiability/bridge_variance.csv", "n_genes", where=where)
        # one row per bridging pool (all genes), sorted by the ratio: the home page and the summary figure plot these
        allg = P.read("16_identifiability/bridge_variance.csv")
        allg = allg[allg["gene_set"] == "all_genes"].sort_values("sum_ratio_batch_over_tissue")
        bridge["pools"] = [
            {"pool_bid": int(r["pool_bid"]), "pool_type": int(r["pool_type"]), "pool_tissue": str(r["pool_tissue"]), "n_vials": int(r["n_vials"]),
             "n_plates": int(r["n_plates"]), "sites": str(r["sites"]), "n_genes": int(r["n_genes"]),
             "sum_ratio_batch_over_tissue": P.val(f"bridge_pool_{int(r['pool_bid'])}_sum_ratio_all_genes", "16_identifiability/bridge_variance.csv",
                                                  "sum_ratio_batch_over_tissue", where={"pool_bid": int(r["pool_bid"]), "gene_set": "all_genes"},
                                                  note=f"{r['pool_tissue']} reference pool {int(r['pool_bid'])} (type {int(r['pool_type'])}) on {int(r['n_plates'])} plates "
                                                       f"({r['sites']}): Σ V_batch / Σ V_tissue over all genes")}
            for _, r in allg.iterrows()]
    w.write("nesting.json", {"nesting": nesting,
                             "estimable_pairs": P.table("estimable_pairs", "16_identifiability/estimable_pairs.csv", "nesting.json", "estimable_pairs"),
                             "batch_counts": P.table("batch_counts", "16_identifiability/batch_counts.csv", "nesting.json", "batch_counts"),
                             "bridge": bridge,
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
    # ---- EDA: variance and batch partition of the leading PCs (03) ------------------------------------
    eda = {}
    for assay in ("TRNSCRPT", "PROT", "METAB"):
        eda[f"variance_{assay}"] = P.table(f"variance_partition_{assay}", f"03_eda/variance_partition_{assay}.csv", "eda.json", f"variance_{assay}")
        eda[f"batch_{assay}"] = P.table(f"batch_partition_{assay}", f"03_eda/batch_partition_{assay}.csv", "eda.json", f"batch_{assay}")
    eda["readout_variance"] = P.table("readout_variance_pc", "03_eda/readout_variance_pc1-3.csv", "eda.json", "readout_variance")
    eda["prot_diagnostic"] = P.table("prot_diagnostic_accuracy", "04_baselines/PROT/diagnostic_accuracy.csv", "eda.json", "prot_diagnostic")
    w.write("eda.json", eda, [rel(RES / f"03_eda/variance_partition_{a}.csv") for a in ("TRNSCRPT", "PROT", "METAB")] + [rel(RES / "04_baselines/PROT/diagnostic_accuracy.csv")])


# ---------------------------------------------------------------------------------------------
# headline: tiles and the transfer ladder
# ---------------------------------------------------------------------------------------------
def export_headline(w: Writer, prov: Prov, rec: pd.DataFrame | None):
    P = prov
    n_vials = int(P.val("n_trnscrpt_vials", "04_baselines/TRNSCRPT/per_fold.csv", "n_test", where={"model": "logreg_l2"}, agg="sum", note="the five test folds partition the vials"))
    n_animals = P.recomputed("n_trnscrpt_animals", int(P.read("04_baselines/TRNSCRPT/per_fold.csv").query("model == 'logreg_l2'").iloc[0][["n_train_animals", "n_test_animals"]].sum()),
                             ["04_baselines/TRNSCRPT/per_fold.csv"], "train + test animals of one fold")
    rc12_ = P.read("12_bodymap/recalibration.csv")
    n_bm_animals = P.recomputed("n_bodymap_adult_animals", int(rc12_["n_recal"].iloc[0] + rc12_["n_test_individuals"].iloc[0]), ["12_bodymap/recalibration.csv"], "recalibration + test individuals")
    n_bm_mapped = P.val("n_bodymap_adult_mapped", "12_bodymap/conformal_transfer.csv", "n_mapped", where={"stage_weeks": 21, "model": "k20", "conformal": "marginal"})
    tiles = [
        {"id": "tile_acc_k20", "value": P.val("tile_acc_k20", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 20}, agg="mean",
                                             note="mean over 5 animal-grouped folds; round-robin selector + logreg_l2"),
         "sd": P.val("tile_acc_k20_sd", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 20}, agg="std"),
         "label": "balanced accuracy of a 20-gene panel", "sub": f"19 tissues · {n_vials} vials · {n_animals} animals · 5 animal-grouped folds, the panel re-selected inside each", "format": "3",
         "source": "results/05_panels/TRNSCRPT/panel_curve.csv (k = 20, mean ± sd over folds)"},
        {"id": "tile_bodymap_k20", "value": P.val("tile_bodymap_k20", "12_bodymap/age_shift_accuracy.csv", "k20", where={"stage_weeks": 21}),
         "ci": (transfer_rung_stats("12_bodymap", "k20", "stage_weeks", 21, "animal_id")["accuracy_ci"] if (REGEN / "12_bodymap" / "scores_target_probs.csv").exists() else None),
         "label": "mapped adult organs named correctly in another lab's rats", "sub": f"rat BodyMap, 21-week adults: 9 of 11 organs mapped, muscle and brain as super-classes; {n_bm_mapped} samples, {n_bm_animals} animals; panel fit on all 50 MoTrPAC animals", "format": "3",
         "source": "results/12_bodymap/age_shift_accuracy.csv (stage 21, k20)"},
        {"id": "tile_bodymap_cov_k20", "value": P.val("tile_bodymap_cov_k20", "12_bodymap/conformal_transfer.csv", "coverage_mapped",
                                                      where={"stage_weeks": 21, "model": "k20", "conformal": "marginal"}),
         "ci": (transfer_rung_stats("12_bodymap", "k20", "stage_weeks", 21, "animal_id")["coverage_ci"] if (REGEN / "12_bodymap" / "scores_target_probs.csv").exists() else None),
         "empty": P.val("tile_bodymap_empty_k20", "12_bodymap/conformal_transfer.csv", "frac_empty_mapped", where={"stage_weeks": 21, "model": "k20", "conformal": "marginal"}),
         "label": "coverage of the 90 % guarantee there", "sub": "sets calibrated on MoTrPAC animals, α = 0.10; the shortfall is empty sets", "format": "3",
         "source": "results/12_bodymap/conformal_transfer.csv (stage 21, k20, marginal)"},
    ]
    cov_tile = tiles[-1]
    rt3 = P.read("31_site_regen/12_bodymap/recal_thresholds.csv")
    rt3 = rt3[(rt3["model"] == "k20") & (rt3["n_recal"] == 3)]
    recal_where = {"model": "k20", "n_recal": 3}
    recal_draws = P.val("tile_bodymap_recal_draws_k20", "12_bodymap/recalibration.csv", "draws", where=recal_where)
    recal_size = P.val("tile_bodymap_recal_size_k20", "12_bodymap/recalibration.csv", "set_size_recalibrated", where=recal_where)
    recal_n_test = P.val("tile_bodymap_recal_n_test_k20", "12_bodymap/recalibration.csv", "n_test_individuals", where=recal_where)
    recal_n_samples = P.recomputed("tile_bodymap_recal_n3_samples", float(rt3["n_cal_scores"].mean()), ["31_site_regen/12_bodymap/recal_thresholds.csv"],
                                   "mean number of calibration samples over the three-animal draws (model k20)")
    ci_txt = f" [{cov_tile['ci'][0]:.2f}, {cov_tile['ci'][1]:.2f}]" if cov_tile.get("ci") else ""
    tiles += [
        {"id": "tile_bodymap_recal_k20", "value": P.val("tile_bodymap_recal_k20", "12_bodymap/recalibration.csv", "coverage_recalibrated", where=recal_where,
                                                        note="coverage of the α = 0.10 sets on the 21-week BodyMap mapped organs after recalibrating the threshold on 3 of its animals; mean over the draws"),
         "line": f"mean of {int(recal_draws)} draws · {recal_size:.2f} tissue per set",
         "label": "observed coverage in another lab after recalibrating on three of its animals",
         "sub": f"sets calibrated on MoTrPAC animals cover {cov_tile['value']:.3f} there{ci_txt}, the shortfall empty sets; 3 calibration animals (≈ {round(recal_n_samples)} samples), {int(recal_n_test)} test animals, α = 0.10",
         "format": "3", "source": "results/12_bodymap/recalibration.csv (k20, n_recal 3: coverage_recalibrated, set_size_recalibrated, draws)"},
        {"id": "tile_estimable", "value": P.val("tile_estimable", "16_identifiability/estimable_pairs.csv", "n_pairs_estimable", where={"assay": "TRNSCRPT"}),
         "total": P.val("tile_estimable_total", "16_identifiability/estimable_pairs.csv", "n_pairs_total", where={"assay": "TRNSCRPT"}),
         "pairs": P.val("tile_estimable_pairs", "16_identifiability/estimable_pairs.csv", "estimable_pairs", where={"assay": "TRNSCRPT"}),
         "pair_text": " and ".join(t.lower() for t in str(P.val("tile_estimable_pairs_text", "16_identifiability/estimable_pairs.csv", "estimable_pairs", where={"assay": "TRNSCRPT"})).split("|")),
         "label": "tissue pairs whose contrast exists inside one processing batch (" + " and ".join(t.lower() for t in str(P.entries[-1]["value"]).split("|")) + ")",
         "sub": "RNA-seq plate, library batch and flowcell; the one pair is also the sex contrast", "format": "of",
         "source": "results/16_identifiability/estimable_pairs.csv (TRNSCRPT)"},
    ]
    # in-distribution accuracy per model
    acc = {"k20": (P.val("acc_k20", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 20}, agg="mean"),
                   P.val("acc_k20_sd", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 20}, agg="std")),
           "k50": (P.val("acc_k50", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 50}, agg="mean"),
                   P.val("acc_k50_sd", "05_panels/TRNSCRPT/panel_curve.csv", "balanced_accuracy", where={"k": 50}, agg="std")),
           "full": (P.val("acc_full", "04_baselines/TRNSCRPT/summary.csv", "balanced_accuracy_mean", where={"model": "logreg_l2"}),
                    P.val("acc_full_sd", "04_baselines/TRNSCRPT/summary.csv", "balanced_accuracy_std", where={"model": "logreg_l2"}))}
    acc_fclassif_k20 = P.val("acc_fclassif_k20", "05_panels/TRNSCRPT/panel_curve_fclassif.csv", "balanced_accuracy", where={"k": 20}, agg="mean", note="F-test selector at the same k")
    # the "1 of 171" tile belongs to the Identifiability page; the home page shows the bridge measurement instead
    tiles_identifiability = [t for t in tiles if t["id"] == "tile_estimable"]
    tiles = [t for t in tiles if t["id"] != "tile_estimable"]
    bridge_file = "16_identifiability/bridge_variance.csv"
    if (RES / bridge_file).exists():
        tiles.append({"id": "tile_bridge",
                      "value": P.val("tile_bridge", bridge_file, "sum_ratio_batch_over_tissue", where={"gene_set": "all_genes", "pool_bid": 80001},
                                     note="batch as a share of the tissue-separating variance: Σ V_batch / Σ V_tissue over all genes, reference pool 99 run on 6 plates at both sites"),
                      "label": "of the tissue signal is batch, measured directly on MoTrPAC's bridging reference pools", "format": "pct1",
                      "sub": (lambda ba: f"between-plate variance of gastrocnemius reference pool 99 over the variance separating the 19 tissues, all genes; 6 extraction plates, both sequencing sites; "
                                         f"{100 * P.recomputed('tile_bridge_pools_min', float(ba.min()), [bridge_file], 'smallest Σ V_batch / Σ V_tissue (all genes) over the bridging pools'):.1f}–"
                                         f"{100 * P.recomputed('tile_bridge_pools_max', float(ba.max()), [bridge_file], 'largest Σ V_batch / Σ V_tissue (all genes) over the bridging pools'):.1f} % across all "
                                         f"{P.recomputed('tile_bridge_n_pools', int(len(ba)), [bridge_file], 'bridging pools with a row in bridge_variance.csv')} pools")
                             (P.read(bridge_file).query("gene_set == 'all_genes'")["sum_ratio_batch_over_tissue"]),
                      "source": "results/16_identifiability/bridge_variance.csv (gene_set all_genes, pool 80001, Σ V_batch / Σ V_tissue)"})
    else:
        tiles.append({"id": "tile_bridge", "label": "batch, measured on MoTrPAC's bridging reference pools", "format": "pct1",
                      **P.pending("tile_bridge", "bridge measurement not computed (scripts/16_identifiability.py --bridge with the portal files)")})
    # the home page renders these four, in this order; the coverage tile stays in `tiles` for the Identifiability page
    home_tiles = ["tile_acc_k20", "tile_bodymap_k20", "tile_bodymap_recal_k20", "tile_bridge"]
    ladder = []

    def rung(**kw):
        ladder.append(kw)

    # 1. in-distribution: accuracy AND coverage of the same phase-06 models (fit 18 / calibrate 22 / test 10 animals per fold);
    #    the full model's coverage rows come from coverage.csv, the k20/k50 rows and every accuracy from the regenerated scores
    cov06 = P.read("06_conformal/TRNSCRPT/coverage.csv")
    acc06 = {}
    if rec is not None:
        a6 = rec.attrs["accuracy"]
        for model in MODELS:
            d = a6[a6["model"] == model]
            acc06[model] = (P.recomputed(f"acc06_{model}", float(d["balanced_accuracy"].mean()), ["31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv"],
                                         f"balanced accuracy of the phase-06 {model} models (18 fit animals per fold) on their 10 test animals, mean over 5 folds"),
                            float(d["balanced_accuracy"].std()))
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
                a_same = acc06.get(model, acc[model])
                rung(rung_id="in_distribution", label="in-distribution (held-out animals)", model=model, variant=variant, calibration=calib,
                     accuracy=a_same[0], accuracy_sd=a_same[1], accuracy_cv40=acc[model][0], accuracy_cv40_sd=acc[model][1],
                     coverage=cov, coverage_sd=sd, empty=empty, wrong_non_empty=max(0.0, 1 - cov - empty), set_size=size,
                     n_samples=n_vials, n_individuals=n_animals, n_calibration_animals=22,
                     source=[src, "results/31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv (accuracy of the same models)"])
    # 2. held-out sex (phase 08). k20 and full rows come from results/08; the k50 rows from the k = 50 rerun under
    #    results/31_site_regen/08_shift_k50 (same script, --k 50). Empty / wrong rates over SEEN-class vials and the
    #    cluster-bootstrap intervals come from the per-vial sets the --save-scores reruns export.
    st = P.read("08_shift/TRNSCRPT/shift_table.csv")
    k50_file = "31_site_regen/08_shift_k50/shift_table.csv"
    have_k50 = (RES / k50_file).exists()
    vials = {}
    for kk in ("k20", "k50"):
        f = REGEN / f"08_shift_{kk}" / "scores_target_vials.csv"
        if f.exists():
            vials[kk] = pd.read_csv(f, dtype={"viallabel": str, "pid": str})
    for split, label in (("train_control_test_trained", "training state: fit on sedentary control animals, tested on trained animals"),
                         ("train_male_test_female", "held-out sex: trained on males, tested on females"), ("train_female_test_male", "held-out sex: trained on females, tested on males")):
        for model in MODELS:
            arm = {"k20": "panel_k20", "k50": "panel_k50", "full": "full"}[model]
            tfile = k50_file if model == "k50" else "08_shift/TRNSCRPT/shift_table.csv"
            if model == "k50" and not have_k50:
                for variant in VARIANTS:
                    rung(rung_id=split, label=label, model=model, variant=variant, calibration="pooled",
                         **P.pending(f"cov_{split}_{model}_{variant}", "phase 08 has not been rerun with --k 50 (make regen-scores)"))
                continue
            vdf = vials.get("k50" if model == "k50" else "k20")
            v = vdf[(vdf["split"] == split) & (vdf["arm"] == arm)] if vdf is not None else None
            seen = v[v["seen"]] if v is not None else None
            srow0 = P.read(tfile)
            srow0 = srow0[(srow0["split"] == split) & (srow0["arm"] == arm)].iloc[0]
            # classes the source model knows: 19 minus the tissues absent from the source (one sex-specific tissue on the held-out sex,
            # none on the training-state split); the per-vial export carries one filled p_ column per model class, which must agree
            n_unseen = 0 if pd.isna(srow0["unseen_classes"]) else len([c for c in str(srow0["unseen_classes"]).replace(";", ",").split(",") if c.strip()])
            n_model_classes = 19 - n_unseen
            if v is not None and len(v):
                pcols = [c for c in v.columns if c.startswith("p_")]
                assert int(v[pcols].notna().all().sum()) == n_model_classes, ("model classes in the per-vial export disagree with unseen_classes", split, model)
            for variant in VARIANTS:
                col = {"marginal": "coverage_target_seen", "mondrian": "coverage_target_seen_mondrian", "floored": "coverage_target_seen_floored"}[variant]
                where = {"split": split, "arm": arm}
                cov = P.val(f"cov_{split}_{model}_{variant}", tfile, col, where=where)
                size_col = {"marginal": "avg_set_size_target", "mondrian": "avg_set_size_target_mondrian", "floored": "avg_set_size_target_floored"}[variant]
                srow = P.read(tfile)
                srow = srow[(srow["split"] == split) & (srow["arm"] == arm)].iloc[0]
                extra = {}
                if seen is not None and len(seen):
                    assert abs(float(seen[f"covered_{variant}"].mean()) - float(cov)) < 1e-9, ("per-vial sets disagree with shift_table", split, model, variant)
                    sz = seen[f"size_{variant}"].to_numpy()
                    cv = seen[f"covered_{variant}"].to_numpy(dtype=bool)
                    extra = {"n_samples_coverage": int(len(seen)), "empty": float((sz == 0).mean()), "wrong_non_empty": float(((~cv) & (sz > 0)).mean()),
                             "coverage_ci": cluster_boot(cv, seen["pid"].to_numpy()),
                             "accuracy_ci": cluster_boot((v["y_pred"] == v["tissue"]).to_numpy(), v["pid"].to_numpy())}
                    if variant == "marginal":
                        P.recomputed(f"empty_seen_{split}_{model}", extra["empty"], [f"31_site_regen/08_shift_{'k50' if model == 'k50' else 'k20'}/scores_target_vials.csv"],
                                     "empty-set rate over the target vials of seen classes (the unseen sex-specific tissue's vials excluded)")
                else:
                    n_seen = int(round(srow["n_test"] * srow["coverage_target_all"] / srow["coverage_target_seen"]))
                    extra = {"n_samples_coverage": n_seen, "empty": P.val(f"empty_{split}_{model}", tfile, "lac_frac_empty_target", where=where) if variant == "marginal" else None,
                             "empty_note": "over all target vials, including the unseen tissue"}
                rung(rung_id=split, label=label, model=model, variant=variant, calibration="pooled",
                     accuracy=P.val(f"acc_{split}_{model}", tfile, "accuracy_all", where=where) if variant == "marginal" else acc_cache[(split, model)],
                     accuracy_seen=P.val(f"accseen_{split}_{model}", tfile, "bal_acc_seen", where=where, note="balanced accuracy over the seen classes") if variant == "marginal" else accseen_cache[(split, model)],
                     n_classes_seen=n_model_classes, accuracy_sd=None, coverage=cov, coverage_sd=None,
                     set_size=P.val(f"size_{split}_{model}_{variant}", tfile, size_col, where=where),
                     coverage_source=P.val(f"covsrc_{split}_{model}", tfile, "coverage_source_id", where=where) if variant == "marginal" else None,
                     recal_n3=P.val(f"recal3_{split}_{model}", tfile, "cov_target_recal_N3", where=where) if variant == "marginal" else None,
                     unseen=P.val(f"unseen_{split}_{model}", tfile, "unseen_classes", where=where) if (variant == "marginal" and n_unseen) else None,
                     n_samples=int(srow["n_test"]), n_individuals=int(srow["n_test_animals"]), n_calibration_animals=int(srow["n_cal_animals"]),
                     n_source_animals=int(srow["n_train_animals"]),
                     source=[f"results/{tfile}"] + ([f"results/31_site_regen/08_shift_{'k50' if model == 'k50' else 'k20'}/scores_target_vials.csv"] if seen is not None else []),
                     **extra)
                if variant == "marginal":
                    acc_cache[(split, model)] = ladder[-1]["accuracy"]
                    accseen_cache[(split, model)] = ladder[-1]["accuracy_seen"]
    # 3. different lab: BodyMap adults
    rc12 = P.read("12_bodymap/recalibration.csv")
    n_adult_animals = int(rc12["n_recal"].iloc[0] + rc12["n_test_individuals"].iloc[0])
    for model in MODELS:
        a = P.val(f"acc_bodymap_{model}", "12_bodymap/age_shift_accuracy.csv", model, where={"stage_weeks": 21})
        bstat = transfer_rung_stats("12_bodymap", model, "stage_weeks", 21, "animal_id") if (REGEN / "12_bodymap" / "scores_target_probs.csv").exists() else None
        if bstat:
            assert abs(bstat["accuracy"] - a) < 1e-9, ("all-animal accuracy from the exported calls disagrees with age_shift_accuracy.csv", model)
            P.recomputed(f"acc_bodymap_{model}_ci", bstat["accuracy_ci"], ["31_site_regen/12_bodymap/scores_target_probs.csv"], "95 % cluster bootstrap over the 8 adult animals")
        for variant in VARIANTS:
            where = {"stage_weeks": 21, "model": model, "conformal": variant}
            rung(rung_id="different_lab", label="different laboratory: rat BodyMap adults (same strain)", model=model, variant=variant, calibration="pooled",
                 accuracy=a, accuracy_sd=None, accuracy_ci=bstat["accuracy_ci"] if bstat else None, accuracy_refit=bstat["accuracy_refit"] if bstat else None,
                 coverage_ci=bstat["coverage_ci"] if (bstat and variant == "marginal") else None, wrong_non_empty=bstat["wrong_non_empty"] if (bstat and variant == "marginal") else None,
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
        gstat = transfer_rung_stats("13_gtex", model, "stage", "adult", "donor") if (REGEN / "13_gtex" / "scores_target_probs.csv").exists() else None
        if gstat:
            assert abs(gstat["accuracy"] - a) < 1e-9, ("all-animal accuracy from the exported calls disagrees with accuracy_overall.csv", model)
            P.recomputed(f"acc_gtex_{model}_ci", gstat["accuracy_ci"], ["31_site_regen/13_gtex/scores_target_probs.csv"], "95 % cluster bootstrap over the 862 donors")
        for variant in VARIANTS:
            where = {"stage": "adult", "model": model, "conformal": variant}
            rung(rung_id="different_species", label="different species: human GTEx v8 (1:1 orthologs)", model=model, variant=variant, calibration="pooled",
                 accuracy=a, accuracy_sd=None, accuracy_ci=gstat["accuracy_ci"] if gstat else None, accuracy_refit=gstat["accuracy_refit"] if gstat else None,
                 coverage_ci=gstat["coverage_ci"] if (gstat and variant == "marginal") else None, wrong_non_empty=gstat["wrong_non_empty"] if (gstat and variant == "marginal") else None,
                 coverage=P.val(f"cov_gtex_{model}_{variant}", "13_gtex/conformal_transfer.csv", "coverage_mapped", where=where), coverage_sd=None,
                 empty=P.val(f"empty_gtex_{model}_{variant}", "13_gtex/conformal_transfer.csv", "frac_empty_mapped", where=where),
                 set_size=P.val(f"size_gtex_{model}_{variant}", "13_gtex/conformal_transfer.csv", "avg_set_size_mapped", where=where),
                 recal_n3=P.val(f"recal3_gtex_{model}", "13_gtex/recalibration.csv", "coverage_recalibrated", where={"model": model, "n_recal": 3}) if variant == "marginal" else None,
                 recal_n3_size=P.val(f"recal3size_gtex_{model}", "13_gtex/recalibration.csv", "set_size_recalibrated", where={"model": model, "n_recal": 3}) if variant == "marginal" else None,
                 n_samples=P.val(f"n_gtex_{model}_{variant}", "13_gtex/conformal_transfer.csv", "n_mapped", where=where), n_individuals=int(go["gtex_donors"].iloc[0]), n_calibration_animals=15,
                 source=["results/13_gtex/accuracy_overall.csv", "results/13_gtex/conformal_transfer.csv", "results/13_gtex/recalibration.csv"])
    # a few more headline numbers used in prose
    extras = {
        "acc_fclassif_k20": acc_fclassif_k20,
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
        "qc_technical": P.val("qc_technical", "16_identifiability/qc_only_summary.csv", "bal_acc_mean", where={"features": "technical"}, note="balanced accuracy, the metric of every other accuracy on the site"),
        "qc_technical_sd": P.val("qc_technical_sd", "16_identifiability/qc_only_summary.csv", "bal_acc_sd", where={"features": "technical"}),
        "qc_composition": P.val("qc_composition", "16_identifiability/qc_only_summary.csv", "bal_acc_mean", where={"features": "composition"}, note="balanced accuracy, the metric of every other accuracy on the site"),
        "qc_composition_sd": P.val("qc_composition_sd", "16_identifiability/qc_only_summary.csv", "bal_acc_sd", where={"features": "composition"}),
        "qc_all": P.val("qc_all", "16_identifiability/qc_only_summary.csv", "bal_acc_mean", where={"features": "all"}, note="balanced accuracy, the metric of every other accuracy on the site"),
        "qc_all_sd": P.val("qc_all_sd", "16_identifiability/qc_only_summary.csv", "bal_acc_sd", where={"features": "all"}),
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
    bv = RES / "16_identifiability" / "bridge_variance.csv"
    if bv.exists():
        b = pd.read_csv(bv)
        sel = b[(b["pool_bid"] == 80001) & (b["gene_set"] == "all_genes")]
        if len(sel) == 1:
            extras["bridge_sum_ratio_all_genes_pool99"] = float(sel["sum_ratio_batch_over_tissue"].iloc[0])
            extras["bridge_n_plates_pool99"] = int(sel["n_plates"].iloc[0])
        allg = b[b["gene_set"] == "all_genes"]
        extras["bridge_pools_min_sum_ratio_all_genes"] = P.recomputed("bridge_pools_min_sum_ratio_all_genes", float(allg["sum_ratio_batch_over_tissue"].min()),
                                                                      ["16_identifiability/bridge_variance.csv"], "smallest Σ V_batch / Σ V_tissue (all genes) over the bridging pools")
        extras["bridge_pools_max_sum_ratio_all_genes"] = P.recomputed("bridge_pools_max_sum_ratio_all_genes", float(allg["sum_ratio_batch_over_tissue"].max()),
                                                                      ["16_identifiability/bridge_variance.csv"], "largest Σ V_batch / Σ V_tissue (all genes) over the bridging pools")
        extras["bridge_n_pools"] = int(len(allg))
        sel = b[(b["pool_bid"] == 80001) & (b["gene_set"] == "panel_k20_expressed_in_pool")]
        if len(sel) == 1:
            extras["bridge_sum_ratio_panel_expressed_pool99"] = float(sel["sum_ratio_batch_over_tissue"].iloc[0])
            extras["bridge_n_panel_expressed_pool99"] = int(sel["n_genes"].iloc[0])
    else:
        extras["bridge_sum_ratio_all_genes_pool99"] = P.pending("bridge_sum_ratio_all_genes_pool99", "phase 16 not run with --bridge")
    rt = REGEN / "13_gtex" / "recal_thresholds.csv"
    if rt.exists():
        r = pd.read_csv(rt)
        d = r[(r["model"] == "k20") & (r["n_recal"] == 3)]
        extras["gtex_recal_k20_n3_frac_inf"] = P.recomputed("gtex_recal_k20_n3_frac_inf", float(np.isinf(d["q_t"]).mean()), ["31_site_regen/13_gtex/recal_thresholds.csv"],
                                                            "fraction of the 20 three-donor recalibration draws (k20) whose threshold is +∞")
    else:
        extras["gtex_recal_k20_n3_frac_inf"] = P.pending("gtex_recal_k20_n3_frac_inf", "per-draw thresholds not regenerated")
    # the design constants the diagrams state (read from the result tables where they are recorded; pipeline constants otherwise)
    from tfp import models as _models
    import inspect as _inspect
    pc = P.read("05_panels/TRNSCRPT/panel_curve.csv")
    pc20 = pc[pc["k"] == 20]
    cov06_ = P.read("06_conformal/TRNSCRPT/coverage.csv")
    design = {"n_outer_folds": int(pc20["fold"].nunique()) if "fold" in pc20 else int(len(pc20)),
              "n_train_animals": int(round(pc20["n_train_animals"].mean())), "n_test_animals": int(round(pc20["n_test_animals"].mean())),
              "n_fit_animals": int(round(cov06_["n_fit_animals"].mean())), "n_cal_animals": int(round(cov06_["n_cal_animals"].mean())),
              "alpha": 0.1, "k_panel": 20, "variance_prefilter": int(_inspect.signature(_models.make_pipeline).parameters["prefilter"].default),
              "C_grid": list(_models.C_GRID), "inner_splits": int(_inspect.signature(_models.fit_tuned).parameters["inner_splits"].default),
              "source": "results/05_panels/TRNSCRPT/panel_curve.csv (folds, animals), results/06_conformal/TRNSCRPT/coverage.csv (fit / calibration animals), src/tfp/models.py (prefilter, C grid, inner splits)"}
    assert design["n_train_animals"] + design["n_test_animals"] == n_animals, design
    extras["n_vials"], extras["n_animals"] = int(n_vials), int(n_animals)
    # vena cava → brown fat calls of the in-distribution 20-gene model, and how many of those vials the consortium flagged
    sp = REGEN / "06_conformal" / "TRNSCRPT" / "scores_test_probs.csv"
    flf = "15_time_course/design/flagged_vials.csv"
    if sp.exists() and (RES / flf).exists():
        sc = pd.read_csv(sp, dtype={"viallabel": str, "pid": str})
        sc = sc[sc["model"] == "k20"]
        pcols_ = [c for c in sc.columns if c.startswith("p_")]
        sc = sc.assign(y_pred=sc[pcols_].idxmax(axis=1).str[2:])
        vb = sc[(sc["tissue"] == "VENACV") & (sc["y_pred"] == "BAT")]
        flagged = set(P.read(flf)["viallabel"].astype(str))
        srcs = ["31_site_regen/06_conformal/TRNSCRPT/scores_test_probs.csv", flf]
        extras["venacv_bat_calls"] = P.recomputed("venacv_bat_calls", int(len(vb)), srcs[:1], "held-out vena cava vials the 20-gene model calls brown fat (phase-06 design, all folds)")
        extras["venacv_bat_calls_flagged"] = P.recomputed("venacv_bat_calls_flagged", int(vb["viallabel"].isin(flagged).sum()), srcs, "of those, vials the consortium flagged as brown-fat contaminated")
        extras["venacv_flagged_total"] = P.val("venacv_flagged_total", flf, "viallabel", agg="count", note="vena cava vials flagged as brown-fat contaminated by the consortium (all female, 1 and 2 weeks)")
        extras["venacv_vials"] = P.recomputed("venacv_vials", int((sc["tissue"] == "VENACV").sum()), srcs[:1], "vena cava vials scored in the phase-06 design")
    else:
        for k in ("venacv_bat_calls", "venacv_bat_calls_flagged", "venacv_flagged_total", "venacv_vials"):
            extras[k] = None
            P.pending(k, "phase 15 design table flagged_vials.csv or the phase-06 regeneration missing")
    extras["n_tissues"] = int(P.val("n_tissues", "16_identifiability/estimable_pairs.csv", "n_tissues", where={"assay": "TRNSCRPT"}))
    w.write("headline.json", {"question": "Can a molecular signature identify a tissue reliably?",
                              "tiles": tiles, "home_tiles": home_tiles, "tiles_identifiability": tiles_identifiability, "accuracy": {m: {"mean": a[0], "sd": a[1]} for m, a in acc.items()}, "ladder": ladder, "extras": extras, "design": design,
                              "rung_order": ["in_distribution", "train_control_test_trained", "train_male_test_female", "train_female_test_male", "different_lab", "different_species"]},
            sorted(prov.sources))


acc_cache: dict = {}
accseen_cache: dict = {}


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
    if (RES / "34_panel_model" / "genes.csv").exists():        # the scoring tool's model (35-animal refit) may select genes outside the all-animal panel
        extra = P.read("34_panel_model/genes.csv")["feature_ID"].astype(str).tolist()
        sets["scoring_model"] = extra
        genes = list(dict.fromkeys(list(genes) + extra))
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
                                  "values": [[round(float(x), 4) for x in Xm[g].to_numpy()] for g in present],
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
                                  "values": [[round(float(x), 4) for x in lb[g].to_numpy()] for g in present_b],
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
# the scoring tool: the 20-gene transfer model, and the panel card
# ---------------------------------------------------------------------------------------------
def export_panel_model(w: Writer, prov: Prov):
    """site/data/panel_model.json — the transfer model that carries the guarantee (phase 34): genes in feature order,
    multinomial logistic coefficients and intercepts, the MoTrPAC z-scoring statistics (fallback for small uploads),
    the 15-animal calibration scores, the all-animal call model, and the validation against the regeneration."""
    P = prov
    before = set(P.sources)
    d = "34_panel_model"
    if not (RES / d / "coefficients.csv").exists():
        print("  panel model: results/34_panel_model absent (scripts/34_panel_model.py); skipped")
        return
    def block(prefix):
        genes = P.table(f"panel_model_{prefix}genes", f"{d}/{prefix}genes.csv", "panel_model.json", f"{prefix or 'scoring.'}genes")
        coef = P.table(f"panel_model_{prefix}coefficients", f"{d}/{prefix}coefficients.csv", "panel_model.json", f"{prefix or 'scoring.'}coefficients")
        sz = P.table(f"panel_model_{prefix}source_z", f"{d}/{prefix}source_z.csv", "panel_model.json", f"{prefix or 'scoring.'}source_z")
        ids = [str(g["feature_ID"]) for g in genes]
        return {"classes": [str(r["class"]) for r in coef],
                "genes": [{"id": str(g["feature_ID"]), "symbol": g["symbol"], "in_all_animal_panel": bool(g["in_all_animal_panel"])} for g in genes],
                "coef": [[float(r[g]) for g in ids] for r in coef], "intercept": [float(r["intercept"]) for r in coef],
                "source_z": {"mean": [float(r["mean"]) for r in sz], "scale": [float(r["scale"]) for r in sz]}}
    model = block("")
    call = block("call_model_")
    cal = P.read("31_site_regen/12_bodymap/scores_calibration.csv")
    cal = cal[cal["model"] == "k20"]
    classes = model["classes"]
    calibration = {"scores": [float(v) for v in cal["score_lac"]], "y_idx": [classes.index(t) for t in cal["tissue"].astype(str)],
                   "n_animals": int(cal["pid"].nunique()), "n_vials": int(len(cal)), "alpha_default": 0.1,
                   "source": "results/31_site_regen/12_bodymap/scores_calibration.csv (model k20; the 15 calibration animals, one vial per tissue)"}
    val = P.table("panel_model_validation", f"{d}/validation.csv", "panel_model.json", "validation")
    n_shared = P.val("panel_model_genes_shared", f"{d}/genes.csv", "feature_ID", where={"in_all_animal_panel": True}, agg="count",
                     note="genes of the scoring model (35-animal refit) that are also in the all-animal 20-gene panel")
    rt = P.read("31_site_regen/12_bodymap/recal_thresholds.csv")
    d3 = rt[(rt["model"] == "k20") & (rt["n_recal"] == 3)]
    context = {"bodymap_coverage_k20": P.val("panel_model_bodymap_coverage_k20", "12_bodymap/conformal_transfer.csv", "coverage_mapped",
                                             where={"stage_weeks": 21, "model": "k20", "conformal": "marginal"},
                                             note="coverage of the α = 0.10 sets on the 21-week BodyMap mapped organs with the MoTrPAC calibration (the home-page tile)"),
               "recal_n3_mean_cal_samples": P.recomputed("panel_model_recal_n3_mean_cal_samples", float(d3["n_cal_scores"].mean()),
                                                         ["31_site_regen/12_bodymap/recal_thresholds.csv"],
                                                         "mean number of calibration samples over the three-animal BodyMap recalibration draws (model k20): what 'three animals' meant in the study"),
               "recal_n3_cal_samples_range": [int(d3["n_cal_scores"].min()), int(d3["n_cal_scores"].max())], "recal_n3_draws": int(len(d3))}
    w.write("panel_model.json", {**model, "calibration": calibration, "alpha_default": 0.1, "call_model": call, "validation": val, "context": context,
                                 "n_genes_shared_with_all_animal_panel": n_shared,
                                 "design": "PanelModels (scripts/12_bodymap_validate.py): z-score representation, multinomial logistic regression C = 0.1, "
                                           "fit on the 35 MoTrPAC animals outside the calibration split; sets calibrated on the other 15 animals; a new "
                                           "sample is z-scored per gene within its own dataset (population sd) and scored by softmax(z · coefᵀ + intercept)"},
            sorted(P.sources - before))
    return n_shared


def export_panel_card(w: Writer):
    """site/data/panel_card.json and .csv — the all-animal 20-gene panel: ids, symbols, human orthologs, marker tissue,
    effect size, bootstrap frequency, risk flags, BodyMap / GTEx status, and the mean log2 CPM per tissue (the reference
    profiles), from the site's own gene and expression exports."""
    G = json.loads((SITE / "genes.json").read_text())
    Ex = json.loads((SITE / "expr_motrpac.json").read_text())
    k20 = list(G["sets"]["k20"])
    tissues = sorted({s["tissue"] for s in Ex["samples"]})
    gi = {g: i for i, g in enumerate(Ex["genes"])}
    by_t = {t: [j for j, s in enumerate(Ex["samples"]) if s["tissue"] == t] for t in tissues}
    def status(fail, weak):
        return "not tested" if fail is None else ("lost" if fail else ("weakened" if weak else "holds"))
    rows = []
    for g in (x for x in G["genes"] if x["id"] in k20):
        prof = {t: (float(np.mean([Ex["values"][gi[g["id"]]][j] for j in by_t[t]])) if g["id"] in gi else None) for t in tissues}
        rows.append({"id": g["id"], "symbol": g["symbol"], "human_gene": g.get("human_gene"), "marker_tissue": g.get("marker_tissue"),
                     "effect_size_log2cpm": g.get("effect_size"), "bootstrap_frequency": g.get("freq"), "training_regulated": g.get("regulated"),
                     "qc_correlated": g.get("qc_flag"), "bodymap": status(g.get("fails_bodymap"), g.get("weakened_bodymap")),
                     "gtex": status(g.get("fails_gtex"), g.get("weakened_gtex")), "profile_mean_log2cpm": prof})
    sources = ["site/data/genes.json (results/05_panels/TRNSCRPT/stability_k20_annotated.csv, candidate_panel_annotated.csv, results/12_bodymap/panel_gene_check.csv, results/13_gtex/panel_gene_check.csv)",
               "site/data/expr_motrpac.json (data/raw/counts, log2 CPM)"]
    w.write("panel_card.json", {"panel": rows, "tissues": tissues, "unit": Ex["unit"], "n_samples_per_tissue": {t: len(v) for t, v in by_t.items()}}, sources)
    flat = []
    for r in rows:
        f = {k: v for k, v in r.items() if k != "profile_mean_log2cpm"}
        f.update({f"mean_log2cpm_{t}": (None if r["profile_mean_log2cpm"][t] is None else round(r["profile_mean_log2cpm"][t], 3)) for t in tissues})
        flat.append(f)
    pd.DataFrame(flat).to_csv(SITE / "panel_card.csv", index=False)
    print(f"  wrote panel_card.csv ({len(flat)} genes × {len(tissues)} tissues)")


# ---------------------------------------------------------------------------------------------
# the Exercise page: what training does and does not do to the fingerprint
# ---------------------------------------------------------------------------------------------
def export_exercise(w: Writer, prov: Prov):
    """site/data/exercise.json — every number of the Exercise page: within-tissue separability of sedentary vs
    trained animals and its batch check (phase 07), where the training signal sits (phase 03), the fingerprint fit
    on controls only tested per training duration (recomputed from the phase-08 per-vial regeneration; enriched from
    phase 15 when present), the study-design dates and physiology (phase 15) and the panel genes' training response
    (phase 05b). `summary` holds every scalar a page reads, each with a provenance id."""
    P = prov
    S: dict = {}
    before = set(P.sources)
    # ---- 1. control vs trained within tissue (phase 07) -----------------------------------------------------
    sep = {"duration": P.table("taskB_duration", "07_fusion/taskB_duration_summary.csv", "exercise.json", "separability.duration"),
           "best_vs_null": P.table("taskB_best_vs_null", "07_fusion/taskB_best_vs_null.csv", "exercise.json", "separability.best_vs_null"),
           "arms_8w": P.table("taskB_summary", "07_fusion/taskB_summary.csv", "exercise.json", "separability.arms_8w"),
           "plex_balance": P.table("taskB_batch_balance", "07_fusion/taskB_batch_balance.csv", "exercise.json", "separability.plex_balance"),
           "design": "control vs 8-week within each of 7 tissues (17–20 animals), animal-grouped sex × group-stratified folds, tuned models; "
                     "13 arms (single omic and fusion); the null is the 95th percentile of the best arm's AUROC under within-sex label permutations"}
    S["fusion_n_beats_single"] = int(P.val("fusion_n_beats_single", "07_fusion/taskB_best_vs_null.csv", "fusion_beats_single_by_gt_sd", agg="sum",
                                           note="tissues where a fusion arm beats the best single omic by more than its fold sd"))
    S["fusion_n_beats_null"] = int(P.val("fusion_n_beats_null", "07_fusion/taskB_best_vs_null.csv", "best_beats_null_p95", agg="sum",
                                         note="tissues where the best arm beats the permutation null (95th percentile of the max over arms)"))
    S["fusion_n_tissues"] = P.val("fusion_n_tissues", "07_fusion/taskB_best_vs_null.csv", "tissue", agg="count")
    for d in ("1w", "2w", "4w", "8w"):
        S[f"taskB_mean_auroc_{d}"] = P.val(f"taskB_mean_auroc_{d}", "07_fusion/taskB_duration_summary.csv", "best_single_auroc", where={"duration": d}, agg="mean",
                                           note="mean over the 7 tissues of the best single-omic AUROC (best of 9 arms, optimistic); 1w/2w/4w animals were sacrificed months apart from the controls")
    # ---- 2. processing covariates alone (phase 07 batch check) -----------------------------------------------
    cov = {"auroc": P.table("batch_covariate_auroc", "07_fusion/batch_covariate_auroc.csv", "exercise.json", "covariates.auroc"),
           "conclusion": P.table("batch_conclusion", "07_fusion/batch_conclusion.csv", "exercise.json", "covariates.conclusion"),
           "feature_sets": {"pheno_collection": "collection (arrival, sacrifice and collection dates, times, staff, cage, freeze times)",
                            "trnscrpt_library": "library preparation (plate, dates, concentrations, well position, sequencing run)",
                            "trnscrpt_depth": "sequencing depth (raw, mapped and splice-junction read counts)",
                            "trnscrpt_qc": "library QC fractions (mapping, mRNA, intronic, rRNA, globin, GC, duplication, RIN, 5′–3′ bias)",
                            "prot_plex_channel": "TMT plex and channel", "all_batch_covariates": "all of the above"}}
    S["verdict_n_training"] = P.val("verdict_n_training", "07_fusion/batch_conclusion.csv", "tissue", where={"verdict": "training"}, agg="count",
                                    note="tissues where no processing covariate separates the arms above its permutation null")
    S["verdict_n_tissues"] = P.val("verdict_n_tissues", "07_fusion/batch_conclusion.csv", "tissue", agg="count")
    S["verdict_overall"] = P.val("verdict_overall", "07_fusion/batch_conclusion.csv", "overall_verdict", where={"tissue": "CORTEX"}, note="the same string in every row")
    for id_, t, fs in (("cov_heart_depth", "HEART", "trnscrpt_depth"), ("cov_kidney_library", "KIDNEY", "trnscrpt_library"), ("cov_skmgn_qc", "SKM-GN", "trnscrpt_qc")):
        where = {"tissue": t, "feature_set": fs, "model": "logreg"}
        S[f"{id_}_logreg"] = P.val(f"{id_}_logreg", "07_fusion/batch_covariate_auroc.csv", "auroc_mean", where=where)
        S[f"{id_}_null95"] = P.val(f"{id_}_null95", "07_fusion/batch_covariate_auroc.csv", "null_p95_auroc", where=where)
        S[f"{id_}_p"] = P.val(f"{id_}_p", "07_fusion/batch_covariate_auroc.csv", "p_perm", where=where)
    S["cov_skmgn_all_rf"] = P.val("cov_skmgn_all_rf", "07_fusion/batch_covariate_auroc.csv", "auroc_mean", where={"tissue": "SKM-GN", "feature_set": "all_batch_covariates", "model": "rf"})
    # ---- 3. where the training signal sits (phase 03, within-tissue PCA) -----------------------------------------
    wt = {a: P.table(f"within_tissue_pca_{a}", f"03_eda/within_tissue_pca_{a}.csv", "exercise.json", f"within_tissue.{a}") for a in ("TRNSCRPT", "PROT", "METAB")}
    wt["readout"] = P.table("readout_within_tissue", "03_eda/readout_within_tissue_summary.csv", "exercise.json", "within_tissue.readout")
    for a in ("TRNSCRPT", "PROT", "METAB"):
        S[f"wt_{a}_n_visible"] = P.val(f"wt_{a}_n_visible", "03_eda/readout_within_tissue_summary.csv", "group_visible", where={"assay": a},
                                       note="tissues whose training-group R² on some PC 1–5 exceeds the 95th percentile of label permutations")
        S[f"wt_{a}_n_tissues"] = P.val(f"wt_{a}_n_tissues", "03_eda/readout_within_tissue_summary.csv", "tissues", where={"assay": a})
        S[f"wt_{a}_median_max_r2"] = P.val(f"wt_{a}_median_max_r2", "03_eda/readout_within_tissue_summary.csv", "median_max_R2_group", where={"assay": a})
    S["wt_heart_trnscrpt_max_r2_group"] = P.val("wt_heart_trnscrpt_max_r2_group", "03_eda/within_tissue_pca_TRNSCRPT.csv", "max_R2_group_PC1-5", where={"tissue": "HEART"})
    S["wt_heart_trnscrpt_pc"] = P.val("wt_heart_trnscrpt_pc", "03_eda/within_tissue_pca_TRNSCRPT.csv", "PC_with_max_group_R2", where={"tissue": "HEART"})
    S["stacked_pc1_r2_group"] = P.val("stacked_pc1_r2_group", "03_eda/variance_partition_TRNSCRPT.csv", "R2_group", where={"PC": "PC1"})
    S["stacked_pc1_r2_tissue"] = P.val("stacked_pc1_r2_tissue", "03_eda/variance_partition_TRNSCRPT.csv", "R2_tissue", where={"PC": "PC1"})
    # ---- 4. the fingerprint fit on controls only, by training duration (phase-08 per-vial regeneration) ------------
    fbd, errors = [], []
    for model, kk, arm in (("k20", "k20", "panel_k20"), ("k50", "k50", "panel_k50"), ("full", "k20", "full")):
        f = REGEN / f"08_shift_{kk}" / "scores_target_vials.csv"
        rf = f"31_site_regen/08_shift_{kk}/scores_target_vials.csv"
        if not f.exists():
            fbd.append({"model": model, **P.pending(f"fbd_acc_{model}_8w", "phase 08 not regenerated with --save-scores (make regen-scores)")})
            continue
        v = pd.read_csv(f, dtype={"viallabel": str, "pid": str})
        v = v[(v["split"] == "train_control_test_trained") & (v["arm"] == arm)]
        for grp in ("1w", "2w", "4w", "8w"):
            d = v[v["group"] == grp]
            ok = (d["y_pred"] == d["tissue"]).to_numpy()
            cv = d["covered_marginal"].to_numpy(dtype=bool)
            sz = d["size_marginal"].to_numpy()
            row = {"model": model, "group": grp, "n_samples": int(len(d)), "n_individuals": int(d["pid"].nunique()),
                   "accuracy": P.recomputed(f"fbd_acc_{model}_{grp}", float(ok.mean()), [rf], f"accuracy of the control-fit {model} model on the {grp} animals' vials (train_control_test_trained split)"),
                   "accuracy_ci": cluster_boot(ok, d["pid"].to_numpy()),
                   "coverage": P.recomputed(f"fbd_cov_{model}_{grp}", float(cv.mean()), [rf], f"coverage of the α = 0.10 marginal sets on the {grp} animals' vials"),
                   "coverage_ci": cluster_boot(cv, d["pid"].to_numpy()),
                   "empty": float((sz == 0).mean()), "wrong_non_empty": float(((~cv) & (sz > 0)).mean()), "set_size": float(sz.mean()),
                   "source": [f"results/{rf}", "results/08_shift/TRNSCRPT/shift_table.csv"]}
            fbd.append(row)
            e = d[d["y_pred"] != d["tissue"]].groupby(["tissue", "y_pred", "sex"]).size().reset_index(name="n")
            errors += [{"model": model, "group": grp, **r} for r in records(e)]
            if model in ("k20", "full") and grp in ("1w", "8w"):
                S[f"fbd_acc_{model}_{grp}"], S[f"fbd_cov_{model}_{grp}"] = row["accuracy"], row["coverage"]
    fbd_block = {"rows": fbd, "errors_by_group": errors,
                 "design": "fit on 7 of the 10 sedentary control animals, calibrate on the other 3 (54 pooled vials, LAC, α = 0.10), test every vial of each "
                           "trained group; the phase-08 train_controls_test_trained split; 95 % cluster-bootstrap intervals over the 10 animals of a group"}
    p15 = RES / "15_time_course" / "fingerprint_by_duration"
    if (p15 / "part5_by_duration.csv").exists():
        p5 = P.read("15_time_course/fingerprint_by_duration/part5_by_duration.csv")
        armname = {"k20": "panel_k20", "full": "full"}
        for r in fbd:
            if r.get("pending") or r["model"] not in armname:
                continue
            sel = p5[(p5["design"] == "control_only_fit7_cal3") & (p5["arm"] == armname[r["model"]]) & (p5["test_group"] == r["group"])]
            assert len(sel) == 1 and abs(float(sel["accuracy"].iloc[0]) - r["accuracy"]) < 1e-9 and abs(float(sel["coverage"].iloc[0]) - r["coverage"]) < 1e-9, \
                ("phase 15 disagrees with the per-vial regeneration", r["model"], r["group"])
            wx = {"design": "control_only_fit7_cal3_excl_BATcontam_VENACV", "arm": armname[r["model"]], "test_group": r["group"]}
            if len(p5[(p5["design"] == wx["design"]) & (p5["arm"] == wx["arm"]) & (p5["test_group"] == wx["test_group"])]) == 1:
                r["accuracy_excl_flagged"] = P.val(f"fbd_acc_excl_{r['model']}_{r['group']}", "15_time_course/fingerprint_by_duration/part5_by_duration.csv", "accuracy", where=wx,
                                                   note="the consortium-flagged brown-fat-contaminated vena cava vials excluded")
                r["coverage_excl_flagged"] = P.val(f"fbd_cov_excl_{r['model']}_{r['group']}", "15_time_course/fingerprint_by_duration/part5_by_duration.csv", "coverage", where=wx)
        ms = P.read("15_time_course/fingerprint_by_duration/part5_matched_summary.csv")
        ref = {}
        for model, arm in armname.items():
            where = {"arm": arm, "test_group": "control_heldout (pooled over folds)"}
            if len(ms[(ms["arm"] == arm) & (ms["test_group"] == where["test_group"])]) == 1:
                ref[model] = {k: P.val(f"fbd_ref_{k}_{model}", "15_time_course/fingerprint_by_duration/part5_matched_summary.csv", k, where=where)
                              for k in ("accuracy", "coverage", "empty_rate") if k in ms.columns}
                for k in ("accuracy_lo", "accuracy_hi", "coverage_lo", "coverage_hi", "n_test_animals"):
                    if k in ms.columns:
                        ref[model][k] = jsonable(ms[(ms["arm"] == arm) & (ms["test_group"] == where["test_group"])][k].iloc[0])
        S["fbd_ref_coverage_k20"] = ref.get("k20", {}).get("coverage")
        if S["fbd_ref_coverage_k20"] is None:
            S["fbd_ref_coverage_k20"] = P.pending("fbd_ref_coverage_k20", "held-out-control reference row not in part5_matched_summary.csv")["pending"] and None
        fbd_block.update({"reference_held_out_controls": ref,
                          "matched_summary": P.table("fbd_matched_summary", "15_time_course/fingerprint_by_duration/part5_matched_summary.csv", "exercise.json", "fingerprint_by_duration.matched_summary")})
        for name, fn in (("differences", "differences.csv"), ("seeds", "part6_seeds_summary.csv"), ("calibration", "calibration_info.csv")):
            if (p15 / fn).exists():
                fbd_block[name] = P.table(f"fbd_{name}", f"15_time_course/fingerprint_by_duration/{fn}", "exercise.json", f"fingerprint_by_duration.{name}")
    else:
        S["fbd_ref_coverage_k20"] = None
        P.pending("fbd_ref_coverage_k20", "phase 15 (scripts/15_fingerprint_by_duration.py, make time-course) not run")
        fbd_block["reference_held_out_controls"] = None
    # ---- 5. study design dates and physiology (phase 15) ---------------------------------------------------------
    d15 = RES / "15_time_course" / "design"
    ph15 = RES / "15_time_course" / "physiology"
    if (d15 / "design_by_group_sex.csv").exists() and (ph15 / "physiology_group_tests.csv").exists():
        design = {"by_group_sex": P.table("design_by_group_sex", "15_time_course/design/design_by_group_sex.csv", "exercise.json", "design.by_group_sex")}
        if (d15 / "design_contrasts.csv").exists():
            design["contrasts"] = P.table("design_contrasts", "15_time_course/design/design_contrasts.csv", "exercise.json", "design.contrasts")
            dc = P.read("15_time_course/design/design_contrasts.csv")
            for d in ("1w", "2w", "4w", "8w"):
                if len(dc[(dc["duration"] == d) & (dc["sex"] == "pooled")]) == 1:
                    S[f"design_gap_days_{d}"] = P.val(f"design_gap_days_{d}", "15_time_course/design/design_contrasts.csv", "sacrifice_gap_days", where={"duration": d, "sex": "pooled"},
                                                      note="median sacrifice date of the group minus the controls' median, days")
        phys = {"group_tests": P.table("physiology_group_tests", "15_time_course/physiology/physiology_group_tests.csv", "exercise.json", "physiology.group_tests")}
        gt = P.read("15_time_course/physiology/physiology_group_tests.csv")
        for var in ("vo2max_change", "fat_change"):
            for sex in ("female", "male"):
                where = {"variable": var, "contrast": "8w vs control", "sex": sex}
                if len(gt[(gt["variable"] == var) & (gt["contrast"] == where["contrast"]) & (gt["sex"] == sex)]) == 1:
                    for col in ("median_diff", "mw_p", "n_trained", "n_control"):
                        S[f"physio_{var}_8w_{sex}_{col}"] = P.val(f"physio_{var}_8w_{sex}_{col}", "15_time_course/physiology/physiology_group_tests.csv", col, where=where)
    else:
        design = phys = None
        for d in ("1w", "2w", "4w", "8w"):
            S[f"design_gap_days_{d}"] = None
            P.pending(f"design_gap_days_{d}", "phase 15 design script (scripts/15_time_course_design.py, make time-course) not run")
        S["physio_vo2max_change_8w_female_median_diff"] = None
        P.pending("physio_vo2max_change_8w_female_median_diff", "phase 15 physiology (make time-course) not run")
    # ---- 6. the panel genes' training response (phase 05b) ---------------------------------------------------
    f05 = RES / "05_panels" / "TRNSCRPT" / "panel_training_summary.csv"
    if f05.exists():
        ptr = {"genes": P.table("panel_training_response", "05_panels/TRNSCRPT/panel_training_response.csv", "exercise.json", "panel_training.genes"),
               "background": P.table("panel_training_background", "05_panels/TRNSCRPT/panel_training_background.csv", "exercise.json", "panel_training.background"),
               "summary": P.table("panel_training_summary", "05_panels/TRNSCRPT/panel_training_summary.csv", "exercise.json", "panel_training.summary")}
        ps = P.read("05_panels/TRNSCRPT/panel_training_summary.csv")
        for col in ("n_genes", "n_regulated_marker_5pct", "max_abs_logfc_marker", "min_tissue_effect", "median_ratio_tissue_over_training", "min_ratio_tissue_over_training"):
            if col in ps.columns and len(ps[ps["gene_set"] == "k20"]) == 1:
                S[f"ptr_k20_{col}"] = P.val(f"ptr_k20_{col}", "05_panels/TRNSCRPT/panel_training_summary.csv", col, where={"gene_set": "k20"})
    else:
        ptr = None
        for col in ("n_genes", "n_regulated_marker_5pct", "max_abs_logfc_marker", "min_tissue_effect", "median_ratio_tissue_over_training", "min_ratio_tissue_over_training"):
            S[f"ptr_k20_{col}"] = None
            P.pending(f"ptr_k20_{col}", "scripts/05_panel_training_response.py (make panel-training) not run")
    w.write("exercise.json", {"separability": sep, "covariates": cov, "within_tissue": wt, "fingerprint_by_duration": fbd_block,
                              "design": design, "physiology": phys, "panel_training": ptr, "summary": S},
            sorted(P.sources - before))

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
    ("recal3size_gtex_k20", 11.70, 0.005), ("gtex_recal_k20_n3_frac_inf", 0.45, 0.005), # QC-only baseline as balanced accuracy since 2026-09-27 (the brief's 0.873 / 0.949 / 0.975 were plain accuracy; bal_acc_mean in the same file)
    ("qc_technical", 0.874, 0.0005), ("qc_technical_sd", 0.025, 0.0005),
    ("qc_composition", 0.952, 0.0005), ("qc_composition_sd", 0.010, 0.0005), ("qc_all", 0.976, 0.0005), ("qc_all_sd", 0.020, 0.0005),
    # the Exercise page (results/SUMMARY.md and results/07_fusion, 03_eda; per-duration rows from the phase-08 regeneration)
    ("fusion_n_beats_single", 0, 0), ("fusion_n_beats_null", 7, 0), ("fusion_n_tissues", 7, 0), ("verdict_n_training", 4, 0), ("verdict_n_tissues", 7, 0),
    ("taskB_mean_auroc_1w", 0.96, 0.005), ("taskB_mean_auroc_2w", 0.99, 0.005), ("taskB_mean_auroc_4w", 0.94, 0.005), ("taskB_mean_auroc_8w", 0.99, 0.005),
    ("cov_heart_depth_logreg", 0.96, 0.005), ("cov_kidney_library_logreg", 0.89, 0.005), ("cov_skmgn_qc_logreg", 0.96, 0.005), ("cov_skmgn_all_rf", 1.0, 1e-9),
    ("wt_heart_trnscrpt_max_r2_group", 0.70, 0.005), ("wt_TRNSCRPT_n_visible", 13, 0), ("wt_PROT_n_visible", 6, 0), ("wt_METAB_n_visible", 18, 0),
    ("fbd_acc_k20_8w", 0.972, 0.0005), ("fbd_cov_k20_8w", 0.917, 0.0005), ("fbd_acc_k20_1w", 0.944, 0.0005), ("fbd_acc_full_8w", 1.000, 1e-9), ("fbd_cov_full_8w", 0.889, 0.0005),
    ("physio_vo2max_change_8w_female_median_diff", 19.25, 0.005), ("physio_vo2max_change_8w_male_median_diff", 17.12, 0.005), ("design_gap_days_8w", 12, 0.5),
    ("fbd_ref_coverage_k20", 0.939, 0.0005), ("ptr_k20_n_regulated_marker_5pct", 6, 0),
    ("n_plates", 17, 0), ("n_lib_batches", 17, 0), ("n_flowcells", 4, 0), ("shared_genes_bodymap", 21040, 0), ("motrpac_genes", 21193, 0),
    ("orthologs_1to1", 14609, 0), ("orthologs_in_gtex", 14569, 0), ("tile_estimable", 1, 0), ("tile_estimable_total", 171, 0), ("cov_id_full_marginal_one_per_animal", 0.916, 0.0005), ("cov_id_full_marginal_pooled", 0.908, 0.0005),
    ("bridge_sum_ratio_all_genes_pool99", 0.017, 0.002),
    ("bridge_pools_max_sum_ratio_all_genes", 0.053, 0.005),  # hippocampus pool 88 on 3 plates: the largest of the six bridging pools
    ("panel_model_recal_n3_mean_cal_samples", 25.4, 1.0),   # "three animals" in the BodyMap recalibration ≈ 25 calibration samples   # the brief's "~1.7 % of the variance that separates tissues"
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


def _prefix_key(e: dict) -> tuple:
    return (e["file"], json.dumps(e.get("where", {}), sort_keys=True), e["column"], e["agg"])


def _pre_fix_value(e: dict, backup: Path):
    """The value of a provenance entry read from a pre-fix results tree ("" when the file or row is absent)."""
    bp = backup / e["file"].replace("results/", "", 1)
    if not bp.exists():
        return ""
    sel = pd.read_csv(bp)
    for k, v in e["where"].items():
        sel = sel[sel[k].astype(str) == str(v)]
    if e["agg"] == "value" and len(sel) == 1:
        return jsonable(sel[e["column"]].iloc[0])
    if e["agg"] == "mean":
        return jsonable(sel[e["column"]].astype(float).mean())
    if e["agg"] == "std":
        return jsonable(sel[e["column"]].astype(float).std())
    if e["agg"] == "sum":
        return jsonable(sel[e["column"]].astype(float).sum())
    if e["agg"] == "count":
        return len(sel)
    return ""


def write_prefix_table(prov_entries: list[dict], backup: Path, out_csv: Path):
    """Generate the committed pre-fix value table from a pre-quantile-fix results tree (run once; the tree is not in the repo)."""
    rows = []
    for e in prov_entries:
        if e.get("pending") or e.get("agg") == "recomputed" or not e.get("file"):
            continue
        k = _prefix_key(e)
        rows.append({"id": e["id"], "file": k[0], "where": k[1], "column": k[2], "agg": k[3], "pre_fix": _pre_fix_value(e, backup)})
    out_csv.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(rows).to_csv(out_csv, index=False)
    n = sum(1 for r in rows if r["pre_fix"] != "")
    print(f"  wrote {out_csv}: {len(rows)} entries, {n} with a pre-fix value")


def reconciliation(prov_entries: list[dict], out_md: Path, prefix_table: Path = None):
    """docs/NUMBERS_RECONCILIATION.md: every headline number, its source, post-fix value and, where the committed
    pre-fix table holds the same entry, the pre-fix value."""
    prefix_table = prefix_table or PREFIX_TABLE
    pre_by_key: dict[tuple, object] = {}
    if prefix_table.exists():
        pt = pd.read_csv(prefix_table, dtype=str, keep_default_na=False)
        for _, r in pt.iterrows():
            v = r["pre_fix"]
            try:
                v = float(v) if v != "" else ""
            except ValueError:
                pass
            pre_by_key[(r["file"], r["where"], r["column"], r["agg"])] = v
    else:
        print(f"  WARNING: {prefix_table} missing; reconciliation has no pre-fix values")
    rows = []
    for e in prov_entries:
        if e.get("pending") or e.get("agg") == "recomputed" or not e.get("file"):
            rows.append({"id": e["id"], "value": e.get("value"), "source": e.get("file") or ", ".join(e.get("files", [])),
                         "where": json.dumps(e.get("where", {})), "column": e.get("column", ""), "agg": e.get("agg", "pending"), "pre_fix": "", "changed": ""})
            continue
        pre = pre_by_key.get(_prefix_key(e), "")
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
             "`results/` (post-fix, 2026-09-25) is the truth; the pre-fix values come from `docs/reconciliation/pre_quantile_fix_values.csv`, "
             "generated once from the pre-fix results of 2026-09-25 (not in the repository). Values are shown to 4 decimals; the JSON holds them unrounded.", "",
             "## 1. Headline numbers (the home-page tiles and the transfer ladder)", ""]
    head_ids = ["tile_acc_k20", "tile_acc_k20_sd", "acc_k50", "acc_full", "acc_fclassif_k20", "tile_bodymap_k20", "tile_bodymap_cov_k20", "tile_bodymap_empty_k20", "tile_bodymap_recal_k20", "tile_bridge",
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
    lines += ["", "## 3. Recomputed and pending entries", "", "| id | value | files / reason |", "|---|---|---|"]
    for e in prov_entries:
        if e.get("pending"):
            lines.append(f"| `{e['id']}` | pending | {e['reason']} |")
        elif e.get("agg") == "recomputed":
            lines.append(f"| `{e['id']}` | {fmt(e['value']) if not isinstance(e['value'], dict) else 'table'} | {', '.join('`' + f + '`' for f in e['files'])}: {e['note']} |")
    out_md.write_text("\n".join(lines) + "\n")
    print(f"  wrote {out_md} ({len(ch)} changed entries)")
    return {"changed": int(len(ch)), "comparable": int((df["changed"] != "").sum()), "document": "docs/NUMBERS_RECONCILIATION.md"}


def readme_table(prov_entries: list[dict]) -> str:
    byid = {e["id"]: e.get("value") for e in prov_entries}

    def f(i, d=3):
        v = byid.get(i)
        return "pending" if v is None else (f"{v:.{d}f}" if isinstance(v, float) else str(v))
    rows = [
        ("20-gene panel, balanced accuracy (19 tissues, 5 animal-grouped folds)", f"{f('acc_k20')} ± {f('acc_k20_sd')}", "results/05_panels/TRNSCRPT/panel_curve.csv"),
        ("50-gene panel / all genes", f"{f('acc_k50')} / {f('acc_full')}", "results/05_panels/TRNSCRPT/panel_curve.csv, results/04_baselines/TRNSCRPT/summary.csv"),
        ("F-test selector at k = 20 (why the selector matters)", f('acc_fclassif_k20'), "results/05_panels/TRNSCRPT/panel_curve_fclassif.csv"),
        ("Coverage of 90 % sets in-distribution, all genes (pooled / one vial per animal)", f"{f('cov_id_full_marginal_pooled')} / {f('cov_id_full_marginal_one_per_animal')}", "results/06_conformal/TRNSCRPT/coverage.csv"),
        ("Coverage of 90 % sets in-distribution, 20-gene panel (pooled / one vial per animal)", f"{f('cov_id_k20_marginal_pooled')} / {f('cov_id_k20_marginal_one_per_animal')}", "results/31_site_regen/06_conformal/TRNSCRPT/scores_*.csv (recomputed)"),
        ("Trained animals, panel fit on the sedentary controls only: accuracy k20 / coverage", f"{f('acc_train_control_test_trained_k20')} / {f('cov_train_control_test_trained_k20_marginal')}", "results/08_shift/TRNSCRPT/shift_table.csv"),
        ("BodyMap adults (another lab): accuracy k20 / coverage / empty sets", f"{f('acc_bodymap_k20')} / {f('cov_bodymap_k20_marginal')} / {f('empty_bodymap_k20_marginal')}", "results/12_bodymap/"),
        ("BodyMap recalibrated on 3 animals: coverage at set size", f"{f('recal3_bodymap_k20')} at {f('recal3size_bodymap_k20', 2)}", "results/12_bodymap/recalibration.csv"),
        ("GTEx (human): accuracy k20 / k50 / full", f"{f('acc_gtex_k20')} / {f('acc_gtex_k50')} / {f('acc_gtex_full')}", "results/13_gtex/accuracy_overall.csv"),
        ("GTEx coverage k20 / empty; recalibrated on 3 donors (all 20 draws, 9 with no finite threshold): coverage at set size", f"{f('cov_gtex_k20_marginal')} / {f('empty_gtex_k20_marginal')}; {f('recal3_gtex_k20')} at {f('recal3size_gtex_k20', 2)}", "results/13_gtex/"),
        ("Estimable tissue pairs within study (RNA-seq)", f"{f('tile_estimable')} of {f('tile_estimable_total')} ({' and '.join(t.lower() for t in str(byid.get('tile_estimable_pairs', '')).split('|'))})", "results/16_identifiability/estimable_pairs.csv"),
        ("Sedentary vs 8-week-trained within tissue: mean best single-omic AUROC / fusion beats single / attributable to training", f"{f('taskB_mean_auroc_8w')} / {f('fusion_n_beats_single', 0)} of {f('fusion_n_tissues', 0)} / {f('verdict_n_training', 0)} of {f('verdict_n_tissues', 0)}", "results/07_fusion/"),
        ("Trained animals, 8-week group only (cohort-matched with the controls): accuracy k20 / coverage", f"{f('fbd_acc_k20_8w')} / {f('fbd_cov_k20_8w')}", "results/31_site_regen/08_shift_k20/scores_target_vials.csv"),
        ("Batch measured on a bridging reference pool run on 6 plates (Σ V_batch / Σ V_tissue, all genes)", f('bridge_sum_ratio_all_genes_pool99'), "results/16_identifiability/bridge_variance.csv"),
        ("QC covariates alone, balanced accuracy: technical / composition / all", f"{f('qc_technical')} / {f('qc_composition')} / {f('qc_all')}", "results/16_identifiability/qc_only_summary.csv"),
    ]
    return "\n".join(["| result | value | source |", "|---|---|---|"] + [f"| {a} | {b} | `{c}` |" for a, b, c in rows])


def abstract(prov_entries: list[dict]) -> str:
    byid = {e["id"]: e.get("value") for e in prov_entries}
    g = lambda i, d=3: ("pending" if byid.get(i) is None else (f"{byid[i]:.{d}f}" if isinstance(byid[i], float) else str(byid[i])))
    pairs = " and ".join(t.lower() for t in str(byid.get("tile_estimable_pairs", "")).split("|"))
    return (
        f"In the MoTrPAC rat endurance-training transcriptomes (19 tissues, 899 vials, 50 animals), a 20-gene panel selected "
        f"inside each fold by a class-aware round-robin rule reaches {g('acc_k20')} balanced accuracy over animal-grouped folds "
        f"(50 genes {g('acc_k50')}; all genes {g('acc_full')}; a univariate F-test at the same size {g('acc_fclassif_k20')}). "
        f"A split-conformal prediction set promising the true tissue 90 % of the time covers {g('cov_id_k20_marginal_one_per_animal')} of held-out "
        f"animals and, with the panel fit on sedentary controls only, {g('cov_train_control_test_trained_k20_marginal')} of trained animals. "
        f"Beyond the study it abstains rather than errs: {g('cov_bodymap_k20_marginal')} of rat BodyMap organs from another laboratory, "
        f"where the panel names every mapped organ correctly, and {g('cov_gtex_k20_marginal')} of human GTEx samples (accuracy {g('acc_gtex_k20')}). "
        f"Three target animals restore observed coverage within species ({g('recal3_bodymap_k20')}; {g('pv_recal_bodymap_3_min_coverage', 2)}–{g('pv_recal_bodymap_3_max_coverage', 2)} per draw); "
        f"across species five donors restore the number, not the information ({g('pv_recal_gtex_5_mean_coverage_all_draws')} at {g('pv_recal_gtex_5_mean_set_size_finite', 1)} tissues per set). "
        f"As in any multi-tissue design, each tissue was processed as a unit: {g('tile_estimable')} of {g('tile_estimable_total')} tissue pairs "
        f"({pairs}) can be contrasted inside one batch, and library QC numbers alone classify tissue at {g('qc_all')}. "
        f"Two external facts establish the fingerprint as biology: the independent laboratory, and batch measured directly on the consortium's "
        f"bridging reference pools, {g('bridge_sum_ratio_all_genes_pool99')} of the variance that separates tissues."
    )


# ---------------------------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser(description="export site data with provenance")
    ap.add_argument("--results", help="results root to read (default: tfp.config.results_root())")
    ap.add_argument("--out", help="output directory (default: site/data)")
    ap.add_argument("--write-prefix-table", action="store_true", help="write docs/reconciliation/pre_quantile_fix_values.csv from --prefix-backup and exit")
    ap.add_argument("--prefix-backup", help="a pre-quantile-fix results tree (only for --write-prefix-table)")
    ap.add_argument("--prefix-table", help="pre-fix value table for --reconciliation (default docs/reconciliation/pre_quantile_fix_values.csv)")
    ap.add_argument("--skip-expr", action="store_true", help="skip the per-sample expression exports (slow GTEx read)")
    ap.add_argument("--check-anchors", action="store_true")
    ap.add_argument("--reconciliation", action="store_true", help="also write docs/NUMBERS_RECONCILIATION.md")
    ap.add_argument("--readme-table", action="store_true", help="print the README key-results table and exit")
    ap.add_argument("--abstract", action="store_true", help="print the abstract and exit")
    args = ap.parse_args()
    configure(Path(args.results) if args.results else None, Path(args.out) if args.out else None)
    SITE.mkdir(parents=True, exist_ok=True)
    if args.readme_table or args.abstract:
        entries = json.loads((SITE / "provenance.json").read_text())["entries"]
        print(readme_table(entries) if args.readme_table else abstract(entries))
        return
    if args.write_prefix_table:
        assert args.prefix_backup, "--write-prefix-table needs --prefix-backup DIR"
        entries = json.loads((SITE / "provenance.json").read_text())["entries"]
        write_prefix_table(entries, Path(args.prefix_backup), Path(args.prefix_table) if args.prefix_table else PREFIX_TABLE)
        return
    print(f"== results root: {RES} ({'live run' if RES == C.RESULTS_DIR else 'snapshot' if RES == C.FROZEN_DIR else 'explicit'}); output: {SITE}")
    generated = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    ghash = git_hash()
    prov = Prov()
    w = Writer(prov, generated, ghash)
    regen06 = REGEN / "06_conformal" / "TRNSCRPT"
    rec = probs = cal = classes = None
    if (regen06 / "scores_test_probs.csv").exists():
        print("== recompute in-distribution coverage by model from the phase-06 scores")
        rec, probs, cal, classes = coverage_from_scores(regen06)
    print("== aggregates")
    export_aggregates(w, prov, rec)
    print("== headline")
    export_headline(w, prov, rec)
    print("== exercise")
    export_exercise(w, prov)
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
    if (C.RAW_DIR / "pheno.csv").exists():
        print("== genes")
        export_genes(w, prov, args.skip_expr)
    else:
        committed = ROOT / "site" / "data"
        kept = []
        for name in ("genes.json", "expr_motrpac.json", "expr_bodymap.json", "expr_gtex.json"):
            if (committed / name).exists():
                if (SITE / name).resolve() != (committed / name).resolve():
                    (SITE / name).write_bytes((committed / name).read_bytes())
                kept.append(name)
        print(f"== genes: data/raw absent; kept the committed exports ({', '.join(kept)})")
    print("== panel model and card")
    export_panel_model(w, prov)
    export_panel_card(w)
    (SITE / "provenance.json").write_text(json.dumps(jsonable({"_meta": {"generated": generated, "git_hash": ghash, "sources": sorted(prov.sources)},
                                                              "entries": prov.entries, "tables": prov.tables}), indent=0))
    print("== manifest")
    export_manifest(w, prov)
    total = sum(p.stat().st_size for p in SITE.glob("*.json"))
    print(f"== {len(prov.entries)} provenance entries, {len(prov.tables)} tables; site/data total {total / 1e6:.2f} MB; largest "
          f"{max(SITE.glob('*.json'), key=lambda p: p.stat().st_size).name}")
    if args.reconciliation:
        counts = reconciliation(prov.entries, ROOT / "docs" / "NUMBERS_RECONCILIATION.md", Path(args.prefix_table) if args.prefix_table else None)
        pp = SITE / "provenance.json"
        pj = json.loads(pp.read_text())
        pj["_meta"]["reconciliation"] = counts
        pp.write_text(json.dumps(pj, indent=0))
    if args.check_anchors:
        ok = check_anchors(prov.entries)
        if not ok:
            sys.exit(2)


if __name__ == "__main__":
    main()
