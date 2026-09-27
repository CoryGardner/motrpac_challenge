#!/usr/bin/env python
"""Export the numbers the "Check samples" page (site/index.html) shows into site/data/product.json, with a provenance
entry per number appended to site/data/provenance.json (ids prefixed `pv_`; re-running replaces them).

Sources: results_product/40_product/ (scripts/40_product_validation.py) and results_frozen/ (the BodyMap unseen-organ
sets for the "fails safe" trust number). Nothing is typed.
Usage: python scripts/41_export_product_data.py [--site DIR]
"""
from __future__ import annotations

import argparse
import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from tfp import config as C

ROOT = C.ROOT
PV = ROOT / "results_product" / "40_product"
FZ = ROOT / "results_frozen"


def jsonable(x):
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer, int)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return None if math.isnan(v) else v
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if x is None or isinstance(x, str):
        return x
    return None if pd.isna(x) else str(x)


def rel(p):
    return str(Path(p).relative_to(ROOT))


class Prov:
    def __init__(self):
        self.entries, self.tables, self.sources, self.cache = [], [], set(), {}

    def read(self, p):
        p = Path(p)
        if p not in self.cache:
            self.cache[p] = pd.read_csv(p)
        self.sources.add(rel(p))
        return self.cache[p]

    def val(self, id_, p, column, where=None, tol=1e-9):
        sel = self.read(p)
        for k, v in (where or {}).items():
            sel = sel[sel[k].astype(str) == str(v)]
        assert len(sel) == 1, f"{id_}: {where} matched {len(sel)} rows in {p}"
        v = jsonable(sel[column].iloc[0])
        self.entries.append({"id": id_, "value": v, "file": rel(p), "where": where or {}, "column": column, "agg": "value", "tol": tol})
        return v

    def table(self, id_, p, json_path, rows, n_rows=None, quantised=False):
        self.read(p)
        t = {"id": id_, "file": rel(p), "json_file": "product.json", "json_path": json_path, "n_rows": n_rows if n_rows is not None else len(rows)}
        if quantised:
            t["quantised"] = True
        self.tables.append(t)
        return rows


def records(df, cols=None):
    d = df[cols] if cols else df
    return [jsonable(r) for r in d.to_dict("records")]


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--site", default=str(ROOT / "site" / "data"), help="site data directory (CI exports to a scratch copy)")
    args = ap.parse_args()
    site = Path(args.site)
    P = Prov()
    out = {}
    # ---- scaling guard evidence ------------------------------------------------------------------------------------
    sc = PV / "scaling.csv"
    df = P.read(sc)
    out["scaling"] = P.table("pv_scaling", sc, "scaling", records(df, ["subset", "mode", "organ", "n", "mapped", "rat_classes", "accuracy", "coverage", "avg_set_size", "frac_empty"]))
    key = {}
    for sub in ("adult_21wk", "all_ages"):
        for mode in ("within_all", "reference", "within_organ_alone"):
            for col in ("accuracy", "coverage", "frac_empty", "n"):
                key[f"{sub}.{mode}.{col}"] = P.val(f"pv_scaling_{sub}_{mode}_{col}", sc, col, {"subset": sub, "mode": mode, "organ": "ALL_MAPPED"})
        for mode in ("within_all", "reference"):
            for organ in ("Thymus", "Uterus"):
                key[f"{sub}.{mode}.{organ}.frac_empty"] = P.val(f"pv_scaling_{sub}_{mode}_{organ}_frac_empty", sc, "frac_empty", {"subset": sub, "mode": mode, "organ": organ})
    out["scaling_key"] = key
    # ---- composition sensitivity (200 random subsets of the 21-week adults) ------------------------------------------
    cs = PV / "composition_summary.csv"
    comp = {}
    for mode in ("within", "reference"):
        for metric in ("accuracy", "coverage", "frac_empty"):
            for col in ("p05", "p50", "p95"):
                comp[f"{mode}.{metric}.{col}"] = P.val(f"pv_comp_{mode}_{metric}_{col}", cs, col, {"mode": mode, "metric": metric})
    comp["n_subsets"] = P.val("pv_comp_n_subsets", cs, "n_subsets", {"mode": "within", "metric": "coverage"})
    comp["size_min"] = P.val("pv_comp_size_min", cs, "subset_size_min", {"mode": "within", "metric": "coverage"})
    comp["size_max"] = P.val("pv_comp_size_max", cs, "subset_size_max", {"mode": "within", "metric": "coverage"})
    out["composition"] = comp
    # ---- per-draw recalibration coverage (replayed from results_frozen) ----------------------------------------------------
    rds = PV / "recal_draws_summary.csv"
    rd = {}
    for ds in ("bodymap", "gtex"):
        for n in (3, 5):
            w = {"dataset": ds, "model": "k20", "n_recal": n}
            rd[f"{ds}.{n}"] = {c: P.val(f"pv_recal_{ds}_{n}_{c}", rds, c, w) for c in
                              ("n_recal", "draws", "mean_coverage_all_draws", "min_coverage", "max_coverage", "n_draws_below_0.90", "n_finite", "n_infinite",
                               "mean_coverage_finite", "min_coverage_finite", "mean_set_size_finite", "min_cal_scores", "max_cal_scores", "n_classes")}
    out["recal_draws"] = rd
    # ---- the flag ---------------------------------------------------------------------------------------------------
    fr = PV / "flag_rates.csv"
    fdf = P.read(fr)
    out["flag"] = P.table("pv_flag_rates", fr, "flag", records(fdf))
    fk = {}
    for setting in fdf["setting"].unique():
        for a in (0.05, 0.1, 0.2):
            for col in ("false_mismatch_rate", "cant_confirm_rate", "swap_vial_detection_rate", "swap_pair_detection_rate", "swap_vial_missed_rate", "n_samples", "n_swaps"):
                fk[f"{setting}.{a}.{col}"] = P.val(f"pv_flag_{setting}_{a}_{col}", fr, col, {"setting": setting, "alpha": a})
    out["flag_key"] = fk
    # ---- vena cava → brown fat -------------------------------------------------------------------------------------
    vc = PV / "venacv_cases.csv"
    vdf = P.read(vc)
    out["venacv"] = P.table("pv_venacv_cases", vc, "venacv", records(vdf.astype({"viallabel": str})))
    vs = PV / "venacv_summary.csv"
    out["venacv_all"] = {("flagged" if fl else "unflagged"): {c: P.val(f"pv_venacv_{'flagged' if fl else 'unflagged'}_{c}", vs, c, {"consortium_flagged": fl})
                                                             for c in ("n_vials", "n_called_bat", "n_consistent", "n_mismatch", "n_cant_confirm")} for fl in (True, False)}
    # ---- reference map ---------------------------------------------------------------------------------------------
    ld = P.read(PV / "pca_loadings.csv")
    rc = P.read(PV / "pca_reference_coords.csv")
    ps = PV / "pca_summary.csv"
    out["pca"] = {"genes": list(ld["feature_ID"]), "centre": [float(v) for v in ld["centre"]], "loadings": [[float(v) for v in ld["PC1"]], [float(v) for v in ld["PC2"]]],
                  "explained": [P.val("pv_pca_explained_PC1", ps, "explained", {"PC": "PC1"}), P.val("pv_pca_explained_PC2", ps, "explained", {"PC": "PC2"})],
                  "reference": {"id": [str(v) for v in rc["viallabel"]], "tissue": list(rc["tissue"]), "x": [round(float(v), 6) for v in rc["PC1"]], "y": [round(float(v), 6) for v in rc["PC2"]]}}
    P.tables.append({"id": "pv_pca_loadings", "file": rel(PV / "pca_loadings.csv"), "json_file": "product.json", "json_path": "pca.loadings", "n_rows": int(len(ld)), "matrix": True})
    P.tables.append({"id": "pv_pca_reference", "file": rel(PV / "pca_reference_coords.csv"), "json_file": "product.json", "json_path": "pca.reference", "n_rows": int(len(rc)), "matrix": True})
    # ---- trust strip: unseen organs in another lab (k20, marginal) ----------------------------------------------------------------------
    ood = FZ / "12_bodymap" / "ood_sets.csv"
    tr = {}
    for organ in ("Thymus", "Uterus"):
        tr[organ] = {"n": P.val(f"pv_ood_{organ}_n", ood, "n", {"model": "k20", "conformal": "marginal", "organ": organ}),
                     "frac_empty": P.val(f"pv_ood_{organ}_frac_empty", ood, "frac_empty", {"model": "k20", "conformal": "marginal", "organ": organ})}
    out["ood"] = tr
    # ---- write -----------------------------------------------------------------------------------------------------
    try:
        gh = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        gh = "n/a"
    out["_meta"] = {"generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), "git_hash": gh, "sources": sorted(P.sources),
                    "note": "every value traces to results_product/40_product/ or results_frozen/ through the pv_* entries of provenance.json"}
    (site / "product.json").write_text(json.dumps(out, indent=0))
    pp = site / "provenance.json"
    prov = json.loads(pp.read_text())
    prov["entries"] = [e for e in prov["entries"] if not str(e.get("id", "")).startswith("pv_")] + P.entries
    prov["tables"] = [t for t in prov.get("tables", []) if not str(t.get("id", "")).startswith("pv_")] + P.tables
    prov.setdefault("_meta", {})["product"] = {"generated": out["_meta"]["generated"], "git_hash": gh, "n_entries": len(P.entries), "n_tables": len(P.tables), "sources": sorted(P.sources)}
    pp.write_text(json.dumps(prov, indent=0))
    print(f"wrote {site / 'product.json'} ({(site / 'product.json').stat().st_size / 1e3:.0f} kB); {len(P.entries)} provenance entries, {len(P.tables)} tables")


if __name__ == "__main__":
    main()
