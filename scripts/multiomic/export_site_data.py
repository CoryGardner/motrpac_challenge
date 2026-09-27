#!/usr/bin/env python
"""Export the numbers the multiomic page shows into site/data/multiomic.json, with a provenance entry per number appended
to site/data/provenance.json (ids prefixed `mo_`; re-running replaces the previous `mo_` entries).

Every value is read from a result file under results_multiomic/ or results_frozen/ (the RNA comparison rows) and the
provenance entry records file, row selector and column, in the format tests/test_site_data.py checks. Nothing is typed.
Usage: python scripts/multiomic/export_site_data.py
"""
from __future__ import annotations

import json
import math
import subprocess
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

from tfp import config as C

ROOT = C.ROOT
R = ROOT / "results_multiomic"
FZ = ROOT / "results_frozen"
SITE = ROOT / "site" / "data"


def jsonable(x):
    if isinstance(x, (np.bool_, bool)):
        return bool(x)
    if isinstance(x, (np.integer, int)):
        return int(x)
    if isinstance(x, (np.floating, float)):
        v = float(x)
        return None if math.isnan(v) else ("inf" if math.isinf(v) else v)
    if isinstance(x, (list, tuple, np.ndarray)):
        return [jsonable(v) for v in x]
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if x is None or isinstance(x, str):
        return x
    if pd.isna(x):
        return None
    return str(x)


def rel(p: Path) -> str:
    return str(Path(p).relative_to(ROOT))


class Prov:
    def __init__(self):
        self.entries, self.tables, self.sources, self.cache = [], [], set(), {}

    def read(self, p: Path) -> pd.DataFrame:
        p = Path(p)
        if p not in self.cache:
            self.cache[p] = pd.read_csv(p)
        self.sources.add(rel(p))
        return self.cache[p]

    def val(self, id_, p: Path, column: str, where: dict | None = None, tol: float = 1e-6):
        df = self.read(p)
        sel = df
        for k, v in (where or {}).items():
            sel = sel[sel[k].astype(str) == str(v)]
        assert len(sel) == 1, f"{id_}: selector {where} matched {len(sel)} rows in {p}"
        v = sel[column].iloc[0]
        v = jsonable(v)
        entry = {"id": id_, "value": v, "file": rel(p), "where": where or {}, "column": column, "agg": "value", "tol": tol}
        if v is None:   # an empty cell (e.g. no genes shared): recorded as pending so the ledger test reads the reason instead of comparing NaN to null
            entry.update({"pending": True, "reason": f"empty cell in {rel(p)} ({column})"})
        self.entries.append(entry)
        return v

    def table(self, id_, p: Path, json_path: str, rows: list[dict]):
        self.sources.add(rel(p))
        self.tables.append({"id": id_, "file": rel(p), "json_file": "multiomic.json", "json_path": json_path, "n_rows": len(rows)})
        return rows


def records(df: pd.DataFrame, cols: list[str] | None = None) -> list[dict]:
    d = df[cols] if cols else df
    return [jsonable(r) for r in d.to_dict("records")]


def main():
    prov = Prov()
    out: dict = {}
    # ---- scales (Phase 1) ---------------------------------------------------------------------------------------------
    vp, vpr, vpc = R / "01_rii" / "variance_partition.csv", R / "01_rii" / "variance_partition_ratio.csv", R / "01_rii" / "variance_partition_complete.csv"
    js = R / "01_rii" / "join_summary.csv"
    out["scales"] = {"rows": [], "n_vials": prov.val("mo_rii_n_vials", js, "n_vials", {"assay": "prot-pr"}), "n_animals": prov.val("mo_rii_n_animals", js, "n_animals", {"assay": "prot-pr"}),
                     "n_proteins_inner": prov.val("mo_rii_n_proteins_inner", js, "n_proteins_inner", {"assay": "prot-pr"}), "n_proteins_complete": prov.val("mo_rii_n_proteins_complete", js, "n_proteins_inner_complete", {"assay": "prot-pr"}),
                     "n_proteins_union": prov.val("mo_rii_n_proteins_union", js, "n_proteins_union", {"assay": "prot-pr"})}
    for pc in ("PC1", "PC2", "PC3"):
        out["scales"]["rows"].append({"PC": pc, "r2_rii": prov.val(f"mo_r2_rii_{pc}", vp, "R2_tissue", {"PC": pc}), "r2_rii_null95": prov.val(f"mo_r2_rii_null95_{pc}", vp, "R2_tissue_null95", {"PC": pc}),
                                      "r2_ratio": prov.val(f"mo_r2_ratio_{pc}", vpr, "R2_tissue", {"PC": pc}), "r2_complete": prov.val(f"mo_r2_complete_{pc}", vpc, "R2_tissue", {"PC": pc}),
                                      "explained_rii": prov.val(f"mo_explained_rii_{pc}", vp, "explained", {"PC": pc}), "explained_ratio": prov.val(f"mo_explained_ratio_{pc}", vpr, "explained", {"PC": pc})})
    ds = R / "01_rii" / "diagnostic_accuracy_summary.csv"
    out["scales"]["missingness_outer_acc"] = prov.val("mo_missingness_outer_acc", ds, "mean", {"quantity": "missingness_outer_acc"})
    out["scales"]["means_removed_acc"] = prov.val("mo_means_removed_acc", ds, "mean", {"quantity": "per_tissue_means_removed_acc"})
    pcs = R / "01_rii" / "panel_curve_summary.csv"
    out["scales"]["panel_k20_bal_acc"] = prov.val("mo_rii_panel_k20_bal_acc", pcs, "bal_acc_mean", {"k": 20})
    out["scales"]["panel_k20_null95"] = prov.val("mo_rii_panel_k20_null95", R / "01_rii" / "panel_null_k20.csv", "null_q95", {"k": 20})
    st = R / "01_rii" / "stability_k20_summary.csv"
    out["scales"]["stability_ge_0_8"] = prov.val("mo_rii_stability_ge08", st, "n_selected_ge_0.8", {"k": 20})
    # ---- ladder: protein beside RNA ---------------------------------------------------------------------------------------
    side = prov.read(R / "09_extensions" / "ladder_side_by_side.csv")
    ladder = []
    for r in side.itertuples():
        ladder.append({"rung_id": r.split, "layer": "protein" if str(r.layer).startswith("protein") else "RNA", "model": str(r.arm), "accuracy": jsonable(r.accuracy), "accuracy_sd": jsonable(getattr(r, "accuracy_sd", np.nan)),
                       "coverage": jsonable(r.coverage), "coverage_sd": jsonable(getattr(r, "coverage_sd", np.nan)), "set_size": jsonable(r.set_size), "empty": jsonable(r.frac_empty), "n_samples": jsonable(r.n_test_samples),
                       "n_individuals": jsonable(r.n_test_animals), "n_classes": jsonable(r.n_classes), "source": str(r.source), "where": "MoTrPAC"})
    prov.table("mo_ladder_within", R / "09_extensions" / "ladder_side_by_side.csv", "ladder_within", ladder)
    # other species: protein → Jiang (cleaned relative and raw ppm) with donor-bootstrap intervals; RNA → GTEx (frozen); RNA 7-class on the same 42 Jiang samples (Phase 5)
    xs = []
    for tag, d, label in (("jiang_relative", R / "03_prot_transfer", "protein → Jiang 2020 (cleaned relative)"), ("jiang_rawppm", R / "03_prot_transfer" / "rawppm", "protein → Jiang 2020 (raw ppm, same processing both sides)")):
        for m in ("k20", "k50", "full"):
            cov_id = f"mo_{tag}_{m}_coverage"
            xs.append({"rung_id": "different_species", "target": label, "layer": "protein", "model": m,
                       "accuracy": prov.val(f"mo_{tag}_{m}_accuracy", d / "accuracy_overall.csv", "accuracy_sample_weighted", {"model": m}),
                       "accuracy_ci": [prov.val(f"mo_{tag}_{m}_acc_lo", d / "accuracy_overall.csv", "acc_ci95_low_donor_boot", {"model": m}), prov.val(f"mo_{tag}_{m}_acc_hi", d / "accuracy_overall.csv", "acc_ci95_high_donor_boot", {"model": m})],
                       "n_samples": prov.val(f"mo_{tag}_{m}_n", d / "accuracy_overall.csv", "n_samples_mapped", {"model": m}), "n_individuals": prov.val(f"mo_{tag}_{m}_donors", d / "accuracy_overall.csv", "n_donors_mapped", {"model": m}),
                       "coverage": prov.val(cov_id, d / "conformal_transfer.csv", "coverage_mapped", {"model": m, "conformal": "marginal"}),
                       "coverage_ci": [prov.val(f"mo_{tag}_{m}_cov_lo", d / "coverage_ci.csv", "ci95_low_donor_boot", {"model": m}), prov.val(f"mo_{tag}_{m}_cov_hi", d / "coverage_ci.csv", "ci95_high_donor_boot", {"model": m})],
                       "empty": prov.val(f"mo_{tag}_{m}_empty", d / "conformal_transfer.csv", "frac_empty_mapped", {"model": m, "conformal": "marginal"}),
                       "set_size": prov.val(f"mo_{tag}_{m}_setsize", d / "conformal_transfer.csv", "avg_set_size_mapped", {"model": m, "conformal": "marginal"}),
                       "ood_empty": prov.val(f"mo_{tag}_{m}_ood_empty", d / "conformal_transfer.csv", "ood_frac_empty", {"model": m, "conformal": "marginal"}), "n_classes": 7, "source": rel(d)})
    g = FZ / "13_gtex"
    for m in ("k20", "k50", "full"):
        xs.append({"rung_id": "different_species", "target": "RNA → GTEx v8 (frozen phase 13)", "layer": "RNA", "model": m,
                   "accuracy": prov.val(f"mo_gtex_{m}_accuracy", g / "accuracy_overall.csv", "accuracy_sample_weighted", {"model": m}), "accuracy_ci": None,
                   "n_samples": prov.val(f"mo_gtex_{m}_n", g / "accuracy_overall.csv", "n_samples", {"model": m}), "n_individuals": 862,
                   "coverage": prov.val(f"mo_gtex_{m}_coverage", g / "conformal_transfer.csv", "coverage_mapped", {"model": m, "conformal": "marginal"}), "coverage_ci": None,
                   "empty": prov.val(f"mo_gtex_{m}_empty", g / "conformal_transfer.csv", "frac_empty_mapped", {"model": m, "conformal": "marginal"}),
                   "set_size": prov.val(f"mo_gtex_{m}_setsize", g / "conformal_transfer.csv", "avg_set_size_mapped", {"model": m, "conformal": "marginal"}), "ood_empty": None, "n_classes": 19, "source": rel(g)})
    fs = R / "05_fusion_transfer" / "fusion_transfer_summary.csv"
    for layer in ("RNA", "protein"):
        for m in ("k20", "full"):
            xs.append({"rung_id": "different_species", "target": f"{layer}, 7-class fingerprint on the 42 Jiang samples with both layers (Phase 5)", "layer": layer, "model": m,
                       "accuracy": prov.val(f"mo_same42_{layer}_{m}_accuracy", fs, "accuracy", {"model": m, "layer": layer}),
                       "accuracy_ci": [prov.val(f"mo_same42_{layer}_{m}_acc_lo", fs, "acc_ci95_low_donor_boot", {"model": m, "layer": layer}), prov.val(f"mo_same42_{layer}_{m}_acc_hi", fs, "acc_ci95_high_donor_boot", {"model": m, "layer": layer})],
                       "n_samples": prov.val(f"mo_same42_{layer}_{m}_n", fs, "n_mapped", {"model": m, "layer": layer}), "n_individuals": prov.val(f"mo_same42_{layer}_{m}_donors", fs, "n_donors", {"model": m, "layer": layer}),
                       "coverage": prov.val(f"mo_same42_{layer}_{m}_coverage", fs, "coverage_motrpac_cal", {"model": m, "layer": layer}),
                       "coverage_ci": [prov.val(f"mo_same42_{layer}_{m}_cov_lo", fs, "cov_ci95_low_donor_boot", {"model": m, "layer": layer}), prov.val(f"mo_same42_{layer}_{m}_cov_hi", fs, "cov_ci95_high_donor_boot", {"model": m, "layer": layer})],
                       "empty": prov.val(f"mo_same42_{layer}_{m}_empty", fs, "frac_empty", {"model": m, "layer": layer}), "set_size": prov.val(f"mo_same42_{layer}_{m}_setsize", fs, "avg_set_size", {"model": m, "layer": layer}),
                       "ood_empty": None, "n_classes": 7, "source": rel(fs)})
    out["ladder_within"] = ladder
    out["ladder_species"] = xs
    # ---- recalibration: coverage AND set size vs donors, with the label-space size ------------------------------------------------
    rs = prov.read(R / "08_verification" / "recalibration_set_sizes.csv")
    rec = records(rs, ["source", "model", "n_recal", "n_test_individuals", "draws", "coverage_recalibrated", "coverage_source_cal_same_test", "set_size_recalibrated", "n_classes_label_space", "set_size_frac_of_classes", "frac_empty_recalibrated", "file"])
    out["recalibration"] = prov.table("mo_recalibration", R / "08_verification" / "recalibration_set_sizes.csv", "recalibration", rec)
    # the n = 0 points (source calibration) for the same runs
    n0 = []
    for label, d, ncls in (("Phase 3 protein → Jiang 2020 (cleaned relative)", R / "03_prot_transfer", 7), ("Phase 3 protein → Jiang 2020 (raw ppm)", R / "03_prot_transfer" / "rawppm", 7),
                           ("Phase 4 hilic_sato", R / "04_metab_transfer" / "hilic_sato", 19), ("Phase 4 deep_sato", R / "04_metab_transfer" / "deep_sato", 9), ("Phase 4 deep_mw", R / "04_metab_transfer" / "deep_mw", 9),
                           ("frozen phase 13 RNA → GTEx", FZ / "13_gtex", 19)):
        ct = prov.read(d / "conformal_transfer.csv")
        primary = {"Phase 4 hilic_sato": "Sedentary", "Phase 4 deep_sato": "Sedentary", "Phase 4 deep_mw": "3"}.get(label, "adult")
        for m in ("k20", "full"):
            sel = ct[(ct["model"] == m) & (ct["conformal"] == "marginal") & (ct["stage"].astype(str) == primary)]
            if len(sel) == 1:
                n0.append({"source": label, "model": m, "n_recal": 0, "coverage": jsonable(sel["coverage_mapped"].iloc[0]), "set_size": jsonable(sel["avg_set_size_mapped"].iloc[0]), "n_classes_label_space": ncls, "stage": primary, "file": rel(d / "conformal_transfer.csv")})
    out["recalibration_n0"] = n0
    # ---- fusion under shift ------------------------------------------------------------------------------------------------------
    fu = prov.read(fs)
    out["fusion"] = prov.table("mo_fusion", fs, "fusion", records(fu, ["model", "layer", "n_mapped", "n_donors", "accuracy", "acc_ci95_low_donor_boot", "acc_ci95_high_donor_boot", "coverage_motrpac_cal", "cov_ci95_low_donor_boot", "cov_ci95_high_donor_boot", "frac_empty", "avg_set_size", "ood_frac_empty"]))
    cs = prov.read(R / "05_fusion_transfer" / "confusion_structure.csv")
    out["confusion_structure"] = prov.table("mo_confusion_structure", R / "05_fusion_transfer" / "confusion_structure.csv", "confusion_structure", records(cs))
    # ---- design comparison ---------------------------------------------------------------------------------------------------------
    dc = prov.read(R / "06_external_identifiability" / "design_comparison.csv")
    out["design"] = prov.table("mo_design", R / "06_external_identifiability" / "design_comparison.csv", "design", records(dc, ["dataset", "batch_variable", "n_samples", "n_tissues", "n_levels", "max_tissues_per_level", "cramers_v", "n_pairs_estimable", "n_pairs_total", "frac_pairs_estimable", "source"]))
    out["design_jiang"] = {"cramers_v": prov.val("mo_jiang_cramers_v", R / "06_external_identifiability" / "design_comparison.csv", "cramers_v", {"dataset": "Jiang2020"}),
                           "n_pairs_estimable": prov.val("mo_jiang_pairs_est", R / "06_external_identifiability" / "design_comparison.csv", "n_pairs_estimable", {"dataset": "Jiang2020"}),
                           "n_pairs_total": prov.val("mo_jiang_pairs_total", R / "06_external_identifiability" / "design_comparison.csv", "n_pairs_total", {"dataset": "Jiang2020"}),
                           "max_tissues_per_run": prov.val("mo_jiang_max_tissues_per_run", R / "06_external_identifiability" / "design_comparison.csv", "max_tissues_per_level", {"dataset": "Jiang2020"}),
                           "n_runs": prov.val("mo_jiang_n_runs", R / "06_external_identifiability" / "design_comparison.csv", "n_levels", {"dataset": "Jiang2020"}),
                           "motrpac_prot_cramers_v": prov.val("mo_prot_cramers_v", R / "06_external_identifiability" / "design_comparison.csv", "cramers_v", {"dataset": "MoTrPAC PROT"}),
                           "motrpac_prot_pairs_total": prov.val("mo_prot_pairs_total", R / "06_external_identifiability" / "design_comparison.csv", "n_pairs_total", {"dataset": "MoTrPAC PROT"})}
    # ---- RNA markers at the protein level -------------------------------------------------------------------------------------------
    rp = prov.read(R / "01_rii" / "rna_protein_correlation.csv")
    edges = np.linspace(-1, 1, 21)
    counts, _ = np.histogram(rp["spearman"].dropna().to_numpy(), bins=edges)
    rps = R / "01_rii" / "rna_protein_correlation_summary.csv"
    out["rna_protein"] = {"hist_edges": [float(e) for e in edges], "hist_counts": [int(c) for c in counts], "n_genes": prov.val("mo_rp_n_genes", rps, "n_genes"), "spearman_median": prov.val("mo_rp_spearman_median", rps, "spearman_median"),
                          "null_median": prov.val("mo_rp_null_median", rps, "null_mismatched_median"), "null_q95": prov.val("mo_rp_null_q95", rps, "null_mismatched_q95"), "frac_above_null_q95": prov.val("mo_rp_frac_above_q95", rps, "frac_genes_above_null_q95"),
                          "frac_same_marker": prov.val("mo_rp_frac_same_marker", rps, "frac_same_marker_tissue")}
    prov.tables.append({"id": "mo_rp_hist", "file": rel(R / "01_rii" / "rna_protein_correlation.csv"), "json_file": "multiomic.json", "json_path": "rna_protein.hist_counts", "n_rows": int(len(rp)), "quantised": True})
    pg = prov.read(R / "01_rii" / "rna_protein_panel_genes.csv")
    pgt = pg[pg["testable_c"] == True]
    out["rna_protein"]["panel_genes"] = prov.table("mo_rp_panel_genes", R / "01_rii" / "rna_protein_panel_genes.csv", "rna_protein.panel_genes",
                                                    [{"gene_symbol": r.gene_symbol, "lists": r.lists, "rna_marker_tissue": r.marker_tissue_19, "protein_marker_tissue": r.prot_marker_tissue, "agree": bool(r.protein_marker_equals_rna_marker19)} for r in pgt.itertuples()])
    prov.tables[-1]["n_rows"] = int(len(pg)); prov.tables[-1]["quantised"] = True   # subset of the CSV (testable rows)
    pss = R / "01_rii" / "rna_protein_panel_summary.csv"
    out["rna_protein"]["n_testable"] = prov.val("mo_rp_n_testable", pss, "n_testable_c"); out["rna_protein"]["n_same_marker"] = prov.val("mo_rp_n_same_marker", pss, "n_same_marker")
    # ---- per-tissue verification, Sato design, B1 overlap ---------------------------------------------------------------------------------------
    pt = prov.read(R / "08_verification" / "per_tissue_k20.csv")
    out["per_tissue_k20"] = prov.table("mo_per_tissue_k20", R / "08_verification" / "per_tissue_k20.csv", "per_tissue_k20", records(pt, ["run", "layer", "jiang_tissue", "rat_class", "n_samples", "n_donors", "accuracy", "top_prediction"]))
    sv = R / "08_verification" / "sato_invariance_design.csv"
    out["sato"] = {k: prov.val(f"mo_sato_{k}", sv, k) for k in ("n_sedentary_mice_union", "n_exercised_mice_union", "n_sedentary_samples_total", "n_exercised_samples_total", "invariance_fit_animals", "invariance_cal_animals",
                                                                 "invariance_test_exercised_samples", "invariance_test_exercised_animals", "invariance_k20_accuracy", "invariance_k20_coverage", "invariance_k20_set_size",
                                                                 "invariance_full_accuracy", "invariance_full_coverage", "invariance_full_set_size", "invariance_n_classes")}
    out["sato"]["exercise_design"] = prov.val("mo_sato_exercise_design", sv, "exercise_design")
    L4 = R / "04_metab_transfer" / "legs_summary.csv"
    out["metabolites"] = []
    for leg in ("hilic_sato", "deep_sato", "deep_mw"):
        out["metabolites"].append({"leg": leg, "matched": prov.val(f"mo_{leg}_matched", L4, "matched_refmet", {"leg": leg}), "n_mapped": prov.val(f"mo_{leg}_n_mapped", L4, "n_mapped", {"leg": leg}),
                                   "n_individuals": prov.val(f"mo_{leg}_n_ind", L4, "n_individuals_mapped", {"leg": leg}), "acc_k20": prov.val(f"mo_{leg}_acc_k20", L4, "acc_k20", {"leg": leg}),
                                   "acc_k20_ci": [prov.val(f"mo_{leg}_acc_k20_lo", L4, "acc_k0_ci_low", {"leg": leg}), prov.val(f"mo_{leg}_acc_k20_hi", L4, "acc_k0_ci_high", {"leg": leg})],
                                   "acc_full": prov.val(f"mo_{leg}_acc_full", L4, "acc_full", {"leg": leg}), "chance": prov.val(f"mo_{leg}_chance", L4, "chance", {"leg": leg}),
                                   "coverage_k20": prov.val(f"mo_{leg}_cov_k20", L4, "coverage_k0_source_cal_primary", {"leg": leg}), "coverage_k20_recal5": prov.val(f"mo_{leg}_cov_k20_recal5", L4, "coverage_k0_recal5", {"leg": leg}),
                                   "set_size_k20_recal5": prov.val(f"mo_{leg}_size_k20_recal5", L4, "set_size_k0_recal5", {"leg": leg}), "n_classes": 19 if leg.startswith("hilic") else 9, "source_tissues": prov.val(f"mo_{leg}_src_tissues", L4, "source_tissues", {"leg": leg})})
    out["metabolites_stopped"] = {"leg": "hilic_mw", "matched": prov.val("mo_hilic_mw_matched", L4, "matched_refmet", {"leg": "hilic_mw"}), "min_overlap": prov.val("mo_hilic_mw_min_overlap", L4, "min_overlap", {"leg": "hilic_mw"})}
    out["metabolites_stopped"]["source_metabolites"] = prov.val("mo_hilic_mw_src_metabolites", L4, "source_metabolites", {"leg": "hilic_mw"})
    # ---- v9: "how it was done" (design counts, matrix build, the 20-protein transfer panel) --------------------------------------
    ts = R / "01_rii" / "tissue_summary.csv"
    dcp = R / "06_external_identifiability" / "design_comparison.csv"
    out["how"] = {"n_tissues_prot": prov.val("mo_how_n_tissues_prot", js, "n_tissues", {"assay": "prot-pr"}),
                  "n_tissues_rna": prov.val("mo_how_n_tissues_rna", dcp, "n_tissues", {"dataset": "MoTrPAC TRNSCRPT"}),
                  "min_peptides": prov.val("mo_how_min_peptides", js, "min_peptides", {"assay": "prot-pr"}),
                  "norm": prov.val("mo_how_norm", js, "norm", {"assay": "prot-pr"}),
                  "release": prov.val("mo_how_release", js, "release", {"assay": "prot-pr"})}
    tsd = prov.read(ts)
    out["how"]["prot_tissues"] = prov.table("mo_how_prot_tissues", ts, "how.prot_tissues", records(tsd, ["tissue", "n_vials", "n_animals", "n_plexes"]))
    pgc = prov.read(R / "03_prot_transfer" / "panel_gene_check.csv")
    panel = [{"protein": r.gene_symbol, "marker_tissue": r.marker_tissue, "direction": "higher" if r.source_effect_z > 0 else "lower",
              "effect_rat_z": jsonable(r.source_effect_z), "effect_human_z": jsonable(r.target_effect_z)} for r in pgc.itertuples()]
    out["how"]["transfer_panel"] = prov.table("mo_how_transfer_panel", R / "03_prot_transfer" / "panel_gene_check.csv", "how.transfer_panel", panel)
    # ---- v9: metabolites — within-MoTrPAC tissue R² of the first metabolite PC, deep_mw per tissue at k20 ------------------------------
    vpm = FZ / "03_eda" / "variance_partition_METAB.csv"
    out["metab_within"] = {"r2_tissue_pc1": prov.val("mo_metab_r2_tissue_PC1", vpm, "R2_tissue", {"PC": "PC1"}),
                           "explained_pc1": prov.val("mo_metab_explained_PC1", vpm, "explained", {"PC": "PC1"})}
    abt = R / "04_metab_transfer" / "deep_mw" / "accuracy_by_tissue.csv"
    ab = prov.read(abt)
    out["deep_mw_by_tissue_k20"] = prov.table("mo_deep_mw_by_tissue_k20", abt, "deep_mw_by_tissue_k20",
                                              records(ab[ab["model"] == "k20"], ["target_tissue", "rat_classes", "n", "n_individuals", "accuracy", "top_prediction", "top_prediction_frac"]))
    prov.tables[-1]["n_rows"] = int(len(ab)); prov.tables[-1]["quantised"] = True   # the k20 rows of the CSV
    ovs = R / "09_extensions" / "overlap_summary.csv"
    out["b1"] = {k: prov.val(f"mo_b1_{k}", ovs, k) for k in ("n_rna_k20", "n_protein_k20", "n_protein_stable_core_ge_0.8", "k20_intersection", "rna_k20_x_protein_core", "genes_in_both_candidate_lists", "genes_in_both_with_rna_marker_in_prot7", "n_markers_agree", "prereg_B1_pass")}
    out["b1"]["genes_in_both"] = prov.val("mo_b1_genes_in_both", ovs, "genes_in_both")
    # ---- write --------------------------------------------------------------------------------------------------------
    try:
        gh = subprocess.check_output(["git", "rev-parse", "--short", "HEAD"], cwd=ROOT, text=True).strip()
    except Exception:
        gh = "n/a"
    out["_meta"] = {"generated": datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"), "git_hash": gh, "branch": "main", "sources": sorted(prov.sources),
                    "note": "follow-up work, merged into main (developed on branch multiomic-overnight); every value traces to a results_multiomic/ or results_frozen/ file through the mo_* entries of provenance.json"}
    (SITE / "multiomic.json").write_text(json.dumps(out, indent=0))
    P = json.loads((SITE / "provenance.json").read_text())
    P["entries"] = [e for e in P["entries"] if not str(e.get("id", "")).startswith("mo_")] + prov.entries
    P["tables"] = [t for t in P.get("tables", []) if not str(t.get("id", "")).startswith("mo_")] + prov.tables
    P.setdefault("_meta", {})["multiomic"] = {"generated": out["_meta"]["generated"], "git_hash": gh, "n_entries": len(prov.entries), "n_tables": len(prov.tables), "sources": sorted(prov.sources)}
    (SITE / "provenance.json").write_text(json.dumps(P, indent=0))
    print(f"wrote {SITE / 'multiomic.json'} ({(SITE / 'multiomic.json').stat().st_size / 1e3:.0f} kB); {len(prov.entries)} provenance entries, {len(prov.tables)} tables, {len(prov.sources)} source files")


if __name__ == "__main__":
    main()
