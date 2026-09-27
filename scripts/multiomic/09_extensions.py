#!/usr/bin/env python
"""Run 2, part B — two cheap extensions on local data (tfp code unchanged; pre-registered in docs/PREREGISTRATION_MULTIOMIC.md).

B1  Overlap between the RNA k20 panel (19 tissues; results_frozen/12_bodymap/panel_gene_check.csv = the all-animal k20 panel,
    results_frozen/05_panels/TRNSCRPT/stability_k20_annotated.csv = the 51-gene stability list) and the protein k20 panel
    (results_multiomic/01_rii/panel_k20_all_animals.csv) plus the protein stable core (stability_k20.csv, ≥ 0.80 over 50
    animal-bootstraps): genes chosen in both layers, marker-tissue agreement.
B2  The protein ladder inside MoTrPAC on the RII scale, through the phase-08 code path (scripts/08_shift_tests.py, imported):
    held-out animals (5-fold, animal-grouped, calibration on 30 % of the training animals), trained animals with the panel fit
    on sedentary controls only, held-out sex (both directions); k = 20 and the all-protein model; accuracy, coverage, set size,
    recalibration on 3 and 5 target animals; beside the RNA rungs from results_frozen/08_shift/TRNSCRPT/shift_table.csv and the
    RNA in-distribution numbers.
Outputs: results_multiomic/09_extensions/.
"""
from __future__ import annotations

import importlib
import json
import sys
import time
import warnings

import numpy as np
import pandas as pd

from tfp import config as C, conformal as cp, io, models, report
from tfp.splits import fit_calibration_split, grouped_kfold, leave_one_sex_out, train_controls_test_trained

ROOT = C.ROOT
R = ROOT / "results_multiomic"
OUT = R / "09_extensions"
FROZEN = C.FROZEN_DIR
sys.path.insert(0, str(ROOT / "scripts"))
shift08 = importlib.import_module("08_shift_tests")


def f(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{float(v):.{nd}f}"


def b1(out):
    feats = pd.read_csv(R / "01_rii" / "rii_features.csv", index_col=0, dtype=str)
    p20 = pd.read_csv(R / "01_rii" / "panel_k20_all_animals.csv", dtype={"protein_id": str})
    pst = pd.read_csv(R / "01_rii" / "stability_k20.csv", dtype={"feature_ID": str})
    pst["ensembl_gene"] = feats["ensembl_gene"].reindex(pst["feature_ID"]).to_numpy()
    pst["gene_symbol"] = feats["gene_symbol"].reindex(pst["feature_ID"]).to_numpy()
    # protein marker tissue for stability proteins: from the all-animal panel where present, else from the RII means (recomputed)
    X = pd.read_parquet(R / "01_rii" / "rii_inner_log2ppm.parquet"); meta = pd.read_csv(R / "01_rii" / "rii_meta.csv", index_col=0, dtype=str).loc[X.index]
    Z = (X - X.mean()) / X.std(ddof=0)
    mu = Z.groupby(meta["tissue"].to_numpy()).mean()
    prot_marker = {pid: mu[pid].idxmax() for pid in pst["feature_ID"] if pid in mu.columns}
    r20 = pd.read_csv(FROZEN / "12_bodymap" / "panel_gene_check.csv")
    r51 = pd.read_csv(FROZEN / "05_panels" / "TRNSCRPT" / "stability_k20_annotated.csv")
    rna = {}
    for _, r in r20.iterrows():
        rna.setdefault(r["feature_ID"], {"ensembl_gene": r["feature_ID"], "gene_symbol": r["gene_symbol"], "rna_marker_tissue": r["marker_tissue"], "in_rna_k20": True, "rna_stability_frequency": np.nan})
    for _, r in r51.iterrows():
        e = rna.setdefault(r["feature_ID"], {"ensembl_gene": r["feature_ID"], "gene_symbol": r["gene_symbol"], "rna_marker_tissue": r["marker_tissue"], "in_rna_k20": False, "rna_stability_frequency": np.nan})
        e["rna_stability_frequency"] = float(r["selection_frequency"]); e["rna_marker_tissue"] = r["marker_tissue"]
    prot = {}
    for _, r in p20.iterrows():
        if pd.isna(r["ensembl_gene"]):
            continue
        prot.setdefault(r["ensembl_gene"], {"protein_ids": set(), "in_protein_k20": False, "protein_stability_frequency": 0.0, "protein_marker_tissue": None})
        prot[r["ensembl_gene"]]["protein_ids"].add(r["protein_id"]); prot[r["ensembl_gene"]]["in_protein_k20"] = True; prot[r["ensembl_gene"]]["protein_marker_tissue"] = r["marker_tissue"]
    for _, r in pst.iterrows():
        if pd.isna(r["ensembl_gene"]):
            continue
        e = prot.setdefault(r["ensembl_gene"], {"protein_ids": set(), "in_protein_k20": False, "protein_stability_frequency": 0.0, "protein_marker_tissue": None})
        e["protein_ids"].add(r["feature_ID"]); e["protein_stability_frequency"] = max(e["protein_stability_frequency"], float(r["selection_frequency"]))
        if e["protein_marker_tissue"] is None:
            e["protein_marker_tissue"] = prot_marker.get(r["feature_ID"])
    PROT7 = {"CORTEX", "HEART", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC"}
    rows = []
    for g in sorted(set(rna) | set(prot)):
        a, b = rna.get(g, {}), prot.get(g, {})
        rows.append({"ensembl_gene": g, "gene_symbol": a.get("gene_symbol") or feats.loc[feats["ensembl_gene"] == g, "gene_symbol"].dropna().iloc[0] if (a.get("gene_symbol") or (feats["ensembl_gene"] == g).any()) else g,
                     "in_rna_k20": bool(a.get("in_rna_k20", False)), "rna_stability_frequency": a.get("rna_stability_frequency", np.nan), "rna_marker_tissue": a.get("rna_marker_tissue"),
                     "in_protein_k20": bool(b.get("in_protein_k20", False)), "protein_stability_frequency": b.get("protein_stability_frequency", np.nan), "protein_in_stable_core": bool(b.get("protein_stability_frequency", 0) >= 0.8),
                     "protein_marker_tissue": b.get("protein_marker_tissue"), "protein_ids": ";".join(sorted(b.get("protein_ids", ()))),
                     "in_both_layers": bool(a) and bool(b), "rna_marker_in_prot7": a.get("rna_marker_tissue") in PROT7 if a else None})
    ov = pd.DataFrame(rows)
    ov["markers_agree"] = np.where(ov["in_both_layers"] & ov["rna_marker_in_prot7"].fillna(False).astype(bool), ov["rna_marker_tissue"] == ov["protein_marker_tissue"], None)
    ov = ov.sort_values(["in_both_layers", "in_rna_k20", "in_protein_k20"], ascending=False)
    ov.to_csv(out / "panel_overlap.csv", index=False)
    both = ov[ov["in_both_layers"]]
    testable = both[both["rna_marker_in_prot7"].fillna(False).astype(bool)]
    summ = pd.DataFrame([{"n_rna_candidates": int(len(rna)), "n_rna_k20": int(ov["in_rna_k20"].sum()), "n_protein_candidates": int(len(prot)), "n_protein_k20": int(ov["in_protein_k20"].sum()), "n_protein_stable_core_ge_0.8": int(ov["protein_in_stable_core"].sum()),
                          "k20_intersection": int((ov["in_rna_k20"] & ov["in_protein_k20"]).sum()), "rna_k20_x_protein_core": int((ov["in_rna_k20"] & ov["protein_in_stable_core"]).sum()),
                          "genes_in_both_candidate_lists": int(len(both)), "genes_in_both_with_rna_marker_in_prot7": int(len(testable)), "n_markers_agree": int(testable["markers_agree"].fillna(False).astype(bool).sum()),
                          "frac_markers_agree": float(testable["markers_agree"].fillna(False).astype(bool).mean()) if len(testable) else np.nan, "prereg_B1_threshold": 0.70,
                          "prereg_B1_pass": bool(len(testable) and testable["markers_agree"].fillna(False).astype(bool).mean() >= 0.70),
                          "genes_in_both": ";".join(both["gene_symbol"].astype(str))}])
    summ.to_csv(out / "overlap_summary.csv", index=False)
    print(summ.T.to_string())
    return ov, summ


def b2(out, seed, k=20, alpha=0.10, recal_ns=(3, 5), recal_repeats=10):
    X = pd.read_parquet(R / "01_rii" / "rii_inner_log2ppm.parquet")
    meta = pd.read_csv(R / "01_rii" / "rii_meta.csv", index_col=0, dtype=str).loc[X.index]
    om = io.OmicsMatrix(X, meta, pd.DataFrame(index=X.columns), "PROT-RII", name="RII")
    Xn = X.to_numpy(dtype=float); y = meta["tissue"].astype(str).to_numpy(); g = meta["pid"].astype(str).to_numpy()
    rows, recals = [], []
    # in-distribution: 5 animal-grouped folds; fit on 70 % of the training animals, calibrate on 30 %, test on the held-out animals
    for fold, (tr, te) in enumerate(grouped_kfold(meta, "tissue", 5, seed)):
        fit_idx, cal_idx = fit_calibration_split(meta, tr, 0.3, seed)
        for arm, kk in (("full", None), (f"panel_k{k}", k)):
            est = models.fit_tuned("logreg_l2", Xn[fit_idx], y[fit_idx], g[fit_idx], k=kk, prefilter=5000, quick=False, seed=seed)
            cls = list(est.classes_)
            scores, _ = cp.calibration_scores(est.predict_proba(Xn[cal_idx]), y[cal_idx], cls, "lac")
            q = cp.conformal_quantile(scores, alpha)
            p = est.predict_proba(Xn[te]); sets = cp.predict_sets(p, q, "lac")
            yi = np.array([cls.index(v) for v in y[te]])
            rows.append({"split": "in_distribution", "fold": fold, "arm": arm, "n_train_animals": len(np.unique(g[fit_idx])), "n_cal_animals": len(np.unique(g[cal_idx])), "n_test_animals": len(np.unique(g[te])), "n_test": len(te),
                         "accuracy_all": float(np.mean(est.predict(Xn[te]) == y[te])), "coverage_target_all": float(sets[np.arange(len(yi)), yi].mean()), "avg_set_size_target": float(sets.sum(axis=1).mean()),
                         "lac_frac_empty_target": float((sets.sum(axis=1) == 0).mean()), "n_classes": len(cls)})
    # the phase-08 shifts through the phase-08 code
    splits = list(leave_one_sex_out(meta)) + list(train_controls_test_trained(meta))
    for name, (tr, te) in splits:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            r, pc, rc, vr = shift08.eval_split(name, tr, te, Xn, y, g, meta, k, alpha, "logreg_l2", 5000, False, seed, recal_ns=recal_ns, recal_repeats=recal_repeats)
        for rr in r:
            rr["fold"] = np.nan; rr["n_classes"] = 7
        rows.extend(r); recals.extend(rc)
        print(f"  {name}: " + "; ".join(f"{rr['arm']} acc {rr['accuracy_all']:.3f} cov {rr['coverage_target_all']:.3f} size {rr['avg_set_size_target']:.2f}" for rr in r))
    tab = pd.DataFrame(rows)
    tab.to_csv(out / "protein_ladder_raw.csv", index=False)
    recal = pd.DataFrame(recals)
    recal["n_classes"] = 7
    recal.to_csv(out / "protein_ladder_recalibration.csv", index=False)
    # summary per split × arm (mean over folds for in-distribution)
    agg_cols = ["accuracy_all", "coverage_target_all", "avg_set_size_target", "lac_frac_empty_target"]
    summ = tab.groupby(["split", "arm"], sort=False).agg(accuracy=("accuracy_all", "mean"), accuracy_sd=("accuracy_all", "std"), coverage=("coverage_target_all", "mean"), coverage_sd=("coverage_target_all", "std"),
                                                          set_size=("avg_set_size_target", "mean"), frac_empty=("lac_frac_empty_target", "mean"), n_test_samples=("n_test", "sum"), n_test_animals=("n_test_animals", "sum"),
                                                          n_folds=("fold", "nunique"), n_train_animals=("n_train_animals", "first"), n_cal_animals=("n_cal_animals", "first"), n_classes=("n_classes", "first")).reset_index()
    # for the shift rows n_test_animals is the split's count (one row), for in-distribution the sum over folds = 60
    summ["layer"] = "protein (RII)"
    summ.to_csv(out / "protein_ladder.csv", index=False)
    # RNA rows from the frozen results
    rna_rows = []
    st = pd.read_csv(FROZEN / "08_shift" / "TRNSCRPT" / "shift_table.csv")
    for r in st[st["split"].isin(["train_male_test_female", "train_female_test_male", "train_control_test_trained"])].itertuples():
        rna_rows.append({"split": r.split, "arm": r.arm, "layer": "RNA (counts)", "accuracy": float(r.accuracy_all), "coverage": float(r.coverage_target_all), "set_size": float(r.avg_set_size_target), "frac_empty": float(getattr(r, "lac_frac_empty_target", np.nan)) if hasattr(r, "lac_frac_empty_target") else np.nan,
                         "n_test_samples": int(r.n_test), "n_test_animals": int(r.n_test_animals), "n_classes": 19, "source": "results_frozen/08_shift/TRNSCRPT/shift_table.csv"})
    pc5 = pd.read_csv(FROZEN / "05_panels" / "TRNSCRPT" / "panel_curve.csv"); cov6 = pd.read_csv(FROZEN / "06_conformal" / "TRNSCRPT" / "coverage.csv")
    c6 = cov6[(cov6["calibration"] == "pooled") & (cov6["conformal"] == "marginal") & (cov6["method"] == "lac") & (cov6["alpha"] == 0.1)]
    hl = json.loads((ROOT / "site" / "data" / "headline.json").read_text())
    k20_id = next(x for x in hl["ladder"] if x["rung_id"] == "in_distribution" and x["model"] == "k20" and x["variant"] == "marginal" and (x.get("calibration") or "pooled") == "pooled")
    rna_rows.append({"split": "in_distribution", "arm": "panel_k20", "layer": "RNA (counts)", "accuracy": float(pc5.loc[pc5["k"] == 20, "balanced_accuracy"].mean()), "coverage": float(k20_id["coverage"]), "set_size": float(k20_id["set_size"]), "frac_empty": float(k20_id["empty"]),
                     "n_test_samples": int(k20_id["n_samples"]), "n_test_animals": int(k20_id["n_individuals"]), "n_classes": 19, "source": "results_frozen/05_panels/TRNSCRPT/panel_curve.csv (accuracy, mean over folds); site/data/headline.json in_distribution k20 marginal pooled (coverage, recomputed from results/31_site_regen/06_conformal)"})
    rna_rows.append({"split": "in_distribution", "arm": "full", "layer": "RNA (counts)", "accuracy": np.nan, "coverage": float(c6["coverage"].mean()), "set_size": float(c6["avg_set_size"].mean()), "frac_empty": float(c6["frac_empty"].mean()),
                     "n_test_samples": int(c6["n"].sum()), "n_test_animals": int(c6["n_test_animals"].sum()), "n_classes": 19, "source": "results_frozen/06_conformal/TRNSCRPT/coverage.csv (pooled, marginal, lac, α = 0.1; mean over folds)"})
    rna = pd.DataFrame(rna_rows)
    side = pd.concat([summ.assign(source="results_multiomic/09_extensions/protein_ladder.csv"), rna], ignore_index=True)
    side["arm"] = side["arm"].replace({"panel_k20": "k20"})
    side.to_csv(out / "ladder_side_by_side.csv", index=False)
    print(side[["split", "arm", "layer", "accuracy", "coverage", "set_size", "n_test_samples", "n_test_animals"]].round(3).to_string(index=False))
    return summ, recal, side


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    ov, s1 = b1(OUT)
    summ, recal, side = b2(OUT, C.SEED)
    s = s1.iloc[0]
    pk = lambda split, arm: summ[(summ["split"] == split) & (summ["arm"] == arm)].iloc[0]
    idk, ctk, mfk, fmk = pk("in_distribution", "panel_k20"), pk("train_control_test_trained", "panel_k20"), pk("train_male_test_female", "panel_k20"), pk("train_female_test_male", "panel_k20")
    rules = {"B2_i_held_out_animals": bool(idk.accuracy >= 0.95 and abs(idk.coverage - 0.90) <= 0.02), "B2_ii_controls_to_trained": bool(ctk.accuracy >= 0.90 and ctk.coverage >= 0.85),
             "B2_iii_held_out_sex_either": bool((mfk.accuracy >= 0.85 and mfk.coverage >= 0.80) or (fmk.accuracy >= 0.85 and fmk.coverage >= 0.80))}
    rec = lambda split, n: recal[(recal["split"] == split) & (recal["arm"] == "panel_k20") & (recal["n_recal_animals"] == n)].iloc[0]
    rna = side[side["layer"].str.startswith("RNA")].set_index(["split", "arm"])
    lines = [f"# Run 2 B — two extensions", "", f"Built by `scripts/multiomic/09_extensions.py` on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}; pre-registered in `docs/PREREGISTRATION_MULTIOMIC.md` (run 2).", "",
             "## B1 — RNA vs protein panels — `panel_overlap.csv`, `overlap_summary.csv`", "", report.df_to_md(s1.T.reset_index().rename(columns={"index": "quantity", 0: "value"})), "",
             report.df_to_md(ov[ov["in_both_layers"] | ov["in_protein_k20"]].drop(columns=["protein_ids"]), floatfmt=".2f", max_rows=60), "",
             "## B2 — the protein ladder inside MoTrPAC (RII) beside the RNA ladder — `protein_ladder.csv`, `ladder_side_by_side.csv`, `protein_ladder_recalibration.csv`", "",
             report.df_to_md(side[["split", "arm", "layer", "accuracy", "accuracy_sd", "coverage", "coverage_sd", "set_size", "frac_empty", "n_test_samples", "n_test_animals", "n_classes"]], floatfmt=".3f"), "",
             report.df_to_md(recal, floatfmt=".3f"), "", "Pre-registered rules: " + ", ".join(f"{k} {'PASS' if v else 'FAIL'}" for k, v in rules.items()), ""]
    (OUT / "README.md").write_text("\n".join(lines))
    sec = [f"- **B1, panels across layers** (`results_multiomic/09_extensions/overlap_summary.csv`): the RNA k20 panel ({int(s.n_rna_k20)} genes, 19 tissues) and the protein k20 panel ({int(s.n_protein_k20)} proteins, 7 tissues) share **{int(s.k20_intersection)}** gene(s); "
           f"RNA k20 ∩ protein stable core (≥ 0.80 in 50 bootstraps, {int(s['n_protein_stable_core_ge_0.8'])} proteins): {int(s.rna_k20_x_protein_core)}. Over the wider candidate lists ({int(s.n_rna_candidates)} RNA, {int(s.n_protein_candidates)} protein genes) {int(s.genes_in_both_candidate_lists)} genes appear in both"
           + (f" ({s.genes_in_both})" if s.genes_in_both else "") + f"; {int(s.genes_in_both_with_rna_marker_in_prot7)} have their RNA marker among the 7 tissues and {int(s.n_markers_agree)} of those agree on the marker tissue → pre-registered rule B1 {'PASS' if s.prereg_B1_pass else 'FAIL (or untestable)'}. "
           "Read: the layers pick different genes for the same tissues; the Phase 1 result (RNA markers keep their marker tissue as proteins, 17 of 20) is about the RNA genes' behaviour at the protein level, not about the protein selector choosing them.",
           f"- **B2, the protein ladder inside MoTrPAC** (`ladder_side_by_side.csv`; RII scale, phase-08 code path, k20, 7 classes): held-out animals accuracy {f(idk.accuracy)} ± {f(idk.accuracy_sd)}, coverage {f(idk.coverage)} ± {f(idk.coverage_sd)} at set size {f(idk.set_size, 2)} (RNA: {f(rna.loc[('in_distribution', 'k20'), 'accuracy'])} / {f(rna.loc[('in_distribution', 'k20'), 'coverage'])}, 19 classes); "
           f"fit on the sedentary controls, tested on the trained animals: accuracy {f(ctk.accuracy)}, coverage {f(ctk.coverage)} at {f(ctk.set_size, 2)} ({int(ctk.n_test_samples)} vials, {int(ctk.n_test_animals)} animals; RNA {f(rna.loc[('train_control_test_trained', 'k20'), 'accuracy'])} / {f(rna.loc[('train_control_test_trained', 'k20'), 'coverage'])}); "
           f"held-out sex M → F {f(mfk.accuracy)} / {f(mfk.coverage)} at {f(mfk.set_size, 2)}, F → M {f(fmk.accuracy)} / {f(fmk.coverage)} at {f(fmk.set_size, 2)} (RNA {f(rna.loc[('train_male_test_female', 'k20'), 'accuracy'])} / {f(rna.loc[('train_male_test_female', 'k20'), 'coverage'])}, {f(rna.loc[('train_female_test_male', 'k20'), 'accuracy'])} / {f(rna.loc[('train_female_test_male', 'k20'), 'coverage'])}). "
           f"Recalibration on 3 trained animals: coverage {f(rec('train_control_test_trained', 3).coverage_target_recalibrated)} at set size {f(rec('train_control_test_trained', 3).set_size_recalibrated, 2)} of 7. Pre-registered rules: "
           + ", ".join(f"{k.replace('_', ' ')} {'PASS' if v else 'FAIL'}" for k, v in rules.items()) + ". Within-study: plex is nested in tissue; context for the external rungs, not a finding."]
    (OUT / "REPORT_SECTION.md").write_text("\n".join(sec) + "\n")
    json.dump({"rules": rules, "minutes": (time.time() - t0) / 60}, open(OUT / "rules.json", "w"), indent=1)
    (OUT / "STATUS.json").write_text(json.dumps({"status": "DONE", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}))
    print(rules); print(f"wrote {OUT} in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
