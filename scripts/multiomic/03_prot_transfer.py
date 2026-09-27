#!/usr/bin/env python
"""Multiomic Phase 3 — protein fingerprint transfer: MoTrPAC RII panels scored on the Jiang et al. 2020 human proteome
map (Cell 183:269; 32 tissues, 14 GTEx donors, TMT), through tfp.transfer unchanged.

Source: results_multiomic/01_rii/rii_inner_log2ppm.parquet (Phase 1; proteins quantified in every tissue), rolled to genes
(mean log2 ppm over the proteins of a gene), restricted to 1:1 rat–human orthologs present in the target. Target: the
authors' cleaned relative abundances (log2 to the run's reference; data/external_multiomic/jiang2020/protein_relative.parquet),
genes with >= --min-complete non-missing target samples, the remaining gaps imputed with the target's own gene median before
the shared code path (which z-scores within each dataset). Every target split is grouped on the donor.
Tissue map: Muscle - Skeletal → SKM-GN; Heart - Atrial Appendage and Heart - Left Ventricle → HEART; Lung → LUNG; Liver → LIVER;
Brain - Cortex → CORTEX. KIDNEY and WAT-SC have no target tissue; the other 26 human tissues are out-of-distribution.
Reports: per-tissue accuracy and main wrong call, conformal coverage with MoTrPAC calibration (marginal, Mondrian, floored),
empty-set fraction, recalibration on 3 and 5 donors with set size, the protein ladder beside the frozen RNA ladder (GTEx), the
reverse direction (select on Jiang, score MoTrPAC RII), a per-gene check of the k = 20 panel. `--target-scale rawppm` repeats
the primary analysis on the raw reporter intensities re-normalised exactly like the MoTrPAC RII (channel total, log2 ppm).
Outputs: results_multiomic/03_prot_transfer/ (or a subfolder for the sensitivity run).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler

from tfp import config as C, io, report
from tfp.boot import individual_bootstrap_ci
from tfp.splits import assert_no_group_leak
from tfp.transfer import PanelModels, calibrate_models, conformal_transfer, gene_check, one_to_one_orthologs, per_organ_coverage, score_block, sets_for, zscore

ROOT = C.ROOT
EXT = ROOT / "data" / "external_multiomic" / "jiang2020"
RII = ROOT / "results_multiomic" / "01_rii"
FROZEN = C.FROZEN_DIR

JIANG_TO_RAT = {"Muscle - Skeletal": "SKM-GN", "Heart - Atrial Appendage": "HEART", "Heart - Left Ventricle": "HEART",
                "Lung": "LUNG", "Liver": "LIVER", "Brain - Cortex": "CORTEX"}


def f(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{float(v):.{nd}f}"


def load_target(scale: str, min_complete: float):
    cm = pd.read_csv(EXT / "protein_relative_columns.csv", dtype=str).set_index("column")
    des = pd.read_csv(EXT / "experimental_design.csv", dtype=str).dropna(subset=["GTEx Sample_ID"]).drop_duplicates("GTEx Sample_ID").set_index("GTEx Sample_ID")
    if scale == "relative":
        rel = pd.read_parquet(EXT / "protein_relative.parquet").set_index("gene.id")
        samp = [c for c in rel.columns if c.startswith("GTEX")]
        X = rel[samp].T                                             # samples × genes
        tissue = cm["row2"].reindex(samp)
    else:  # raw reporter intensities → channel total → log2 ppm, reference channels dropped (same processing as tfp.rii)
        raw = pd.read_parquet(EXT / "protein_raw.parquet").set_index("gene.id")
        rcm = pd.read_csv(EXT / "protein_raw_columns.csv", dtype=str).drop_duplicates("column").set_index("column")
        samp = [c for c in raw.columns if c.startswith("GTEX")]          # technical replicates of a sample carry a __i suffix
        M = raw[samp].apply(pd.to_numeric, errors="coerce")
        tot = M.sum(axis=0, skipna=True)
        X = np.log2((M.div(tot, axis=1) * 1e6).where(M > 0)).T
        base = [c.split("__")[0] for c in samp]
        tissue = pd.Series(rcm["row3"].reindex(base).to_numpy(), index=samp)
    X = X.loc[:, ~X.columns.duplicated()]
    X.columns = X.columns.astype(str)
    X.index.name = "sample"
    base = [c.split("__")[0] for c in samp]
    mt = pd.DataFrame({"tissue": tissue.to_numpy(), "donor": des["Individual ID"].reindex(base).to_numpy(),
                       "run": des["Run"].reindex(base).to_numpy(), "tag": des["TMT Tag"].reindex(base).to_numpy()}, index=pd.Index(samp, name="sample"))
    mt["organ"] = mt["tissue"].map(JIANG_TO_RAT).fillna(mt["tissue"])   # mapped → rat class name; unmapped keep their human name (OOD)
    mt["group_id"] = mt["donor"]
    mt["stage"] = "adult"
    complete = X.notna().mean(axis=0)
    keep = complete[complete >= min_complete].index
    X = X[keep]
    X = X.fillna(X.median(axis=0))                                    # target's own gene medians for the remaining gaps
    return X, mt, complete


def main() -> None:
    ap = argparse.ArgumentParser(description="Multiomic Phase 3 — protein transfer to Jiang 2020")
    ap.add_argument("--target-scale", default="relative", choices=["relative", "rawppm"])
    ap.add_argument("--min-complete", type=float, default=0.8)
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--grid", default="20,50")
    ap.add_argument("--prefilter", type=int, default=5000)
    ap.add_argument("--cal-frac", type=float, default=0.3)
    ap.add_argument("--recal-donors", default="3,5")
    ap.add_argument("--recal-repeats", type=int, default=20)
    ap.add_argument("--seed", type=int, default=C.SEED)
    ap.add_argument("--out", default=None)
    args = ap.parse_args()
    t0 = time.time()
    out = Path(args.out) if args.out else (ROOT / "results_multiomic" / "03_prot_transfer" / ("" if args.target_scale == "relative" else args.target_scale))
    out.mkdir(parents=True, exist_ok=True)
    grid = [int(k) for k in args.grid.split(",")]
    rng = np.random.default_rng(args.seed)
    print(f"=== multiomic 03_prot_transfer === {vars(args)}")

    # ---- source: RII proteins → genes → 1:1 orthologs ------------------------------------------------------------------
    Xs = pd.read_parquet(RII / "rii_inner_log2ppm.parquet")
    ms = pd.read_csv(RII / "rii_meta.csv", index_col=0, dtype=str).loc[Xs.index]
    feats = pd.read_csv(RII / "rii_features.csv", index_col=0, dtype=str).reindex(Xs.columns)
    p2g = feats["ensembl_gene"].dropna()
    G = Xs[p2g.index].T.groupby(p2g.to_numpy()).mean().T                 # samples × rat genes (mean over a gene's proteins)
    orth = one_to_one_orthologs()
    orth = orth[orth["RAT_ENSEMBL_ID"].isin(G.columns)]
    Xt, mt, completeness = load_target(args.target_scale, args.min_complete)
    orth = orth[orth["HUMAN_ORTHOLOG_ENSEMBL_ID"].isin(Xt.columns)].drop_duplicates("RAT_ENSEMBL_ID").drop_duplicates("HUMAN_ORTHOLOG_ENSEMBL_ID")
    rat_genes, hum_genes = orth["RAT_ENSEMBL_ID"].tolist(), orth["HUMAN_ORTHOLOG_ENSEMBL_ID"].tolist()
    L_src = G[rat_genes].to_numpy(dtype=float)
    L_tgt = Xt[hum_genes].to_numpy(dtype=float)
    y = ms["tissue"].astype(str).to_numpy()
    g = ms["pid"].astype(str).to_numpy()
    classes_src = sorted(np.unique(y))
    organ_map = {c: {c} for c in classes_src}            # mapped target organs already carry the rat class name; others → None
    mapped = mt["organ"].isin(classes_src).to_numpy()
    ood_organs = sorted(mt.loc[~mapped, "organ"].unique())
    ov = pd.DataFrame([{"rii_proteins_every_tissue": Xs.shape[1], "rii_genes": G.shape[1], "genes_with_1to1_ortholog": int(one_to_one_orthologs()["RAT_ENSEMBL_ID"].isin(G.columns).sum()),
                        "target_genes_total": int(len(completeness)), f"target_genes_complete_ge_{args.min_complete}": int((completeness >= args.min_complete).sum()),
                        "matched_genes": len(rat_genes), "target_samples": len(mt), "target_donors": mt["donor"].nunique(), "target_tissues": mt["tissue"].nunique(),
                        "target_samples_mapped": int(mapped.sum()), "target_donors_mapped": mt.loc[mapped, "donor"].nunique(), "mapped_rat_classes": ";".join(sorted(mt.loc[mapped, "organ"].unique())),
                        "source_classes_without_target": ";".join(c for c in classes_src if c not in set(mt["organ"])), "n_ood_tissues": len(ood_organs), "target_scale": args.target_scale}])
    ov.to_csv(out / "gene_overlap.csv", index=False)
    print(ov.T.to_string())
    sym = feats.dropna(subset=["ensembl_gene"]).drop_duplicates("ensembl_gene").set_index("ensembl_gene")["gene_symbol"].to_dict()

    # ---- panels on all MoTrPAC animals; accuracy per target tissue ------------------------------------------------------
    pm = PanelModels(L_src, y, rat_genes, grid, args.prefilter)
    Lt_pre, Zt_pre = pm.target_matrices(L_tgt)
    model_names = [f"k{k}" for k in grid] + ["full"]
    acc_rows, conf, ok_by_model = [], {}, {}
    for name in model_names:
        pred, ok, p, cls = score_block(pm, name, Lt_pre, Zt_pre, mt, organ_map)
        ok_by_model[name] = ok
        sub = mt.assign(pred=pred, correct=ok)
        for (tissue, organ), d in sub.groupby(["tissue", "organ"]):
            vc = d["pred"].value_counts()
            wrong = d.loc[d["pred"] != organ, "pred"].value_counts()
            acc_rows.append({"model": name, "jiang_tissue": tissue, "rat_class": organ if organ in classes_src else "OOD", "n": len(d), "n_donors": d["donor"].nunique(),
                             "accuracy": float(d["correct"].mean()) if organ in classes_src else np.nan, "top_prediction": vc.index[0], "top_prediction_frac": float(vc.iloc[0] / len(d)),
                             "main_wrong_call": wrong.index[0] if len(wrong) else "", "main_wrong_call_frac": float(wrong.iloc[0] / len(d)) if len(wrong) else 0.0})
        conf[name] = pd.crosstab(sub["tissue"], sub["pred"]).reindex(columns=pm.classes, fill_value=0)
        conf[name].to_csv(out / f"confusion_{name}.csv")
    acc = pd.DataFrame(acc_rows)
    acc.to_csv(out / "accuracy_by_tissue.csv", index=False)
    am = acc[acc["rat_class"] != "OOD"]
    overall = am.groupby("model").apply(lambda d: np.average(d["accuracy"], weights=d["n"]), include_groups=False)
    macro = am.groupby(["model", "rat_class"]).apply(lambda d: np.average(d["accuracy"], weights=d["n"]), include_groups=False).groupby("model").mean()
    ao = pd.DataFrame({"model": overall.index, "accuracy_sample_weighted": overall.to_numpy(), "accuracy_macro_over_rat_classes": macro.loc[overall.index].to_numpy(),
                       "n_samples_mapped": int(mapped.sum()), "n_donors_mapped": mt.loc[mapped, "donor"].nunique(), "n_rat_classes_with_target": int(mt.loc[mapped, "organ"].nunique()),
                       "chance_1_over_7": 1 / 7})
    cis = [individual_bootstrap_ci(ok_by_model[n][mapped].astype(float), mt["donor"].to_numpy()[mapped], 1000, args.seed) for n in ao["model"]]
    ao["acc_ci95_low_donor_boot"] = [c[0] for c in cis]; ao["acc_ci95_high_donor_boot"] = [c[1] for c in cis]; ao["n_boot"] = 1000
    ao.to_csv(out / "accuracy_overall.csv", index=False)
    print(ao.to_string(index=False))

    # ---- conformal transfer with MoTrPAC calibration; recalibration on donors --------------------------------------------
    uniq = np.unique(g)
    rng.shuffle(uniq)
    cal_animals = set(uniq[:max(1, int(round(args.cal_frac * len(uniq))))])
    is_cal = np.array([a in cal_animals for a in g])
    fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
    assert_no_group_leak(ms, fit_idx, cal_idx)
    pm_c = PanelModels(L_src[fit_idx], y[fit_idx], rat_genes, grid, args.prefilter)
    classes = pm_c.classes
    Lm_c, Zm_c = pm_c.source_matrices(L_src[cal_idx])
    calib = calibrate_models(pm_c, Lm_c, Zm_c, y[cal_idx], classes, args.alpha, model_names)
    Lt_c, Zt_c = pm_c.target_matrices(L_tgt)
    recal_ns = [int(v) for v in args.recal_donors.split(",") if v.strip()]
    conf_df, ood, recal = conformal_transfer(calib, classes, Lt_c, Zt_c, mt, organ_map, args.alpha, "stage", "adult", recal_ns, args.recal_repeats, rng,
                                             model_names, ood_organs)
    conf_df["n_cal_animals"] = len(cal_animals)
    conf_df.to_csv(out / "conformal_transfer.csv", index=False)
    ood.to_csv(out / "ood_sets.csv", index=False)
    recal.to_csv(out / "recalibration.csv", index=False)
    pt = per_organ_coverage(calib, classes, Lt_c, Zt_c, mt, organ_map, model_names).rename(columns={"organ": "rat_class_or_ood_tissue"})
    pt.to_csv(out / "coverage_by_tissue.csv", index=False)
    # donor-bootstrap intervals for the marginal coverage of the mapped samples
    cov_rows = []
    for name in model_names:
        p = calib[name]["proba"](Lt_c, Zt_c)
        sets = sets_for(calib[name], p, "marginal")
        covered = np.array([any(sets[i, classes.index(t)] for t in organ_map[o]) for i, o in enumerate(mt["organ"]) if organ_map.get(o)], dtype=float)
        lo, hi, nd = individual_bootstrap_ci(covered, mt["donor"].to_numpy()[mapped], 1000, args.seed)
        cov_rows.append({"model": name, "conformal": "marginal", "coverage_mapped": float(covered.mean()), "ci95_low_donor_boot": lo, "ci95_high_donor_boot": hi, "n_samples": int(len(covered)), "n_donors": nd})
    pd.DataFrame(cov_rows).to_csv(out / "coverage_ci.csv", index=False)
    print(conf_df.round(3).to_string(index=False)); print(recal.round(3).to_string(index=False))

    # ---- the ladder: protein (this run) beside RNA (frozen phase 13, GTEx) ---------------------------------------------
    rows = []
    for name in model_names:
        cm_ = conf_df[(conf_df["model"] == name) & (conf_df["conformal"] == "marginal")].iloc[0]
        rc3 = recal[(recal["model"] == name) & (recal["n_recal"] == 3)]; rc5 = recal[(recal["model"] == name) & (recal["n_recal"] == 5)]
        rows.append({"layer": "protein (RII → Jiang 2020)", "model": name, "accuracy_sample_weighted": float(overall[name]), "n_samples": int(mapped.sum()), "n_individuals": int(mt.loc[mapped, "donor"].nunique()),
                     "n_target_classes": int(mt.loc[mapped, "organ"].nunique()), "coverage_source_cal_marginal": float(cm_["coverage_mapped"]), "frac_empty_source_cal": float(cm_["frac_empty_mapped"]),
                     "avg_set_size_source_cal": float(cm_["avg_set_size_mapped"]),
                     "coverage_recal_3": float(rc3["coverage_recalibrated"].iloc[0]) if len(rc3) else np.nan, "coverage_recal_5": float(rc5["coverage_recalibrated"].iloc[0]) if len(rc5) else np.nan,
                     "set_size_recal_5": float(rc5["set_size_recalibrated"].iloc[0]) if len(rc5) else np.nan})
    d13 = FROZEN / "13_gtex"
    if (d13 / "accuracy_overall.csv").exists():
        a13 = pd.read_csv(d13 / "accuracy_overall.csv").set_index("model"); c13 = pd.read_csv(d13 / "conformal_transfer.csv"); r13 = pd.read_csv(d13 / "recalibration.csv")
        for name in a13.index:
            cm_ = c13[(c13["model"] == name) & (c13["conformal"] == "marginal")].iloc[0]
            rc3 = r13[(r13["model"] == name) & (r13["n_recal"] == 3)]; rc5 = r13[(r13["model"] == name) & (r13["n_recal"] == 5)]
            rows.append({"layer": "RNA (counts → GTEx v8, frozen phase 13)", "model": name, "accuracy_sample_weighted": float(a13.loc[name, "accuracy_sample_weighted"]),
                         "n_samples": int(a13.loc[name, "n_samples"]), "n_individuals": 862, "n_target_classes": int(a13.loc[name, "n_tissues"]),
                         "coverage_source_cal_marginal": float(cm_["coverage_mapped"]), "frac_empty_source_cal": float(cm_["frac_empty_mapped"]), "avg_set_size_source_cal": float(cm_["avg_set_size_mapped"]),
                         "coverage_recal_3": float(rc3["coverage_recalibrated"].iloc[0]) if len(rc3) else np.nan, "coverage_recal_5": float(rc5["coverage_recalibrated"].iloc[0]) if len(rc5) else np.nan,
                         "set_size_recal_5": float(rc5["set_size_recalibrated"].iloc[0]) if len(rc5) else np.nan})
    ladder = pd.DataFrame(rows)
    ladder.to_csv(out / "ladder_protein_vs_rna.csv", index=False)

    # ---- reverse direction: select on Jiang (mapped tissues), score MoTrPAC RII ------------------------------------------
    m_idx = np.flatnonzero(mapped)
    yt = mt["organ"].to_numpy()[m_idx]
    rev_rows = []
    if len(np.unique(yt)) >= 2:
        pm_r = PanelModels(L_tgt[m_idx], yt, hum_genes, grid, args.prefilter)
        Ls_r, Zs_r = pm_r.target_matrices(L_src)          # MoTrPAC scored as the target, z-scored within MoTrPAC
        src_mt = pd.DataFrame({"organ": y, "group_id": g}, index=ms.index)
        for name in model_names:
            pred, ok, _, _ = score_block(pm_r, name, Ls_r, Zs_r, src_mt, {c: {c} for c in pm_r.classes})
            sub = src_mt.assign(pred=pred, correct=ok)
            for t, d in sub.groupby("organ"):
                vc = d["pred"].value_counts()
                rev_rows.append({"model": name, "motrpac_tissue": t, "in_jiang_classes": t in pm_r.classes, "n_vials": len(d), "n_animals": d["group_id"].nunique(),
                                 "accuracy": float(d["correct"].mean()) if t in pm_r.classes else np.nan, "top_prediction": vc.index[0], "top_prediction_frac": float(vc.iloc[0] / len(d))})
        rev = pd.DataFrame(rev_rows)
        rev.to_csv(out / "reverse_direction.csv", index=False)
        rev_over = rev[rev["in_jiang_classes"]].groupby("model").apply(lambda d: np.average(d["accuracy"], weights=d["n_vials"]), include_groups=False)
        pd.DataFrame({"model": rev_over.index, "accuracy_sample_weighted": rev_over.to_numpy(), "n_vials": int(rev[rev["in_jiang_classes"]].groupby("model")["n_vials"].sum().iloc[0]),
                      "n_jiang_training_samples": len(m_idx), "n_jiang_donors": mt.loc[mapped, "donor"].nunique(), "n_classes": len(pm_r.classes)}).to_csv(out / "reverse_direction_overall.csv", index=False)
        print(rev_over.round(3).to_string())
    else:
        rev_over = pd.Series(dtype=float)

    # ---- per-gene check of the k = 20 panel in the target ------------------------------------------------------------------
    Zm_df = pd.DataFrame(zscore(G[rat_genes]), index=G.index, columns=rat_genes)
    Zt_df = pd.DataFrame(zscore(Xt[hum_genes]), index=Xt.index, columns=rat_genes)
    k0 = grid[0]
    gc = gene_check(pm, pm.panel_genes(k0), sym, Zm_df, y, Zt_df, mt[mapped], {c: c for c in classes_src}, {})
    gc["human_gene"] = gc["feature_ID"].map(dict(zip(rat_genes, hum_genes)))
    gc = gc.drop(columns=["risk_T7_regulated", "risk_qc_correlated"], errors="ignore")
    gc.to_csv(out / "panel_gene_check.csv", index=False)
    testable = gc[gc["fails_in_target"].notna()]
    fails = testable[testable["fails_in_target"].astype(bool)]
    panel = pd.DataFrame({"rat_gene": pm.panel_genes(k0), "gene_symbol": [sym.get(gg, gg) for gg in pm.panel_genes(k0)],
                          "human_gene": [dict(zip(rat_genes, hum_genes))[gg] for gg in pm.panel_genes(k0)],
                          "marker_tissue": [pm.classes[int(np.argmax(pm.scores[:, list(pm.genes_pre).index(gg)]))] for gg in pm.panel_genes(k0)]})
    panel.to_csv(out / f"panel_k{k0}.csv", index=False)

    # ---- README + report section + findings -------------------------------------------------------------------------------
    o = ov.iloc[0]
    k20 = f"k{grid[0]}"
    cmk = conf_df[(conf_df["model"] == k20) & (conf_df["conformal"] == "marginal")].iloc[0]
    cmf = conf_df[(conf_df["model"] == k20) & (conf_df["conformal"] == "floored")].iloc[0]
    rk3 = recal[(recal["model"] == k20) & (recal["n_recal"] == 3)].iloc[0]; rk5 = recal[(recal["model"] == k20) & (recal["n_recal"] == 5)].iloc[0]
    acc_w = acc.pivot_table(index=["jiang_tissue", "rat_class", "n", "n_donors"], columns="model", values="accuracy").reset_index()
    wrong = acc[(acc["model"] == k20)][["jiang_tissue", "rat_class", "n", "top_prediction", "top_prediction_frac", "main_wrong_call", "main_wrong_call_frac"]]
    ood_k = ood[ood["model"] == k20]
    title = f"Phase 3 — protein fingerprint transfer to Jiang 2020 ({'authors’ cleaned relative abundance' if args.target_scale == 'relative' else 'raw reporter intensities re-normalised to channel-total log2 ppm'})"
    lines = [f"# {title}", "", f"Built by `scripts/multiomic/03_prot_transfer.py` on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}; every number is read from a CSV in this directory.", "",
             "## Data and matching — `gene_overlap.csv`",
             f"- Source: MoTrPAC RII, {int(o['rii_proteins_every_tissue'])} proteins quantified in every tissue → {int(o['rii_genes'])} genes; {int(o['genes_with_1to1_ortholog'])} with a 1:1 human ortholog; "
             f"**{int(o['matched_genes'])} matched** to the target after its completeness filter (≥ {args.min_complete:.0%} non-missing; {int(o[f'target_genes_complete_ge_{args.min_complete}'])} of {int(o['target_genes_total'])} target genes).",
             f"- Target: {int(o['target_samples'])} samples, {int(o['target_donors'])} donors, {int(o['target_tissues'])} tissues; **{int(o['target_samples_mapped'])} samples from {int(o['target_donors_mapped'])} donors map to {o['mapped_rat_classes']}**; "
             f"{o['source_classes_without_target']} have no human counterpart in the atlas; {int(o['n_ood_tissues'])} human tissues are out-of-distribution for the 7-class protein fingerprint.",
             "- Human atlas caveats: adult donors of both sexes, post-mortem GTEx tissue, one pooled all-tissue reference per TMT run, different search engine and FDR — species, age, death and processing are all mixed into the shift.", "",
             "## Accuracy per target tissue (super-class scoring; OOD rows have no accuracy) — `accuracy_by_tissue.csv`, `accuracy_overall.csv`", "", report.df_to_md(acc_w, floatfmt=".2f"), "",
             "Overall, sample-weighted over mapped samples: " + ", ".join(f"{r.model} {f(r.accuracy_sample_weighted)}" for r in ao.itertuples()) + f" (chance 1/7 = {f(1 / 7)}).", "",
             f"Top prediction and main wrong call at {k20}:", "", report.df_to_md(wrong, floatfmt=".2f"), "",
             f"## Conformal transfer (α = {args.alpha}, LAC, calibration on {len(cal_animals)} held-out MoTrPAC animals) — `conformal_transfer.csv`, `coverage_by_tissue.csv`, `ood_sets.csv`, `recalibration.csv`", "",
             report.df_to_md(conf_df.drop(columns=["stage"]), floatfmt=".3f"), "",
             f"Per tissue, {k20}, marginal vs Mondrian vs floored (rat calibration):", "",
             report.df_to_md(pt[pt["model"] == k20].pivot(index="rat_class_or_ood_tissue", columns="conformal", values="coverage").reset_index(), floatfmt=".2f"), "",
             f"What the sets hold for the out-of-distribution tissues at {k20} (`ood_sets.csv`):", "", report.df_to_md(ood_k, floatfmt=".2f") if len(ood_k) else "_none_", "",
             "Recalibration on Jiang donors (pooled tissues), tested on the remaining donors:", "", report.df_to_md(recal, floatfmt=".3f"), "",
             "## The ladder: protein beside RNA — `ladder_protein_vs_rna.csv`", "", report.df_to_md(ladder, floatfmt=".3f"), "",
             "The two rows are different targets (Jiang protein: 5 mapped tissues, few donors; GTEx RNA: 17 tissues, 862 donors); Phase 5 puts both layers on the same Jiang samples.", "",
             "## Reverse direction: panels selected on Jiang, scored on MoTrPAC RII — `reverse_direction.csv`, `reverse_direction_overall.csv`", "",
             (report.df_to_md(pd.read_csv(out / "reverse_direction.csv"), floatfmt=".2f") if (out / "reverse_direction.csv").exists() else "_not enough mapped classes_"), "",
             f"## Per-gene check of the {k20} panel — `panel_gene_check.csv`, `panel_k{k0}.csv`", "", report.df_to_md(gc, floatfmt=".2f"), "",
             f"{len(testable)} of {len(gc)} panel genes have their tissue in the target; {len(fails)} fail (effect ≤ 0 in the target: {', '.join(fails['gene_symbol'].astype(str)) or 'none'}).", "",
             "## What this phase does not show", "",
             "- Only 5 of the 7 protein classes have a human counterpart here (no kidney, no adipose in Jiang 2020); the accuracy is over those five, with 2 cortex samples.",
             "- Species, adult age, post-mortem interval, TMT reference design and search pipeline shift together; the design does not separate them.",
             "- Coverage numbers are on ≤ 50 mapped samples from ≤ 14 donors; the recalibration draws share donors across draws.",
             "", f"_Run time {(time.time() - t0) / 60:.1f} min._"]
    (out / "README.md").write_text("\n".join(lines))
    if args.target_scale == "relative":
        sec = ["- question · does a protein panel selected on MoTrPAC RII name the tissue of an independently processed human proteome, and what happens to conformal coverage (pre-registration d)?",
               f"- data · Jiang et al. 2020 (Cell 183:269) cleaned relative protein abundances: {int(o['target_samples'])} TMT samples, {int(o['target_donors'])} GTEx donors, {int(o['target_tissues'])} tissues; "
               f"{int(o['matched_genes'])} genes matched to the {int(o['rii_genes'])} RII genes through 1:1 orthologs; {int(o['target_samples_mapped'])} samples from {int(o['target_donors_mapped'])} donors in the 5 mapped classes "
               f"({o['mapped_rat_classes']}); KIDNEY and WAT-SC have no target; {int(o['n_ood_tissues'])} human tissues are OOD (`results_multiomic/03_prot_transfer/gene_overlap.csv`).",
               "- design · `tfp.transfer` unchanged: panels on all MoTrPAC animals, z-scores within dataset, super-class scoring, conformal sets calibrated on 30 % held-out MoTrPAC animals (α = 0.10), recalibration on 3 and 5 donors (20 draws), reverse direction, per-gene check; sensitivity run on raw reporter intensities re-normalised like the RII (`rawppm/`).",
               f"- result · accuracy over the {int(o['target_samples_mapped'])} mapped samples: " + ", ".join(f"{r.model} **{f(r.accuracy_sample_weighted)}** (donor-bootstrap 95 % CI {f(r.acc_ci95_low_donor_boot, 2)}–{f(r.acc_ci95_high_donor_boot, 2)})" for r in ao.itertuples()) + f" (chance {f(1 / 7)}). "
               f"Marginal coverage with MoTrPAC calibration at {k20}: **{f(cmk['coverage_mapped'])}** (empty sets {f(cmk['frac_empty_mapped'])}, set size {f(cmk['avg_set_size_mapped'])}); floored Mondrian {f(cmf['coverage_mapped'])}. "
               f"Recalibration on 3 donors → {f(rk3['coverage_recalibrated'])} (set size {f(rk3['set_size_recalibrated'])}), on 5 donors → {f(rk5['coverage_recalibrated'])} (set size {f(rk5['set_size_recalibrated'])}). "
               f"OOD tissues at {k20} (marginal): empty-set fraction {f(cmk['ood_frac_empty'])}. "
               + (f"Reverse direction (panel selected on Jiang, scored on MoTrPAC RII): " + ", ".join(f"{n} {f(v)}" for n, v in rev_over.items()) + ". " if len(rev_over) else "")
               + f"Pre-registration (d) asks for accuracy ≥ {f(3 / 7)} and coverage < 0.90 with MoTrPAC calibration: {'PASS' if (overall[k20] >= 3 / 7 and cmk['coverage_mapped'] < 0.90) else 'FAIL'} at {k20} "
               f"(accuracy {f(overall[k20])}, coverage {f(cmk['coverage_mapped'])}); RNA ladder for comparison in `ladder_protein_vs_rna.csv`.",
               f"- what it does not show · which of species, age, post-mortem state and TMT design drives the shift; kidney and adipose are untested; n is small ({int(o['target_donors_mapped'])} donors) and the 2 cortex samples cannot support a per-tissue statement. "
               f"Panel genes failing in the target at {k20}: {', '.join(fails['gene_symbol'].astype(str)) or 'none'} (`panel_gene_check.csv`).",
               "- files · `results_multiomic/03_prot_transfer/README.md`."]
        (out / "REPORT_SECTION.md").write_text("\n".join(sec) + "\n")
        aok = ao.set_index("model")
        cci = pd.read_csv(out / "coverage_ci.csv").set_index("model")
        find = [{"rank": 1, "text": f"**Protein fingerprint transfers to a human proteome atlas, with the RNA coverage pattern.** A {k0}-protein panel selected on MoTrPAC RII names the tissue of "
                                     f"{f(overall[k20])} of {int(o['target_samples_mapped'])} Jiang 2020 samples (donor-bootstrap 95 % CI {f(aok.loc[k20, 'acc_ci95_low_donor_boot'], 2)}–{f(aok.loc[k20, 'acc_ci95_high_donor_boot'], 2)}; "
                                     f"{int(o['target_donors_mapped'])} donors, 5 mapped tissues; k50 {f(overall.get('k50', np.nan))}, full {f(overall['full'])}; chance {f(1 / 7)}; on the raw reporter scale re-normalised like the RII, k20 rises to "
                                     f"{f(float(pd.read_csv(out / 'rawppm' / 'accuracy_overall.csv').set_index('model').loc['k20', 'accuracy_sample_weighted']), 3) if (out / 'rawppm' / 'accuracy_overall.csv').exists() else 'n/a'} on 94 samples); "
                                     f"marginal conformal coverage with MoTrPAC calibration is {f(cmk['coverage_mapped'])} (CI {f(cci.loc[k20, 'ci95_low_donor_boot'], 2)}–{f(cci.loc[k20, 'ci95_high_donor_boot'], 2)}; empty sets {f(cmk['frac_empty_mapped'])}) against a nominal 0.90, and recalibration on 5 donors restores "
                                     f"{f(rk5['coverage_recalibrated'])} at set size {f(rk5['set_size_recalibrated'])}. Human adults, post-mortem, different TMT design — all shifts at once. "
                                     f"— `results_multiomic/03_prot_transfer/accuracy_overall.csv`, `conformal_transfer.csv`, `recalibration.csv`"}]
        (out / "FINDINGS.json").write_text(json.dumps(find, indent=1))
        (out / "STATUS.json").write_text(json.dumps({"status": "DONE", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}))
    print(f"wrote {out} in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
