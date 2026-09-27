#!/usr/bin/env python
"""Multiomic Phase 3b — the MoTrPAC RII protein panels scored on two one-sample-per-tissue atlases:
  wang2019   Wang et al. 2019 MSB 15:e8503, 29 human tissues, label-free gene-level intensities (Table EV1 sheet C), matched
             to RII genes through 1:1 rat–human orthologs (Ensembl ids);
  geiger2013 Geiger et al. 2013 MCP 12:1709, 28 mouse tissues, SILAC H/L ratios to one common SILAC-mouse standard
             (Table S1), matched by gene symbol (rat ↔ mouse; caveat: symbol match, not an orthology table).
Each atlas sample is one tissue from one (or one pooled) individual, so accuracy is over samples and every split is
grouped on the sample; conformal coverage with MoTrPAC calibration is reported, recalibration on 3 and 5 target samples
(each a different individual/tissue) is reported for completeness. Same code path as Phase 3 (tfp.transfer).
Outputs: results_multiomic/03_prot_transfer/<atlas>/.
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

from tfp import config as C, report
from tfp.splits import assert_no_group_leak
from tfp.transfer import PanelModels, calibrate_models, conformal_transfer, one_to_one_orthologs, per_organ_coverage, score_block

ROOT = C.ROOT
EXT = ROOT / "data" / "external_multiomic"
RII = ROOT / "results_multiomic" / "01_rii"

MAPS = {
    "wang2019": {"Brain": "CORTEX", "Heart": "HEART", "Kidney": "KIDNEY", "Liver": "LIVER", "Lung": "LUNG", "Fat": "WAT-SC"},
    "geiger2013": {"Brain cortex": "CORTEX", "Heart": "HEART", "Kidney cortex": "KIDNEY", "Kidney medulla": "KIDNEY", "Liver": "LIVER", "Lung": "LUNG",
                   "Muscle": "SKM-GN", "Diaphragm": "SKM-GN", "White fat": "WAT-SC"},
}
CAVEAT = {"wang2019": "whole brain → CORTEX; subcutaneous 'Fat' → WAT-SC; no skeletal muscle in the atlas (SKM-GN has no target); one adult donor per tissue, label-free intensities.",
          "geiger2013": "pooled mice per tissue; diaphragm and limb muscle both → SKM-GN, kidney cortex and medulla both → KIDNEY; matched by gene symbol across species."}


def f(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{float(v):.{nd}f}"


def load_source():
    Xs = pd.read_parquet(RII / "rii_inner_log2ppm.parquet")
    ms = pd.read_csv(RII / "rii_meta.csv", index_col=0, dtype=str).loc[Xs.index]
    feats = pd.read_csv(RII / "rii_features.csv", index_col=0, dtype=str).reindex(Xs.columns)
    return Xs, ms, feats


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--atlas", required=True, choices=list(MAPS))
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--grid", default="20,50")
    ap.add_argument("--prefilter", type=int, default=5000)
    ap.add_argument("--cal-frac", type=float, default=0.3)
    ap.add_argument("--recal", default="3,5")
    ap.add_argument("--recal-repeats", type=int, default=20)
    ap.add_argument("--min-complete", type=float, default=0.8)
    ap.add_argument("--seed", type=int, default=C.SEED)
    args = ap.parse_args()
    t0 = time.time()
    out = ROOT / "results_multiomic" / "03_prot_transfer" / args.atlas
    out.mkdir(parents=True, exist_ok=True)
    grid = [int(k) for k in args.grid.split(",")]
    rng = np.random.default_rng(args.seed)
    Xs, ms, feats = load_source()
    y = ms["tissue"].astype(str).to_numpy(); g = ms["pid"].astype(str).to_numpy()
    classes_src = sorted(np.unique(y))

    if args.atlas == "wang2019":
        T = pd.read_csv(EXT / "wang2019" / "wang2019_tissue_ibaq.csv", index_col=0)   # genes × tissues (intensity)
        Xt = np.log2(T.where(T > 0)).T                                                      # samples × genes, 0 → NaN
        p2g = feats["ensembl_gene"].dropna()
        G = Xs[p2g.index].T.groupby(p2g.to_numpy()).mean().T
        orth = one_to_one_orthologs()
        orth = orth[orth["RAT_ENSEMBL_ID"].isin(G.columns) & orth["HUMAN_ORTHOLOG_ENSEMBL_ID"].isin(Xt.columns)].drop_duplicates("RAT_ENSEMBL_ID").drop_duplicates("HUMAN_ORTHOLOG_ENSEMBL_ID")
        src_keys, tgt_keys = orth["RAT_ENSEMBL_ID"].tolist(), orth["HUMAN_ORTHOLOG_ENSEMBL_ID"].tolist()
        sym = feats.dropna(subset=["ensembl_gene"]).drop_duplicates("ensembl_gene").set_index("ensembl_gene")["gene_symbol"].to_dict()
        match_note = "1:1 rat–human orthologs by Ensembl id"
    else:
        T = pd.read_csv(EXT / "geiger2013" / "tableS1_ratio_HL_normalized.csv", index_col=0)
        ann = pd.read_csv(EXT / "geiger2013" / "tableS1_annotation.csv", index_col=0, dtype=str)
        T = T.drop(columns=[c for c in T.columns if c.strip() == ""], errors="ignore")
        gsym = ann["Gene names"].astype(str).str.split(";").str[0].str.strip().str.upper()
        T = T.groupby(gsym.reindex(T.index).to_numpy()).mean()                             # gene-symbol level (mean over protein groups)
        T = T[T.index != "NAN"]
        Xt = np.log2(T.where(T > 0)).T
        s2 = feats["gene_symbol"].dropna().str.upper()
        G = Xs[s2.index].T.groupby(s2.to_numpy()).mean().T
        common = sorted(set(G.columns) & set(Xt.columns))
        src_keys, tgt_keys = common, common
        sym = {s: s for s in common}
        match_note = "gene symbol, upper-cased (rat ↔ mouse)"
    completeness = Xt[tgt_keys].notna().mean(axis=0)
    keep = set(completeness[completeness >= args.min_complete].index)
    pairs = [(s, t) for s, t in zip(src_keys, tgt_keys) if t in keep]
    src_keys, tgt_keys = [s for s, _ in pairs], [t for _, t in pairs]
    Xt_k = Xt[tgt_keys]
    Xt_k = Xt_k.fillna(Xt_k.median(axis=0))
    L_src = G[src_keys].to_numpy(dtype=float); L_tgt = Xt_k.to_numpy(dtype=float)
    mt = pd.DataFrame({"tissue": Xt_k.index, "organ": [MAPS[args.atlas].get(t, t) for t in Xt_k.index], "group_id": Xt_k.index, "stage": "adult"}, index=Xt_k.index)
    organ_map = {c: {c} for c in classes_src}
    mapped = mt["organ"].isin(classes_src).to_numpy()
    ood_organs = sorted(mt.loc[~mapped, "organ"].unique())
    ov = pd.DataFrame([{"atlas": args.atlas, "matching": match_note, "rii_genes": G.shape[1], "target_genes": int(Xt.shape[1]), "matched_genes_before_completeness": len(completeness),
                        f"matched_genes_complete_ge_{args.min_complete}": len(tgt_keys), "target_samples": len(mt), "target_samples_mapped": int(mapped.sum()),
                        "mapped_rat_classes": ";".join(sorted(mt.loc[mapped, "organ"].unique())), "source_classes_without_target": ";".join(c for c in classes_src if c not in set(mt["organ"])),
                        "n_ood_samples": int((~mapped).sum())}])
    ov.to_csv(out / "gene_overlap.csv", index=False); print(ov.T.to_string())

    pm = PanelModels(L_src, y, src_keys, grid, args.prefilter)
    Lt_pre, Zt_pre = pm.target_matrices(L_tgt)
    model_names = [f"k{k}" for k in grid] + ["full"]
    rows = []
    for name in model_names:
        pred, ok, p, cls = score_block(pm, name, Lt_pre, Zt_pre, mt, organ_map)
        for i, s in enumerate(mt.index):
            rows.append({"model": name, "atlas_tissue": mt.loc[s, "tissue"], "rat_class": mt.loc[s, "organ"] if mapped[i] else "OOD", "prediction": pred[i],
                         "correct": bool(ok[i]) if mapped[i] else np.nan, "p_max": float(p[i].max())})
    acc = pd.DataFrame(rows)
    acc.to_csv(out / "predictions_by_sample.csv", index=False)
    am = acc[acc["rat_class"] != "OOD"]
    ao = am.groupby("model").agg(accuracy=("correct", "mean"), n_mapped=("correct", "size")).reset_index()
    ao["n_ood"] = int((~mapped).sum()); ao["chance_1_over_7"] = 1 / 7
    ao.to_csv(out / "accuracy_overall.csv", index=False); print(ao.to_string(index=False))
    # conformal with MoTrPAC calibration
    uniq = np.unique(g); rng.shuffle(uniq)
    cal_animals = set(uniq[:max(1, int(round(args.cal_frac * len(uniq))))])
    is_cal = np.array([a in cal_animals for a in g])
    fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
    assert_no_group_leak(ms, fit_idx, cal_idx)
    pm_c = PanelModels(L_src[fit_idx], y[fit_idx], src_keys, grid, args.prefilter)
    classes = pm_c.classes
    Lm_c, Zm_c = pm_c.source_matrices(L_src[cal_idx])
    calib = calibrate_models(pm_c, Lm_c, Zm_c, y[cal_idx], classes, args.alpha, model_names)
    Lt_c, Zt_c = pm_c.target_matrices(L_tgt)
    recal_ns = [int(v) for v in args.recal.split(",") if v.strip()]
    conf_df, ood, recal = conformal_transfer(calib, classes, Lt_c, Zt_c, mt, organ_map, args.alpha, "stage", "adult", recal_ns, args.recal_repeats, rng, model_names, ood_organs)
    conf_df.to_csv(out / "conformal_transfer.csv", index=False); ood.to_csv(out / "ood_sets.csv", index=False); recal.to_csv(out / "recalibration.csv", index=False)
    pt = per_organ_coverage(calib, classes, Lt_c, Zt_c, mt, organ_map, model_names)
    pt.to_csv(out / "coverage_by_tissue.csv", index=False)
    print(conf_df[conf_df["conformal"] == "marginal"].round(3).to_string(index=False))
    k20 = f"k{grid[0]}"
    cmk = conf_df[(conf_df["model"] == k20) & (conf_df["conformal"] == "marginal")].iloc[0]
    o = ov.iloc[0]
    tab = acc[acc["model"] == k20][["atlas_tissue", "rat_class", "prediction", "correct", "p_max"]]
    lines = [f"# Phase 3b — protein transfer to {args.atlas} (one sample per tissue)", "", f"Built by `scripts/multiomic/03b_prot_transfer_atlases.py` on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}; numbers from the CSVs here.", "",
             f"- Matching: {o['matching']}; {int(o[f'matched_genes_complete_ge_{args.min_complete}'])} genes after the ≥ {args.min_complete:.0%} completeness filter (`gene_overlap.csv`).",
             f"- Target: {int(o['target_samples'])} tissue samples; {int(o['target_samples_mapped'])} map to {o['mapped_rat_classes']}; without target: {o['source_classes_without_target'] or 'none'}; {int(o['n_ood_samples'])} OOD samples.",
             f"- Caveats: {CAVEAT[args.atlas]}", "",
             "## Accuracy over mapped samples — `accuracy_overall.csv`, `predictions_by_sample.csv`", "", report.df_to_md(ao, floatfmt=".3f"), "",
             f"Per sample at {k20}:", "", report.df_to_md(tab, floatfmt=".2f", max_rows=40), "",
             f"## Conformal sets with MoTrPAC calibration (α = {args.alpha}) — `conformal_transfer.csv`, `recalibration.csv`, `ood_sets.csv`", "",
             report.df_to_md(conf_df.drop(columns=["stage"]), floatfmt=".3f"), "", report.df_to_md(recal, floatfmt=".3f"), "",
             "## What this does not show", "", "- One sample per tissue: no within-tissue spread, no donor-level statement; recalibration 'individuals' are single tissue samples.", "",
             f"_Run time {(time.time() - t0) / 60:.1f} min._"]
    (out / "README.md").write_text("\n".join(lines))
    summary = {"atlas": args.atlas, "n_mapped": int(mapped.sum()), "n_ood": int((~mapped).sum()), "matched_genes": len(tgt_keys),
               **{f"acc_{r.model}": float(r.accuracy) for r in ao.itertuples()}, "coverage_k20_marginal_source_cal": float(cmk["coverage_mapped"]), "frac_empty_k20_source_cal": float(cmk["frac_empty_mapped"]),
               "ood_frac_empty_k20": float(cmk["ood_frac_empty"]) if not np.isnan(cmk["ood_frac_empty"]) else None}
    (out / "summary.json").write_text(json.dumps(summary, indent=1))
    print(f"wrote {out} in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
