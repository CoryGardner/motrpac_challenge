#!/usr/bin/env python
"""Phase 13 — human transfer: score GTEx with the MoTrPAC rat fingerprint (docs/GTEX_TRANSFER.md),
through the same code path as the rat BodyMap validation (tfp.transfer).

Inputs: data/external/gtex_tpm_subset.csv (samples × human genes, log2(TPM+1); scripts/11_gtex_prepare.py)
and gtex_meta.csv (SAMPID, donor, SMTSD, rat_tissue). Orthology: data/raw/rat_to_human_gene.csv,
1:1 pairs only (one human gene per rat gene and one rat gene per human gene), matched to GTEx by
human Ensembl id, symbol as fallback. Every GTEx split is grouped on the donor.

Mapping: Muscle - Skeletal → {SKM-GN, SKM-VL} and each brain region to its rat region; Artery - Aorta →
VENACV (imperfect, flagged); all other GTEx tissues one-to-one. Genes are z-scored within species.
Outputs (results/13_gtex/): gene_overlap.csv, panel_survival.csv, accuracy_by_tissue.csv,
confusion_<model>.csv, conformal_transfer.csv, recalibration.csv, native_panel.csv, panel_gene_check.csv;
one REPORT.md section.
"""
from __future__ import annotations

import numpy as np
import pandas as pd

from tfp import cli, config as C, io, models, report
from tfp.splits import assert_no_group_leak
from tfp.transfer import PanelModels, calibrate_models, conformal_transfer, gene_check, gtex_symbols, match_gtex_orthologs, \
    one_to_one_orthologs, per_organ_coverage, save_transfer_scores, score_block, zscore

GTEX_TO_RAT = {
    "Whole Blood": {"BLOOD"}, "Muscle - Skeletal": {"SKM-GN", "SKM-VL"}, "Adipose - Subcutaneous": {"WAT-SC"},
    "Heart - Left Ventricle": {"HEART"}, "Liver": {"LIVER"}, "Kidney - Cortex": {"KIDNEY"}, "Lung": {"LUNG"},
    "Brain - Cortex": {"CORTEX"}, "Brain - Hippocampus": {"HIPPOC"}, "Brain - Hypothalamus": {"HYPOTH"},
    "Colon - Transverse": {"COLON"}, "Small Intestine - Terminal Ileum": {"SMLINT"}, "Spleen": {"SPLEEN"},
    "Adrenal Gland": {"ADRNL"}, "Ovary": {"OVARY"}, "Testis": {"TESTES"}, "Artery - Aorta": {"VENACV"},
}
TISSUE_TO_ORGAN = {t: o for o, ts in GTEX_TO_RAT.items() for t in ts}


def main() -> None:
    ap = cli.common_parser("Human transfer on GTEx")
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--grid", default="20,50")
    ap.add_argument("--cal-frac", type=float, default=0.3)
    ap.add_argument("--recal-donors", default="3,5")
    ap.add_argument("--recal-repeats", type=int, default=20)
    ap.add_argument("--save-scores", action="store_true",
                    help="also write per-sample probabilities of the calibrated models, the calibration scores, the "
                         "per-draw recalibration thresholds, classes.json and organ_map.json")
    args = ap.parse_args()
    cli.banner("13_gtex_transfer", args)
    out = cli.outdir("13_gtex", args.out)
    grid = [int(k) for k in args.grid.split(",")]
    rng = np.random.default_rng(args.seed)

    # ---- data ------------------------------------------------------------------------------
    pheno = io.load_pheno()
    om = io.stack_tissues("TRNSCRPT", source="counts", pheno=pheno, verbose=False)
    Xg = pd.read_csv(C.EXTERNAL_DIR / "gtex_tpm_subset.csv", index_col=0)      # samples × human genes (log2 TPM+1)
    mg = pd.read_csv(C.EXTERNAL_DIR / "gtex_meta.csv", dtype=str).set_index("SAMPID").loc[Xg.index]
    mg["organ"] = mg["SMTSD"]
    mg["group_id"] = mg["donor"]
    mg["stage"] = "adult"
    # ---- orthology (1:1; shared with phase 14) ------------------------------------------------
    orth = one_to_one_orthologs()
    orth = orth[orth["RAT_ENSEMBL_ID"].isin(om.X.columns)]
    by_ens, n_by_ensembl, used_symbol = match_gtex_orthologs(orth, Xg.columns, gtex_symbols())
    rat_genes = by_ens["RAT_ENSEMBL_ID"].tolist()
    hum_genes = by_ens["HUMAN_ORTHOLOG_ENSEMBL_ID"].tolist()
    ov = pd.DataFrame([{"motrpac_genes": om.n_features, "rat_genes_with_1to1_human_ortholog": int(orth.shape[0]),
                        "orthologs_present_in_gtex": len(rat_genes), "matched_by_human_ensembl_id": n_by_ensembl,
                        "matched_by_symbol_fallback": used_symbol,
                        "lost_to_orthology": om.n_features - len(rat_genes), "gtex_genes": Xg.shape[1], "gtex_samples": len(mg),
                        "gtex_donors": mg["donor"].nunique()}])
    ov.to_csv(out / "gene_overlap.csv", index=False)
    sym = io.map_to_gene_symbols(om.X.columns)
    Zm = zscore(om.X[rat_genes])
    Zt = zscore(Xg[hum_genes])
    L_src = om.X[rat_genes].to_numpy(dtype=float)
    L_tgt = Xg[hum_genes].to_numpy(dtype=float)
    y = om.meta["tissue"].astype(str).to_numpy()
    g = om.groups()
    print(f"  1:1 orthologs in GTEx: {len(rat_genes)} of {om.n_features} MoTrPAC genes; GTEx {len(mg)} samples, {mg['donor'].nunique()} donors, "
          f"{mg['organ'].nunique()} tissues")

    # ---- panels on all 50 MoTrPAC animals (ortholog space) ---------------------------------
    pm = PanelModels(L_src, y, rat_genes, grid, args.prefilter)
    Lt_pre, Zt_pre = pm.target_matrices(L_tgt)
    model_names = [f"k{k}" for k in grid] + ["full"]
    # panel survival: the ortholog-space panels vs the all-gene panels of phase 12
    pm_all = PanelModels(om.X.to_numpy(dtype=float), y, list(om.X.columns), grid, args.prefilter, full=False)
    surv = []
    for k in grid:
        full_panel = pm_all.panel_genes(k)
        surv.append({"panel": f"k={k} (all-gene fit)", "n_genes": k, "with_1to1_ortholog_in_gtex": sum(gg in set(rat_genes) for gg in full_panel),
                     "lost": ";".join(str(sym.get(gg, gg)) for gg in full_panel if gg not in set(rat_genes))})
    core_p = C.RESULTS_DIR / "05_panels" / "TRNSCRPT" / "candidate_panel_annotated.csv"
    if core_p.exists():
        core = pd.read_csv(core_p)
        surv.append({"panel": "T5 stable core", "n_genes": len(core), "with_1to1_ortholog_in_gtex": int(core["feature_ID"].isin(rat_genes).sum()),
                     "lost": ";".join(core.loc[~core["feature_ID"].isin(rat_genes), "gene_symbol"].astype(str))})
    surv = pd.DataFrame(surv)
    surv.to_csv(out / "panel_survival.csv", index=False)

    # ---- accuracy per GTEx tissue -----------------------------------------------------------
    acc_rows, conf = [], {}
    for name in model_names:
        pred, ok, _, _ = score_block(pm, name, Lt_pre, Zt_pre, mg, GTEX_TO_RAT)
        sub = mg.assign(pred=pred, correct=ok)
        for organ, d in sub.groupby("organ"):
            acc_rows.append({"model": name, "gtex_tissue": organ, "rat_classes": "/".join(sorted(GTEX_TO_RAT[organ])), "n": len(d),
                             "n_donors": d["donor"].nunique(), "accuracy": float(d["correct"].mean()),
                             "top_prediction": d["pred"].value_counts().index[0], "top_prediction_frac": float(d["pred"].value_counts().iloc[0] / len(d))})
        conf[name] = pd.crosstab(sub["organ"], sub["pred"]).reindex(columns=pm.classes, fill_value=0)
        conf[name].to_csv(out / f"confusion_{name}.csv")
    acc = pd.DataFrame(acc_rows)
    acc.to_csv(out / "accuracy_by_tissue.csv", index=False)
    acc_wide = acc.pivot(index="gtex_tissue", columns="model", values="accuracy").reset_index()
    # convention used everywhere: sample-weighted accuracy (every mapped sample counts once); macro = mean over tissues
    overall = acc.groupby("model").apply(lambda d: np.average(d["accuracy"], weights=d["n"]), include_groups=False)
    macro = acc.groupby("model")["accuracy"].mean()
    pd.DataFrame({"model": overall.index, "accuracy_sample_weighted": overall.to_numpy(), "accuracy_macro_over_tissues": macro.loc[overall.index].to_numpy(),
                  "n_samples": int(len(mg)), "n_tissues": int(mg["organ"].nunique())}).to_csv(out / "accuracy_overall.csv", index=False)

    # ---- conformal transfer -----------------------------------------------------------------
    uniq = np.unique(g)
    rng.shuffle(uniq)
    cal_animals = set(uniq[:max(1, int(round(args.cal_frac * len(uniq))))])
    is_cal = np.array([a in cal_animals for a in g])
    fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
    assert_no_group_leak(om.meta, fit_idx, cal_idx)
    pm_c = PanelModels(L_src[fit_idx], y[fit_idx], rat_genes, grid, args.prefilter)
    classes = pm_c.classes
    Lm_c, Zm_c = pm_c.source_matrices(L_src[cal_idx])
    calib = calibrate_models(pm_c, Lm_c, Zm_c, y[cal_idx], classes, args.alpha, model_names)
    Lt_c, Zt_c = pm_c.target_matrices(L_tgt)
    recal_ns = [int(v) for v in args.recal_donors.split(",") if v.strip()]
    collect = [] if args.save_scores else None
    conf_df, ood, recal = conformal_transfer(calib, classes, Lt_c, Zt_c, mg, GTEX_TO_RAT, args.alpha, "stage", "adult",
                                             recal_ns, args.recal_repeats, rng, model_names, [], collect=collect)
    if args.save_scores:
        save_transfer_scores(out, calib, classes, mg, ["organ", "donor", "stage"], Lt_c, Zt_c, om.meta.iloc[cal_idx],
                             collect, GTEX_TO_RAT, model_names, pm_all=pm, Lp_all=Lt_pre, Zp_all=Zt_pre)
    conf_df.to_csv(out / "conformal_transfer.csv", index=False)
    recal.to_csv(out / "recalibration.csv", index=False)
    # per-tissue coverage with rat calibration: marginal vs Mondrian vs marginal-floor Mondrian
    pt = per_organ_coverage(calib, classes, Lt_c, Zt_c, mg, GTEX_TO_RAT, model_names).rename(columns={"organ": "gtex_tissue", "coverage": "coverage_rat_cal"})
    pt.to_csv(out / "coverage_by_tissue.csv", index=False)
    pt_marg = pt[pt["conformal"] == "marginal"]

    # ---- native GTEx panel under donor-grouped CV -------------------------------------------
    Zt_df = pd.DataFrame(Zt, index=mg.index, columns=hum_genes)
    meta_g = pd.DataFrame({"pid": mg["donor"].to_numpy(), "bid": mg.index, "sex": "unknown", "group": "adult",
                           "tissue": mg["organ"].to_numpy(), "assay": "GTEX"}, index=mg.index)
    om_g = io.OmicsMatrix(Zt_df, meta_g, pd.DataFrame(index=hum_genes), "GTEX", name="gtex")
    curve, sel, preds = models.panel_curve(om_g, label="tissue", grid=[20], kind="logreg_l2", n_splits=5, prefilter=args.prefilter,
                                           quick=True, seed=args.seed, verbose=False, return_predictions=True)
    native_acc = float((preds["y_true"] == preds["y_pred"]).mean())
    imported = {n: float(overall[n]) for n in model_names}
    native = pd.DataFrame([{"native_gtex_panel_k20_accuracy": native_acc, "native_panel_balanced_accuracy": float(curve["balanced_accuracy"].mean()),
                            **{f"imported_{n}_accuracy": v for n, v in imported.items()},
                            "gap_native_minus_imported_k20": native_acc - imported["k20"], "n_donors": mg["donor"].nunique(),
                            "native_panel_genes_fold0": ";".join(sel[sel["fold"] == 0]["feature_ID"].tolist())}])
    native.to_csv(out / "native_panel.csv", index=False)

    # ---- per-gene check of the k=20 panel in GTEx --------------------------------------------
    flags = {}
    for fn in ("stability_k20_annotated.csv", "candidate_panel_annotated.csv"):
        p = C.RESULTS_DIR / "05_panels" / "TRNSCRPT" / fn
        if p.exists():
            for _, r in pd.read_csv(p).iterrows():
                flags.setdefault(r["feature_ID"], {"risk_T7_regulated": bool(r.get("risk_T7_regulated", False)),
                                                  "risk_qc_correlated": bool(r.get("risk_qc_correlated", False))})
    Zm_df = pd.DataFrame(Zm, index=om.X.index, columns=rat_genes)
    rat_to_hum = dict(zip(rat_genes, hum_genes))
    Zt_rat_cols = pd.DataFrame(Zt, index=mg.index, columns=rat_genes)   # same columns, named by the rat gene
    gc = gene_check(pm, pm.panel_genes(20 if 20 in grid else grid[0]), sym, Zm_df, y, Zt_rat_cols, mg, TISSUE_TO_ORGAN, flags)
    gc["human_gene"] = gc["feature_ID"].map(rat_to_hum)
    gc.to_csv(out / "panel_gene_check.csv", index=False)
    testable = gc[gc["fails_in_target"].notna()]
    fails = testable[testable["fails_in_target"].astype(bool)]
    weak = testable[testable["weakened"].astype(bool)]

    # ---- report ---------------------------------------------------------------------------------
    body = [
        f"{om.notes[0]}\n\nGTEx v8 subset: {len(mg)} samples from {mg['donor'].nunique()} donors in {mg['organ'].nunique()} tissues "
        f"(≤ 150 donors per tissue); every GTEx split is grouped on the donor. Orthology: {int(ov['rat_genes_with_1to1_human_ortholog'].iloc[0])} "
        f"MoTrPAC genes have a 1:1 human ortholog, {len(rat_genes)} of those are in GTEx ({n_by_ensembl} by human Ensembl id, {used_symbol} by symbol fallback); "
        f"{int(ov['lost_to_orthology'].iloc[0])} of {om.n_features} genes are lost to orthology. Genes z-scored within species; panels "
        "re-selected in the ortholog space on all 50 MoTrPAC animals.",
        "Panel survival (the all-gene panels of phase 12, and the T5 stable core):\n\n" + report.df_to_md(surv),
        "Tissue mapping: Muscle - Skeletal → {SKM-GN, SKM-VL} super-class; brain regions one-to-one; Artery - Aorta → VENACV is imperfect "
        "and read as a caveat; no out-of-distribution tissue was included.",
        "**Per-tissue accuracy** (super-class scoring):\n\n" + report.df_to_md(acc_wide, floatfmt=".2f")
        + "\n\nOverall, sample-weighted (the convention used throughout; every mapped sample counts once): "
        + ", ".join(f"{n} {v:.3f}" for n, v in imported.items())
        + "; macro over tissues: " + ", ".join(f"{n} {macro[n]:.3f}" for n in model_names),
    ]
    for name in model_names:
        body.append(f"Confusion, {name} (rows = GTEx tissue, columns = predicted rat tissue; zero columns dropped):\n\n"
                    + report.df_to_md(conf[name].loc[:, (conf[name] > 0).any(axis=0)].reset_index()))
    body += [
        f"**Conformal transfer (α = {args.alpha}, LAC; rat calibration on {len(cal_animals)} held-out MoTrPAC animals):**\n\n"
        + report.df_to_md(conf_df, floatfmt=".3f") + "\n\nPer GTEx tissue, rat calibration (marginal):\n\n"
        + report.df_to_md(pt_marg.pivot(index="gtex_tissue", columns="model", values="coverage_rat_cal").reset_index(), floatfmt=".2f")
        + "\n\nPer GTEx tissue, k = 20, marginal vs Mondrian vs marginal-floor Mondrian (rat calibration):\n\n"
        + report.df_to_md(pt[pt["model"] == "k20"].pivot(index="gtex_tissue", columns="conformal", values="coverage_rat_cal").reset_index(), floatfmt=".2f"),
        "Recalibration on GTEx donors (pooled tissues), tested on the remaining donors, vs the rat calibration on the same donors:\n\n"
        + report.df_to_md(recal, floatfmt=".3f"),
        "**Native GTEx 20-gene panel** (donor-grouped 5-fold CV, all 17 tissues) vs the imported rat panels on the same samples:\n\n"
        + report.df_to_md(native.drop(columns=["native_panel_genes_fold0"]).T.reset_index().rename(columns={"index": "quantity", 0: "value"}), floatfmt=".3f"),
        "**Per-gene check of the k = 20 rat panel in GTEx** (effect = mean z in the mapped tissue minus the highest other-tissue mean):\n\n"
        + report.df_to_md(gc, floatfmt=".2f")
        + f"\n\n{len(testable)} of {len(gc)} panel genes have their tissue in GTEx; {len(fails)} fail "
        f"({', '.join(fails['gene_symbol'].astype(str)) or 'none'}) and {len(weak)} are weakened "
        f"({', '.join(weak['gene_symbol'].astype(str)) or 'none'}); of these, "
        f"{int((pd.concat([fails, weak])['risk_T7_regulated'] == True).sum())} carry the training-regulated flag and "
        f"{int((pd.concat([fails, weak])['risk_qc_correlated'] == True).sum())} the QC-correlation flag.",
    ]
    report.add_section("13 · Human transfer on GTEx v8", "\n\n".join(body),
                       params={"alpha": args.alpha, "grid": args.grid, "orthologs": len(rat_genes), "donors": int(mg["donor"].nunique())})
    print(acc_wide.round(2).to_string(index=False)); print(overall.round(3).to_string())
    print(conf_df.round(3).to_string(index=False)); print(recal.round(3).to_string(index=False))
    print(native.drop(columns=["native_panel_genes_fold0"]).T.round(3).to_string()); print(f"wrote {out}")


if __name__ == "__main__":
    main()
