#!/usr/bin/env python
"""Phase 12 — external validation of the rat tissue fingerprint on the rat BodyMap (GSE53960 via
the Bioconductor `bodymapRat` package; see docs/EXTERNAL_VALIDATION.md and R/export_bodymap.R).

Inputs: data/external/bodymap_counts.csv (Ensembl genes × biological samples, technical runs summed),
data/external/bodymap_meta.csv. The BodyMap carries no animal identifier: the animal is the replicate
index within organ × stage × sex (animal_id = sex_stage_replicate), and every BodyMap split is
grouped on it.

Steps: gene matching (Ensembl ids; symbol fallback only if needed) and panel survival; organ mapping
(muscle and brain scored as super-classes, thymus/uterus out-of-distribution); primary test on the
21-week adults with panels of k = 20 and 50 and the full model fit on all 50 MoTrPAC animals, each
gene z-scored within its own dataset; conformal transfer (MoTrPAC calibration, pooled and Mondrian;
recalibration on 3 and 5 BodyMap animals); age shift (2, 6, 104 weeks vs 21); a native BodyMap
panel under animal-grouped CV; and per-gene failure analysis against the T5 risk flags.

Outputs (results/12_bodymap/): gene_overlap.csv, panel_survival.csv, accuracy_by_organ.csv,
confusion_<model>_adult.csv, conformal_transfer.csv, ood_sets.csv, recalibration.csv, age_shift.csv,
native_panel.csv, panel_gene_check.csv; one REPORT.md section.
"""
from __future__ import annotations

import itertools

import numpy as np
import pandas as pd
from tfp import cli, config as C, conformal as cp, io, models, report
from tfp.splits import assert_no_group_leak
from tfp.transfer import PanelModels, calibrate_models, conformal_transfer, per_organ_coverage, save_transfer_scores, \
    score_block, zscore

ORGAN_MAP = {"Adrenal": {"ADRNL"}, "Brain": {"CORTEX", "HIPPOC", "HYPOTH"}, "Heart": {"HEART"}, "Kidney": {"KIDNEY"},
             "Liver": {"LIVER"}, "Lung": {"LUNG"}, "Muscle": {"SKM-GN", "SKM-VL"}, "Spleen": {"SPLEEN"},
             "Testes": {"TESTES"}, "Thymus": None, "Uterus": None}
TISSUE_TO_ORGAN = {t: o for o, ts in ORGAN_MAP.items() if ts for t in ts}


def load_bodymap():
    X = pd.read_csv(C.EXTERNAL_DIR / "bodymap_counts.csv", index_col=0)
    X.index = X.index.astype(str)
    meta = pd.read_csv(C.EXTERNAL_DIR / "bodymap_meta.csv", dtype=str).set_index("sample")
    meta = meta.loc[X.columns]
    meta["stage_weeks"] = meta["stage_weeks"].astype(int)
    meta["sex"] = meta["sex"].map({"F": "female", "M": "male"}).fillna(meta["sex"])
    meta["replicate_index"] = meta["replicate"].astype(str).str.rsplit("_", n=1).str[-1]
    meta["animal_id"] = meta["sex"].str[0] + "_" + meta["stage_weeks"].astype(str) + "_" + meta["replicate_index"]
    return X.T, meta  # samples × genes


def main() -> None:
    ap = cli.common_parser("External validation on the rat BodyMap")
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--grid", default="20,50")
    ap.add_argument("--cal-frac", type=float, default=0.3)
    ap.add_argument("--recal-animals", default="3,5")
    ap.add_argument("--recal-repeats", type=int, default=20)
    ap.add_argument("--save-scores", action="store_true",
                    help="also write per-sample probabilities of the calibrated models, the calibration scores, the "
                         "per-draw recalibration thresholds, classes.json and organ_map.json")
    args = ap.parse_args()
    cli.banner("12_bodymap_validate", args)
    out = cli.outdir("12_bodymap", args.out)
    grid = [int(k) for k in args.grid.split(",")]
    body = []

    # ---- data ------------------------------------------------------------------------------
    pheno = io.load_pheno()
    om = io.stack_tissues("TRNSCRPT", source="counts", pheno=pheno, verbose=False)
    Xb_counts, mb = load_bodymap()
    shared = [g for g in om.X.columns if g in set(Xb_counts.columns)]
    f2g = io.load_feature_to_gene()
    sym = io.map_to_gene_symbols(om.X.columns, f2g)
    ov = {"motrpac_genes": om.n_features, "bodymap_genes": Xb_counts.shape[1], "shared_ensembl": len(shared),
          "motrpac_only": om.n_features - len(shared), "shared_frac_of_motrpac": len(shared) / om.n_features,
          "symbol_fallback_needed": len(shared) / om.n_features < 0.9}
    pd.DataFrame([ov]).to_csv(out / "gene_overlap.csv", index=False)
    Zm = zscore(om.X[shared])                                   # MoTrPAC log2 CPM → z within MoTrPAC
    lb = io.log_cpm(Xb_counts, log=True)                        # BodyMap total-library log2 CPM
    Zb_all = zscore(lb[shared])                                 # z within BodyMap (all 316 samples)
    y = om.meta["tissue"].astype(str).to_numpy()
    g = om.groups()
    mb["group_id"] = mb["animal_id"]
    n_animals_bm = mb["animal_id"].nunique()
    print(f"  shared genes {len(shared)}/{om.n_features}; BodyMap {len(mb)} samples, {n_animals_bm} animal ids, "
          f"organs {mb['organ'].value_counts().to_dict()}")

    # ---- panels on all 50 MoTrPAC animals ---------------------------------------------------
    L_src = om.X[shared].to_numpy(dtype=float)
    pm = PanelModels(L_src, y, shared, grid, args.prefilter)
    Lb_pre, Zb_pre = pm.target_matrices(lb[shared].to_numpy(dtype=float))   # z within BodyMap, all 316 samples
    surv = []
    core = pd.read_csv(C.RESULTS_DIR / "05_panels" / "TRNSCRPT" / "candidate_panel_annotated.csv") if (C.RESULTS_DIR / "05_panels" / "TRNSCRPT" / "candidate_panel_annotated.csv").exists() else None
    for k in grid:
        genes = pm.panel_genes(k)
        surv.append({"panel": f"k={k} (all-animal fit)", "n_genes": k, "present_in_bodymap": sum(gg in set(Xb_counts.columns) for gg in genes),
                     "genes": ";".join(str(sym.get(gg, gg)) for gg in genes)})
    if core is not None:
        surv.append({"panel": "T5 stable core", "n_genes": len(core), "present_in_bodymap": int(core["feature_ID"].isin(Xb_counts.columns).sum()),
                     "genes": ";".join(core["gene_symbol"].astype(str))})
    surv = pd.DataFrame(surv)
    surv.to_csv(out / "panel_survival.csv", index=False)
    model_names = [f"k{k}" for k in grid] + ["full"]

    # ---- primary test: adults, and every stage --------------------------------------------
    acc_rows, conf = [], {}
    for stage in sorted(mb["stage_weeks"].unique()):
        m_s = (mb["stage_weeks"] == stage).to_numpy()
        for name in model_names:
            pred, ok, _, _ = score_block(pm, name, Lb_pre[m_s], Zb_pre[m_s], mb[m_s], ORGAN_MAP)
            sub = mb[m_s].assign(pred=pred, correct=ok)
            for organ, d in sub.groupby("organ"):
                acc_rows.append({"stage_weeks": stage, "model": name, "organ": organ, "n": len(d),
                                 "accuracy": float(d["correct"].mean()) if ORGAN_MAP.get(organ) else np.nan,
                                 "top_prediction": d["pred"].value_counts().index[0],
                                 "top_prediction_frac": float(d["pred"].value_counts().iloc[0] / len(d))})
            if stage == 21:
                conf[name] = pd.crosstab(sub["organ"], sub["pred"]).reindex(columns=pm.classes, fill_value=0)
                conf[name].to_csv(out / f"confusion_{name}_adult.csv")
    acc = pd.DataFrame(acc_rows)
    acc.to_csv(out / "accuracy_by_organ.csv", index=False)
    mapped = acc[acc["accuracy"].notna()]
    acc_wide = mapped[mapped["stage_weeks"] == 21].pivot(index="organ", columns="model", values="accuracy").reset_index()
    overall = mapped.groupby(["stage_weeks", "model"]).apply(lambda d: np.average(d["accuracy"], weights=d["n"]), include_groups=False).unstack("model").reset_index()
    overall.to_csv(out / "age_shift_accuracy.csv", index=False)

    # ---- conformal transfer (shared code path) ----------------------------------------------
    rng = np.random.default_rng(args.seed)
    uniq = np.unique(g)
    rng.shuffle(uniq)
    n_cal = max(1, int(round(args.cal_frac * len(uniq))))
    cal_animals = set(uniq[:n_cal])
    is_cal = np.array([a in cal_animals for a in g])
    fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
    assert_no_group_leak(om.meta, fit_idx, cal_idx)
    pm_c = PanelModels(L_src[fit_idx], y[fit_idx], shared, grid, args.prefilter)
    classes = pm_c.classes
    Lm_c, Zm_c = pm_c.source_matrices(L_src[cal_idx])
    calib = calibrate_models(pm_c, Lm_c, Zm_c, y[cal_idx], classes, args.alpha, model_names)
    Lt_c, Zt_c = pm_c.target_matrices(lb[shared].to_numpy(dtype=float))
    adults = (mb["stage_weeks"] == 21).to_numpy()
    collect = [] if args.save_scores else None
    conf_df, ood, recal = conformal_transfer(calib, classes, Lt_c, Zt_c, mb, ORGAN_MAP, args.alpha, "stage_weeks", 21,
                                             [int(v) for v in args.recal_animals.split(",") if v.strip()], args.recal_repeats, rng,
                                             model_names, ["Thymus", "Uterus"], collect=collect)
    if args.save_scores:
        save_transfer_scores(out, calib, classes, mb, ["organ", "stage_weeks", "sex", "animal_id"], Lt_c, Zt_c,
                             om.meta.iloc[cal_idx], collect, ORGAN_MAP, model_names, pm_all=pm, Lp_all=Lb_pre, Zp_all=Zb_pre)
    conf_df.to_csv(out / "conformal_transfer.csv", index=False)
    ood.to_csv(out / "ood_sets.csv", index=False)
    recal.to_csv(out / "recalibration.csv", index=False)
    # per-organ coverage on the adults: marginal vs Mondrian vs marginal-floor Mondrian (MoTrPAC calibration)
    po = per_organ_coverage(calib, classes, Lt_c[adults], Zt_c[adults], mb[adults], ORGAN_MAP, model_names)
    po.to_csv(out / "coverage_by_organ.csv", index=False)
    po_w = po[po["model"] == "k20"].pivot(index="organ", columns="conformal", values="coverage").reset_index()

    # ---- age shift table ----------------------------------------------------------------------
    age = overall.copy()
    cov_age = conf_df[conf_df["conformal"] == "marginal"].pivot(index="stage_weeks", columns="model", values="coverage_mapped").add_prefix("cov_").reset_index()
    age = age.merge(cov_age, on="stage_weeks")
    for name in model_names:
        ref = age.loc[age["stage_weeks"] == 21, name].iloc[0]
        age[f"acc_drop_{name}"] = ref - age[name]
        refc = age.loc[age["stage_weeks"] == 21, f"cov_{name}"].iloc[0]
        age[f"cov_drop_{name}"] = refc - age[f"cov_{name}"]
    age.to_csv(out / "age_shift.csv", index=False)

    # ---- native BodyMap panel under animal-grouped CV -----------------------------------------
    Zb_df = pd.DataFrame(Zb_all, index=mb.index, columns=shared)
    meta_bm = pd.DataFrame({"pid": mb["animal_id"].to_numpy(), "bid": mb.index, "sex": mb["sex"].to_numpy(),
                            "group": mb["stage_weeks"].astype(str).to_numpy(), "tissue": mb["organ"].to_numpy(), "assay": "BODYMAP"}, index=mb.index)
    om_bm = io.OmicsMatrix(Zb_df, meta_bm, pd.DataFrame(index=shared), "BODYMAP", name="bodymap")
    curve, sel = models.panel_curve(om_bm, label="tissue", grid=[20], kind="logreg_l2", n_splits=5, prefilter=args.prefilter,
                                    quick=True, seed=args.seed, verbose=False)
    native_acc = float(curve["balanced_accuracy"].mean())
    # imported panel on all BodyMap mapped-organ samples (super-class), for the same comparison
    imported = {}
    for name in model_names:
        _, ok, _, _ = score_block(pm, name, Lb_pre, Zb_pre, mb, ORGAN_MAP)
        m_map = np.array([bool(ORGAN_MAP.get(o)) for o in mb["organ"]])
        imported[name] = float(ok[m_map].mean())
    native_map = None
    # native panel accuracy restricted to the 9 mapped organs (from the same CV predictions)
    _, _, preds = models.panel_curve(om_bm, label="tissue", grid=[20], kind="logreg_l2", n_splits=5, prefilter=args.prefilter,
                                     quick=True, seed=args.seed, verbose=False, return_predictions=True)
    pm_ok = preds["y_true"].map(lambda o: bool(ORGAN_MAP.get(o)))
    native_map = float((preds.loc[pm_ok, "y_true"] == preds.loc[pm_ok, "y_pred"]).mean())
    native = pd.DataFrame([{"native_panel_k20_balanced_accuracy_11_organs": native_acc,
                            "native_panel_k20_accuracy_9_mapped_organs": native_map,
                            **{f"imported_{n}_accuracy_9_mapped_organs_all_ages": v for n, v in imported.items()},
                            "gap_native_minus_imported_k20": native_map - imported["k20"], "n_bodymap_animals": n_animals_bm,
                            "native_panel_genes": ";".join(str(sym.get(gg, gg)) for gg in sel[sel["fold"] == 0]["feature_ID"].tolist())}])
    native.to_csv(out / "native_panel.csv", index=False)

    # ---- per-gene failure analysis -------------------------------------------------------------
    genes20 = pm.panel_genes(max(grid[0], 20)) if 20 in grid else pm.panel_genes(grid[0])
    flags = {}
    for fn in ("stability_k20_annotated.csv", "candidate_panel_annotated.csv"):
        p = C.RESULTS_DIR / "05_panels" / "TRNSCRPT" / fn
        if p.exists():
            a = pd.read_csv(p)
            for _, r in a.iterrows():
                flags.setdefault(r["feature_ID"], {"risk_T7_regulated": bool(r.get("risk_T7_regulated", False)),
                                                  "risk_qc_correlated": bool(r.get("risk_qc_correlated", False))})
    sub_ad = mb[adults]
    Zb_ad = pd.DataFrame(Zb_all[adults], index=sub_ad.index, columns=shared)
    Zm_df = pd.DataFrame(Zm, index=om.X.index, columns=shared)
    rows = []
    for gg in genes20:
        j = list(pm.genes_pre).index(gg)
        i_t = int(np.argmax(pm.scores[:, j]))
        marker = pm.classes[i_t]
        organ = TISSUE_TO_ORGAN.get(marker)
        mu_m = Zm_df[gg].groupby(y).mean()
        eff_m = float(mu_m[marker] - mu_m.drop(marker).max())
        row = {"feature_ID": gg, "gene_symbol": sym.get(gg, gg), "marker_tissue": marker, "bodymap_organ": organ or "not in BodyMap",
               "motrpac_effect_z": eff_m}
        if organ:
            mu_b = Zb_ad[gg].groupby(sub_ad["organ"]).mean()
            eff_b = float(mu_b[organ] - mu_b.drop(organ).max())
            row.update({"bodymap_effect_z": eff_b, "bodymap_top_organ": str(mu_b.idxmax()),
                        "fails_in_bodymap": bool(eff_b <= 0), "weakened": bool(0 < eff_b < 0.5 * eff_m)})
        else:
            row.update({"bodymap_effect_z": np.nan, "bodymap_top_organ": None, "fails_in_bodymap": None, "weakened": None})
        row.update(flags.get(gg, {"risk_T7_regulated": None, "risk_qc_correlated": None}))
        rows.append(row)
    gene_check = pd.DataFrame(rows)
    gene_check.to_csv(out / "panel_gene_check.csv", index=False)
    # every panel gene's mean log2 CPM per BodyMap organ × stage (reads the juvenile confusions: Pgk2, Mybph, Hbq1b, Fcrl5)
    jv = []
    for gg in genes20:
        j = list(pm.genes_pre).index(gg)
        marker = pm.classes[int(np.argmax(pm.scores[:, j]))]
        m = pd.DataFrame({"v": lb[gg].to_numpy(dtype=float), "organ": mb["organ"].to_numpy(), "stage_weeks": mb["stage_weeks"].to_numpy()})
        for (organ, stage), dd in m.groupby(["organ", "stage_weeks"]):
            jv.append({"feature_ID": gg, "gene_symbol": sym.get(gg, gg), "marker_tissue": marker, "organ": organ, "stage_weeks": int(stage),
                       "mean_log2_cpm": float(dd["v"].mean()), "n": int(len(dd))})
    pd.DataFrame(jv).to_csv(out / "juvenile_marker_check.csv", index=False)
    testable = gene_check[gene_check["fails_in_bodymap"].notna()]
    fails = testable[testable["fails_in_bodymap"].astype(bool)]
    weak = testable[testable["weakened"].astype(bool)]

    # ---- report ---------------------------------------------------------------------------------
    cov_adult = conf_df[conf_df["stage_weeks"] == 21]
    body = [
        f"{om.notes[0]}\n\nBodyMap: {len(mb)} biological samples from {n_animals_bm} animal ids (technical runs summed; the animal "
        "is the replicate index within organ × stage × sex, GEO/SRA carry no animal id, so every BodyMap split is grouped on "
        "sex_stage_replicate — this is an assumption), organs " + ", ".join(f"{o} {n}" for o, n in mb["organ"].value_counts().items())
        + f"; stages 2/6/21/104 weeks. Gene matching by Ensembl id: {len(shared)} of {om.n_features} MoTrPAC genes present "
        f"({ov['shared_frac_of_motrpac']:.1%}); no symbol fallback needed.",
        "Panel survival:\n\n" + report.df_to_md(surv.drop(columns=["genes"])),
        "Organ mapping: Muscle → {SKM-GN, SKM-VL} and Brain → {CORTEX, HIPPOC, HYPOTH} scored as super-classes; Adrenal, Heart, "
        "Kidney, Liver, Lung, Spleen, Testes one-to-one; Thymus and Uterus out-of-distribution. Each gene is z-scored within its "
        "own dataset (MoTrPAC on 899 samples, BodyMap on 316); classifiers are round-robin panels + L2 logistic regression fit "
        "on all 50 MoTrPAC animals ('full' = the 5,000-gene variance prefilter).",
        "**Primary test, 21-week adults — per-organ accuracy (super-class scoring):**\n\n" + report.df_to_md(acc_wide, floatfmt=".2f")
        + "\n\nOverall (sample-weighted) by stage and model:\n\n" + report.df_to_md(overall, floatfmt=".3f"),
    ]
    for name in model_names:
        body.append(f"Confusion, adults, {name} (rows = BodyMap organ, columns = predicted MoTrPAC tissue; zero columns dropped):\n\n"
                    + report.df_to_md(conf[name].loc[:, (conf[name] > 0).any(axis=0)].reset_index()))
    body += [
        f"**Conformal transfer (α = {args.alpha}, LAC; calibration = {len(cal_animals)} held-out MoTrPAC animals, models refit on the "
        f"other {len(uniq) - len(cal_animals)}):** coverage on the mapped BodyMap adults, empty-set rate, and what the sets hold for "
        "the out-of-distribution organs:\n\n" + report.df_to_md(cov_adult, floatfmt=".3f") + "\n\n" + report.df_to_md(ood, floatfmt=".2f"),
        "Per-organ coverage on the adults with MoTrPAC calibration, k = 20, marginal vs Mondrian vs marginal-floor Mondrian "
        "(each class quantile floored at the marginal one; thymus and uterus have no mapped class, so no coverage):\n\n"
        + report.df_to_md(po_w, floatfmt=".2f"),
        "Recalibration on BodyMap adult animals (pooled organs; the calibration label is the organ's best-scoring mapped tissue), "
        "tested on the remaining adult animals, vs the MoTrPAC calibration on the same animals:\n\n" + report.df_to_md(recal, floatfmt=".3f"),
        "**Age shift** (accuracy and marginal-LAC coverage per stage, drops relative to 21 weeks):\n\n" + report.df_to_md(age, floatfmt=".3f"),
        "**Native BodyMap panel** (k = 20, round-robin + L2 logistic, 5 animal-grouped folds on all 316 samples) vs the imported "
        "MoTrPAC panels scored on the same mapped organs:\n\n" + report.df_to_md(native.drop(columns=["native_panel_genes"]).T.reset_index().rename(columns={"index": "quantity", 0: "value"}), floatfmt=".3f")
        + "\n\nNative panel genes (fold 0): " + native["native_panel_genes"].iloc[0].replace(";", ", "),
        "**Per-gene check of the k = 20 MoTrPAC panel in BodyMap adults** (effect = mean z in the marker organ minus the highest "
        "other-organ mean; 'fails' = effect ≤ 0):\n\n" + report.df_to_md(gene_check, floatfmt=".2f")
        + f"\n\n{len(testable)} of {len(gene_check)} panel genes have their marker tissue in BodyMap; {len(fails)} fail "
        f"({', '.join(fails['gene_symbol'].astype(str)) or 'none'}) and {len(weak)} are weakened to < half their MoTrPAC effect "
        f"({', '.join(weak['gene_symbol'].astype(str)) or 'none'}). Of the failing/weakened genes, "
        f"{int((pd.concat([fails, weak])['risk_T7_regulated'] == True).sum())} carry the training-regulated flag and "
        f"{int((pd.concat([fails, weak])['risk_qc_correlated'] == True).sum())} the QC-correlation flag from T5.",
    ]
    report.add_section("12 · External validation on the rat BodyMap (GSE53960)", "\n\n".join(body),
                       params={"alpha": args.alpha, "grid": args.grid, "shared_genes": len(shared), "source": "bodymapRat"})
    print(acc_wide.round(2).to_string(index=False))
    print(overall.round(3).to_string(index=False))
    print(cov_adult.round(3).to_string(index=False))
    print(po_w.round(2).to_string(index=False))
    print(recal.round(3).to_string(index=False))
    print(native.drop(columns=["native_panel_genes"]).T.round(3).to_string())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
