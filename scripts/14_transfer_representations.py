#!/usr/bin/env python
"""Phase 14 — representation and transferability-aware selection, on both transfer targets.

Three representations of a k-gene panel (tfp.transfer): per-gene z-scores within each dataset,
within-sample ranks of the panel genes, and top-scoring pairs. Two selectors: the standard round
robin, and a transferability-aware one that forbids a tissue from picking a gene that is
training-regulated in that tissue or correlated (|r| > --r-thresh) with the library mRNA fraction
in that tissue (both computed on MoTrPAC only; inside the training fold for the CV cost).
For every selector × representation × k: MoTrPAC animal-grouped CV balanced accuracy (the cost),
per-tissue accuracy on BodyMap adults and on GTEx (super-class scoring), source-calibrated coverage
and recalibrated coverage / set size on 3 and 5 target individuals (what it buys).

Accuracy convention on the targets: sample-weighted over the mapped samples (every sample counts
once), the same convention as phases 12 and 13.

--gtex-matrix cpm scores GTEx on log2 CPM from the read counts (scripts/11_gtex_prepare.py --reads)
instead of log2 TPM: the units test for the rank / pair representations (MoTrPAC and BodyMap are
CPM). --targets, --representations, --selectors and --skip-cv restrict a run; --out redirects it.

Outputs (results/14_transfer/ or --out): motrpac_cv_cost.csv, target_accuracy.csv,
target_summary.csv, recalibration.csv, panels.csv; one REPORT.md section.
"""
from __future__ import annotations

import numpy as np
import pandas as pd
from sklearn.metrics import balanced_accuracy_score

from tfp import cli, config as C, io, report
from tfp.splits import assert_no_group_leak, grouped_kfold
from tfp.transfer import REPRESENTATIONS, PanelModels, calibrate_models, conformal_transfer, exclusion_mask, \
    gtex_symbols, match_gtex_orthologs, one_to_one_orthologs, score_block

BODYMAP_MAP = {"Adrenal": {"ADRNL"}, "Brain": {"CORTEX", "HIPPOC", "HYPOTH"}, "Heart": {"HEART"}, "Kidney": {"KIDNEY"},
               "Liver": {"LIVER"}, "Lung": {"LUNG"}, "Muscle": {"SKM-GN", "SKM-VL"}, "Spleen": {"SPLEEN"},
               "Testes": {"TESTES"}, "Thymus": None, "Uterus": None}
GTEX_MAP = {"Whole Blood": {"BLOOD"}, "Muscle - Skeletal": {"SKM-GN", "SKM-VL"}, "Adipose - Subcutaneous": {"WAT-SC"},
            "Heart - Left Ventricle": {"HEART"}, "Liver": {"LIVER"}, "Kidney - Cortex": {"KIDNEY"}, "Lung": {"LUNG"},
            "Brain - Cortex": {"CORTEX"}, "Brain - Hippocampus": {"HIPPOC"}, "Brain - Hypothalamus": {"HYPOTH"},
            "Colon - Transverse": {"COLON"}, "Small Intestine - Terminal Ileum": {"SMLINT"}, "Spleen": {"SPLEEN"},
            "Adrenal Gland": {"ADRNL"}, "Ovary": {"OVARY"}, "Testis": {"TESTES"}, "Artery - Aorta": {"VENACV"}}
SELECTORS = ("standard", "transfer_aware")
TARGET_NAMES = {"bodymap": "BodyMap adults", "gtex": "GTEx"}


def load_bodymap(shared_with):
    X = pd.read_csv(C.EXTERNAL_DIR / "bodymap_counts.csv", index_col=0)
    X.index = X.index.astype(str)
    meta = pd.read_csv(C.EXTERNAL_DIR / "bodymap_meta.csv", dtype=str).set_index("sample").loc[X.columns]
    meta["stage_weeks"] = meta["stage_weeks"].astype(int)
    meta["replicate_index"] = meta["replicate"].astype(str).str.rsplit("_", n=1).str[-1]
    meta["group_id"] = meta["sex"].astype(str).str[0] + "_" + meta["stage_weeks"].astype(str) + "_" + meta["replicate_index"]
    lb = io.log_cpm(X.T, log=True)
    shared = [g for g in shared_with if g in set(lb.columns)]
    return lb[shared].to_numpy(dtype=float), shared, meta


def load_gtex(motrpac_genes, matrix="tpm"):
    """GTEx samples × 1:1-ortholog genes (rat ids), on log2 TPM or log2 CPM; matching as in phase 13."""
    path = C.EXTERNAL_DIR / ("gtex_cpm_subset.csv" if matrix == "cpm" else "gtex_tpm_subset.csv")
    if not path.exists():
        raise SystemExit(f"{path} is missing — run scripts/11_gtex_prepare.py" + (" --reads <gene_reads.gct.gz>" if matrix == "cpm" else ""))
    Xg = pd.read_csv(path, index_col=0)
    mg = pd.read_csv(C.EXTERNAL_DIR / "gtex_meta.csv", dtype=str).set_index("SAMPID").loc[Xg.index]
    mg["organ"] = mg["SMTSD"]
    mg["group_id"] = mg["donor"]
    mg["stage"] = "adult"
    orth = one_to_one_orthologs()
    orth = orth[orth["RAT_ENSEMBL_ID"].isin(motrpac_genes)]
    by, _, _ = match_gtex_orthologs(orth, Xg.columns, gtex_symbols())
    return Xg[by["HUMAN_ORTHOLOG_ENSEMBL_ID"].tolist()].to_numpy(dtype=float), by["RAT_ENSEMBL_ID"].tolist(), mg


def main() -> None:
    ap = cli.common_parser("Representations and transferability-aware selection")
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--grid", default="20,50")
    ap.add_argument("--cal-frac", type=float, default=0.3)
    ap.add_argument("--recal", default="3,5")
    ap.add_argument("--recal-repeats", type=int, default=20)
    ap.add_argument("--r-thresh", type=float, default=0.5)
    ap.add_argument("--targets", default="bodymap,gtex")
    ap.add_argument("--representations", default=",".join(REPRESENTATIONS))
    ap.add_argument("--selectors", default=",".join(SELECTORS))
    ap.add_argument("--gtex-matrix", default="tpm", choices=["tpm", "cpm"], help="GTEx units: log2 TPM (portal) or log2 CPM from read counts")
    ap.add_argument("--skip-cv", action="store_true", help="skip the MoTrPAC CV cost (targets only)")
    args = ap.parse_args()
    cli.banner("14_transfer_representations", args)
    out = cli.outdir("14_transfer", args.out)
    grid = [int(k) for k in args.grid.split(",")]
    recal_ns = [int(v) for v in args.recal.split(",") if v.strip()]
    n_splits = 3 if args.quick else args.n_splits
    reps = [r.strip() for r in args.representations.split(",") if r.strip()]
    selectors = [s.strip() for s in args.selectors.split(",") if s.strip()]
    target_keys = [t.strip().lower() for t in args.targets.split(",") if t.strip()]
    units = {"bodymap": "log2 CPM", "gtex": "log2 CPM (read counts)" if args.gtex_matrix == "cpm" else "log2 TPM"}

    pheno = io.load_pheno()
    om = io.stack_tissues("TRNSCRPT", source="counts", pheno=pheno, verbose=False)
    y = om.meta["tissue"].astype(str).to_numpy()
    g = om.groups()
    classes = sorted(np.unique(y))
    sym = io.map_to_gene_symbols(om.X.columns)
    # transferability covariates (MoTrPAC only): library mRNA fraction per vial, training-regulated genes per tissue
    tm = pd.read_csv(C.META_DIR / "TRNSCRPT.csv", dtype=str, low_memory=False).drop_duplicates("viallabel").set_index("viallabel")
    qc = pd.to_numeric(tm["pct_mrna"], errors="coerce").reindex(om.X.index).to_numpy()
    reg = pd.read_csv(C.RAW_DIR / "training_regulated_features.csv", dtype=str, usecols=["feature_ID", "assay", "tissue"])
    reg = reg[reg["assay"] == "TRNSCRPT"]
    reg_map = reg.groupby("tissue")["feature_ID"].agg(set).to_dict()
    all_genes = list(om.X.columns)
    L_all = om.X.to_numpy(dtype=float)

    # ---- 1. MoTrPAC CV cost (full gene set; exclusion mask computed on the training fold) ---------
    cv_s, excl = None, pd.DataFrame()
    if not args.skip_cv:
        cv_rows, excl_rows = [], []
        for fold, (tr, te) in enumerate(grouped_kfold(om.meta, "tissue", n_splits, args.seed)):
            assert_no_group_leak(om.meta, tr, te)
            for selector in selectors:
                ex = None if selector == "standard" else exclusion_mask(L_all[tr], y[tr], all_genes, classes, qc=qc[tr], regulated=reg_map, r_thresh=args.r_thresh)
                if ex is not None and fold == 0:
                    excl_rows.append({"selector": selector, "n_gene_tissue_pairs_excluded": int(ex.sum()),
                                      "frac_of_gene_tissue_pairs": float(ex.mean()), "genes_excluded_for_all_tissues": int(ex.all(axis=0).sum())})
                for rep in reps:
                    pm = PanelModels(L_all[tr], y[tr], all_genes, grid, args.prefilter, representation=rep, exclude=ex, full=False)
                    Lp_te, Zp_te = pm.source_matrices(L_all[te])
                    for k in grid:
                        p, cls = pm.proba(f"k{k}", Lp_te, Zp_te)
                        pred = np.asarray(cls)[p.argmax(axis=1)]
                        cv_rows.append({"selector": selector, "representation": rep, "k": k, "fold": fold,
                                        "balanced_accuracy": balanced_accuracy_score(y[te], pred), "n_test_animals": int(om.meta.iloc[te]["pid"].nunique())})
            print(f"  MoTrPAC CV fold {fold} done")
        cv = pd.DataFrame(cv_rows)
        cv.to_csv(out / "motrpac_cv_cost.csv", index=False)
        cv_s = cv.groupby(["representation", "k", "selector"])["balanced_accuracy"].agg(["mean", "std"]).unstack("selector")
        cv_s.columns = [f"{a}_{b}" for a, b in cv_s.columns]
        cv_s = cv_s.reset_index()
        if {"mean_standard", "mean_transfer_aware"} <= set(cv_s.columns):
            cv_s["cost_of_exclusion"] = cv_s["mean_standard"] - cv_s["mean_transfer_aware"]
        excl = pd.DataFrame(excl_rows)

    # ---- 2. targets ------------------------------------------------------------------------------
    targets = {}
    if "bodymap" in target_keys:
        Lb, shared_b, mb = load_bodymap(all_genes)
        targets[TARGET_NAMES["bodymap"]] = (Lb, shared_b, mb, BODYMAP_MAP, "stage_weeks", 21, ["Thymus", "Uterus"])
    if "gtex" in target_keys:
        Lg, rat_g, mg = load_gtex(set(all_genes), args.gtex_matrix)
        targets[TARGET_NAMES["gtex"]] = (Lg, rat_g, mg, GTEX_MAP, "stage", "adult", [])
    acc_rows, sum_rows, recal_rows, panel_rows = [], [], [], []
    for tname, (L_tgt, genes_t, mt, omap, stage_col, primary, ood) in targets.items():
        L_src = om.X[genes_t].to_numpy(dtype=float)
        col = {gg: i for i, gg in enumerate(genes_t)}
        primary_m = (mt[stage_col] == primary).to_numpy()
        # conformal split of the source animals (same for every configuration)
        uniq = np.unique(g)
        rng_split = np.random.default_rng(args.seed)
        rng_split.shuffle(uniq)
        cal_animals = set(uniq[:max(1, int(round(args.cal_frac * len(uniq))))])
        is_cal = np.array([a in cal_animals for a in g])
        fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
        assert_no_group_leak(om.meta, fit_idx, cal_idx)
        ex_ta_all = exclusion_mask(L_src, y, genes_t, classes, qc=qc, regulated=reg_map, r_thresh=args.r_thresh)
        for selector in selectors:
            ex_all = None if selector == "standard" else ex_ta_all
            ex_fit = None if selector == "standard" else exclusion_mask(L_src[fit_idx], y[fit_idx], genes_t, classes, qc=qc[fit_idx], regulated=reg_map, r_thresh=args.r_thresh)
            for rep in reps:
                pm = PanelModels(L_src, y, genes_t, grid, args.prefilter, representation=rep, exclude=ex_all, full=False)
                Lt_pre, Zt_pre = pm.target_matrices(L_tgt)
                for k in grid:
                    pred, ok, _, _ = score_block(pm, f"k{k}", Lt_pre[primary_m], Zt_pre[primary_m], mt[primary_m], omap)
                    sub = mt[primary_m].assign(correct=ok, pred=pred)
                    for organ, d in sub.groupby("organ"):
                        if not omap.get(organ):
                            continue
                        acc_rows.append({"target": tname, "selector": selector, "representation": rep, "k": k, "organ": organ,
                                         "n": len(d), "accuracy": float(d["correct"].mean()), "top_prediction": d["pred"].value_counts().index[0]})
                    mapped = np.array([bool(omap.get(o)) for o in sub["organ"]])
                    per_organ = [float(d["correct"].mean()) for o, d in sub.groupby("organ") if omap.get(o)]
                    sum_rows.append({"target": tname, "selector": selector, "representation": rep, "k": k,
                                     "accuracy": float(ok[mapped].mean()), "accuracy_macro_over_tissues": float(np.mean(per_organ)),
                                     "n": int(mapped.sum())})
                    if rep == "zscore":
                        genes_k = pm.panel_genes(k)
                        refused = []
                        if selector == "standard":
                            for gg in genes_k:
                                j = list(pm.genes_pre).index(gg)
                                marker = pm.classes[int(np.argmax(pm.scores[:, j]))]
                                if ex_ta_all[classes.index(marker), col[gg]]:
                                    refused.append(str(sym.get(gg, gg)))
                        panel_rows.append({"target": tname, "selector": selector, "k": k,
                                           "genes": ";".join(str(sym.get(gg, gg)) for gg in genes_k),
                                           "excluded_under_transfer_aware": ";".join(refused)})
                # conformal: models on the fit animals, calibration on the held-out source animals
                pm_c = PanelModels(L_src[fit_idx], y[fit_idx], genes_t, grid, args.prefilter, representation=rep, exclude=ex_fit, full=False)
                Lm_c, Zm_c = pm_c.source_matrices(L_src[cal_idx])
                names = [f"k{k}" for k in grid]
                calib = calibrate_models(pm_c, Lm_c, Zm_c, y[cal_idx], pm_c.classes, args.alpha, names)
                Lt_c, Zt_c = pm_c.target_matrices(L_tgt)
                conf_df, _, recal = conformal_transfer(calib, pm_c.classes, Lt_c, Zt_c, mt, omap, args.alpha, stage_col, primary,
                                                       recal_ns, args.recal_repeats, np.random.default_rng(args.seed), names, ood)
                src_cov = conf_df[(conf_df[stage_col] == primary) & (conf_df["conformal"] == "marginal")].set_index("model")
                for _, r in recal.iterrows():
                    recal_rows.append({"target": tname, "selector": selector, "representation": rep, "k": int(r["model"][1:]),
                                       "n_recal": int(r["n_recal"]), "coverage_source_cal": float(src_cov.loc[r["model"], "coverage_mapped"]),
                                       "frac_empty_source_cal": float(src_cov.loc[r["model"], "frac_empty_mapped"]),
                                       "coverage_recalibrated": r["coverage_recalibrated"], "set_size_recalibrated": r["set_size_recalibrated"],
                                       "frac_empty_recalibrated": r["frac_empty_recalibrated"]})
                print(f"  {tname:14s} {selector:14s} {rep:6s} acc " + ", ".join(f"k{k}={[s for s in sum_rows if s['target']==tname and s['selector']==selector and s['representation']==rep and s['k']==k][0]['accuracy']:.3f}" for k in grid))
    acc = pd.DataFrame(acc_rows)
    acc.to_csv(out / "target_accuracy.csv", index=False)
    summ = pd.DataFrame(sum_rows)
    summ.to_csv(out / "target_summary.csv", index=False)
    rc = pd.DataFrame(recal_rows)
    rc.to_csv(out / "recalibration.csv", index=False)
    panels = pd.DataFrame(panel_rows)
    panels.to_csv(out / "panels.csv", index=False)

    # ---- report ----------------------------------------------------------------------------------
    summ_w = summ.pivot_table(index=["target", "k", "representation"], columns="selector", values="accuracy").reset_index()
    if {"standard", "transfer_aware"} <= set(summ_w.columns):
        summ_w["gain_from_exclusion"] = summ_w["transfer_aware"] - summ_w["standard"]
    rc_w = rc.pivot_table(index=["target", "k", "representation", "n_recal"], columns="selector",
                          values=["coverage_source_cal", "coverage_recalibrated", "set_size_recalibrated"]).reset_index()
    rc_w.columns = ["_".join(c).rstrip("_") if isinstance(c, tuple) else c for c in rc_w.columns]
    acc20 = acc[acc["k"] == 20].pivot_table(index=["target", "organ"], columns=["selector", "representation"], values="accuracy").reset_index()
    acc20.columns = ["_".join(c).rstrip("_") if isinstance(c, tuple) else c for c in acc20.columns]
    tgt_txt = "; ".join(f"{TARGET_NAMES[k]} on {units[k]}" for k in target_keys if k in TARGET_NAMES)
    body = [
        f"{om.notes[0]}\n\nRepresentations: `zscore` (per-gene z within each dataset), `rank` (within-sample ranks of the k panel genes), "
        "`pairs` (top-scoring pairs: which of two panel genes is higher within the sample, k(k−1)/2 binary features). Selectors: "
        f"`standard` round robin, and `transfer_aware`, which forbids a tissue from picking a gene training-regulated in that tissue "
        f"or with |r| > {args.r_thresh} to the library mRNA fraction in that tissue (MoTrPAC covariates only; computed inside the "
        "training fold for the CV cost, on the fit animals for the calibrated models). Classifier: L2 logistic regression on the "
        f"representation; targets z-scored within themselves; muscle and brain scored as super-classes. Targets and units: {tgt_txt}. "
        "Target accuracy is sample-weighted over the mapped samples (the convention of phases 12 and 13); the macro mean over "
        "tissues is in `target_summary.csv`.",
        ("Exclusion mask on fold 0 of the MoTrPAC CV:\n\n" + (report.df_to_md(excl, floatfmt=".3f") if len(excl) else "_none_")) if not args.skip_cv else "MoTrPAC CV cost skipped (--skip-cv).",
        ("**Cost in MoTrPAC: animal-grouped CV balanced accuracy (mean, sd over folds):**\n\n" + report.df_to_md(cv_s, floatfmt=".3f")) if cv_s is not None else "",
        "**What it buys on the targets: overall accuracy on the mapped tissues (BodyMap 21-week adults; GTEx all samples):**\n\n"
        + report.df_to_md(summ_w, floatfmt=".3f"),
        "Per-tissue accuracy at k = 20:\n\n" + report.df_to_md(acc20, floatfmt=".2f"),
        "**Coverage (α = 0.1, LAC): source-calibrated on held-out MoTrPAC animals vs recalibrated on 3 / 5 target individuals "
        "(pooled tissues), with the recalibrated set size:**\n\n" + report.df_to_md(rc_w, floatfmt=".3f"),
        "Panels (z-score selection is shared by all representations); `excluded_under_transfer_aware` lists the standard panel's "
        "genes that the transfer-aware selector would have refused for their marker tissue:\n\n" + report.df_to_md(panels),
    ]
    title = "14 · Representations and transferability-aware selection (" + ", ".join(TARGET_NAMES[k] for k in target_keys if k in TARGET_NAMES) + ")"
    if "gtex" in target_keys and args.gtex_matrix == "cpm":
        title += " — units test: GTEx on log2 CPM from read counts"
    report.add_section(title, "\n\n".join(b for b in body if b),
                       params={"grid": args.grid, "alpha": args.alpha, "r_thresh": args.r_thresh, "recal": args.recal,
                               "targets": args.targets, "gtex_matrix": args.gtex_matrix, "representations": args.representations,
                               "selectors": args.selectors, "skip_cv": args.skip_cv})
    if cv_s is not None:
        print(cv_s.round(3).to_string(index=False))
    print(summ_w.round(3).to_string(index=False))
    print(rc_w.round(3).to_string(index=False))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
