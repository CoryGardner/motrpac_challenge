#!/usr/bin/env python
"""Phase 10 — collect the headline numbers from every phase into results/SUMMARY.md.

REPORT.md is the append-only log (every run, in order). SUMMARY.md is the one-page view built
from the latest CSVs, plus an "open questions" checklist to carry to the hackathon.

Conventions used throughout SUMMARY.md: target accuracy on BodyMap and GTEx is SAMPLE-WEIGHTED
(every mapped sample counts once); the macro mean over tissues is kept in the CSVs only.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from motrpac import cli, config as C, report

DURATIONS = ("1w", "2w", "4w", "8w")
VARIANTS = ("marginal", "mondrian", "floored")


def read(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path) if path.exists() else None


def flat(df: pd.DataFrame) -> pd.DataFrame:
    df.columns = ["_".join(str(x) for x in c).rstrip("_") if isinstance(c, tuple) else c for c in df.columns]
    return df


def main() -> None:
    ap = cli.common_parser("Build SUMMARY.md")
    args = ap.parse_args()
    cli.banner("10_make_report", args)
    R = C.RESULTS_DIR
    parts = ["# MoTrPAC hackathon — pre-hackathon summary\n",
             "_Built from results/ CSVs; see REPORT.md for the full log with parameters and figures. Accuracy on the external "
             "targets (BodyMap, GTEx) is sample-weighted throughout: every mapped sample counts once._"]
    if cli.is_synthetic():
        parts.append("\n> **All numbers below come from SYNTHETIC data.** Re-run after the real export.\n")

    inv = read(R / "02_inventory" / "inventory.csv")
    if inv is not None:
        norm = inv[inv["kind"] == "norm"]
        parts.append("## Data\n\n" + report.df_to_md(norm.groupby("assay").agg(
            tissues=("tissue", "nunique"), samples=("n_samples", "sum"), animals=("n_animals", "max"),
            features_median=("n_features", "median")).reset_index())
            + "\n\nUnmatched samples: " + str(int(inv["n_unmatched"].sum())) + ". Fusion tissues: CORTEX, HEART, KIDNEY, LIVER, LUNG, SKM-GN, WAT-SC.")

    eda = read(R / "03_eda" / "readout_variance_pc1-3.csv")
    if eda is not None:
        parts.append("## EDA: what the leading PCs are\n\n" + report.df_to_md(eda[["assay", "PC", "explained", "R2_tissue", "R2_sex", "R2_group"]], floatfmt=".3f")
                     + "\n\nProteomics carries no tissue signal after its per-tissue normalization; transcript batch covariates "
                     "are nested in tissue; TMT channel explains up to 0.91 of a within-tissue proteomics PC after the design is removed.")

    for label, d in (("TRNSCRPT counts, 19 tissues", "TRNSCRPT"), ("METAB, 19 tissues", "METAB"), ("METAB, 9 core tissues", "METAB_core9"),
                     ("PROT, 7 tissues — DIAGNOSTIC ONLY, not a fingerprint", "PROT")):
        s = read(R / "04_baselines" / d / "summary.csv")
        if s is not None:
            cols = [c for c in ("model", "balanced_accuracy_mean", "balanced_accuracy_std", "macro_f1_mean", "fit_seconds_mean") if c in s.columns]
            parts.append(f"## Tissue-fingerprint baselines — {label}\n\n" + report.df_to_md(s[cols], floatfmt=".3f"))
    diag = read(R / "04_baselines" / "PROT" / "diagnostic_accuracy.csv")
    if diag is not None:
        parts.append("PROT diagnostic (fold 0): missingness indicators alone " + f"{diag['missingness_indicators_only'].iloc[0]:.2f}, "
                     f"after per-tissue mean removal {diag['per_tissue_means_removed'].iloc[0]:.2f}, chance {diag['chance_balanced'].iloc[0]:.2f}. "
                     "Plex-structured missingness is a tissue label; proteomics is usable within tissue only.")

    for label, d in (("TRNSCRPT", "TRNSCRPT"), ("METAB core 9", "METAB_core9")):
        c = read(R / "05_panels" / d / "panel_curve.csv")
        if c is not None:
            agg = c.groupby("k")["balanced_accuracy"].agg(["mean", "std"]).reset_index()
            parts.append(f"## Panel size — {label}\n\n" + report.df_to_md(agg, floatfmt=".3f"))
            ann = read(R / "05_panels" / d / "candidate_panel_annotated.csv")
            if ann is not None and len(ann):
                cols = [c for c in ("feature_ID", "gene_symbol", "selection_frequency", "marker_tissue", "effect_size", "risk_T7_regulated", "risk_qc_correlated") if c in ann.columns]
                parts.append("Stable candidate panel (≥ 80% animal-bootstrap selection), with risk flags (annotation, not exclusion):\n\n"
                             + report.df_to_md(ann[cols], floatfmt=".2f")
                             + ("\n\nA `__dup1` feature is the second platform's measurement of a metabolite measured on several platforms "
                                "(provenance in `data/raw/meta/METAB_FEATURES.csv`)." if ann["feature_ID"].astype(str).str.contains("__dup").any() else ""))

    for label, d in (("primary: one vial per animal, α = δ = 0.10, 22 calibration animals", "TRNSCRPT"),
                     ("secondary: pooled vials, α = δ = 0.05", "TRNSCRPT_pooled_a05")):
        cov = read(R / "06_conformal" / d / "coverage.csv")
        dist = read(R / "06_conformal" / d / "certificate_distribution.csv")
        val = read(R / "06_conformal" / d / "certificate_validity.csv")
        if cov is None:
            continue
        cv = cov[(cov.get("conformal", "marginal") == "marginal")] if "conformal" in cov.columns else cov
        cs = cv.groupby(["calibration", "method", "alpha"]).agg(coverage=("coverage", "mean"), avg_set_size=("avg_set_size", "mean")).reset_index()
        parts.append(f"## Conformal — {label}\n\n" + report.df_to_md(cs, floatfmt=".3f")
                     + "\n\nWith one vial per animal the threshold is one of n calibration scores, so coverage moves in steps of "
                     "1/(n+1): with 22 points α = 0.05 takes the largest score; with 12 points α = 0.05 has no finite threshold "
                     "(textbook rank 13 > 12), so every set holds every tissue. The pooled rows sit on 1 − α but treat vials of "
                     "one animal as independent. (Quantile fixed 2026-09-25; results/QUANTILE_FIX_CHANGES.md.)")
        if dist is not None and val is not None:
            k_half = dist.loc[dist["frac_certified"] >= 0.5, "k"].min() if (dist["frac_certified"] >= 0.5).any() else None
            parts.append("Certificate distribution (fraction of fold × repeat splits certifying each k): "
                         + ", ".join(f"k={int(r.k)}: {r.frac_certified:.2f}" for r in dist.itertuples())
                         + f". Smallest k certified in ≥ half of the splits: **{k_half}**; no k certified in "
                         f"{int(val['certified_k'].isna().sum())} of {len(val)} splits; test error ≤ α in "
                         f"{val['test_err_le_alpha'].dropna().astype(bool).mean():.2f} of certified splits.")
    mv = read(R / "06_conformal" / "TRNSCRPT" / "per_tissue_marginal_vs_mondrian_alpha0.1.csv")
    if mv is not None:
        low = mv[mv["coverage_marginal"] < 0.9]
        parts.append("Conditional coverage fix (Mondrian, α = 0.1, pooled): " + "; ".join(
            f"{r.y_true} {r.coverage_marginal:.2f} → {r.coverage_mondrian:.2f} (set size {r.set_size_mondrian:.2f})" for r in low.itertuples()) + "."
            + (" With the marginal floor: " + "; ".join(f"{r.y_true} {r.coverage_floored:.2f} (set size {r.set_size_floored:.2f})" for r in low.itertuples()) + "."
               if "coverage_floored" in mv.columns else ""))
    ad = read(R / "06_conformal" / "TRNSCRPT" / "certificate_alpha_delta_grid.csv")
    if ad is not None:
        parts.append("(α, δ) grid, fraction of splits certifying k = 20: " + ", ".join(
            f"({r.alpha:.2f}, {r.delta:.2f}): {r[3]:.2f}" for r in ad.itertuples(index=False)) + ". Zero-error sizing: 59 animals for (0.05, 0.05), 22 for (0.10, 0.10).")

    # ---- conditional coverage: marginal vs Mondrian vs marginal-floor Mondrian --------------------------------------
    po = read(R / "12_bodymap" / "coverage_by_organ.csv")
    vs = read(R / "08_shift" / "TRNSCRPT" / "shift_per_class_variants.csv")
    st = read(R / "08_shift" / "TRNSCRPT" / "shift_table.csv")
    if mv is not None and "coverage_floored" in mv.columns:
        cc = ["## Conditional coverage — marginal vs Mondrian vs marginal-floor Mondrian (α = 0.1, LAC)\n\n"
              "Mondrian calibrates one quantile per class; the marginal-floor variant raises every class quantile to at least the "
              "marginal quantile, so its sets contain both the marginal and the Mondrian sets and keep both guarantees, paying in set "
              "size on the easy classes.\n\nIn-distribution per tissue (TRNSCRPT full model, pooled calibration, mean over folds):\n\n"
              + report.df_to_md(mv[["y_true", "coverage_marginal", "coverage_mondrian", "coverage_floored", "set_size_marginal", "set_size_mondrian", "set_size_floored"]], floatfmt=".2f")
              + "\n\nLowest per-tissue coverage: " + ", ".join(f"{v} {mv[f'coverage_{v}'].min():.2f} ({mv.loc[mv[f'coverage_{v}'].idxmin(), 'y_true']})" for v in VARIANTS)
              + "; mean set size: " + ", ".join(f"{v} {mv[f'set_size_{v}'].mean():.2f}" for v in VARIANTS) + "."]
        if po is not None:
            p20 = po[po["model"] == "k20"]
            w = p20.pivot(index="organ", columns="conformal", values="coverage")[list(VARIANTS)]
            w2 = p20.pivot(index="organ", columns="conformal", values="avg_set_size")[list(VARIANTS)].add_prefix("set_size_")
            w = w.add_prefix("coverage_").join(w2).reset_index()
            mapped = p20[p20["coverage"].notna()]
            pooled = mapped.groupby("conformal").apply(lambda d: np.average(d["coverage"], weights=d["n"]), include_groups=False)
            cc.append("BodyMap adults, MoTrPAC calibration, k = 20 panel (per organ; thymus and uterus have no mapped class):\n\n"
                      + report.df_to_md(w, floatfmt=".2f") + "\n\nPooled over the mapped adults: "
                      + ", ".join(f"{v} {pooled[v]:.2f}" for v in VARIANTS if v in pooled) + ".")
        if vs is not None:
            sex_splits = [sp for sp in ("train_male_test_female", "train_female_test_male") if sp in set(vs["split"])]
            sub = vs[vs["split"].isin(sex_splits) & (vs["arm"] == "panel_k20") & vs["seen"]]
            if len(sub):
                w = sub.pivot_table(index="y_true", columns=["split", "conformal"], values="coverage")
                w = w[[(sp, v) for sp in sex_splits for v in VARIANTS if (sp, v) in w.columns]]
                w.columns = [f"{'M→F' if sp.startswith('train_male') else 'F→M'} {v}" for sp, v in w.columns]
                txt = "Held-out-sex shift, k = 20 panel calibrated on the source sex (per tissue, seen classes; M→F = trained on males, tested on females):\n\n" + report.df_to_md(w.reset_index(), floatfmt=".2f")
                if st is not None and "coverage_target_seen_floored" in st.columns:
                    ss = st[st["split"].isin(sex_splits) & (st["arm"] == "panel_k20")]
                    txt += "\n\nPooled over the target sex: " + "; ".join(
                        f"{'M→F' if r.split.startswith('train_male') else 'F→M'} marginal {r.coverage_target_seen:.2f} (set size {r.avg_set_size_target:.2f}), "
                        f"Mondrian {r.coverage_target_seen_mondrian:.2f} ({r.avg_set_size_target_mondrian:.2f}), floored {r.coverage_target_seen_floored:.2f} ({r.avg_set_size_target_floored:.2f})"
                        for r in ss.itertuples()) + "."
                cc.append(txt)
        parts.append("\n\n".join(cc))

    fus = read(R / "07_fusion" / "taskA_summary.csv")
    if fus is not None:
        parts.append("## Fusion vs single-omic (tissue identity, TRNSCRPT + METAB; PROT excluded: missingness = tissue label)\n\n" + report.df_to_md(fus, floatfmt=".3f"))
    bv = read(R / "07_fusion" / "taskB_best_vs_null.csv")
    if bv is not None:
        cols = [c for c in ("tissue", "best_arm", "best_auroc", "best_sd", "best_single", "single_auroc", "fusion_minus_single",
                            "fusion_beats_single_by_gt_sd", "null_p95_max_auroc", "p_perm", "best_beats_null_p95") if c in bv.columns]
        parts.append("## Fusion vs single-omic (control vs 8w within tissue, best arm vs the best-of-13-arms permutation null)\n\n"
                     + report.df_to_md(bv[cols], floatfmt=".3f")
                     + f"\n\nFusion beats the best single-omic baseline by more than its fold sd in {int(bv['fusion_beats_single_by_gt_sd'].sum())} of {len(bv)} tissues; "
                     f"the best arm beats the permutation null (95th percentile of the max over arms) in {int(bv['best_beats_null_p95'].sum())} of {len(bv)}. "
                     "Where a fusion arm is listed as best it ties a single-omic arm at the same AUROC; the null's 95th percentile is "
                     f"{bv['null_p95_max_auroc'].min():.2f}–{bv['null_p95_max_auroc'].max():.2f}, so any arm below 0.9 would be indistinguishable from noise. "
                     "Control vs 8w within tissue is separable by each single omic alone (AUROC ≈ 1); plex/channel balance for these "
                     "20-animal subsets is in `results/07_fusion/taskB_batch_balance.csv`.")
    ds = read(R / "07_fusion" / "taskB_duration_summary.csv")
    if ds is not None and ds["duration"].nunique() > 1:
        order = [d for d in DURATIONS if d in set(ds["duration"])]
        wide = ds.pivot(index="tissue", columns="duration", values="best_single_auroc").reset_index()[["tissue"] + order]
        means = ds.groupby("duration")["best_single_auroc"].mean().reindex(order)
        parts.append("**Task B by training duration** (best single-omic AUROC, control vs 1w / 2w / 4w / 8w within tissue, same tuned "
                     "protocol; the 1w, 2w and 4w animals were sacrificed months apart from the controls, the 8w animals on the same days):\n\n"
                     + report.df_to_md(wide, floatfmt=".3f") + "\n\nMean over tissues: " + ", ".join(f"{d} {v:.2f}" for d, v in means.items()) + ".")
    ca = read(R / "07_fusion" / "batch_covariate_auroc.csv")
    if ca is not None:
        w = flat(ca.pivot_table(index="tissue", columns=["feature_set", "model"], values="auroc_mean")).reset_index()
        null = ca[ca["model"] == "logreg"].groupby("tissue")["null_p95_auroc"].max().rename("null_p95_logreg")
        w = w.merge(null.reset_index(), on="tissue")
        parts.append("**Covariates only, no omics** (control vs 8w within tissue; animal-grouped CV AUROC of logistic regression / random "
                     "forest on collection covariates, library-preparation variables, sequencing depth, library QC fractions and TMT "
                     "plex/channel; `null_p95_logreg` = 95th percentile of the within-sex label-permutation null):\n\n"
                     + report.df_to_md(w, floatfmt=".2f"))
    ss_ = read(R / "07_fusion" / "batch_variable_screen_summary.csv")
    if ss_ is not None:
        sep = ss_[ss_["n_tissues_separating"] > 0]
        by = {cat: sep[sep["category"] == cat] for cat in ("design", "outcome", "baseline", "collection", "library", "depth", "qc", "processing")}
        parts.append(f"**Every PHENO / meta variable screened within tissue** ({len(ss_)} variables; a variable 'separates' at AUROC ≥ 0.9 or "
                     "Cramér's V ≥ 0.8 with p < 0.01, or when it is missing in one group only): "
                     + "; ".join(f"{cat} {len(d)}" + (" (" + ", ".join(d["variable"].head(8)) + (", …" if len(d) > 8 else "") + ")" if len(d) and cat not in ("design", "outcome") else "")
                                 for cat, d in by.items())
                     + ". Design variables separate by construction, outcome variables (training log, VO2max, body composition) are the "
                     "training effect itself; only the collection, library, depth, QC and processing rows would indicate a batch.")
    concl = R / "07_fusion" / "batch_conclusion.md"
    if concl.exists():
        paras = [p for p in concl.read_text().split("\n\n") if p.strip()]
        keep = [p for p in paras if p.startswith("**Verdict") or p.startswith("(a)") or p.startswith("**Read")]
        parts.append("**Training, batch, or unresolvable?**\n\n" + "\n\n".join(keep))

    for label, d in (("TRNSCRPT, k = 20", "TRNSCRPT"), ("METAB core 9, k = 10", "METAB_core9")):
        t = read(R / "08_shift" / d / "shift_table.csv")
        if t is not None:
            cols = [c for c in ("split", "arm", "accuracy_all", "coverage_source_id", "coverage_target_seen", "coverage_drop", "cov_target_recal_N3", "unseen_classes") if c in t.columns]
            parts.append(f"## Shift — {label} (α = 0.1, pooled source calibration; `cov_target_recal_N3` = after calibrating on 3 target animals)\n\n"
                         + report.df_to_md(t[cols], floatfmt=".3f"))

    dis = read(R / "09_discordance" / "summary.csv")
    if dis is not None:
        parts.append("## Discordance (PROT restricted to complete features)\n\n" + report.df_to_md(dis, floatfmt=".2f")
                     + "\n\nAt FDR 0.10 no gene is significant in both layers with opposite signs in any tissue except male LIVER (29%) and "
                     "male HEART (3%): discordance here is 'significant in one layer only' (0.5–5.8% of genes). HEART and LIVER, where TMT "
                     "channel encodes sex, show 2–4× more significant genes in males than females; the other tissues are balanced (Q29).")
    ce = read(R / "09_discordance" / "channel_effect.csv")
    if ce is not None:
        parts.append("TMT channel residualized (unsupervised): median ρ before → after: " + "; ".join(
            f"{r.tissue} {r.median_rho_before:.2f} → {r.median_rho_after:.2f}" for r in ce.itertuples()) + ". HEART and LIVER move because channel encodes sex there.")
    au = read(R / "09_discordance" / "prediction_auroc.csv")
    if au is not None and "auroc" in au.columns:
        parts.append(report.df_to_md(au.groupby(["target", "model"])["auroc"].agg(["mean", "std"]).reset_index(), floatfmt=".3f"))
    im = read(R / "09_discordance" / "covariate_importance.csv")
    if im is not None and len(im):
        top = im.groupby(["target", "covariate"])["logreg_drop_one_dAUROC"].mean().reset_index().sort_values(["target", "logreg_drop_one_dAUROC"], ascending=[True, False]).groupby("target").head(3)
        parts.append("Covariates carrying the discordance signal (largest drop-one ΔAUROC): " + "; ".join(
            f"{r.target}: {r.covariate} ({r.logreg_drop_one_dAUROC:+.3f})" for r in top.itertuples()) + ".")

    bm = read(R / "12_bodymap" / "age_shift_accuracy.csv")
    if bm is not None:
        ct = read(R / "12_bodymap" / "conformal_transfer.csv")
        rc = read(R / "12_bodymap" / "recalibration.csv")
        nat = read(R / "12_bodymap" / "native_panel.csv")
        gc = read(R / "12_bodymap" / "panel_gene_check.csv")
        ao = read(R / "12_bodymap" / "accuracy_by_organ.csv")
        jv = read(R / "12_bodymap" / "juvenile_marker_check.csv")
        parts.append("## External validation — rat BodyMap (GSE53960 via bodymapRat; 316 samples, 11 organs, 4 ages)\n\n"
                     "Accuracy on the 9 mapped organs by stage (sample-weighted; super-class scoring for muscle and brain; panels fit on all 50 MoTrPAC "
                     "animals, genes z-scored within each dataset):\n\n" + report.df_to_md(bm, floatfmt=".3f"))
        if ct is not None:
            ad = ct[(ct["stage_weeks"] == 21)]
            parts.append("Conformal transfer on adults with MoTrPAC calibration (α = 0.1, LAC; `floored` = marginal-floor Mondrian):\n\n"
                         + report.df_to_md(ad[["model", "conformal", "coverage_mapped", "frac_empty_mapped", "ood_frac_empty"]], floatfmt=".2f"))
        if rc is not None:
            parts.append("Recalibration on BodyMap adult animals:\n\n" + report.df_to_md(rc[["model", "n_recal", "coverage_recalibrated",
                         "coverage_source_cal_same_test", "set_size_recalibrated"]], floatfmt=".2f"))
        if nat is not None:
            parts.append(f"Native BodyMap 20-gene panel: {nat['native_panel_k20_accuracy_9_mapped_organs'].iloc[0]:.3f} vs imported k = 20 "
                         f"{nat['imported_k20_accuracy_9_mapped_organs_all_ages'].iloc[0]:.3f} (k = 50 {nat['imported_k50_accuracy_9_mapped_organs_all_ages'].iloc[0]:.3f}, "
                         f"full {nat['imported_full_accuracy_9_mapped_organs_all_ages'].iloc[0]:.3f}) on the same mapped organs, all ages.")
        if gc is not None:
            t = gc[gc["fails_in_bodymap"].notna()]
            bad = t[t["fails_in_bodymap"].astype(bool) | t["weakened"].astype(bool)]
            untested = gc[gc["fails_in_bodymap"].isna()]
            parts.append("Panel genes that fail or weaken in BodyMap adults: " + (", ".join(f"{r.gene_symbol} ({r.marker_tissue}"
                         + (", regulated" if r.risk_T7_regulated is True else "") + (", QC-correlated" if r.risk_qc_correlated is True else "") + ")"
                         for r in bad.itertuples()) or "none") + f". The {len(untested)} panel markers whose tissue has no BodyMap organ — "
                         + ", ".join(f"{r.gene_symbol} ({r.marker_tissue})" for r in untested.itertuples())
                         + " — i.e. ovary, adipose, blood, BAT, small intestine, colon and vena cava, are **untested within species**; their only "
                         "external check is the human transfer.")
        if ao is not None:
            j2 = ao[(ao["stage_weeks"] == 2) & (ao["model"] == "k20")].set_index("organ")
            def lvl(gene, organ, stage):
                if jv is None:
                    return np.nan
                m = jv[(jv["gene_symbol"] == gene) & (jv["organ"] == organ) & (jv["stage_weeks"] == stage)]
                return float(m["mean_log2_cpm"].iloc[0]) if len(m) else np.nan
            if {"Testes", "Spleen"} <= set(j2.index):
                parts.append(f"Juveniles (2 weeks): the 20-gene panel calls testes **{j2.loc['Testes', 'top_prediction']}** "
                             f"({j2.loc['Testes', 'accuracy']:.0%} correct) and spleen **{j2.loc['Spleen', 'top_prediction']}** ({j2.loc['Spleen', 'accuracy']:.0%}). "
                             f"The panel's only testis marker, Pgk2, is a spermatid-specific glycolytic isozyme that is not expressed before puberty "
                             f"({lvl('Pgk2', 'Testes', 2):.1f} log2 CPM in 2-week testes vs {lvl('Pgk2', 'Testes', 21):.1f} in adults), while the SKM-VL "
                             f"marker Mybph is higher in the juvenile testis ({lvl('Mybph', 'Testes', 2):.1f} vs {lvl('Mybph', 'Testes', 21):.1f}), so the "
                             f"juvenile testis carries muscle evidence and no testis evidence. The juvenile spleen, still erythropoietic, expresses more of the "
                             f"blood marker Hbq1b ({lvl('Hbq1b', 'Spleen', 2):.1f} vs {lvl('Hbq1b', 'Spleen', 21):.1f}) and less of the spleen marker Fcrl5 "
                             f"({lvl('Fcrl5', 'Spleen', 2):.1f} vs {lvl('Fcrl5', 'Spleen', 21):.1f}), so it is called blood. Both are developmental "
                             "stage effects on single markers: the panel is an adult fingerprint.")
        parts.append("Read: accuracy transfers, the guarantee does not — MoTrPAC-calibrated sets are empty for a third of the BodyMap "
                     "adults because the model is less confident on the new library chemistry; three BodyMap animals recalibrate it. "
                     "Thymus and uterus, absent from MoTrPAC, mostly receive empty sets, the desired out-of-distribution behaviour.")
    gx = read(R / "13_gtex" / "accuracy_by_tissue.csv")
    if gx is not None:
        ov = read(R / "13_gtex" / "gene_overlap.csv")
        ct = read(R / "13_gtex" / "conformal_transfer.csv")
        rc = read(R / "13_gtex" / "recalibration.csv")
        nat = read(R / "13_gtex" / "native_panel.csv")
        gc = read(R / "13_gtex" / "panel_gene_check.csv")
        ovr = read(R / "13_gtex" / "accuracy_overall.csv")
        wide = gx.pivot(index="gtex_tissue", columns="model", values="accuracy").reset_index()
        overall = gx.groupby("model").apply(lambda d: np.average(d["accuracy"], weights=d["n"]), include_groups=False)
        parts.append("## Human transfer — GTEx v8 (" + (f"{int(ov['gtex_samples'].iloc[0])} samples, {int(ov['gtex_donors'].iloc[0])} donors, " if ov is not None else "")
                     + "≤ 150 donors per tissue; 1:1 orthologs only" + (f", {int(ov['orthologs_present_in_gtex'].iloc[0])} genes, {int(ov['lost_to_orthology'].iloc[0])} lost to orthology" if ov is not None else "") + ")\n\n"
                     "**Accuracy convention:** sample-weighted over all mapped GTEx samples (every sample counts once; tissues have 60–150 "
                     "samples, so this is close to the macro mean, which is kept in `results/13_gtex/accuracy_overall.csv`). The same "
                     "convention is used in the questions document and the abstract.\n\n"
                     "Per-tissue accuracy (super-class scoring for muscle; aorta → VENACV is imperfect):\n\n" + report.df_to_md(wide, floatfmt=".2f")
                     + "\n\nOverall, sample-weighted: " + ", ".join(f"{n} {v:.3f}" for n, v in overall.items())
                     + ("; macro over tissues: " + ", ".join(f"{r.model} {r.accuracy_macro_over_tissues:.3f}" for r in ovr.itertuples()) if ovr is not None else ""))
        if ct is not None:
            parts.append("Conformal transfer with rat calibration (α = 0.1, LAC):\n\n" + report.df_to_md(ct[["model", "conformal", "coverage_mapped", "frac_empty_mapped"]], floatfmt=".2f"))
        if rc is not None:
            parts.append("Recalibration on GTEx donors:\n\n" + report.df_to_md(rc[["model", "n_recal", "coverage_recalibrated", "coverage_source_cal_same_test", "set_size_recalibrated"]], floatfmt=".2f"))
        if nat is not None:
            parts.append(f"Native GTEx 20-gene panel (donor-grouped CV, sample-weighted): {nat['native_gtex_panel_k20_accuracy'].iloc[0]:.3f} vs imported k = 20 "
                         f"{nat['imported_k20_accuracy'].iloc[0]:.3f}, k = 50 {nat['imported_k50_accuracy'].iloc[0]:.3f}, full {nat['imported_full_accuracy'].iloc[0]:.3f}.")
        if gc is not None:
            t = gc[gc["fails_in_target"].notna()]
            bad = t[t["fails_in_target"].astype(bool) | t["weakened"].astype(bool)]
            parts.append("Rat panel genes that fail or weaken in GTEx: " + (", ".join(f"{r.gene_symbol} ({r.marker_tissue}"
                         + (", regulated" if r.risk_T7_regulated is True else "") + (", QC-correlated" if r.risk_qc_correlated is True else "") + ")"
                         for r in bad.itertuples()) or "none") + f"; {len(gc) - len(t)} of {len(gc)} have no GTEx tissue to test.")
    cv14 = read(R / "14_transfer" / "motrpac_cv_cost.csv")
    ts = read(R / "14_transfer" / "target_summary.csv")
    if cv14 is not None:
        rc14 = read(R / "14_transfer" / "recalibration.csv")
        cvs = cv14.groupby(["representation", "k", "selector"])["balanced_accuracy"].mean().unstack("selector").reset_index()
        cvs["cost"] = cvs["standard"] - cvs["transfer_aware"]
        parts.append("## Representations and transferability-aware selection\n\n"
                     "MoTrPAC animal-grouped CV balanced accuracy by representation, panel size and selector "
                     "(transfer_aware refuses genes that are training-regulated or |r| > 0.5 with mRNA fraction in their marker tissue):\n\n"
                     + report.df_to_md(cvs, floatfmt=".3f"))
        if ts is not None:
            tw = ts.pivot_table(index=["target", "k", "representation"], columns="selector", values="accuracy").reset_index()
            tw["gain"] = tw["transfer_aware"] - tw["standard"]
            parts.append("Target accuracy on mapped tissues (sample-weighted; BodyMap adults, GTEx on log2 TPM):\n\n" + report.df_to_md(tw, floatfmt=".3f"))
        if rc14 is not None:
            r3 = rc14[rc14["n_recal"] == 3].pivot_table(index=["target", "k", "representation"], columns="selector",
                                                        values=["coverage_source_cal", "coverage_recalibrated", "set_size_recalibrated"]).reset_index()
            parts.append("Coverage at α = 0.1: source-calibrated vs recalibrated on 3 target individuals, with recalibrated set size:\n\n"
                         + report.df_to_md(flat(r3), floatfmt=".2f"))
    # ---- units test: GTEx on log2 CPM from read counts vs log2 TPM ----------------------------------------------------
    ts_c = read(R / "14_transfer_cpm" / "target_summary.csv")
    if ts is not None and ts_c is not None:
        rc_t, rc_c = read(R / "14_transfer" / "recalibration.csv"), read(R / "14_transfer_cpm" / "recalibration.csv")
        ta_t, ta_c = read(R / "14_transfer" / "target_accuracy.csv"), read(R / "14_transfer_cpm" / "target_accuracy.csv")
        a = ts[(ts["target"] == "GTEx") & (ts["selector"] == "standard")][["k", "representation", "accuracy"]].rename(columns={"accuracy": "accuracy_tpm"})
        b = ts_c[(ts_c["target"] == "GTEx") & (ts_c["selector"] == "standard")][["k", "representation", "accuracy"]].rename(columns={"accuracy": "accuracy_cpm"})
        m = a.merge(b, on=["k", "representation"])
        if rc_t is not None and rc_c is not None:
            for name, rcx in (("coverage_source_cal_tpm", rc_t), ("coverage_source_cal_cpm", rc_c)):
                q = rcx[(rcx["target"] == "GTEx") & (rcx["selector"] == "standard") & (rcx["n_recal"] == 3)][["k", "representation", "coverage_source_cal"]]
                m = m.merge(q.rename(columns={"coverage_source_cal": name}), on=["k", "representation"], how="left")
        m = m.sort_values(["representation", "k"])
        txt = ["**Units test — GTEx scored on log2 CPM from read counts (the MoTrPAC / BodyMap unit) instead of log2 TPM** "
               "(standard selector; accuracy sample-weighted; coverage = source-calibrated, α = 0.1):\n\n" + report.df_to_md(m, floatfmt=".3f")]
        if ta_t is not None and ta_c is not None:
            pt = []
            for lab, ta in (("tpm", ta_t), ("cpm", ta_c)):
                q = ta[(ta["target"] == "GTEx") & (ta["selector"] == "standard") & (ta["k"] == 20)].pivot(index="organ", columns="representation", values="accuracy")
                pt.append(q.add_suffix(f"_{lab}"))
            pt = pd.concat(pt, axis=1)
            pt = pt[[c for r in ("rank", "pairs", "zscore") for c in (f"{r}_tpm", f"{r}_cpm") if c in pt.columns]].reset_index()
            txt.append("Per GTEx tissue at k = 20:\n\n" + report.df_to_md(pt, floatfmt=".2f"))
        rp = m[m["representation"].isin(["rank", "pairs"])]
        zs = m[m["representation"] == "zscore"]
        delta = float((rp["accuracy_cpm"] - rp["accuracy_tpm"]).mean())
        gap_cpm = float((zs.set_index("k")["accuracy_cpm"].reindex(rp["k"]).to_numpy() - rp["accuracy_cpm"].to_numpy()).mean()) if len(zs) else np.nan
        if abs(delta) < 0.1 and (np.isnan(gap_cpm) or gap_cpm > 0.1):
            verdict = ("**biology, not units**: matching the unit changes rank / pair accuracy by " + f"{delta:+.2f} on average and leaves it "
                       f"{gap_cpm:.2f} below the z-score representation on the same CPM matrix; the within-sample ordering of the panel genes differs "
                       "between rat and human tissues (and, since the rank representation transfers within species to BodyMap, the difference is "
                       "species, not library chemistry).")
        elif not np.isnan(gap_cpm) and gap_cpm <= 0.05:
            verdict = f"**units**: on CPM the rank / pair representations reach the z-score representation (Δ vs TPM {delta:+.2f})."
        else:
            verdict = f"**partly units**: matching the unit moves rank / pair accuracy by {delta:+.2f} but leaves a {gap_cpm:.2f} gap to the z-score representation."
        txt.append("Read: the cross-species failure of ranks and pairs is " + verdict)
        parts.append("\n\n".join(txt))
    parts.append("""## Open questions to carry to the hackathon

- [ ] The certificate is a distribution, not a number (k = 15–20 in half the splits, none in 30%): do we present the (α, δ) grid or pick (0.10, 0.10)?
- [ ] Conditional coverage: Mondrian fixes the muscle pair and hippocampus at α = 0.1, the marginal-floor variant keeps the marginal guarantee as well — is class-conditional validity the claim, or do we merge SKM-GN/SKM-VL?
- [ ] Shift: the sex shift breaks coverage (transcripts and, more sharply, metabolomics; section on shift above) while accuracy holds; three target animals recalibrate it — is a recalibration protocol part of the deliverable?
- [ ] VENACV: 24% of samples land in BAT with a 20-gene panel; dissection question for the organizers (Q26).
- [ ] Panel risks: Akr1c3 and Ccl25 are training-regulated, Cyp21a1 and Fcrl5 track library mRNA fraction — keep, replace, or report?
- [ ] External checks done (BodyMap within species, GTEx across species; docs/EXTERNAL_VALIDATION.md, docs/GTEX_TRANSFER.md): seven panel markers have no BodyMap organ and are untested within species; is a second rat cohort available?
- [ ] Proteomics is within-tissue only (plex missingness = tissue label; channel = sex in HEART/LIVER): confirm with the organizers (Q21, Q22, Q25).
- [ ] Discordance is one-layer-only significance at FDR 0.10; is that the track's definition (Q30)?
- [ ] Control vs 8w within tissue: batch check done (fusion section) — where it is unresolvable, can the organizers supply library pooling / loading records?
""")
    out = R / "SUMMARY.md"
    out.write_text("\n\n".join(parts))
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
