#!/usr/bin/env python
"""Multiomic Phase 5 — fusion judged by transfer, on matched RNA + protein (Jiang 2020: 182 GTEx samples with both layers).

Both layers are fit on MoTrPAC restricted to the 7 proteomics tissues (the fusion set), so their class probabilities share one
label space: RNA = raw counts → log2 CPM (`tfp.io.stack_tissues(source="counts")`), protein = RII log2 ppm (Phase 1), each in
the 1:1-ortholog space of the corresponding Jiang table (RNA log TPM; cleaned relative protein abundance). Models: the pipeline's
PanelModels (k = 20, 50, full) per layer; late fusion = mean of the two layers' class probabilities; stacked fusion = multinomial
logistic regression on the two layers' out-of-fold MoTrPAC probabilities (animal-grouped folds), fit on MoTrPAC only.
Conformal sets (LAC, α = 0.10) are calibrated on the same 30 % held-out MoTrPAC animals for every model. Every Jiang split is
grouped on the donor. Mapped Jiang tissues: Muscle - Skeletal → SKM-GN, Heart (both sites) → HEART, Lung, Liver, Brain - Cortex;
KIDNEY and WAT-SC have no target; the other tissues are out-of-distribution.
Also: the 19-class RNA fingerprint scored on the same Jiang RNA (the RNA ladder on these samples), the confusion structure of
the two layers on Jiang (hard and probability-weighted off-diagonal correlation), and the local MoTrPAC version.
Outputs: results_multiomic/05_fusion_transfer/.
"""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd
from sklearn.linear_model import LogisticRegression

from tfp import config as C, conformal as cp, io, report
from tfp.splits import assert_no_group_leak, grouped_kfold
from tfp.transfer import PanelModels, one_to_one_orthologs, zscore

ROOT = C.ROOT
EXT = ROOT / "data" / "external_multiomic" / "jiang2020"
RII = ROOT / "results_multiomic" / "01_rii"
OUT = ROOT / "results_multiomic" / "05_fusion_transfer"
PROT7 = ["CORTEX", "HEART", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC"]
JIANG_TO_RAT7 = {"Muscle - Skeletal": "SKM-GN", "Heart - Atrial Appendage": "HEART", "Heart - Left Ventricle": "HEART", "Lung": "LUNG", "Liver": "LIVER", "Brain - Cortex": "CORTEX"}
JIANG_TO_RAT19 = {"Muscle - Skeletal": {"SKM-GN", "SKM-VL"}, "Heart - Atrial Appendage": {"HEART"}, "Heart - Left Ventricle": {"HEART"}, "Lung": {"LUNG"}, "Liver": {"LIVER"}, "Brain - Cortex": {"CORTEX"},
                  "Adrenal Gland": {"ADRNL"}, "Colon - Transverse": {"COLON"}, "Small Intestine - Terminal Ileum": {"SMLINT"}, "Spleen": {"SPLEEN"}, "Testis": {"TESTES"}, "Ovary": {"OVARY"}, "Artery - Aorta": {"VENACV"}}
ALPHA, GRID, SEED, CAL_FRAC = 0.10, [20, 50], C.SEED, 0.3


def f(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{float(v):.{nd}f}"


def key(s):
    return "-".join(str(s).split("-")[:3])


def probs(pm, name, L):
    Lp, Zp = pm.target_matrices(L)
    p, cls = pm.proba(name, Lp, Zp)
    return p, list(cls)


def oof_probs(L, y, g, meta, keys, name, k_grid, n_splits=5):
    """Out-of-fold class probabilities on MoTrPAC (animal-grouped) for one layer."""
    P = np.zeros((len(y), len(np.unique(y))))
    classes = sorted(np.unique(y))
    for tr, te in grouped_kfold(meta, "tissue", n_splits, SEED):
        pm = PanelModels(L[tr], y[tr], keys, k_grid, 5000)
        Lp, Zp = pm.source_matrices(L[te])
        p, cls = pm.proba(name, Lp, Zp)
        P[te] = p[:, [list(cls).index(c) for c in classes]]
    return P, classes


def sets_stats(p, q, classes, y_true_sets):
    sets = cp.predict_sets(p, q, "lac")
    cov = np.array([any(sets[i, classes.index(t)] for t in ts if t in classes) for i, ts in enumerate(y_true_sets)])
    size = sets.sum(axis=1)
    return float(cov.mean()), float((size == 0).mean()), float(size.mean())


def main():
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(SEED)
    pheno = io.load_pheno()
    # ---- MoTrPAC, 7 tissues, both layers ---------------------------------------------------------------------------------
    rna = io.stack_tissues("TRNSCRPT", tissues=PROT7, source="counts", pheno=pheno, verbose=False)
    Xp = pd.read_parquet(RII / "rii_inner_log2ppm.parquet")
    mp = pd.read_csv(RII / "rii_meta.csv", index_col=0, dtype=str).loc[Xp.index]
    feats = pd.read_csv(RII / "rii_features.csv", index_col=0, dtype=str).reindex(Xp.columns)
    p2g = feats["ensembl_gene"].dropna()
    Gp = Xp[p2g.index].T.groupby(p2g.to_numpy()).mean().T
    # ---- Jiang, matched RNA + protein ---------------------------------------------------------------------------------
    rel = pd.read_parquet(EXT / "protein_relative.parquet").set_index("gene.id")
    rel = rel[~rel.index.duplicated(keep="first")]
    cm = pd.read_csv(EXT / "protein_relative_columns.csv", dtype=str).drop_duplicates("column").set_index("column")
    des = pd.read_csv(EXT / "experimental_design.csv", dtype=str).dropna(subset=["GTEx Sample_ID"]).drop_duplicates("GTEx Sample_ID").set_index("GTEx Sample_ID")
    psamp = [c for c in rel.columns if c.startswith("GTEX")]
    rnaJ = pd.read_parquet(EXT / "rna_log_tpm.parquet").set_index("gene.id")
    rnaJ = rnaJ[~rnaJ.index.duplicated(keep="first")]
    rsamp = [c for c in rnaJ.columns if c.startswith("GTEX")]
    rk = {key(s): s for s in rsamp}
    both = [(s, rk[key(s)]) for s in psamp if key(s) in rk]
    mt = pd.DataFrame({"protein_sample": [a for a, _ in both], "rna_sample": [b for _, b in both]}, index=[key(a) for a, _ in both])
    mt["tissue"] = cm["row2"].reindex(mt["protein_sample"]).to_numpy()
    mt["donor"] = des["Individual ID"].reindex(mt["protein_sample"]).to_numpy()
    mt["organ7"] = mt["tissue"].map(JIANG_TO_RAT7)
    mt["organ19"] = [JIANG_TO_RAT19.get(t) for t in mt["tissue"]]
    mapped = mt["organ7"].notna().to_numpy()
    # ---- orthology and matched feature spaces -------------------------------------------------------------------------
    orth = one_to_one_orthologs()
    o_r = orth[orth["RAT_ENSEMBL_ID"].isin(rna.X.columns) & orth["HUMAN_ORTHOLOG_ENSEMBL_ID"].isin(rnaJ.index)].drop_duplicates("RAT_ENSEMBL_ID").drop_duplicates("HUMAN_ORTHOLOG_ENSEMBL_ID")
    comp_p = rel[psamp].notna().mean(axis=1)
    o_p = orth[orth["RAT_ENSEMBL_ID"].isin(Gp.columns) & orth["HUMAN_ORTHOLOG_ENSEMBL_ID"].isin(comp_p[comp_p >= 0.8].index)].drop_duplicates("RAT_ENSEMBL_ID").drop_duplicates("HUMAN_ORTHOLOG_ENSEMBL_ID")
    Lr_src = rna.X[o_r["RAT_ENSEMBL_ID"]].to_numpy(float); yr = rna.meta["tissue"].astype(str).to_numpy(); gr = rna.meta["pid"].astype(str).to_numpy()
    Lp_src = Gp[o_p["RAT_ENSEMBL_ID"]].to_numpy(float); yp = mp["tissue"].astype(str).to_numpy(); gp = mp["pid"].astype(str).to_numpy()
    RJ = rnaJ.loc[o_r["HUMAN_ORTHOLOG_ENSEMBL_ID"], mt["rna_sample"]].T
    PJ = rel.loc[o_p["HUMAN_ORTHOLOG_ENSEMBL_ID"], mt["protein_sample"]].T
    PJ = PJ.fillna(PJ.median(axis=0))
    Lr_tgt, Lp_tgt = RJ.to_numpy(float), PJ.to_numpy(float)
    ov = pd.DataFrame([{"jiang_samples_with_both_layers": len(mt), "jiang_donors": mt["donor"].nunique(), "jiang_tissues": mt["tissue"].nunique(), "mapped_samples_7class": int(mapped.sum()),
                        "mapped_donors_7class": mt.loc[mapped, "donor"].nunique(), "rna_genes_matched": len(o_r), "protein_genes_matched": len(o_p),
                        "motrpac_rna_vials": rna.n_samples, "motrpac_rna_animals": rna.n_animals, "motrpac_prot_vials": len(mp), "motrpac_prot_animals": mp["pid"].nunique()}])
    ov.to_csv(OUT / "overlap.csv", index=False); print(ov.T.to_string())
    classes = PROT7
    y_sets = [{o} if isinstance(o, str) else set() for o in mt["organ7"]]
    # ---- all-animal models per layer ---------------------------------------------------------------------------------
    pm_r = PanelModels(Lr_src, yr, list(o_r["RAT_ENSEMBL_ID"]), GRID, 5000)
    pm_p = PanelModels(Lp_src, yp, list(o_p["RAT_ENSEMBL_ID"]), GRID, 5000)
    model_names = [f"k{k}" for k in GRID] + ["full"]
    # ---- calibration animals: 30 % of the animals present in BOTH layers, held out of both fits -------------------------
    shared_animals = np.array(sorted(set(gr) & set(gp)))
    rng.shuffle(shared_animals)
    cal_animals = set(shared_animals[:max(1, int(round(CAL_FRAC * len(shared_animals))))])
    r_fit, r_cal = np.flatnonzero(~np.isin(gr, list(cal_animals))), np.flatnonzero(np.isin(gr, list(cal_animals)))
    p_fit, p_cal = np.flatnonzero(~np.isin(gp, list(cal_animals))), np.flatnonzero(np.isin(gp, list(cal_animals)))
    assert_no_group_leak(rna.meta, r_fit, r_cal); assert_no_group_leak(mp, p_fit, p_cal)
    # calibration rows must be the same (animal, tissue) in both layers for fusion: align on pid + tissue
    kr = pd.Series([f"{a}|{t}" for a, t in zip(gr[r_cal], yr[r_cal])], index=r_cal)
    kp = pd.Series([f"{a}|{t}" for a, t in zip(gp[p_cal], yp[p_cal])], index=p_cal)
    common_keys = sorted(set(kr) & set(kp))
    r_cal_a = np.array([kr[kr == k].index[0] for k in common_keys]); p_cal_a = np.array([kp[kp == k].index[0] for k in common_keys])
    y_cal = yr[r_cal_a]
    assert (yp[p_cal_a] == y_cal).all()
    pm_rc = PanelModels(Lr_src[r_fit], yr[r_fit], list(o_r["RAT_ENSEMBL_ID"]), GRID, 5000)
    pm_pc = PanelModels(Lp_src[p_fit], yp[p_fit], list(o_p["RAT_ENSEMBL_ID"]), GRID, 5000)
    # stacker training data: out-of-fold probabilities on the FIT animals, aligned on (animal, tissue)
    meta_rf = rna.meta.iloc[r_fit]; meta_pf = mp.iloc[p_fit]
    rows_summary, rows_tissue, conf_store = [], [], {}
    for name in model_names:
        Pr_oof, cls = oof_probs(Lr_src[r_fit], yr[r_fit], gr[r_fit], meta_rf, list(o_r["RAT_ENSEMBL_ID"]), name, GRID)
        Pp_oof, _ = oof_probs(Lp_src[p_fit], yp[p_fit], gp[p_fit], meta_pf, list(o_p["RAT_ENSEMBL_ID"]), name, GRID)
        kr_f = pd.Series([f"{a}|{t}" for a, t in zip(gr[r_fit], yr[r_fit])]); kp_f = pd.Series([f"{a}|{t}" for a, t in zip(gp[p_fit], yp[p_fit])])
        ck = sorted(set(kr_f) & set(kp_f))
        ir = np.array([kr_f[kr_f == k].index[0] for k in ck]); ip = np.array([kp_f[kp_f == k].index[0] for k in ck])
        Xst = np.hstack([Pr_oof[ir], Pp_oof[ip]]); yst = yr[r_fit][ir]
        stacker = LogisticRegression(C=1.0, max_iter=3000).fit(Xst, yst)
        # probabilities: calibration rows (models without cal animals), target rows (models without cal animals → for sets; all-animal → for accuracy)
        def layer_probs(pm_rr, pm_pp, Lr, Lp, is_source):
            mats_r = pm_rr.source_matrices(Lr) if is_source else pm_rr.target_matrices(Lr)
            mats_p = pm_pp.source_matrices(Lp) if is_source else pm_pp.target_matrices(Lp)
            pr, cr = pm_rr.proba(name, *mats_r)
            pp, cpp = pm_pp.proba(name, *mats_p)
            pr = pr[:, [list(cr).index(c) for c in classes]]; pp = pp[:, [list(cpp).index(c) for c in classes]]
            return pr, pp
        pr_cal, pp_cal = layer_probs(pm_rc, pm_pc, Lr_src[r_cal_a], Lp_src[p_cal_a], True)
        pr_t, pp_t = layer_probs(pm_rc, pm_pc, Lr_tgt, Lp_tgt, False)
        pr_t_all, pp_t_all = layer_probs(pm_r, pm_p, Lr_tgt, Lp_tgt, False)
        fusions = {"RNA": (pr_cal, pr_t, pr_t_all), "protein": (pp_cal, pp_t, pp_t_all),
                   "late_mean": ((pr_cal + pp_cal) / 2, (pr_t + pp_t) / 2, (pr_t_all + pp_t_all) / 2)}
        st_cls = list(stacker.classes_)
        def stk(a, b):
            q = stacker.predict_proba(np.hstack([a, b])); return q[:, [st_cls.index(c) for c in classes]]
        fusions["stacked_LR"] = (stk(pr_cal, pp_cal), stk(pr_t, pp_t), stk(pr_t_all, pp_t_all))
        yci = np.array([classes.index(v) for v in y_cal])
        for fname, (pc_, pt_, pt_all) in fusions.items():
            q = cp.conformal_quantile(cp.lac_scores(pc_, yci), ALPHA)
            pred = np.asarray(classes)[pt_all.argmax(axis=1)]
            correct = np.array([p in s for p, s in zip(pred, y_sets)])
            cov, emp, size = sets_stats(pt_[mapped], q, classes, [y_sets[i] for i in np.flatnonzero(mapped)])
            sets_ood = cp.predict_sets(pt_[~mapped], q, "lac")
            rows_summary.append({"model": name, "layer": fname, "n_mapped": int(mapped.sum()), "n_donors": mt.loc[mapped, "donor"].nunique(), "accuracy": float(correct[mapped].mean()),
                                 "coverage_motrpac_cal": cov, "frac_empty": emp, "avg_set_size": size, "ood_frac_empty": float((sets_ood.sum(axis=1) == 0).mean()), "n_cal_rows": len(y_cal),
                                 "cal_accuracy": float(np.mean(np.asarray(classes)[pc_.argmax(axis=1)] == y_cal))})
            sub = mt.assign(pred=pred, correct=correct)[mapped]
            for organ, d in sub.groupby("organ7"):
                rows_tissue.append({"model": name, "layer": fname, "rat_class": organ, "n": len(d), "accuracy": float(d["correct"].mean()), "top_prediction": d["pred"].value_counts().index[0]})
            if name == "k20" or name == "full":
                conf_store[(name, fname)] = (pd.crosstab(pd.Series(mt.loc[mapped, "organ7"].to_numpy(), name="true"), pd.Series(pred[mapped], name="pred")).reindex(index=sorted(set(mt.loc[mapped, "organ7"])), columns=classes, fill_value=0),
                                             pd.DataFrame(pt_all[mapped], columns=classes).groupby(mt.loc[mapped, "organ7"].to_numpy()).mean())
    summ = pd.DataFrame(rows_summary); summ.to_csv(OUT / "fusion_transfer_summary.csv", index=False)
    pd.DataFrame(rows_tissue).to_csv(OUT / "fusion_accuracy_by_tissue.csv", index=False)
    print(summ.round(3).to_string(index=False))
    # ---- confusion structure across layers (Jiang) ----------------------------------------------------------------------
    cs_rows = []
    for name in ("k20", "full"):
        hr, sr = conf_store[(name, "RNA")]; hp, sp = conf_store[(name, "protein")]
        hr.to_csv(OUT / f"confusion_rna_{name}.csv"); hp.to_csv(OUT / f"confusion_protein_{name}.csv"); sr.to_csv(OUT / f"soft_confusion_rna_{name}.csv"); sp.to_csv(OUT / f"soft_confusion_protein_{name}.csv")
        def offdiag(m):
            m = m.div(m.sum(axis=1).replace(0, np.nan), axis=0)
            return np.array([m.loc[i, j] for i in m.index for j in m.columns if i != j])
        a, b = offdiag(hr), offdiag(hp); sa, sb = offdiag(sr), offdiag(sp)
        okh = ~(np.isnan(a) | np.isnan(b)); oks = ~(np.isnan(sa) | np.isnan(sb))
        cs_rows.append({"model": name, "where": "Jiang 2020 (transfer)", "n_offdiag_cells": int(okh.sum()),
                        "pearson_hard_offdiag": float(np.corrcoef(a[okh], b[okh])[0, 1]) if a[okh].std() > 0 and b[okh].std() > 0 else np.nan,
                        "pearson_soft_offdiag": float(np.corrcoef(sa[oks], sb[oks])[0, 1]) if sa[oks].std() > 0 and sb[oks].std() > 0 else np.nan,
                        "per_tissue_accuracy_pearson": float(np.corrcoef(pd.DataFrame(rows_tissue).query("model == @name and layer == 'RNA'").set_index("rat_class")["accuracy"],
                                                                         pd.DataFrame(rows_tissue).query("model == @name and layer == 'protein'").set_index("rat_class")["accuracy"])[0, 1])})
    # local version: MoTrPAC 7 tissues, out-of-fold, soft confusion of RNA vs RII
    for name in ("k20", "full"):
        Pr, cls = oof_probs(Lr_src, yr, gr, rna.meta, list(o_r["RAT_ENSEMBL_ID"]), name, GRID)
        Pp, _ = oof_probs(Lp_src, yp, gp, mp, list(o_p["RAT_ENSEMBL_ID"]), name, GRID)
        sr = pd.DataFrame(Pr, columns=cls).groupby(yr).mean(); sp = pd.DataFrame(Pp, columns=cls).groupby(yp).mean()
        hr = pd.crosstab(pd.Series(yr), pd.Series(np.asarray(cls)[Pr.argmax(axis=1)])).reindex(index=cls, columns=cls, fill_value=0)
        hp = pd.crosstab(pd.Series(yp), pd.Series(np.asarray(cls)[Pp.argmax(axis=1)])).reindex(index=cls, columns=cls, fill_value=0)
        sr.to_csv(OUT / f"local_soft_confusion_rna_{name}.csv"); sp.to_csv(OUT / f"local_soft_confusion_protein_{name}.csv")
        off = lambda m: np.array([m.loc[i, j] for i in m.index for j in m.columns if i != j], dtype=float)
        a, b, sa, sb = off(hr), off(hp), off(sr), off(sp)
        cs_rows.append({"model": name, "where": "MoTrPAC 7 tissues (animal-grouped out-of-fold)", "n_offdiag_cells": len(a),
                        "pearson_hard_offdiag": float(np.corrcoef(a, b)[0, 1]) if a.std() > 0 and b.std() > 0 else np.nan,
                        "pearson_soft_offdiag": float(np.corrcoef(sa, sb)[0, 1]) if sa.std() > 0 and sb.std() > 0 else np.nan,
                        "rna_oof_accuracy": float(np.mean(np.asarray(cls)[Pr.argmax(axis=1)] == yr)), "protein_oof_accuracy": float(np.mean(np.asarray(cls)[Pp.argmax(axis=1)] == yp))})
    cs = pd.DataFrame(cs_rows); cs.to_csv(OUT / "confusion_structure.csv", index=False); print(cs.round(3).to_string(index=False))
    # ---- the 19-class RNA fingerprint on the same Jiang RNA samples (RNA ladder on Jiang) ---------------------------------
    rna19 = io.stack_tissues("TRNSCRPT", source="counts", pheno=pheno, verbose=False)
    o19 = orth[orth["RAT_ENSEMBL_ID"].isin(rna19.X.columns) & orth["HUMAN_ORTHOLOG_ENSEMBL_ID"].isin(rnaJ.index)].drop_duplicates("RAT_ENSEMBL_ID").drop_duplicates("HUMAN_ORTHOLOG_ENSEMBL_ID")
    pm19 = PanelModels(rna19.X[o19["RAT_ENSEMBL_ID"]].to_numpy(float), rna19.meta["tissue"].astype(str).to_numpy(), list(o19["RAT_ENSEMBL_ID"]), GRID, 5000)
    L19 = rnaJ.loc[o19["HUMAN_ORTHOLOG_ENSEMBL_ID"], mt["rna_sample"]].T.to_numpy(float)
    Lp19, Zp19 = pm19.target_matrices(L19)
    m19 = mt["organ19"].notna().to_numpy()
    r19 = []
    for name in model_names:
        p, cls = pm19.proba(name, Lp19, Zp19)
        pred = np.asarray(cls)[p.argmax(axis=1)]
        ok = np.array([bool(s) and pr in s for pr, s in zip(pred, mt["organ19"])])
        r19.append({"model": name, "n_mapped_19class": int(m19.sum()), "n_donors": mt.loc[m19, "donor"].nunique(), "accuracy": float(ok[m19].mean()), "n_tissues_mapped": int(mt.loc[m19, "tissue"].nunique()), "rna_genes": len(o19)})
    r19 = pd.DataFrame(r19); r19.to_csv(OUT / "rna19_on_jiang.csv", index=False); print(r19.round(3).to_string(index=False))
    # ---- README, section, findings ------------------------------------------------------------------------------------
    o = ov.iloc[0]
    k20 = summ[summ["model"] == "k20"].set_index("layer"); full = summ[summ["model"] == "full"].set_index("layer")
    best_single = max(["RNA", "protein"], key=lambda l: k20.loc[l, "accuracy"])
    def robust(fl):
        return all(fl.loc[fu, "accuracy"] >= fl.loc[best_single, "accuracy"] and fl.loc[fu, "coverage_motrpac_cal"] >= fl.loc[best_single, "coverage_motrpac_cal"] and fl.loc[fu, "frac_empty"] <= fl.loc[best_single, "frac_empty"] for fu in ["late_mean", "stacked_LR"] for _ in [0])
    verdict_k20 = {fu: bool(k20.loc[fu, "accuracy"] >= k20.loc[best_single, "accuracy"] and k20.loc[fu, "coverage_motrpac_cal"] >= k20.loc[best_single, "coverage_motrpac_cal"] and k20.loc[fu, "frac_empty"] <= k20.loc[best_single, "frac_empty"]) for fu in ["late_mean", "stacked_LR"]}
    lines = [f"# Phase 5 — fusion judged by transfer (Jiang 2020, matched RNA + protein)", "", f"Built by `scripts/multiomic/05_fusion_transfer.py` on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}; numbers from the CSVs here.", "",
             f"- {int(o['jiang_samples_with_both_layers'])} Jiang samples ({int(o['jiang_donors'])} donors, {int(o['jiang_tissues'])} tissues) have both layers; **{int(o['mapped_samples_7class'])} from {int(o['mapped_donors_7class'])} donors map to the 7-class fingerprint** (CORTEX, HEART, LIVER, LUNG, SKM-GN present; KIDNEY, WAT-SC absent). RNA: {int(o['rna_genes_matched'])} matched genes; protein: {int(o['protein_genes_matched'])} (`overlap.csv`).",
             f"- MoTrPAC 7-tissue source: RNA {int(o['motrpac_rna_vials'])} vials / {int(o['motrpac_rna_animals'])} animals; RII {int(o['motrpac_prot_vials'])} vials / {int(o['motrpac_prot_animals'])} animals; calibration on {len(cal_animals)} animals held out of both layers ({int(summ['n_cal_rows'].iloc[0])} animal × tissue rows).", "",
             "## Accuracy, coverage and empty sets on the same mapped Jiang samples — `fusion_transfer_summary.csv`", "", report.df_to_md(summ, floatfmt=".3f"), "",
             f"Rule (pre-registration): a fusion is 'more robust' only if accuracy AND coverage ≥ the better single layer ({best_single} at k20) and empty-set fraction ≤ it. At k20: late_mean {'YES' if verdict_k20['late_mean'] else 'NO'}, stacked_LR {'YES' if verdict_k20['stacked_LR'] else 'NO'}.", "",
             "Per tissue — `fusion_accuracy_by_tissue.csv`:", "", report.df_to_md(pd.DataFrame(rows_tissue).pivot_table(index=["model", "rat_class", "n"], columns="layer", values="accuracy").reset_index(), floatfmt=".2f"), "",
             "## Confusion structure across layers — `confusion_structure.csv`", "", report.df_to_md(cs, floatfmt=".3f"), "",
             "Off-diagonal correlation of the row-normalised confusion matrices (hard = argmax calls; soft = mean class probability per true tissue). On Jiang the true classes are the 5 mapped tissues; locally all 7.", "",
             "## The 19-class RNA fingerprint on the same Jiang RNA — `rna19_on_jiang.csv`", "", report.df_to_md(r19, floatfmt=".3f"), "",
             "## What this does not show", "", f"- {int(o['mapped_samples_7class'])} samples from {int(o['mapped_donors_7class'])} donors, five classes, two of them tiny (cortex, liver): the fusion comparison is directional, not a precise estimate.",
             "- The protein layer here is the authors' cleaned relative abundance; Phase 3's raw-scale sensitivity gave higher protein accuracy, so the protein arm is conservative.", "", f"_Run time {(time.time() - t0) / 60:.1f} min._"]
    (OUT / "README.md").write_text("\n".join(lines))
    sec = ["- question · is a two-layer fingerprint more robust under the species shift than either layer alone — accuracy AND coverage AND empty-set fraction on the same human samples — and do the layers confuse the same pairs?",
           f"- data · Jiang 2020: {int(o['jiang_samples_with_both_layers'])} samples with matched RNA and protein, {int(o['mapped_samples_7class'])} from {int(o['mapped_donors_7class'])} donors in the 5 mapped classes; MoTrPAC 7-tissue source, RNA {int(o['rna_genes_matched'])} / protein {int(o['protein_genes_matched'])} matched genes (`results_multiomic/05_fusion_transfer/overlap.csv`).",
           "- design · both layers fit on the 7 proteomics tissues (one label space); late fusion = mean probability; stacked LR fit on animal-grouped out-of-fold MoTrPAC probabilities; one calibration set (30 % of animals, held out of both layers) for every model; every Jiang statement on donors.",
           f"- result · at k20 on the same {int(o['mapped_samples_7class'])} samples: RNA accuracy {f(k20.loc['RNA', 'accuracy'])} / coverage {f(k20.loc['RNA', 'coverage_motrpac_cal'])} / empty {f(k20.loc['RNA', 'frac_empty'])}; protein {f(k20.loc['protein', 'accuracy'])} / {f(k20.loc['protein', 'coverage_motrpac_cal'])} / {f(k20.loc['protein', 'frac_empty'])}; "
           f"late mean {f(k20.loc['late_mean', 'accuracy'])} / {f(k20.loc['late_mean', 'coverage_motrpac_cal'])} / {f(k20.loc['late_mean', 'frac_empty'])}; stacked LR {f(k20.loc['stacked_LR', 'accuracy'])} / {f(k20.loc['stacked_LR', 'coverage_motrpac_cal'])} / {f(k20.loc['stacked_LR', 'frac_empty'])}. "
           f"Full models: RNA {f(full.loc['RNA', 'accuracy'])} / {f(full.loc['RNA', 'coverage_motrpac_cal'])}, protein {f(full.loc['protein', 'accuracy'])} / {f(full.loc['protein', 'coverage_motrpac_cal'])}, late mean {f(full.loc['late_mean', 'accuracy'])} / {f(full.loc['late_mean', 'coverage_motrpac_cal'])}, stacked {f(full.loc['stacked_LR', 'accuracy'])} / {f(full.loc['stacked_LR', 'coverage_motrpac_cal'])}. "
           f"Pre-registered robustness rule at k20 (≥ {best_single} on accuracy and coverage, ≤ on empty sets): late mean {'PASS' if verdict_k20['late_mean'] else 'FAIL'}, stacked {'PASS' if verdict_k20['stacked_LR'] else 'FAIL'}. "
           f"Confusion structure on Jiang (k20): off-diagonal Pearson hard {f(cs.iloc[0]['pearson_hard_offdiag'])}, soft {f(cs.iloc[0]['pearson_soft_offdiag'])}; per-tissue accuracy correlation across layers {f(cs.iloc[0]['per_tissue_accuracy_pearson'])}; locally (MoTrPAC out-of-fold, k20) soft {f(cs.iloc[2]['pearson_soft_offdiag'])}. "
           f"The 19-class RNA fingerprint on the same Jiang RNA: k20 {f(r19.set_index('model').loc['k20', 'accuracy'])}, full {f(r19.set_index('model').loc['full', 'accuracy'])} over {int(r19['n_mapped_19class'].iloc[0])} mapped samples in {int(r19['n_tissues_mapped'].iloc[0])} tissues (`rna19_on_jiang.csv`).",
           f"- what it does not show · a precise fusion benefit: {int(o['mapped_samples_7class'])} samples, {int(o['mapped_donors_7class'])} donors, two classes with ≤ 5 samples; coverage under MoTrPAC calibration is near zero for every model, so 'coverage ≥' is a comparison of collapses. The protein arm uses the cleaned relative scale (conservative, see Phase 3).",
           "- files · `results_multiomic/05_fusion_transfer/README.md`."]
    (OUT / "REPORT_SECTION.md").write_text("\n".join(sec) + "\n")
    finds = [{"rank": 5, "text": f"**Fusion under the species shift (same {int(o['mapped_samples_7class'])} Jiang samples, {int(o['mapped_donors_7class'])} donors, 5 tissues).** At k20 accuracy / MoTrPAC-calibrated coverage / empty sets: RNA {f(k20.loc['RNA', 'accuracy'])} / {f(k20.loc['RNA', 'coverage_motrpac_cal'])} / {f(k20.loc['RNA', 'frac_empty'])}, "
                                  f"protein {f(k20.loc['protein', 'accuracy'])} / {f(k20.loc['protein', 'coverage_motrpac_cal'])} / {f(k20.loc['protein', 'frac_empty'])}, late-mean fusion {f(k20.loc['late_mean', 'accuracy'])} / {f(k20.loc['late_mean', 'coverage_motrpac_cal'])} / {f(k20.loc['late_mean', 'frac_empty'])}, "
                                  f"stacked {f(k20.loc['stacked_LR', 'accuracy'])} / {f(k20.loc['stacked_LR', 'coverage_motrpac_cal'])} / {f(k20.loc['stacked_LR', 'frac_empty'])}; pre-registered robustness rule: late mean {'PASS' if verdict_k20['late_mean'] else 'FAIL'}, stacked {'PASS' if verdict_k20['stacked_LR'] else 'FAIL'}. "
                                  f"The two layers' confusion structures on Jiang correlate at {f(cs.iloc[0]['pearson_soft_offdiag'])} (soft off-diagonal, k20). — `results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv`, `confusion_structure.csv`"}]
    (OUT / "FINDINGS.json").write_text(json.dumps(finds, indent=1))
    (OUT / "STATUS.json").write_text(json.dumps({"status": "DONE", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}))
    print(f"wrote {OUT} in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
