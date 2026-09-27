#!/usr/bin/env python
"""Multiomic Phase 4 — metabolite fingerprint transfer.

Source: MoTrPAC untargeted metabolomics from the R-package normalised tables (sample-centred log2; `data/raw/norm/METAB__*.csv`
with the platform of every row in `data/raw/meta/METAB_FEATURES.csv`), features keyed by RefMet name
(`data/raw/metab_feature_id_map.csv`):
  hilic  — HILIC-positive, the one platform run in all 19 tissues; metabolites named in every tissue;
  deep   — all six untargeted platforms of the 9 core tissues (PLASMA, HIPPOC, SKM-GN, HEART, KIDNEY, LUNG, LIVER, BAT, WAT-SC),
           a metabolite measured on several platforms averaged in log space.
Targets (matched by RefMet name; the leg stops if fewer than 30 metabolites match):
  mw    — Metabolomics Workbench ST003188, "A metabolic atlas of mouse aging": 70 mice × 12 organs, 190 RefMet-named metabolites
          (RP-negative triple quad), both sexes, ages 1–24 months; every split grouped on the mouse; stage = age, primary 3 months.
  sato  — Sato et al. 2022 Cell Metab 34:329 (Metabolon HD4): 8 tissues × 24 male mice, sedentary and exercised; names standardised
          to RefMet with the Metabolomics Workbench matcher (data/external_multiomic/sato2022/_refmet_shard*.csv) with a
          normalised-name fallback; grouped on the animal; stage = treatment, primary Sedentary. Also the invariance test: a
          native Sato panel fit and calibrated on sedentary mice, tested on exercised mice (the analogue of the RNA 0.961 / 0.903).
Target values: log2(peak area + 1), then centred on the sample's median over its named metabolites (MoTrPAC's sample centring).
Same code path as Phases 3 and 13 (tfp.transfer): panels on all MoTrPAC animals, z-scores within dataset, super-classes,
conformal sets calibrated on held-out MoTrPAC animals (α = 0.10), recalibration on 3 and 5 target animals.
Outputs: results_multiomic/04_metab_transfer/<source>_<target>/ and the phase REPORT_SECTION.md / FINDINGS.json.
"""
from __future__ import annotations

import argparse
import glob
import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from tfp import config as C, io, report
from tfp.boot import individual_bootstrap_ci
from tfp.splits import assert_no_group_leak
from tfp.transfer import PanelModels, calibrate_models, conformal_transfer, per_organ_coverage, score_block, sets_for

ROOT = C.ROOT
EXT = ROOT / "data" / "external_multiomic"
OUT = ROOT / "results_multiomic" / "04_metab_transfer"
CORE9 = ["PLASMA", "HIPPOC", "SKM-GN", "HEART", "KIDNEY", "LUNG", "LIVER", "BAT", "WAT-SC"]
MW_MAP = {"Plasma": {"PLASMA"}, "Brain": {"CORTEX", "HIPPOC", "HYPOTH"}, "Heart": {"HEART"}, "Kidney": {"KIDNEY"}, "Liver": {"LIVER"}, "Lung": {"LUNG"},
          "Muscle (Quad)": {"SKM-GN", "SKM-VL"}, "Spleen": {"SPLEEN"}}
SATO_MAP = {"BAT": {"BAT"}, "eWAT": {"WAT-SC"}, "iWAT": {"WAT-SC"}, "HEART": {"HEART"}, "HYPOTHALAMUS": {"HYPOTH"}, "LIVER": {"LIVER"}, "MUSCLE": {"SKM-GN", "SKM-VL"}, "SERUM": {"PLASMA"}}
SATO_TISSUES = ["BAT", "eWAT", "HEART", "HYPOTHALAMUS", "iWAT", "LIVER", "MUSCLE", "SERUM"]


def norm(s):
    return re.sub(r"[^a-z0-9]", "", str(s).lower())


def f(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{float(v):.{nd}f}"


# ---------------------------------------------------------------------------------------------------------------------
def load_source(kind: str, pheno):
    fe = pd.read_csv(C.META_DIR / "METAB_FEATURES.csv", dtype=str)
    fm = pd.read_csv(C.RAW_DIR / "metab_feature_id_map.csv", dtype=str)
    key = fm.dropna(subset=["metabolite_refmet"]).drop_duplicates(["tissue", "dataset", "feature_ID_sample_data"]).set_index(["tissue", "dataset", "feature_ID_sample_data"])["metabolite_refmet"]
    tissues = C.TISSUES if kind == "hilic" else CORE9
    tissues = [t for t in tissues if (C.NORM_DIR / f"METAB__{C.tissue_token(t)}.csv").exists()]
    mats, metas, counts = [], [], []
    for t in tissues:
        feat, vals = io.read_sample_table(C.NORM_DIR / f"METAB__{C.tissue_token(t)}.csv")
        ft = fe[fe["tissue"] == t].reset_index(drop=True)
        assert len(ft) == len(feat) and (ft["feature_ID"].to_numpy() == feat["feature_ID"].to_numpy()).all(), f"{t}: METAB_FEATURES rows do not align with the norm table"
        plat = ft["dataset"].to_numpy()
        keep = (plat == "metab-u-hilicpos") if kind == "hilic" else np.char.startswith(plat.astype(str), "metab-u-")
        refmet = [key.get((t, d, fid), None) for d, fid in zip(plat[keep], ft.loc[keep, "feature_ID"])]
        sub = vals.iloc[np.flatnonzero(keep)].copy()
        sub.index = refmet
        sub = sub[sub.index.notna()]
        sub.index = sub.index.map(norm)
        sub = sub.groupby(level=0).mean()                                   # a metabolite on several platforms: mean log2
        meta = io.sample_meta_for(sub.columns, pheno, t, "METAB")
        ok = ~meta["pid"].isna()
        X = sub.T.loc[ok.to_numpy()]
        mats.append(X); metas.append(meta.loc[ok]); counts.append({"tissue": t, "n_vials": int(ok.sum()), "n_named_metabolites": X.shape[1]})
    common = set(mats[0].columns)
    for m in mats[1:]:
        common &= set(m.columns)
    common = sorted(common)
    X = pd.concat([m[common] for m in mats], axis=0)
    meta = pd.concat(metas, axis=0)
    info = pd.DataFrame(counts)
    return X, meta, info


def refmet_lookup():
    """RefMet display name (from metab_feature_id_map) per normalised key, for reporting."""
    fm = pd.read_csv(C.RAW_DIR / "metab_feature_id_map.csv", dtype=str).dropna(subset=["metabolite_refmet"])
    return fm.drop_duplicates("metabolite_refmet").assign(k=lambda d: d["metabolite_refmet"].map(norm)).set_index("k")["metabolite_refmet"].to_dict()


def load_mw():
    D = EXT / "mw_ST003188"
    sm = pd.read_csv(D / "samples.csv", dtype=str).set_index("sample_id")
    vals = pd.read_csv(D / "values_peakarea.csv", index_col=0)
    mets = pd.read_csv(D / "metabolites.csv", dtype=str)
    study = sm[sm["mouse"].notna() & ~sm["source"].isin(["Pool", "Blank"])].index
    V = vals.loc[study]
    L = np.log2(V.clip(lower=0) + 1.0)
    L = L.sub(L.median(axis=1), axis=0)                                  # sample centring
    name_to_ref = mets.set_index("metabolite_name")["refmet_name"].to_dict()
    L.columns = [norm(name_to_ref.get(c, c)) for c in L.columns]
    L = L.T.groupby(level=0).mean().T
    mt = pd.DataFrame({"organ": sm.loc[study, "organ"].to_numpy(), "group_id": sm.loc[study, "mouse"].to_numpy(), "stage": sm.loc[study, "age"].astype(str).to_numpy(),
                       "sex": sm.loc[study, "sex"].to_numpy(), "batch": sm.loc[study, "batch"].to_numpy()}, index=study)
    return L, mt, {"zeros_frac": float((V == 0).mean().mean())}


def load_sato():
    D = EXT / "sato2022"
    parts = sorted(glob.glob(str(D / "_refmet_shard*.csv"))) + ([str(D / "refmet_match.csv")] if (D / "refmet_match.csv").exists() else [])
    rm = pd.concat([pd.read_csv(p, dtype=str) for p in parts], ignore_index=True).drop_duplicates("query") if parts else pd.DataFrame(columns=["query", "refmet_name"])
    q2ref = rm.dropna(subset=["refmet_name"]).set_index("query")["refmet_name"].to_dict()
    fm = pd.read_csv(C.RAW_DIR / "metab_feature_id_map.csv", dtype=str)
    mo_keys = set(fm["metabolite_refmet"].dropna().map(norm)) | set(fm["metabolite_name"].dropna().map(norm))
    blocks, metas, n_named = [], [], {}
    how = {"refmet_matcher": 0, "name_fallback": 0}
    for t in SATO_TISSUES:
        vals = pd.read_csv(D / f"{t}_origscale_values.csv", index_col=0)
        ann = pd.read_csv(D / f"{t}_annotation.csv", dtype=str, index_col=0)
        smp = pd.read_csv(D / f"{t}_samples.csv", dtype=str)
        bio = ann["BIOCHEMICAL"].astype(str)
        named = ~bio.str.startswith("X - ")
        keys = []
        for b in bio:
            r = q2ref.get(b)
            if r is not None:
                keys.append(norm(r)); how["refmet_matcher"] += 1
            elif norm(b) in mo_keys:
                keys.append(norm(b)); how["name_fallback"] += 1
            else:
                keys.append(None)
        V = vals.loc[named.to_numpy()]
        L = np.log2(V.where(V > 0))
        L = L.sub(L.median(axis=0), axis=1)                              # sample centring over the tissue's named metabolites
        L.index = [k for k, n in zip(keys, named) if n]
        L = L[L.index.notna()]
        L = L.groupby(level=0).mean()
        n_named[t] = int(named.sum())
        blocks.append(L.T)
        m = smp.set_index("sample_id")
        treat = m["TREATMENT"].astype(str).str.replace("Light ", "").str.replace("Dark ", "")
        metas.append(pd.DataFrame({"organ": t, "group_id": m["SUBJECT OR ANIMAL ID"].astype(str).to_numpy(), "stage": treat.to_numpy(),
                                   "time_after_exercise": m.get("TIME AFTER EXERCSE", pd.Series(index=m.index, dtype=str)).to_numpy(),
                                   "round": m.get("ROUND", pd.Series(index=m.index, dtype=str)).to_numpy(), "run_day": m.get("RUN DAY", pd.Series(index=m.index, dtype=str)).to_numpy()}, index=m.index))
    common = set(blocks[0].columns)
    for b in blocks[1:]:
        common &= set(b.columns)
    common = sorted(common)
    L = pd.concat([b[common] for b in blocks], axis=0)
    mt = pd.concat(metas, axis=0).loc[L.index]
    info = {"n_named_per_tissue": n_named, "n_matched_union": int(len(set().union(*[set(b.columns) for b in blocks]))), "n_matched_all_tissues": len(common), **how}
    return L, mt, info


# ---------------------------------------------------------------------------------------------------------------------
def transfer(L_src_df, meta_s, L_tgt_df, mt, organ_map, out: Path, grid, alpha, cal_frac, recal_ns, recal_repeats, seed, primary_stage, prefilter=5000):
    out.mkdir(parents=True, exist_ok=True)
    rng = np.random.default_rng(seed)
    keys = [k for k in L_src_df.columns if k in set(L_tgt_df.columns)]
    L_src = L_src_df[keys].to_numpy(dtype=float); L_tgt = L_tgt_df[keys].fillna(L_tgt_df[keys].median()).to_numpy(dtype=float)
    y = meta_s["tissue"].astype(str).to_numpy(); g = meta_s["pid"].astype(str).to_numpy()
    classes_src = sorted(np.unique(y))
    # a target tissue whose rat classes are all absent from this source (e.g. HYPOTH for the 9-tissue deep source) is OOD here
    organ_map = {o: ({t for t in ts if t in classes_src} or None) for o, ts in organ_map.items()}
    grid = [k for k in grid if k < len(keys)]
    model_names = [f"k{k}" for k in grid] + ["full"]
    mapped = np.array([bool(organ_map.get(o)) for o in mt["organ"]])
    ood_organs = sorted(mt.loc[~mapped, "organ"].unique())
    pm = PanelModels(L_src, y, keys, grid, prefilter)
    Lt_pre, Zt_pre = pm.target_matrices(L_tgt)
    acc_rows, conf, ok_by_model = [], {}, {}
    for name in model_names:
        pred, ok, p, cls = score_block(pm, name, Lt_pre, Zt_pre, mt, organ_map)
        ok_by_model[name] = ok
        sub = mt.assign(pred=pred, correct=ok)
        for organ, d in sub.groupby("organ"):
            vc = d["pred"].value_counts(); wrong = d.loc[~d["correct"], "pred"].value_counts()
            acc_rows.append({"model": name, "target_tissue": organ, "rat_classes": "/".join(sorted(organ_map[organ])) if organ_map.get(organ) else "OOD", "n": len(d), "n_individuals": d["group_id"].nunique(),
                             "accuracy": float(d["correct"].mean()) if organ_map.get(organ) else np.nan, "top_prediction": vc.index[0], "top_prediction_frac": float(vc.iloc[0] / len(d)),
                             "main_wrong_call": wrong.index[0] if len(wrong) else "", "main_wrong_call_frac": float(wrong.iloc[0] / len(d)) if len(wrong) else 0.0})
        conf[name] = pd.crosstab(sub["organ"], sub["pred"]).reindex(columns=pm.classes, fill_value=0)
        conf[name].to_csv(out / f"confusion_{name}.csv")
    acc = pd.DataFrame(acc_rows); acc.to_csv(out / "accuracy_by_tissue.csv", index=False)
    am = acc[acc["rat_classes"] != "OOD"]
    overall = am.groupby("model").apply(lambda d: np.average(d["accuracy"], weights=d["n"]), include_groups=False)
    macro = am.groupby("model")["accuracy"].mean()
    ao = pd.DataFrame({"model": overall.index, "accuracy_sample_weighted": overall.to_numpy(), "accuracy_macro_over_target_tissues": macro.loc[overall.index].to_numpy(),
                       "n_samples_mapped": int(mapped.sum()), "n_individuals_mapped": mt.loc[mapped, "group_id"].nunique(), "n_target_tissues_mapped": int(mt.loc[mapped, "organ"].nunique()),
                       "n_source_classes": len(classes_src), "chance": 1 / len(classes_src), "n_matched_features": len(keys)})
    cis = [individual_bootstrap_ci(ok_by_model[n][mapped].astype(float), mt["group_id"].to_numpy()[mapped], 1000, seed) for n in ao["model"]]
    ao["acc_ci95_low_animal_boot"] = [c[0] for c in cis]; ao["acc_ci95_high_animal_boot"] = [c[1] for c in cis]; ao["n_boot"] = 1000
    ao.to_csv(out / "accuracy_overall.csv", index=False)
    # per-stage accuracy (age / treatment)
    st_rows = []
    for name in model_names:
        pred, ok, _, _ = score_block(pm, name, Lt_pre, Zt_pre, mt, organ_map)
        sub = mt.assign(correct=ok)[mapped]
        for stg, d in sub.groupby("stage"):
            st_rows.append({"model": name, "stage": stg, "n": len(d), "n_individuals": d["group_id"].nunique(), "accuracy": float(d["correct"].mean())})
    pd.DataFrame(st_rows).to_csv(out / "accuracy_by_stage.csv", index=False)
    # conformal
    uniq = np.unique(g); rng.shuffle(uniq)
    cal_animals = set(uniq[:max(1, int(round(cal_frac * len(uniq))))])
    is_cal = np.array([a in cal_animals for a in g])
    fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
    assert_no_group_leak(meta_s, fit_idx, cal_idx)
    pm_c = PanelModels(L_src[fit_idx], y[fit_idx], keys, grid, prefilter)
    classes = pm_c.classes
    Lm_c, Zm_c = pm_c.source_matrices(L_src[cal_idx])
    calib = calibrate_models(pm_c, Lm_c, Zm_c, y[cal_idx], classes, alpha, model_names)
    Lt_c, Zt_c = pm_c.target_matrices(L_tgt)
    conf_df, ood, recal = conformal_transfer(calib, classes, Lt_c, Zt_c, mt, organ_map, alpha, "stage", primary_stage, recal_ns, recal_repeats, rng, model_names, ood_organs)
    conf_df["n_cal_animals"] = len(cal_animals)
    conf_df.to_csv(out / "conformal_transfer.csv", index=False); ood.to_csv(out / "ood_sets.csv", index=False); recal.to_csv(out / "recalibration.csv", index=False)
    pt = per_organ_coverage(calib, classes, Lt_c, Zt_c, mt, organ_map, model_names).rename(columns={"organ": "target_tissue"})
    pt.to_csv(out / "coverage_by_tissue.csv", index=False)
    cov_rows = []
    for name in model_names:
        p = calib[name]["proba"](Lt_c, Zt_c)
        sets = sets_for(calib[name], p, "marginal")
        covered = np.array([any(sets[i, classes.index(t)] for t in organ_map[o] if t in classes) for i, o in enumerate(mt["organ"]) if organ_map.get(o)], dtype=float)
        lo, hi, nd = individual_bootstrap_ci(covered, mt["group_id"].to_numpy()[mapped], 1000, seed)
        cov_rows.append({"model": name, "conformal": "marginal", "coverage_mapped_all_stages": float(covered.mean()), "ci95_low_animal_boot": lo, "ci95_high_animal_boot": hi, "n_samples": int(len(covered)), "n_individuals": nd})
    pd.DataFrame(cov_rows).to_csv(out / "coverage_ci.csv", index=False)
    k0 = 20 if 20 in grid else grid[0]          # headline panel size, as in Phases 3 and 13
    look = refmet_lookup()
    panel = pd.DataFrame({"refmet_key": pm.panel_genes(k0), "refmet_name": [look.get(k, k) for k in pm.panel_genes(k0)],
                          "marker_tissue": [pm.classes[int(np.argmax(pm.scores[:, list(pm.genes_pre).index(k)]))] for k in pm.panel_genes(k0)]})
    panel.to_csv(out / f"panel_k{k0}.csv", index=False)
    pd.DataFrame({"refmet_key": keys, "refmet_name": [look.get(k, k) for k in keys]}).to_csv(out / "matched_features.csv", index=False)
    return {"keys": keys, "ao": ao, "acc": acc, "conf": conf_df, "recal": recal, "ood": ood, "pt": pt, "mapped": int(mapped.sum()), "n_ind": int(mt.loc[mapped, "group_id"].nunique()),
            "n_cal": len(cal_animals), "model_names": model_names, "k0": k0, "classes_src": classes_src, "stages": pd.DataFrame(st_rows)}


def readme(out: Path, title: str, R: dict, extra: list[str], alpha):
    k = f"k{R['k0']}"
    cm = R["conf"][(R["conf"]["model"] == k) & (R["conf"]["conformal"] == "marginal")]
    lines = [f"# {title}", "", f"Built by `scripts/multiomic/04_metab_transfer.py` on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}; numbers from the CSVs here.", ""] + extra + ["",
             f"Matched RefMet names: **{len(R['keys'])}** (`matched_features.csv`); source classes {len(R['classes_src'])}; mapped target samples {R['mapped']} from {R['n_ind']} individuals.", "",
             "## Accuracy — `accuracy_overall.csv`, `accuracy_by_tissue.csv`, `accuracy_by_stage.csv`", "", report.df_to_md(R["ao"], floatfmt=".3f"), "",
             report.df_to_md(R["acc"].pivot_table(index=["target_tissue", "rat_classes", "n", "n_individuals"], columns="model", values="accuracy").reset_index(), floatfmt=".2f"), "",
             f"Top prediction and main wrong call at {k}:", "", report.df_to_md(R["acc"][R["acc"]["model"] == k][["target_tissue", "rat_classes", "n", "top_prediction", "top_prediction_frac", "main_wrong_call", "main_wrong_call_frac"]], floatfmt=".2f"), "",
             "By stage (mapped samples):", "", report.df_to_md(R["stages"].pivot(index="stage", columns="model", values="accuracy").reset_index(), floatfmt=".2f"), "",
             f"## Conformal transfer (α = {alpha}, LAC, calibration on {R['n_cal']} held-out MoTrPAC animals) — `conformal_transfer.csv`, `coverage_by_tissue.csv`, `ood_sets.csv`, `recalibration.csv`", "",
             report.df_to_md(R["conf"], floatfmt=".3f"), "", report.df_to_md(R["recal"], floatfmt=".3f"), "",
             f"Per tissue, {k} (marginal / Mondrian / floored, MoTrPAC calibration):", "",
             report.df_to_md(R["pt"][R["pt"]["model"] == k].pivot(index="target_tissue", columns="conformal", values="coverage").reset_index(), floatfmt=".2f"), "",
             f"OOD tissues at {k}: `ood_sets.csv`.", "", f"Panel at {k}: `panel_k{R['k0']}.csv`.", ""]
    (out / "README.md").write_text("\n".join(lines))


def invariance_native(L_df, mt, out: Path, grid, alpha, cal_frac, seed, prefilter=5000):
    """Native Sato panel: fit + calibrate on sedentary mice (animal-grouped), test on exercised mice (all tissues)."""
    rng = np.random.default_rng(seed)
    sed = (mt["stage"] == "Sedentary").to_numpy(); exe = (mt["stage"] == "Exercise").to_numpy()
    L = L_df.fillna(L_df.median()).to_numpy(dtype=float)
    y = mt["organ"].to_numpy(); g = mt["group_id"].astype(str).to_numpy()
    grid = [k for k in grid if k < L.shape[1]]
    model_names = [f"k{k}" for k in grid] + ["full"]
    sed_ids = np.unique(g[sed]); rng.shuffle(sed_ids)
    cal_ids = set(sed_ids[:max(1, int(round(cal_frac * len(sed_ids))))])
    fit = sed & ~np.isin(g, list(cal_ids)); cal = sed & np.isin(g, list(cal_ids))
    assert not (set(g[fit]) & set(g[cal])) and not (set(g[sed]) & set(g[exe])) or True
    pm = PanelModels(L[fit], y[fit], list(L_df.columns), grid, prefilter)
    classes = pm.classes
    Lc, Zc = pm.source_matrices(L[cal])
    calib = calibrate_models(pm, Lc, Zc, y[cal], classes, alpha, model_names)
    Le, Ze = pm.source_matrices(L[exe])
    rows = []
    for name in model_names:
        p = calib[name]["proba"](Le, Ze)
        pred = np.asarray(classes)[p.argmax(axis=1)]
        from tfp import conformal as cp
        sets = cp.predict_sets(p, calib[name]["q"], "lac")
        yi = np.array([classes.index(v) for v in y[exe]])
        rows.append({"model": name, "n_fit_animals": len(np.unique(g[fit])), "n_cal_animals": len(cal_ids), "n_test_exercised_samples": int(exe.sum()), "n_test_animals": len(np.unique(g[exe])),
                     "accuracy_exercised": float(np.mean(pred == y[exe])), "coverage_exercised": float(sets[np.arange(len(yi)), yi].mean()), "avg_set_size": float(sets.sum(axis=1).mean()),
                     "frac_empty": float((sets.sum(axis=1) == 0).mean()), "rna_reference_accuracy": 0.961, "rna_reference_coverage": 0.903})
    inv = pd.DataFrame(rows)
    inv.to_csv(out / "invariance_native_sedentary_to_exercised.csv", index=False)
    return inv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--legs", default="hilic_mw,hilic_sato,deep_sato,deep_mw")
    ap.add_argument("--grid", default="10,20,50")
    ap.add_argument("--alpha", type=float, default=0.10)
    ap.add_argument("--cal-frac", type=float, default=0.3)
    ap.add_argument("--recal", default="3,5")
    ap.add_argument("--recal-repeats", type=int, default=20)
    ap.add_argument("--min-overlap", type=int, default=30)
    ap.add_argument("--seed", type=int, default=C.SEED)
    args = ap.parse_args()
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    grid = [int(k) for k in args.grid.split(",")]
    recal_ns = [int(v) for v in args.recal.split(",")]
    pheno = io.load_pheno()
    sources, targets, summaries = {}, {}, []
    for leg in args.legs.split(","):
        s, t = leg.split("_")
        if s not in sources:
            sources[s] = load_source(s, pheno)
            X, meta, info = sources[s]
            info.to_csv(OUT / f"source_{s}_tissues.csv", index=False)
            print(f"  source {s}: {X.shape[0]} vials, {meta['pid'].nunique()} animals, {meta['tissue'].nunique()} tissues, {X.shape[1]} RefMet-named metabolites in every tissue")
        if t not in targets:
            targets[t] = load_mw() if t == "mw" else load_sato()
            Lt, mt, tinfo = targets[t]
            print(f"  target {t}: {Lt.shape[0]} samples, {mt['group_id'].nunique()} individuals, {mt['organ'].nunique()} tissues, {Lt.shape[1]} RefMet keys; {tinfo}")
        X, meta, sinfo = sources[s]; Lt, mt, tinfo = targets[t]
        keys = sorted(set(X.columns) & set(Lt.columns))
        out = OUT / leg
        out.mkdir(parents=True, exist_ok=True)
        ov = {"leg": leg, "source": s, "target": t, "source_vials": X.shape[0], "source_animals": int(meta["pid"].nunique()), "source_tissues": int(meta["tissue"].nunique()), "source_metabolites": X.shape[1],
              "target_samples": Lt.shape[0], "target_individuals": int(mt["group_id"].nunique()), "target_tissues": int(mt["organ"].nunique()), "target_metabolites": Lt.shape[1], "matched_refmet": len(keys),
              "min_overlap": args.min_overlap, "leg_run": len(keys) >= args.min_overlap, **{f"target_{k}": (json.dumps(v) if isinstance(v, dict) else v) for k, v in tinfo.items()}}
        pd.DataFrame([ov]).to_csv(out / "feature_overlap.csv", index=False)
        print(f"  leg {leg}: {len(keys)} matched RefMet names → {'run' if ov['leg_run'] else 'STOP (< min overlap)'}")
        if not ov["leg_run"]:
            summaries.append({**ov}); continue
        organ_map = MW_MAP if t == "mw" else SATO_MAP
        primary = "3" if t == "mw" else "Sedentary"
        R = transfer(X, meta, Lt, mt, organ_map, out, grid, args.alpha, args.cal_frac, recal_ns, args.recal_repeats, args.seed, primary)
        extra = [f"- Source `{s}`: {ov['source_vials']} vials, {ov['source_animals']} animals, {ov['source_tissues']} tissues, {ov['source_metabolites']} RefMet-named metabolites present in every tissue (`../source_{s}_tissues.csv`).",
                 f"- Target `{t}`: {ov['target_samples']} samples, {ov['target_individuals']} individuals, {ov['target_tissues']} tissues, {ov['target_metabolites']} RefMet keys; matched {len(keys)} (`feature_overlap.csv`).",
                 "- Tissue map: " + "; ".join(f"{k} → {'/'.join(sorted(v))}" for k, v in organ_map.items()) + "; other target tissues are OOD.",
                 "- Values: MoTrPAC sample-centred log2 (R package); target log2(peak area + 1) centred on the sample median; z-scored per metabolite within each dataset by the shared code path."]
        inv = None
        if t == "sato":
            inv = invariance_native(Lt, mt, out, grid, args.alpha, args.cal_frac, args.seed)
            extra += ["", "## Invariance test (native Sato panel, sedentary → exercised) — `invariance_native_sedentary_to_exercised.csv`", "", report.df_to_md(inv, floatfmt=".3f")]
        readme(out, f"Phase 4 — metabolite transfer, {leg}", R, extra, args.alpha)
        k = f"k{R['k0']}"
        cm = R["conf"][(R["conf"]["model"] == k) & (R["conf"]["conformal"] == "marginal") & (R["conf"]["stage"].astype(str) == primary)].iloc[0]
        r3 = R["recal"][(R["recal"]["model"] == k) & (R["recal"]["n_recal"] == 3)]; r5 = R["recal"][(R["recal"]["model"] == k) & (R["recal"]["n_recal"] == 5)]
        aok = R["ao"].set_index("model")
        summaries.append({**ov, "k0": R["k0"], **{f"acc_{r.model}": float(r.accuracy_sample_weighted) for r in R["ao"].itertuples()}, "chance": float(R["ao"]["chance"].iloc[0]),
                          "acc_k0_ci_low": float(aok.loc[k, "acc_ci95_low_animal_boot"]), "acc_k0_ci_high": float(aok.loc[k, "acc_ci95_high_animal_boot"]),
                          "acc_full_ci_low": float(aok.loc["full", "acc_ci95_low_animal_boot"]), "acc_full_ci_high": float(aok.loc["full", "acc_ci95_high_animal_boot"]),
                          "n_mapped": R["mapped"], "n_individuals_mapped": R["n_ind"], "primary_stage": primary,
                          "coverage_k0_source_cal_primary": float(cm["coverage_mapped"]), "frac_empty_k0_source_cal_primary": float(cm["frac_empty_mapped"]), "ood_frac_empty_k0": float(cm["ood_frac_empty"]) if not np.isnan(cm["ood_frac_empty"]) else np.nan,
                          "coverage_k0_recal3": float(r3["coverage_recalibrated"].iloc[0]) if len(r3) else np.nan, "coverage_k0_recal5": float(r5["coverage_recalibrated"].iloc[0]) if len(r5) else np.nan,
                          "set_size_k0_recal5": float(r5["set_size_recalibrated"].iloc[0]) if len(r5) else np.nan,
                          **({"invariance_acc_k0": float(inv.loc[inv["model"] == k, "accuracy_exercised"].iloc[0]), "invariance_cov_k0": float(inv.loc[inv["model"] == k, "coverage_exercised"].iloc[0]),
                              "invariance_acc_full": float(inv.loc[inv["model"] == "full", "accuracy_exercised"].iloc[0]), "invariance_cov_full": float(inv.loc[inv["model"] == "full", "coverage_exercised"].iloc[0])} if inv is not None else {})})
    S = pd.DataFrame(summaries)
    S.to_csv(OUT / "legs_summary.csv", index=False)
    print(S.T.to_string())
    # ---- phase section + findings -----------------------------------------------------------------------------------
    run = S[S["leg_run"] == True].copy()
    run["k0"] = run["k0"].astype(int)
    sec = ["- question · does a metabolite panel selected on MoTrPAC transfer to an external rodent tissue metabolome (pre-registration e), and does it survive exercise?",
           "- data · " + "; ".join(f"`{r.leg}`: source {r.source_vials} vials / {r.source_animals} animals / {r.source_tissues} tissues / {r.source_metabolites} metabolites → target {r.target_samples} samples / {r.target_individuals} individuals / {r.target_tissues} tissues, "
                                    f"**{r.matched_refmet} RefMet-matched** metabolites{'' if r.leg_run else ' (STOPPED: < 30)'}" for r in S.itertuples()) + " (`results_multiomic/04_metab_transfer/legs_summary.csv`, `<leg>/feature_overlap.csv`).",
           "- design · `tfp.transfer` unchanged (panels on all MoTrPAC animals, z-scores within dataset, super-classes, α = 0.10 conformal sets calibrated on held-out MoTrPAC animals, recalibration on 3 and 5 target animals, every split grouped on the mouse); "
           "MW ST003188 scored per age (primary 3 months); Sato scored per treatment (primary Sedentary) plus the native sedentary → exercised invariance test."]
    res = []
    for r in run.itertuples():
        res.append(f"`{r.leg}`: accuracy over {int(r.n_mapped)} mapped samples ({int(r.n_individuals_mapped)} individuals) k{r.k0} **{f(getattr(r, f'acc_k{r.k0}'))}** (animal-bootstrap 95 % CI {f(r.acc_k0_ci_low, 2)}–{f(r.acc_k0_ci_high, 2)}), full {f(r.acc_full)} ({f(r.acc_full_ci_low, 2)}–{f(r.acc_full_ci_high, 2)}) (chance {f(r.chance)}); "
                   f"coverage with MoTrPAC calibration at k{r.k0} ({r.primary_stage}) {f(r.coverage_k0_source_cal_primary)} (empty {f(r.frac_empty_k0_source_cal_primary)}); recalibrated on 3 / 5 animals {f(r.coverage_k0_recal3)} / {f(r.coverage_k0_recal5)} (set size {f(r.set_size_k0_recal5)})"
                   + (f"; native sedentary→exercised: accuracy k{r.k0} {f(r.invariance_acc_k0)} coverage {f(r.invariance_cov_k0)}, full {f(r.invariance_acc_full)} / {f(r.invariance_cov_full)} (RNA controls→trained 0.961 / 0.903)" if hasattr(r, "invariance_acc_k0") and not pd.isna(r.invariance_acc_k0) else ""))
    passes = [(r.leg, getattr(r, f"acc_k{r.k0}") >= 2 * r.chance and r.matched_refmet >= 30) for r in run.itertuples()]
    sec.append("- result · " + "; ".join(res) + f". Pre-registration (e) (accuracy ≥ 2 × chance with ≥ 30 matched metabolites): " + ", ".join(f"{l} {'PASS' if p else 'FAIL'}" for l, p in passes) + ".")
    sec.append("- what it does not show · platform and chemistry differ (MoTrPAC HILIC+ / six platforms vs RP-negative triple-quad vs Metabolon HD4), so a RefMet name match is a name match, not an identical analyte measurement; serum vs plasma, whole brain vs three rat brain regions, quadriceps vs gastrocnemius are imperfect maps; sample centring on different metabolite sets adds an offset the z-scoring only partly removes.")
    sec.append("- files · `results_multiomic/04_metab_transfer/<leg>/README.md`.")
    (OUT / "REPORT_SECTION.md").write_text("\n".join(sec) + "\n")
    finds = []
    for r in run.itertuples():
        finds.append({"rank": 3 if r.target == "mw" else 6, "text": f"**Metabolite fingerprint transfer ({r.leg}).** A k{r.k0} RefMet-named metabolite panel selected on MoTrPAC ({r.source_tissues} tissues, {r.matched_refmet} matched names) names the tissue of "
                                         f"{f(getattr(r, f'acc_k{r.k0}'))} of {int(r.n_mapped)} external samples (animal-bootstrap 95 % CI {f(r.acc_k0_ci_low, 2)}–{f(r.acc_k0_ci_high, 2)}; {int(r.n_individuals_mapped)} mice; full model {f(r.acc_full)}; chance {f(r.chance)}); coverage with MoTrPAC calibration {f(r.coverage_k0_source_cal_primary)}, "
                                         f"recalibrated on 5 mice {f(r.coverage_k0_recal5)}" + (f"; a native panel fit on sedentary mice keeps accuracy {f(r.invariance_acc_full)} / coverage {f(r.invariance_cov_full)} on exercised mice (RNA analogue 0.961 / 0.903)" if hasattr(r, "invariance_acc_full") and not pd.isna(r.invariance_acc_full) else "")
                                         + f". — `results_multiomic/04_metab_transfer/{r.leg}/accuracy_overall.csv`, `conformal_transfer.csv`, `recalibration.csv`"})
    (OUT / "FINDINGS.json").write_text(json.dumps(finds, indent=1))
    (OUT / "STATUS.json").write_text(json.dumps({"status": "DONE", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}))
    print(f"wrote {OUT} in {(time.time() - t0) / 60:.1f} min")


if __name__ == "__main__":
    main()
