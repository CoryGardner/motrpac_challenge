#!/usr/bin/env python
"""Time course, part 3: covariate-only gradient, control vs 1w / 2w / 4w / 8w within each fusion tissue, using ONLY RNA
processing covariates (no omics).

Feature sets reuse scripts/07_batch_check.py's definitions (categorize / covariate_frame on data/raw/meta/TRNSCRPT.csv,
plate row / column derived from Lib_barcode_well exactly as there):
  trnscrpt_depth, trnscrpt_qc, trnscrpt_library (as defined: includes plate row / column),
  trnscrpt_library_noposition (same minus Lib_barcode_well_row / _col, the plate-position axis),
  rna_depth_qc (= depth ∪ qc; the PRE-SPECIFIED decisive set, NOTES.md §0).
Models: logistic (C = 0.1, with within-sex permutation null) and RF (500 trees, no null), as in 07_batch_check.
Animals: Task B animals (all three omics, same as part 2), control + the duration; folds grouped_kfold(meta, 'sex_group',
min(4, n // 3), SEED). Reproduction check: 8w on 07_batch_check's own animal set vs results/07_fusion/batch_covariate_auroc.csv.

Run from code/pipeline:  PYTHONPATH=src python investigations/time_course/tc3_covariate_gradient.py [--n-perm 200]
Outputs: results/15_time_course/2_3_gradients/{part3_covariate_auroc.csv, part3_repro_check.csv, part3_profiles.csv,
part3_vs_part2.csv, part3_gradient.png}
"""
from __future__ import annotations

import os

os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import argparse
import warnings

from joblib import Parallel, delayed

import numpy as np
import pandas as pd
from sklearn.compose import ColumnTransformer
from sklearn.ensemble import RandomForestClassifier
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import roc_auc_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import OneHotEncoder, StandardScaler

from tfp import config as C, io
from tfp.splits import assert_no_group_leak, grouped_kfold

import sys
sys.path.insert(0, os.path.dirname(__file__))
from tc2_molecular_gradient import DURS, OUT, TISSUES, classify_profile, load_script  # noqa: E402

NULL_SETS = {"trnscrpt_depth", "trnscrpt_qc", "rna_depth_qc"}   # logistic permutation null only for these (cost)
SETS = ["trnscrpt_depth", "trnscrpt_qc", "trnscrpt_library", "trnscrpt_library_noposition", "rna_depth_qc"]
DECISIVE = ("rna_depth_qc", "logreg")


def cv_auroc(X, num, cat, y, meta, n_splits, seed, model, n_perm, perm_seed):
    warnings.simplefilter("ignore")
    rng = np.random.default_rng(perm_seed)
    pre = ColumnTransformer([("num", Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler())]), num),
                             ("cat", OneHotEncoder(handle_unknown="ignore"), cat)])
    clf = (LogisticRegression(C=0.1, max_iter=3000) if model == "logreg"
           else RandomForestClassifier(n_estimators=500, class_weight="balanced", random_state=seed, n_jobs=1))
    pipe = Pipeline([("pre", pre), ("clf", clf)])
    folds = list(grouped_kfold(meta, "sex_group", n_splits, seed))

    def run(yy):
        a = []
        for tr, te in folds:
            assert_no_group_leak(meta, tr, te)
            if len(np.unique(yy[te])) < 2 or len(np.unique(yy[tr])) < 2:
                continue
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pipe.fit(X.iloc[tr], yy[tr])
                a.append(roc_auc_score(yy[te], pipe.predict_proba(X.iloc[te])[:, 1]))
        return a
    a = run(y)
    out = {"auroc_mean": float(np.mean(a)), "auroc_sd": float(np.std(a)), "auroc_sd_ddof1": float(np.std(a, ddof=1)), "n_folds": len(a)}
    if n_perm:
        sex = meta["sex"].to_numpy()
        null = []
        for _ in range(n_perm):
            yp = y.copy()
            for sx in np.unique(sex):
                idx = np.flatnonzero(sex == sx)
                yp[idx] = rng.permutation(yp[idx])
            null.append(float(np.mean(run(yp))))
        out.update({"null_mean": float(np.mean(null)), "null_p95_auroc": float(np.quantile(null, 0.95)),
                    "p_perm": float((np.sum(np.asarray(null) >= out["auroc_mean"]) + 1) / (n_perm + 1)), "n_perm": n_perm})
    return out


def frames_for(bc, d, tv):
    trn_cols = [c for c in tv.columns if bc.categorize(c)]
    base = {"trnscrpt_library": [c for c in trn_cols if bc.categorize(c) == "library"],
            "trnscrpt_depth": [c for c in trn_cols if bc.categorize(c) == "depth"],
            "trnscrpt_qc": [c for c in trn_cols if bc.categorize(c) == "qc"]}
    base["trnscrpt_library_noposition"] = [c for c in base["trnscrpt_library"] if c not in ("Lib_barcode_well_row", "Lib_barcode_well_col")]
    base["rna_depth_qc"] = base["trnscrpt_depth"] + base["trnscrpt_qc"]
    return {fs: bc.covariate_frame(tv, cols) for fs, cols in base.items()}


def vial_table(tm, d):
    tv = tm.set_index("viallabel").reindex(d["viallabel"]).set_axis(d.index).copy()
    w = tv["Lib_barcode_well"].astype("string").str.upper().str.extract(r"^([A-H])(\d{1,2})$")
    tv["Lib_barcode_well_row"] = w[0].map({c: i + 1 for i, c in enumerate("ABCDEFGH")}).astype("string")
    tv["Lib_barcode_well_col"] = w[1]
    return tv


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--n-perm", type=int, default=200)
    ap.add_argument("--jobs", type=int, default=8)
    args = ap.parse_args()
    seed = C.SEED
    warnings.simplefilter("ignore")
    bc = load_script("07_batch_check.py")
    fus = load_script("07_fusion_vs_baselines.py")
    pheno = io.load_pheno()
    ph = pd.read_csv(C.RAW_DIR / "pheno.csv", dtype=str, low_memory=False)
    tm = pd.read_csv(C.META_DIR / "TRNSCRPT.csv", dtype=str, low_memory=False).drop_duplicates("viallabel")

    tasks = []
    for t in TISSUES:
        _, meta_t, _ = fus.build_blocks([t], ["TRNSCRPT", "PROT", "METAB"], pheno, "counts", complete=True)
        taskb_pids = set(meta_t["pid"].astype(str))
        for dur in DURS:
            for animals in (["taskB"] + (["batchcheck"] if dur == "8w" else [])):
                d = ph[(ph["tissue"] == t) & ph["group"].isin(["control", dur]) & ph["viallabel"].isin(tm["viallabel"])].drop_duplicates("pid").copy()
                if animals == "taskB":
                    d = d[d["pid"].isin(taskb_pids)]
                d = d.set_index("pid")
                y = (d["group"] == dur).astype(int).to_numpy()
                meta = pd.DataFrame({"pid": d.index, "sex": d["sex"].to_numpy(), "group": d["group"].to_numpy()}, index=d.index)
                meta["sex_group"] = meta["sex"] + "/" + meta["group"]
                n_splits = min(4, len(d) // 3)
                frames = frames_for(bc, d, vial_table(tm, d))
                for fs in SETS:
                    X, num, cat = frames[fs]
                    if X.shape[1] == 0:
                        continue
                    for model in ("logreg", "rf"):
                        if animals == "batchcheck" and fs in ("trnscrpt_library_noposition", "rna_depth_qc"):
                            continue
                        npm = args.n_perm if (model == "logreg" and animals == "taskB" and fs in NULL_SETS) else 0
                        row = {"tissue": t, "duration": dur, "weeks": DURS[dur], "animals": animals, "feature_set": fs, "model": model,
                               "n_animals": len(d), "n_control": int((y == 0).sum()), "n_trained": int(y.sum()), "n_features": X.shape[1]}
                        tasks.append((row, (X, num, cat, y, meta, n_splits, seed, model, npm, seed + len(tasks))))
    print(f"{len(tasks)} covariate classifiers, n_perm={args.n_perm}", flush=True)
    outs = Parallel(n_jobs=args.jobs, verbose=5)(delayed(cv_auroc)(*a) for _, a in tasks)
    rows, repro = [], []
    for (row, _), r in zip(tasks, outs):
        (repro if row["animals"] == "batchcheck" else rows).append({**row, **r})
    res = pd.DataFrame(rows)
    res.to_csv(OUT / "part3_covariate_auroc.csv", index=False)
    rp = pd.DataFrame(repro)
    old = pd.read_csv(C.RESULTS_DIR / "07_fusion" / "batch_covariate_auroc.csv")
    rp = rp.merge(old[["tissue", "feature_set", "model", "auroc_mean", "n_animals"]].rename(columns={"auroc_mean": "auroc_published", "n_animals": "n_published"}),
                  on=["tissue", "feature_set", "model"], how="left")
    rp["abs_diff"] = (rp["auroc_mean"] - rp["auroc_published"]).abs()
    rp.to_csv(OUT / "part3_repro_check.csv", index=False)
    print(f"reproduction of batch_covariate_auroc.csv (8w): max |diff| = {rp['abs_diff'].max():.4f} over {len(rp)} rows")

    # profiles: every set × model, decisive flagged
    prow = []
    for (t, fs, model), d in res.groupby(["tissue", "feature_set", "model"]):
        d = d.set_index("weeks")
        noise = float(np.mean([d.loc[1, "auroc_sd_ddof1"], d.loc[8, "auroc_sd_ddof1"]]))
        cls, rho = classify_profile(d["auroc_mean"].to_dict(), noise)
        prow.append({"tissue": t, "feature_set": fs, "model": model, "decisive": (fs, model) == DECISIVE, "profile": cls, "rho": rho, "noise": noise,
                     **{f"auroc_{w}w": d.loc[w, "auroc_mean"] for w in (1, 2, 4, 8)},
                     **{f"p_{w}w": (d.loc[w, "p_perm"] if "p_perm" in d else np.nan) for w in (1, 2, 4, 8)}})
    prof = pd.DataFrame(prow)
    prof.to_csv(OUT / "part3_profiles.csv", index=False)
    print(prof[prof["model"] == "logreg"].round(2).to_string(index=False))

    # side by side with part 2
    p2 = OUT / "part2_profiles.csv"
    if p2.exists():
        m2 = pd.read_csv(p2)
        dec = prof[prof["decisive"]][["tissue", "profile", "rho"] + [f"auroc_{w}w" for w in (1, 2, 4, 8)]].rename(
            columns={"profile": "covariate_profile", "rho": "covariate_rho"})
        both = m2[["tissue", "TRNSCRPT_profile", "TRNSCRPT_rho", "fixed_neglogloss_profile"]].merge(dec, on="tissue")
        both.to_csv(OUT / "part3_vs_part2.csv", index=False)
        print(both.round(2).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(13, 4), sharey=True)
    for ax, fs in zip(axs, ["rna_depth_qc", "trnscrpt_depth", "trnscrpt_library"]):
        for t, d in res[(res["feature_set"] == fs) & (res["model"] == "logreg")].groupby("tissue"):
            d = d.sort_values("weeks")
            ax.errorbar(d["weeks"], d["auroc_mean"], yerr=d["auroc_sd_ddof1"], marker="o", capsize=2, label=t)
        ax.set_xscale("log", base=2)
        ax.set_xticks([1, 2, 4, 8])
        ax.set_xticklabels(["1w*", "2w*", "4w*", "8w"])
        ax.axhline(0.5, color="grey", ls=":", lw=1)
        ax.set_title(f"{fs} (logistic)", fontsize=9)
        ax.grid(alpha=0.3)
    axs[0].set_ylabel("covariate-only AUROC, control vs Nw (±fold sd)")
    axs[0].legend(fontsize=7, frameon=False)
    fig.suptitle("RNA processing covariates only (Task B animals). *1w/2w/4w: controls from a different cohort/season", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "part3_gradient.png", dpi=150)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
