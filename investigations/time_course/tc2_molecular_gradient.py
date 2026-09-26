#!/usr/bin/env python
"""Time course, part 2: molecular duration gradient, control vs 1w / 2w / 4w / 8w within each fusion tissue.

(a) best-single-omic AUROC ± fold sd, read from results/07_fusion/taskB_duration_summary.csv (as published);
(b) the PRE-SPECIFIED arm single:TRNSCRPT/logreg_l2, read from results/07_fusion/taskB_by_duration.csv (per-fold AUROC and
    held-out log-loss), plus optionally (--perm-arm N) a within-sex label-permutation null of that arm, recomputed with the
    same folds, prefilter and tuning as Task B (the observed arm is recomputed too and checked against the file);
(c) unsupervised separation D (see NOTES.md §0) for TRNSCRPT / PROT / METAB on the Task B animals, 1000 within-sex
    permutations;
(d) profile classification by the pre-specified rule (NOTES.md §0).

Animals: those with all three omics in the tissue (Task B's build_blocks), control + the duration.
Run from code/pipeline:  PYTHONPATH=src python investigations/time_course/tc2_molecular_gradient.py [--perm-arm 200]
Outputs: results/15_time_course/2_3_gradients/{part2_fixed_arm_per_fold.csv, part2_pc_separation.csv,
part2_fixed_arm_null.csv, part2_matrix.csv, part2_profiles.csv, part2_gradient.png}
"""
from __future__ import annotations

import os

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "8")
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import argparse
import importlib.util
import time
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from scipy import stats

from motrpac import config as C, io, models
from motrpac.splits import assert_no_group_leak, grouped_kfold

PIPE = Path(__file__).resolve().parents[2]
OUT = PIPE / "results" / "15_time_course" / "2_3_gradients"
FUS = PIPE / "results" / "07_fusion"
TISSUES = ["CORTEX", "HEART", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC"]
DURS = {"1w": 1, "2w": 2, "4w": 4, "8w": 8}
ASSAYS = ["TRNSCRPT", "PROT", "METAB"]
FIXED_ARM, FIXED_MODEL = "single:TRNSCRPT", "logreg_l2"
N_PC, N_TOP, N_PERM_D = 5, 2000, 1000


def load_script(name):
    spec = importlib.util.spec_from_file_location(name.replace(".py", "").lstrip("0123456789_"), PIPE / "scripts" / name)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def classify_profile(vals: dict, noise: float) -> tuple[str, float]:
    """vals: {weeks: metric}. Rule pre-specified in NOTES.md §0."""
    w = sorted(vals)
    v = np.array([vals[k] for k in w], float)
    rho = round(float(stats.spearmanr(w, v)[0]), 9) if np.ptp(v) > 0 else 0.0   # round: 0.8 comes back as 0.7999…
    if rho >= 0.8 and v[-1] - v[0] > noise:
        return "INCREASING", rho
    if v[-1] == v.min() and v[0] - v[-1] > noise:
        return "DECREASING", rho
    if np.ptp(v) <= noise:
        return "FLAT", rho
    return "NON-MONOTONE", rho


def pc_separation(X: np.ndarray, sex: np.ndarray, grp: np.ndarray, rng, n_perm=N_PERM_D):
    X = X.copy()
    for s in np.unique(sex):                       # remove sex (label-free, groups are sex-balanced)
        X[sex == s] -= X[sex == s].mean(axis=0)
    top = np.argsort(X.var(axis=0))[::-1][:N_TOP]
    X = X[:, top]
    sd = X.std(axis=0)
    X = X[:, sd > 0] / sd[sd > 0]
    U, S, _ = np.linalg.svd(X - X.mean(axis=0), full_matrices=False)
    Z = U[:, :N_PC] * S[:N_PC]

    def D(g):
        a, b = Z[g == 1], Z[g == 0]
        ma, mb = a.mean(0), b.mean(0)
        w = (((a - ma) ** 2).sum() + ((b - mb) ** 2).sum()) / (len(g) - 2)
        return float(np.linalg.norm(ma - mb) / np.sqrt(w))
    obs = D(grp)
    null = []
    for _ in range(n_perm):
        gp = grp.copy()
        for s in np.unique(sex):
            idx = np.flatnonzero(sex == s)
            gp[idx] = rng.permutation(gp[idx])
        null.append(D(gp))
    null = np.asarray(null)
    return {"D_obs": obs, "D_null_mean": float(null.mean()), "D_null_sd": float(null.std()),
            "e": obs - float(null.mean()), "z": (obs - null.mean()) / null.std(),
            "p_perm": float((np.sum(null >= obs) + 1) / (n_perm + 1)), "n_features_in": int(X.shape[1])}


def fixed_arm_cv(X, meta, n_splits, seed, y):
    """single:TRNSCRPT/logreg_l2 exactly as Task B: sex column appended, prefilter 2000, fit_tuned."""
    g = meta["pid"].astype(str).to_numpy()
    rows = []
    for fold, (tr, te) in enumerate(grouped_kfold(meta, "sex_group", n_splits, seed)):
        assert_no_group_leak(meta, tr, te)
        est = models.fit_tuned(FIXED_MODEL, X[tr], y[tr], g[tr], prefilter=2000, quick=False, seed=seed)
        P = est.predict_proba(X[te])
        m = models.classification_metrics(y[te], est.predict(X[te]), P, list(est.classes_))
        rows.append({"fold": fold, **m})
    return pd.DataFrame(rows)


def perm_one(X, meta, n_splits, seed, perm_seed):
    rng = np.random.default_rng(perm_seed)
    m = meta.copy()
    grp = m["group"].astype(str).to_numpy().copy()
    for sx in np.unique(m["sex"]):
        idx = np.flatnonzero((m["sex"] == sx).to_numpy())
        grp[idx] = rng.permutation(grp[idx])
    m["group"] = grp
    m["sex_group"] = m["sex"].astype(str) + "/" + m["group"].astype(str)
    pf = fixed_arm_cv(X, m, n_splits, seed, grp)
    return float(pf["auroc"].mean()), float(pf["log_loss"].mean())


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--perm-arm", type=int, default=0, help="permutations of the fixed arm per tissue × duration (0 = none)")
    ap.add_argument("--perm-jobs", type=int, default=8)
    ap.add_argument("--time-one", action="store_true", help="time the fixed-arm CV + 2 permutations on HEART 8w, then stop")
    ap.add_argument("--tissues", default=",".join(TISSUES))
    ap.add_argument("--profiles-only", action="store_true", help="recompute profiles + figure from part2_matrix.csv")
    args = ap.parse_args()
    if args.profiles_only:
        profiles_and_figure(pd.read_csv(OUT / "part2_matrix.csv"))
        return
    OUT.mkdir(parents=True, exist_ok=True)
    seed = C.SEED
    fus = load_script("07_fusion_vs_baselines.py")
    pheno = io.load_pheno()
    tissues = args.tissues.split(",")

    # (a) + (b) from the published Task B files
    summ = pd.read_csv(FUS / "taskB_duration_summary.csv")
    pf_all = pd.read_csv(FUS / "taskB_by_duration.csv")
    fa = pf_all[(pf_all["arm"] == FIXED_ARM) & (pf_all["model"] == FIXED_MODEL)].copy()
    fa.to_csv(OUT / "part2_fixed_arm_per_fold.csv", index=False)

    rng = np.random.default_rng(seed)
    sep_rows, null_rows, check_rows = [], [], []
    for t in tissues:
        blocks, meta_t, _ = fus.build_blocks([t], ASSAYS, pheno, "counts", complete=True)
        for dur in DURS:
            msk = meta_t["group"].isin(["control", dur]).to_numpy()
            sm = meta_t.loc[msk].copy()
            sm["sex_group"] = sm["sex"].astype(str) + "/" + sm["group"].astype(str)
            y01 = (sm["group"] == dur).astype(int).to_numpy()
            sex = sm["sex"].astype(str).to_numpy()
            n_c, n_t = int((y01 == 0).sum()), int((y01 == 1).sum())
            n_splits = min(4, sm["pid"].nunique() // 3)
            if args.time_one:
                if (t, dur) != ("HEART", "8w"):
                    continue
            else:
                for a in ASSAYS:
                    r = pc_separation(blocks[a].loc[msk].to_numpy(float), sex, y01, rng)
                    sep_rows.append({"tissue": t, "duration": dur, "weeks": DURS[dur], "assay": a, "n_control": n_c, "n_trained": n_t,
                                     "n_male": int((sex == "male").sum()), "n_female": int((sex == "female").sum()), **r})
                print(f"{t} {dur}: n={n_c}+{n_t}  D(TRN) e={sep_rows[-3]['e']:.3f} z={sep_rows[-3]['z']:.2f}", flush=True)
            if args.perm_arm or args.time_one:
                X = fus.with_sex(blocks["TRNSCRPT"].loc[msk], sm).to_numpy(float)
                y = sm["group"].astype(str).to_numpy()
                t0 = time.time()
                obs = fixed_arm_cv(X, sm, n_splits, seed, y)
                t_obs = time.time() - t0
                pub = fa[(fa["tissue"] == t) & (fa["duration"] == dur)]
                check_rows.append({"tissue": t, "duration": dur, "auroc_recomputed": obs["auroc"].mean(), "auroc_file": pub["auroc"].mean(),
                                   "logloss_recomputed": obs["log_loss"].mean(), "logloss_file": pub["log_loss"].mean()})
                n_p = 2 if args.time_one else args.perm_arm
                t0 = time.time()
                res = Parallel(n_jobs=1 if args.time_one else args.perm_jobs)(
                    delayed(perm_one)(X, sm, n_splits, seed, seed + 7919 * i + 13 * DURS[dur]) for i in range(n_p))
                t_perm = time.time() - t0
                for i, (au, ll) in enumerate(res):
                    null_rows.append({"tissue": t, "duration": dur, "perm": i, "auroc": au, "log_loss": ll})
                print(f"   fixed arm {t} {dur}: obs {t_obs:.1f}s (auroc {obs['auroc'].mean():.3f} vs file {pub['auroc'].mean():.3f}); "
                      f"{n_p} perms {t_perm:.1f}s", flush=True)
                if args.time_one:
                    per = t_perm / n_p
                    print(f"projected fixed-arm null, 200 perms × 28 comparisons: sequential {per * 200 * 28 / 60:.0f} min; "
                          f"with 8 jobs ≈ {per * 200 * 28 / 60 / 8:.0f} min")
                    return
    if check_rows:
        pd.DataFrame(check_rows).to_csv(OUT / "part2_fixed_arm_recompute_check.csv", index=False)
    if null_rows:
        pd.DataFrame(null_rows).to_csv(OUT / "part2_fixed_arm_null.csv", index=False)
    sep = pd.DataFrame(sep_rows)
    sep.to_csv(OUT / "part2_pc_separation.csv", index=False)

    # matrix
    fsum = fa.groupby(["tissue", "duration"]).agg(fixed_auroc=("auroc", "mean"), fixed_auroc_sd=("auroc", "std"),
                                                  fixed_logloss=("log_loss", "mean"), fixed_logloss_sd=("log_loss", "std"),
                                                  n_folds=("fold", "nunique"), n_animals=("n_animals", "first")).reset_index()
    mat = summ[["tissue", "duration", "weeks", "n_animals", "best_single_arm", "best_single_auroc", "best_single_sd"]].merge(
        fsum.drop(columns="n_animals"), on=["tissue", "duration"])
    trn = sep[sep["assay"] == "TRNSCRPT"].set_index(["tissue", "duration"])
    for a in ASSAYS:
        s = sep[sep["assay"] == a].set_index(["tissue", "duration"])
        for c in ("e", "z", "p_perm", "D_null_sd"):
            mat[f"{a}_{c}"] = [s.loc[(r.tissue, r.duration), c] if (r.tissue, r.duration) in s.index else np.nan for r in mat.itertuples()]
    mat["n_control"] = [trn.loc[(r.tissue, r.duration), "n_control"] for r in mat.itertuples()]
    mat["n_trained"] = [trn.loc[(r.tissue, r.duration), "n_trained"] for r in mat.itertuples()]
    if null_rows:
        nl = pd.DataFrame(null_rows)
        pv = []
        for r in mat.itertuples():
            q = nl[(nl["tissue"] == r.tissue) & (nl["duration"] == r.duration)]
            pv.append({"tissue": r.tissue, "duration": r.duration,
                       "fixed_auroc_p": float((np.sum(q["auroc"] >= r.fixed_auroc) + 1) / (len(q) + 1)),
                       "fixed_logloss_p": float((np.sum(q["log_loss"] <= r.fixed_logloss) + 1) / (len(q) + 1)),
                       "fixed_logloss_null_median": float(q["log_loss"].median()), "n_perm_arm": len(q)})
        mat = mat.merge(pd.DataFrame(pv), on=["tissue", "duration"])
    mat = mat.sort_values(["tissue", "weeks"])
    mat.to_csv(OUT / "part2_matrix.csv", index=False)
    profiles_and_figure(mat)


def profiles_and_figure(mat):

    # profiles
    prow = []
    for t, d in mat.groupby("tissue"):
        d = d.set_index("weeks")
        out = {"tissue": t}
        for a in ASSAYS:
            noise = float(np.sqrt(d.loc[1, f"{a}_D_null_sd"] ** 2 + d.loc[8, f"{a}_D_null_sd"] ** 2))
            cls, rho = classify_profile(d[f"{a}_e"].to_dict(), noise)
            out.update({f"{a}_profile": cls, f"{a}_rho": rho, f"{a}_noise": noise})
        noise = float(np.mean([d.loc[1, "fixed_logloss_sd"], d.loc[8, "fixed_logloss_sd"]]))
        cls, rho = classify_profile((-d["fixed_logloss"]).to_dict(), noise)   # higher = better separation
        out.update({"fixed_neglogloss_profile": cls, "fixed_neglogloss_rho": rho, "fixed_logloss_noise": noise})
        noise = float(np.mean([d.loc[1, "fixed_auroc_sd"], d.loc[8, "fixed_auroc_sd"]]))
        cls, rho = classify_profile(d["fixed_auroc"].to_dict(), noise)
        out.update({"fixed_auroc_profile": cls, "fixed_auroc_rho": rho})
        prow.append(out)
    prof = pd.DataFrame(prow)
    prof.to_csv(OUT / "part2_profiles.csv", index=False)
    print(mat.round(3).to_string(index=False))
    print(prof.round(2).to_string(index=False))

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axs = plt.subplots(1, 3, figsize=(13, 4))
    for t, d in mat.groupby("tissue"):
        d = d.sort_values("weeks")
        axs[0].errorbar(d["weeks"], d["TRNSCRPT_e"], yerr=d["TRNSCRPT_D_null_sd"], marker="o", capsize=2, label=t)
        axs[1].errorbar(d["weeks"], d["fixed_logloss"], yerr=d["fixed_logloss_sd"], marker="o", capsize=2, label=t)
        axs[2].errorbar(d["weeks"], d["fixed_auroc"], yerr=d["fixed_auroc_sd"], marker="o", capsize=2, label=t)
    for ax, yl in zip(axs, ["TRNSCRPT PC separation, D − null mean (±null sd)", "held-out log-loss, TRNSCRPT/logreg_l2 (±fold sd)",
                            "AUROC, TRNSCRPT/logreg_l2 (±fold sd)"]):
        ax.set_xscale("log", base=2)
        ax.set_xticks([1, 2, 4, 8])
        ax.set_xticklabels(["1w*", "2w*", "4w*", "8w"])
        ax.set_ylabel(yl, fontsize=8)
        ax.grid(alpha=0.3)
    axs[0].axhline(0, color="grey", lw=1, ls=":")
    axs[0].legend(fontsize=7, frameon=False)
    fig.suptitle("Control vs Nw within tissue (Task B animals, n ≈ 17–20). *1w/2w/4w: controls from a different cohort/season", fontsize=9)
    fig.tight_layout()
    fig.savefig(OUT / "part2_gradient.png", dpi=150)
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
