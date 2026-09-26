#!/usr/bin/env python
"""Time course, parts 5 and 6: tissue fingerprint (19-way TRNSCRPT tissue classifier) under training duration.

PART 5  Train + calibrate on sedentary controls only (the `make shift` controls->trained split: controls from
        splits.train_controls_test_trained, fit/calibration by fit_calibration_split(meta, controls, 0.3, SEED),
        pooled-vial LAC calibration, alpha 0.10, logreg_l2, prefilter 5000, k = 20 panel and full model).
        One fitted model per arm, evaluated SEPARATELY on the 1w, 2w, 4w and 8w animals.
        Control arm ("matched"): 5 grouped folds over the 10 control animals; per fold 2 control animals are held
        out, the other 8 are split 70/30 by animal into fit/calibration (6/2), and the same fitted model is scored
        on the 2 held-out controls and on every trained duration. Held-out controls are pooled over the folds
        (each control animal tested once); trained durations are averaged over folds.
PART 6  Hold out 8w; fit + calibrate on (a) control+1w, (b) control+1w+2w+4w (= make shift's holdout_group_8w),
        (c) control+4w (size-matched to (a), temporally adjacent), (d) control only (part-5 primary);
        same 10 8w test animals, same settings.

Intervals: percentile cluster bootstrap over test animals (B = 2000), conditional on the fitted model and the
calibration set (they do NOT include calibration-set randomness; the matched-fold spread shows that part).

Outputs: results/15_time_course/5_6_fingerprint/{part5_by_duration.csv, part5_matched_folds.csv,
part5_matched_summary.csv, part6_extrapolation.csv, calibration_info.csv, per_vial.csv.gz, fingerprint_time_course.png}
Run from code/pipeline:  PYTHONPATH=src python investigations/time_course/tc56_fingerprint_shift.py
"""
from __future__ import annotations

import os

os.environ.setdefault("LOKY_MAX_CPU_COUNT", "8")   # GridSearchCV uses n_jobs=-1; cap it (shared machine)
os.environ.setdefault("OMP_NUM_THREADS", "2")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "2")

import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "scripts"))
from motrpac import cli, config as C, conformal as cp, io, models  # noqa: E402
from motrpac.splits import assert_no_group_leak, fit_calibration_split, grouped_kfold, \
    train_controls_test_trained  # noqa: E402

ALPHA, K, MODEL, PREFILTER, SEED, B = 0.10, 20, "logreg_l2", 5000, C.SEED, 2000
DURS = ["1w", "2w", "4w", "8w"]
OUT = C.RESULTS_DIR / "15_time_course" / "5_6_fingerprint"
# consortium outlier list (read-only); rows with reason "BAT contamination" are female 1w/2w vena cava
OUTLIERS_TXT = Path(__file__).resolve().parents[4] / "data/analysis/rat-training-06/c1.0/resources/OUTLIERS.txt"


def bat_contaminated_pids() -> set[str]:
    o = pd.read_csv(OUTLIERS_TXT, sep="\t", dtype=str)
    o = o[o["reason"].fillna("").str.contains("BAT contamination") & (o["tissue"] == "VENACV")]
    return set(o["pid"].astype(str))


def drop_contam(d: pd.DataFrame, pids: set[str]) -> pd.DataFrame:
    return d[~((d["y_true"] == "VENACV") & d["pid"].isin(pids))]


def boot_diff(d1: pd.DataFrame, d2: pd.DataFrame, rng, paired: bool) -> dict:
    """metric(d1) - metric(d2); paired = same animals in both (resample animals jointly)."""
    point = {k: metrics(d1)[k] - metrics(d2)[k] for k in ("accuracy", "coverage", "empty_rate", "mean_set_size")}
    b1, b2 = {p: x for p, x in d1.groupby("pid")}, {p: x for p, x in d2.groupby("pid")}
    p1, p2 = np.array(sorted(b1)), np.array(sorted(b2))
    boots = []
    for _ in range(B):
        s1 = rng.choice(p1, size=len(p1), replace=True)
        s2 = s1 if paired else rng.choice(p2, size=len(p2), replace=True)
        m1 = metrics(pd.concat([b1[p] for p in s1], ignore_index=True))
        m2 = metrics(pd.concat([b2[p] for p in s2], ignore_index=True))
        boots.append({k: m1[k] - m2[k] for k in point})
    boots = pd.DataFrame(boots)
    row = {"n_animals_1": len(p1), "n_animals_2": len(p2), "paired": paired}
    for k, v in point.items():
        row[f"d_{k}"] = v
        row[f"d_{k}_lo"], row[f"d_{k}_hi"] = float(boots[k].quantile(0.025)), float(boots[k].quantile(0.975))
    return row


def fit_and_calibrate(X, y, g, meta, fit_idx, cal_idx, k):
    assert_no_group_leak(meta, fit_idx, cal_idx)
    est = models.fit_tuned(MODEL, X[fit_idx], y[fit_idx], g[fit_idx], k=k, prefilter=PREFILTER, quick=False, seed=SEED)
    cls = list(est.classes_)
    scores, n_drop = cp.calibration_scores(est.predict_proba(X[cal_idx]), y[cal_idx], cls, "lac")
    qhat = cp.conformal_quantile(scores, ALPHA)
    n = len(scores)
    rank = int(np.ceil((n + 1) * (1 - ALPHA) - 1e-9))
    info = {"n_fit_animals": int(meta.iloc[fit_idx]["pid"].nunique()), "n_fit_vials": len(fit_idx),
            "n_cal_animals": int(meta.iloc[cal_idx]["pid"].nunique()), "n_cal_vials": n, "n_cal_dropped": n_drop,
            "conformal_rank": rank, "qhat": qhat, "qhat_is_inf": bool(np.isinf(qhat)),
            "cal_score_max": float(scores.max()), "qhat_is_max_score": bool(np.isfinite(qhat) and qhat >= scores.max()),
            # what an animal-level calibration (one score per animal) would give with the same animals
            "animal_level_rank": int(np.ceil((meta.iloc[cal_idx]['pid'].nunique() + 1) * (1 - ALPHA) - 1e-9)),
            "animal_level_qhat_inf": bool(np.ceil((meta.iloc[cal_idx]['pid'].nunique() + 1) * (1 - ALPHA) - 1e-9)
                                          > meta.iloc[cal_idx]['pid'].nunique())}
    return est, cls, qhat, info


def score_vials(est, cls, qhat, X, y, meta, idx):
    p = est.predict_proba(X[idx])
    sets = cp.predict_sets(p, qhat, "lac")
    yi, seen = cp.class_indices(y[idx], cls)
    covered = np.zeros(len(idx), dtype=bool)
    covered[seen] = sets[np.flatnonzero(seen), yi[seen]]
    m = meta.iloc[idx]
    return pd.DataFrame({"pid": m["pid"].astype(str).to_numpy(), "group": m["group"].astype(str).to_numpy(),
                         "sex": m["sex"].astype(str).to_numpy(), "y_true": y[idx], "y_pred": est.predict(X[idx]),
                         "covered": covered, "seen": seen, "set_size": sets.sum(axis=1)})


def metrics(d: pd.DataFrame) -> dict:
    s = d[d["seen"]]
    return {"accuracy": float((d["y_true"] == d["y_pred"]).mean()),
            "bal_acc": float(models.balanced_accuracy_score(s["y_true"], s["y_pred"])) if len(s) else np.nan,
            "coverage": float(s["covered"].mean()) if len(s) else np.nan,
            "empty_rate": float((d["set_size"] == 0).mean()), "mean_set_size": float(d["set_size"].mean()),
            "singleton_rate": float((d["set_size"] == 1).mean())}


def summarize(d: pd.DataFrame, rng) -> dict:
    point = metrics(d)
    pids = d["pid"].unique()
    by = {p: x for p, x in d.groupby("pid")}
    boots = []
    for _ in range(B):
        draw = rng.choice(pids, size=len(pids), replace=True)
        boots.append(metrics(pd.concat([by[p] for p in draw], ignore_index=True)))
    boots = pd.DataFrame(boots)
    row = {"n_test_animals": len(pids), "n_test_vials": len(d), "n_unseen_vials": int((~d["seen"]).sum())}
    for key, v in point.items():
        row[key] = v
        row[f"{key}_lo"], row[f"{key}_hi"] = (float(np.nanquantile(boots[key], 0.025)), float(np.nanquantile(boots[key], 0.975)))
    # Clopper-Pearson on vials (ignores within-animal clustering; shown for reference only)
    s = d[d["seen"]]
    from scipy.stats import beta
    x, n = int(s["covered"].sum()), len(s)
    row["coverage_cp_lo"] = float(beta.ppf(0.025, x, n - x + 1)) if x > 0 else 0.0
    row["coverage_cp_hi"] = float(beta.ppf(0.975, x + 1, n - x)) if x < n else 1.0
    # animals whose every vial is covered / fraction of animals with coverage < 1-alpha
    per_animal = s.groupby("pid")["covered"].mean()
    row["animals_cov_below_1ma"] = int((per_animal < 1 - ALPHA).sum())
    return row


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    pheno = io.load_pheno()
    om = io.stack_tissues("TRNSCRPT", source=cli.resolve_source("TRNSCRPT", "auto"), join="inner", pheno=pheno)
    X, y, g, meta = om.X.to_numpy(dtype=float), om.meta["tissue"].astype(str).to_numpy(), om.groups(), om.meta
    grp = meta["group"].astype(str).to_numpy()
    print(f"loaded {X.shape}, {meta['tissue'].nunique()} tissues, {meta['pid'].nunique()} animals "
          f"({time.time() - t0:.0f}s)", flush=True)
    print(meta.groupby("group")["pid"].nunique().to_dict(), flush=True)
    rng = np.random.default_rng(SEED)
    contam = bat_contaminated_pids()
    print(f"consortium 'BAT contamination' VENACV animals: {len(contam)}", flush=True)
    arms = (("full", None), (f"panel_k{K}", K))
    cal_rows, vial_frames = [], []

    # ---------------- PART 5 primary: controls -> each duration --------------------------------
    [(name, (tr, te))] = list(train_controls_test_trained(meta))
    fit_idx, cal_idx = fit_calibration_split(meta, tr, 0.3, SEED)      # exactly as eval_split does it
    p5, diffs = [], []
    for arm, k in arms:
        est, cls, qhat, info = fit_and_calibrate(X, y, g, meta, fit_idx, cal_idx, k)
        cal_rows.append({"design": "p5_primary_control_only", "arm": arm, **info})
        dv = score_vials(est, cls, qhat, X, y, meta, te)
        vial_frames.append(dv.assign(design="p5_primary_control_only", arm=arm))
        pooled = metrics(dv)
        print(f"  [check vs results/08_shift train_control_test_trained] {arm}: acc={pooled['accuracy']:.6f} "
              f"cov={pooled['coverage']:.6f} size={pooled['mean_set_size']:.6f} qhat={qhat:.6f}", flush=True)
        for dur in DURS + ["pooled_trained"]:
            d = dv if dur == "pooled_trained" else dv[dv["group"] == dur]
            p5.append({"design": "control_only_fit7_cal3", "arm": arm, "test_group": dur, **summarize(d, rng),
                       "qhat": qhat, "n_cal_vials": info["n_cal_vials"], "n_cal_animals": info["n_cal_animals"]})
            p5.append({"design": "control_only_fit7_cal3_excl_BATcontam_VENACV", "arm": arm, "test_group": dur,
                       **summarize(drop_contam(d, contam), rng), "qhat": qhat, "n_cal_vials": info["n_cal_vials"],
                       "n_cal_animals": info["n_cal_animals"]})
        for excl in (False, True):
            d8 = dv[dv["group"] == "8w"]
            for dur in ("1w", "2w", "4w"):
                d1 = dv[dv["group"] == dur]
                if excl:
                    d1 = drop_contam(d1, contam)
                diffs.append({"part": 5, "arm": arm, "comparison": f"{dur} - 8w (control-only model)",
                              "excl_BATcontam_VENACV": excl, **boot_diff(d1, d8, rng, paired=False)})
    pd.DataFrame(p5).to_csv(OUT / "part5_by_duration.csv", index=False)
    print(f"part 5 primary done ({time.time() - t0:.0f}s)", flush=True)

    # ---------------- PART 5 control arm: matched folds over controls --------------------------
    ctrl = np.flatnonzero(grp == "control")
    mc = meta.iloc[ctrl].reset_index(drop=True)
    mf_rows, held_ctrl = [], {a: [] for a, _ in arms}
    trained_by_fold = {a: [] for a, _ in arms}
    for fold, (ftr, fte) in enumerate(grouped_kfold(mc, "sex", n_splits=5, seed=SEED)):
        tr_c, te_c = ctrl[ftr], ctrl[fte]
        assert_no_group_leak(meta, tr_c, te_c)
        f_idx, c_idx = fit_calibration_split(meta, tr_c, 0.3, SEED + fold)
        for arm, k in arms:
            est, cls, qhat, info = fit_and_calibrate(X, y, g, meta, f_idx, c_idx, k)
            cal_rows.append({"design": f"p5_matched_fold{fold}", "arm": arm, **info})
            test_idx = np.concatenate([te_c, np.flatnonzero(grp != "control")])
            dv = score_vials(est, cls, qhat, X, y, meta, test_idx).assign(fold=fold)
            dv["test_group"] = np.where(np.isin(np.arange(len(test_idx)), np.arange(len(te_c))), "control_heldout", dv["group"])
            held_ctrl[arm].append(dv[dv["test_group"] == "control_heldout"])
            trained_by_fold[arm].append(dv[dv["test_group"] != "control_heldout"])
            vial_frames.append(dv.assign(design=f"p5_matched_fold{fold}", arm=arm))
            for tg, d in dv.groupby("test_group"):
                mf_rows.append({"fold": fold, "arm": arm, "test_group": tg, "n_test_animals": d["pid"].nunique(),
                                "n_test_vials": len(d), "n_fit_animals": info["n_fit_animals"],
                                "n_cal_animals": info["n_cal_animals"], "n_cal_vials": info["n_cal_vials"],
                                "qhat": qhat, **metrics(d)})
        print(f"  matched fold {fold} done ({time.time() - t0:.0f}s)", flush=True)
    mf = pd.DataFrame(mf_rows)
    mf.to_csv(OUT / "part5_matched_folds.csv", index=False)
    ms = []
    for arm, _ in arms:
        hc = pd.concat(held_ctrl[arm], ignore_index=True)       # every control animal tested once
        ms.append({"arm": arm, "test_group": "control_heldout (pooled over folds)", "n_folds": 5, **summarize(hc, rng)})
        sub = mf[(mf["arm"] == arm) & (mf["test_group"] != "control_heldout")]
        for tg, s in sub.groupby("test_group"):
            r = {"arm": arm, "test_group": f"{tg} (mean over folds)", "n_folds": len(s),
                 "n_test_animals": int(s["n_test_animals"].iloc[0]), "n_test_vials": int(s["n_test_vials"].iloc[0])}
            for key in ("accuracy", "bal_acc", "coverage", "empty_rate", "mean_set_size", "singleton_rate"):
                r[key], r[f"{key}_fold_sd"] = float(s[key].mean()), float(s[key].std(ddof=1))
                r[f"{key}_fold_min"], r[f"{key}_fold_max"] = float(s[key].min()), float(s[key].max())
            ms.append(r)
        # per-fold paired difference: coverage(8w) - coverage(held-out controls), same model and calibration
        s8 = mf[(mf["arm"] == arm)].pivot(index="fold", columns="test_group", values="coverage")
        for tg in DURS:
            dlt = s8[tg] - s8["control_heldout"]
            ms.append({"arm": arm, "test_group": f"coverage[{tg}] - coverage[control_heldout], per fold",
                       "n_folds": len(dlt), "coverage": float(dlt.mean()), "coverage_fold_sd": float(dlt.std(ddof=1)),
                       "coverage_fold_min": float(dlt.min()), "coverage_fold_max": float(dlt.max())})
    pd.DataFrame(ms).to_csv(OUT / "part5_matched_summary.csv", index=False)
    print(f"part 5 matched done ({time.time() - t0:.0f}s)", flush=True)

    # ---------------- PART 6: temporal extrapolation to 8w -------------------------------------
    te8 = np.flatnonzero(grp == "8w")
    designs = {"a_control+1w": ["control", "1w"], "b_control+1w+2w+4w": ["control", "1w", "2w", "4w"],
               "c_control+4w": ["control", "4w"], "d_control_only": ["control"]}
    p6 = []
    for dname, groups in designs.items():
        trg = np.flatnonzero(np.isin(grp, groups))
        assert_no_group_leak(meta, trg, te8)
        f_idx, c_idx = fit_calibration_split(meta, trg, 0.3, SEED)
        for arm, k in arms:
            est, cls, qhat, info = fit_and_calibrate(X, y, g, meta, f_idx, c_idx, k)
            cal_rows.append({"design": f"p6_{dname}", "arm": arm, **info})
            dv = score_vials(est, cls, qhat, X, y, meta, te8)
            vial_frames.append(dv.assign(design=f"p6_{dname}", arm=arm))
            cal_groups = meta.iloc[c_idx].groupby("group")["pid"].nunique().to_dict()
            p6.append({"design": dname, "arm": arm, "train_groups": "+".join(groups),
                       "n_train_animals": int(meta.iloc[trg]["pid"].nunique()), "n_fit_animals": info["n_fit_animals"],
                       "n_cal_animals": info["n_cal_animals"], "n_cal_vials": info["n_cal_vials"],
                       "cal_animals_by_group": str(cal_groups), "qhat": qhat, **summarize(dv, rng)})
            print(f"  {dname} {arm}: acc={p6[-1]['accuracy']:.4f} cov={p6[-1]['coverage']:.4f} "
                  f"size={p6[-1]['mean_set_size']:.3f} ({time.time() - t0:.0f}s)", flush=True)
    pd.DataFrame(p6).to_csv(OUT / "part6_extrapolation.csv", index=False)
    vp6 = pd.concat([f for f in vial_frames if f["design"].iloc[0].startswith("p6_")], ignore_index=True)
    for arm, _ in arms:
        ref = vp6[(vp6["design"] == "p6_b_control+1w+2w+4w") & (vp6["arm"] == arm)]
        for dname in ("a_control+1w", "c_control+4w", "d_control_only"):
            cur = vp6[(vp6["design"] == f"p6_{dname}") & (vp6["arm"] == arm)]
            diffs.append({"part": 6, "arm": arm, "comparison": f"{dname} - b_control+1w+2w+4w (same 8w animals)",
                          "excl_BATcontam_VENACV": False, **boot_diff(cur, ref, rng, paired=True)})
    pd.DataFrame(diffs).to_csv(OUT / "differences.csv", index=False)
    pd.DataFrame(cal_rows).to_csv(OUT / "calibration_info.csv", index=False)
    pd.concat(vial_frames, ignore_index=True).to_csv(OUT / "per_vial.csv.gz", index=False)
    print(f"done ({time.time() - t0:.0f}s) -> {OUT}", flush=True)


if __name__ == "__main__":
    main()
