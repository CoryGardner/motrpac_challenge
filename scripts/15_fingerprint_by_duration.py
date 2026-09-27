#!/usr/bin/env python
"""Phase 15b — is the tissue fingerprint invariant to training duration?

The task is TISSUE classification: the 19-tissue TRNSCRPT fingerprint with the phase 04–08 settings. It is not
classification of training state. One model is fitted and calibrated on sedentary controls only and scored
separately on the 1w, 2w, 4w and 8w animals; a matched-fold control arm gives the in-distribution reference; and
the calibration set is moved back in time to ask whether temporal distance costs coverage on held-out 8w animals.

Settings (they match `make shift` for TRNSCRPT). Counts source, log2 CPM with the stacked gene filter, all study
vials of the 19 tissues (50 animals, 10 per group), `logreg_l2` tuned by inner GroupKFold, variance prefilter 5000
inside the fit, two arms: the k = 20 round-robin panel and the full model; LAC scores, α = 0.10, calibration on the
pooled vials of the calibration animals, the textbook conformal quantile (post-2026-09-25). `--quick` only reduces
the number of fit/calibration splits and bootstrap draws; the fits are never reduced, so the reproduction check
below holds in every mode. Known deviation, inherited: the transcript gene filter runs before the split.

PART 5, primary design `control_only_fit7_cal3`. The phase-08 controls→trained split (`train_controls_test_trained`,
then `fit_calibration_split(meta, controls, 0.3, SEED)`): 7 fit and 3 calibration control animals (54 calibration
vials), one fitted model per arm, evaluated per duration (10 animals, ~180 vials each). Rows
`control_only_fit7_cal3_excl_BATcontam_VENACV` drop the vena-cava vials of the 10 female 1w/2w animals the consortium
flags "BAT contamination" (`io.load_outliers()`, 30 vials; they are still in the count matrices).
Reference arm (`part5_matched_*`): 5 grouped folds over the 10 control animals; per fold 2 controls are held out,
the other 8 are split 6 fit / 2 calibration by animal, and the same model is scored on the held-out controls and on
every trained duration. Held-out controls are pooled over folds (each control tested once); trained durations are
averaged over folds, and coverage[duration] − coverage[held-out controls] is taken within fold.

PART 6, temporal extrapolation to held-out 8w. Fit + calibrate (70/30 by animal) on (a) control+1w, (b)
control+1w+2w+4w (= `make shift`'s holdout_group_8w), (c) control+4w (size-matched to (a), adjacent in time),
(d) control only (the part-5 model); the same 10 8w animals are the test set. `part6_seeds*`: the four designs
repeated over `--seeds` random fit/calibration splits (`fit_calibration_split` seeds SEED+1..SEED+N, model seed
fixed), so the (a) − (b) comparison includes calibration-set randomness.

Intervals. 95 % percentile cluster bootstrap over test animals (`--bootstraps` draws), conditional on the fitted
model and the calibration set. `differences.csv`: duration − 8w for the control-only model (unpaired animal
bootstrap) and design − (b) on the same 8w animals (paired). Clopper–Pearson intervals on vials are given for
reference only; they ignore the clustering of vials within animals.

Reproduction check. Pooled over the 40 trained animals, the primary design must reproduce
results/08_shift/TRNSCRPT/shift_table.csv (`train_control_test_trained`: accuracy, coverage, set size, empty rate),
and per duration the numbers recomputed from results/31_site_regen/08_shift_k20/scores_target_vials.csv (accuracy =
mean(y_pred == tissue), coverage = mean(covered_marginal), empty = mean(size_marginal == 0)); design (b) must
reproduce `holdout_group_8w`. The script raises if any accuracy or coverage differs by more than 1e-6. qhat is
compared but only reported: a 4th-decimal wobble from BLAS threading was seen before while the sets were identical,
so the per-vial sets (y_pred, covered, set size) are compared vial by vial and the counts reported. The check is
skipped, with a message, when the reference files are absent (results/ first, then results_frozen/).

Cohort caveat. Only 8w shares the controls' arrival cohort and sacrifice dates (phase 15a). 1w/2w were sacrificed
Jan–Mar 2019, 4w Nov–Dec 2018, controls and 8w Sep–Oct 2018, so every 1w/2w/4w-vs-control difference is also a
cohort, season and date difference.

Claim labels. [measured] = read from the CSVs here; [interpreted] = inference.

Reading of the 2026-09-26 run [measured unless marked]. Every set has size 0 or 1, so 1 − coverage is the empty-set
(abstention) rate plus a few wrong singletons. The control-trained model is best on 8w (full model accuracy 1.000,
coverage 0.889) and worst on 1w/2w (0.961/0.967, 0.867/0.878): the ordering runs opposite to training duration and
with cohort distance. At 1w and 2w, 5 of the 7 and 5 of the 6 full-model errors are VENACV → BAT, all in the
consortium-flagged females; without those vials 1w/2w sit level with 4w/8w. Held-out controls cover 0.92–0.94, all
trained groups 0.87–0.93; the only visible gap is between held-out controls and every trained group alike (−0.01 to
−0.05, fold sd ≈ 0.03, resting on 2 held-out controls per fold), not graded by duration. Calibrating on control+1w
does not lose coverage on 8w relative to control+1w+2w+4w (+0.01 to +0.02 over 10 splits, within the split-to-split
sd). With 10 animals per group a coverage change below about 0.05 is not detectable [interpreted]. Pooled-vial
calibration keeps the threshold finite with 3 calibration animals, but the 18 vials of one animal are not
exchangeable units, so the 1 − α guarantee is nominal; animal-level calibration (one score per animal) would need
at least 9 animals for a finite threshold at α = 0.10 [interpreted].

Outputs (results/15_time_course/fingerprint_by_duration/): part5_by_duration.csv, part5_matched_folds.csv,
part5_matched_summary.csv, part6_extrapolation.csv, part6_seeds.csv, part6_seeds_summary.csv, differences.csv,
calibration_info.csv, reproduction_check.csv, per_vial.csv.gz, fingerprint_time_course.png, NOTES.md (this text plus
the run's tables). About 3 min on 20 cores (100 tuned fits; the bootstrap is a sum over per-animal counts, so the
draws are cheap); `--quick` (2 splits, 200 draws) about 1 min.
"""
from __future__ import annotations

import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import beta

from tfp import cli, config as C, conformal as cp, io, models, report
from tfp.splits import assert_no_group_leak, fit_calibration_split, grouped_kfold, train_controls_test_trained

DURS = ["1w", "2w", "4w", "8w"]
PART6_DESIGNS = {"a_control+1w": ["control", "1w"], "b_control+1w+2w+4w": ["control", "1w", "2w", "4w"],
                 "c_control+4w": ["control", "4w"], "d_control_only": ["control"]}
METRICS = ["accuracy", "bal_acc", "coverage", "empty_rate", "mean_set_size", "singleton_rate"]
DIFF_METRICS = ["accuracy", "coverage", "empty_rate", "mean_set_size"]
TOL = 1e-6


# ------------------------------------------------------------------------------------------ consortium flags
def bat_contaminated_pids() -> set[str]:
    """Animals whose vena-cava vials the consortium flags 'BAT contamination' (female 1w/2w)."""
    o = io.load_outliers()
    if o is None:
        print("  WARNING: data/raw/outliers.csv absent; the *_excl_BATcontam_VENACV rows equal the unfiltered rows")
        return set()
    o = o[o["reason"].fillna("").str.contains("BAT contamination") & (o["tissue"] == "VENACV")]
    pids = set(o["pid"].astype(str))
    assert len(pids) == 10 and len(o) == 30, f"expected 10 animals / 30 vials flagged, found {len(pids)} / {len(o)}"
    return pids


def drop_contam(d: pd.DataFrame, pids: set[str]) -> pd.DataFrame:
    return d[~((d["y_true"] == "VENACV") & d["pid"].isin(pids))]


# ------------------------------------------------------------------------------------------ metrics, bootstrap
def animal_stats(d: pd.DataFrame, pids: np.ndarray, classes: list[str]) -> dict[str, np.ndarray]:
    """Per-animal sufficient statistics of a scored frame, rows in the order of `pids`; every metric of `metrics()`
    is a function of their sums, so a cluster bootstrap over animals is a sum over drawn rows."""
    pos = pd.Index(pids).get_indexer(d["pid"].to_numpy())
    assert (pos >= 0).all()
    n_p = len(pids)
    seen = d["seen"].to_numpy(dtype=bool)
    correct = d["y_true"].to_numpy() == d["y_pred"].to_numpy()
    size = d["set_size"].to_numpy(dtype=float)
    covered = d["covered"].to_numpy(dtype=bool)
    A = {k: np.bincount(pos, weights=w.astype(float), minlength=n_p) for k, w in (
        ("n", np.ones(len(d))), ("correct", correct), ("seen", seen), ("cov_seen", covered & seen),
        ("empty", size == 0), ("size", size), ("single", size == 1))}
    ci = pd.Index(classes).get_indexer(d["y_true"].to_numpy())
    m = seen & (ci >= 0)
    T, Cc = np.zeros((n_p, len(classes))), np.zeros((n_p, len(classes)))
    np.add.at(T, (pos[m], ci[m]), 1.0)
    np.add.at(Cc, (pos[m], ci[m]), correct[m].astype(float))
    A["T"], A["Cc"] = T, Cc
    return A


def metrics_from(A: dict[str, np.ndarray], idx: np.ndarray | None = None) -> dict:
    tot = {k: (v[idx] if idx is not None else v).sum(axis=0) for k, v in A.items()}
    present = tot["T"] > 0
    return {"accuracy": float(tot["correct"] / tot["n"]),
            "bal_acc": float((tot["Cc"][present] / tot["T"][present]).mean()) if present.any() else np.nan,
            "coverage": float(tot["cov_seen"] / tot["seen"]) if tot["seen"] else np.nan,
            "empty_rate": float(tot["empty"] / tot["n"]), "mean_set_size": float(tot["size"] / tot["n"]),
            "singleton_rate": float(tot["single"] / tot["n"])}


def metrics(d: pd.DataFrame, classes: list[str]) -> dict:
    return metrics_from(animal_stats(d, d["pid"].unique(), classes))


def summarize(d: pd.DataFrame, classes: list[str], rng: np.random.Generator, B: int, alpha: float) -> dict:
    """Point estimates with 95 % percentile cluster-bootstrap intervals over animals (B draws)."""
    pids = d["pid"].unique()
    A = animal_stats(d, pids, classes)
    point = metrics_from(A)
    boots = {k: np.empty(B) for k in METRICS}
    for b in range(B):
        m = metrics_from(A, rng.choice(len(pids), size=len(pids), replace=True))
        for k in METRICS:
            boots[k][b] = m[k]
    row = {"n_test_animals": len(pids), "n_test_vials": len(d), "n_unseen_vials": int((~d["seen"]).sum())}
    for k in METRICS:
        row[k] = point[k]
        row[f"{k}_lo"], row[f"{k}_hi"] = float(np.nanquantile(boots[k], 0.025)), float(np.nanquantile(boots[k], 0.975))
    x, n = int(A["cov_seen"].sum()), int(A["seen"].sum())     # Clopper–Pearson on vials (ignores clustering)
    row["coverage_cp_lo"] = float(beta.ppf(0.025, x, n - x + 1)) if x > 0 else 0.0
    row["coverage_cp_hi"] = float(beta.ppf(0.975, x + 1, n - x)) if x < n else 1.0
    has = A["seen"] > 0
    row["animals_cov_below_1ma"] = int((A["cov_seen"][has] / A["seen"][has] < 1 - alpha).sum())
    return row


def boot_diff(d1: pd.DataFrame, d2: pd.DataFrame, classes: list[str], rng: np.random.Generator, B: int, paired: bool) -> dict:
    """metric(d1) − metric(d2) with a cluster bootstrap over animals; paired = the same animals in both frames."""
    p1, p2 = np.array(sorted(d1["pid"].unique())), np.array(sorted(d2["pid"].unique()))
    if paired:
        assert np.array_equal(p1, p2), "paired difference needs the same animals in both frames"
    A1, A2 = animal_stats(d1, p1, classes), animal_stats(d2, p2, classes)
    m1, m2 = metrics_from(A1), metrics_from(A2)
    point = {k: m1[k] - m2[k] for k in DIFF_METRICS}
    boots = {k: np.empty(B) for k in DIFF_METRICS}
    for b in range(B):
        s1 = rng.choice(len(p1), size=len(p1), replace=True)
        s2 = s1 if paired else rng.choice(len(p2), size=len(p2), replace=True)
        b1, b2 = metrics_from(A1, s1), metrics_from(A2, s2)
        for k in DIFF_METRICS:
            boots[k][b] = b1[k] - b2[k]
    row = {"n_animals_1": len(p1), "n_animals_2": len(p2), "paired": paired}
    for k, v in point.items():
        row[f"d_{k}"] = v
        row[f"d_{k}_lo"], row[f"d_{k}_hi"] = float(np.nanquantile(boots[k], 0.025)), float(np.nanquantile(boots[k], 0.975))
    return row


# ------------------------------------------------------------------------------------------ fit, calibrate, score
def fit_and_calibrate(X, y, g, meta, fit_idx, cal_idx, k, alpha, model, prefilter, seed):
    assert_no_group_leak(meta, fit_idx, cal_idx)
    est = models.fit_tuned(model, X[fit_idx], y[fit_idx], g[fit_idx], k=k, prefilter=prefilter, quick=False, seed=seed)
    cls = list(est.classes_)
    scores, n_drop = cp.calibration_scores(est.predict_proba(X[cal_idx]), y[cal_idx], cls, "lac")
    qhat = cp.conformal_quantile(scores, alpha)
    n, n_cal_animals = len(scores), int(meta.iloc[cal_idx]["pid"].nunique())
    animal_rank = int(np.ceil((n_cal_animals + 1) * (1 - alpha) - 1e-9))
    info = {"n_fit_animals": int(meta.iloc[fit_idx]["pid"].nunique()), "n_fit_vials": len(fit_idx),
            "n_cal_animals": n_cal_animals, "n_cal_vials": n, "n_cal_dropped": n_drop,
            "conformal_rank": int(np.ceil((n + 1) * (1 - alpha) - 1e-9)), "qhat": qhat, "qhat_is_inf": bool(np.isinf(qhat)),
            "cal_score_max": float(scores.max()), "qhat_is_max_score": bool(np.isfinite(qhat) and qhat >= scores.max()),
            # what an animal-level calibration (one score per animal) would give with the same animals
            "animal_level_rank": animal_rank, "animal_level_qhat_inf": bool(animal_rank > n_cal_animals)}
    return est, cls, qhat, info


def score_vials(est, cls, qhat, X, y, meta, idx) -> pd.DataFrame:
    p = est.predict_proba(X[idx])
    sets = cp.predict_sets(p, qhat, "lac")
    yi, seen = cp.class_indices(y[idx], cls)
    covered = np.zeros(len(idx), dtype=bool)
    covered[seen] = sets[np.flatnonzero(seen), yi[seen]]
    m = meta.iloc[idx]
    return pd.DataFrame({"viallabel": m.index.astype(str).to_numpy(), "pid": m["pid"].astype(str).to_numpy(),
                         "group": m["group"].astype(str).to_numpy(), "sex": m["sex"].astype(str).to_numpy(),
                         "y_true": y[idx], "y_pred": est.predict(X[idx]), "covered": covered, "seen": seen,
                         "set_size": sets.sum(axis=1)})


# ------------------------------------------------------------------------------------------ reproduction check
def _reference(rel: str) -> Path | None:
    for root in (C.RESULTS_DIR, C.FROZEN_DIR):
        if (root / rel).exists():
            return root / rel
    return None


def reproduction_check(k: int, classes: list[str], primary_vials: pd.DataFrame, primary_qhat: dict, primary_ncal: dict,
                       part6b_vials: pd.DataFrame, part6b_qhat: dict, part6b_ncal: dict) -> pd.DataFrame:
    """Compare the primary design (pooled and per duration) and part-6 design (b) with phase 08 and its regeneration."""
    rows = []

    def add(ref, split, arm, group, metric, ours, theirs, checked=True):
        rows.append({"reference": ref, "split": split, "arm": arm, "test_group": group, "metric": metric, "ours": ours,
                     "reference_value": theirs, "abs_diff": abs(ours - theirs) if pd.notna(ours) and pd.notna(theirs) else np.nan,
                     "checked": checked})

    st_path = _reference("08_shift/TRNSCRPT/shift_table.csv")
    if st_path is None:
        print("  reproduction check vs 08_shift/TRNSCRPT/shift_table.csv skipped: file not found")
    else:
        st = pd.read_csv(st_path).set_index(["split", "arm"])
        for split, vials, qh, nc in (("train_control_test_trained", primary_vials, primary_qhat, primary_ncal),
                                     ("holdout_group_8w", part6b_vials, part6b_qhat, part6b_ncal)):
            for arm in ("full", f"panel_k{k}"):
                if (split, arm) not in st.index:
                    continue
                r, m = st.loc[(split, arm)], metrics(vials[vials["arm"] == arm], classes)
                add(st_path.name, split, arm, "pooled", "accuracy", m["accuracy"], float(r["accuracy_all"]))
                add(st_path.name, split, arm, "pooled", "coverage", m["coverage"], float(r["coverage_target_seen"]))
                add(st_path.name, split, arm, "pooled", "mean_set_size", m["mean_set_size"], float(r["avg_set_size_target"]))
                add(st_path.name, split, arm, "pooled", "empty_rate", m["empty_rate"], float(r["lac_frac_empty_target"]))
                add(st_path.name, split, arm, "pooled", "qhat", qh[arm], float(r["qhat"]), checked=False)
                add(st_path.name, split, arm, "pooled", "n_cal_vials", float(nc[arm]), float(r["n_cal"]))
    rg_path = _reference(f"31_site_regen/08_shift_k{k}/scores_target_vials.csv")
    if rg_path is None:
        print(f"  reproduction check vs 31_site_regen/08_shift_k{k}/scores_target_vials.csv skipped: file not found")
    else:
        rg = pd.read_csv(rg_path, dtype={"viallabel": str, "pid": str})
        rg = rg[rg["split"] == "train_control_test_trained"]
        for arm in ("full", f"panel_k{k}"):
            ours = primary_vials[primary_vials["arm"] == arm]
            theirs = rg[rg["arm"] == arm]
            for grp in DURS:
                o, t = ours[ours["group"] == grp], theirs[theirs["group"] == grp]
                add(rg_path.name, "train_control_test_trained", arm, grp, "accuracy", float((o["y_true"] == o["y_pred"]).mean()), float((t["y_pred"] == t["tissue"]).mean()))
                add(rg_path.name, "train_control_test_trained", arm, grp, "coverage", float(o.loc[o["seen"], "covered"].mean()), float(t.loc[t["seen"], "covered_marginal"].mean()))
                add(rg_path.name, "train_control_test_trained", arm, grp, "empty_rate", float((o["set_size"] == 0).mean()), float((t["size_marginal"] == 0).mean()))
            j = ours.merge(theirs[["viallabel", "y_pred", "covered_marginal", "size_marginal"]], on="viallabel", suffixes=("", "_ref"))
            add(rg_path.name, "train_control_test_trained", arm, "per_vial", "n_vials_compared", float(len(j)), float(len(ours)), checked=False)
            add(rg_path.name, "train_control_test_trained", arm, "per_vial", "n_y_pred_differ", float((j["y_pred"] != j["y_pred_ref"]).sum()), 0.0, checked=False)
            add(rg_path.name, "train_control_test_trained", arm, "per_vial", "n_covered_differ", float((j["covered"] != j["covered_marginal"]).sum()), 0.0, checked=False)
            add(rg_path.name, "train_control_test_trained", arm, "per_vial", "n_set_size_differ", float((j["set_size"] != j["size_marginal"]).sum()), 0.0, checked=False)
    chk = pd.DataFrame(rows)
    if len(chk):
        with pd.option_context("display.width", 200, "display.max_rows", 200):
            print("  reproduction check:")
            print(chk.to_string(index=False, float_format=lambda v: f"{v:.9g}"))
        core = chk[chk["checked"] & chk["metric"].isin(["accuracy", "coverage"])]
        bad = core[core["abs_diff"] > TOL]
        if len(bad):
            raise AssertionError(f"reproduction check failed (|Δ| > {TOL}):\n{bad.to_string(index=False)}")
        sets = chk[chk["metric"].str.startswith("n_") & chk["metric"].str.endswith("_differ")]
        n_diff = int(sets["ours"].sum()) if len(sets) else 0
        q = chk[chk["metric"] == "qhat"]
        print(f"  accuracy/coverage reproduce to {TOL:g}; per-vial sets differing: {n_diff}; "
              f"max |Δ qhat| = {q['abs_diff'].max():.2e}" + (" (sets identical, thresholds wobble)" if n_diff == 0 and q["abs_diff"].max() > 1e-9 else ""))
    return chk


# ------------------------------------------------------------------------------------------ figure
def figure(out: Path, p5: pd.DataFrame, ms: pd.DataFrame, s6: pd.DataFrame, k: int, n_seeds: int, alpha: float) -> Path:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    arms = ("full", f"panel_k{k}")
    col = {"full": "#2a78d6", f"panel_k{k}": "#eb6834"}
    lab = {"full": "full model", f"panel_k{k}": f"k = {k} panel"}
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.8), gridspec_kw={"width_ratios": [1.2, 1.2, 1]})
    order = ["control"] + DURS
    x = np.arange(len(order))
    for ax, metric, title in ((axes[0], "coverage", f"A  Conformal coverage (α = {alpha:g}), control-trained model"),
                              (axes[1], "accuracy", "B  Tissue accuracy, control-trained model")):
        for j, arm in enumerate(arms):
            off, c = (j - 0.5) * 0.22, col[arm]
            hc = ms[(ms["arm"] == arm) & ms["test_group"].str.startswith("control_heldout")].iloc[0]
            ax.errorbar(0 + off, hc[metric], yerr=[[hc[metric] - hc[f"{metric}_lo"]], [hc[f"{metric}_hi"] - hc[metric]]],
                        fmt="s", mfc="white", color=c, ms=8, lw=2, capsize=3)
            d = p5[(p5["design"] == "control_only_fit7_cal3") & (p5["arm"] == arm) & p5["test_group"].isin(DURS)].set_index("test_group").loc[DURS]
            ax.errorbar(x[1:] + off, d[metric], yerr=[d[metric] - d[f"{metric}_lo"], d[f"{metric}_hi"] - d[metric]],
                        fmt="o-", color=c, ms=8, lw=2, capsize=3, label=f"{lab[arm]} (fit 7 / cal 3 controls)")
            dx = p5[p5["design"].str.contains("excl") & (p5["arm"] == arm) & p5["test_group"].isin(["1w", "2w"])].set_index("test_group").loc[["1w", "2w"]]
            ax.plot(x[1:3] + off, dx[metric], "x", color=c, ms=9, mew=2,
                    label=f"{lab[arm]}, excl. consortium-flagged VENACV vials" if metric == "accuracy" else None)
        if metric == "coverage":
            ax.axhline(1 - alpha, color="#888", lw=1, ls="--")
            ax.text(4.35, 1 - alpha + 0.001, "1 − α", va="bottom", ha="right", fontsize=9, color="#555")
        ax.set_xticks(x, ["held-out\ncontrols*", "1w†", "2w†", "4w†", "8w"])
        ax.set_xlabel("test group (10 animals each, ~180 vials)")
        ax.set_ylabel(metric + " (vials)")
        ax.set_title(title, fontsize=10, loc="left")
        ax.grid(axis="y", color="#e5e5e5")
        ax.set_axisbelow(True)
        for s in ("top", "right"):
            ax.spines[s].set_visible(False)
    axes[0].set_ylim(0.78, 1.0)
    axes[1].set_ylim(0.84, 1.005)
    axes[0].legend(fontsize=8, loc="lower left", frameon=False)
    axes[1].legend(fontsize=8, loc="lower right", frameon=False)
    ax = axes[2]
    des = ["d_control_only", "a_control+1w", "c_control+4w", "b_control+1w+2w+4w"]
    lbl = ["control\n(10)", "control+1w\n(20)", "control+4w\n(20)", "control+1w\n+2w+4w (40)"]
    for j, arm in enumerate(arms):
        off = (j - 0.5) * 0.22
        d = s6[s6["arm"] == arm]
        for i, dn in enumerate(des):
            v = d[d["design"] == dn]["coverage"].to_numpy()
            ax.scatter(np.full(len(v), i + off) + np.linspace(-0.05, 0.05, len(v)), v, s=14, color=col[arm], alpha=0.35, lw=0)
        m = d.groupby("design")["coverage"].agg(["mean", "std"]).loc[des]
        ax.errorbar(np.arange(4) + off, m["mean"], yerr=m["std"].fillna(0), fmt="D", color=col[arm], ms=7, lw=2, capsize=3,
                    label=f"{lab[arm]}: mean ± sd over {n_seeds} fit/cal splits")
    ax.axhline(1 - alpha, color="#888", lw=1, ls="--")
    ax.set_xticks(np.arange(4), lbl, fontsize=8)
    ax.set_xlabel("fit+cal groups (n animals); test: 10 held-out 8w animals")
    ax.set_ylabel("coverage on 8w vials")
    ax.set_ylim(0.84, 0.98)
    ax.set_title("C  Held-out 8w: calibration source vs coverage", fontsize=10, loc="left")
    ax.legend(fontsize=8, loc="upper right", frameon=False)
    ax.grid(axis="y", color="#e5e5e5")
    ax.set_axisbelow(True)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    fig.text(0.01, 0.005, "Bars: 95% cluster-bootstrap over test animals (conditional on one fit/cal split). "
             "* held-out controls: 5 matched folds (fit 6 / cal 2), pooled. "
             "† 1w/2w/4w differ from controls in arrival cohort, sacrifice season and date as well as training; only 8w is cohort-matched.",
             fontsize=7.5, color="#444")
    fig.tight_layout(rect=(0, 0.03, 1, 1))
    path = out / "fingerprint_time_course.png"
    fig.savefig(path, dpi=150)
    plt.close(fig)
    return path


# ------------------------------------------------------------------------------------------ main
def main() -> None:
    ap = cli.common_parser("The tissue fingerprint by training duration (parts 5 and 6 of the time-course investigation)")
    ap.add_argument("--k", type=int, default=20, help="panel size (the certified k of phase 06)")
    ap.add_argument("--alpha", type=float, default=0.1)
    ap.add_argument("--model", default="logreg_l2")
    ap.add_argument("--seeds", type=int, default=10, help="fit/calibration splits for the part-6 robustness run (--quick: 2)")
    ap.add_argument("--bootstraps", type=int, default=2000, help="cluster-bootstrap draws (--quick: 200)")
    ap.add_argument("--no-figure", action="store_true")
    args = ap.parse_args()
    if args.quick:
        args.seeds, args.bootstraps = min(args.seeds, 2), min(args.bootstraps, 200)
    cli.banner("15_fingerprint_by_duration", args)
    out = cli.outdir("15_time_course/fingerprint_by_duration", args.out)
    k, alpha, B, seed = args.k, args.alpha, args.bootstraps, args.seed
    arms = (("full", None), (f"panel_k{k}", k))
    t0 = time.time()

    pheno = io.load_pheno()
    om = io.stack_tissues(args.assay, tissues=cli.parse_tissues(args.tissues), source=cli.resolve_source(args.assay, args.source),
                          join="inner", pheno=pheno)
    X, y, g, meta = om.X.to_numpy(dtype=float), om.meta["tissue"].astype(str).to_numpy(), om.groups(), om.meta
    grp = meta["group"].astype(str).to_numpy()
    classes = sorted(np.unique(y))
    animals_by_group = meta.groupby("group")["pid"].nunique().reindex(C.GROUP_ORDER).to_dict()
    print(f"loaded {X.shape[0]} vials × {X.shape[1]} genes, {len(classes)} tissues, {meta['pid'].nunique()} animals "
          f"{animals_by_group} ({time.time() - t0:.0f}s)", flush=True)
    contam = bat_contaminated_pids()
    print(f"consortium 'BAT contamination' VENACV animals: {len(contam)}", flush=True)
    rng = np.random.default_rng(seed)
    cal_rows, vial_frames, diffs = [], [], []

    # ---------------- PART 5 primary: controls -> each duration --------------------------------------------
    [(_, (tr, te))] = list(train_controls_test_trained(meta))
    fit_idx, cal_idx = fit_calibration_split(meta, tr, 0.3, seed)          # exactly as phase 08's eval_split
    p5, primary_qhat, primary_ncal = [], {}, {}
    for arm, kk in arms:
        est, cls, qhat, info = fit_and_calibrate(X, y, g, meta, fit_idx, cal_idx, kk, alpha, args.model, args.prefilter, seed)
        cal_rows.append({"design": "p5_primary_control_only", "arm": arm, **info})
        primary_qhat[arm], primary_ncal[arm] = qhat, info["n_cal_vials"]
        dv = score_vials(est, cls, qhat, X, y, meta, te)
        vial_frames.append(dv.assign(design="p5_primary_control_only", arm=arm))
        pooled = metrics(dv, classes)
        print(f"  {arm}: fit {info['n_fit_animals']} / cal {info['n_cal_animals']} control animals ({info['n_cal_vials']} vials), "
              f"pooled trained acc={pooled['accuracy']:.4f} cov={pooled['coverage']:.4f} size={pooled['mean_set_size']:.4f} qhat={qhat:.6f}", flush=True)
        for dur in DURS + ["pooled_trained"]:
            d = dv if dur == "pooled_trained" else dv[dv["group"] == dur]
            common = {"qhat": qhat, "n_cal_vials": info["n_cal_vials"], "n_cal_animals": info["n_cal_animals"]}
            p5.append({"design": "control_only_fit7_cal3", "arm": arm, "test_group": dur, **summarize(d, classes, rng, B, alpha), **common})
            p5.append({"design": "control_only_fit7_cal3_excl_BATcontam_VENACV", "arm": arm, "test_group": dur,
                       **summarize(drop_contam(d, contam), classes, rng, B, alpha), **common})
        d8 = dv[dv["group"] == "8w"]
        for excl in (False, True):
            for dur in ("1w", "2w", "4w"):
                d1 = dv[dv["group"] == dur]
                diffs.append({"part": 5, "arm": arm, "comparison": f"{dur} - 8w (control-only model)", "excl_BATcontam_VENACV": excl,
                              **boot_diff(drop_contam(d1, contam) if excl else d1, d8, classes, rng, B, paired=False)})
    p5 = pd.DataFrame(p5)
    p5.to_csv(out / "part5_by_duration.csv", index=False)
    primary_vials = pd.concat([f for f in vial_frames if f["design"].iloc[0] == "p5_primary_control_only"], ignore_index=True)
    print(f"part 5 primary done ({time.time() - t0:.0f}s)", flush=True)

    # ---------------- PART 5 reference arm: matched folds over the controls -------------------------------
    ctrl = np.flatnonzero(grp == "control")
    mc = meta.iloc[ctrl].reset_index(drop=True)
    mf_rows, held_ctrl = [], {a: [] for a, _ in arms}
    for fold, (ftr, fte) in enumerate(grouped_kfold(mc, "sex", n_splits=5, seed=seed)):
        tr_c, te_c = ctrl[ftr], ctrl[fte]
        assert_no_group_leak(meta, tr_c, te_c)
        f_idx, c_idx = fit_calibration_split(meta, tr_c, 0.3, seed + fold)
        test_idx = np.concatenate([te_c, np.flatnonzero(grp != "control")])
        for arm, kk in arms:
            est, cls, qhat, info = fit_and_calibrate(X, y, g, meta, f_idx, c_idx, kk, alpha, args.model, args.prefilter, seed)
            cal_rows.append({"design": f"p5_matched_fold{fold}", "arm": arm, **info})
            dv = score_vials(est, cls, qhat, X, y, meta, test_idx).assign(fold=fold)
            dv["test_group"] = np.where(np.arange(len(test_idx)) < len(te_c), "control_heldout", dv["group"])
            held_ctrl[arm].append(dv[dv["test_group"] == "control_heldout"])
            vial_frames.append(dv.assign(design=f"p5_matched_fold{fold}", arm=arm))
            for tg, d in dv.groupby("test_group"):
                mf_rows.append({"fold": fold, "arm": arm, "test_group": tg, "n_test_animals": d["pid"].nunique(), "n_test_vials": len(d),
                                "n_fit_animals": info["n_fit_animals"], "n_cal_animals": info["n_cal_animals"],
                                "n_cal_vials": info["n_cal_vials"], "qhat": qhat, **metrics(d, classes)})
        print(f"  matched fold {fold} done ({time.time() - t0:.0f}s)", flush=True)
    mf = pd.DataFrame(mf_rows)
    mf.to_csv(out / "part5_matched_folds.csv", index=False)
    ms = []
    for arm, _ in arms:
        hc = pd.concat(held_ctrl[arm], ignore_index=True)                  # every control animal tested once
        ms.append({"arm": arm, "test_group": "control_heldout (pooled over folds)", "n_folds": mf["fold"].nunique(),
                   **summarize(hc, classes, rng, B, alpha)})
        sub = mf[(mf["arm"] == arm) & (mf["test_group"] != "control_heldout")]
        for tg, s in sub.groupby("test_group"):
            r = {"arm": arm, "test_group": f"{tg} (mean over folds)", "n_folds": len(s),
                 "n_test_animals": int(s["n_test_animals"].iloc[0]), "n_test_vials": int(s["n_test_vials"].iloc[0])}
            for key in METRICS:
                r[key], r[f"{key}_fold_sd"] = float(s[key].mean()), float(s[key].std(ddof=1))
                r[f"{key}_fold_min"], r[f"{key}_fold_max"] = float(s[key].min()), float(s[key].max())
            ms.append(r)
        piv = mf[mf["arm"] == arm].pivot(index="fold", columns="test_group", values="coverage")
        for tg in DURS:                                                     # paired within fold: same model, same calibration
            dlt = piv[tg] - piv["control_heldout"]
            ms.append({"arm": arm, "test_group": f"coverage[{tg}] - coverage[control_heldout], per fold", "n_folds": len(dlt),
                       "coverage": float(dlt.mean()), "coverage_fold_sd": float(dlt.std(ddof=1)),
                       "coverage_fold_min": float(dlt.min()), "coverage_fold_max": float(dlt.max())})
    ms = pd.DataFrame(ms)
    ms.to_csv(out / "part5_matched_summary.csv", index=False)
    print(f"part 5 matched done ({time.time() - t0:.0f}s)", flush=True)

    # ---------------- PART 6: temporal extrapolation of the calibration to held-out 8w --------------------
    te8 = np.flatnonzero(grp == "8w")
    p6, p6_qhat_b, p6_ncal_b = [], {}, {}
    for dname, groups in PART6_DESIGNS.items():
        trg = np.flatnonzero(np.isin(grp, groups))
        assert_no_group_leak(meta, trg, te8)
        f_idx, c_idx = fit_calibration_split(meta, trg, 0.3, seed)
        for arm, kk in arms:
            est, cls, qhat, info = fit_and_calibrate(X, y, g, meta, f_idx, c_idx, kk, alpha, args.model, args.prefilter, seed)
            cal_rows.append({"design": f"p6_{dname}", "arm": arm, **info})
            dv = score_vials(est, cls, qhat, X, y, meta, te8)
            vial_frames.append(dv.assign(design=f"p6_{dname}", arm=arm))
            if dname == "b_control+1w+2w+4w":
                p6_qhat_b[arm], p6_ncal_b[arm] = qhat, info["n_cal_vials"]
            cal_groups = meta.iloc[c_idx].groupby("group")["pid"].nunique().to_dict()
            p6.append({"design": dname, "arm": arm, "train_groups": "+".join(groups), "n_train_animals": int(meta.iloc[trg]["pid"].nunique()),
                       "n_fit_animals": info["n_fit_animals"], "n_cal_animals": info["n_cal_animals"], "n_cal_vials": info["n_cal_vials"],
                       "cal_animals_by_group": str(cal_groups), "qhat": qhat, **summarize(dv, classes, rng, B, alpha)})
            print(f"  {dname} {arm}: acc={p6[-1]['accuracy']:.4f} cov={p6[-1]['coverage']:.4f} size={p6[-1]['mean_set_size']:.3f} "
                  f"({time.time() - t0:.0f}s)", flush=True)
    p6 = pd.DataFrame(p6)
    p6.to_csv(out / "part6_extrapolation.csv", index=False)
    vp6 = pd.concat([f for f in vial_frames if f["design"].iloc[0].startswith("p6_")], ignore_index=True)
    for arm, _ in arms:
        ref = vp6[(vp6["design"] == "p6_b_control+1w+2w+4w") & (vp6["arm"] == arm)]
        for dname in ("a_control+1w", "c_control+4w", "d_control_only"):
            cur = vp6[(vp6["design"] == f"p6_{dname}") & (vp6["arm"] == arm)]
            diffs.append({"part": 6, "arm": arm, "comparison": f"{dname} - b_control+1w+2w+4w (same 8w animals)",
                          "excl_BATcontam_VENACV": False, **boot_diff(cur, ref, classes, rng, B, paired=True)})
    diffs = pd.DataFrame(diffs)
    diffs.to_csv(out / "differences.csv", index=False)
    cal_info = pd.DataFrame(cal_rows)
    cal_info.to_csv(out / "calibration_info.csv", index=False)
    pd.concat(vial_frames, ignore_index=True).to_csv(out / "per_vial.csv.gz", index=False)
    print(f"part 6 done ({time.time() - t0:.0f}s)", flush=True)

    # ---------------- reproduction check against phase 08 and its regeneration --------------------------------
    part6b_vials = vp6[vp6["design"] == "p6_b_control+1w+2w+4w"]
    chk = reproduction_check(k, classes, primary_vials, primary_qhat, primary_ncal, part6b_vials, p6_qhat_b, p6_ncal_b)
    chk.to_csv(out / "reproduction_check.csv", index=False)

    # ---------------- PART 6 robustness: random fit/calibration splits ----------------------------------------
    rows = []
    for s in range(1, args.seeds + 1):
        for dname, groups in PART6_DESIGNS.items():
            trg = np.flatnonzero(np.isin(grp, groups))
            f_idx, c_idx = fit_calibration_split(meta, trg, 0.3, seed + s)
            for arm, kk in arms:
                est, cls, qhat, info = fit_and_calibrate(X, y, g, meta, f_idx, c_idx, kk, alpha, args.model, args.prefilter, seed)
                dv = score_vials(est, cls, qhat, X, y, meta, te8)
                rows.append({"split_seed": seed + s, "design": dname, "arm": arm, "n_fit_animals": info["n_fit_animals"],
                             "n_cal_animals": info["n_cal_animals"], "n_cal_vials": info["n_cal_vials"], "qhat": qhat,
                             "n_test_animals": dv["pid"].nunique(), "n_test_vials": len(dv), **metrics(dv, classes)})
        print(f"  split {s}/{args.seeds} done ({time.time() - t0:.0f}s)", flush=True)
    s6 = pd.DataFrame(rows)
    s6.to_csv(out / "part6_seeds.csv", index=False)
    summ = []
    for (dname, arm), d in s6.groupby(["design", "arm"]):
        r = {"design": dname, "arm": arm, "n_splits": len(d), "n_test_animals": int(d["n_test_animals"].iloc[0]),
             "n_cal_animals": int(d["n_cal_animals"].iloc[0])}
        for key in ("accuracy", "coverage", "empty_rate", "mean_set_size", "qhat"):
            r[f"{key}_mean"], r[f"{key}_sd"] = float(d[key].mean()), float(d[key].std(ddof=1))
            r[f"{key}_min"], r[f"{key}_max"] = float(d[key].min()), float(d[key].max())
        summ.append(r)
    for arm, _ in arms:                                                     # paired over splits: design minus (b)
        piv = s6[s6["arm"] == arm].pivot(index="split_seed", columns="design", values="coverage")
        for dname in ("a_control+1w", "c_control+4w", "d_control_only"):
            dd = piv[dname] - piv["b_control+1w+2w+4w"]
            summ.append({"design": f"coverage[{dname}] - coverage[b]", "arm": arm, "n_splits": len(dd),
                         "coverage_mean": float(dd.mean()), "coverage_sd": float(dd.std(ddof=1)),
                         "coverage_min": float(dd.min()), "coverage_max": float(dd.max()), "n_splits_positive": int((dd > 0).sum())})
    s6s = pd.DataFrame(summ)
    s6s.to_csv(out / "part6_seeds_summary.csv", index=False)

    fig_path = None if args.no_figure else figure(out, p5, ms, s6, k, args.seeds, alpha)
    elapsed = time.time() - t0

    # ---------------- NOTES.md: the docstring plus the run's tables --------------------------------------------
    show5 = ["design", "arm", "test_group", "n_test_animals", "n_test_vials", "accuracy", "accuracy_lo", "accuracy_hi", "coverage",
             "coverage_lo", "coverage_hi", "empty_rate", "mean_set_size", "qhat"]
    showm = ["arm", "test_group", "n_folds", "n_test_animals", "accuracy", "coverage", "coverage_lo", "coverage_hi", "coverage_fold_sd",
             "coverage_fold_min", "coverage_fold_max", "empty_rate"]
    show6 = ["design", "arm", "n_train_animals", "n_fit_animals", "n_cal_animals", "n_cal_vials", "cal_animals_by_group", "qhat", "accuracy",
             "accuracy_lo", "accuracy_hi", "coverage", "coverage_lo", "coverage_hi", "empty_rate"]
    show6s = ["design", "arm", "n_splits", "accuracy_mean", "coverage_mean", "coverage_sd", "coverage_min", "coverage_max", "qhat_mean", "n_splits_positive"]
    showd = ["part", "arm", "comparison", "excl_BATcontam_VENACV", "d_accuracy", "d_accuracy_lo", "d_accuracy_hi", "d_coverage", "d_coverage_lo", "d_coverage_hi"]
    showc = ["design", "arm", "n_fit_animals", "n_cal_animals", "n_cal_vials", "conformal_rank", "qhat", "cal_score_max", "animal_level_rank", "animal_level_qhat_inf"]
    per_vial_sizes = pd.concat(vial_frames)["set_size"].value_counts().sort_index().to_dict()
    run = [
        "## Run [measured]",
        f"`{Path(__file__).name}` on {time.strftime('%Y-%m-%d %H:%M')}: assay {args.assay}, source {cli.resolve_source(args.assay, args.source)}, "
        f"{X.shape[0]} vials × {X.shape[1]} genes, {len(classes)} tissues, {meta['pid'].nunique()} animals {animals_by_group}; model `{args.model}`, "
        f"prefilter {args.prefilter}, k = {k}, α = {alpha:g}, {B} bootstrap draws, {args.seeds} fit/calibration splits in part 6, seed {seed}"
        f"{' (--quick)' if args.quick else ''}; {elapsed / 60:.1f} min. Set sizes over every scored vial: {per_vial_sizes}. "
        f"Consortium-flagged VENACV animals excluded in the `_excl_` rows: {len(contam)}.",
        "### Reproduction check against phase 08 and its regeneration\n\n" + (report.df_to_md(chk.assign(**{c: chk[c].map(lambda v: f'{v:.9g}' if pd.notna(v) else '') for c in ('ours', 'reference_value')}), floatfmt=".2e", max_rows=100) if len(chk) else "_skipped: reference files not found_"),
        "### Calibration sets (`calibration_info.csv`)\n\n" + report.df_to_md(cal_info[showc], floatfmt=".4f"),
        "### Part 5, control-trained model by duration (`part5_by_duration.csv`)\n\n" + report.df_to_md(p5[show5], floatfmt=".4f"),
        "### Part 5, matched folds (`part5_matched_summary.csv`)\n\n" + report.df_to_md(ms[showm], floatfmt=".4f"),
        "### Part 6, temporal extrapolation, single split (`part6_extrapolation.csv`)\n\n" + report.df_to_md(p6[show6], floatfmt=".4f"),
        f"### Part 6, {args.seeds} fit/calibration splits (`part6_seeds_summary.csv`)\n\n" + report.df_to_md(s6s[show6s], floatfmt=".4f"),
        "### Bootstrap differences (`differences.csv`)\n\n" + report.df_to_md(diffs[showd], floatfmt=".4f"),
    ]
    if fig_path is not None:
        run.append(report.figure_md(fig_path, "A, B: coverage and accuracy by duration; C: calibration source vs coverage on held-out 8w"))
    notes = "# Phase 15b: the tissue fingerprint by training duration\n\n" + __doc__.split("\n", 1)[1].strip() + "\n\n" + "\n\n".join(run) + "\n"
    (out / "NOTES.md").write_text(notes)
    report.add_section("15b · Tissue fingerprint by training duration", "\n\n".join(run[1:4]),
                       params={"k": k, "alpha": alpha, "seeds": args.seeds, "bootstraps": B, "quick": args.quick})
    print(f"done ({elapsed:.0f}s) -> {out}", flush=True)


if __name__ == "__main__":
    main()
