"""PART 4 — physiology anchor: does an animal's molecular "trained-ness" track its own physiological change?

Run from code/pipeline:  PYTHONPATH=src python investigations/time_course/tc4_physiology.py [--quick]
Outputs: results/15_time_course/4_physiology/

PRE-SPECIFIED (written before any correlation was computed)
- Phenotype source: data/raw/pheno.csv (dotted column names, portal/R-package c1.0 PHENO; NOT the c4.0 package).
  Index _1 = pre-training test, _2 = post-training test (verified from days_vo2_1 < training.day1_days < days_vo2_2).
  calculated.variables.vo2_max_change == vo2_max_2 - vo2_max_1 and pct_body_fat_change == nmr_fat_2 - nmr_fat_1
  (checked in this script). 1w / 2w animals have NO post-training VO2max or NMR measurement.
- Molecular score (primary): out-of-fold logit of P(trained) from logreg_l2 (C = 0.1 fixed, no tuning) on
  median-impute -> top-2000-variance -> z-score features + a sex indicator (as Task B, add_sex=True), per omic
  (TRNSCRPT, PROT, METAB), trained control vs Nw on the Task B animals (all three omics; build_blocks from
  scripts/07_fusion_vs_baselines.py), folds grouped_kfold(meta, "sex_group", min(4, n // 3), seed); the OOF score
  is averaged over 10 fold seeds (SEED + r). Secondary scores: centroid-axis projection and Euclidean distance
  from the control centroid (both on the same preprocessing, control centroid from training-fold animals only,
  so never includes the scored animal).
- Primary test family: WITHIN the 8w group, sex-adjusted (within-sex ranks, permutation within sex) correlation
  between the primary score and {vo2max_change, fat_change}: 7 tissues x 3 omics x 2 variables = 42 tests;
  BH-FDR and Bonferroni over those 42. Everything else is secondary / descriptive.
- Pooled (control + 8w) correlations are reported but are NOT evidence about batch vs training: group membership
  drives both the score and the physiology, so any group separator correlates.
- Negative control: the same within-group correlation in the CONTROL animals (their score from the same OOF model).
"""
from __future__ import annotations

import argparse
import importlib.util
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd
from scipy import stats
from sklearn.metrics import roc_auc_score

from tfp import cli, config as C, io, models
from tfp.splits import assert_no_group_leak, grouped_kfold

HERE = Path(__file__).resolve()
PIPE = HERE.parents[2]
_spec = importlib.util.spec_from_file_location("fus07", PIPE / "scripts" / "07_fusion_vs_baselines.py")
fus07 = importlib.util.module_from_spec(_spec)
sys.modules["fus07"] = fus07
_spec.loader.exec_module(fus07)

TISSUES = ["CORTEX", "HEART", "KIDNEY", "LIVER", "LUNG", "SKM-GN", "WAT-SC"]
ASSAYS = ["TRNSCRPT", "PROT", "METAB"]
GROUP_MAP = {"Eight-week program Control Group": "control", "One-week program": "1w", "Two-week program": "2w",
             "Four-week program": "4w", "Eight-week program Training Group": "8w"}
PHYS_PRIMARY = ["vo2max_change", "fat_change"]
PHYS_ALL = PHYS_PRIMARY + ["lean_change", "weight_change", "vo2max_pre", "fat_pre"]
PREFILTER = 2000
N_FOLD_SEEDS = 10
N_PERM = 20000


# ------------------------------------------------------------------------------------------------ physiology
def physiology() -> pd.DataFrame:
    p = pd.read_csv(C.RAW_DIR / "pheno.csv", dtype=str, low_memory=False).drop_duplicates("pid")
    num = lambda c: pd.to_numeric(p[c], errors="coerce")
    d = pd.DataFrame({"pid": p["pid"].str.strip(), "bid": p["bid"].str.strip(),
                      "sex": p["registration.sex"].str.strip().str.lower(),
                      "group": p["key.anirandgroup"].map(GROUP_MAP)})
    d["days_vo2_pre"], d["days_vo2_post"] = num("vo2.max.test.days_vo2_1"), num("vo2.max.test.days_vo2_2")
    d["days_train_start"] = num("training.day1_days")
    d["days_nmr_pre"], d["days_nmr_post"] = num("nmr.testing.days_nmr_1"), num("nmr.testing.days_nmr_2")
    d["vo2max_pre"], d["vo2max_post"] = num("vo2.max.test.vo2_max_1"), num("vo2.max.test.vo2_max_2")
    d["vo2max_change"] = d["vo2max_post"] - d["vo2max_pre"]
    d["fat_pre"], d["fat_post"] = num("nmr.testing.nmr_fat_1"), num("nmr.testing.nmr_fat_2")
    d["fat_change"] = d["fat_post"] - d["fat_pre"]
    d["lean_change"] = num("nmr.testing.nmr_lean_2") - num("nmr.testing.nmr_lean_1")
    d["weight_change"] = num("nmr.testing.nmr_weight_2") - num("nmr.testing.nmr_weight_1")
    # consistency with the consortium's calculated variables
    for mine, theirs in [("vo2max_change", "calculated.variables.vo2_max_change"),
                         ("fat_change", "calculated.variables.pct_body_fat_change"),
                         ("lean_change", "calculated.variables.pct_body_lean_change")]:
        diff = (d[mine] - num(theirs)).abs().max()
        assert not (diff > 1e-6), f"{mine} disagrees with {theirs} (max |diff| {diff})"
    # pre/post ordering relative to training start (only where training exists; controls have no training days)
    ok = d.dropna(subset=["days_vo2_post"])
    assert (ok["days_vo2_pre"] < ok["days_vo2_post"]).all()
    tr = d.dropna(subset=["days_train_start", "days_vo2_post"])
    assert ((tr["days_vo2_pre"] < tr["days_train_start"]) & (tr["days_train_start"] < tr["days_vo2_post"])).all()
    return d


def iqr_str(x):
    x = x.dropna()
    if len(x) == 0:
        return np.nan, np.nan, np.nan
    return x.median(), x.quantile(0.25), x.quantile(0.75)


def physiology_tables(ph: pd.DataFrame, out: Path) -> pd.DataFrame:
    rows = []
    for var in ["vo2max_pre", "vo2max_post", "vo2max_change", "fat_pre", "fat_post", "fat_change", "lean_change",
                "weight_change", "days_vo2_pre", "days_vo2_post", "days_train_start"]:
        for grp in ["control", "1w", "2w", "4w", "8w"]:
            for sex in ["all", "female", "male"]:
                s = ph[(ph["group"] == grp) & ((ph["sex"] == sex) | (sex == "all"))][var]
                med, q1, q3 = iqr_str(s)
                rows.append({"variable": var, "group": grp, "sex": sex, "n_animals": int(s.notna().sum()),
                             "n_in_group": int(len(s)), "median": med, "q25": q1, "q75": q3})
    tab = pd.DataFrame(rows)
    tab.to_csv(out / "physiology_by_group.csv", index=False)
    tests = []
    for var in ["vo2max_change", "fat_change", "lean_change", "weight_change", "vo2max_pre", "fat_pre"]:
        for grp in ["8w", "4w"]:
            for sex in ["female", "male"]:
                a = ph[(ph["group"] == grp) & (ph["sex"] == sex)][var].dropna()
                b = ph[(ph["group"] == "control") & (ph["sex"] == sex)][var].dropna()
                u = stats.mannwhitneyu(a, b, alternative="two-sided")
                tests.append({"variable": var, "contrast": f"{grp} vs control", "sex": sex, "n_trained": len(a),
                              "n_control": len(b), "median_diff": a.median() - b.median(), "mw_p": u.pvalue,
                              "cohort_matched": grp == "8w"})
    tt = pd.DataFrame(tests)
    tt.to_csv(out / "physiology_group_tests.csv", index=False)
    return tab


# ------------------------------------------------------------------------------------------------ OOF scores
def oof_scores(X: pd.DataFrame, meta: pd.DataFrame, dur: str, n_seeds: int) -> tuple[pd.DataFrame, dict]:
    Xs = fus07.with_sex(X, meta).to_numpy(dtype=float)
    y = (meta["group"].astype(str) == dur).astype(int).to_numpy()
    n_splits = min(4, meta["pid"].nunique() // 3)
    acc = {k: np.zeros((n_seeds, len(y))) for k in ["logit", "proj", "dist"]}
    for r in range(n_seeds):
        filled = np.zeros(len(y), bool)
        for tr, te in grouped_kfold(meta, "sex_group", n_splits, C.SEED + r):
            assert_no_group_leak(meta, tr, te)
            pipe = models.make_pipeline("logreg_l2", prefilter=PREFILTER, seed=C.SEED)
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                pipe.fit(Xs[tr], y[tr])
            acc["logit"][r, te] = pipe.decision_function(Xs[te])  # log-odds of the trained class
            pre = pipe[:-1]  # impute -> prefilter -> scale (-> passthrough), fit on the training fold
            Ztr, Zte = pre.transform(Xs[tr]), pre.transform(Xs[te])
            mc, mt = Ztr[y[tr] == 0].mean(0), Ztr[y[tr] == 1].mean(0)
            w = mt - mc
            acc["proj"][r, te] = (Zte - (mc + mt) / 2) @ w / np.linalg.norm(w)
            acc["dist"][r, te] = np.linalg.norm(Zte - mc, axis=1) / np.sqrt(Zte.shape[1])
            filled[te] = True
        assert filled.all()
    res = pd.DataFrame({"pid": meta["pid"].astype(str).to_numpy(), "bid": (meta["bid"] if "bid" in meta.columns else meta.index.to_series()).astype(str).to_numpy(),
                        "sex": meta["sex"].astype(str).to_numpy(), "group": meta["group"].astype(str).to_numpy()})
    info = {"n_splits": n_splits, "n_features_in": X.shape[1]}
    for k, v in acc.items():
        res[f"score_{k}"] = v.mean(0)
        res[f"score_{k}_seed_sd"] = v.std(0)
        info[f"oof_auroc_{k}"] = float(np.mean([roc_auc_score(y, v[r]) for r in range(n_seeds)]))
    return res, info


# ------------------------------------------------------------------------------------------------ correlation
def within_sex_ranks(v: np.ndarray, sex: np.ndarray) -> np.ndarray:
    out = np.empty(len(v))
    for s in np.unique(sex):
        m = sex == s
        out[m] = stats.rankdata(v[m]) / (m.sum() + 1) - 0.5
    return out


def perm_corr(x, yv, sex, adjust: bool, rng, n_perm=N_PERM):
    """Spearman (adjust=False) or sex-stratified rank correlation (adjust=True; ranks within sex, permute within sex).
    Two-sided permutation p with the +1 correction."""
    if adjust:
        rx, ry = within_sex_ranks(x, sex), within_sex_ranks(yv, sex)
    else:
        rx, ry = stats.rankdata(x), stats.rankdata(yv)
    rx = (rx - rx.mean()) / np.linalg.norm(rx - rx.mean())
    ry = (ry - ry.mean()) / np.linalg.norm(ry - ry.mean())
    obs = float(rx @ ry)
    P = np.tile(np.arange(len(ry)), (n_perm, 1))
    if adjust:
        for s in np.unique(sex):
            idx = np.where(sex == s)[0]
            P[:, idx] = idx[rng.permuted(np.tile(np.arange(len(idx)), (n_perm, 1)), axis=1)]
    else:
        P = rng.permuted(P, axis=1)
    null = ry[P] @ rx
    p = (1 + np.sum(np.abs(null) >= abs(obs) - 1e-12)) / (n_perm + 1)
    return obs, p


def critical_rho(n: int, alpha=0.05, n_perm=200000, seed=0) -> float:
    """Exact-by-simulation two-sided critical |Spearman rho| under H0 at sample size n."""
    rng = np.random.default_rng(seed)
    r = np.arange(n, dtype=float)
    r = (r - r.mean()) / np.linalg.norm(r - r.mean())
    null = np.abs(r[rng.permuted(np.tile(np.arange(n), (n_perm, 1)), axis=1)] @ r)
    return float(np.quantile(null, 1 - alpha))


def rho_for_power(n: int, power=0.8, alpha=0.05) -> float:
    """Fisher-z approximation for Spearman (variance 1.06/(n-3)): true |rho| needed for the given power."""
    z = (stats.norm.ppf(1 - alpha / 2) + stats.norm.ppf(power)) * np.sqrt(1.06 / (n - 3))
    return float(np.tanh(z))


def bh(p):
    p = np.asarray(p, float)
    o = np.argsort(p)
    q = p[o] * len(p) / np.arange(1, len(p) + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty_like(q)
    out[o] = np.minimum(q, 1)
    return out


# ------------------------------------------------------------------------------------------------ figure
def figure(scores: pd.DataFrame, corr: pd.DataFrame, path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    col = {"female": "#eb6834", "male": "#2a78d6"}
    mk = {"female": "o", "male": "s"}
    rows = [(a, v) for v in PHYS_PRIMARY for a in ASSAYS]
    fig, axes = plt.subplots(len(rows), len(TISSUES), figsize=(2.3 * len(TISSUES), 2.0 * len(rows)), squeeze=False)
    lab = {"vo2max_change": "VO2max change", "fat_change": "% fat change (pts)"}
    for i, (a, v) in enumerate(rows):
        for j, t in enumerate(TISSUES):
            ax = axes[i, j]
            s = scores[(scores["tissue"] == t) & (scores["assay"] == a) & (scores["duration"] == "8w") & (scores["group"] == "8w")]
            for sx in ["female", "male"]:
                ss = s[s["sex"] == sx]
                ax.scatter(ss["score_logit"], ss[v], s=22, c=col[sx], marker=mk[sx], edgecolors="white", linewidths=0.8,
                           label=sx if (i == 0 and j == 0) else None)
            c = corr[(corr["tissue"] == t) & (corr["assay"] == a) & (corr["physio"] == v) & (corr["set"] == "within_8w")
                     & (corr["score"] == "logit") & (corr["duration"] == "8w") & (corr["adjust"] == "sex_ranks")]
            if len(c):
                c = c.iloc[0]
                ax.set_title(f"{t} {a}\nrho_sex={c['rho']:+.2f} p={c['p_perm']:.2f} n={int(c['n'])}", fontsize=7)
            ax.tick_params(labelsize=6)
            ax.grid(alpha=0.2, lw=0.5)
            for sp in ["top", "right"]:
                ax.spines[sp].set_visible(False)
            if j == 0:
                ax.set_ylabel(lab[v], fontsize=7)
            if i == len(rows) - 1:
                ax.set_xlabel("OOF logit P(8w)", fontsize=7)
    fig.legend(loc="upper right", fontsize=7, frameon=False)
    fig.suptitle("8w animals only: out-of-fold molecular score (logreg_l2, control vs 8w) vs own physiological change",
                 fontsize=9)
    fig.tight_layout(rect=(0, 0, 1, 0.98))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def forest(corr: pd.DataFrame, rcrit: float, path: Path):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig, axes = plt.subplots(1, 2, figsize=(9, 6), sharey=True)
    for ax, v in zip(axes, PHYS_PRIMARY):
        d = corr[(corr["physio"] == v) & (corr["score"] == "logit") & (corr["duration"] == "8w") & (corr["adjust"] == "sex_ranks")]
        labels = [f"{t} {a}" for t in TISSUES for a in ASSAYS]
        ypos = np.arange(len(labels))[::-1]
        for setname, c, off in [("within_8w", "#2a78d6", 0.15), ("within_control", "#8a8a85", -0.15)]:
            dd = d[d["set"] == setname].set_index(["tissue", "assay"])
            vals = [dd.loc[(t, a), "rho"] if (t, a) in dd.index else np.nan for t in TISSUES for a in ASSAYS]
            ax.scatter(vals, ypos + off, s=20, c=c, label=setname.replace("_", " "))
        ax.axvspan(-rcrit, rcrit, color="#e6e6e3", zorder=0, label=f"|rho| < {rcrit:.2f} (critical, 5F+5M, two-sided 0.05)")
        ax.axvline(0, color="#555", lw=0.8)
        ax.set_yticks(ypos, labels, fontsize=7)
        ax.set_xlim(-1, 1)
        ax.set_xlabel("sex-stratified rank correlation with " + v, fontsize=8)
        for sp in ["top", "right"]:
            ax.spines[sp].set_visible(False)
    axes[0].legend(fontsize=7, frameon=False, loc="lower left")
    fig.suptitle("Within-group score-physiology correlation (primary logit score, control-vs-8w model)", fontsize=9)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ main
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true", help="2 tissues, 2 fold seeds, 2000 permutations")
    args = ap.parse_args()
    global N_PERM
    tissues = TISSUES[:2] if args.quick else TISSUES
    n_seeds = 2 if args.quick else N_FOLD_SEEDS
    if args.quick:
        N_PERM = 2000
    out = PIPE / "results" / "15_time_course" / "4_physiology"
    out.mkdir(parents=True, exist_ok=True)

    ph = physiology()
    physiology_tables(ph, out)
    ph.to_csv(out / "physiology_per_animal.csv", index=False)

    pheno = io.load_pheno()
    source_trn = cli.resolve_source("TRNSCRPT", "auto")
    all_scores, infos = [], []
    for t in tissues:
        blocks, meta_t, notes = fus07.build_blocks([t], ASSAYS, pheno, source_trn, complete=True)
        for dur in ["8w", "4w"]:
            m = meta_t["group"].isin(["control", dur]).to_numpy()
            sub_meta = meta_t.loc[m].copy()
            sub_meta["sex_group"] = sub_meta["sex"].astype(str) + "/" + sub_meta["group"].astype(str)
            for a in ASSAYS:
                sc, info = oof_scores(blocks[a].loc[m], sub_meta, dur, n_seeds)
                sc["tissue"], sc["assay"], sc["duration"] = t, a, dur
                all_scores.append(sc)
                infos.append({"tissue": t, "assay": a, "duration": dur, "n_animals": int(sub_meta["pid"].nunique()),
                              "n_control": int((sub_meta["group"] == "control").sum()),
                              "n_trained": int((sub_meta["group"] == dur).sum()), **info})
                print(f"{t} {a} control vs {dur}: n={sub_meta['pid'].nunique()} AUROC(logit)={info['oof_auroc_logit']:.3f}")
    scores = pd.concat(all_scores, ignore_index=True)
    scores = scores.merge(ph.drop(columns=["bid", "sex", "group"]), on="pid", how="left")
    scores.to_csv(out / "oof_scores_with_physiology.csv", index=False)
    pd.DataFrame(infos).to_csv(out / "oof_model_auroc.csv", index=False)

    rng = np.random.default_rng(C.SEED)
    rows = []
    for (t, a, dur), s in scores.groupby(["tissue", "assay", "duration"]):
        sets = {"pooled": s, f"within_{dur}": s[s["group"] == dur], "within_control": s[s["group"] == "control"]}
        for setname, d in sets.items():
            for sc in ["logit", "proj", "dist"]:
                for v in PHYS_ALL:
                    dd = d.dropna(subset=[v, f"score_{sc}"])
                    if len(dd) < 5:
                        continue
                    for adj in ["none", "sex_ranks"]:
                        rho, p = perm_corr(dd[f"score_{sc}"].to_numpy(), dd[v].to_numpy(), dd["sex"].to_numpy(),
                                           adj == "sex_ranks", rng)
                        rows.append({"tissue": t, "assay": a, "duration": dur, "set": setname, "score": sc, "physio": v,
                                     "adjust": adj, "n": len(dd), "n_female": int((dd["sex"] == "female").sum()),
                                     "n_male": int((dd["sex"] == "male").sum()), "rho": rho, "p_perm": p})
    corr = pd.DataFrame(rows)
    # sign a training response would give: higher trained-ness <-> larger VO2max gain / larger fat loss (negative change)
    corr["expected_sign"] = corr["physio"].map({"vo2max_change": 1, "fat_change": -1, "lean_change": 1, "weight_change": -1})
    corr["rho_in_expected_direction"] = corr["rho"] * corr["expected_sign"]
    prim = ((corr["set"] == "within_8w") & (corr["duration"] == "8w") & (corr["score"] == "logit")
            & (corr["physio"].isin(PHYS_PRIMARY)) & (corr["adjust"] == "sex_ranks"))
    corr["primary"] = prim
    corr.loc[prim, "q_bh_primary"] = bh(corr.loc[prim, "p_perm"])
    corr.loc[prim, "p_bonf_primary"] = np.minimum(corr.loc[prim, "p_perm"] * prim.sum(), 1)
    # secondary families: BH within each (duration, set, score, adjust) family over tissue x assay x physio
    corr["q_bh_family"] = np.nan
    for _, idx in corr.groupby(["duration", "set", "score", "adjust"]).groups.items():
        corr.loc[idx, "q_bh_family"] = bh(corr.loc[idx, "p_perm"])
    corr.to_csv(out / "score_physiology_correlations.csv", index=False)

    n_typ = int(corr.loc[prim, "n"].median()) if prim.any() else 10
    power = pd.DataFrame([{"n": n, "critical_abs_rho_two_sided_0.05": critical_rho(n),
                           "true_abs_rho_for_80pct_power": rho_for_power(n)} for n in [5, 8, 9, 10, 11, 12, 18, 20]])
    strat = []
    for nf, nm in [(5, 5), (4, 5), (5, 6), (6, 6), (9, 9), (10, 10)]:
        sx = np.array(["female"] * nf + ["male"] * nm)
        r0 = np.random.default_rng(1)
        x = r0.normal(size=nf + nm)
        null = [abs(perm_corr(x, r0.normal(size=nf + nm), sx, True, r0, n_perm=1)[0]) for _ in range(20000)]
        strat.append({"n": nf + nm, "n_female": nf, "n_male": nm,
                      "critical_abs_rho_sex_stratified_0.05": float(np.quantile(null, 0.95))})
    power = power.merge(pd.DataFrame(strat), on="n", how="outer")
    power.to_csv(out / "power_reference.csv", index=False)
    rcrit = float(power.loc[power["n"] == 10, "critical_abs_rho_sex_stratified_0.05"].iloc[0])  # 5 F + 5 M
    figure(scores, corr, out / "score_vs_physiology_8w.png")
    forest(corr, rcrit, out / "within_group_rho_forest.png")
    print(corr[prim].sort_values("p_perm").head(10).to_string())
    print(f"primary tests: {int(prim.sum())}; median n {n_typ}")


if __name__ == "__main__":
    main()
