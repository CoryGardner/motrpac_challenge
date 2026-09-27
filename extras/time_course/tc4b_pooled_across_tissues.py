"""Part 4b (pre-specified 2026-09-26, before looking): one pooled within-group test of the physiology anchor.
Among 8w animals only, each animal's molecular 'trained-ness' = mean over all tissue x omic scores of its
within-sex, within-(tissue, omic) percentile rank of the out-of-fold logit (from tc4_physiology.py; control-vs-8w
logreg_l2, Task B folds). Correlate with the within-sex percentile rank of its VO2max change and of its fat change
(Spearman), permutation p by shuffling physiology within sex (10,000). Expected direction: + for VO2max, - for fat.
Negative control: the same in controls. Run from code/pipeline: python investigations/time_course/tc4b_pooled_across_tissues.py"""
from pathlib import Path
import numpy as np
import pandas as pd
from scipy import stats

D = Path("results/15_time_course/4_physiology")
s = pd.read_csv(D / "oof_scores_with_physiology.csv", dtype={"pid": str})
rows = []
for grp in ("8w", "control"):
    g = s[s["group"] == grp].copy()
    g["r"] = g.groupby(["tissue", "assay", "sex"])["score_logit"].rank(pct=True)
    a = g.groupby(["pid", "sex"]).agg(score=("r", "mean"), n_scores=("r", "size"),
                                      vo2=("vo2max_change", "first"), fat=("fat_change", "first")).reset_index()
    for var in ("vo2", "fat"):
        b = a.dropna(subset=[var]).copy()
        b["pr"] = b.groupby("sex")[var].rank(pct=True)
        rho = stats.spearmanr(b["score"], b["pr"]).statistic
        rng = np.random.default_rng(20260925)
        null = []
        for _ in range(10000):
            perm = b.groupby("sex")["pr"].transform(lambda x: rng.permutation(x.to_numpy()))
            null.append(stats.spearmanr(b["score"], perm).statistic)
        null = np.array(null)
        rows.append({"group": grp, "physiology": {"vo2": "VO2max change", "fat": "fat change (pts)"}[var],
                     "n_animals": len(b), "n_female": int((b.sex == "female").sum()), "n_male": int((b.sex == "male").sum()),
                     "median_scores_per_animal": float(b["n_scores"].median()), "spearman_rho": rho,
                     "p_perm_two_sided": float((np.abs(null) >= abs(rho)).mean()),
                     "critical_abs_rho_0.05": float(np.quantile(np.abs(null), 0.95))})
out = pd.DataFrame(rows)
out.to_csv(D / "pooled_across_tissues.csv", index=False)
print(out.round(3).to_string(index=False))
