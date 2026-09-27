# %% [markdown]
# ## 4. Tissue-fingerprint baselines (transcriptome, 19 tissues)
#
# **Question.** How well do *simple, tuned* classifiers tell the 19 tissues apart from RNA-seq, when
# every sample of a test animal is kept out of training?
#
# **Why it matters.** Every later claim — compact panels (section 5), certified panel sizes, shift
# tests, fusion — has to beat this reference on the same folds. If the simple models are already near
# the ceiling, "tissue identity from the transcriptome" is not a result by itself; the content is in
# how *small* and how *robust* the fingerprint can be.
#
# **What to look for.** (1) The leakage guard: a split over *samples* puts the same animal on both
# sides, and the pipeline's assertion refuses it; the animal-grouped split passes. (2) The balanced
# accuracy of the three baselines, recomputed here, next to the published values. (3) Which tissue
# pairs remain confused even with all genes.
#
# Design (from the Makefile target `baselines`, with the published `MODELS=centroid,logreg_l2,rf`):
# raw counts stacked across tissues → log2 CPM, 5 animal-grouped outer folds stratified on tissue,
# variance prefilter to 5,000 genes inside each fold, a small grid search inside each training fold.
# The L1 logistic model is in the library's default list but was **never run** in the published
# pipeline; this notebook reproduces the three models that were.

# %%
section("4 baselines")


# from src/tfp/cli.py::resolve_source (verbatim) — cli.py is not pasted into this notebook
def resolve_source(assay: str, source: str) -> str:
    if source != "auto":
        return source
    if assay == "TRNSCRPT" and any((C.COUNTS_DIR).glob(f"{assay}__*.csv")):
        return "counts"
    return "norm"


# from scripts/04_fingerprint_baselines.py::main (trimmed: argparse defaults inlined —
# --assay TRNSCRPT, --source auto, no --tissues, no --complete-features, no --drop-incomplete-samples)
S04_SOURCE = resolve_source("TRNSCRPT", "auto")
pheno = cached("pheno", io.load_pheno)
S04_OM = cached("TRNSCRPT_counts", lambda: io.stack_tissues(
    "TRNSCRPT", tissues=None, source=resolve_source("TRNSCRPT", "auto"), join="inner", pheno=pheno,
    complete=False, drop_incomplete_samples=None))
om = S04_OM
print("source:", S04_SOURCE)
print(om)
# the stacking note (genes kept by the pre-split filter). Picked by content: the script's report used
# om.notes[-1], which for counts is the long per-gene "expressed in n tissues" list, not this note.
print(next(n for n in om.notes if n.startswith("stack ")))
record("s04.n_samples", om.n_samples, "04")
record("s04.n_features", om.n_features, "04")
record("s04.n_animals", om.n_animals, "04")

# %% [markdown]
# ### 4a. The leakage guard, fired on purpose
#
# Each rat contributes up to 19 tissue samples. A plain `KFold` over *samples* scatters one animal's
# tissues across train and test, so a model could learn that animal's individual signature (sex,
# library batch, genotype) and get credit for it on the test side. `splits.assert_no_group_leak` checks
# every split the pipeline makes. Below it is run first on a deliberately wrong split, then on the
# grouped split the pipeline actually uses.

# %%
from sklearn.model_selection import KFold

S04_N_SPLITS = 5          # --n-splits default
# (1) the wrong split: 5-fold over samples, ignoring which animal each sample came from
bad_tr, bad_te = next(KFold(n_splits=S04_N_SPLITS, shuffle=True, random_state=C.SEED).split(om.X))
try:
    splits.assert_no_group_leak(om.meta, bad_tr, bad_te)          # group = "pid" by default
    print("sample-level split: no leak detected (unexpected)")
except AssertionError as e:
    print("sample-level split REJECTED ->", e)
    n_shared = len(set(om.meta.iloc[bad_tr]["pid"]) & set(om.meta.iloc[bad_te]["pid"]))
    print(f"  animals on both sides: {n_shared} of {om.n_animals}")
    record("s04.leak_shared_animals_samplekfold", n_shared, "04", "this run; plain KFold(5, shuffle, SEED)")

# (2) the pipeline's split: stratified on tissue, grouped by animal (pid). grouped_kfold itself calls the
# assertion on every fold; here it is called once more, explicitly, so the pass is visible.
# These are exactly the folds models.evaluate_cv makes (same meta, label, n_splits, seed) — section 13 reuses them.
S04_FOLDS = list(grouped_kfold(om.meta, "tissue", S04_N_SPLITS, C.SEED))
# the fold of each vial, written so notebook 02 (section 13) can assert it rebuilt exactly these folds
(OUT / "s04").mkdir(parents=True, exist_ok=True)
pd.DataFrame({"viallabel": S04_OM.meta.index[np.concatenate([te for _, te in S04_FOLDS])],
              "fold": np.concatenate([[k] * len(te) for k, (_, te) in enumerate(S04_FOLDS)])}).to_csv(OUT / "s04" / "folds.csv", index=False)
rows = []
for f, (tr, te) in enumerate(S04_FOLDS):
    splits.assert_no_group_leak(om.meta, tr, te)
    rows.append(splits.describe_split(om.meta, tr, te) | {"fold": f})
print("animal-grouped split: all", len(S04_FOLDS), "folds pass the assertion")
display(pd.DataFrame(rows).set_index("fold"))

# %% [markdown]
# **What this shows.** The assertion is live: a sample-level split is refused with the list of animals
# that sit on both sides, and the animal-grouped folds pass. The table gives the animals and samples per
# fold. **What it does not show.** Grouping by animal removes animal-level leakage only. It does not
# remove *batch* leakage: each tissue was extracted and sequenced as one batch (section 2), so "tissue"
# and "library batch" cannot be separated by any split of these data.

# %%
# from scripts/04_fingerprint_baselines.py::main (trimmed: dropped --quick and --time-one-fold branches,
# report.add_section, writing per_fold/predictions CSVs to results/).
# Runtime: the three models took 150-200 s here on the Spark (random forest is about half of it). To keep
# the default run short, RECOMPUTE=False fits centroid and logreg_l2 live and LOADS the random-forest
# folds from results/; RECOMPUTE=True fits all three.
S04_MODELS = ["centroid", "logreg_l2", "rf"]      # Makefile MODELS (published run); logreg_l1 was never run
S04_LIVE = S04_MODELS if RECOMPUTE else ["centroid", "logreg_l2"]
S04_PREFILTER = 5000                              # --prefilter default
s04_res = models.evaluate_cv(om, label="tissue", kinds=S04_LIVE, n_splits=S04_N_SPLITS,
                             prefilter=S04_PREFILTER, quick=False, seed=C.SEED)
s04_per_fold = s04_res.per_fold.assign(source="recomputed")
s04_conf = dict(s04_res.confusion)
for kind in [k for k in S04_MODELS if k not in S04_LIVE]:          # loaded, not recomputed
    pf = pd.read_csv(res("04_baselines", "TRNSCRPT", "per_fold.csv"))
    s04_per_fold = pd.concat([s04_per_fold, pf[pf["model"] == kind].assign(source="loaded")], ignore_index=True)
    s04_conf[kind] = pd.read_csv(res("04_baselines", "TRNSCRPT", f"confusion_{kind}.csv"), index_col=0)
summ = CVResult(s04_per_fold.drop(columns="source"), s04_res.predictions, s04_conf).summary()
fit_time = s04_per_fold.groupby("model")["fit_seconds"].agg(fit_seconds_mean="mean", fit_seconds_total="sum").reset_index()
summ = summ.merge(fit_time, on="model").merge(s04_per_fold.groupby("model")["source"].first().reset_index(), on="model")
(OUT / "s04").mkdir(parents=True, exist_ok=True)
summ.to_csv(OUT / "s04" / "summary.csv", index=False)

# recomputed vs published, side by side
pub = pd.read_csv(res("04_baselines", "TRNSCRPT", "summary.csv"))
cmp = summ[["model", "source", "balanced_accuracy_mean", "balanced_accuracy_std", "macro_f1_mean", "log_loss_mean"]].merge(
    pub[["model", "balanced_accuracy_mean", "balanced_accuracy_std"]], on="model", suffixes=("", "_published"))
cmp["diff"] = cmp["balanced_accuracy_mean"] - cmp["balanced_accuracy_mean_published"]
display(cmp.round(4))
for _, r in cmp.iterrows():
    record(f"s04.bal_acc_{r['model']}", r["balanced_accuracy_mean"], "04", r["source"])
best = summ.sort_values("balanced_accuracy_mean", ascending=False).iloc[0]
print(f"best baseline: {best['model']}  balanced accuracy {best['balanced_accuracy_mean']:.3f} ± {best['balanced_accuracy_std']:.3f} (sd over folds)")
print("per-fold balanced accuracy:")
display(s04_per_fold.pivot(index="fold", columns="model", values="balanced_accuracy").round(3))

# %%
# from scripts/04_fingerprint_baselines.py::confusable_pairs (verbatim)
def confusable_pairs(conf: pd.DataFrame, top: int = 10) -> pd.DataFrame:
    rows = []
    for a in conf.index:
        for b in conf.columns:
            if a != b and conf.loc[a, b] > 0:
                rows.append({"true": a, "predicted": b, "count": int(conf.loc[a, b]),
                             "frac_of_true": float(conf.loc[a, b] / max(conf.loc[a].sum(), 1))})
    cols = ["true", "predicted", "count", "frac_of_true"]
    df = pd.DataFrame(rows, columns=cols)  # keeps the header even when nothing is confused
    return df.sort_values("count", ascending=False).head(top)


# remaining confusions with all genes, per model (pooled over the 5 test folds), vs the published lists
for kind in S04_MODELS:
    conf = s04_conf[kind]
    live = confusable_pairs(conf).reset_index(drop=True)
    pubp = pd.read_csv(res("04_baselines", "TRNSCRPT", f"confusable_pairs_{kind}.csv"))  # top 10 pairs
    key = ["true", "predicted", "count"]
    same = live[key].sort_values(key).reset_index(drop=True).equals(pubp[key].sort_values(key).reset_index(drop=True))
    print(f"{kind} ({'recomputed' if kind in S04_LIVE else 'loaded'}): {int(live['count'].sum())} misclassified samples of {int(conf.to_numpy().sum())}; "
          f"confused pairs identical to published: {same}")
    if kind == "logreg_l2":
        display(live)
        record("s04.n_errors_logreg_l2", int(live["count"].sum()), "04")
show_png(plots.confusion_heatmap(s04_conf["logreg_l2"],
                                 "TRNSCRPT logreg_l2: confusion (row-normalized)", OUT / "s04" / "confusion_logreg_l2.png"),
         width=600)

# %% [markdown]
# **What this shows.** All three tuned baselines separate the 19 tissues almost perfectly under
# animal-grouped cross-validation, and the recomputed numbers match the published ones (the `diff`
# column; the `source` column says which rows were fitted in this run and which were loaded). The fold-to-fold spread is the same size as the gaps between models, so the three are
# effectively tied; the L2 logistic regression is the reference later sections must beat. The few
# errors that remain are mostly between tissues that are anatomically adjacent or similar (listed above).
#
# **What it does not show.** High accuracy here is not evidence of a biological "fingerprint" in any
# strong sense: tissue is confounded with processing batch, and the test animals come from the same
# cohorts, lab and pipeline as the training animals. It says nothing about transfer to another lab,
# another library chemistry, or trained animals held out as a group (sections 8–10). The gene filter
# (≥ 10 counts in ≥ 20 % of one tissue's samples) is applied to the full matrix before splitting — a
# known deviation, reproduced here as the pipeline ran it. That filter sees the tissue labels of the
# test samples too, so any bias it adds is optimistic; the variance prefilter and all tuning are inside the folds.
