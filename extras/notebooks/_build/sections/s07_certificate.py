# %% [markdown]
# ## 7. The certified panel size (Learn-then-Test with Clopper–Pearson)
#
# **Question.** Section 5 found the smallest panel whose *average* cross-validated balanced accuracy
# reaches 1 − α. A certificate asks for more: a panel size k such that, with probability at least
# 1 − δ over the draw of calibration animals, the panel's error on new animals is at most α. The
# pipeline certifies k by Learn-then-Test: for each k on a fixed grid (10, 15, 20, 30, 50, 100), fit
# the k-gene panel on the fit animals, count its errors on the calibration animals, compute the
# one-sided Clopper–Pearson upper bound on the error rate, and test from the largest k downward,
# stopping at the first k whose bound exceeds α. The primary certificate (`make conformal`, first
# command) uses α = δ = 0.10, 22 calibration animals and **one vial per animal** as the loss unit;
# the fit/calibration split is repeated 20 times in each of 5 outer folds (100 splits), so the
# certified k comes with a distribution.
#
# **Why it matters.** "20 genes identify the tissue" is a point estimate; the certificate says how
# many genes one can *promise* with this many animals — and whether the study is large enough to
# promise anything.
#
# **What to look for.** (1) The zero-error sizing table: how many calibration animals a certificate
# needs when the panel makes no calibration error at all. (2) How often each k is certified over the
# 100 splits, and how often no k is. (3) Validity: at the certified k, is the test error ≤ α in at
# least 1 − δ of the splits? (4) Why the certified k is larger than the point-estimate k.

# %% [markdown]
# ### 7a. Zero-error sizing: how many calibration units does a certificate need?
#
# With zero errors in n calibration units, the one-sided Clopper–Pearson upper bound is
# 1 − δ^(1/n). It falls to α only when n ≥ ln δ / ln(1 − α). The cell computes that bound for the
# pipeline's four (α, δ) pairs and checks it against the pipeline's `sizing_table.csv` (post-fix and
# pre-fix copies), and against the library's `clopper_pearson_upper`.

# %%
section("7 certificate")
# from scripts/06_conformal_certify.py (verbatim constant and the sizing block of main)
SIZING_PAIRS = [(0.05, 0.05), (0.10, 0.05), (0.10, 0.10), (0.20, 0.10)]
s07_sizing = pd.DataFrame([{"alpha": a, "delta": d, "n_zero_error_needed": int(np.ceil(np.log(d) / np.log(1 - a)))}
                           for a, d in SIZING_PAIRS])
# check with the library bound: the smallest n with UCB(0 errors, n, δ) ≤ α
s07_sizing["n_by_cp_search"] = [next(n for n in range(1, 500) if cp.clopper_pearson_upper(0, n, d) <= a) for a, d in SIZING_PAIRS]
S07_PRE = REPO / "backup" / "pipeline_history" / "results_pre_quantile_fix_2026-09-25" / "06_conformal"
pub = pd.read_csv(res("06_conformal", "TRNSCRPT", "sizing_table.csv"))
pre = pd.read_csv(S07_PRE / "TRNSCRPT" / "sizing_table.csv")
s07_sizing["published_post_fix"] = pub["n_zero_error_needed"].to_numpy()
s07_sizing["published_pre_fix"] = pre["n_zero_error_needed"].to_numpy()
display(s07_sizing)
for (a, d), n in zip(SIZING_PAIRS, s07_sizing["n_zero_error_needed"]):
    record(f"s07.sizing_a{a:.2f}_d{d:.2f}", n, "07", "ceil(ln δ / ln(1 − α))")

# the bound at the primary design (n = 22 calibration animals, δ = 0.10) for 0, 1, 2 errors
s07_n_cal, s07_alpha, s07_delta = 22, 0.10, 0.10
s07_ucb = pd.DataFrame({"n_err": [0, 1, 2]})
s07_ucb["point_error"] = s07_ucb["n_err"] / s07_n_cal
s07_ucb["cp_upper_bound"] = [cp.clopper_pearson_upper(e, s07_n_cal, s07_delta) for e in s07_ucb["n_err"]]
s07_ucb["passes_alpha"] = s07_ucb["cp_upper_bound"] <= s07_alpha
print(f"Clopper–Pearson (1 − δ) upper bound with n = {s07_n_cal} calibration animals, δ = {s07_delta}:")
display(s07_ucb.round(4))
record("s07.ucb_0err_n22", s07_ucb.loc[0, "cp_upper_bound"], "07", "CP upper bound, 0 errors of 22, δ=0.10")
record("s07.ucb_1err_n22", s07_ucb.loc[1, "cp_upper_bound"], "07", "CP upper bound, 1 error of 22, δ=0.10")

# %% [markdown]
# **What this shows.** The table's first column is the whole design problem: to certify a tissue call
# at (0.05, 0.05) one needs more zero-error calibration animals than the study has per tissue, and at
# the primary (0.10, 0.10) the 22 calibration animals are *exactly* the minimum — the bound with zero
# errors sits just under α, and a single wrong calibration animal pushes it far above. So the primary
# certificate is a zero-error certificate: a panel size is certified only if the panel classifies the
# one random vial of every one of the 22 calibration animals correctly.
#
# **What it does not show.** Nothing here depends on the data; it is arithmetic on the binomial. It
# also does not depend on the conformal quantile (no quantile is involved).

# %% [markdown]
# ### 7b. The certificate over 100 fit/calibration splits
#
# With `RECOMPUTE = False` the 100-split certificate is loaded from the pipeline's
# `results/06_conformal/TRNSCRPT/` (it fits 100 + 20 panel families, a few minutes). With
# `RECOMPUTE = True` the cell below re-runs it (copied script logic) and writes under `OUT/s07/`.
# Either way, the distribution, the validity check and the (α, δ) grid are then computed *in this
# notebook* from the per-split tables and compared with the pipeline's summary files.

# %%
# from scripts/06_conformal_certify.py::panel_family (verbatim)
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


def panel_family(X_fit, y_fit, grid, prefilter, model):
    """Preprocessing (impute → variance prefilter → scale → round-robin scores) fit once on the fit
    animals; then one k-panel classifier per k. Same model as the quick pipeline of phase 05
    (logreg_l2: C = 0.1). Returns predict(k, X)."""
    if model != "logreg_l2":
        raise ValueError("the certificate is defined for --model logreg_l2 (round-robin panel + L2 logistic regression)")
    imp = SimpleImputer(strategy="median").fit(X_fit)
    pre = models.VarianceTopK(prefilter).fit(imp.transform(X_fit))
    sc = StandardScaler().fit(pre.transform(imp.transform(X_fit)))
    Xs = sc.transform(pre.transform(imp.transform(X_fit)))
    sel = models.RoundRobinSelector(k=max(grid)).fit(Xs, y_fit)
    order = sel.order_
    fam = {}
    for k in grid:
        idx = np.sort(order[:k])
        fam[k] = (idx, LogisticRegression(solver="lbfgs", C=0.1, max_iter=3000).fit(Xs[:, idx], y_fit))

    def predict(k, X):
        Z = sc.transform(pre.transform(imp.transform(X)))
        idx, clf = fam[k]
        return clf.predict(Z[:, idx])
    return predict


# from scripts/06_conformal_certify.py::certify_losses (verbatim)
def certify_losses(predict, X, y, g, cal_idx, cal_opa, te, grid, alpha, delta, one_per_animal):
    """Certified k under the chosen loss and under the per-animal any-wrong loss; test error per k."""
    errors, errors_animal, test_err = {}, {}, {}
    m_opa = np.isin(cal_idx, cal_opa)
    for k in grid:
        wrong = predict(k, X[cal_idx]) != y[cal_idx]
        errors[k] = (int(wrong[m_opa].sum()), int(m_opa.sum())) if one_per_animal else (int(wrong.sum()), len(cal_idx))
        any_wrong = pd.Series(wrong).groupby(g[cal_idx]).any()
        errors_animal[k] = (int(any_wrong.sum()), int(len(any_wrong)))
        test_err[k] = float((predict(k, X[te]) != y[te]).mean()) if len(te) else np.nan
    k1, t1 = cp.certify_panel_size(errors, alpha, delta)
    k2, t2 = cp.certify_panel_size(errors_animal, alpha, delta)
    t1["loss"] = "one_vial_per_animal" if one_per_animal else "pooled_vials"
    t2["loss"] = "animal_any_tissue_wrong"
    return k1, k2, pd.concat([t1, t2], ignore_index=True), test_err, errors

# %%
# the flags of `make conformal` (first command)
s07_args = SimpleNamespace(label="tissue", model="logreg_l2", n_splits=5, prefilter=5000, seed=C.SEED, cal_frac=0.3,
                           cal_animals=22, n_repeats=20, one_per_animal=True, alpha=0.1, delta=0.1,
                           grid=[10, 15, 20, 30, 50, 100], k_of_interest=20)
LOSS = "one_vial_per_animal"
if RECOMPUTE:
    # from scripts/06_conformal_certify.py::main — the repeat loop, certificate leg (b), and the deployment certificate
    # (trimmed: conformal-set leg (a) is section 6; `args.` → `a.`; `om` → `om_trn`)
    a = s07_args
    om_trn = cached("TRNSCRPT_counts", lambda: io.stack_tissues(
        "TRNSCRPT", tissues=None, source=resolve_source("TRNSCRPT", "auto"), join="inner", pheno=cached("pheno", io.load_pheno),
        complete=False, drop_incomplete_samples=None))
    X = om_trn.X.to_numpy(dtype=float)
    y = om_trn.meta[a.label].astype(str).to_numpy()
    g = om_trn.groups()
    grid = [k for k in a.grid if k <= X.shape[1]]
    cert_rows, valid_rows = [], []
    for fold, (tr, te) in enumerate(grouped_kfold(om_trn.meta, a.label, a.n_splits, a.seed)):
        for rep in range(a.n_repeats):
            seed_r = a.seed + 1000 * fold + rep
            fit_idx, cal_idx = split_train(om_trn.meta, tr, a.cal_frac, a.cal_animals, seed_r)
            cal_opa = cp.one_per_group(cal_idx, g, seed_r)
            predict = panel_family(X[fit_idx], y[fit_idx], grid, a.prefilter, a.model)
            k1, k2, table, test_err, errs = certify_losses(predict, X, y, g, cal_idx, cal_opa, te, grid, a.alpha, a.delta,
                                                           a.one_per_animal)
            table["fold"], table["repeat"] = fold, rep
            table["test_err"] = table["k"].map(test_err)
            cert_rows.append(table)
            valid_rows.append({"fold": fold, "repeat": rep, "n_fit_animals": int(om_trn.meta.iloc[fit_idx]["pid"].nunique()),
                               "n_cal_animals": int(om_trn.meta.iloc[cal_idx]["pid"].nunique()),
                               "n_cal_loss": int(len(cal_opa)) if a.one_per_animal else int(len(cal_idx)),
                               "certified_k": k1, "certified_k_animal_any_wrong": k2,
                               "test_err_at_certified_k": test_err.get(k1, np.nan) if k1 else np.nan,
                               "test_err_le_alpha": (test_err.get(k1, 1.0) <= a.alpha) if k1 else None})
        print(f"  fold {fold}: certified k over {a.n_repeats} repeats: {[v['certified_k'] for v in valid_rows if v['fold'] == fold]}")
    s07_cert = pd.concat(cert_rows, ignore_index=True)
    s07_valid = pd.DataFrame(valid_rows)
    all_idx = np.arange(len(y))
    final_rows = []
    for rep in range(a.n_repeats):
        fit_idx, cal_idx = split_train(om_trn.meta, all_idx, a.cal_frac, a.cal_animals, a.seed + 999 + rep)
        cal_opa = cp.one_per_group(cal_idx, g, a.seed + 999 + rep)
        predict = panel_family(X[fit_idx], y[fit_idx], grid, a.prefilter, a.model)
        k1, k2, table, _, _ = certify_losses(predict, X, y, g, cal_idx, cal_opa, np.array([], dtype=int), grid,
                                             a.alpha, a.delta, a.one_per_animal)
        final_rows.append({"repeat": rep, "n_fit_animals": int(om_trn.meta.iloc[fit_idx]["pid"].nunique()),
                           "n_cal_animals": int(om_trn.meta.iloc[cal_idx]["pid"].nunique()), "certified_k": k1,
                           "certified_k_animal_any_wrong": k2})
    s07_final = pd.DataFrame(final_rows)
    del a
    (OUT / "s07").mkdir(parents=True, exist_ok=True)
    s07_cert.to_csv(OUT / "s07" / "certificate_by_fold.csv", index=False)
    s07_valid.to_csv(OUT / "s07" / "certificate_validity.csv", index=False)
    s07_final.to_csv(OUT / "s07" / "certificate_final_repeats.csv", index=False)
else:
    s07_cert = pd.read_csv(res("06_conformal", "TRNSCRPT", "certificate_by_fold.csv"))
    s07_valid = pd.read_csv(res("06_conformal", "TRNSCRPT", "certificate_validity.csv"))
    s07_final = pd.read_csv(res("06_conformal", "TRNSCRPT", "certificate_final_repeats.csv"))
print(f"{len(s07_valid)} fit/calibration splits; {s07_valid['n_cal_loss'].iloc[0]} calibration units per split "
      f"({s07_valid['n_fit_animals'].iloc[0]} fit / {s07_valid['n_cal_animals'].iloc[0]} calibration animals)")

# %%
# from scripts/06_conformal_certify.py::main (verbatim: distribution of the certified k over folds × repeats)
grid = s07_args.grid
valid = s07_valid
dist = []
for k in grid:
    dist.append({"k": k, "frac_certified": float(valid["certified_k"].apply(lambda v: v is not None and not pd.isna(v) and v <= k).mean()),
                 "frac_certified_animal_any_wrong": float(valid["certified_k_animal_any_wrong"].apply(lambda v: v is not None and not pd.isna(v) and v <= k).mean())})
dist = pd.DataFrame(dist)
k_half = dist.loc[dist["frac_certified"] >= 0.5, "k"].min() if (dist["frac_certified"] >= 0.5).any() else None
ks = valid["certified_k"].dropna()
n_none = int(valid["certified_k"].isna().sum())
s07_dist = dist.merge(pd.read_csv(res("06_conformal", "TRNSCRPT", "certificate_distribution.csv")), on="k", suffixes=("", "_published"))
print("fraction of the splits whose certified k is ≤ k (cumulative), this notebook vs the pipeline's file:")
display(s07_dist)
print("certified k over the splits (value: number of splits):", ks.astype(int).value_counts().sort_index().to_dict(),
      f"| none certified: {n_none}")
print(f"certified k: min {int(ks.min())}, median {ks.median()}, max {int(ks.max())}; smallest k certified in ≥ half the splits: {k_half}")
s07_validity = float(valid["test_err_le_alpha"].dropna().astype(bool).mean())
print(f"validity: fraction of certified splits with test error ≤ α at the certified k = {s07_validity:.2f} "
      f"(target ≥ 1 − δ = {1 - s07_args.delta:.2f}); worst test error at a certified k = {valid['test_err_at_certified_k'].max():.3f}")
print("deployment certificate on all animals (20 repeats), certified k:",
      s07_final["certified_k"].value_counts(dropna=False).sort_index().to_dict())
record("s07.frac_certified_k20", float(dist.loc[dist.k == 20, "frac_certified"].iloc[0]), "07", "fraction of 100 splits certifying k ≤ 20")
record("s07.frac_certified_k15", float(dist.loc[dist.k == 15, "frac_certified"].iloc[0]), "07", "fraction of 100 splits certifying k ≤ 15")
record("s07.n_none_certified", n_none, "07", "splits certifying no k")
record("s07.k_half", k_half, "07", "smallest k certified in ≥ half the splits")
record("s07.validity", s07_validity, "07", "fraction of certified splits with test error ≤ α")
record("s07.frac_any_animal_any_wrong", float(dist["frac_certified_animal_any_wrong"].max()), "07", "per-animal any-wrong loss: fraction certifying any k")
del valid, dist, grid

# %%
# (α, δ) grid on the same calibration errors (from scripts/06_conformal_certify.py::main, verbatim logic;
# the per-split error counts `err_store` are rebuilt from certificate_by_fold.csv, loss = one vial per animal)
cb = s07_cert[s07_cert["loss"] == LOSS]
err_store = [{int(r.k): (int(r.n_err), int(r.n_cal)) for r in grp.itertuples()} for _, grp in cb.groupby(["fold", "repeat"], sort=True)]
ad_grid = [(0.05, 0.05), (0.10, 0.05), (0.10, 0.10), (0.20, 0.05), (0.20, 0.10), (0.05, 0.10)]
ad_rows = []
for a_, d_ in ad_grid:
    ks_ad = [cp.certify_panel_size(e, a_, d_)[0] for e in err_store]
    ks_ok = [k for k in ks_ad if k is not None]
    ad_rows.append({"alpha": a_, "delta": d_, "n_splits": len(ks_ad),
                    f"frac_certified_k{s07_args.k_of_interest}": float(np.mean([(k is not None and k <= s07_args.k_of_interest) for k in ks_ad])),
                    "frac_any_certified": float(np.mean([k is not None for k in ks_ad])),
                    "median_certified_k": float(np.median(ks_ok)) if ks_ok else np.nan})
s07_ad = pd.DataFrame(ad_rows).merge(pd.read_csv(res("06_conformal", "TRNSCRPT", "certificate_alpha_delta_grid.csv")),
                                     on=["alpha", "delta"], suffixes=("", "_published"))
display(s07_ad)
record("s07.ad_a20_d10_frac_k20", float(s07_ad.loc[(s07_ad.alpha == 0.2) & (s07_ad.delta == 0.1), "frac_certified_k20"].iloc[0]), "07",
       "(α, δ) = (0.20, 0.10): fraction certifying k ≤ 20")
del cb

# %% [markdown]
# **What this shows.** The certified k is a random variable, not a number: across the 100 splits it
# lands most often on the smallest certifiable sizes, but it is spread over the whole grid and no k
# is certified at all in a sizeable share of splits (the "none certified" count). Where it certifies,
# it is valid — the test error at the certified k stays at or below α in the validity line above.
# The per-animal "any tissue of this animal wrong" loss never certifies anything (its column in the
# distribution table): an animal counts as wrong if any of its ~18 vials is wrong, and at every k
# enough calibration animals have one wrong vial to break a zero-error certificate. Loosening the target to
# α = 0.20 (the grid) makes the certificate routine; tightening δ or α to 0.05 makes it impossible —
# consistent with the sizing table.
#
# **What it does not show.** 100 splits of one dataset are not 100 independent studies: the splits
# share animals, so the distribution is narrower than a new study's would be. The deployment
# certificate re-uses all animals and has no test set, so it has no validity check of its own.

# %% [markdown]
# ### 7c. Why the certified k exceeds the point-estimate k
#
# **Question.** Section 5's panel curve reaches mean balanced accuracy ≥ 1 − α at a smaller k than
# the certificate usually certifies. Where does the gap come from? The table below is derived from the
# per-split certificate table: for every k, the mean test error over the splits (what a point estimate
# sees), the fraction of splits in which the panel made **zero** errors on the 22 calibration animals,
# the zero-error probability that the mean test error predicts if the 22 calibration vials were
# independent draws, `(1 − err)^22`, and the fraction of splits certifying *this very* k (fixed-sequence
# testing also needs every larger k to pass).

# %%
# point-estimate k: from the phase-05 panel curve (verbatim logic of scripts/06_conformal_certify.py::main)
curve = pd.read_csv(res("05_panels", "TRNSCRPT", "panel_curve.csv")).groupby("k")["balanced_accuracy"].mean()
s07_point_k = int(curve[curve >= 1 - s07_args.alpha].index.min()) if (curve >= 1 - s07_args.alpha).any() else None
print(f"point-estimate k (smallest k with mean 5-fold balanced accuracy ≥ {1 - s07_args.alpha:.2f}): {s07_point_k}")
record("s07.point_k", s07_point_k, "07", "phase-05 curve, mean balanced accuracy ≥ 0.90")

cb = s07_cert[s07_cert["loss"] == LOSS]
n = int(cb["n_cal"].iloc[0])
s07_why = cb.groupby("k").agg(mean_test_err=("test_err", "mean"), mean_cal_err=("err_hat", "mean"),
                              frac_zero_cal_errors=("n_err", lambda v: float((v == 0).mean())),
                              mean_ucb=("ucb", "mean"), frac_k_passes_alone=("ucb", lambda v: float((v <= s07_args.alpha).mean()))).reset_index()
s07_why["pred_zero_err_prob"] = (1 - s07_why["mean_test_err"]) ** n
s07_why["frac_certified_this_k"] = s07_why["k"].map(s07_valid["certified_k"].value_counts() / len(s07_valid)).fillna(0.0)
s07_why["frac_certified_le_k"] = s07_why["k"].map(lambda k: float((s07_valid["certified_k"] <= k).mean()))
s07_why["curve_bal_acc"] = s07_why["k"].map(curve)
display(s07_why.round(3))
for k in (15, 20):
    r = s07_why.set_index("k").loc[k]
    record(f"s07.frac_zero_err_k{k}", r["frac_zero_cal_errors"], "07", f"fraction of splits with 0 calibration errors at k={k}")
    record(f"s07.mean_test_err_k{k}", r["mean_test_err"], "07", f"mean test error at k={k} over 100 splits")
del cb, n, curve

# %% [markdown]
# **What this shows — the derivation.** Read the table from left to right.
#
# 1. *The point estimate is an average.* The point-estimate k is the first grid value whose mean
#    balanced accuracy (`curve_bal_acc`) clears 1 − α; its mean test error is already well below α.
#    On average, the panel is good enough there — and that is all the phase-05 curve says.
# 2. *The certificate is a zero-error test at this design.* From 7a, with 22 calibration units and
#    δ = 0.10, zero errors give a bound just below α and one error gives a bound far above it. So a
#    k passes only in the splits where the panel got all 22 calibration vials right
#    (`frac_k_passes_alone` equals `frac_zero_cal_errors`).
# 3. *Zero errors in 22 draws is rare unless the error rate is far below α.* If the error rate is p,
#    the chance of 22 correct calls is about (1 − p)^22; the `pred_zero_err_prob` column shows this
#    tracks the observed zero-error fraction. Even at the point-estimate k, whose error is well below
#    α, 22 correct calls in a row happen in only a minority of splits; a k is certified in about half
#    of the splits only once its error is several times smaller than α.
# 4. *Fixed-sequence testing compounds it.* Testing starts at the largest k and stops at the first
#    failure, so a k is certified only if it *and every larger k* made zero calibration errors — one
#    unlucky vial at k = 30 ends the certificate at 50, whatever smaller panels do. That is why the
#    cumulative certified fraction at a k (`frac_certified_le_k`) is at or below that k's own pass rate.
#
# Together: the certified k is the smallest k whose error is small enough that 22 consecutive correct
# calls are likely — a stricter target than "average error ≤ α", so it sits above the point-estimate
# k. The gap is the price of a guarantee over an estimate, and with this study's animal count it
# shrinks only if α or δ are relaxed or more calibration animals are added.
#
# **What it does not show.** The calibration vials are not independent draws of one error rate (the
# error rate differs by tissue, and one vial per animal is a random tissue), so `(1 − err)^22` is an
# approximation used only to explain the pattern; the certificate itself uses the exact counts.

# %% [markdown]
# ### 7d. Did the conformal-quantile fix move the certificate?
#
# **What to look for.** It should not: the certificate never calls the conformal quantile. It counts
# calibration errors and bounds them with Clopper–Pearson. The cell compares every certificate file
# of both certificate runs (primary; and the secondary pooled-vial run at α = δ = 0.05) between the
# pre-fix backup and the post-fix pipeline results.

# %%
rows = []
for run in ("TRNSCRPT", "TRNSCRPT_pooled_a05"):
    for f in ("certificate_by_fold.csv", "certificate_validity.csv", "certificate_distribution.csv", "certificate_final.csv",
              "certificate_final_repeats.csv", "certificate_alpha_delta_grid.csv", "sizing_table.csv"):
        a_, b_ = pd.read_csv(S07_PRE / run / f), pd.read_csv(res("06_conformal", run, f))
        same_shape = a_.shape == b_.shape and list(a_.columns) == list(b_.columns)
        ints = [c for c in a_.columns if c in ("k", "n_cal", "n_err", "certified", "certified_k", "certified_k_animal_any_wrong",
                                               "n_zero_error_needed", "frac_certified", "frac_any_certified", "median_certified_k",
                                               "frac_certified_k20", "test_err_le_alpha", "loss")]
        counts_equal = same_shape and all(a_[c].astype(str).equals(b_[c].astype(str)) for c in ints)
        flo = [c for c in a_.columns if c not in ints and pd.api.types.is_float_dtype(a_[c])]
        max_diff = max([float(np.nanmax(np.abs(a_[c].to_numpy() - b_[c].to_numpy()))) if len(a_) else 0.0 for c in flo] or [0.0]) if same_shape else np.nan
        rows.append({"run": run, "file": f, "identical_bytes": (S07_PRE / run / f).read_bytes() == res("06_conformal", run, f).read_bytes(),
                     "counts_and_ks_equal": counts_equal, "max_abs_float_diff": max_diff})
s07_fixcmp = pd.DataFrame(rows)
display(s07_fixcmp)
print("certified k (primary, 100 splits) pre-fix == post-fix:",
      pd.read_csv(S07_PRE / "TRNSCRPT" / "certificate_validity.csv")["certified_k"].equals(
          pd.read_csv(res("06_conformal", "TRNSCRPT", "certificate_validity.csv"))["certified_k"]))
record("s07.n_cert_files_counts_changed", int((~s07_fixcmp["counts_and_ks_equal"]).sum()), "07",
       "certificate files whose error counts / certified k changed with the quantile fix")
del rows

# %% [markdown]
# **What this shows.** Every error count and every certified k is unchanged by the quantile fix, in
# both certificate runs. Where a file is not byte-identical, the only differences are in the last
# floating-point digits of the Clopper–Pearson bound (the `max_abs_float_diff` column), i.e. numerical
# noise in `scipy.stats.beta.ppf` between runs — not a change in any count or decision.
#
# **What it does not show.** That the certificate is *right*; only that it is independent of the
# conformal quantile, as it should be.
