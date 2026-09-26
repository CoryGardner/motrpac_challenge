# %% [markdown]
# ## 6. Conformal prediction sets: does "90 % coverage" hold on new animals?
#
# **Question.** A split-conformal wrapper turns the tissue classifier into a *set*-valued predictor with
# a promise: on new animals, the true tissue is in the set with probability at least 1 − α. Does that
# promise hold here, marginally (over all vials) and per tissue? Two scores are compared: **LAC**
# (the set holds every class with p̂ ≥ 1 − q̂; smallest sets, can be empty) and **APS** (classes are
# added in order of p̂ until the cumulative mass passes q̂; never empty, larger sets). Two
# calibrations use the *same* calibration animals: **pooled** (every vial of the 22 calibration
# animals) and **one vial per animal** (22 points, exchangeable by construction). On the pooled
# calibration two class-conditional variants are added: **Mondrian** (one quantile per tissue) and
# **floored Mondrian** (each tissue's quantile never below the marginal one).
#
# **Why it matters.** Coverage is the number a user of a tissue-ID tool would rely on; a marginal
# guarantee can hold while one tissue is almost never covered.
#
# **What to look for.** Marginal LAC coverage close to 1 − α; APS well above it (over-coverage) with
# sets a little larger than one class; per-tissue coverage far below 1 − α for the confusable tissues
# under marginal calibration, restored by Mondrian at the price of larger sets for those tissues.
#
# The settings are those of the pipeline's `make conformal` (first command): 5 outer animal-grouped
# folds, 22 calibration animals taken from each training fold, L2 logistic regression tuned on the
# remaining fit animals, α ∈ {0.05, 0.10, 0.20}. Only the conformal-set leg is run here (it is cheap);
# the 100-split certificate is section 7.
#
# **Note on the quantile.** On 2026-09-25 the library's `conformal_quantile` was changed to the
# textbook rule (the ⌈(n+1)(1−α)⌉-th smallest score, and an infinite threshold — the full set — when
# that rank exceeds n). The library pasted above is the fixed one; the cells below also print the
# pipeline's pre-fix numbers so the effect of the fix is visible.

# %%
section("6 conformal sets")

# -- the data: the same stacked TRNSCRPT matrix phases 04, 06 and 08 use (shared cache key) --
# from src/motrpac/cli.py::resolve_source (verbatim)
def resolve_source(assay: str, source: str) -> str:
    if source != "auto":
        return source
    if assay == "TRNSCRPT" and any((C.COUNTS_DIR).glob(f"{assay}__*.csv")):
        return "counts"
    return "norm"


pheno = cached("pheno", io.load_pheno)
om_trn = cached("TRNSCRPT_counts", lambda: io.stack_tissues(
    "TRNSCRPT", tissues=None, source=resolve_source("TRNSCRPT", "auto"), join="inner", pheno=pheno,
    complete=False, drop_incomplete_samples=None))

# -- the flags of `make conformal` (first command) --
s06_args = SimpleNamespace(label="tissue", model="logreg_l2", n_splits=5, prefilter=5000, seed=C.SEED, quick=False,
                           cal_frac=0.3, cal_animals=22, conditional=True)
X = om_trn.X.to_numpy(dtype=float)
y = om_trn.meta[s06_args.label].astype(str).to_numpy()
classes = sorted(np.unique(y))
y_idx = np.array([classes.index(v) for v in y])
g = om_trn.groups()
print(f"{X.shape[0]} vials × {X.shape[1]} genes, {len(np.unique(g))} animals, {len(classes)} tissues")


# from scripts/06_conformal_certify.py::proba_on (verbatim)
def proba_on(est, X, classes):
    p = est.predict_proba(X)
    order = [list(est.classes_).index(c) for c in classes]
    return p[:, order]


# from scripts/06_conformal_certify.py::split_train (verbatim)
def split_train(meta, tr, cal_frac, cal_animals, seed):
    """(fit, calibration) by animal: an absolute number of calibration animals if given, else a fraction."""
    if cal_animals is None:
        return fit_calibration_split(meta, tr, cal_frac, seed)
    tr = np.asarray(tr)
    rng = np.random.default_rng(seed)
    pids = meta.iloc[tr]["pid"].astype(str).to_numpy()
    uniq = np.unique(pids)
    rng.shuffle(uniq)
    n_cal = min(max(1, int(cal_animals)), len(uniq) - 1)
    cal = set(uniq[:n_cal])
    is_cal = np.array([p in cal for p in pids])
    fit_idx, cal_idx = tr[~is_cal], tr[is_cal]
    assert_no_group_leak(meta, fit_idx, cal_idx)
    return fit_idx, cal_idx

# %%
# from scripts/06_conformal_certify.py::main, the fold loop, conformal-set leg (a) only
# (trimmed: the repeat loop runs only rep == 0, which is the only repeat that computes sets;
#  the certificate leg (b) is section 7; `args.` → `s06_args.`)
t0 = time.perf_counter()
cov_rows, pc_rows = [], []
args = s06_args
for fold, (tr, te) in enumerate(grouped_kfold(om_trn.meta, args.label, args.n_splits, args.seed)):
    rep = 0
    seed_r = args.seed + 1000 * fold + rep
    fit_idx, cal_idx = split_train(om_trn.meta, tr, args.cal_frac, args.cal_animals, seed_r)
    cal_opa = cp.one_per_group(cal_idx, g, seed_r)          # one random vial per calibration animal
    # (a) full-feature model → conformal sets, calibrated pooled and one-per-animal on the same animals
    est = models.fit_tuned(args.model, X[fit_idx], y[fit_idx], g[fit_idx], k=None, prefilter=args.prefilter,
                           quick=args.quick, seed=args.seed)
    p_te = proba_on(est, X[te], classes)
    for mode, cidx in (("pooled", cal_idx), ("one_per_animal", cal_opa)):
        p_cal = proba_on(est, X[cidx], classes)
        for method in ("lac", "aps"):
            for alpha in C.ALPHAS:
                variants = [("marginal", False)] + ([("mondrian", True), ("floored", "floored")] if (args.conditional and mode == "pooled") else [])
                for cname, cond in variants:
                    r = cp.calibrate_and_evaluate(p_cal, y_idx[cidx], p_te, y_idx[te], classes, alpha, method,
                                                  extra_test=om_trn.meta.iloc[te][["sex", "group"]], conditional=cond)
                    cov_rows.append({"fold": fold, "calibration": mode, "conformal": cname, "method": method, "alpha": alpha,
                                     "qhat": r.qhat, "n_cal": r.n_cal, "n_cal_animals": int(om_trn.meta.iloc[cidx]["pid"].nunique()),
                                     "n_fit_animals": int(om_trn.meta.iloc[fit_idx]["pid"].nunique()),
                                     "n_test_animals": int(om_trn.meta.iloc[te]["pid"].nunique()), **r.overall})
                    pc_rows.append(r.per_class.assign(fold=fold, calibration=mode, conformal=cname, method=method, alpha=alpha))
    lac = [c for c in cov_rows if c["fold"] == fold and c["method"] == "lac" and c["alpha"] == 0.1 and c["conformal"] == "marginal"]
    print(f"  fold {fold}: {lac[0]['n_fit_animals']} fit / {lac[0]['n_cal_animals']} calibration / {lac[0]['n_test_animals']} test animals; "
          f"LAC α=0.1 coverage pooled={lac[0]['coverage']:.3f} one-per-animal={lac[1]['coverage']:.3f}")
del args
s06_cov = pd.DataFrame(cov_rows)
s06_pc = pd.concat(pc_rows, ignore_index=True)
(OUT / "s06").mkdir(parents=True, exist_ok=True)
s06_cov.to_csv(OUT / "s06" / "coverage.csv", index=False)
s06_pc.to_csv(OUT / "s06" / "per_class_coverage.csv", index=False)
print(f"conformal-set leg: {time.perf_counter() - t0:.0f} s")

# %% [markdown]
# ### 6a. Marginal coverage and set size
#
# The table below averages the five folds (the pipeline's own summary). Next to each recomputed value
# it prints the pipeline's value from `results/06_conformal/TRNSCRPT/coverage.csv` (post-fix; same
# seeds, so the two should agree to rounding) and the pre-fix value from the backup copy.

# %%
# from scripts/06_conformal_certify.py::main (verbatim: the cov_summary aggregation)
def s06_summary(cov):
    cov = cov[cov["conformal"] == "marginal"]
    return cov.groupby(["calibration", "method", "alpha"]).agg(
        coverage=("coverage", "mean"), coverage_sd=("coverage", "std"), avg_set_size=("avg_set_size", "mean"),
        frac_singleton=("frac_singleton", "mean"), frac_empty=("frac_empty", "mean"), n_cal=("n_cal", "mean")).reset_index()


S06_PRE = REPO / "backup" / "pipeline_history" / "results_pre_quantile_fix_2026-09-25" / "06_conformal" / "TRNSCRPT"
s06_pub = pd.read_csv(res("06_conformal", "TRNSCRPT", "coverage.csv"))     # post-fix pipeline run
s06_pre = pd.read_csv(S06_PRE / "coverage.csv")                              # pre-fix pipeline run
keys = ["calibration", "method", "alpha"]
s06_sum = s06_summary(s06_cov)
s06_cmp = (s06_sum[keys + ["coverage", "coverage_sd", "avg_set_size", "frac_empty", "n_cal"]]
           .merge(s06_summary(s06_pub)[keys + ["coverage", "avg_set_size"]].rename(columns=lambda c: c if c in keys else c + "_pub"), on=keys)
           .merge(s06_summary(s06_pre)[keys + ["coverage", "avg_set_size"]].rename(columns=lambda c: c if c in keys else c + "_prefix"), on=keys))
display(s06_cmp.round(3))
# the recomputation must match the post-fix pipeline run row by row (same folds, seeds and fits)
num = ["coverage", "avg_set_size", "qhat", "frac_empty"]
mrg = s06_cov.merge(s06_pub, on=["fold", "calibration", "conformal", "method", "alpha"], suffixes=("", "_pub"))
print(f"rows matched: {len(mrg)} of {len(s06_cov)}; max |recomputed − published| over coverage, set size, q̂, frac_empty: "
      f"{max(np.nanmax(np.abs(np.where(np.isinf(mrg[c]) & np.isinf(mrg[c + '_pub']), 0, mrg[c] - mrg[c + '_pub']))) for c in num):.2e}")


def s06_get(tbl, cal, method, alpha, col="coverage"):
    return float(tbl.loc[(tbl.calibration == cal) & (tbl.method == method) & (tbl.alpha == alpha), col].iloc[0])


for cal in ("pooled", "one_per_animal"):
    for method in ("lac", "aps"):
        tag = f"{'pooled' if cal == 'pooled' else 'opa'}_{method}_a10"
        record(f"s06.cov_{tag}", s06_get(s06_sum, cal, method, 0.1), "06", f"mean over 5 folds, {cal}, {method.upper()}, α=0.1")
        record(f"s06.size_{tag}", s06_get(s06_sum, cal, method, 0.1, "avg_set_size"), "06", f"avg set size, {cal}, {method.upper()}, α=0.1")
record("s06.cov_opa_lac_a05", s06_get(s06_sum, "one_per_animal", "lac", 0.05), "06", "one vial per animal, LAC, α=0.05")

# %% [markdown]
# **What this shows.** Pooled LAC coverage sits close to 1 − α at every α, with an empty-set rate that
# is essentially the miss rate: LAC misses a vial by returning no tissue at all rather than a wrong
# one. APS over-covers at every α, because it can never return an empty set, and pays in set size.
# With one vial per animal the calibration has only 22 points, so its coverage swings more across
# folds (larger sd). Compare the `_prefix` columns: where the fix moved a number, it moved the
# one-vial-per-animal rows (with 22 calibration points the old rule sat one rank higher) — the pooled
# rows, with hundreds of points, barely change.
#
# **What it does not show.** Marginal coverage says nothing about any single tissue (next cell), and
# it is averaged over five folds of about ten test animals each, so a fold-level value is noisy.

# %% [markdown]
# ### 6b. Per-tissue (conditional) coverage: marginal vs Mondrian vs floored
#
# **What to look for.** Tissues far below 1 − α under marginal calibration, and whether the
# class-conditional quantiles bring them back and at what set size.
#
# **These values differ from the published report.** They come from the pipeline after the
# 2026-09-25 conformal-quantile fix. The pre-fix code took a threshold one rank too high, so the
# report's per-tissue Mondrian and floored coverages were higher than the textbook rule supports. The
# report's claim that flooring keeps every tissue at or above its old minimum no longer holds. Every
# moved value, old and new, with the document sentence that quotes it, is in
# `docs/QUANTILE_FIX_DELTA.md`. The table below is compared with the *post-fix* results files.

# %%
# per-tissue coverage at LAC α = 0.1 (pooled calibration), mean over folds, for the three variants
# (after scripts/06_conformal_certify.py::main, the `mpt` table)
mp = s06_pc[(s06_pc["calibration"] == "pooled") & (s06_pc["method"] == "lac") & (s06_pc["alpha"] == 0.1)]
s06_mpt = mp.groupby(["y_true", "conformal"]).agg(coverage=("coverage", "mean"), set_size=("avg_set_size", "mean")).unstack("conformal")
s06_mpt.columns = [f"{a}_{b}" for a, b in s06_mpt.columns]
s06_mpt = s06_mpt.reset_index().sort_values("coverage_marginal")
# one-vial-per-animal marginal coverage for the same tissues
opa = s06_pc[(s06_pc["calibration"] == "one_per_animal") & (s06_pc["method"] == "lac") & (s06_pc["alpha"] == 0.1)]
s06_mpt = s06_mpt.merge(opa.groupby("y_true")["coverage"].mean().rename("coverage_opa_marginal").reset_index(), on="y_true")
display(s06_mpt.round(2))
pub_mpt = pd.read_csv(res("06_conformal", "TRNSCRPT", "per_tissue_marginal_vs_mondrian_alpha0.1.csv"))
chk = s06_mpt.merge(pub_mpt, on="y_true", suffixes=("", "_pub"))
print("max |recomputed − published| per-tissue coverage:",
      f"{max(np.abs(chk[c] - chk[c + '_pub']).max() for c in ['coverage_marginal', 'coverage_mondrian', 'coverage_floored']):.2e}")
under = s06_mpt[s06_mpt["coverage_marginal"] < 0.85]
print(f"tissues under 0.85 (marginal, pooled): {', '.join(under['y_true'])}")
# Mondrian / floored overall (pooled, LAC, α = 0.1): coverage and set size over all test vials
s06_var = (s06_cov[s06_cov.calibration == "pooled"].groupby(["conformal", "method", "alpha"])
           .agg(coverage=("coverage", "mean"), avg_set_size=("avg_set_size", "mean"), frac_empty=("frac_empty", "mean")).reset_index())
pre_var = pd.read_csv(S06_PRE / "coverage_marginal_vs_mondrian.csv")
display(s06_var.merge(pre_var, on=["conformal", "method", "alpha"], suffixes=("", "_prefix")).round(3))
record("s06.skmgn_cov_marginal", float(s06_mpt.set_index("y_true").loc["SKM-GN", "coverage_marginal"]), "06", "SKM-GN, pooled LAC α=0.1")
record("s06.skmgn_cov_mondrian", float(s06_mpt.set_index("y_true").loc["SKM-GN", "coverage_mondrian"]), "06", "SKM-GN, Mondrian LAC α=0.1")
record("s06.n_tissues_under085", len(under), "06", "tissues with marginal pooled LAC α=0.1 coverage < 0.85")
for v in ("mondrian", "floored"):
    row = s06_var[(s06_var.conformal == v) & (s06_var.method == "lac") & (s06_var.alpha == 0.1)].iloc[0]
    record(f"s06.cov_{v}_lac_a10", row["coverage"], "06", f"{v}, pooled, LAC, α=0.1")
    record(f"s06.size_{v}_lac_a10", row["avg_set_size"], "06", f"{v}, pooled, LAC, α=0.1")
# why the α = 0.05 Mondrian sets grew with the quantile fix: a class needs n ≥ n_min calibration vials for the
# ⌈(n+1)(1−α)⌉-th score to exist; below that its threshold is +inf and the class enters every set
for a_ in C.ALPHAS:
    n_min = next(n for n in range(1, 200) if np.ceil((n + 1) * (1 - a_) - 1e-9) <= n)
    print(f"α = {a_:.2f}: a class needs ≥ {n_min} calibration vials for a finite Mondrian threshold")
cal_counts = pd.Series(y[cal_idx]).value_counts()   # pooled calibration vials per tissue, last fold
print("pooled calibration vials per tissue (last fold): min", cal_counts.min(), "for", ", ".join(cal_counts[cal_counts == cal_counts.min()].index),
      "| others", cal_counts[cal_counts > cal_counts.min()].min(), "–", cal_counts.max())
v05 = s06_var[(s06_var.conformal == "mondrian") & (s06_var.method == "lac") & (s06_var.alpha == 0.05)].iloc[0]
record("s06.size_mondrian_lac_a05", v05["avg_set_size"], "06", "Mondrian, pooled, LAC, α=0.05 (moved by the quantile fix)")

# %%
# One figure: (A) coverage vs nominal 1 − α; (B) average set size; (C) per-tissue coverage at LAC α = 0.1.
fig, axes = plt.subplots(1, 3, figsize=(16, 4.6), gridspec_kw={"width_ratios": [1, 1, 2.2]})
styles = {("pooled", "marginal"): ("pooled, marginal", "-o"), ("one_per_animal", "marginal"): ("one vial/animal, marginal", "--o"),
          ("pooled", "mondrian"): ("pooled, Mondrian", "-s"), ("pooled", "floored"): ("pooled, floored Mondrian", ":^")}
colors = {"lac": "tab:blue", "aps": "tab:orange"}
agg = s06_cov.groupby(["calibration", "conformal", "method", "alpha"]).agg(coverage=("coverage", "mean"),
                                                                          size=("avg_set_size", "mean")).reset_index()
for (cal, cname), (lab, ls) in styles.items():
    for method in ("lac", "aps"):
        d = agg[(agg.calibration == cal) & (agg.conformal == cname) & (agg.method == method)].sort_values("alpha")
        axes[0].plot(1 - d["alpha"], d["coverage"], ls, color=colors[method], label=f"{method.upper()}, {lab}", ms=4, alpha=0.8)
        axes[1].plot(1 - d["alpha"], d["size"], ls, color=colors[method], ms=4, alpha=0.8)
axes[0].plot([0.78, 1.0], [0.78, 1.0], color="grey", lw=0.8)
axes[0].set(xlabel="nominal coverage 1 − α", ylabel="empirical coverage (mean of 5 folds)", title="A. Coverage on held-out animals")
axes[1].set(xlabel="nominal coverage 1 − α", ylabel="average set size (tissues)", title="B. Set size")
axes[0].legend(fontsize=7, frameon=False)
d = s06_mpt.sort_values("coverage_marginal")
xx = np.arange(len(d))
for off, col, lab in ((-0.27, "coverage_marginal", "marginal"), (0, "coverage_mondrian", "Mondrian"), (0.27, "coverage_floored", "floored")):
    axes[2].bar(xx + off, d[col], width=0.27, label=lab)
axes[2].axhline(0.9, color="red", ls="--", lw=1, label="1 − α = 0.90")
axes[2].set_xticks(xx, d["y_true"], rotation=90, fontsize=8)
axes[2].set(ylabel="coverage (pooled calibration)", ylim=(0, 1.05), title="C. Per-tissue coverage, LAC, α = 0.10")
axes[2].legend(fontsize=7, frameon=False, ncol=4, loc="lower center", bbox_to_anchor=(0.5, 1.06))
axes[2].set_title("C. Per-tissue coverage, LAC, α = 0.10", pad=24)
fig.tight_layout()
fig.savefig(OUT / "s06" / "conformal_summary.png", dpi=150)
plt.show()

# %% [markdown]
# **What this shows.** The marginal guarantee holds on average while a handful of tissues — the two
# skeletal muscles first, then the brain regions and vena cava, i.e. the tissues the classifier
# confuses with a neighbour — are covered far less often than 1 − α (panel C and the table above).
# Mondrian calibration pulls every tissue towards 1 − α: *up* for the confusable tissues (their own,
# looser threshold lets the true tissue into the set instead of leaving it empty) and *down* for the
# easy ones (their own threshold, from a couple of dozen well-classified vials, is tighter than the
# marginal one). At α = 0.10 the set sizes equal the coverages, i.e. the sets are singletons or empty
# under all three variants. The floored variant keeps the easy tissues at the marginal threshold, so
# it covers at least as well as both. After the quantile fix the α = 0.05 Mondrian and floored sets
# are much larger than before (the `_prefix` columns): the two sex-specific tissues have too few
# calibration vials (the printed minimum) for a finite threshold at α = 0.05, so under the textbook
# rule OVARY and TESTES enter *every* set — the honest answer "this calibration set cannot support
# α = 0.05 for these classes"; the old rule silently used their largest score instead.
#
# **What it does not show.** Per-tissue coverage here comes from about ten test animals per tissue
# per fold, so a single value can move by 0.1 with one animal. Pooled calibration treats the vials of
# one animal as exchangeable with other animals' vials; its test-animal coverage is the empirical
# check, not a proof. None of this is a guarantee under shift (section 8).
