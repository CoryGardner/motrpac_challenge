# %% [markdown]
# ## 8. Shift tests: does the fingerprint — and its coverage guarantee — survive a shift?
#
# **Question.** Sections 6–7 evaluated on held-out animals drawn from the same design. Here the test
# animals differ systematically from the training animals: train on one sex and test on the other
# (both ways); hold out a whole time point (8 w, 1 w); train on sedentary controls and test on every
# trained animal. For each split the pipeline fits the full-feature classifier and the k = 20 panel,
# calibrates LAC sets (α = 0.10) on 30 % of the *source* animals (pooled vials), and reports accuracy,
# in-distribution coverage on held-out source animals, and coverage on the target. It then asks
# whether calibrating on a few *target* animals (3 or 5, 10 random draws) buys the guarantee back.
#
# **Why it matters.** Conformal coverage is guaranteed only under exchangeability. A shift can break
# the guarantee while accuracy still looks fine — and a tool that reports "90 % sets" would then be
# silently wrong.
#
# **What to look for.** Accuracy that stays high while target coverage falls below source coverage
# (the guarantee breaks before the classifier does); unseen classes under the sex shift (OVARY and
# TESTES exist in one sex only; they are listed and excluded from the "seen" columns); whether a few
# target animals restore coverage.
#
# The TRNSCRPT splits are recomputed live (the pipeline phase took about a minute and a half). The
# metabolomics arm (9 core tissues, complete features only, k = 10) is small and is recomputed too.

# %%
section("8 shift tests")


# from scripts/08_shift_tests.py::coverage_with (verbatim)
def coverage_with(est, cls, X_te, y_te, qhat, method):
    sets = cp.predict_sets(est.predict_proba(X_te), qhat, method)
    cov, n_dropped = cp.coverage_seen(sets, y_te, cls)
    return cov, float(sets.sum(axis=1).mean()), n_dropped


# from scripts/08_shift_tests.py::eval_split (verbatim except one name: `models.balanced_accuracy_score` →
# `balanced_accuracy_score`, because the notebook's `models` namespace holds only the module's own definitions,
# not its imports; the sklearn function is the same object)
def eval_split(name, tr, te, X, y, g, meta, k, alpha, model, prefilter, quick, seed, method="lac",
               one_per_animal=False, recal_ns=(), recal_repeats=10):
    fit_idx, cal_idx = fit_calibration_split(meta, tr, 0.3, seed)
    cal_used = cp.one_per_group(cal_idx, g, seed) if one_per_animal else cal_idx
    rows, per_class, recal_rows, variant_rows = [], [], [], []
    for arm, kk in (("full", None), (f"panel_k{k}", k)):
        est = models.fit_tuned(model, X[fit_idx], y[fit_idx], g[fit_idx], k=kk, prefilter=prefilter, quick=quick, seed=seed)
        cls = list(est.classes_)
        unseen = sorted(set(y[te]) - set(cls))
        scores, n_cal_dropped = cp.calibration_scores(est.predict_proba(X[cal_used]), y[cal_used], cls, method)
        qhat = cp.conformal_quantile(scores, alpha)
        p_te = est.predict_proba(X[te])
        sets_te = cp.predict_sets(p_te, qhat, method)
        yi_te, seen_mask = cp.class_indices(y[te], cls)
        covered = np.zeros(len(te), dtype=bool)
        covered[seen_mask] = sets_te[np.flatnonzero(seen_mask), yi_te[seen_mask]]
        size = sets_te.sum(axis=1)
        # class-conditional variants from the same calibration scores: Mondrian and marginal-floor Mondrian
        yi_cal, seen_cal = cp.class_indices(y[cal_used], cls)
        q_variants = {"mondrian": cp.conformal_quantile_per_class(scores, yi_cal[seen_cal], alpha, len(cls), fallback=qhat),
                      "floored": cp.conformal_quantile_per_class(scores, yi_cal[seen_cal], alpha, len(cls), fallback=qhat, floor=qhat)}
        variant_sets = {"marginal": sets_te}
        variant_cov = {}
        for vname, qv in q_variants.items():
            sv = cp.predict_sets_conditional(p_te, qv, method)
            cv_ = np.zeros(len(te), dtype=bool)
            cv_[seen_mask] = sv[np.flatnonzero(seen_mask), yi_te[seen_mask]]
            variant_sets[vname], variant_cov[vname] = sv, cv_
        variant_cov["marginal"] = covered
        y_pred = est.predict(X[te])
        acc_all = float(np.mean(y_pred == y[te]))
        bal_seen = balanced_accuracy_score(y[te][seen_mask], y_pred[seen_mask]) if seen_mask.any() else np.nan
        # source in-distribution reference: calibrate on half of the calibration animals, test on the other half
        s_fit, s_cal = fit_calibration_split(meta, cal_idx, 0.5, seed + 1)
        s_fit_used = cp.one_per_group(s_fit, g, seed + 1) if one_per_animal else s_fit
        src_scores, n_src_fit_dropped = cp.calibration_scores(est.predict_proba(X[s_fit_used]), y[s_fit_used], cls, method)
        q_src = cp.conformal_quantile(src_scores, alpha)
        cov_src, size_src, n_src_cal_dropped = coverage_with(est, cls, X[s_cal], y[s_cal], q_src, method)
        n_src_dropped = n_src_fit_dropped + n_src_cal_dropped
        # empty-set rate (LAC) and APS sets under the same shift, same calibration animals
        frac_empty_te = float((sets_te.sum(axis=1) == 0).mean())
        sets_src_lac = cp.predict_sets(est.predict_proba(X[s_cal]), q_src, method)
        frac_empty_src = float((sets_src_lac.sum(axis=1) == 0).mean())
        aps_scores_cal, _ = cp.calibration_scores(est.predict_proba(X[cal_used]), y[cal_used], cls, "aps")
        q_aps = cp.conformal_quantile(aps_scores_cal, alpha)
        cov_aps_te, size_aps_te, _ = coverage_with(est, cls, X[te], y[te], q_aps, "aps")
        aps_src_scores, _ = cp.calibration_scores(est.predict_proba(X[s_fit_used]), y[s_fit_used], cls, "aps")
        _, size_aps_src, _ = coverage_with(est, cls, X[s_cal], y[s_cal], cp.conformal_quantile(aps_src_scores, alpha), "aps")
        if n_cal_dropped or n_src_dropped:
            print(f"  {name} {arm}: dropped {n_cal_dropped} calibration and {n_src_dropped} source-check samples "
                  f"whose class is absent from the fitted estimator ({', '.join(sorted(set(y[cal_idx]) - set(cls))) or '-'})")
        row = {"split": name, "arm": arm, "n_train_animals": int(meta.iloc[tr]["pid"].nunique()),
               "n_test_animals": int(meta.iloc[te]["pid"].nunique()), "n_test": len(te),
               "n_cal": int(len(scores)), "n_cal_animals": int(meta.iloc[cal_idx]["pid"].nunique()),
               "n_cal_dropped_unseen": int(n_cal_dropped), "n_src_dropped_unseen": int(n_src_dropped),
               "unseen_classes": ",".join(unseen), "accuracy_all": acc_all, "bal_acc_seen": bal_seen,
               "coverage_source_id": cov_src, "avg_set_size_source": size_src,
               "coverage_target_all": float(covered.mean()),
               "coverage_target_seen": float(covered[seen_mask].mean()) if seen_mask.any() else np.nan,
               "avg_set_size_target": float(size.mean()), "qhat": qhat, "alpha": alpha,
               "lac_frac_empty_source": frac_empty_src, "lac_frac_empty_target": frac_empty_te,
               "aps_coverage_target_seen": cov_aps_te, "aps_set_size_source": size_aps_src, "aps_set_size_target": size_aps_te}
        for vname in ("mondrian", "floored"):
            row[f"coverage_target_seen_{vname}"] = float(variant_cov[vname][seen_mask].mean()) if seen_mask.any() else np.nan
            row[f"avg_set_size_target_{vname}"] = float(variant_sets[vname].sum(axis=1).mean())
        for vname, sv in variant_sets.items():
            variant_rows.append(pd.DataFrame({"split": name, "arm": arm, "conformal": vname, "y_true": y[te], "covered": variant_cov[vname],
                                              "size": sv.sum(axis=1), "seen": seen_mask})
                                .groupby(["split", "arm", "conformal", "y_true"]).agg(coverage=("covered", "mean"), avg_set_size=("size", "mean"),
                                                                                       n=("covered", "size"), seen=("seen", "all")).reset_index())
        # recalibration on N target animals, test on the remaining target animals (same animals for both quantiles)
        rng = np.random.default_rng(seed + 7)
        t_pids = np.unique(g[te])
        for n_recal in recal_ns:
            if n_recal >= len(t_pids):
                continue
            cov_t, cov_s, sz_t, sz_s, dropped = [], [], [], [], []
            for _ in range(recal_repeats):
                chosen = set(rng.choice(t_pids, size=n_recal, replace=False))
                m_cal = np.array([p in chosen for p in g[te]])
                tcal, ttest = te[m_cal], te[~m_cal]
                assert_no_group_leak(meta, tcal, ttest)
                t_scores, n_drop = cp.calibration_scores(est.predict_proba(X[tcal]), y[tcal], cls, method)
                q_t = cp.conformal_quantile(t_scores, alpha)
                c_t, s_t, _ = coverage_with(est, cls, X[ttest], y[ttest], q_t, method)
                c_s, s_s, _ = coverage_with(est, cls, X[ttest], y[ttest], qhat, method)
                cov_t.append(c_t); cov_s.append(c_s); sz_t.append(s_t); sz_s.append(s_s); dropped.append(n_drop)
            row[f"cov_target_recal_N{n_recal}"] = float(np.nanmean(cov_t))
            row[f"cov_target_srccal_N{n_recal}"] = float(np.nanmean(cov_s))
            recal_rows.append({"split": name, "arm": arm, "n_recal_animals": n_recal, "repeats": recal_repeats,
                               "coverage_target_recalibrated": float(np.nanmean(cov_t)), "coverage_target_recal_sd": float(np.nanstd(cov_t)),
                               "coverage_target_source_cal_same_test": float(np.nanmean(cov_s)),
                               "set_size_recalibrated": float(np.mean(sz_t)), "set_size_source_cal": float(np.mean(sz_s)),
                               "n_recal_dropped_unseen_mean": float(np.mean(dropped))})
        rows.append(row)
        per_class.append(pd.DataFrame({"split": name, "arm": arm, "y_true": y[te], "covered": covered, "size": size})
                         .groupby(["split", "arm", "y_true"]).agg(coverage=("covered", "mean"), avg_set_size=("size", "mean"),
                                                                   n=("covered", "size")).reset_index())
    return rows, pd.concat(per_class), recal_rows, pd.concat(variant_rows)


# from scripts/08_shift_tests.py::main (trimmed: argparse → explicit flags of `make shift`; report and plot dropped;
# the table/per-class assembly is kept verbatim)
def s08_run(om, k, alpha=0.1, model="logreg_l2", prefilter=5000, seed=C.SEED, holdout_groups="8w,1w",
            recal_ns=(3, 5), recal_repeats=10, one_per_animal=False, label="tissue"):
    X = om.X.to_numpy(dtype=float)
    y = om.meta[label].astype(str).to_numpy()
    g = om.groups()
    splits = list(leave_one_sex_out(om.meta))
    splits += list(leave_one_group_out(om.meta, "group", [h.strip() for h in holdout_groups.split(",")]))
    splits += list(train_controls_test_trained(om.meta))
    rows, pcs, recals, variants = [], [], [], []
    for name, (tr, te) in splits:
        r, pc, rc, vr = eval_split(name, tr, te, X, y, g, om.meta, k, alpha, model, prefilter, False,
                                   seed, one_per_animal=one_per_animal, recal_ns=recal_ns, recal_repeats=recal_repeats)
        rows.extend(r); pcs.append(pc); recals.extend(rc); variants.append(vr)
        for rr in r:
            print(f"  {name:28s} {rr['arm']:10s} acc={rr['accuracy_all']:.3f} cov_src={rr['coverage_source_id']:.3f} "
                  f"cov_target={rr['coverage_target_all']:.3f} (seen: {rr['coverage_target_seen']:.3f}) "
                  f"size={rr['avg_set_size_target']:.2f} unseen={rr['unseen_classes'] or '-'}")
    table = pd.DataFrame(rows)
    table["coverage_drop"] = table["coverage_source_id"] - table["coverage_target_seen"]
    table["set_size_change"] = table["avg_set_size_target"] - table["avg_set_size_source"]
    return table, pd.concat(pcs, ignore_index=True), pd.DataFrame(recals), pd.concat(variants, ignore_index=True)

# %%
# TRNSCRPT, `make shift` first command: --assay TRNSCRPT --k 20 --recalibrate-target 3,5
pheno = cached("pheno", io.load_pheno)
om_trn = cached("TRNSCRPT_counts", lambda: io.stack_tissues(
    "TRNSCRPT", tissues=None, source=resolve_source("TRNSCRPT", "auto"), join="inner", pheno=pheno,
    complete=False, drop_incomplete_samples=None))
t0 = time.perf_counter()
s08_trn, s08_trn_pc, s08_trn_recal, s08_trn_var = s08_run(om_trn, k=20)
(OUT / "s08").mkdir(parents=True, exist_ok=True)
s08_trn.to_csv(OUT / "s08" / "TRNSCRPT_shift_table.csv", index=False)
print(f"TRNSCRPT shift splits: {time.perf_counter() - t0:.0f} s")

# %% [markdown]
# ### 8a. TRNSCRPT, k = 20: accuracy vs coverage under each shift
#
# The first table is the recomputed shift table. The second puts the target coverage next to the
# pipeline's post-fix file (`results/08_shift/TRNSCRPT/shift_table.csv`, which this run should
# reproduce) and the pre-fix backup, so the effect of the quantile fix is visible split by split.

# %%
cols = ["split", "arm", "n_train_animals", "n_test_animals", "n_cal", "unseen_classes", "accuracy_all", "bal_acc_seen",
        "coverage_source_id", "coverage_target_seen", "coverage_drop", "avg_set_size_target", "lac_frac_empty_target",
        "cov_target_recal_N3", "cov_target_recal_N5"]
display(s08_trn[cols].round(3))

S08_PRE = REPO / "backup" / "pipeline_history" / "results_pre_quantile_fix_2026-09-25" / "08_shift"


# notebook-only helper (not in the pipeline): presentation of tables the pipeline wrote
def s08_compare(table, run):
    """Recomputed vs the pipeline's post-fix shift_table.csv vs the pre-fix backup, on the headline columns."""
    pub = pd.read_csv(res("08_shift", run, "shift_table.csv"))
    pre = pd.read_csv(S08_PRE / run / "shift_table.csv")
    heads = ["accuracy_all", "coverage_source_id", "coverage_target_seen", "cov_target_recal_N3", "qhat"]
    out = table[["split", "arm"] + heads]
    out = out.merge(pub[["split", "arm"] + heads], on=["split", "arm"], suffixes=("", "_pub"))
    out = out.merge(pre[["split", "arm"] + heads].rename(columns={h: h + "_prefix" for h in heads}), on=["split", "arm"])
    diff = max(float(np.nanmax(np.abs(out[h] - out[h + "_pub"]))) for h in heads)
    print(f"{run}: max |recomputed − published (post-fix)| over {', '.join(heads)}: {diff:.2e}")
    return out[["split", "arm"] + [c for h in heads for c in (h, h + "_pub", h + "_prefix")]]


s08_trn_cmp = s08_compare(s08_trn, "TRNSCRPT")
display(s08_trn_cmp.drop(columns=[c for c in s08_trn_cmp.columns if c.startswith("qhat")]).round(3))
short = {"train_male_test_female": "m2f", "train_female_test_male": "f2m", "holdout_group_8w": "8w",
         "holdout_group_1w": "1w", "train_control_test_trained": "ctrl2trained"}
for r in s08_trn.itertuples():
    tag = f"s08.trn_{short[r.split]}_{'k20' if r.arm == 'panel_k20' else 'full'}"
    record(tag + "_acc", r.accuracy_all, "08", f"TRNSCRPT {r.split} {r.arm}: accuracy (all test vials)")
    record(tag + "_cov", r.coverage_target_seen, "08", f"TRNSCRPT {r.split} {r.arm}: target coverage, seen classes")
    record(tag + "_covsrc", r.coverage_source_id, "08", f"TRNSCRPT {r.split} {r.arm}: in-distribution source coverage")
    record(tag + "_recalN3", r.cov_target_recal_N3, "08", f"TRNSCRPT {r.split} {r.arm}: coverage after recalibrating on 3 target animals")

# %% [markdown]
# **What this shows.** Under the sex shift, balanced accuracy on the seen classes stays high while
# target coverage falls below the source's in-distribution coverage — the guarantee breaks before
# the classifier does, most clearly male→female. Each sex split also has an unseen class (OVARY or
# TESTES), which no source-calibrated set can cover. Held-out 8 w animals look like an ordinary
# cross-validation fold; held-out 1 w and controls→trained lose less than male→female. Under every
# shift the LAC misses are mostly empty sets (`lac_frac_empty_target`), rarely wrong singletons. The `_prefix`
# columns show how the quantile fix moved each value: the source calibration sets hold from a few
# dozen to a couple of hundred vials (`n_cal`), so the fix moves the threshold by one rank and most
# target coverages by a point or two (a few by more); the recalibration on 3 target animals, with
# far fewer vials, moves more (next table).
#
# **What it does not show.** Each split is one draw of calibration animals (the recalibration rows
# are averages over 10 draws); the sex shift also mixes sex with arrival cohort (the design confound
# described in the project notes), so "sex shift" here means "sex and cohort shift".

# %% [markdown]
# ### 8b. Recalibration on 3 and 5 target animals
#
# **What to look for.** On the *same* remaining target animals: coverage with the source-calibrated
# threshold vs with a threshold calibrated on N target animals.

# %%
recal_cols = ["split", "arm", "n_recal_animals", "coverage_target_source_cal_same_test", "coverage_target_recalibrated",
              "coverage_target_recal_sd", "set_size_source_cal", "set_size_recalibrated"]
pre_recal = pd.read_csv(S08_PRE / "TRNSCRPT" / "shift_recalibration.csv")
s08_recal = s08_trn_recal[recal_cols].merge(
    pre_recal[["split", "arm", "n_recal_animals", "coverage_target_recalibrated"]].rename(columns={"coverage_target_recalibrated": "recalibrated_prefix"}),
    on=["split", "arm", "n_recal_animals"])
display(s08_recal.round(3))
print("range of recalibrated coverage (N = 3):",
      f"{s08_trn_recal.loc[s08_trn_recal.n_recal_animals == 3, 'coverage_target_recalibrated'].min():.3f} – "
      f"{s08_trn_recal.loc[s08_trn_recal.n_recal_animals == 3, 'coverage_target_recalibrated'].max():.3f}")

# %% [markdown]
# **What this shows.** On the same remaining target animals, a threshold calibrated on only 3 target
# animals brings coverage back to about 1 − α in every split and arm, where the source-calibrated
# threshold fell short; 5 animals do no better on average. Before the fix the recalibrated coverage
# was a few points *above* 1 − α (`recalibrated_prefix`): with a few dozen calibration vials, the
# old rule's extra rank over-covered. After the fix it sits close to the nominal level, as it should.
#
# **What it does not show.** Each value is a mean over 10 random draws of target animals (the sd
# column is the draw-to-draw spread), not a guarantee for one draw.

# %% [markdown]
# ### 8c. Metabolomics, 9 core tissues (k = 10): the sharpest example
#
# `make shift`, second command: `--assay METAB --k 10 --tissues <core 9> --drop-incomplete-samples 0.2
# --complete-features --recalibrate-target 3,5`. Recomputed live (small matrix).

# %%
CORE9 = "SKM-GN,HEART,KIDNEY,LIVER,LUNG,BAT,WAT-SC,HIPPOC,PLASMA"


# from src/motrpac/cli.py::parse_tissues (verbatim)
def parse_tissues(arg: str | None) -> list[str] | None:
    if not arg:
        return None
    return [C.tissue_from_token(C.tissue_token(t.strip())) for t in arg.split(",") if t.strip()]


om_met9 = cached("METAB_core9_complete", lambda: io.stack_tissues(
    "METAB", tissues=parse_tissues(CORE9), source=resolve_source("METAB", "auto"), join="inner", pheno=pheno,
    complete=True, drop_incomplete_samples=0.2))
print(f"METAB core 9: {om_met9.X.shape[0]} vials × {om_met9.X.shape[1]} features, {om_met9.meta['pid'].nunique()} animals")
s08_met, s08_met_pc, s08_met_recal, _ = s08_run(om_met9, k=10)
s08_met.to_csv(OUT / "s08" / "METAB_core9_shift_table.csv", index=False)
s08_met_cmp = s08_compare(s08_met, "METAB_core9")
display(s08_met_cmp.drop(columns=[c for c in s08_met_cmp.columns if c.startswith("qhat")]).round(3))
m2f = s08_met.set_index(["split", "arm"]).loc[("train_male_test_female", "full")]
print(f"METAB male→female, full model: accuracy {m2f['accuracy_all']:.3f}, target coverage {m2f['coverage_target_seen']:.3f}, "
      f"after recalibration on 3 female animals {m2f['cov_target_recal_N3']:.3f}")
low = s08_met_pc[(s08_met_pc.split == "train_male_test_female") & (s08_met_pc.arm == "full")].sort_values("coverage")
print("per-tissue coverage, male→female, full model:", ", ".join(f"{r.y_true} {r.coverage:.2f}" for r in low.itertuples()))
for r in s08_met.itertuples():
    tag = f"s08.met_{short[r.split]}_{'k10' if r.arm == 'panel_k10' else 'full'}"
    record(tag + "_acc", r.accuracy_all, "08", f"METAB core9 {r.split} {r.arm}: accuracy")
    record(tag + "_cov", r.coverage_target_seen, "08", f"METAB core9 {r.split} {r.arm}: target coverage, seen classes")
record("s08.met_m2f_full_recalN3", m2f["cov_target_recal_N3"], "08", "METAB core9 male→female full: recalibrated on 3 animals")

# %%
# One figure: target coverage per split — source-calibrated vs in-distribution vs recalibrated on 3 target animals.
fig, axes = plt.subplots(1, 2, figsize=(15, 4.6), sharey=True)
for ax, (tbl, title) in zip(axes, ((s08_trn, "TRNSCRPT (panel k = 20 and full)"), (s08_met, "METAB core 9 (panel k = 10 and full)"))):
    lab = tbl["split"].map(short).str.replace("ctrl2trained", "ctrl→tr") + "\n" + tbl["arm"].str.replace("panel_", "")
    xx = np.arange(len(tbl))
    ax.bar(xx - 0.27, tbl["coverage_source_id"], 0.27, label="source animals (in distribution)", color="grey")
    ax.bar(xx, tbl["coverage_target_seen"], 0.27, label="target, source-calibrated", color="tab:red")
    ax.bar(xx + 0.27, tbl["cov_target_recal_N3"], 0.27, label="target, recalibrated on 3 target animals", color="tab:blue")
    ax.plot(xx, tbl["accuracy_all"], "k^", label="accuracy (all target vials)")
    ax.axhline(0.9, color="red", ls="--", lw=1)
    ax.set_xticks(xx, lab, fontsize=8)
    ax.set(title=title, ylim=(0.5, 1.02))
axes[0].set_ylabel("coverage (LAC, α = 0.10) / accuracy")
axes[0].legend(fontsize=7, frameon=False, loc="lower left")
fig.tight_layout()
fig.savefig(OUT / "s08" / "shift_summary.png", dpi=150)
plt.show()

# %% [markdown]
# **What this shows.** Metabolomics makes the point most sharply: a model trained on males classifies
# every female vial correctly, yet the source-calibrated sets cover the female vials far less often
# than 1 − α, with some tissues almost never covered (the per-tissue line above). The reason is that
# the source-calibrated threshold encodes the source's *confidence*, and the classifier is less
# confident — though still right — on the other sex; LAC then returns empty sets. Recalibrating on
# 3 target animals restores coverage to about 1 − α or above in every row, for both assays, and 5
# animals do not do better — a few target animals buy the guarantee back. Compare the post-fix and
# pre-fix columns: the fix lowered several source-calibrated target coverages (most visibly
# metabolomics male→female, full model), which sharpens the headline rather than changing it.
#
# **What it does not show.** Why the other sex is less confident (sex biology, arrival cohort, or
# batch — they are confounded in this design). Coverage after recalibration is an average over 10
# draws of 3 animals from a small target, not a guarantee for any one draw; and one-vial-per-animal
# calibration, the exchangeable version, is not usable at these sizes (too few calibration points).
