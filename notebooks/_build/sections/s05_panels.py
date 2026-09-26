# %% [markdown]
# ## 5. Compact panels: how few genes identify a tissue?
#
# **Question.** If the classifier may use only *k* genes, chosen inside each training fold, how does
# accuracy grow with *k*? And does *how* the genes are chosen matter?
#
# **Why it matters.** A small panel is what a hackathon project could actually explain, transfer to
# another dataset (sections 9–10) or certify (section 6). The selector is the design choice that
# matters most: the obvious one, `SelectKBest(f_classif)`, ranks genes by a single multiclass F-score.
# The pipeline's round-robin selector instead ranks genes per tissue (one-vs-rest effect size) and takes
# the best remaining gene of each tissue in turn.
#
# **What to look for.** (1) The two selectors at k = 20 on the same folds. (2) Where the round-robin
# curve flattens. (3) Which tissues are still confused at k = 20 — derived from the predictions, not
# typed in. (4) Whether the genes chosen are stable across animal bootstraps, and what risk flags they
# carry.
#
# Design (Makefile target `panels`, TRNSCRPT line): classifier `logreg_l2` on the *k* selected genes,
# the same 5 animal-grouped folds as section 4, variance prefilter 5,000 → scaling → selector, all inside
# each fold; `--compare-selectors --stability-k 20`, 50 animal bootstraps for stability.
# **Default mode** runs a reduced grid live (k = 5, 10, 20, 50 round-robin; k = 20 F-test) and loads
# the full published grid for the figure; `RECOMPUTE=True` runs the full grid for both selectors.

# %%
section("5 panels")
om = S04_OM                                   # the TRNSCRPT matrix loaded in section 4 (same cached object)
S05_OUT = OUT / "s05"
S05_OUT.mkdir(parents=True, exist_ok=True)
S05_GRID_FULL = C.PANEL_GRID                  # --grid default
S05_GRID = S05_GRID_FULL if RECOMPUTE else [5, 10, 20, 50]
S05_FC_GRID = S05_GRID_FULL if RECOMPUTE else [20]

# from scripts/05_compact_panels.py::main (trimmed: argparse defaults inlined — --model logreg_l2,
# --selector roundrobin, n_splits 5, prefilter 5000; grid reduced in default mode; no REPORT.md)
s05_curve, s05_selected, s05_preds = models.panel_curve(
    om, label="tissue", grid=S05_GRID, kind="logreg_l2", n_splits=5, prefilter=5000, quick=False,
    seed=C.SEED, selector="roundrobin", return_predictions=True, verbose=False)
# --compare-selectors: the same folds and classifier, F-test selector
s05_curve_fc, _ = models.panel_curve(om, label="tissue", grid=S05_FC_GRID, kind="logreg_l2", n_splits=5,
                                     prefilter=5000, quick=False, seed=C.SEED, selector="fclassif", verbose=False)
s05_curve.to_csv(S05_OUT / "panel_curve.csv", index=False)
s05_curve_fc.to_csv(S05_OUT / "panel_curve_fclassif.csv", index=False)

# recomputed vs published, per k and selector (mean balanced accuracy over the 5 folds)
pub_rr = pd.read_csv(res("05_panels", "TRNSCRPT", "panel_curve.csv"))
pub_fc = pd.read_csv(res("05_panels", "TRNSCRPT", "panel_curve_fclassif.csv"))
agg = lambda d: d.groupby("k")["balanced_accuracy"].agg(["mean", "std"])
s05_cmp = pd.concat({
    ("roundrobin", "recomputed"): agg(s05_curve)["mean"], ("roundrobin", "published"): agg(pub_rr)["mean"],
    ("roundrobin", "published_sd"): agg(pub_rr)["std"],
    ("fclassif", "recomputed"): agg(s05_curve_fc)["mean"], ("fclassif", "published"): agg(pub_fc)["mean"],
}, axis=1).sort_index()
display(s05_cmp.round(3))
rr20, fc20 = agg(s05_curve).loc[20, "mean"], agg(s05_curve_fc).loc[20, "mean"]
print(f"k = 20: round-robin {rr20:.3f} vs F-test {fc20:.3f} (published {agg(pub_rr).loc[20, 'mean']:.3f} vs {agg(pub_fc).loc[20, 'mean']:.3f})")
record("s05.bal_acc_rr_k20", rr20, "05")
record("s05.bal_acc_fclassif_k20", fc20, "05")
record("s05.bal_acc_rr_k50", agg(s05_curve).loc[50, "mean"], "05")
record("s05.bal_acc_rr_k10", agg(s05_curve).loc[10, "mean"], "05")

# how many different tissues do the 20 F-test genes of fold 0 "belong to"? (their top one-vs-rest tissue)
# — the mechanism behind the gap. Recomputed on fold 0's training data, the same steps as the pipeline.
from sklearn.impute import SimpleImputer
from sklearn.preprocessing import StandardScaler
tr0 = S04_FOLDS[0][0]
X0 = om.X.to_numpy(dtype=float)[tr0]
y0 = om.meta["tissue"].astype(str).to_numpy()[tr0]
Z0 = StandardScaler().fit_transform(models.VarianceTopK(5000).fit_transform(SimpleImputer(strategy="median").fit_transform(X0)))
fc_sel = models.make_selector("fclassif", 20).fit(Z0, y0)
rr_sel = models.make_selector("roundrobin", 20).fit(Z0, y0)
owner = rr_sel.scores_.argmax(axis=0)          # tissue with the largest one-vs-rest score, per gene
for name, idx in (("F-test", np.flatnonzero(fc_sel.get_support())), ("round-robin", rr_sel.idx_)):
    tissues_hit = pd.Series(rr_sel.classes_[owner[idx]]).value_counts()
    print(f"fold 0, {name:11s}: 20 genes are top markers of {tissues_hit.size} of {len(rr_sel.classes_)} tissues; "
          f"most common: {tissues_hit.index[0]} ({tissues_hit.iloc[0]} genes)")

# %%
# The figure: the full published grid (default) or this run's full grid (RECOMPUTE), both selectors.
full_rr, full_fc = (s05_curve, s05_curve_fc) if RECOMPUTE else (pub_rr, pub_fc)
fig, ax = plt.subplots(figsize=(7, 4.3))
for d, lab, c in ((full_rr, "round-robin (per-tissue) selector", "C0"), (full_fc, "F-test (SelectKBest f_classif)", "C3")):
    a = agg(d)
    ax.errorbar(a.index, a["mean"], yerr=a["std"].fillna(0), marker="o", capsize=3, color=c, label=lab)
if not RECOMPUTE:   # overlay this run's live points on the published curves
    ax.scatter(agg(s05_curve).index, agg(s05_curve)["mean"], marker="x", s=70, color="k", zorder=5, label="recomputed in this run")
    ax.scatter(agg(s05_curve_fc).index, agg(s05_curve_fc)["mean"], marker="x", s=70, color="k", zorder=5)
ax.axhline(1 / om.meta["tissue"].nunique(), color="grey", ls=":", label="chance")
ax.set_xscale("log"); ax.set_xlabel("panel size k (genes)"); ax.set_ylabel("balanced accuracy (mean ± sd over folds)")
ax.set_title("TRNSCRPT, 19 tissues: accuracy vs panel size (logreg_l2)"); ax.grid(alpha=0.3); ax.legend(frameon=False, fontsize=8)
fig.tight_layout(); fig.savefig(S05_OUT / "panel_curve_both_selectors.png", dpi=120); plt.show()
k99 = agg(full_rr).query("mean >= 0.99").index.min()
print("smallest k with mean balanced accuracy >= 0.99 (round-robin, point estimate):", k99)
record("s05.smallest_k_bal_acc_ge_099", k99, "05")

# %% [markdown]
# **What this shows.** With the round-robin selector a panel of a couple of dozen genes already comes
# close to the all-gene baseline of section 4, and the curve is flat from there on. With the F-test
# selector the same 20 genes fail badly: the fold-0 diagnosis above shows why — the top F-scores
# concentrate on a few tissues with huge, easy markers, leaving most tissues without a gene of their own.
# The live points sit on the published curve.
#
# **What it does not show.** These are point estimates over 5 folds of 10 test animals each; the sd bars
# are fold-to-fold spread, not a guarantee. Which *k* is enough with a stated error rate is section 6's
# job (conformal certification). And "smallest k ≥ 0.99" depends on the grid points that were tried.

# %%
# Confusion at k = 20, from this run's pooled out-of-fold predictions
# from scripts/05_compact_panels.py::main (verbatim confusion block, --confusion-k 20 only)
classes = sorted(om.meta["tissue"].astype(str).unique())
pk = s05_preds[s05_preds["k"] == 20]
s05_conf20 = pd.crosstab(pk["y_true"], pk["y_pred"]).reindex(index=classes, columns=classes, fill_value=0)
s05_conf20.to_csv(S05_OUT / "confusion_k20.csv")
rows = [{"k": 20, "true": a, "predicted": b, "count": int(s05_conf20.loc[a, b]),
         "frac_of_true": float(s05_conf20.loc[a, b] / max(s05_conf20.loc[a].sum(), 1))}
        for a in classes for b in classes if a != b and s05_conf20.loc[a, b] > 0]
s05_pairs20 = pd.DataFrame(rows, columns=["k", "true", "predicted", "count", "frac_of_true"]).sort_values("count", ascending=False)

pub_conf20 = pd.read_csv(res("05_panels", "TRNSCRPT", "confusion_k20.csv"), index_col=0)
print("confusion matrix at k = 20 identical to published confusion_k20.csv:",
      s05_conf20.equals(pub_conf20.reindex(index=classes, columns=classes)))
n_err = int(s05_pairs20["count"].sum())
print(f"{n_err} errors of {len(pk)} test samples (pooled over folds); top confused pairs:")
display(s05_pairs20.head(8).reset_index(drop=True).round(2))
top = s05_pairs20.iloc[0]
share = top["count"] / n_err
print(f"largest single error: {top['true']} -> {top['predicted']}: {top['count']} samples "
      f"({top['frac_of_true']:.0%} of all {top['true']} samples; {share:.0%} of all errors)")
record("s05.n_errors_k20", n_err, "05")
record("s05.top_error_count_k20", top["count"], "05", f"{top['true']}->{top['predicted']}")
show_png(plots.confusion_heatmap(s05_conf20, "TRNSCRPT k = 20 round-robin panel: confusion (row-normalized)",
                                 S05_OUT / "confusion_k20.png"), width=600)

# %% [markdown]
# **What this shows.** At k = 20 the errors are few and concentrated: the table is sorted by count, so
# its first rows are the pairs that dominate. The largest one sends vena cava samples to brown adipose
# tissue — consistent with the brown-fat contamination of some vena cava vials that the consortium
# flagged as outliers (see the workspace notes; many of those vials remain in the shipped matrices).
# Most of the rest are between the two skeletal muscles and the two brain regions, which are genuinely similar.
#
# **What it does not show.** This is one set of 5 folds; with 10 test animals per fold, a single
# contaminated animal moves several counts. Whether the vena-cava errors are exactly the flagged vials
# is not checked here.

# %% [markdown]
# ### 5b. Is the panel stable, and what is in it?
#
# Stability selection refits the selection step (impute → prefilter → scale → round-robin, k = 20) on 50
# bootstrap resamples of *animals* and counts how often each gene is chosen. Genes chosen in at least
# 80 % of resamples form the candidate panel. Each is then annotated (script `05_annotate_panel.py`):
# its marker tissue, effect size, whether it is training-regulated in that tissue (a risk for the
# controls → trained shift, section 8) and whether it correlates with a library QC metric (`pct_mrna`)
# within its tissue (a risk for transfer to another library chemistry). Flags are annotation, not
# exclusion.

# %%
if RECOMPUTE:
    # from scripts/05_compact_panels.py::main (verbatim call; --stability-k 20, --n-boot 50 default)
    s05_stab = models.stability_selection(om, label="tissue", k=20, n_boot=50, prefilter=5000,
                                          seed=C.SEED, selector="roundrobin")
    s05_stab.to_csv(S05_OUT / "stability_k20.csv", index=False)
else:
    s05_stab = pd.read_csv(res("05_panels", "TRNSCRPT", "stability_k20.csv"), dtype={"feature_ID": str})
s05_stab["gene_symbol"] = io.map_to_gene_symbols(s05_stab["feature_ID"]).values
s05_panel = s05_stab[s05_stab["selection_frequency"] >= 0.8].reset_index(drop=True)   # candidate_panel.csv
pub_panel = pd.read_csv(res("05_panels", "TRNSCRPT", "candidate_panel.csv"), dtype={"feature_ID": str})
print(f"{len(s05_stab)} genes selected at least once in the bootstraps; {len(s05_panel)} in >= 80 % "
      f"(published candidate panel: {len(pub_panel)}; same genes: {set(s05_panel['feature_ID']) == set(pub_panel['feature_ID'])})")
record("s05.n_candidate_panel", len(s05_panel), "05")

# %%
# from scripts/05_annotate_panel.py::main (trimmed: argparse defaults inlined — qc_col pct_mrna,
# qc_risk 0.5, label tissue; panel = the candidate panel above; no REPORT.md). Annotation only.
qc_col, qc_risk = "pct_mrna", 0.5
panel = s05_panel.drop_duplicates("feature_ID").reset_index(drop=True)
y = om.meta["tissue"].astype(str).to_numpy()
classes = sorted(np.unique(y))
ids = [f for f in panel["feature_ID"] if f in om.X.columns]
full = SimpleImputer(strategy="median").fit_transform(om.X.to_numpy(dtype=float))   # full matrix: annotation only
sel = models.RoundRobinSelector(k=len(ids)).fit(full, y)
scores = sel.scores_
col = {f: j for j, f in enumerate(om.X.columns)}
rank_pos = {i: {j: r + 1 for r, j in enumerate(np.argsort(-scores[i]))} for i in range(len(classes))}
top_full = set(np.asarray(om.X.columns)[sel.idx_])
mu = om.X[ids].groupby(y).mean()
reg = pd.read_csv(C.RAW_DIR / "training_regulated_features.csv", dtype=str, usecols=["feature_ID", "assay", "tissue"])
reg_map = reg[reg["assay"] == "TRNSCRPT"].groupby("feature_ID")["tissue"].agg(lambda s: sorted(set(s))).to_dict()
m = pd.read_csv(C.META_DIR / "TRNSCRPT.csv", dtype=str, low_memory=False)
qc = pd.to_numeric(m.drop_duplicates("viallabel").set_index("viallabel")[qc_col], errors="coerce").reindex(om.X.index)
symbols = io.map_to_gene_symbols(ids)
rows = []
for f in ids:
    j = col[f]; sc = scores[:, j]; order = np.argsort(-sc)
    i1, i2 = int(order[0]), int(order[1])
    t1 = classes[i1]; means = mu[f]; others = means.drop(t1); x = om.X[f]
    r_in, r_max, t_max = np.nan, np.nan, None
    for t in classes:
        mk = (y == t) & x.notna().to_numpy() & qc.notna().to_numpy()
        if mk.sum() >= 8 and x[mk].std() > 0 and qc[mk].std() > 0:
            r = float(np.corrcoef(x[mk], qc[mk])[0, 1])
            if t == t1:
                r_in = r
            if np.isnan(r_max) or abs(r) > abs(r_max):
                r_max, t_max = r, t
    reg_t = reg_map.get(f, [])
    rows.append({"feature_ID": f, "gene_symbol": symbols.get(f),
                 "selection_frequency": panel.set_index("feature_ID").loc[f, "selection_frequency"],
                 "marker_tissue": t1, "ovr_score": float(sc[i1]), "rank_in_tissue": int(rank_pos[i1][j]),
                 "runner_up_tissue": classes[i2], "runner_up_score": float(sc[i2]),
                 "in_full_data_roundrobin_top_k": f in top_full, "mean_in_marker_tissue": float(means[t1]),
                 "next_highest_tissue": str(others.idxmax()), "next_highest_mean": float(others.max()),
                 "effect_size": float(means[t1] - others.max()),
                 "regulated_in_marker_tissue": t1 in reg_t, "regulated_tissues": ";".join(reg_t),
                 f"r_{qc_col}_in_marker_tissue": r_in, f"max_abs_r_{qc_col}": r_max, "tissue_of_max_r": t_max,
                 "risk_T7_regulated": t1 in reg_t, "risk_qc_correlated": bool(not np.isnan(r_in) and abs(r_in) >= qc_risk)})
s05_ann = pd.DataFrame(rows).sort_values(["marker_tissue", "ovr_score"], ascending=[True, False]).reset_index(drop=True)
s05_ann.to_csv(S05_OUT / "candidate_panel_annotated.csv", index=False)

show = ["gene_symbol", "selection_frequency", "marker_tissue", "rank_in_tissue", "effect_size", "next_highest_tissue",
        "regulated_tissues", f"r_{qc_col}_in_marker_tissue", "risk_T7_regulated", "risk_qc_correlated"]
display(s05_ann[show].round(2))
pub_ann = pd.read_csv(res("05_panels", "TRNSCRPT", "candidate_panel_annotated.csv"), dtype={"feature_ID": str})
num = ["ovr_score", "effect_size", f"r_{qc_col}_in_marker_tissue"]
a, b = s05_ann.set_index("feature_ID"), pub_ann.set_index("feature_ID").reindex(s05_ann["feature_ID"])
print("annotation vs published: marker tissues identical:", (a["marker_tissue"] == b["marker_tissue"]).all(),
      "| max |diff| in", num, "=", f"{np.nanmax(np.abs(a[num].to_numpy() - b[num].to_numpy())):.2e}")
uncovered = [c for c in classes if c not in set(s05_ann["marker_tissue"])]
print(f"tissues with a stable marker: {s05_ann['marker_tissue'].nunique()} of {len(classes)}; without: {', '.join(uncovered)}")
print(f"flags: training-regulated in marker tissue {int(s05_ann['risk_T7_regulated'].sum())}, "
      f"|r| >= {qc_risk} with {qc_col} in marker tissue {int(s05_ann['risk_qc_correlated'].sum())}")
record("s05.n_panel_tissues_covered", s05_ann["marker_tissue"].nunique(), "05")
record("s05.n_flag_regulated", int(s05_ann["risk_T7_regulated"].sum()), "05")
record("s05.n_flag_qc", int(s05_ann["risk_qc_correlated"].sum()), "05")

# %% [markdown]
# **What this shows.** Only part of the k = 20 selection is stable across animal bootstraps: the stable
# genes are textbook tissue markers (adrenal steroidogenesis, kidney uromodulin, lung claudin, the
# hypothalamic neuropeptide, and so on), most with a large effect size. The other chosen genes rotate
# between near-equivalent candidates, which is why accuracy is stable even though the gene list is not.
# Many tissues have no stable marker of their own (printed above). A few stable genes carry a flag:
# training-regulated in their own tissue, or correlated with library mRNA content.
#
# **What it does not show.** The one-vs-rest scores, means and QC correlations in the annotation are
# computed on the full matrix — they describe the panel, they are not an evaluation. A flag is a risk,
# not a demonstrated failure: whether a flagged gene actually breaks under the controls → trained shift
# or in another dataset is tested in sections 8–10. Stability selection here uses one seed and 50
# resamples, so frequencies near the threshold are uncertain.
