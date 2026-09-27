# %% [markdown]
# ## 10. Across species: scoring human GTEx tissues with the rat fingerprint
#
# **Question.** Can the rat tissue fingerprint name *human* tissues? GTEx v8 (17 tissues that have a
# MoTrPAC counterpart, at most 150 donors each) is scored through the same code path as the BodyMap
# (section 9), after an orthology step: only strict 1:1 rat–human ortholog pairs are kept.
#
# **Why it matters.** A human user of the fingerprint would meet three losses at once: genes without a
# 1:1 ortholog, a different species' tissue biology, and a different library protocol and unit (the GTEx
# portal gives TPM; MoTrPAC is analysed as CPM). The section measures what survives.
#
# **Design (from `scripts/13_gtex_transfer.py` and `scripts/14_transfer_representations.py`, flags of
# `make gtex` / `make transfer`).** Panels are re-selected in the ortholog space on all 50 MoTrPAC
# animals; genes are z-scored within species; every GTEx split is grouped on the donor. Mapping: Muscle
# - Skeletal → {SKM-GN, SKM-VL} as a super-class, each brain region to its rat region, Artery - Aorta →
# VENACV (imperfect: an artery for a vein; read as a caveat), all others one-to-one. α = 0.10 (LAC).
#
# **What to look for.** How many genes are lost, and at which step; which tissues transfer and which
# collapse into a wrong rat tissue (read off the confusion matrix, not asserted); whether the rat-
# calibrated conformal sets cover human samples, and what recalibrating on a few donors costs in set
# size; and whether representations that do not depend on the gene's scale (within-sample ranks,
# top-scoring pairs) transfer better than z-scores.
#
# **Mode.** Phase 13 and phase 14 take minutes to tens of minutes (the GTEx matrix alone is a
# multi-gigabyte CSV), so with `RECOMPUTE = False` their tables are loaded from `results/13_gtex/`, `results/14_transfer/`
# and `results/14_transfer_cpm/` (post-fix). The orthology bookkeeping is recomputed live. With
# `RECOMPUTE = True` the copied script logic below re-runs into `_outputs/s10/`.

# %%
section("10 GTEx")

# from scripts/13_gtex_transfer.py (verbatim): GTEx tissue → MoTrPAC tissue(s)
GTEX_TO_RAT = {
    "Whole Blood": {"BLOOD"}, "Muscle - Skeletal": {"SKM-GN", "SKM-VL"}, "Adipose - Subcutaneous": {"WAT-SC"},
    "Heart - Left Ventricle": {"HEART"}, "Liver": {"LIVER"}, "Kidney - Cortex": {"KIDNEY"}, "Lung": {"LUNG"},
    "Brain - Cortex": {"CORTEX"}, "Brain - Hippocampus": {"HIPPOC"}, "Brain - Hypothalamus": {"HYPOTH"},
    "Colon - Transverse": {"COLON"}, "Small Intestine - Terminal Ileum": {"SMLINT"}, "Spleen": {"SPLEEN"},
    "Adrenal Gland": {"ADRNL"}, "Ovary": {"OVARY"}, "Testis": {"TESTES"}, "Artery - Aorta": {"VENACV"},
}
s10_args = SimpleNamespace(alpha=0.10, grid="20,50", cal_frac=0.3, recal_donors="3,5", recal_repeats=20,
                           prefilter=5000, seed=C.SEED, n_splits=5, r_thresh=0.5)
S10_OUT = OUT / "s10"
S10_OUT.mkdir(parents=True, exist_ok=True)
PREFIX = REPO / "backup" / "pipeline_history" / "results_pre_quantile_fix_2026-09-25"   # pre-fix copies (read-only)

pheno = cached("pheno", io.load_pheno)
om_trn = cached("TRNSCRPT_counts", lambda: io.stack_tissues(
    "TRNSCRPT", tissues=None, source="counts", join="inner", pheno=pheno, complete=False, drop_incomplete_samples=None))

# ---- orthology losses, step by step (live). The last two steps are exactly the pipeline's
# transfer.one_to_one_orthologs + transfer.match_gtex_orthologs; the intermediate counts re-trace the
# filters inside one_to_one_orthologs so each loss can be seen. ----
_mot = set(om_trn.X.columns)
_r = pd.read_csv(C.RAW_DIR / "rat_to_human_gene.csv", dtype=str).dropna(subset=["RAT_ENSEMBL_ID"])
_r_m = _r[_r["RAT_ENSEMBL_ID"].isin(_mot)]
_has_h = _r_m[_r_m["HUMAN_ORTHOLOG_ENSEMBL_ID"].notna() | _r_m["HUMAN_ORTHOLOG_SYMBOL"].notna()]
s10_orth = one_to_one_orthologs()
s10_orth = s10_orth[s10_orth["RAT_ENSEMBL_ID"].isin(_mot)]
_gtex_cols = pd.read_csv(C.EXTERNAL_DIR / "gtex_tpm_subset.csv", index_col=0, nrows=0).columns   # header only
s10_pairs, _n_ens, _n_sym = match_gtex_orthologs(s10_orth, _gtex_cols, gtex_symbols())
s10_steps = pd.DataFrame([
    ("MoTrPAC genes (stacked TRNSCRPT)", om_trn.n_features),
    ("… with a row in the RGD rat–human table", _r_m["RAT_ENSEMBL_ID"].nunique()),
    ("… with any human ortholog (id or symbol)", _has_h["RAT_ENSEMBL_ID"].nunique()),
    ("… whose ortholog is strictly 1:1 (both directions)", len(s10_orth)),
    ("… whose human gene is a GTEx column (Ensembl id + symbol fallback, deduplicated)", len(s10_pairs)),
], columns=["step", "genes"])
s10_steps["lost_at_step"] = (-s10_steps["genes"].diff()).fillna(0).astype(int)
print(s10_steps.to_string(index=False))
print(f"\nGTEx matches: {_n_ens} by human Ensembl id, {_n_sym} by symbol fallback; GTEx matrix has {len(_gtex_cols)} genes")
_pub = pd.read_csv(res("13_gtex", "gene_overlap.csv")).iloc[0]
_chk = {"rat_genes_with_1to1_human_ortholog": len(s10_orth), "orthologs_present_in_gtex": len(s10_pairs),
        "matched_by_human_ensembl_id": _n_ens, "matched_by_symbol_fallback": _n_sym,
        "lost_to_orthology": om_trn.n_features - len(s10_pairs)}
print(pd.DataFrame({"recomputed": pd.Series(_chk), "published": _pub[list(_chk)]}).to_string())
_mg = pd.read_csv(C.EXTERNAL_DIR / "gtex_meta.csv", dtype=str)
print(f"\nGTEx metadata: {len(_mg)} samples, {_mg['donor'].nunique()} donors, {_mg['SMTSD'].nunique()} tissues "
      f"(published: {_pub['gtex_samples']} samples, {_pub['gtex_donors']} donors)")
record("s10.orthologs_1to1", len(s10_orth), "10")
record("s10.orthologs_in_gtex", len(s10_pairs), "10")
record("s10.lost_to_orthology", om_trn.n_features - len(s10_pairs), "10")

# %%
# from scripts/13_gtex_transfer.py (verbatim; renamed from TISSUE_TO_ORGAN, which section 9 uses for the BodyMap)
GTEX_TO_RAT_T2O = {t: o for o, ts in GTEX_TO_RAT.items() for t in ts}

# ---- phase 13: scripts/13_gtex_transfer.py::main, copied (trimmed: argparse, report.add_section and printing
# dropped; results written to `out`) ----
def s10_run_phase13(out: Path) -> dict:
    args = s10_args
    grid = [int(k) for k in args.grid.split(",")]
    rng = np.random.default_rng(args.seed)
    om = om_trn
    Xg = pd.read_csv(C.EXTERNAL_DIR / "gtex_tpm_subset.csv", index_col=0)      # samples × human genes (log2 TPM+1)
    mg = pd.read_csv(C.EXTERNAL_DIR / "gtex_meta.csv", dtype=str).set_index("SAMPID").loc[Xg.index]
    mg["organ"] = mg["SMTSD"]
    mg["group_id"] = mg["donor"]
    mg["stage"] = "adult"
    orth = one_to_one_orthologs()
    orth = orth[orth["RAT_ENSEMBL_ID"].isin(om.X.columns)]
    by_ens, n_by_ensembl, used_symbol = match_gtex_orthologs(orth, Xg.columns, gtex_symbols())
    rat_genes = by_ens["RAT_ENSEMBL_ID"].tolist()
    hum_genes = by_ens["HUMAN_ORTHOLOG_ENSEMBL_ID"].tolist()
    ov = pd.DataFrame([{"motrpac_genes": om.n_features, "rat_genes_with_1to1_human_ortholog": int(orth.shape[0]),
                        "orthologs_present_in_gtex": len(rat_genes), "matched_by_human_ensembl_id": n_by_ensembl,
                        "matched_by_symbol_fallback": used_symbol,
                        "lost_to_orthology": om.n_features - len(rat_genes), "gtex_genes": Xg.shape[1], "gtex_samples": len(mg),
                        "gtex_donors": mg["donor"].nunique()}])
    sym = io.map_to_gene_symbols(om.X.columns)
    Zm = zscore(om.X[rat_genes])
    Zt = zscore(Xg[hum_genes])
    L_src = om.X[rat_genes].to_numpy(dtype=float)
    L_tgt = Xg[hum_genes].to_numpy(dtype=float)
    y = om.meta["tissue"].astype(str).to_numpy()
    g = om.groups()
    pm = PanelModels(L_src, y, rat_genes, grid, args.prefilter)
    Lt_pre, Zt_pre = pm.target_matrices(L_tgt)
    model_names = [f"k{k}" for k in grid] + ["full"]
    pm_all = PanelModels(om.X.to_numpy(dtype=float), y, list(om.X.columns), grid, args.prefilter, full=False)
    surv = []
    for k in grid:
        full_panel = pm_all.panel_genes(k)
        surv.append({"panel": f"k={k} (all-gene fit)", "n_genes": k, "with_1to1_ortholog_in_gtex": sum(gg in set(rat_genes) for gg in full_panel),
                     "lost": ";".join(str(sym.get(gg, gg)) for gg in full_panel if gg not in set(rat_genes))})
    core_p = C.RESULTS_DIR / "05_panels" / "TRNSCRPT" / "candidate_panel_annotated.csv"
    if core_p.exists():
        core = pd.read_csv(core_p)
        surv.append({"panel": "T5 stable core", "n_genes": len(core), "with_1to1_ortholog_in_gtex": int(core["feature_ID"].isin(rat_genes).sum()),
                     "lost": ";".join(core.loc[~core["feature_ID"].isin(rat_genes), "gene_symbol"].astype(str))})
    surv = pd.DataFrame(surv)
    acc_rows, conf = [], {}
    for name in model_names:
        pred, ok, _, _ = score_block(pm, name, Lt_pre, Zt_pre, mg, GTEX_TO_RAT)
        sub = mg.assign(pred=pred, correct=ok)
        for organ, d in sub.groupby("organ"):
            acc_rows.append({"model": name, "gtex_tissue": organ, "rat_classes": "/".join(sorted(GTEX_TO_RAT[organ])), "n": len(d),
                             "n_donors": d["donor"].nunique(), "accuracy": float(d["correct"].mean()),
                             "top_prediction": d["pred"].value_counts().index[0], "top_prediction_frac": float(d["pred"].value_counts().iloc[0] / len(d))})
        conf[name] = pd.crosstab(sub["organ"], sub["pred"]).reindex(columns=pm.classes, fill_value=0)
    acc = pd.DataFrame(acc_rows)
    overall = acc.groupby("model").apply(lambda d: np.average(d["accuracy"], weights=d["n"]), include_groups=False)
    macro = acc.groupby("model")["accuracy"].mean()
    acc_overall = pd.DataFrame({"model": overall.index, "accuracy_sample_weighted": overall.to_numpy(), "accuracy_macro_over_tissues": macro.loc[overall.index].to_numpy(),
                                "n_samples": int(len(mg)), "n_tissues": int(mg["organ"].nunique())})
    uniq = np.unique(g)
    rng.shuffle(uniq)
    cal_animals = set(uniq[:max(1, int(round(args.cal_frac * len(uniq))))])
    is_cal = np.array([a in cal_animals for a in g])
    fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
    assert_no_group_leak(om.meta, fit_idx, cal_idx)
    pm_c = PanelModels(L_src[fit_idx], y[fit_idx], rat_genes, grid, args.prefilter)
    classes = pm_c.classes
    Lm_c, Zm_c = pm_c.source_matrices(L_src[cal_idx])
    calib = calibrate_models(pm_c, Lm_c, Zm_c, y[cal_idx], classes, args.alpha, model_names)
    Lt_c, Zt_c = pm_c.target_matrices(L_tgt)
    recal_ns = [int(v) for v in args.recal_donors.split(",") if v.strip()]
    conf_df, ood, recal = conformal_transfer(calib, classes, Lt_c, Zt_c, mg, GTEX_TO_RAT, args.alpha, "stage", "adult",
                                             recal_ns, args.recal_repeats, rng, model_names, [])
    pt = per_organ_coverage(calib, classes, Lt_c, Zt_c, mg, GTEX_TO_RAT, model_names).rename(columns={"organ": "gtex_tissue", "coverage": "coverage_rat_cal"})
    Zt_df = pd.DataFrame(Zt, index=mg.index, columns=hum_genes)
    meta_g = pd.DataFrame({"pid": mg["donor"].to_numpy(), "bid": mg.index, "sex": "unknown", "group": "adult",
                           "tissue": mg["organ"].to_numpy(), "assay": "GTEX"}, index=mg.index)
    om_g = io.OmicsMatrix(Zt_df, meta_g, pd.DataFrame(index=hum_genes), "GTEX", name="gtex")
    curve, sel, preds = models.panel_curve(om_g, label="tissue", grid=[20], kind="logreg_l2", n_splits=5, prefilter=args.prefilter,
                                           quick=True, seed=args.seed, verbose=False, return_predictions=True)
    native_acc = float((preds["y_true"] == preds["y_pred"]).mean())
    imported = {n: float(overall[n]) for n in model_names}
    native = pd.DataFrame([{"native_gtex_panel_k20_accuracy": native_acc, "native_panel_balanced_accuracy": float(curve["balanced_accuracy"].mean()),
                            **{f"imported_{n}_accuracy": v for n, v in imported.items()},
                            "gap_native_minus_imported_k20": native_acc - imported["k20"], "n_donors": mg["donor"].nunique(),
                            "native_panel_genes_fold0": ";".join(sel[sel["fold"] == 0]["feature_ID"].tolist())}])
    flags = {}
    for fn in ("stability_k20_annotated.csv", "candidate_panel_annotated.csv"):
        p = C.RESULTS_DIR / "05_panels" / "TRNSCRPT" / fn
        if p.exists():
            for _, r in pd.read_csv(p).iterrows():
                flags.setdefault(r["feature_ID"], {"risk_T7_regulated": bool(r.get("risk_T7_regulated", False)),
                                                  "risk_qc_correlated": bool(r.get("risk_qc_correlated", False))})
    Zm_df = pd.DataFrame(Zm, index=om.X.index, columns=rat_genes)
    rat_to_hum = dict(zip(rat_genes, hum_genes))
    Zt_rat_cols = pd.DataFrame(Zt, index=mg.index, columns=rat_genes)
    gc = gene_check(pm, pm.panel_genes(20 if 20 in grid else grid[0]), sym, Zm_df, y, Zt_rat_cols, mg, GTEX_TO_RAT_T2O, flags)
    gc["human_gene"] = gc["feature_ID"].map(rat_to_hum)
    tables = {"gene_overlap": ov, "panel_survival": surv, "accuracy_by_tissue": acc, "accuracy_overall": acc_overall,
              "conformal_transfer": conf_df, "recalibration": recal, "coverage_by_tissue": pt, "native_panel": native,
              "panel_gene_check": gc, **{f"confusion_{n}": conf[n].reset_index() for n in model_names}}
    for name, t in tables.items():
        t.to_csv(out / f"{name}.csv", index=False)
    return tables


S10_T13 = ["panel_survival", "accuracy_by_tissue", "accuracy_overall", "conformal_transfer", "recalibration",
           "coverage_by_tissue", "native_panel", "panel_gene_check", "confusion_k20", "confusion_k50", "confusion_full"]
if RECOMPUTE:
    t0 = time.perf_counter()
    (S10_OUT / "13_gtex").mkdir(exist_ok=True)
    s10_t = s10_run_phase13(S10_OUT / "13_gtex")
    print(f"phase 13 recomputed in {time.perf_counter() - t0:.0f} s")
else:
    s10_t = {n: pd.read_csv(res("13_gtex", f"{n}.csv")) for n in S10_T13}
    print("phase 13 tables loaded from", RES / "13_gtex")
s10_pre = {n: pd.read_csv(PREFIX / "13_gtex" / f"{n}.csv") for n in S10_T13}

# %% [markdown]
# **What the orthology table shows.** Each row is a filter; `lost_at_step` is how many MoTrPAC genes it
# removes. The recomputed counts are set beside the pipeline's `gene_overlap.csv`. The panel-survival
# table below then asks the sharper question: of the genes the *all-gene* MoTrPAC panels would pick,
# how many survive orthology (the ortholog-space panels used for GTEx are re-selected, so a lost marker
# is replaced by the next-best gene, not left as a hole).
#
# ### 10.1 Per-tissue accuracy at k = 20, k = 50 and the full model

# %%
print(s10_t["panel_survival"].to_string(index=False))
s10_acc = s10_t["accuracy_by_tissue"].pivot(index="gtex_tissue", columns="model", values="accuracy")[["k20", "k50", "full"]]
_n = s10_t["accuracy_by_tissue"].drop_duplicates("gtex_tissue").set_index("gtex_tissue")["n"]
print("\nPer-tissue accuracy (super-class scoring for muscle):")
print(s10_acc.assign(n=_n).round(2).to_string())
# overall: sample-weighted (every GTEx sample counts once — the convention used throughout) and macro over tissues,
# recomputed here from the per-tissue table and set beside the pipeline's accuracy_overall.csv
_ov = pd.DataFrame({"sample_weighted": s10_acc.mul(_n, axis=0).sum() / _n.sum(), "macro": s10_acc.mean()})
_pov = s10_t["accuracy_overall"].set_index("model")
_ov["published_sample_weighted"] = _pov["accuracy_sample_weighted"]
_ov["published_macro"] = _pov["accuracy_macro_over_tissues"]
print("\nOverall:")
print(_ov.round(3).to_string())
for m in ("k20", "k50", "full"):
    record(f"s10.acc_{m}", _ov.loc[m, "sample_weighted"], "10")
record("s10.acc_heart_k20", s10_acc.loc["Heart - Left Ventricle", "k20"], "10")
record("s10.acc_ovary_k20", s10_acc.loc["Ovary", "k20"], "10")

# %% [markdown]
# ### 10.2 Where the misses go: the confusion matrix, derived
#
# The pipeline writes one confusion table per model (rows = GTEx tissue, columns = predicted rat
# tissue). The cell below *derives* from it, by code: the row-normalised matrix (shown as a heat
# map), a check that its super-class diagonal reproduces the per-tissue accuracy above exactly, and for
# every tissue the rat tissue most of its errors go to.

# %%
# notebook-only helper (not in the pipeline): presentation of tables the pipeline wrote
def s10_confusion_summary(cf: pd.DataFrame) -> pd.DataFrame:
    """Per GTEx tissue: n, correct (super-class diagonal), accuracy, the main wrong destination and its share."""
    rows = []
    for tissue, r in cf.iterrows():
        ok_cols = [c for c in GTEX_TO_RAT[tissue] if c in cf.columns]
        wrong = r.drop(ok_cols)
        rows.append({"gtex_tissue": tissue, "n": int(r.sum()), "correct": int(r[ok_cols].sum()),
                     "accuracy": r[ok_cols].sum() / r.sum(), "main_wrong_call": wrong.idxmax() if wrong.sum() else "—",
                     "main_wrong_frac": wrong.max() / r.sum()})
    return pd.DataFrame(rows).set_index("gtex_tissue")


s10_cf = {m: s10_t[f"confusion_{m}"].set_index("organ") for m in ("k20", "k50", "full")}
s10_cs = s10_confusion_summary(s10_cf["k20"])
print("k = 20, derived from confusion_k20.csv:")
print(s10_cs.round(2).to_string())
_d = max(float((s10_confusion_summary(s10_cf[m])["accuracy"] - s10_acc[m]).abs().max()) for m in s10_cf)
print("max |confusion-derived − published per-tissue accuracy| over the three models:", _d)
_cfn = s10_cf["k20"].div(s10_cf["k20"].sum(axis=1), axis=0)
_cfn = _cfn.loc[:, (_cfn > 0).any(axis=0)]
fig, ax = plt.subplots(figsize=(10, 6))
im = ax.imshow(_cfn.to_numpy(), cmap="viridis", vmin=0, vmax=1, aspect="auto")
ax.set_xticks(range(_cfn.shape[1]), _cfn.columns, rotation=60, ha="right", fontsize=8)
ax.set_yticks(range(_cfn.shape[0]), _cfn.index, fontsize=8)
ax.set_xlabel("predicted rat tissue"); ax.set_title("GTEx → rat, k = 20: row-normalised confusion")
fig.colorbar(im, ax=ax, label="fraction of the GTEx tissue's samples")
plt.tight_layout(); plt.show()
record("s10.heart_to_skm_frac_k20", float(s10_cf["k20"].loc["Heart - Left Ventricle", ["SKM-GN", "SKM-VL"]].sum()
                                           / s10_cf["k20"].loc["Heart - Left Ventricle"].sum()), "10")
record("s10.ovary_to_heart_frac_k20", float(s10_cf["k20"].loc["Ovary", "HEART"] / s10_cf["k20"].loc["Ovary"].sum()), "10")

# %% [markdown]
# **What this shows.** The fingerprint crosses species unevenly: tissues with a strong, conserved
# marker program (kidney, liver, lung, spleen, skeletal muscle) transfer at every panel size, while at
# k = 20 some tissues collapse almost entirely into one wrong rat tissue — the `main_wrong_call` column
# names where (e.g. the human left ventricle into the rat skeletal-muscle classes, the human ovary into
# rat heart). Larger models recover most of these; the full model is the best overall. The
# derived-vs-published check confirms the confusion tables and the accuracy table describe the same
# predictions.
#
# **What it does not show.** *Why* a tissue collapses. A 20-gene panel has one or two markers per
# tissue, so a single marker that is not conserved (or not expressed in the human sample) can move a
# whole tissue; the panel-gene check in `results/13_gtex/panel_gene_check.csv` lists which markers fail
# in GTEx, but which failure causes which collapse is not tested here. Aorta → VENACV is an
# artery-for-vein mapping and its accuracy is not a fair test of either tissue.
#
# ### 10.3 Conformal coverage with rat calibration, and recalibration on GTEx donors
#
# As in 9.3: models refit on the MoTrPAC fit animals, calibrated on held-out MoTrPAC animals, sets
# formed on every GTEx sample; then the threshold is recalibrated on 3 or 5 GTEx donors (pooled
# tissues, 20 draws) and tested on the remaining donors. Every number is printed beside its pre-fix
# value.

# %%
_key = ["model", "conformal"]
_ct = s10_t["conformal_transfer"].merge(s10_pre["conformal_transfer"], on=["stage"] + _key, suffixes=("", "_prefix"))
s10_cov = _ct[_key + ["coverage_mapped", "coverage_mapped_prefix", "frac_empty_mapped", "frac_empty_mapped_prefix",
                      "avg_set_size_mapped", "avg_set_size_mapped_prefix"]]
print("All GTEx samples, rat calibration, α = 0.10:")
print(s10_cov.round(3).to_string(index=False))
_pt = s10_t["coverage_by_tissue"]
_ptw = _pt[(_pt["model"] == "k20")].pivot(index="gtex_tissue", columns="conformal", values="coverage_rat_cal")
print("\nPer GTEx tissue, k = 20, rat calibration:")
print(_ptw[["marginal", "mondrian", "floored"]].assign(accuracy_k20=s10_acc["k20"]).round(2).to_string())
_rc = s10_t["recalibration"].merge(s10_pre["recalibration"], on=["model", "n_recal"], suffixes=("", "_prefix"))
s10_recal = _rc[["model", "n_recal", "draws", "coverage_recalibrated", "coverage_recalibrated_prefix",
                 "coverage_source_cal_same_test", "coverage_source_cal_same_test_prefix",
                 "set_size_recalibrated", "set_size_recalibrated_prefix", "frac_empty_recalibrated"]]
print("\nRecalibration on GTEx donors:")
print(s10_recal.round(3).to_string(index=False))
_c = s10_cov.set_index(_key)
for m in ("k20", "k50", "full"):
    for v in ("marginal", "mondrian", "floored"):
        record(f"s10.cov_{m}_{v}", _c.loc[(m, v), "coverage_mapped"], "10")
_r = s10_recal.set_index(["model", "n_recal"])
for m in ("k20", "k50", "full"):
    for n in (3, 5):
        record(f"s10.recal_cov_{m}_n{n}", _r.loc[(m, n), "coverage_recalibrated"], "10")
        record(f"s10.recal_size_{m}_n{n}", _r.loc[(m, n), "set_size_recalibrated"], "10")

# %% [markdown]
# Two of the post-fix numbers above need a mechanism check before they are read: the Mondrian / floored
# sets are never empty, and recalibration on 3 donors reaches high coverage with very large sets. The
# next cell derives both from the script's own seeded draws: (i) which tissues have too few rat
# calibration vials for a finite per-tissue threshold (phase 13 uses the same seed and the same
# calibration split as phase 12); (ii) for each recalibration draw, how many GTEx samples the chosen
# donors contribute, and whether that is enough for a finite textbook threshold, ⌈(n+1)(1 − α)⌉ ≤ n.

# %%
# (i) rat calibration animals, as in scripts/13_gtex_transfer.py: rng = default_rng(seed); rng.shuffle(unique pids)
_g = om_trn.groups()
_uniq = np.unique(_g)
_rng = np.random.default_rng(s10_args.seed)
_rng.shuffle(_uniq)
_cal = set(_uniq[:max(1, int(round(s10_args.cal_frac * len(_uniq))))])
_n_cal = om_trn.meta.loc[[a in _cal for a in _g], "tissue"].value_counts()
_inf = sorted(t for t, n in _n_cal.items() if np.ceil((n + 1) * (1 - s10_args.alpha) - 1e-9) > n)
print("rat tissues whose Mondrian threshold is infinite (in every Mondrian / floored set):", _inf,
      {t: int(_n_cal[t]) for t in _inf}, "calibration vials")
# (ii) the recalibration draws of transfer.conformal_transfer, replayed with the same generator (after the shuffle
# above): per model, per n_recal, 20 draws of donors via rng.choice(ids, n, replace=False); every GTEx tissue is
# mapped, so the calibration set is every sample of the chosen donors
_mg = pd.read_csv(C.EXTERNAL_DIR / "gtex_meta.csv", dtype=str)
_per_donor = _mg["donor"].value_counts()
_ids = sorted(_mg["donor"].unique())
_rows = []
for m in ("k20", "k50", "full"):
    for n_recal in (3, 5):
        for _ in range(s10_args.recal_repeats):
            chosen = tuple(_rng.choice(_ids, size=n_recal, replace=False))
            n = int(_per_donor[list(chosen)].sum())
            _rows.append({"model": m, "n_recal": n_recal, "n_cal_samples": n,
                          "threshold_infinite": bool(np.ceil((n + 1) * (1 - s10_args.alpha) - 1e-9) > n)})
s10_draws = pd.DataFrame(_rows)
print(f"\nGTEx samples per donor: median {int(_per_donor.median())}, range {_per_donor.min()}–{_per_donor.max()}")
s10_draw_sum = s10_draws.groupby(["model", "n_recal"]).agg(median_cal_samples=("n_cal_samples", "median"),
                                                           frac_draws_infinite=("threshold_infinite", "mean"))
s10_draw_sum["set_size_recalibrated"] = s10_recal.set_index(["model", "n_recal"])["set_size_recalibrated"]
print(s10_draw_sum.round(3).to_string())
record("s10.frac_recal_draws_infinite_k20_n3", s10_draw_sum.loc[("k20", 3), "frac_draws_infinite"], "10")

# %% [markdown]
# **What this shows.** Across species the rat-calibrated guarantee fails badly: marginal coverage on
# GTEx is far below the 90 % target for every model and lowest for the *full* model, the most accurate
# one — the shortfall is empty sets: on human data the classifiers' top probabilities rarely reach the
# threshold set on rat animals, so accuracy (argmax) and coverage (thresholded probability) come
# apart. Per tissue, coverage can be far below accuracy for tissues the argmax gets right.
#
# The Mondrian and floored rows look better after the fix only through the always-included tissues of
# the cell above: GTEx ovary and testis are "covered" by every set (per-tissue table) although the
# argmax never or rarely names them, and every other tissue pays with larger sets.
#
# Recalibration on donors: with 3 donors a large fraction of the draws (printed above) contribute too
# few samples for a finite threshold, so the set is every tissue and "coverage" is trivial — read the huge set size, not the
# coverage. With 5 donors more draws have a finite threshold, and coverage near the target comes with
# sets of several tissues for the panels — a weak answer. Before the fix the same small calibration
# sets returned their largest score instead of an infinite threshold, which is why the pre-fix
# 3-donor sets were smaller and their coverage lower.
#
# **What it does not show.** The replayed draws reproduce the calibration *sizes* only; the per-draw
# coverage is not recomputed here (that needs the fitted models — `RECOMPUTE = True`). A finite
# threshold from so few samples is still a very noisy one.
#
# ### 10.4 Native GTEx panel versus the imported rat panels
#
# The ceiling, as in 9.2: a 20-gene panel selected and cross-validated on GTEx itself (donor-grouped
# 5-fold CV, all 17 tissues, same samples as the imported scores).

# %%
_nat = s10_t["native_panel"].iloc[0]
print(pd.Series({"native GTEx k20 (accuracy)": _nat["native_gtex_panel_k20_accuracy"],
                 "native GTEx k20 (balanced accuracy)": _nat["native_panel_balanced_accuracy"],
                 "imported rat k20": _nat["imported_k20_accuracy"], "imported rat k50": _nat["imported_k50_accuracy"],
                 "imported rat full": _nat["imported_full_accuracy"],
                 "gap native − imported k20": _nat["gap_native_minus_imported_k20"]}).round(3).to_string())
record("s10.native_acc", _nat["native_gtex_panel_k20_accuracy"], "10")
_gc = s10_t["panel_gene_check"]
_bad = _gc[(_gc["fails_in_target"] == True) | (_gc["weakened"] == True)]
print(f"\nk = 20 rat panel genes that fail or weaken in GTEx ({len(_bad)} of {len(_gc)}):")
print(_bad[["gene_symbol", "marker_tissue", "source_effect_z", "target_effect_z", "target_top_organ", "fails_in_target",
            "risk_T7_regulated", "risk_qc_correlated"]].round(2).to_string(index=False))

# %% [markdown]
# **What this shows.** Human tissues are easy to tell apart with human markers — the native panel is
# near the top of the scale — so the imported panels' shortfall is a transfer loss, not a hard task. A
# large share of the rat panel's markers fail or weaken in GTEx (the count is printed above), several of
# them markers the MoTrPAC-side risk flags (training-regulated, correlated with library mRNA fraction)
# had already marked.
#
# **What it does not show.** The native panel is a CV estimate on GTEx; the imported number is an
# external test — different kinds of number, set side by side as a ceiling, not as a fair contest.
#
# ### 10.5 Do scale-free representations transfer better? (phase 14)
#
# Three representations of a k-gene panel: per-gene **z-scores** within each dataset (used above),
# within-sample **ranks** of the panel genes, and **top-scoring pairs** (for every pair of panel genes,
# which one is higher in the sample). Ranks and pairs do not depend on the gene's scale, so in principle
# they should survive a change of units or library protocol. Two selectors: `standard` round robin and
# `transfer_aware` (a tissue may not pick a gene that is training-regulated in it, or correlated with the
# library mRNA fraction at |r| > 0.5, computed on MoTrPAC only). The cost side is MoTrPAC animal-grouped
# CV; the benefit side is the two external targets.

# %%
# ---- phase 14: scripts/14_transfer_representations.py, copied (trimmed: argparse, report.add_section and printing
# dropped; `targets` / `gtex_matrix` / `skip_cv` arguments replace --targets / --gtex-matrix / --skip-cv) ----
from sklearn.metrics import balanced_accuracy_score   # the script's own import (phase-14 CV)

# from scripts/14_transfer_representations.py (verbatim): the organ maps and names
BODYMAP_MAP = {"Adrenal": {"ADRNL"}, "Brain": {"CORTEX", "HIPPOC", "HYPOTH"}, "Heart": {"HEART"}, "Kidney": {"KIDNEY"},
               "Liver": {"LIVER"}, "Lung": {"LUNG"}, "Muscle": {"SKM-GN", "SKM-VL"}, "Spleen": {"SPLEEN"},
               "Testes": {"TESTES"}, "Thymus": None, "Uterus": None}
GTEX_MAP = dict(GTEX_TO_RAT)   # identical to phase 13's GTEX_TO_RAT in the script
SELECTORS = ("standard", "transfer_aware")
TARGET_NAMES = {"bodymap": "BodyMap adults", "gtex": "GTEx"}


# from scripts/14_transfer_representations.py::load_bodymap (verbatim; renamed — this is NOT section 9's
# loader: no sex relabelling, group_id built from the raw sex code, returns log2 CPM on the shared genes)
def s10_load_bodymap_p14(shared_with):
    X = pd.read_csv(C.EXTERNAL_DIR / "bodymap_counts.csv", index_col=0)
    X.index = X.index.astype(str)
    meta = pd.read_csv(C.EXTERNAL_DIR / "bodymap_meta.csv", dtype=str).set_index("sample").loc[X.columns]
    meta["stage_weeks"] = meta["stage_weeks"].astype(int)
    meta["replicate_index"] = meta["replicate"].astype(str).str.rsplit("_", n=1).str[-1]
    meta["group_id"] = meta["sex"].astype(str).str[0] + "_" + meta["stage_weeks"].astype(str) + "_" + meta["replicate_index"]
    lb = io.log_cpm(X.T, log=True)
    shared = [g for g in shared_with if g in set(lb.columns)]
    return lb[shared].to_numpy(dtype=float), shared, meta


# from scripts/14_transfer_representations.py::load_gtex (verbatim)
def load_gtex(motrpac_genes, matrix="tpm"):
    """GTEx samples × 1:1-ortholog genes (rat ids), on log2 TPM or log2 CPM; matching as in phase 13."""
    path = C.EXTERNAL_DIR / ("gtex_cpm_subset.csv" if matrix == "cpm" else "gtex_tpm_subset.csv")
    if not path.exists():
        raise SystemExit(f"{path} is missing — run scripts/11_gtex_prepare.py" + (" --reads <gene_reads.gct.gz>" if matrix == "cpm" else ""))
    Xg = pd.read_csv(path, index_col=0)
    mg = pd.read_csv(C.EXTERNAL_DIR / "gtex_meta.csv", dtype=str).set_index("SAMPID").loc[Xg.index]
    mg["organ"] = mg["SMTSD"]
    mg["group_id"] = mg["donor"]
    mg["stage"] = "adult"
    orth = one_to_one_orthologs()
    orth = orth[orth["RAT_ENSEMBL_ID"].isin(motrpac_genes)]
    by, _, _ = match_gtex_orthologs(orth, Xg.columns, gtex_symbols())
    return Xg[by["HUMAN_ORTHOLOG_ENSEMBL_ID"].tolist()].to_numpy(dtype=float), by["RAT_ENSEMBL_ID"].tolist(), mg


def s10_run_phase14(out: Path, targets=("bodymap", "gtex"), gtex_matrix="tpm", skip_cv=False) -> dict:
    args = s10_args
    out.mkdir(parents=True, exist_ok=True)
    grid = [int(k) for k in args.grid.split(",")]
    recal_ns = [3, 5]
    n_splits = args.n_splits
    reps = list(REPRESENTATIONS)
    selectors = list(SELECTORS)
    target_keys = list(targets)
    om = om_trn
    y = om.meta["tissue"].astype(str).to_numpy()
    g = om.groups()
    classes = sorted(np.unique(y))
    sym = io.map_to_gene_symbols(om.X.columns)
    tm = pd.read_csv(C.META_DIR / "TRNSCRPT.csv", dtype=str, low_memory=False).drop_duplicates("viallabel").set_index("viallabel")
    qc = pd.to_numeric(tm["pct_mrna"], errors="coerce").reindex(om.X.index).to_numpy()
    reg = pd.read_csv(C.RAW_DIR / "training_regulated_features.csv", dtype=str, usecols=["feature_ID", "assay", "tissue"])
    reg = reg[reg["assay"] == "TRNSCRPT"]
    reg_map = reg.groupby("tissue")["feature_ID"].agg(set).to_dict()
    all_genes = list(om.X.columns)
    L_all = om.X.to_numpy(dtype=float)
    tables = {}
    if not skip_cv:
        cv_rows = []
        for fold, (tr, te) in enumerate(grouped_kfold(om.meta, "tissue", n_splits, args.seed)):
            assert_no_group_leak(om.meta, tr, te)
            for selector in selectors:
                ex = None if selector == "standard" else exclusion_mask(L_all[tr], y[tr], all_genes, classes, qc=qc[tr], regulated=reg_map, r_thresh=args.r_thresh)
                for rep in reps:
                    pm = PanelModels(L_all[tr], y[tr], all_genes, grid, args.prefilter, representation=rep, exclude=ex, full=False)
                    Lp_te, Zp_te = pm.source_matrices(L_all[te])
                    for k in grid:
                        p, cls = pm.proba(f"k{k}", Lp_te, Zp_te)
                        pred = np.asarray(cls)[p.argmax(axis=1)]
                        cv_rows.append({"selector": selector, "representation": rep, "k": k, "fold": fold,
                                        "balanced_accuracy": balanced_accuracy_score(y[te], pred), "n_test_animals": int(om.meta.iloc[te]["pid"].nunique())})
        tables["motrpac_cv_cost"] = pd.DataFrame(cv_rows)
    tg = {}
    if "bodymap" in target_keys:
        Lb, shared_b, mb = s10_load_bodymap_p14(all_genes)
        tg[TARGET_NAMES["bodymap"]] = (Lb, shared_b, mb, BODYMAP_MAP, "stage_weeks", 21, ["Thymus", "Uterus"])
    if "gtex" in target_keys:
        Lg, rat_g, mg = load_gtex(set(all_genes), gtex_matrix)
        tg[TARGET_NAMES["gtex"]] = (Lg, rat_g, mg, GTEX_MAP, "stage", "adult", [])
    acc_rows, sum_rows, recal_rows = [], [], []
    for tname, (L_tgt, genes_t, mt, omap, stage_col, primary, ood) in tg.items():
        L_src = om.X[genes_t].to_numpy(dtype=float)
        primary_m = (mt[stage_col] == primary).to_numpy()
        uniq = np.unique(g)
        rng_split = np.random.default_rng(args.seed)
        rng_split.shuffle(uniq)
        cal_animals = set(uniq[:max(1, int(round(args.cal_frac * len(uniq))))])
        is_cal = np.array([a in cal_animals for a in g])
        fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
        assert_no_group_leak(om.meta, fit_idx, cal_idx)
        ex_ta_all = exclusion_mask(L_src, y, genes_t, classes, qc=qc, regulated=reg_map, r_thresh=args.r_thresh)
        for selector in selectors:
            ex_all = None if selector == "standard" else ex_ta_all
            ex_fit = None if selector == "standard" else exclusion_mask(L_src[fit_idx], y[fit_idx], genes_t, classes, qc=qc[fit_idx], regulated=reg_map, r_thresh=args.r_thresh)
            for rep in reps:
                pm = PanelModels(L_src, y, genes_t, grid, args.prefilter, representation=rep, exclude=ex_all, full=False)
                Lt_pre, Zt_pre = pm.target_matrices(L_tgt)
                for k in grid:
                    pred, ok, _, _ = score_block(pm, f"k{k}", Lt_pre[primary_m], Zt_pre[primary_m], mt[primary_m], omap)
                    sub = mt[primary_m].assign(correct=ok, pred=pred)
                    for organ, d in sub.groupby("organ"):
                        if not omap.get(organ):
                            continue
                        acc_rows.append({"target": tname, "selector": selector, "representation": rep, "k": k, "organ": organ,
                                         "n": len(d), "accuracy": float(d["correct"].mean()), "top_prediction": d["pred"].value_counts().index[0]})
                    mapped = np.array([bool(omap.get(o)) for o in sub["organ"]])
                    per_organ = [float(d["correct"].mean()) for o, d in sub.groupby("organ") if omap.get(o)]
                    sum_rows.append({"target": tname, "selector": selector, "representation": rep, "k": k,
                                     "accuracy": float(ok[mapped].mean()), "accuracy_macro_over_tissues": float(np.mean(per_organ)),
                                     "n": int(mapped.sum())})
                pm_c = PanelModels(L_src[fit_idx], y[fit_idx], genes_t, grid, args.prefilter, representation=rep, exclude=ex_fit, full=False)
                Lm_c, Zm_c = pm_c.source_matrices(L_src[cal_idx])
                names = [f"k{k}" for k in grid]
                calib = calibrate_models(pm_c, Lm_c, Zm_c, y[cal_idx], pm_c.classes, args.alpha, names)
                Lt_c, Zt_c = pm_c.target_matrices(L_tgt)
                conf_df, _, recal = conformal_transfer(calib, pm_c.classes, Lt_c, Zt_c, mt, omap, args.alpha, stage_col, primary,
                                                       recal_ns, args.recal_repeats, np.random.default_rng(args.seed), names, ood)
                src_cov = conf_df[(conf_df[stage_col] == primary) & (conf_df["conformal"] == "marginal")].set_index("model")
                for _, r in recal.iterrows():
                    recal_rows.append({"target": tname, "selector": selector, "representation": rep, "k": int(r["model"][1:]),
                                       "n_recal": int(r["n_recal"]), "coverage_source_cal": float(src_cov.loc[r["model"], "coverage_mapped"]),
                                       "frac_empty_source_cal": float(src_cov.loc[r["model"], "frac_empty_mapped"]),
                                       "coverage_recalibrated": r["coverage_recalibrated"], "set_size_recalibrated": r["set_size_recalibrated"],
                                       "frac_empty_recalibrated": r["frac_empty_recalibrated"]})
    tables.update({"target_accuracy": pd.DataFrame(acc_rows), "target_summary": pd.DataFrame(sum_rows),
                   "recalibration": pd.DataFrame(recal_rows)})
    for name, t in tables.items():
        t.to_csv(out / f"{name}.csv", index=False)
    return tables

S10_T14 = ["target_accuracy", "target_summary", "recalibration"]
if RECOMPUTE:
    s10_p14 = s10_run_phase14(S10_OUT / "14_transfer")                                    # `make transfer`, first command
    s10_p14cpm = s10_run_phase14(S10_OUT / "14_transfer_cpm", targets=("gtex",), gtex_matrix="cpm", skip_cv=True)  # second
else:
    s10_p14 = {n: pd.read_csv(res("14_transfer", f"{n}.csv")) for n in S10_T14 + ["motrpac_cv_cost"]}
    s10_p14cpm = {n: pd.read_csv(res("14_transfer_cpm", f"{n}.csv")) for n in S10_T14}
    print("phase 14 tables loaded from", RES / "14_transfer", "and", RES / "14_transfer_cpm")
s10_p14_pre = {n: pd.read_csv(PREFIX / "14_transfer" / f"{n}.csv") for n in S10_T14}

# cost: MoTrPAC CV balanced accuracy (mean ± sd over the 5 folds)
_cv = s10_p14["motrpac_cv_cost"].groupby(["representation", "k", "selector"])["balanced_accuracy"].agg(["mean", "std"]).unstack("selector")
print("MoTrPAC animal-grouped CV balanced accuracy (the cost):")
print(_cv.round(3).to_string())
# benefit: accuracy on the external targets (sample-weighted over mapped samples; BodyMap adults; GTEx on log2 TPM)
s10_rep = s10_p14["target_summary"].pivot_table(index=["target", "k", "representation"], columns="selector", values="accuracy")
print("\nTarget accuracy:")
print(s10_rep.round(3).to_string())
_rc14 = s10_p14["recalibration"].merge(s10_p14_pre["recalibration"], on=["target", "selector", "representation", "k", "n_recal"],
                                      suffixes=("", "_prefix"))
_rc14 = _rc14[(_rc14["n_recal"] == 3) & (_rc14["selector"] == "standard")]
print("\nCoverage, standard selector, source-calibrated vs recalibrated on 3 target individuals (with pre-fix values):")
print(_rc14[["target", "k", "representation", "coverage_source_cal", "coverage_source_cal_prefix", "coverage_recalibrated",
             "coverage_recalibrated_prefix", "set_size_recalibrated", "set_size_recalibrated_prefix"]].round(3).to_string(index=False))
for (t, k, r), v in s10_rep["standard"].items():
    if t == "GTEx":
        record(f"s10.p14_gtex_{r}_k{k}_acc", v, "10")
_rcs = _rc14.set_index(["target", "k", "representation"])
for (t, k, r) in _rcs.index:
    tag = "gtex" if t == "GTEx" else "bodymap"
    record(f"s10.p14_{tag}_{r}_k{k}_srccov", _rcs.loc[(t, k, r), "coverage_source_cal"], "10")
    record(f"s10.p14_{tag}_{r}_k{k}_recalcov_n3", _rcs.loc[(t, k, r), "coverage_recalibrated"], "10")

# %% [markdown]
# **What this shows.** The cost side is small: in MoTrPAC cross-validation every representation and
# both selectors are close to one another (the transfer-aware selector even helps the rank panels), so
# the choice is decided on the targets. Within species (BodyMap adults) all three representations
# transfer, z-scores best at k = 20. Across species z-scores are clearly best, and ranks and pairs fall
# far behind at k = 20 — the opposite of the "scale-free transfers better" hope. The transfer-aware
# selector does not help the z-score panels on either target. The coverage table repeats the pattern
# of 10.3: source-calibrated coverage far below target on GTEx for every representation; recalibrated
# coverage on 3 target individuals looks restored, but on GTEx it comes with very large sets, for the
# same small-calibration-set reason shown in 10.3 (phase 14 draws its own donors, so the replay above
# does not apply draw-for-draw).
#
# **What it does not show.** One seed, one calibration split, five CV folds; the differences between
# representations in MoTrPAC CV are of the order of the fold-to-fold spread printed above.
#
# %% [markdown]
# ### 10.6 Is the rank / pair failure a units problem? (the CPM test)
#
# On GTEx the scale-free representations do *worse* than z-scores at k = 20, the opposite of the
# hope. One obvious suspect is units: GTEx was scored on log2 TPM (length-normalised), MoTrPAC and the
# BodyMap are log2 CPM, and ranks *within* a sample are changed by a per-gene length factor while
# z-scores within a dataset are not. The pipeline tested it by re-scoring GTEx on log2 CPM computed from
# the GTEx read counts (`make transfer`, second command). If units were the cause, rank and pair accuracy
# should rise towards z-score accuracy on CPM.

# %%
_a = s10_p14["target_summary"]
_a = _a[_a["target"] == "GTEx"]
s10_units = _a.merge(s10_p14cpm["target_summary"], on=["target", "selector", "representation", "k"], suffixes=("_tpm", "_cpm"))
s10_units["cpm_minus_tpm"] = s10_units["accuracy_cpm"] - s10_units["accuracy_tpm"]
print(s10_units[["selector", "representation", "k", "accuracy_tpm", "accuracy_cpm", "cpm_minus_tpm"]].round(3).to_string(index=False))
_std = s10_units[s10_units["selector"] == "standard"]
_rp = _std[_std["representation"].isin(["rank", "pairs"])]
s10_units_mean = float(_rp["cpm_minus_tpm"].mean())
_z = _std[_std["representation"] == "zscore"].set_index("k")["accuracy_cpm"]
s10_gap_cpm = float((_rp.set_index("k")["accuracy_cpm"] - _z).mean())
print(f"\nstandard selector, rank + pairs, mean (CPM − TPM) accuracy: {s10_units_mean:+.4f}")
print(f"standard selector, on CPM, mean (rank/pairs − z-score) accuracy: {s10_gap_cpm:+.3f}")
_ta = s10_p14["target_accuracy"]
_tc = s10_p14cpm["target_accuracy"]
_pt20 = (pd.concat([_ta[_ta["target"] == "GTEx"].assign(units="tpm"), _tc.assign(units="cpm")])
         .query("k == 20 and selector == 'standard'")
         .pivot_table(index="organ", columns=["representation", "units"], values="accuracy"))
print("\nPer GTEx tissue, k = 20, standard selector:")
print(_pt20.round(2).to_string())
record("s10.units_rankpairs_cpm_minus_tpm", s10_units_mean, "10")
record("s10.units_rankpairs_minus_z_on_cpm", s10_gap_cpm, "10")

# %% [markdown]
# **What this shows (measured).** Matching the units barely moves rank and pair accuracy on GTEx — the
# mean change is printed above, and it is small against the gap between rank/pairs and z-scores, which
# remains on the CPM matrix. Units are therefore ruled out as the main cause. Within species (BodyMap
# adults) the rank representation does transfer (target-accuracy table in 10.5). The per-tissue table
# shows the failure is concentrated: some tissues (kidney, liver, skeletal muscle, ovary) are called
# well by ranks, others (adrenal, aorta, brain regions, heart, spleen, testis, blood) almost never.
#
# **What it does not show — the failure is unexplained.** The test rules out one cause (TPM vs CPM
# units); it does not identify the actual cause. Candidate explanations — species differences in the
# within-sample ordering of the panel genes, differences in library protocol, the 5,000-gene prefilter
# being chosen on rat variance — have not been tested here and are not claimed. Note also that the
# transfer-aware selector *raises* rank/pair accuracy on GTEx at k = 20 while lowering it at k = 50
# (10.5), which no single one of those candidates obviously predicts.
#
# **What this section as a whole does not show.** GTEx donors are post-mortem adults with their own
# ischaemic-time and cause-of-death effects; the transfer numbers mix those with the species effect.
# Only one GTEx subset (≤ 150 donors per tissue) and one seed for the calibration split were run.
