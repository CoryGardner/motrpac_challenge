# %% [markdown]
# ## 9. External validation within species: the rat BodyMap
#
# **Question.** The tissue fingerprint was fit on MoTrPAC's 50 animals, all sequenced by one
# consortium with one protocol. Does it still name the tissue when it is applied to a different rat
# dataset — the rat BodyMap (GSE53960, via the Bioconductor `bodymapRat` package), a different lab
# and library protocol, with animals at four ages (2, 6, 21 and 104 weeks)? And does the
# conformal *guarantee* (90 % coverage at α = 0.10) travel with the accuracy?
#
# **Why it matters.** A within-cohort cross-validation number says nothing about a new lab. This is
# the cheapest external test available: same species, so no orthology step, and the 21-week adults are
# the closest in age to the MoTrPAC animals.
#
# **Design (from `scripts/12_bodymap_validate.py`, flags of `make bodymap`).** Genes are matched by
# Ensembl id and z-scored *within each dataset*. Muscle and brain are scored as super-classes (BodyMap
# "Muscle" is correct if called SKM-GN or SKM-VL; "Brain" if called CORTEX, HIPPOC or HYPOTH); thymus
# and uterus have no MoTrPAC counterpart and are out-of-distribution. The BodyMap has no animal
# identifier, so the animal is taken to be the replicate index within sex × stage — an assumption the
# pipeline states. Panels of k = 20 and 50 genes and the "full" model (5,000-gene variance prefilter)
# are fit on all 50 MoTrPAC animals.
#
# **What to look for.** (i) Accuracy on the adults, and how it changes with age; (ii) whether a panel
# selected *on the BodyMap itself* does better than the imported one, in the same scope; (iii) the
# conformal coverage with MoTrPAC calibration versus after recalibrating on a few BodyMap animals;
# (iv) *why* the juvenile testis and spleen fail.
#
# **Mode.** The pipeline phase takes a few minutes, so with `RECOMPUTE = False` its result tables are
# loaded from `results/12_bodymap/` (written after the 2026-09-25 conformal-quantile fix). With
# `RECOMPUTE = True` the copied script logic below re-runs and writes to `_outputs/s09/`. Cheap checks
# (gene matching, marker expression) are recomputed live in both modes.

# %%
section("9 BodyMap")

# from scripts/12_bodymap_validate.py (verbatim): the organ mapping
ORGAN_MAP = {"Adrenal": {"ADRNL"}, "Brain": {"CORTEX", "HIPPOC", "HYPOTH"}, "Heart": {"HEART"}, "Kidney": {"KIDNEY"},
             "Liver": {"LIVER"}, "Lung": {"LUNG"}, "Muscle": {"SKM-GN", "SKM-VL"}, "Spleen": {"SPLEEN"},
             "Testes": {"TESTES"}, "Thymus": None, "Uterus": None}
TISSUE_TO_ORGAN = {t: o for o, ts in ORGAN_MAP.items() if ts for t in ts}


# from scripts/12_bodymap_validate.py::load_bodymap (verbatim; renamed s09_load_bodymap because
# section 10 copies phase 14's different load_bodymap)
def s09_load_bodymap():
    X = pd.read_csv(C.EXTERNAL_DIR / "bodymap_counts.csv", index_col=0)
    X.index = X.index.astype(str)
    meta = pd.read_csv(C.EXTERNAL_DIR / "bodymap_meta.csv", dtype=str).set_index("sample")
    meta = meta.loc[X.columns]
    meta["stage_weeks"] = meta["stage_weeks"].astype(int)
    meta["sex"] = meta["sex"].map({"F": "female", "M": "male"}).fillna(meta["sex"])
    meta["replicate_index"] = meta["replicate"].astype(str).str.rsplit("_", n=1).str[-1]
    meta["animal_id"] = meta["sex"].str[0] + "_" + meta["stage_weeks"].astype(str) + "_" + meta["replicate_index"]
    return X.T, meta  # samples × genes


# the flags of `make bodymap` (script defaults; cli.common_parser gives prefilter 5000, seed C.SEED)
s09_args = SimpleNamespace(alpha=0.10, grid="20,50", cal_frac=0.3, recal_animals="3,5", recal_repeats=20,
                           prefilter=5000, seed=C.SEED)
S09_OUT = OUT / "s09"
S09_OUT.mkdir(parents=True, exist_ok=True)
# pre-fix copies of the pipeline's results (for "published" columns); read-only
PREFIX = REPO / "backup" / "pipeline_history" / "results_pre_quantile_fix_2026-09-25"

# -- data: the stacked MoTrPAC TRNSCRPT counts (same call as phases 04/12/13/14; shared cache key) --
pheno = cached("pheno", io.load_pheno)
om_trn = cached("TRNSCRPT_counts", lambda: io.stack_tissues(
    "TRNSCRPT", tissues=None, source="counts", join="inner", pheno=pheno, complete=False, drop_incomplete_samples=None))
s09_Xb_counts, s09_mb = s09_load_bodymap()
s09_lb = io.log_cpm(s09_Xb_counts, log=True)                 # BodyMap total-library log2 CPM (as the script)
s09_sym = io.map_to_gene_symbols(om_trn.X.columns, io.load_feature_to_gene())

# -- gene matching by Ensembl id, recomputed live (script 12, "data" block) --
_bm_genes = set(s09_Xb_counts.columns)
s09_shared = [g for g in om_trn.X.columns if g in _bm_genes]
s09_overlap = {"motrpac_genes": om_trn.n_features, "bodymap_genes": s09_Xb_counts.shape[1],
               "shared_ensembl": len(s09_shared), "motrpac_only": om_trn.n_features - len(s09_shared),
               "shared_frac_of_motrpac": len(s09_shared) / om_trn.n_features}
_pub = pd.read_csv(res("12_bodymap", "gene_overlap.csv")).iloc[0]
print(pd.DataFrame({"recomputed": pd.Series(s09_overlap), "published": _pub[list(s09_overlap)]}).to_string())
print(f"\nBodyMap: {len(s09_mb)} samples, {s09_mb['animal_id'].nunique()} animal ids; samples per organ × stage:")
print(pd.crosstab(s09_mb["organ"], s09_mb["stage_weeks"]).to_string())
record("s09.shared_genes", len(s09_shared), "09")
record("s09.bodymap_samples", len(s09_mb), "09")

# %%
# ---- the pipeline phase: scripts/12_bodymap_validate.py::main, copied (trimmed: argparse, report.add_section and
# the printing dropped; results written to `out`, a folder under OUT; the native-panel CV predictions are
# additionally kept, so the native panel can be scored on the adults alone — see 9.2) ----
def s09_run_phase12(out: Path) -> dict:
    args = s09_args
    grid = [int(k) for k in args.grid.split(",")]
    om = om_trn
    Xb_counts, mb = s09_Xb_counts, s09_mb.copy()
    shared = s09_shared
    sym = s09_sym
    lb = s09_lb
    Zm = zscore(om.X[shared])                                   # MoTrPAC log2 CPM → z within MoTrPAC
    Zb_all = zscore(lb[shared])                                 # z within BodyMap (all 316 samples)
    y = om.meta["tissue"].astype(str).to_numpy()
    g = om.groups()
    mb["group_id"] = mb["animal_id"]
    n_animals_bm = mb["animal_id"].nunique()
    # ---- panels on all 50 MoTrPAC animals ----
    L_src = om.X[shared].to_numpy(dtype=float)
    pm = PanelModels(L_src, y, shared, grid, args.prefilter)
    Lb_pre, Zb_pre = pm.target_matrices(lb[shared].to_numpy(dtype=float))
    surv = []
    core_p = C.RESULTS_DIR / "05_panels" / "TRNSCRPT" / "candidate_panel_annotated.csv"
    core = pd.read_csv(core_p) if core_p.exists() else None
    for k in grid:
        genes = pm.panel_genes(k)
        surv.append({"panel": f"k={k} (all-animal fit)", "n_genes": k, "present_in_bodymap": sum(gg in set(Xb_counts.columns) for gg in genes),
                     "genes": ";".join(str(sym.get(gg, gg)) for gg in genes)})
    if core is not None:
        surv.append({"panel": "T5 stable core", "n_genes": len(core), "present_in_bodymap": int(core["feature_ID"].isin(Xb_counts.columns).sum()),
                     "genes": ";".join(core["gene_symbol"].astype(str))})
    surv = pd.DataFrame(surv)
    model_names = [f"k{k}" for k in grid] + ["full"]
    # ---- primary test: adults, and every stage ----
    acc_rows, conf = [], {}
    for stage in sorted(mb["stage_weeks"].unique()):
        m_s = (mb["stage_weeks"] == stage).to_numpy()
        for name in model_names:
            pred, ok, _, _ = score_block(pm, name, Lb_pre[m_s], Zb_pre[m_s], mb[m_s], ORGAN_MAP)
            sub = mb[m_s].assign(pred=pred, correct=ok)
            for organ, d in sub.groupby("organ"):
                acc_rows.append({"stage_weeks": stage, "model": name, "organ": organ, "n": len(d),
                                 "accuracy": float(d["correct"].mean()) if ORGAN_MAP.get(organ) else np.nan,
                                 "top_prediction": d["pred"].value_counts().index[0],
                                 "top_prediction_frac": float(d["pred"].value_counts().iloc[0] / len(d))})
            if stage == 21:
                conf[name] = pd.crosstab(sub["organ"], sub["pred"]).reindex(columns=pm.classes, fill_value=0)
    acc = pd.DataFrame(acc_rows)
    mapped = acc[acc["accuracy"].notna()]
    overall = mapped.groupby(["stage_weeks", "model"]).apply(lambda d: np.average(d["accuracy"], weights=d["n"]), include_groups=False).unstack("model").reset_index()
    # ---- conformal transfer (shared code path) ----
    rng = np.random.default_rng(args.seed)
    uniq = np.unique(g)
    rng.shuffle(uniq)
    n_cal = max(1, int(round(args.cal_frac * len(uniq))))
    cal_animals = set(uniq[:n_cal])
    is_cal = np.array([a in cal_animals for a in g])
    fit_idx, cal_idx = np.flatnonzero(~is_cal), np.flatnonzero(is_cal)
    assert_no_group_leak(om.meta, fit_idx, cal_idx)
    pm_c = PanelModels(L_src[fit_idx], y[fit_idx], shared, grid, args.prefilter)
    classes = pm_c.classes
    Lm_c, Zm_c = pm_c.source_matrices(L_src[cal_idx])
    calib = calibrate_models(pm_c, Lm_c, Zm_c, y[cal_idx], classes, args.alpha, model_names)
    Lt_c, Zt_c = pm_c.target_matrices(lb[shared].to_numpy(dtype=float))
    adults = (mb["stage_weeks"] == 21).to_numpy()
    conf_df, ood, recal = conformal_transfer(calib, classes, Lt_c, Zt_c, mb, ORGAN_MAP, args.alpha, "stage_weeks", 21,
                                             [int(v) for v in args.recal_animals.split(",") if v.strip()], args.recal_repeats, rng,
                                             model_names, ["Thymus", "Uterus"])
    po = per_organ_coverage(calib, classes, Lt_c[adults], Zt_c[adults], mb[adults], ORGAN_MAP, model_names)
    # ---- age shift table ----
    age = overall.copy()
    cov_age = conf_df[conf_df["conformal"] == "marginal"].pivot(index="stage_weeks", columns="model", values="coverage_mapped").add_prefix("cov_").reset_index()
    age = age.merge(cov_age, on="stage_weeks")
    for name in model_names:
        ref = age.loc[age["stage_weeks"] == 21, name].iloc[0]
        age[f"acc_drop_{name}"] = ref - age[name]
        refc = age.loc[age["stage_weeks"] == 21, f"cov_{name}"].iloc[0]
        age[f"cov_drop_{name}"] = refc - age[f"cov_{name}"]
    # ---- native BodyMap panel under animal-grouped CV ----
    Zb_df = pd.DataFrame(Zb_all, index=mb.index, columns=shared)
    meta_bm = pd.DataFrame({"pid": mb["animal_id"].to_numpy(), "bid": mb.index, "sex": mb["sex"].to_numpy(),
                            "group": mb["stage_weeks"].astype(str).to_numpy(), "tissue": mb["organ"].to_numpy(), "assay": "BODYMAP"}, index=mb.index)
    om_bm = io.OmicsMatrix(Zb_df, meta_bm, pd.DataFrame(index=shared), "BODYMAP", name="bodymap")
    curve, sel = models.panel_curve(om_bm, label="tissue", grid=[20], kind="logreg_l2", n_splits=5, prefilter=args.prefilter,
                                    quick=True, seed=args.seed, verbose=False)
    native_acc = float(curve["balanced_accuracy"].mean())
    imported = {}
    for name in model_names:
        _, ok, _, _ = score_block(pm, name, Lb_pre, Zb_pre, mb, ORGAN_MAP)
        m_map = np.array([bool(ORGAN_MAP.get(o)) for o in mb["organ"]])
        imported[name] = float(ok[m_map].mean())
    _, _, preds = models.panel_curve(om_bm, label="tissue", grid=[20], kind="logreg_l2", n_splits=5, prefilter=args.prefilter,
                                     quick=True, seed=args.seed, verbose=False, return_predictions=True)
    pm_ok = preds["y_true"].map(lambda o: bool(ORGAN_MAP.get(o)))
    native_map = float((preds.loc[pm_ok, "y_true"] == preds.loc[pm_ok, "y_pred"]).mean())
    native = pd.DataFrame([{"native_panel_k20_balanced_accuracy_11_organs": native_acc,
                            "native_panel_k20_accuracy_9_mapped_organs": native_map,
                            **{f"imported_{n}_accuracy_9_mapped_organs_all_ages": v for n, v in imported.items()},
                            "gap_native_minus_imported_k20": native_map - imported["k20"], "n_bodymap_animals": n_animals_bm,
                            "native_panel_genes": ";".join(str(sym.get(gg, gg)) for gg in sel[sel["fold"] == 0]["feature_ID"].tolist())}])
    # ---- per-gene failure analysis ----
    genes20 = pm.panel_genes(max(grid[0], 20)) if 20 in grid else pm.panel_genes(grid[0])
    flags = {}
    for fn in ("stability_k20_annotated.csv", "candidate_panel_annotated.csv"):
        p = C.RESULTS_DIR / "05_panels" / "TRNSCRPT" / fn
        if p.exists():
            a = pd.read_csv(p)
            for _, r in a.iterrows():
                flags.setdefault(r["feature_ID"], {"risk_T7_regulated": bool(r.get("risk_T7_regulated", False)),
                                                  "risk_qc_correlated": bool(r.get("risk_qc_correlated", False))})
    sub_ad = mb[adults]
    Zb_ad = pd.DataFrame(Zb_all[adults], index=sub_ad.index, columns=shared)
    Zm_df = pd.DataFrame(Zm, index=om.X.index, columns=shared)
    rows = []
    for gg in genes20:
        j = list(pm.genes_pre).index(gg)
        i_t = int(np.argmax(pm.scores[:, j]))
        marker = pm.classes[i_t]
        organ = TISSUE_TO_ORGAN.get(marker)
        mu_m = Zm_df[gg].groupby(y).mean()
        eff_m = float(mu_m[marker] - mu_m.drop(marker).max())
        row = {"feature_ID": gg, "gene_symbol": sym.get(gg, gg), "marker_tissue": marker, "bodymap_organ": organ or "not in BodyMap",
               "motrpac_effect_z": eff_m}
        if organ:
            mu_b = Zb_ad[gg].groupby(sub_ad["organ"]).mean()
            eff_b = float(mu_b[organ] - mu_b.drop(organ).max())
            row.update({"bodymap_effect_z": eff_b, "bodymap_top_organ": str(mu_b.idxmax()),
                        "fails_in_bodymap": bool(eff_b <= 0), "weakened": bool(0 < eff_b < 0.5 * eff_m)})
        else:
            row.update({"bodymap_effect_z": np.nan, "bodymap_top_organ": None, "fails_in_bodymap": None, "weakened": None})
        row.update(flags.get(gg, {"risk_T7_regulated": None, "risk_qc_correlated": None}))
        rows.append(row)
    gene_check_df = pd.DataFrame(rows)
    jv = []
    for gg in genes20:
        j = list(pm.genes_pre).index(gg)
        marker = pm.classes[int(np.argmax(pm.scores[:, j]))]
        m = pd.DataFrame({"v": lb[gg].to_numpy(dtype=float), "organ": mb["organ"].to_numpy(), "stage_weeks": mb["stage_weeks"].to_numpy()})
        for (organ, stage), dd in m.groupby(["organ", "stage_weeks"]):
            jv.append({"feature_ID": gg, "gene_symbol": sym.get(gg, gg), "marker_tissue": marker, "organ": organ, "stage_weeks": int(stage),
                       "mean_log2_cpm": float(dd["v"].mean()), "n": int(len(dd))})
    tables = {"panel_survival": surv, "accuracy_by_organ": acc, "age_shift_accuracy": overall, "conformal_transfer": conf_df,
              "ood_sets": ood, "recalibration": recal, "coverage_by_organ": po, "age_shift": age, "native_panel": native,
              "panel_gene_check": gene_check_df, "juvenile_marker_check": pd.DataFrame(jv),
              **{f"confusion_{n}_adult": conf[n].reset_index() for n in model_names}}
    for name, t in tables.items():
        t.to_csv(out / f"{name}.csv", index=False)
    preds.to_csv(out / "native_panel_cv_predictions.csv", index=False)      # [notebook] extra output
    return tables


S09_TABLES = ["panel_survival", "accuracy_by_organ", "age_shift_accuracy", "conformal_transfer", "ood_sets",
              "recalibration", "coverage_by_organ", "age_shift", "native_panel", "panel_gene_check",
              "juvenile_marker_check", "confusion_k20_adult", "confusion_k50_adult", "confusion_full_adult"]
if RECOMPUTE:
    t0 = time.perf_counter()
    s09_t = s09_run_phase12(S09_OUT)
    print(f"phase 12 recomputed in {time.perf_counter() - t0:.0f} s → {S09_OUT}")
else:
    s09_t = {n: pd.read_csv(res("12_bodymap", f"{n}.csv")) for n in S09_TABLES}
    print("phase 12 tables loaded from", RES / "12_bodymap")
# the same tables as the pipeline published them before the 2026-09-25 quantile fix
s09_pre = {n: pd.read_csv(PREFIX / "12_bodymap" / f"{n}.csv") for n in S09_TABLES}

# %% [markdown]
# ### 9.1 Does the imported fingerprint name the BodyMap organs?
#
# The panel-survival table says whether the 20 and 50 MoTrPAC panel genes exist in the BodyMap
# annotation at all. Then: per-organ accuracy on the 21-week adults, the sample-weighted accuracy at
# every age, and the adult confusion matrix of the 20-gene panel (rows = BodyMap organ, columns = the
# predicted MoTrPAC tissue; zero columns dropped). These numbers do not involve conformal prediction,
# so the quantile fix cannot have changed them — the cell checks that too.

# %%
print(s09_t["panel_survival"].drop(columns=["genes"]).to_string(index=False))
_acc = s09_t["accuracy_by_organ"]
s09_acc_adult = _acc[(_acc["stage_weeks"] == 21) & _acc["accuracy"].notna()].pivot(index="organ", columns="model", values="accuracy")
print("\nPer-organ accuracy, 21-week adults (super-class scoring):")
print(s09_acc_adult[["k20", "k50", "full"]].round(2).to_string())
s09_age = s09_t["age_shift_accuracy"].set_index("stage_weeks")[["k20", "k50", "full"]]
print("\nSample-weighted accuracy on the 9 mapped organs, by age (weeks):")
print(s09_age.round(3).to_string())
_pre_age = s09_pre["age_shift_accuracy"].set_index("stage_weeks")[["k20", "k50", "full"]]
print("max |this − pre-fix| over the age table:", float((s09_age - _pre_age).abs().max().max()))
_cf = s09_t["confusion_k20_adult"].set_index("organ")
print("\nAdult confusion, k = 20:")
print(_cf.loc[:, (_cf > 0).any(axis=0)].to_string())
for k in ("k20", "k50", "full"):
    record(f"s09.acc_adult_{k}", s09_age.loc[21, k], "09")
record("s09.acc_2wk_k20", s09_age.loc[2, "k20"], "09")
record("s09.acc_104wk_k20", s09_age.loc[104, "k20"], "09")

# %% [markdown]
# **What this shows.** All panel genes are present in the BodyMap, and on the adults the imported
# panels name almost every mapped organ; a muscle sample called SKM-VL instead of SKM-GN is
# still correct under super-class scoring. Accuracy falls at the two ends of life, most at 2 weeks and
# for the 20-gene panel (9.4 asks why). The out-of-distribution rows (thymus, uterus) are forced into
# some MoTrPAC tissue by the argmax — the classifier has no "none of these" answer; only the conformal
# sets in 9.3 can abstain.
#
# **What it does not show.** The adult test is small (a handful of samples per organ), so a perfect
# per-organ accuracy has wide uncertainty. Several panel markers belong to MoTrPAC tissues the BodyMap does
# not have (ovary, adipose, blood, BAT, small intestine, colon, vena cava) and are untested here.
#
# ### 9.2 Native versus imported panel — in the same scope
#
# A 20-gene panel selected and cross-validated *on the BodyMap itself* (5 animal-grouped folds, all
# samples, all ages, all 11 organs) is the ceiling an imported panel can be compared with. The pipeline
# reports the native panel on two scopes (balanced accuracy over all 11 organs; plain accuracy on the 9
# mapped organs) and the imported panels on the 9 mapped organs at all ages. Scores computed on
# different scopes are not comparable, so the table below puts each number in its scope; the adults-only
# scope for the imported panels is the 21-week row of the age table above.

# %%
_nat = s09_t["native_panel"].iloc[0]
s09_scope = pd.DataFrame({
    "native k20 (BodyMap CV)": [_nat["native_panel_k20_balanced_accuracy_11_organs"], _nat["native_panel_k20_accuracy_9_mapped_organs"], np.nan],
    "imported k20": [np.nan, _nat["imported_k20_accuracy_9_mapped_organs_all_ages"], s09_age.loc[21, "k20"]],
    "imported k50": [np.nan, _nat["imported_k50_accuracy_9_mapped_organs_all_ages"], s09_age.loc[21, "k50"]],
    "imported full": [np.nan, _nat["imported_full_accuracy_9_mapped_organs_all_ages"], s09_age.loc[21, "full"]]},
    index=["11 organs, all ages (balanced acc.)", "9 mapped organs, all ages (accuracy)", "9 mapped organs, adults only (accuracy)"])
# native panel on the adults alone: the pipeline saves only the all-ages summary, so the native CV predictions are
# recomputed here in both modes (cheap: 316 samples) with the same call as scripts/12_bodymap_validate.py::main
# (copied; the `return_predictions=True` call), then checked against the published all-ages number
_zb = pd.DataFrame(zscore(s09_lb[s09_shared]), index=s09_mb.index, columns=s09_shared)
_meta_bm = pd.DataFrame({"pid": s09_mb["animal_id"].to_numpy(), "bid": s09_mb.index, "sex": s09_mb["sex"].to_numpy(),
                         "group": s09_mb["stage_weeks"].astype(str).to_numpy(), "tissue": s09_mb["organ"].to_numpy(),
                         "assay": "BODYMAP"}, index=s09_mb.index)
_om_bm = io.OmicsMatrix(_zb, _meta_bm, pd.DataFrame(index=s09_shared), "BODYMAP", name="bodymap")
_, _, s09_native_pred = models.panel_curve(_om_bm, label="tissue", grid=[20], kind="logreg_l2", n_splits=5,
                                           prefilter=s09_args.prefilter, quick=True, seed=s09_args.seed, verbose=False,
                                           return_predictions=True)
_p = s09_native_pred.merge(s09_mb[["stage_weeks"]], left_on="viallabel", right_index=True)
_p = _p[_p["y_true"].map(lambda o: bool(ORGAN_MAP.get(o)))]
_native_all_live = float((_p["y_true"] == _p["y_pred"]).mean())
print(f"native k20, 9 mapped organs, all ages: recomputed {_native_all_live:.4f}  published "
      f"{_nat['native_panel_k20_accuracy_9_mapped_organs']:.4f}")
_pa = _p[_p["stage_weeks"] == 21]
s09_scope.iloc[2, 0] = float((_pa["y_true"] == _pa["y_pred"]).mean())
print(f"native k20 on the adults alone: {len(_pa)} mapped adult samples\n")
record("s09.native_acc_9_live", _native_all_live, "09")
record("s09.native_acc_9_adults", s09_scope.iloc[2, 0], "09")
record("s09.imported_k20_acc_9_adults", s09_scope.iloc[2, 1], "09")
print(s09_scope.round(3).to_string())
print("\ngap native − imported k20, 9 mapped organs all ages:", round(float(_nat["gap_native_minus_imported_k20"]), 3))
print("native panel genes (fold 0):", _nat["native_panel_genes"].replace(";", ", "))
record("s09.native_bal_acc_11", _nat["native_panel_k20_balanced_accuracy_11_organs"], "09")
record("s09.native_acc_9", _nat["native_panel_k20_accuracy_9_mapped_organs"], "09")
record("s09.imported_k20_acc_9_all_ages", _nat["imported_k20_accuracy_9_mapped_organs_all_ages"], "09")

# %% [markdown]
# **What this shows.** In the one scope where both are measured (9 mapped organs, all ages) the native
# panel is ahead of the imported 20-gene panel, and the gap narrows as the imported model grows to 50
# genes and to the full 5,000-gene model. On the adults alone the imported panels are already at or
# near the top of the scale, as is the native panel, so the all-ages gap comes from the juvenile and aged
# animals: see the age table above, and 9.4.
#
# **What it does not show.** The native number is a cross-validation estimate on assumed animal
# ids (see the design note), not an independent test; the imported number is a true external test. The native panel is
# almost entirely different genes from the MoTrPAC panel (compare the gene lists), so "native beats
# imported" does not say the MoTrPAC markers are wrong — only that other markers separate the BodyMap
# organs as well. The two scopes answer different questions: on the adults alone both panels are at
# the top of the scale (the last row), while across all four ages the native panel leads (the middle
# row), so the all-ages gap is an age effect, not a panel-quality effect. The native adults-only cell
# is recomputed in this notebook: the pipeline saves only the all-ages native summary.
#
# ### 9.3 Does the conformal guarantee travel?
#
# Models are refit on 70 % of the MoTrPAC animals and calibrated (LAC, α = 0.10) on the held-out
# 30 %; the sets are then formed on the BodyMap. Three variants: marginal, Mondrian (one quantile per
# tissue) and floored Mondrian (each tissue's quantile at least the marginal one). Then the pooled
# threshold is *recalibrated* on 3 or 5 BodyMap adult animals and tested on the remaining adults (20
# draws). These numbers depend on the conformal quantile, so each one is printed next to its pre-fix
# value.

# %%
_key = ["stage_weeks", "model", "conformal"]
_ct = s09_t["conformal_transfer"].merge(s09_pre["conformal_transfer"], on=_key, suffixes=("", "_prefix"))
s09_cov = _ct[_ct["stage_weeks"] == 21][_key[1:] + ["coverage_mapped", "coverage_mapped_prefix", "frac_empty_mapped",
                                                     "frac_empty_mapped_prefix", "avg_set_size_mapped", "ood_frac_empty"]]
print("Adults (21 weeks), MoTrPAC calibration, α = 0.10:")
print(s09_cov.round(3).to_string(index=False))
print("\nMarginal coverage (k20) by age, this vs pre-fix:")
print(_ct[(_ct["model"] == "k20") & (_ct["conformal"] == "marginal")][["stage_weeks", "coverage_mapped", "coverage_mapped_prefix"]].round(3).to_string(index=False))
_rc = s09_t["recalibration"].merge(s09_pre["recalibration"], on=["model", "n_recal"], suffixes=("", "_prefix"))
s09_recal = _rc[["model", "n_recal", "n_test_individuals", "coverage_recalibrated", "coverage_recalibrated_prefix",
                 "coverage_source_cal_same_test", "coverage_source_cal_same_test_prefix", "set_size_recalibrated", "set_size_recalibrated_prefix"]]
print("\nRecalibration on BodyMap adult animals (mean of 20 draws):")
print(s09_recal.round(3).to_string(index=False))
print("\nWhat the MoTrPAC-calibrated sets hold for the out-of-distribution organs (adults):")
print(s09_t["ood_sets"][s09_t["ood_sets"]["conformal"] == "marginal"].to_string(index=False))
_c = s09_cov.set_index(["model", "conformal"])
for m in ("k20", "k50", "full"):
    for v in ("marginal", "mondrian", "floored"):
        record(f"s09.cov_adult_{m}_{v}", _c.loc[(m, v), "coverage_mapped"], "09")
_r = s09_recal.set_index(["model", "n_recal"])
for m in ("k20", "full"):
    for n in (3, 5):
        record(f"s09.recal_cov_{m}_n{n}", _r.loc[(m, n), "coverage_recalibrated"], "09")
    record(f"s09.srccal_cov_same_test_{m}", _r.loc[(m, 3), "coverage_source_cal_same_test"], "09")

# %% [markdown]
# After the fix, the Mondrian and floored sets are never empty and hold two or more tissues on average
# (table above). The next cell checks, live, the explanation the textbook rule suggests: the script's
# calibration split (the first draw of its seeded generator) leaves some tissues with so few calibration
# vials that ⌈(n+1)(1 − α)⌉ > n, so their per-tissue threshold is infinite and they enter *every* set.

# %%
# the calibration animals of scripts/12_bodymap_validate.py (same generator, same first draw: the shuffle of the
# unique animal ids comes before any other use of `rng` in the script)
_g = om_trn.groups()
_uniq = np.unique(_g)
np.random.default_rng(s09_args.seed).shuffle(_uniq)
_cal = set(_uniq[:max(1, int(round(s09_args.cal_frac * len(_uniq))))])
_n_cal = om_trn.meta.loc[[a in _cal for a in _g], "tissue"].value_counts().sort_index()
s09_mondrian = pd.DataFrame({"n_cal_vials": _n_cal,
                             "rank_needed": np.ceil((_n_cal + 1) * (1 - s09_args.alpha) - 1e-9).astype(int)})
s09_mondrian["threshold_infinite"] = s09_mondrian["rank_needed"] > s09_mondrian["n_cal_vials"]
print(f"{len(_cal)} calibration animals; per-tissue calibration vials and the textbook Mondrian rank (α = {s09_args.alpha}):")
print(s09_mondrian.T.to_string())
print("tissues in every Mondrian / floored set:", list(s09_mondrian.index[s09_mondrian["threshold_infinite"]]))
_po = s09_t["coverage_by_organ"]
print("\nk = 20 Mondrian, adults: mean set size per BodyMap organ")
print(_po[(_po["model"] == "k20") & (_po["conformal"] == "mondrian")].set_index("organ")[["n", "coverage", "frac_empty", "avg_set_size"]].round(2).to_string())
record("s09.n_mondrian_infinite", int(s09_mondrian["threshold_infinite"].sum()), "09")

# %% [markdown]
# **What this shows.** Accuracy transfers but the guarantee does not: with MoTrPAC calibration the
# marginal coverage on BodyMap adults is far below the 90 % target, and almost all of the marginal
# shortfall is *empty* sets (compare coverage with the empty fraction) — on the new library chemistry the model's top
# probability is lower than on MoTrPAC and falls under the MoTrPAC threshold. Plain Mondrian covers
# less than marginal; the floored variant covers the most. After the quantile fix the Mondrian and
# floored sets are no longer empty, but not because they got better: the tissues with too few
# calibration vials (the sex-specific ones, per the table above) now have an infinite threshold and sit
# in every set, which raises the set size without covering the right organ. Recalibrating on a few
# BodyMap animals restores coverage near the target with sets of about one tissue. Under marginal
# calibration thymus and uterus, which MoTrPAC never saw, mostly receive empty sets — the desired "none
# of these" behaviour; under Mondrian they receive the always-included tissues instead. The pre-fix
# columns show what the quantile fix moved (not the marginal adult coverage, nor the source-calibrated
# coverage on the recalibration test animals).
#
# **What it does not show.** The recalibration draws share test animals (there are few adult animal
# ids — see the sample table in the first cell), so the 20 draws are not independent and their mean has
# no simple error bar. Coverage "restored" on the pooled organs says nothing about per-organ coverage
# after recalibration. The calibration labels used for recalibration are the organ's best-scoring
# mapped tissue (a muscle sample calibrates as whichever of SKM-GN / SKM-VL the model prefers), which
# makes the super-class task slightly easier.
#
# ### 9.4 Why the juvenile testis and spleen fail
#
# At 2 weeks the 20-gene panel gets two organs wrong in a systematic way. The first table reads, from
# the phase's own output, what the panel calls each 2-week organ. The second reads the relevant panel
# markers' mean expression (log2 CPM, BodyMap) in the juvenile and the adult organ: Pgk2 (the panel's
# only testis marker) and Mybph (an SKM-VL marker) in testes; Hbq1b (the blood marker) and Fcrl5 (the
# spleen marker) in spleen.

# %%
_a2 = _acc[(_acc["model"] == "k20") & (_acc["organ"].isin(["Testes", "Spleen"]))]
print(_a2.pivot_table(index="organ", columns="stage_weeks", values="top_prediction", aggfunc="first").to_string())
print(_a2.pivot_table(index="organ", columns="stage_weeks", values="accuracy").round(2).to_string())
_jv = s09_t["juvenile_marker_check"]
_pairs = {"Testes": ["Pgk2", "Mybph"], "Spleen": ["Hbq1b", "Fcrl5"]}
s09_juv = (_jv[[gs in _pairs.get(o, []) for gs, o in zip(_jv["gene_symbol"], _jv["organ"])]]
           .pivot_table(index=["organ", "gene_symbol", "marker_tissue"], columns="stage_weeks", values="mean_log2_cpm"))
print("\nMean log2 CPM in BodyMap, by age (weeks):")
print(s09_juv.round(2).to_string())
record("s09.acc_2wk_testes_k20", float(_a2[(_a2["organ"] == "Testes") & (_a2["stage_weeks"] == 2)]["accuracy"].iloc[0]), "09")
record("s09.acc_2wk_spleen_k20", float(_a2[(_a2["organ"] == "Spleen") & (_a2["stage_weeks"] == 2)]["accuracy"].iloc[0]), "09")
record("s09.pgk2_testes_2wk", s09_juv.loc[("Testes", "Pgk2", "TESTES"), 2], "09")
record("s09.pgk2_testes_21wk", s09_juv.loc[("Testes", "Pgk2", "TESTES"), 21], "09")

# %% [markdown]
# **Supporting evidence beyond the panel (a notebook addition, not a pipeline output).** The two
# explanations proposed in the pipeline's summary each predict something about genes *outside* the
# panel, which can be checked directly in the BodyMap:
# - if the juvenile testis lacks Pgk2 because it has no post-meiotic germ cells yet, other
#   spermatid-specific genes (protamine Prm1, transition proteins Tnp1/Tnp2, acrosomal Acrv1) should
#   also be near zero at 2 weeks and high in adults;
# - if the juvenile spleen looks like blood because it is still haematopoietic, erythroid-lineage
#   genes (Alas2, Klf1, Gypa, the adult β-globin Hbb) should be higher in the 2-week spleen than in the
#   adult spleen.
#
# The gene choice is background knowledge (the textbook markers of each cell type), stated before
# looking; the values are computed live from the BodyMap counts.

# %%
_sym2id = s09_sym.dropna().reset_index().drop_duplicates("gene_symbol").set_index("gene_symbol")["index"]
s09_support_genes = {"Testes": ["Pgk2", "Prm1", "Tnp1", "Tnp2", "Acrv1"],
                     "Spleen": ["Hbq1b", "Alas2", "Klf1", "Gypa", "Hbb", "Fcrl5"]}
_rows = []
for organ, genes in s09_support_genes.items():
    for gs in genes:
        gid = _sym2id.get(gs)
        if gid is None or gid not in s09_lb.columns:
            print(f"  {gs}: not found in the BodyMap / MoTrPAC annotation — skipped")
            continue
        v = s09_lb[gid].groupby([s09_mb["organ"], s09_mb["stage_weeks"]]).mean()
        _rows.append({"organ": organ, "gene": gs, **{f"{s} wk": v.get((organ, s), np.nan) for s in (2, 6, 21, 104)},
                      "adult_minus_2wk": v.get((organ, 21), np.nan) - v.get((organ, 2), np.nan)})
s09_support = pd.DataFrame(_rows)
print("Mean log2 CPM in the BodyMap organ, by age:")
print(s09_support.round(2).to_string(index=False))

# %% [markdown]
# **What this shows (measured).** At 2 weeks the panel calls the testis a muscle tissue and the spleen
# blood (first table). The testis marker Pgk2 is essentially absent from the 2-week testis and high in
# the adult one, while the muscle marker Mybph is higher in the juvenile testis than in the adult; the
# blood marker Hbq1b is higher and the spleen marker Fcrl5 lower in the juvenile spleen than in the
# adult. The last table shows whether the non-panel genes move the way each explanation predicts
# (read the `adult_minus_2wk` column: positive for the spermatid genes in testes, negative for the
# erythroid genes in spleen, if the explanations hold).
# They do, for every such gene found in the annotation. The aged (104-week) testis also
# expresses the spermatid genes, Pgk2 included, well below the adult level, and the panel's accuracy on
# aged testes is lower too (first table) — the same single-marker mechanism may be at work there,
# though it was not the pipeline's question.
#
# **Interpretation (not measured here).** Pgk2 is a spermatid-specific glycolytic isozyme; in the rat,
# spermatids first appear weeks after the age of the youngest BodyMap group, so a 2-week testis has no
# cells that make it. The rodent spleen remains a site of erythropoiesis in early postnatal life, so a juvenile spleen
# carries red-cell-lineage transcripts. Both are developmental stage effects acting on single
# markers: a 20-gene panel with one marker per tissue fails when that marker's cell type is absent,
# and the larger models (k = 50, full) degrade less. The panel is an *adult* fingerprint.
#
# **What it does not show.** Hbq1b (ENSRNOG00000028114) is the rn6 globin gene that this repository's
# data notes flag as an annotation artifact (it vanishes in the rn7 reprocessing), so the blood-marker
# half of the spleen story rests on a gene of doubtful identity; the erythroid genes in the last table
# are the check that does not depend on it. The BodyMap has only a few samples per organ × age × sex
# cell (first cell), so none of these differences carries a test; they are descriptive.
