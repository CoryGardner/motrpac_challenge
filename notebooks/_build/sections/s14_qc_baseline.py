# %% [markdown]
# ## 13. NEW ANALYSIS — can RNA-seq QC covariates alone classify tissue?
#
# > **This section is new analysis, not part of the pipeline.** Nothing here was in `results/`; the
# > values it records become the reference values for later runs. It follows the scoping probe
# > `code/probes/track_3/p4_qc_only_baseline.py` (feature lists copied from it), but evaluates on
# > **section 4's exact animal-grouped folds**, so it is directly comparable with the fingerprint
# > baselines there.
#
# **Question.** If a classifier never sees a gene — only the per-library QC numbers the consortium
# reports (RIN, read counts, duplication, mapping and composition fractions) — how well does it
# identify the tissue?
#
# **Why it matters.** The evaluation rules require a tuned simple baseline. The simplest one here is a
# covariates-only model. Each tissue was extracted, library-prepped and sequenced as one batch, so the
# QC numbers can carry the tissue label without any biology. If they do, the within-study fingerprint
# accuracy of section 4 cannot on its own prove that the fingerprint is biology.
#
# **What to look for.** (1) Accuracy per fold, mean ± sd, animals per fold, for the technical,
# composition and combined feature sets, against chance. (2) How many RNA plates, library batches and
# flowcells each tissue spans, and how many tissues share each. (3) The two-part reading at the end.

# %%
section("13 QC-only baseline (new)")
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import accuracy_score, balanced_accuracy_score
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler

# ---- the matrix and folds of section 4: reuse them when section 4 ran, otherwise rebuild identically ----
# from scripts/04_fingerprint_baselines.py::main (trimmed: Makefile `baselines` flags for TRNSCRPT — --assay TRNSCRPT,
# --source auto (→ counts), no --tissues, no --complete-features, no --drop-incomplete-samples); shared cache key
if "S04_OM" in globals():
    s14_om = S04_OM
else:
    _src = "counts" if (C.COUNTS_DIR.exists() and any(C.COUNTS_DIR.glob("TRNSCRPT__*.csv"))) else "norm"  # cli.resolve_source
    s14_om = cached("TRNSCRPT_counts", lambda: io.stack_tissues(
        "TRNSCRPT", tissues=None, source=_src, join="inner", pheno=cached("pheno", io.load_pheno),
        complete=False, drop_incomplete_samples=None))
# the identical call section 4 makes (models.evaluate_cv uses the same one inside)
s14_folds = list(grouped_kfold(s14_om.meta, "tissue", 5, C.SEED))
if "S04_FOLDS" in globals():
    assert len(S04_FOLDS) == len(s14_folds) and all(
        np.array_equal(a_tr, b_tr) and np.array_equal(a_te, b_te) for (a_tr, a_te), (b_tr, b_te) in zip(S04_FOLDS, s14_folds))
    s14_folds = S04_FOLDS
    print("using section 4's folds (S04_FOLDS); rebuilt folds are identical")
else:
    print("section 4 is in notebook 01: folds rebuilt here with the identical call grouped_kfold(om.meta, 'tissue', 5, C.SEED)")
    _fp = OUT / "s04" / "folds.csv"
    if _fp.exists():   # the fold assignment notebook 01 wrote: assert the rebuilt folds are the same
        _f4 = pd.read_csv(_fp, dtype={"viallabel": str}).set_index("viallabel")["fold"]
        _f14 = pd.Series(np.concatenate([[k] * len(te) for k, (_, te) in enumerate(s14_folds)]),
                         index=s14_om.meta.index[np.concatenate([te for _, te in s14_folds])])
        assert _f14.sort_index().equals(_f4.reindex(_f14.index).sort_index()), "rebuilt folds differ from section 4's"
        print("rebuilt folds are identical to the fold assignment notebook 01 (section 4) wrote to _outputs/s04/folds.csv")
    else:
        print("(notebook 01 has not been run here, so the fold identity is by construction only)")
print(f"matrix: {s14_om.n_samples} vials, {s14_om.n_animals} animals, {s14_om.meta['tissue'].nunique()} tissues")

# ---- the consortium RNA-seq QC table (c1.0; read-only), study vials only, aligned to the section-4 rows ----
S14_QC = REPO / "data/quant-id/rat-training-06/c1.0/transcriptomics/qa-qc/motrpac_pass1b-06_transcript-rna-seq_qa-qc-metrics.csv"
_q = pd.read_csv(S14_QC, dtype=str)
print(f"QC table: {_q.shape[0]} rows × {_q.shape[1]} columns")
_q = _q[_q["vial_label"].str.startswith("9")]                     # study vials (reference standards start with 8)
print(f"study vials: {len(_q)}; duplicated vial labels: {int(_q['vial_label'].duplicated().sum())}")
_q = _q.drop_duplicates("vial_label").set_index("vial_label")
_idx = s14_om.meta.index.astype(str)
_miss = int((~_idx.isin(_q.index)).sum())
print(f"section-4 vials without a QC row: {_miss}; QC vials not in the section-4 matrix: {int((~_q.index.isin(_idx)).sum())}")
assert _miss == 0, "every section-4 vial needs a QC row, or the folds would have to be re-indexed"
s14_qc = _q.loc[_idx]                                               # same row order as S04_OM.meta → the folds apply as is
# the tissue label comes from the phenotype join of section 4 (om.meta), not from the QC table's free-text 'Tissue'
s14_y = s14_om.meta["tissue"].astype(str).to_numpy()
print("QC 'Tissue' text ↔ tissue code is one-to-one:",
      bool(pd.crosstab(s14_qc["Tissue"].to_numpy(), s14_y).gt(0).sum(axis=1).eq(1).all()))

# %%
# ---- feature sets, copied from code/probes/track_3/p4_qc_only_baseline.py (verbatim lists) ----
s14_qc = s14_qc.copy()
s14_qc["reads_log10"] = np.log10(s14_qc["reads"].astype(float))
tech = ["RIN", "r_260_280", "r_260_230", "pct_adapter_detected", "pct_trimmed", "pct_GC", "pct_dup_sequence",
        "pct_umi_dup", "pct_multimapped", "median_5_3_bias", "reads_log10", "avg_input_read_length"]
comp = ["pct_rRNA", "pct_globin", "pct_chrM", "pct_chrX", "pct_chrY", "pct_mrna", "pct_coding", "pct_utr",
        "pct_intronic", "pct_intergenic"]
S14_SETS = {"technical": tech, "composition": comp, "all": tech + comp}
del tech, comp

# Model: multinomial logistic regression (C = 1, as in the probe) in a Pipeline with median imputation and
# standardization — all fit inside the training fold (the probe imputed with full-data medians; here it is in-fold).
def s14_model():
    return Pipeline([("imp", SimpleImputer(strategy="median")), ("sc", StandardScaler()),
                     ("clf", LogisticRegression(C=1.0, max_iter=10000))])


_rows = []
for name, cols in S14_SETS.items():
    X = s14_qc[cols].apply(pd.to_numeric, errors="coerce").to_numpy(float)
    for f, (tr, te) in enumerate(s14_folds):
        assert_no_group_leak(s14_om.meta, tr, te)
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            m = s14_model().fit(X[tr], s14_y[tr])
        p = m.predict(X[te])
        _rows.append({"features": name, "n_features": len(cols), "fold": f, "n_test_vials": len(te),
                      "n_test_animals": int(s14_om.meta.iloc[te]["pid"].nunique()),
                      "accuracy": accuracy_score(s14_y[te], p), "balanced_accuracy": balanced_accuracy_score(s14_y[te], p)})
s14_pf = pd.DataFrame(_rows)
print("per fold:")
print(s14_pf.pivot(index="fold", columns="features", values="accuracy")[list(S14_SETS)]
      .join(s14_pf.groupby("fold")[["n_test_vials", "n_test_animals"]].first()).to_string(float_format="%.3f"))
s14_summ = s14_pf.groupby("features", sort=False).agg(n_features=("n_features", "first"), acc_mean=("accuracy", "mean"),
                                                        acc_sd=("accuracy", "std"), bal_acc_mean=("balanced_accuracy", "mean"),
                                                        bal_acc_sd=("balanced_accuracy", "std"))
print("\nsummary (5 animal-grouped folds, mean ± sd):")
print(s14_summ.to_string(float_format="%.3f"))
print(f"chance (1 / number of tissues): {1 / len(np.unique(s14_y)):.3f}")
for name in S14_SETS:
    record(f"s14.acc_{name}", s14_summ.loc[name, "acc_mean"], "14", note="new analysis; section-4 folds")
    record(f"s14.bal_acc_{name}", s14_summ.loc[name, "bal_acc_mean"], "14", note="new analysis; section-4 folds")
(OUT / "s14").mkdir(exist_ok=True)
s14_pf.to_csv(OUT / "s14" / "qc_only_per_fold.csv", index=False)

# %%
# ---- for comparison: the gene-based fingerprint on the same folds (section 4's published baseline summary) ----
_fp = pd.read_csv(res("04_baselines", "TRNSCRPT", "summary.csv"))
print(_fp[[c for c in _fp.columns if c in ("model", "accuracy_mean", "accuracy_std", "balanced_accuracy_mean", "balanced_accuracy_std", "n_test_animals_mean")]]
      .to_string(index=False, float_format="%.3f"))
print("(the probe p4 reported balanced accuracy on its own sex-balanced animal folds; it wrote no file to the repo, so its values"
      " are not reproduced here)")

# %% [markdown]
# **What (1) shows.** With no gene at all, a dozen purely technical QC numbers identify the tissue far
# above chance, the composition fractions do better, and the combined set comes close to the gene-based
# fingerprint on the very same folds. The consortium's library QC table is, in effect, a tissue label.
#
# **What it does not show.** *Why* the QC numbers carry tissue. The composition fractions (mitochondrial,
# globin, rRNA, intronic reads) differ between tissues for biological reasons too, so the composition set
# is not "pure batch". The technical set is closer to pure processing, but RIN and duplication also
# depend on the tissue's RNA. The next cell looks at the processing batches directly.

# %%
# ---- (2) processing batches per tissue: which QC-table columns encode plate, library batch and flowcell ----
S14_BATCH_COLS = {"RNA_extr_plate_ID": "RNA extraction plate", "Lib_batch_ID": "library batch",
                  "Lib_prep_date": "library prep date", "Seq_flowcell_ID": "flowcell", "Seq_flowcell_lane": "flowcell lane",
                  "Seq_batch": "sequencing batch", "GET_site": "sequencing site"}
_b = s14_qc[list(S14_BATCH_COLS)].assign(tissue=s14_y)
_per_tissue = _b.groupby("tissue")[list(S14_BATCH_COLS)].nunique()
print("levels of each batch column within a tissue (1 = the whole tissue in one batch):")
print(_per_tissue.to_string())
_shared = pd.DataFrame({c: {"levels": _b[c].nunique(),
                            "max tissues per level": int(_b.groupby(c)["tissue"].nunique().max()),
                            "levels shared by >1 tissue": int((_b.groupby(c)["tissue"].nunique() > 1).sum()),
                            "tissues in exactly 1 level": int((_per_tissue[c] == 1).sum())} for c in S14_BATCH_COLS}).T
print("\nper batch column:")
print(_shared.to_string())
record("s14.n_tissues_one_plate", int((_per_tissue["RNA_extr_plate_ID"] == 1).sum()), "14")
record("s14.n_plates", _b["RNA_extr_plate_ID"].nunique(), "14")
record("s14.n_flowcells", _b["Seq_flowcell_ID"].nunique(), "14")

# %% [markdown]
# **What (2) shows.** Within a tissue, every vial shares one RNA extraction plate, one library batch, one
# flowcell and one sequencing site (the library prep date is the only column that splits some tissues).
# Across tissues the nesting is not one-to-one: a plate or library batch holds at most a couple of
# tissues, while a flowcell carries several tissues at once. So "batch is nested in tissue" means *each
# tissue sits in one batch*, not *each batch holds one tissue* — which is why plate / library batch are
# near-perfect tissue proxies while flowcell alone is not.

# %%
# ---- (3) the other half of the argument: an independently processed dataset (rat BodyMap) ----
# The fingerprint was fit on all MoTrPAC animals and applied to BodyMap (GSE53960), a different lab, different
# library prep and sequencing, so none of the MoTrPAC batches exist there. Adults = 21-week animals.
# (loaded from the post-fix results; accuracy did not change with the conformal quantile fix, which touches set sizes only)
_bm = pd.read_csv(res("12_bodymap", "age_shift_accuracy.csv"))
print(_bm.to_string(index=False, float_format="%.3f"))
_ad = _bm.set_index("stage_weeks").loc[21]
print(f"\nBodyMap adults (21 weeks): accuracy full model {_ad['full']:.3f}, k = 20 panel {_ad['k20']:.3f}, k = 50 panel {_ad['k50']:.3f}")
print(f"vs QC-only on MoTrPAC (all QC features): {s14_summ.loc['all', 'acc_mean']:.3f} — a QC-only model cannot transfer, "
      "because the QC table and its batches do not exist for BodyMap")
record("s14.bodymap_adult_k20", _ad["k20"], "14")
record("s14.bodymap_adult_full", _ad["full"], "14")
# ...and what does NOT transfer: the conformal calibration. Coverage of the adult BodyMap organs (α = 0.10, LAC)
# with the threshold calibrated on 15 held-out MoTrPAC animals, then re-calibrated on 3 or 5 BodyMap animals.
_ct = pd.read_csv(res("12_bodymap", "conformal_transfer.csv"))
_ct = _ct[(_ct["stage_weeks"] == 21) & (_ct["conformal"] == "marginal")].set_index("model")
_rc = pd.read_csv(res("12_bodymap", "recalibration.csv")).set_index(["model", "n_recal"])
s14_bm_cov = pd.DataFrame({
    "accuracy (adults)": [_ad[m] for m in ("k20", "k50", "full")],
    "coverage, MoTrPAC-calibrated": [_ct.loc[m, "coverage_mapped"] for m in ("k20", "k50", "full")],
    "empty-set rate, MoTrPAC-calibrated": [_ct.loc[m, "frac_empty_mapped"] for m in ("k20", "k50", "full")],
    "coverage, recalibrated on 3 BodyMap animals": [_rc.loc[(m, 3), "coverage_recalibrated"] for m in ("k20", "k50", "full")],
    "coverage, recalibrated on 5": [_rc.loc[(m, 5), "coverage_recalibrated"] for m in ("k20", "k50", "full")],
    "set size, recalibrated on 3": [_rc.loc[(m, 3), "set_size_recalibrated"] for m in ("k20", "k50", "full")]},
    index=["k20", "k50", "full"])
display(s14_bm_cov.round(3))
record("s14.bodymap_adult_cov_k20", s14_bm_cov.loc["k20", "coverage, MoTrPAC-calibrated"], "14")
record("s14.bodymap_adult_empty_k20", s14_bm_cov.loc["k20", "empty-set rate, MoTrPAC-calibrated"], "14")
record("s14.bodymap_recal3_cov_k20", s14_bm_cov.loc["k20", "coverage, recalibrated on 3 BodyMap animals"], "14")
for _n in ("X", "m", "p", "f", "tr", "te", "name", "cols"):
    globals().pop(_n, None)

# %% [markdown]
# **The honest reading has two parts, and neither is sufficient alone.**
#
# 1. **Within MoTrPAC, tissue accuracy cannot separate biology from processing.** QC covariates alone
#    identify the tissue, because batch is nested in tissue: each tissue was extracted, library-prepped
#    and sequenced as one batch. A classifier that reached the section-4 fingerprint accuracy *could* in
#    principle be reading the batch. Within-study accuracy — however high, however leakage-safe the
#    folds — therefore does not prove the fingerprint is biological.
# 2. **But the fingerprint transfers to data processed independently.** Rat BodyMap was collected,
#    extracted and sequenced by another lab; none of the MoTrPAC plates, library batches or flowcells
#    exist there, yet the MoTrPAC-trained fingerprint identifies the adult BodyMap organs (the accuracy
#    column above). A batch signature cannot do that; tissue biology can. What does *not* travel is the
#    calibration. With thresholds set on MoTrPAC animals, coverage of the same adult samples falls far
#    below 1 − α, almost all of it through empty sets (the coverage and empty-set columns): the model
#    names the right organ but is less confident on another lab's libraries than on its own. A few
#    BodyMap animals restore coverage at about one class per set (the recalibration columns). So the
#    MoTrPAC *confidence scale* is study-specific, as a batch-influenced quantity would be, while the
#    *ranking* of tissues is not.
#
# Part 1 alone would say "the fingerprint may be batch"; part 2 alone would leave open that the
# within-study numbers are inflated by batch. Together: the fingerprint carries transferable biology,
# and the *within-study* accuracy is not the evidence for it — the transfer is. (BodyMap scores only the
# organs it shares with MoTrPAC, with muscle and brain as super-classes; see section 14.)
