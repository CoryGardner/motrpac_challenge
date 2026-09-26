# %% [markdown]
# ## 3. Proteomics diagnostic: why cross-tissue PROT values are not comparable
#
# **Question.** Section 2 found no tissue axis in the stacked proteomics PCA, yet a cross-tissue
# proteomics classifier still predicts tissue well. What is it using?
#
# **Why it matters.** If the classifier lives on how each tissue was *processed* (its TMT plexes, its
# missing-value pattern, residual offsets of the per-tissue normalization) rather than on protein
# abundance, then a "proteomic tissue fingerprint" is an artefact, and proteomics has to be used within
# tissue only.
#
# **What to look for.** Three test accuracies on the same held-out animals, against chance
# (one over the number of tissues):
# 1. the tuned L2 logistic-regression pipeline, fitted as in section 4 (`model_as_fitted`);
# 2. a model that sees **only which values are missing** (`missingness_indicators_only`);
# 3. a model after each tissue's feature means are subtracted (`per_tissue_means_removed`). This uses
#    the tissue labels to centre the data, so it is a diagnostic of what is left once mean offsets are
#    erased — never an evaluation.
#
# If (1) and (2) are high and (3) falls toward chance, the signal is processing, not biology. Then a
# second check: inside HEART and LIVER, is the TMT channel a sample was run in confounded with its sex?
#
# Settings are those of `scripts/04_prot_diagnostic.py --assay PROT` (`make prot-diagnostic`):
# normalized PROT tables, inner feature join, animal-grouped stratified 5-fold split with the pipeline
# seed, **outer fold 0 only** (as published), prefilter 5,000, full tuning grid.

# %%
section("3 Proteomics diagnostic")
from sklearn.impute import SimpleImputer
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler


# from scripts/04_prot_diagnostic.py::r2_cat (verbatim)
def r2_cat(x: np.ndarray, cat: np.ndarray) -> float:
    ok = ~np.isnan(x)
    x, cat = x[ok], cat[ok]
    if len(x) < 3 or x.var() == 0 or len(np.unique(cat)) < 2:
        return np.nan
    between = sum(((x[cat == c].mean() - x.mean()) ** 2) * (cat == c).sum() for c in np.unique(cat))
    return float(between / (x.var() * len(x)))


# from scripts/04_prot_diagnostic.py::main (trimmed: argparse/report dropped; defaults of
# common_parser inlined: n_splits=5, prefilter=5000, seed=C.SEED, quick=False, no --complete-features)
S03_OUT = OUT / "s03"
S03_OUT.mkdir(parents=True, exist_ok=True)
om = io.stack_tissues("PROT", tissues=None, source="norm", join="inner", pheno=pheno,
                      complete=False, drop_incomplete_samples=None)
X = om.X.to_numpy(dtype=float)
y = om.meta["tissue"].astype(str).to_numpy()
g = om.groups()
feat = np.asarray(om.X.columns)
classes = sorted(np.unique(y))
tr, te = next(iter(grouped_kfold(om.meta, "tissue", 5, C.SEED)))
print(f"fold 0: {len(np.unique(g[tr]))} training / {len(np.unique(g[te]))} test animals; "
      f"{len(tr)} / {len(te)} samples; {len(classes)} tissues; overall fraction of NaN values {np.isnan(X).mean():.3f}")

# 1. tissue R² per protein, all samples (descriptive)
tr2 = pd.Series([r2_cat(X[:, j], y) for j in range(X.shape[1])], index=feat, name="tissue_R2")
print("across-tissue R² of tissue per protein, quantiles:", tr2.quantile([0.5, 0.9, 0.99]).round(3).to_dict(),
      f"; proteins with R² > 0.5: {int((tr2 > 0.5).sum())} of {len(tr2)}")

# %%
# 2. the fitted model (tuned inside the training fold only) and its largest coefficients
est = models.fit_tuned("logreg_l2", X[tr], y[tr], g[tr], prefilter=5000, quick=False, seed=C.SEED)
clf = est.named_steps["clf"]
kept = est.named_steps["prefilter"].idx_
coef = pd.DataFrame(clf.coef_, index=list(clf.classes_), columns=feat[kept])
meta_path = C.META_DIR / "PROT.csv"
channel = None
if meta_path.exists():
    meta = pd.read_csv(meta_path, dtype=str, low_memory=False).drop_duplicates("viallabel").set_index("viallabel")
    if "tmt11_channel" in meta.columns:
        channel = meta["tmt11_channel"].reindex(om.X.index).to_numpy()
sex = om.meta["sex"].astype(str).to_numpy()
rows = []
for t in classes:
    top = coef.loc[t].abs().sort_values(ascending=False).head(3)
    for f in top.index:
        j = int(np.flatnonzero(feat == f)[0])
        mk = y == t
        nan_by_t = pd.Series(np.isnan(X[:, j])).groupby(y).mean()
        rows.append({"tissue": t, "feature_ID": f, "coef": float(coef.loc[t, f]),
                     "tissue_R2_all_samples": float(tr2[f]),
                     "R2_sex_within_tissue": r2_cat(X[mk, j], sex[mk]),
                     "R2_channel_within_tissue": r2_cat(X[mk, j], channel[mk]) if channel is not None else np.nan,
                     "frac_nan_in_tissue": float(nan_by_t[t]), "max_frac_nan_other_tissue": float(nan_by_t.drop(t).max())})
top_df = pd.DataFrame(rows)
top_df.to_csv(S03_OUT / "diagnostic_top_coefs.csv", index=False)
print("largest-|coefficient| proteins per tissue (3 each):")
display(top_df.round(3))

# %%
# 3. fold-0 accuracies: as fitted / missingness only / per-tissue means removed (label-using)
# from scripts/04_prot_diagnostic.py::main (verbatim)
acc = {"model_as_fitted": float(np.mean(est.predict(X[te]) == y[te]))}
M = np.isnan(X).astype(float)
if M[tr].std(axis=0).max() > 0:
    miss = LogisticRegression(C=0.1, max_iter=3000).fit(M[tr], y[tr])
    acc["missingness_indicators_only"] = float(np.mean(miss.predict(M[te]) == y[te]))
else:
    acc["missingness_indicators_only"] = np.nan
Xi = SimpleImputer(strategy="median").fit(X[tr]).transform(X)
Xc = Xi.copy()
for t in classes:  # centre every tissue on its TRAINING-fold feature means (uses labels: diagnostic only)
    mu = Xi[tr][y[tr] == t].mean(axis=0)
    Xc[y == t] -= mu
sc = StandardScaler().fit(Xc[tr])
cen = LogisticRegression(C=0.1, max_iter=3000).fit(sc.transform(Xc[tr]), y[tr])
acc["per_tissue_means_removed"] = float(np.mean(cen.predict(sc.transform(Xc[te])) == y[te]))
acc["chance_balanced"] = 1.0 / len(classes)
acc_df = pd.DataFrame([acc])
acc_df.to_csv(S03_OUT / "diagnostic_accuracy.csv", index=False)

s03_pub = pd.read_csv(res("04_baselines", "PROT", "diagnostic_accuracy.csv"))
print(f"fold-0 test accuracy ({len(classes)} tissues; chance = 1/{len(classes)}):")
display(pd.concat([acc_df.assign(run="recomputed"), s03_pub.assign(run="published")]).set_index("run").round(4))
for k in acc:
    record(f"s03.{k}", acc[k], "03")

# %% [markdown]
# **What this shows.** The quantile printout: for most proteins tissue explains only a small share of
# the variance across the stacked samples — the per-tissue reference ratios and per-sample centring
# remove most mean differences — but a minority keep tissue-specific offsets. The accuracy table: the fitted model classifies held-out
# animals' tissues far above chance; the missing-value pattern *alone* does at least as well; and once
# each tissue's mean offsets are removed, accuracy falls to near chance. The classifier therefore rests
# on how each tissue was processed (which proteins each tissue's plexes detected, and residual offsets
# of the per-tissue normalization), not on a proteomic tissue fingerprint. The top-coefficient table
# shows the same from the model's side: its most-used proteins are among the minority with large
# tissue offsets, several are partly missing in their own or in other tissues, and their within-tissue
# variation mostly tracks TMT channel more than sex.
#
# **What it does not show.** This is one fold of five, as in the pipeline's diagnostic, so its
# accuracies have no fold-to-fold spread; the full 5-fold PROT baseline is in section 4's results. The
# mean-removal model uses the labels to centre and is not a valid classifier. It does not say that
# proteomics lacks tissue biology — only that these per-tissue-normalized tables cannot show it.

# %%
# ---- TMT channel vs sex within each proteomics tissue ----
# from scripts/03_eda.py::shared_with (verbatim). Association of a covariate with a design factor on a
# 0–1 scale (Cramér's V² for a categorical covariate); 1 = the covariate is the factor in this tissue.
def shared_with(cov, factor) -> float:
    cov = pd.Series(np.asarray(cov)).reset_index(drop=True)
    fac = pd.Series(np.asarray(factor)).astype(str).reset_index(drop=True)
    ok = cov.notna().to_numpy()
    if ok.sum() < 3 or cov[ok].nunique() < 2 or fac[ok].nunique() < 2:
        return np.nan
    if pd.api.types.is_numeric_dtype(cov):
        return float(r2_categorical(cov[ok].to_numpy(dtype=float)[:, None], fac[ok])[0])
    ct = pd.crosstab(cov[ok].astype(str), fac[ok]).to_numpy(dtype=float)
    n = ct.sum()
    expected = np.outer(ct.sum(1), ct.sum(0)) / n
    chi2 = ((ct - expected) ** 2 / expected).sum()
    return float(chi2 / (n * (min(ct.shape) - 1)))


# recompute R2_tmt11_channel~sex per tissue (as scripts/03_eda.py::main does on each tissue's own table),
# next to the published within-tissue table
s03_wt = pd.read_csv(res("03_eda", "within_tissue_pca_PROT.csv")).set_index("tissue")
s03_rows = []
for t in s03_wt.index:
    sub = io.load_norm("PROT", t, pheno)
    ch = meta["tmt11_channel"].reindex(sub.X.index)
    s03_rows.append({"tissue": t, "n": sub.n_samples, "n_channels": int(ch.nunique()),
                     "R2_channel~sex (recomputed)": shared_with(ch, sub.meta["sex"]),
                     "R2_channel~sex (published)": s03_wt.loc[t, "R2_tmt11_channel~sex"],
                     "max_R2_sex_PC1-5 (published)": s03_wt.loc[t, "max_R2_sex_PC1-5"],
                     "max_R2_channel|design_PC1-5 (published)": s03_wt.loc[t, "max_R2_tmt11_channel|design_PC1-5"]})
s03_ch = pd.DataFrame(s03_rows).set_index("tissue")
display(s03_ch.round(3))
# the sex-by-channel crosstab in the two tissues where they coincide
for t in ("HEART", "LIVER"):
    sub = io.load_norm("PROT", t, pheno)
    print(f"\n{t}: samples per TMT channel x sex")
    display(pd.crosstab(meta["tmt11_channel"].reindex(sub.X.index).to_numpy(), sub.meta["sex"].to_numpy(),
                        rownames=["tmt11_channel"], colnames=["sex"]).T)
    record(f"s03.{t}_R2_channel_sex", s03_ch.loc[t, "R2_channel~sex (recomputed)"], "03")

# %% [markdown]
# **What this shows.** In HEART and LIVER the TMT channel and sex are the same variable (association 1:
# each channel holds one sex only, as the crosstabs show), while in the other proteomics tissues they
# are unrelated. So in those two tissues any sex difference in protein abundance is also a channel
# difference and cannot be attributed to biology from these data alone. The last published column is
# the share of a within-tissue PC that channel still explains after the sex × group design means are
# removed — a purely technical component in several tissues.
#
# **What it does not show.** Channel–sex confounding does not mean the sex effects in HEART and LIVER
# are false; it means this design cannot separate them from a channel effect. Correcting for channel in
# those tissues would remove the sex effect with it.
