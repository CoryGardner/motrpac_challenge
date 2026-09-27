# %% [markdown]
# ## 12. Transcript–protein discordance
#
# **Question.** In the seven tissues with both RNA-seq and proteomics, how often do the two layers
# disagree about a gene's training response — opposite signs, or significant in one layer only — and
# can simple gene covariates predict which genes disagree?
#
# **Why it matters.** Track 2 is about RNA–protein discordance. Before building on it, one needs to know
# what "discordance" mostly *is* in these data. If nearly every discordant gene is "significant in one
# layer, not significant in the other", that is at least partly a statement about statistical power,
# not about biology. And a predictor of discordance is only interesting if it does not use information
# that already defines the label.
#
# **What to look for.** (1) The per-gene RNA–protein correlation across animals, and how much the TMT
# channel (which encodes sex in HEART and LIVER) moves it. (2) The three counts over all
# tissue × sex × time-point contrasts at FDR 0.10: gene × contrast pairs significant in either layer,
# significant in one layer only, and significant in both with opposite signs. (3) The AUROC for
# predicting "one-layer-only" genes **with and without the circular covariate** `is_training_regulated`
# — the consortium's training-regulated flag, derived from the same differential-analysis tables that
# define "significant in one layer" (the script's own caveat).
#
# **Mode.** Correlations, DA counts and the discordance-prediction AUROCs are recomputed live in both modes;
# per-gene tables go to `_outputs/s12/`. Only the per-covariate drop-one importances (a refit per covariate per
# fold, most of this phase's run time) are loaded from `results/09_discordance/` unless `RECOMPUTE = True`;
# the AUROCs themselves are the same fits either way.

# %%
section("12 discordance")
# from scripts/09_discordance.py::main (trimmed: argparse, rho histograms and report.add_section dropped; output → OUT/s12;
# defaults --assay-a TRNSCRPT --assay-b PROT --fdr 0.1, complete protein features; `make discordance` passes no flags)


# [notebook workaround] src/tfp/discordance.py::ptm_site_counts has a function-local relative import
# (`from .io import read_sample_table`) that the builder does not strip (it strips only top-level ones), so the
# pasted version raises ImportError here. Below: the same function, verbatim, minus that one line
# (read_sample_table is already defined bare in this namespace). Reported to the lead as a build.py fix.
def _s12_ptm_site_counts(tissue: str, assay: str, f2g=None, raw_dir=None) -> pd.Series:
    """Number of quantified PTM sites per gene in this tissue (0 when the assay has no table here)."""
    p = (raw_dir or C.RAW_DIR) / "norm" / f"{assay}__{C.tissue_token(tissue)}.csv"
    if not p.exists():
        return pd.Series(dtype=float, name=f"n_{assay.lower()}_sites")
    feat, _ = read_sample_table(p)
    acc = pd.Series([str(f).rsplit("_", 1)[0] for f in feat.index])
    genes = map_to_gene_symbols(acc.unique(), f2g, raw_dir)
    g = acc.map(genes.to_dict())
    return g.dropna().value_counts().rename(f"n_{assay.lower()}_sites")


import inspect as _s12_inspect
D = discordance
if "from .io import" in _s12_inspect.getsource(discordance.ptm_site_counts):
    D = SimpleNamespace(**{**vars(discordance), "ptm_site_counts": _s12_ptm_site_counts})
S12_OUT = OUT / "s12"
S12_OUT.mkdir(parents=True, exist_ok=True)
_A, _B, _FDR, _complete_b, _seed = "TRNSCRPT", "PROT", 0.1, True, C.SEED
pheno = cached("pheno", io.load_pheno)
f2g = io.load_feature_to_gene()
files = io.list_sample_files("norm")
tissues = sorted(set(files[files["assay"] == _A]["tissue"]) & set(files[files["assay"] == _B]["tissue"]))
da_files = {C.tissue_from_token(p.stem.split("__")[1]) for p in C.DA_DIR.glob(f"{_B}__*.csv")} if C.DA_DIR.exists() else set()
ptm_tissues = {a: set(files[files["assay"] == a]["tissue"]) for a in ("PHOSPHO", "ACETYL", "UBIQ")}
try:
    reg = pd.read_csv(C.RAW_DIR / "training_regulated_features.csv", dtype=str)
    regulated = set(io.map_to_gene_symbols(reg["feature_ID"], f2g).dropna())
except Exception:
    regulated = None
chan = D.channel_map(_B)

summary, rates, aucs, imps, chan_rows = [], [], [], [], []
s12_rho = {}
for t in tissues:
    A, B, meta = D.paired_layers(t, _A, _B, pheno, f2g, complete_b=_complete_b)
    corr = D.per_gene_correlation(A, B)
    row = {"tissue": t, "n_animals": len(meta), "n_genes": len(corr), "median_rho": corr["rho"].median(),
           "frac_rho_gt_0.3": float((corr["rho"] > 0.3).mean()), "frac_rho_lt_0": float((corr["rho"] < 0).mean())}
    # channel residualization (unsupervised) → ρ after
    if chan is not None:
        A2, B2, _ = D.paired_layers(t, _A, _B, pheno, f2g, complete_b=_complete_b, residualize_b=chan)
        corr2 = D.per_gene_correlation(A2, B2).rename(columns={"rho": "rho_after_channel"})
        cc = corr.merge(corr2[["gene", "rho_after_channel"]], on="gene")
        cc.to_csv(S12_OUT / f"per_gene_rho_{C.tissue_token(t)}.csv", index=False)
        chan_rows.append({"tissue": t, "median_rho_before": cc["rho"].median(), "median_rho_after": cc["rho_after_channel"].median(),
                          "frac_rho_lt_0_before": float((cc["rho"] < 0).mean()), "frac_rho_lt_0_after": float((cc["rho_after_channel"] < 0).mean()),
                          "frac_genes_sign_flip": float((np.sign(cc["rho"]) != np.sign(cc["rho_after_channel"])).mean())})
        row["median_rho_after_channel"] = cc["rho_after_channel"].median()
        row["frac_rho_lt_0_after_channel"] = float((cc["rho_after_channel"] < 0).mean())
    s12_rho[t] = corr["rho"].to_numpy()
    ptm = {}
    for a in ("PHOSPHO", "ACETYL", "UBIQ"):
        if t in ptm_tissues[a]:
            ptm[f"n_{a.lower()}_sites"] = D.ptm_site_counts(t, a, f2g)
    feat = D.discordance_features(A, B, corr, regulated, ptm)
    if t in da_files:
        dd = D.da_discordance(t, _A, _B, _FDR, f2g)
        dd.to_csv(S12_OUT / f"da_discordance_{C.tissue_token(t)}.csv", index=False)
        r = dd.groupby(["sex", "comparison_group"]).agg(
            n_genes=("gene", "size"), n_sig_either=("sig_a", lambda s: int((s | dd.loc[s.index, "sig_b"]).sum())),
            n_concordant=("concordant", "sum"), n_sign_discordant=("discordant_sign", "sum"),
            n_one_layer_only=("discordant_sig", "sum")).reset_index()
        r["tissue"] = t
        rates.append(r)
        gene_flag = dd.groupby("gene").agg(any_sign_discordant=("discordant_sign", "any"),
                                            any_one_layer=("discordant_sig", "any")).reset_index()
        feat = feat.merge(gene_flag, on="gene", how="left")
        feat["any_sign_discordant"] = feat["any_sign_discordant"].fillna(False).astype(int)
        feat["any_one_layer"] = feat["any_one_layer"].fillna(False).astype(int)
        fcols = [c for c in ("mean_a", "var_a", "mean_b", "var_b", "prot_abundance_rank", "rho", "is_training_regulated")
                 if c in feat.columns] + list(ptm)
        for target in ("any_sign_discordant", "any_one_layer"):
            # the training-regulated flag is derived from the same DA tables that define the target
            # ('significant in one layer'), so it is reported with and without that flag
            for tag, cols in (("", fcols), (" (no regulated flag)", [c for c in fcols if c != "is_training_regulated"])):
                # [notebook] importance only when RECOMPUTE (same folds and fits; importance adds drop-one refits)
                if RECOMPUTE:
                    res_, imp = D.predict_discordance(feat, target, cols, seed=_seed, importance=True)
                else:
                    res_, imp = D.predict_discordance(feat, target, cols, seed=_seed, importance=False), []
                res_["tissue"], res_["target"] = t, target + tag
                aucs.append(res_)
                if len(imp):
                    imp["tissue"], imp["target"] = t, target + tag
                    imps.append(imp)
        row["frac_genes_any_sign_discordant"] = float(feat["any_sign_discordant"].mean())
        row["frac_genes_any_one_layer"] = float(feat["any_one_layer"].mean())
    feat["rho_low"] = (feat["rho"] < 0.1).astype(int)
    res_ = D.predict_discordance(feat, "rho_low", [c for c in ("mean_a", "var_a", "mean_b", "var_b", "prot_abundance_rank",
                                                                "is_training_regulated") if c in feat.columns] + list(ptm), seed=_seed)
    res_["tissue"], res_["target"] = t, "rho_low(<0.1)"
    aucs.append(res_)
    summary.append(row)
    print(f"  {t}: {len(corr)} genes, {len(meta)} animals, median rho={row['median_rho']:.2f}"
          + (f" → {row['median_rho_after_channel']:.2f} after channel" if "median_rho_after_channel" in row else ""))
s12_summ = pd.DataFrame(summary)
s12_rates = pd.concat(rates, ignore_index=True)
s12_au = pd.concat(aucs, ignore_index=True)
s12_imp = (pd.concat(imps, ignore_index=True) if RECOMPUTE and imps
           else pd.read_csv(res("09_discordance", "covariate_importance.csv")))   # default mode: loaded
for _n, _df in (("summary", s12_summ), ("discordance_rates", s12_rates), ("prediction_auroc", s12_au), ("covariate_importance", s12_imp)):
    _df.to_csv(S12_OUT / f"{_n}.csv", index=False)
# keep the namespace clean: the verbatim loop used generic names (the library alias stays `discordance`;
# note `res_` above — the script's `res` would have shadowed the notebook's res() helper)
for _n in ("D", "f2g", "files", "tissues", "da_files", "ptm_tissues", "reg", "regulated", "chan", "summary", "rates",
           "aucs", "imps", "chan_rows", "t", "A", "B", "meta", "corr", "A2", "B2", "corr2", "cc", "row", "ptm", "a", "feat", "dd",
           "r", "gene_flag", "fcols", "target", "tag", "cols", "imp", "res_"):
    globals().pop(_n, None)

# %%
# ---- (1) per-gene RNA–protein Spearman ρ across animals, before and after removing the TMT channel means ----
_pub = pd.read_csv(res("09_discordance", "summary.csv"))
_cols = ["tissue", "n_animals", "n_genes", "median_rho", "median_rho_after_channel", "frac_rho_lt_0", "frac_rho_gt_0.3"]
print(s12_summ[_cols].to_string(index=False, float_format="%.3f"))
_d = (s12_summ.set_index("tissue")[["median_rho", "median_rho_after_channel", "n_genes"]]
      - _pub.set_index("tissue")[["median_rho", "median_rho_after_channel", "n_genes"]]).abs().max()
print("\nrecomputed vs published summary.csv, max |difference|:", _d.round(12).to_dict())
fig, ax = plt.subplots(figsize=(7, 3.2))
for t, r in s12_rho.items():
    ax.hist(r, bins=60, histtype="step", density=True, label=t)
ax.axvline(0, color="grey", lw=0.8, ls=":")
ax.set_xlabel("per-gene Spearman ρ (RNA vs protein, across animals)")
ax.set_ylabel("density")
ax.legend(fontsize=7, ncol=2, frameon=False)
ax.set_title("RNA–protein correlation per gene, by tissue")
fig.tight_layout()
plt.show()
record("s12.median_rho_WATSC", s12_summ.set_index("tissue").loc["WAT-SC", "median_rho"], "12")
record("s12.median_rho_CORTEX", s12_summ.set_index("tissue").loc["CORTEX", "median_rho"], "12")

# %% [markdown]
# **What (1) shows.** Across animals, the RNA and protein levels of a gene are only weakly correlated
# in most tissues (median ρ near zero in CORTEX, highest in WAT-SC), with a large fraction of negative
# correlations. Removing the TMT channel means moves ρ materially only in HEART and LIVER — the two
# tissues where the channel encodes sex, so there the "correction" also removes a real biological
# (sex) source of concordance. Read those two tissues separately.

# %%
# ---- (2) DA-based discordance: every tissue × sex × time-point contrast at FDR 0.10, pooled ----
_tot = s12_rates[["n_sig_either", "n_concordant", "n_sign_discordant", "n_one_layer_only"]].sum()
_pubr = pd.read_csv(res("09_discordance", "discordance_rates.csv"))
_pubt = _pubr[["n_sig_either", "n_concordant", "n_sign_discordant", "n_one_layer_only"]].sum()
print(f"contrasts: {len(s12_rates)} (tissue × sex × time point), genes per tissue {s12_rates['n_genes'].min()}–{s12_rates['n_genes'].max()}")
print(pd.DataFrame({"recomputed": _tot, "published": _pubt}).to_string())
print(f"\none-layer-only share of pairs significant in either layer: {_tot['n_one_layer_only'] / _tot['n_sig_either']:.3f}")
_sd = s12_rates[s12_rates["n_sign_discordant"] > 0][["tissue", "sex", "comparison_group", "n_sig_either", "n_concordant", "n_sign_discordant"]]
print("contrasts holding the sign-discordant pairs:\n" + _sd.to_string(index=False))
print(f"pairs from 8-week contrasts: {int(s12_rates.loc[s12_rates['comparison_group'] == '8w', 'n_sig_either'].sum())} of {int(_tot['n_sig_either'])}")
record("s12.n_sig_either", _tot["n_sig_either"], "12")
record("s12.n_one_layer_only", _tot["n_one_layer_only"], "12")
record("s12.n_sign_discordant", _tot["n_sign_discordant"], "12")
record("s12.n_concordant", _tot["n_concordant"], "12")

# %% [markdown]
# **What (2) shows.** These are the three counts quoted in the findings: gene × contrast pairs
# significant in RNA or protein; of those, the ones significant in one layer only; and the ones
# significant in both with opposite signs. Nearly all "discordance" is one-layer-only significance;
# opposite-sign pairs are a handful, all in male 8-week HEART and LIVER — the channel = sex tissues —
# where the denominators are too small to carry a rate.
#
# **What it does not show.** That one-layer-only genes are biologically discordant. With few animals
# per contrast, "significant in one layer, not in the other" is what two noisy tests of the same true
# effect produce routinely; the findings read it as mostly statistical power (an interpretation, not a
# measurement). The counts also inherit the consortium's DA thresholds and the c1.0 release.

# %%
# ---- (3) can gene covariates predict "one-layer-only" genes? With and without the circular flag ----
# gene-grouped 5-fold CV inside each tissue (library function discordance.predict_discordance); pooled like SUMMARY.md
_pubau = pd.read_csv(res("09_discordance", "prediction_auroc.csv"))
_agg = lambda a: a.groupby(["target", "model"])["auroc"].agg(["mean", "std", "count"])
s12_auc_tab = _agg(s12_au).join(_agg(_pubau)[["mean"]].rename(columns={"mean": "published_mean"}))
print(s12_auc_tab.to_string(float_format="%.3f"))
_w = s12_auc_tab.loc[("any_one_layer", "logreg"), "mean"]
_wo = s12_auc_tab.loc[("any_one_layer (no regulated flag)", "logreg"), "mean"]
print(f"\nlogreg AUROC, any_one_layer: with is_training_regulated {_w:.3f}  →  without it {_wo:.3f}  (drop {_w - _wo:.3f})")
_rw = s12_auc_tab.loc[("any_one_layer", "rf"), "mean"]
_rwo = s12_auc_tab.loc[("any_one_layer (no regulated flag)", "rf"), "mean"]
print(f"rf     AUROC, any_one_layer: with {_rw:.3f}  →  without {_rwo:.3f}")
# which covariates carry the signal (mean drop-one ΔAUROC of the logistic model over tissues and folds)
# (importances: recomputed when RECOMPUTE, else loaded from results/09_discordance/covariate_importance.csv)
_top = (s12_imp.groupby(["target", "covariate"])["logreg_drop_one_dAUROC"].mean().reset_index()
        .sort_values(["target", "logreg_drop_one_dAUROC"], ascending=[True, False]).groupby("target").head(3))
print("\ntop covariates by drop-one ΔAUROC (logreg):\n" + _top.to_string(index=False, float_format="%+.3f"))
record("s12.auroc_one_layer_logreg_with_flag", _w, "12")
record("s12.auroc_one_layer_logreg_without_flag", _wo, "12")
record("s12.auroc_one_layer_rf_with_flag", _rw, "12")
record("s12.auroc_one_layer_rf_without_flag", _rwo, "12")

# %% [markdown]
# **What (3) shows.** Gene covariates predict which genes are one-layer-only well above chance, and the
# circular covariate `is_training_regulated` adds to that — as it must, because the flag comes from the
# same consortium DA tables that define the label. The rows *without* the flag are the honest number;
# the gap between the two is the part of the "predictability" that was circular. Without the flag, the
# transcript variance and mean and the RNA–protein ρ carry most of the signal — covariates of
# measurement power (how noisy the transcript is, how well the layers agree at baseline), which fits the
# "mostly power" reading of (2).
#
# **What it does not show.** A biological mechanism for discordance. The predictor is fit on
# tissue-specific gene sets with gene-grouped folds, but its covariates describe measurability, not
# regulation; the sign-discordant target has too few positives to model at all (the "too few
# positives" rows in `prediction_auroc.csv`).
