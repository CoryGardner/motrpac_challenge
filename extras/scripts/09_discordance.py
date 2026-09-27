#!/usr/bin/env python
"""Phase 09 — transcript–protein discordance (track 2): how much, where, and is it predictable?

Per tissue with both layers (PROT restricted to proteins with no missing value in the tissue):
  1. per-gene Spearman ρ across matched animals, before and after an unsupervised TMT-channel
     residualization of the protein layer (per-channel means subtracted, no labels)
  2. DA-based discordance per sex × time point (sign disagreement / one-layer-only significance)
  3. gene-grouped CV: predict 'DA-discordant gene' from transcript mean/variance, protein
     abundance rank, ρ, the training-regulated flag and PTM site counts (PHOSPHO everywhere,
     ACETYL/UBIQ where measured), with per-covariate importance

Outputs (results/09_discordance/): per_gene_rho_<tissue>.csv, rho_hist_<tissue>.png,
da_discordance_<tissue>.csv, discordance_rates.csv, channel_effect.csv, prediction_auroc.csv,
covariate_importance.csv, summary.csv
"""
from __future__ import annotations

import numpy as np
import pandas as pd

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))  # extras/ holds tfp_extras
from tfp import cli, config as C, io, plots, report  # noqa: E402
from tfp_extras import discordance as D  # noqa: E402


def main() -> None:
    ap = cli.common_parser("Transcript–protein discordance")
    ap.add_argument("--assay-a", default="TRNSCRPT")
    ap.add_argument("--assay-b", default="PROT")
    ap.add_argument("--fdr", type=float, default=0.1)
    ap.add_argument("--no-complete-b", action="store_true", help="keep layer-b features with missing values")
    args = ap.parse_args()
    cli.banner("09_discordance", args)
    out = cli.outdir("09_discordance", args.out)
    pheno = io.load_pheno()
    f2g = io.load_feature_to_gene()
    files = io.list_sample_files("norm")
    tissues = cli.parse_tissues(args.tissues) or sorted(
        set(files[files["assay"] == args.assay_a]["tissue"]) & set(files[files["assay"] == args.assay_b]["tissue"]))
    da_files = {C.tissue_from_token(p.stem.split("__")[1]) for p in C.DA_DIR.glob(f"{args.assay_b}__*.csv")} if C.DA_DIR.exists() else set()
    ptm_tissues = {a: set(files[files["assay"] == a]["tissue"]) for a in ("PHOSPHO", "ACETYL", "UBIQ")}
    try:
        reg = pd.read_csv(C.RAW_DIR / "training_regulated_features.csv", dtype=str)
        regulated = set(io.map_to_gene_symbols(reg["feature_ID"], f2g).dropna())
    except Exception:
        regulated = None
    chan = D.channel_map(args.assay_b)
    complete_b = not args.no_complete_b

    summary, rates, aucs, imps, chan_rows, body = [], [], [], [], [], []
    for t in tissues:
        A, B, meta = D.paired_layers(t, args.assay_a, args.assay_b, pheno, f2g, complete_b=complete_b)
        corr = D.per_gene_correlation(A, B)
        corr.to_csv(out / f"per_gene_rho_{C.tissue_token(t)}.csv", index=False)
        fig = plots.histogram(corr["rho"].to_numpy(), f"{t}: per-gene Spearman ρ ({args.assay_a} vs {args.assay_b}, n={len(meta)} animals)",
                              "Spearman ρ across animals", out / f"rho_hist_{C.tissue_token(t)}.png")
        row = {"tissue": t, "n_animals": len(meta), "n_genes": len(corr), "median_rho": corr["rho"].median(),
               "frac_rho_gt_0.3": float((corr["rho"] > 0.3).mean()), "frac_rho_lt_0": float((corr["rho"] < 0).mean())}
        # channel residualization (unsupervised) → ρ after
        if chan is not None:
            A2, B2, _ = D.paired_layers(t, args.assay_a, args.assay_b, pheno, f2g, complete_b=complete_b, residualize_b=chan)
            corr2 = D.per_gene_correlation(A2, B2).rename(columns={"rho": "rho_after_channel"})
            cc = corr.merge(corr2[["gene", "rho_after_channel"]], on="gene")
            cc.to_csv(out / f"per_gene_rho_{C.tissue_token(t)}.csv", index=False)
            chan_rows.append({"tissue": t, "median_rho_before": cc["rho"].median(), "median_rho_after": cc["rho_after_channel"].median(),
                              "frac_rho_lt_0_before": float((cc["rho"] < 0).mean()), "frac_rho_lt_0_after": float((cc["rho_after_channel"] < 0).mean()),
                              "frac_rho_gt_0.3_before": float((cc["rho"] > 0.3).mean()), "frac_rho_gt_0.3_after": float((cc["rho_after_channel"] > 0.3).mean()),
                              "median_abs_change": float((cc["rho_after_channel"] - cc["rho"]).abs().median()),
                              "frac_genes_sign_flip": float((np.sign(cc["rho"]) != np.sign(cc["rho_after_channel"])).mean())})
            row["median_rho_after_channel"] = cc["rho_after_channel"].median()
            row["frac_rho_lt_0_after_channel"] = float((cc["rho_after_channel"] < 0).mean())
        ptm = {}
        for a in ("PHOSPHO", "ACETYL", "UBIQ"):
            if t in ptm_tissues[a]:
                ptm[f"n_{a.lower()}_sites"] = D.ptm_site_counts(t, a, f2g)
        feat = D.discordance_features(A, B, corr, regulated, ptm)
        if t in da_files:
            dd = D.da_discordance(t, args.assay_a, args.assay_b, args.fdr, f2g)
            dd.to_csv(out / f"da_discordance_{C.tissue_token(t)}.csv", index=False)
            r = dd.groupby(["sex", "comparison_group"]).agg(
                n_genes=("gene", "size"), n_sig_either=("sig_a", lambda s: int((s | dd.loc[s.index, "sig_b"]).sum())),
                n_concordant=("concordant", "sum"), n_sign_discordant=("discordant_sign", "sum"),
                n_one_layer_only=("discordant_sig", "sum")).reset_index()
            r["tissue"] = t
            r["frac_sign_discordant_of_both_sig"] = r["n_sign_discordant"] / (r["n_sign_discordant"] + r["n_concordant"]).replace(0, np.nan)
            r["frac_one_layer_only_of_sig_either"] = r["n_one_layer_only"] / r["n_sig_either"].replace(0, np.nan)
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
                    res, imp = D.predict_discordance(feat, target, cols, seed=args.seed, importance=True)
                    res["tissue"], res["target"] = t, target + tag
                    aucs.append(res)
                    if len(imp):
                        imp["tissue"], imp["target"] = t, target + tag
                        imps.append(imp)
            row["frac_genes_any_sign_discordant"] = float(feat["any_sign_discordant"].mean())
            row["frac_genes_any_one_layer"] = float(feat["any_one_layer"].mean())
        feat.to_csv(out / f"gene_features_{C.tissue_token(t)}.csv", index=False)
        feat["rho_low"] = (feat["rho"] < 0.1).astype(int)
        res = D.predict_discordance(feat, "rho_low", [c for c in ("mean_a", "var_a", "mean_b", "var_b", "prot_abundance_rank", "is_training_regulated") if c in feat.columns] + list(ptm), seed=args.seed)
        res["tissue"], res["target"] = t, "rho_low(<0.1)"
        aucs.append(res)
        summary.append(row)
        body.append(report.figure_md(fig, f"{t}: distribution of per-gene RNA–protein correlation"))
        print(f"  {t}: {len(corr)} genes, {len(meta)} animals, median rho={row['median_rho']:.2f}"
              + (f" → {row['median_rho_after_channel']:.2f} after channel" if "median_rho_after_channel" in row else ""))

    summ = pd.DataFrame(summary)
    summ.to_csv(out / "summary.csv", index=False)
    body.insert(0, f"Per-tissue sample-level discordance ({args.assay_b} restricted to features with no missing value in the tissue):\n\n"
                + report.df_to_md(summ, floatfmt=".2f"))
    if chan_rows:
        ce = pd.DataFrame(chan_rows)
        ce.to_csv(out / "channel_effect.csv", index=False)
        body.append("TMT channel residualized from the protein layer (per-channel means subtracted within tissue, no labels): "
                    "per-gene ρ before vs after:\n\n" + report.df_to_md(ce, floatfmt=".3f")
                    + "\n\nRead: the change in the median ρ and in the fraction of negative ρ is the share of the apparent "
                    "discordance that channel structure explains. In HEART and LIVER channel encodes sex, so residualizing "
                    "channel there also removes the sex difference from the protein layer — a biological source of concordance — "
                    "which is why those two tissues have to be read separately.")
    if rates:
        rt = pd.concat(rates, ignore_index=True)
        rt.to_csv(out / "discordance_rates.csv", index=False)
        body.append(f"DA-based discordance per sex × time point (FDR {args.fdr}):\n\n" + report.df_to_md(rt, floatfmt=".2f"))
        per_sex = rt.groupby(["tissue", "sex"]).agg(sign_discordant_rate=("frac_sign_discordant_of_both_sig", "mean"),
                                                    one_layer_rate=("frac_one_layer_only_of_sig_either", "mean"),
                                                    n_sig_either=("n_sig_either", "sum")).reset_index()
        per_sex["channel_encodes_sex"] = per_sex["tissue"].isin(["HEART", "LIVER"])
        per_sex.to_csv(out / "discordance_rates_per_sex.csv", index=False)
        grp = per_sex.groupby("channel_encodes_sex")[["sign_discordant_rate", "one_layer_rate"]].mean()
        body.append("Per-sex discordance rates, HEART and LIVER (TMT channel = sex) vs the other tissues:\n\n"
                    + report.df_to_md(per_sex, floatfmt=".3f") + "\n\nMeans: " + "; ".join(
                        f"{'channel=sex' if k else 'other'}: sign-discordant {v['sign_discordant_rate']:.3f}, one-layer-only {v['one_layer_rate']:.3f}"
                        for k, v in grp.iterrows())
                    + ". If HEART/LIVER are outliers here, the channel–sex confound (Q22) inflates or deflates their protein DA and the "
                    "discordance rates inherit it.")
    if aucs:
        au = pd.concat(aucs, ignore_index=True)
        au.to_csv(out / "prediction_auroc.csv", index=False)
        au_s = au.groupby(["tissue", "target", "model"])["auroc"].agg(["mean", "std"]).reset_index()
        body.append("Predicting discordance from gene covariates (gene-grouped CV AUROC):\n\n" + report.df_to_md(au_s, floatfmt=".3f"))
    if imps:
        im = pd.concat(imps, ignore_index=True)
        im.to_csv(out / "covariate_importance.csv", index=False)
        im_s = im.groupby(["target", "covariate"]).agg(logreg_abs_coef=("logreg_abs_coef", "mean"),
                                                        logreg_drop_one_dAUROC=("logreg_drop_one_dAUROC", "mean"),
                                                        rf_importance=("rf_importance", "mean")).reset_index().sort_values(["target", "logreg_drop_one_dAUROC"], ascending=[True, False])
        body.append("Which covariates carry the signal (mean over tissues; drop-one ΔAUROC = loss when the covariate is removed "
                    "from the logistic model):\n\n" + report.df_to_md(im_s, floatfmt=".3f"))
    body.append("Caveat: `is_training_regulated` comes from the consortium's DA tables, the same tables that define 'significant in "
                "one layer', so its predictive power for `any_one_layer` is partly circular; the '(no regulated flag)' rows are the "
                "honest version. If AUROC ≈ 0.5 the covariates carry no signal — a real result.")
    report.add_section("09 · Transcript–protein discordance", "\n\n".join(body),
                       params={"assays": f"{args.assay_a}/{args.assay_b}", "fdr": args.fdr, "complete_b": complete_b})
    print(summ.round(3).to_string())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
