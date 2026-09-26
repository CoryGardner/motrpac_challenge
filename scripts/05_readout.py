#!/usr/bin/env python
"""Phase 04/05 read-out: headline numbers of the baselines and compact panels pulled from the
result CSVs (results/04_baselines/*, results/05_panels/*), with the interpretation for the
hackathon. Run after phases 04, 04b (PROT diagnostic), 05 and 05b (annotation). Writes
results/05_panels/readout_*.csv and appends one REPORT.md section.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

from motrpac import cli, config as C, report

R = C.RESULTS_DIR


def read(path: Path) -> pd.DataFrame | None:
    return pd.read_csv(path) if path.exists() else None


def main() -> None:
    ap = cli.common_parser("Baselines + panels read-out")
    args = ap.parse_args()
    cli.banner("05_readout", args)
    out = cli.outdir("05_panels", args.out)
    body = []

    # ---- baselines ------------------------------------------------------------------------
    rows = []
    for arm, d in (("TRNSCRPT counts, 19 tissues", "TRNSCRPT"), ("PROT norm, 7 tissues", "PROT"),
                   ("METAB norm, 19 tissues", "METAB"), ("METAB norm, 9 core tissues", "METAB_core9")):
        s = read(R / "04_baselines" / d / "summary.csv")
        pf = read(R / "04_baselines" / d / "per_fold.csv")
        if s is None:
            continue
        for _, r in s.iterrows():
            rows.append({"arm": arm, "model": r["model"], "bal_acc": r["balanced_accuracy_mean"], "sd": r["balanced_accuracy_std"],
                         "macro_f1": r["macro_f1_mean"], "log_loss": r.get("log_loss_mean", np.nan),
                         "fit_s_per_fold": r.get("fit_seconds_mean", np.nan),
                         "train/test animals": f"{int(pf['n_train_animals'].iloc[0])}/{int(pf['n_test_animals'].iloc[0])}" if pf is not None and "n_train_animals" in pf else ""})
    base = pd.DataFrame(rows)
    base.to_csv(out / "readout_baselines.csv", index=False)
    pairs = read(R / "04_baselines" / "TRNSCRPT" / "confusable_pairs_logreg_l2.csv")
    diag = read(R / "04_baselines" / "PROT" / "diagnostic_accuracy.csv")
    body.append("### Baselines (T4)\n\n" + report.df_to_md(base, floatfmt=".3f") + "\n\n"
                "Transcripts: all four baselines are within one sd of each other at 0.98–0.995; the tuned L2 logistic "
                "regression is the reference every later model has to beat, and the L1 model buys nothing for 40× the "
                "fit time. Remaining confusions with all genes available (`logreg_l2`, pooled over folds):\n\n"
                + (report.df_to_md(pairs, floatfmt=".2f") if pairs is not None and len(pairs) else "_none_") + "\n\n"
                "Metabolomics separates all tissues with 130 hilicpos metabolites (19 tissues) and perfectly with 435 "
                "features on the 9 core tissues. **Proteomics is not a result**: the 0.99–1.00 accuracy is an artefact "
                + (f"— on fold 0 the NaN pattern alone classifies tissue with accuracy {diag['missingness_indicators_only'].iloc[0]:.2f} "
                   f"and removing per-tissue means leaves {diag['per_tissue_means_removed'].iloc[0]:.2f} (chance {diag['chance_balanced'].iloc[0]:.2f}); "
                   if diag is not None else "") +
                "proteins are quantified plex by plex, plexes are nested in tissue, so which proteins are missing is a "
                "tissue label in disguise, and the per-tissue normalization leaves only residual offsets (see the PROT "
                "diagnostic section).")

    # ---- panels ---------------------------------------------------------------------------
    for arm, d in (("TRNSCRPT", "TRNSCRPT"), ("METAB core 9", "METAB_core9")):
        pc = read(R / "05_panels" / d / "panel_curve.csv")
        if pc is None:
            continue
        agg = pc.groupby("k").agg(bal_acc=("balanced_accuracy", "mean"), sd=("balanced_accuracy", "std"),
                                  macro_f1=("macro_f1", "mean")).reset_index()
        other = read(R / "05_panels" / d / "panel_curve_fclassif.csv")
        if other is not None:
            agg = agg.merge(other.groupby("k")["balanced_accuracy"].mean().rename("bal_acc_fclassif").reset_index(), on="k")
        cons = read(R / "05_panels" / d / "selection_consistency.csv")
        if cons is not None:
            agg = agg.merge(cons, on="k", how="left")
        agg.to_csv(out / f"readout_curve_{d}.csv", index=False)
        k99 = agg.loc[agg["bal_acc"] >= 0.99, "k"].min() if (agg["bal_acc"] >= 0.99).any() else None
        ann = read(R / "05_panels" / d / "candidate_panel_annotated.csv")
        pk = {k: read(R / "05_panels" / d / f"confusable_pairs_k{k}.csv") for k in (20, 30)}
        txt = [f"### Compact panels (T5) — {arm}\n\n" + report.df_to_md(agg, floatfmt=".3f"),
               f"Smallest k with mean balanced accuracy ≥ 0.99: **{k99}** (point estimate; the certified k comes in phase 06). "
               "`frac_in_all_folds` is the share of the features chosen at that k that every fold chose."]
        for k, cp in pk.items():
            if cp is not None:
                txt.append(f"Confused pairs at k={k}: " + (", ".join(f"{r['true']}→{r['predicted']} ({int(r['count'])})" for _, r in cp.head(8).iterrows()) if len(cp) else "none")
                           + f"; {int(cp['count'].sum()) if len(cp) else 0} errors in total.")
        if ann is not None:
            marker_t = set(ann["marker_tissue"])
            n_reg, n_qc = int(ann["risk_T7_regulated"].sum()), int(ann["risk_qc_correlated"].sum())
            flagged = ann[ann["risk_T7_regulated"] | ann["risk_qc_correlated"]]
            qc_col = [c for c in ann.columns if c.startswith("r_") and c.endswith("_in_marker_tissue")]
            txt.append(f"Candidate panel: {len(ann)} features stable in ≥ 80% of animal bootstraps, covering "
                       f"{len(marker_t)} tissues as top marker. Flags: {n_reg} training-regulated in their marker tissue "
                       f"(risk for controls→trained, T7), {n_qc} with |r| ≥ 0.5 to the library QC metric in their marker "
                       "tissue (risk for external validation). Flagged features:\n\n"
                       + (report.df_to_md(flagged[["feature_ID", "gene_symbol", "marker_tissue", "effect_size", "regulated_tissues"] + qc_col],
                                          floatfmt=".2f") if len(flagged) else "_none_"))
        body.append("\n\n".join(txt))

    report.add_section("04/05 · Baselines and panels read-out (real data)", "\n\n".join(body),
                       params={"inputs": "results/04_baselines/*, results/05_panels/*"})
    print("\n\n".join(body)[:3000])
    print(f"wrote {out}/readout_*.csv and the REPORT.md section")


if __name__ == "__main__":
    main()
