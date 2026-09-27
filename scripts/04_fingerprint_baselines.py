#!/usr/bin/env python
"""Phase 04 — tissue classification with tuned simple baselines under animal-grouped CV.

Outputs (results/04_baselines/<assay>/): per_fold.csv, summary.csv, predictions.csv,
confusion_<model>.png, confusable_pairs.csv
"""
from __future__ import annotations

import pandas as pd

from tfp import cli, config as C, io, models, plots, report


def confusable_pairs(conf: pd.DataFrame, top: int = 10) -> pd.DataFrame:
    rows = []
    for a in conf.index:
        for b in conf.columns:
            if a != b and conf.loc[a, b] > 0:
                rows.append({"true": a, "predicted": b, "count": int(conf.loc[a, b]),
                             "frac_of_true": float(conf.loc[a, b] / max(conf.loc[a].sum(), 1))})
    cols = ["true", "predicted", "count", "frac_of_true"]
    df = pd.DataFrame(rows, columns=cols)  # keeps the header even when nothing is confused
    return df.sort_values("count", ascending=False).head(top)


def main() -> None:
    ap = cli.common_parser("Tissue fingerprint baselines")
    ap.add_argument("--models", default=",".join(models.BASELINES))
    ap.add_argument("--label", default="tissue")
    ap.add_argument("--time-one-fold", action="store_true",
                    help="run outer fold 0 only, print fit seconds per model and the projected total, write nothing to REPORT.md")
    args = ap.parse_args()
    cli.banner("04_fingerprint_baselines", args)
    source = cli.resolve_source(args.assay, args.source)
    out = cli.outdir(f"04_baselines/{args.assay}", args.out)
    pheno = io.load_pheno()
    om = io.stack_tissues(args.assay, tissues=cli.parse_tissues(args.tissues), source=source, join="inner", pheno=pheno,
                           complete=args.complete_features,
                          drop_incomplete_samples=args.drop_incomplete_samples)
    kinds = [k.strip() for k in args.models.split(",") if k.strip()]
    n_splits = 3 if args.quick else args.n_splits
    prefilter = min(args.prefilter, 1000) if args.quick else args.prefilter

    if args.time_one_fold:
        res = models.evaluate_cv(om, label=args.label, kinds=kinds, n_splits=n_splits, prefilter=prefilter,
                                 quick=args.quick, seed=args.seed, max_folds=1)
        timing = res.per_fold[["model", "fit_seconds", "n_train_animals", "n_test_animals", "balanced_accuracy"]].copy()
        timing["projected_total_seconds"] = timing["fit_seconds"] * n_splits
        timing.to_csv(out / "timing_fold0.csv", index=False)
        print(timing.to_string(index=False))
        print(f"projected total for {n_splits} folds: {timing['projected_total_seconds'].sum() / 60:.1f} min")
        return
    res = models.evaluate_cv(om, label=args.label, kinds=kinds, n_splits=n_splits, prefilter=prefilter,
                             quick=args.quick, seed=args.seed)
    res.per_fold.to_csv(out / "per_fold.csv", index=False)
    summ = res.summary()
    fit_time = res.per_fold.groupby("model")["fit_seconds"].agg(fit_seconds_mean="mean", fit_seconds_total="sum").reset_index()
    summ = summ.merge(fit_time, on="model")
    summ.to_csv(out / "summary.csv", index=False)
    res.predictions.to_csv(out / "predictions.csv", index=False)
    folds = res.per_fold[res.per_fold["model"] == kinds[0]][["fold", "n_train_animals", "n_test_animals", "n_test"]]
    print(summ.to_string())

    figs, pairs_md = [], []
    for kind, conf in res.confusion.items():
        conf.to_csv(out / f"confusion_{kind}.csv")
        figs.append(plots.confusion_heatmap(conf, f"{args.assay} {kind}: confusion (row-normalized)", out / f"confusion_{kind}.png"))
        cp = confusable_pairs(conf)
        cp.to_csv(out / f"confusable_pairs_{kind}.csv", index=False)
        if len(cp):
            pairs_md.append(f"**{kind}** most confused pairs:\n\n" + report.df_to_md(cp))
    best = summ.sort_values("balanced_accuracy_mean", ascending=False).iloc[0]
    body = [
        f"{om.notes[-1]}; label = `{args.label}` ({om.meta[args.label].nunique()} classes); "
        f"{n_splits}-fold animal-grouped CV, tuning inside folds, variance prefilter {prefilter}.",
        report.df_to_md(summ[["model", "balanced_accuracy_mean", "balanced_accuracy_std", "macro_f1_mean",
                              "log_loss_mean", "n_folds", "fit_seconds_mean"]].rename(columns=lambda c: c.replace("_mean", "")
                                                                  .replace("balanced_accuracy_std", "bal.acc sd")), floatfmt=".3f"),
        "Animals per outer fold (train / test) — the same folds for every model:\n\n" + report.df_to_md(folds),
        f"Best baseline: **{best['model']}** (balanced accuracy {best['balanced_accuracy_mean']:.3f} ± "
        f"{best['balanced_accuracy_std']:.3f}). Anything fancier has to beat this on the same folds.",
        report.figure_md(figs[0], f"{args.assay}: confusion matrix, {kinds[0]}") if figs else "",
        "\n\n".join(pairs_md) if pairs_md else "No confusions at all — with all features, tissue identity is trivial; the panel-size and shift phases are where the content is.",
    ]
    subset = f", {om.meta['tissue'].nunique()} tissues" if args.tissues else ""
    report.add_section(f"04 · Tissue-fingerprint baselines ({args.assay}/{source}{subset})", "\n\n".join(b for b in body if b),
                       params={"assay": args.assay, "source": source, "n_splits": n_splits, "quick": args.quick,
                               "tissues": args.tissues, "models": args.models})
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
