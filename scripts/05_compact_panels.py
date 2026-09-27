#!/usr/bin/env python
"""Phase 05 — how small can a tissue fingerprint be, and is the selection stable?

Outputs (results/05_panels/<assay>/): panel_curve.csv, panel_curve.png, selected_by_fold.csv,
stability_k<k>.csv, candidate_panel.csv (features chosen in >= 80% of bootstraps, annotated)
"""
from __future__ import annotations

import pandas as pd

from tfp import cli, config as C, io, models, plots, report


def main() -> None:
    ap = cli.common_parser("Compact panel curves and stability selection")
    ap.add_argument("--label", default="tissue")
    ap.add_argument("--model", default="logreg_l2", help="classifier used on top of the k selected features")
    ap.add_argument("--grid", default=",".join(map(str, C.PANEL_GRID)))
    ap.add_argument("--selector", default="roundrobin", choices=["roundrobin", "fclassif"],
                    help="panel selector; fclassif is kept to demonstrate why univariate F-scores fail for multiclass panels")
    ap.add_argument("--compare-selectors", action="store_true", help="also run the F-test selector for comparison")
    ap.add_argument("--stability-k", type=int, default=10)
    ap.add_argument("--n-boot", type=int, default=50)
    ap.add_argument("--confusion-k", default="20,30", help="panel sizes at which to report confused tissue pairs")
    args = ap.parse_args()
    cli.banner("05_compact_panels", args)
    source = cli.resolve_source(args.assay, args.source)
    out = cli.outdir(f"05_panels/{args.assay}", args.out)
    pheno = io.load_pheno()
    om = io.stack_tissues(args.assay, tissues=cli.parse_tissues(args.tissues), source=source, join="inner", pheno=pheno,
                           complete=args.complete_features,
                          drop_incomplete_samples=args.drop_incomplete_samples)
    grid = [int(k) for k in args.grid.split(",")]
    n_splits = 3 if args.quick else args.n_splits
    n_boot = 10 if args.quick else args.n_boot
    prefilter = min(args.prefilter, 1000) if args.quick else args.prefilter

    curve, selected, preds = models.panel_curve(om, label=args.label, grid=grid, kind=args.model, n_splits=n_splits,
                                                prefilter=prefilter, quick=args.quick, seed=args.seed,
                                                selector=args.selector, return_predictions=True)
    curve.to_csv(out / "panel_curve.csv", index=False)
    preds.to_csv(out / "panel_predictions.csv", index=False)
    classes = sorted(om.meta[args.label].astype(str).unique())
    pairs_md = []
    for kk in [int(v) for v in args.confusion_k.split(",") if v.strip()]:
        pk = preds[preds["k"] == kk]
        if pk.empty:
            continue
        conf = pd.crosstab(pk["y_true"], pk["y_pred"]).reindex(index=classes, columns=classes, fill_value=0)
        conf.to_csv(out / f"confusion_k{kk}.csv")
        rows = [{"k": kk, "true": a, "predicted": b, "count": int(conf.loc[a, b]),
                 "frac_of_true": float(conf.loc[a, b] / max(conf.loc[a].sum(), 1))}
                for a in classes for b in classes if a != b and conf.loc[a, b] > 0]
        cp = pd.DataFrame(rows, columns=["k", "true", "predicted", "count", "frac_of_true"]).sort_values("count", ascending=False)
        cp.to_csv(out / f"confusable_pairs_k{kk}.csv", index=False)
        pairs_md.append(f"Confused pairs at k={kk} (pooled over folds; {int(cp['count'].sum())} errors of {len(pk)} test samples):\n\n"
                        + (report.df_to_md(cp.head(12), floatfmt=".2f") if len(cp) else "_none_"))
    compare_md = ""
    if args.compare_selectors:
        other = "fclassif" if args.selector == "roundrobin" else "roundrobin"
        curve2, _ = models.panel_curve(om, label=args.label, grid=grid, kind=args.model, n_splits=n_splits,
                                       prefilter=prefilter, quick=args.quick, seed=args.seed, selector=other, verbose=False)
        curve2.to_csv(out / f"panel_curve_{other}.csv", index=False)
        cmp = pd.concat([curve.assign(selector=args.selector), curve2.assign(selector=other)])
        cmp_agg = cmp.groupby(["selector", "k"])["balanced_accuracy"].mean().unstack("selector").reset_index()
        compare_md = "Selector comparison (mean balanced accuracy):\n\n" + report.df_to_md(cmp_agg)
    selected.to_csv(out / "selected_by_fold.csv", index=False)
    agg = curve.groupby("k").agg(bal_acc_mean=("balanced_accuracy", "mean"), bal_acc_sd=("balanced_accuracy", "std"),
                                 macro_f1_mean=("macro_f1", "mean")).reset_index()
    fig = plots.panel_curve_plot(curve, "balanced_accuracy", f"{args.assay}: accuracy vs panel size ({args.model})",
                                 out / "panel_curve.png")

    # how consistent is the selection across folds at each k?
    cons = (selected.groupby(["k", "feature_ID"]).size().reset_index(name="n_folds")
            .assign(frac_folds=lambda d: d["n_folds"] / n_splits))
    cons_summary = cons.groupby("k").apply(lambda d: pd.Series({"n_unique_features": len(d),
                                                                "frac_in_all_folds": float((d["frac_folds"] == 1).mean())}),
                                           include_groups=False).reset_index()
    cons_summary.to_csv(out / "selection_consistency.csv", index=False)

    stab = models.stability_selection(om, label=args.label, k=args.stability_k, n_boot=n_boot,
                                      prefilter=prefilter, seed=args.seed, selector=args.selector)
    stab.to_csv(out / f"stability_k{args.stability_k}.csv", index=False)
    try:
        stab["gene_symbol"] = io.map_to_gene_symbols(stab["feature_ID"]).values
    except Exception as e:  # annotation is optional
        stab["gene_symbol"] = None
        print(f"  (gene annotation skipped: {e})")
    panel = stab[stab["selection_frequency"] >= 0.8]
    panel.to_csv(out / "candidate_panel.csv", index=False)

    smallest_perfect = agg.loc[agg["bal_acc_mean"] >= 0.99, "k"].min() if (agg["bal_acc_mean"] >= 0.99).any() else None
    body = [
        f"{om.notes[-1]}; classifier on selected features: `{args.model}`; {n_splits} animal-grouped folds; "
        f"selection inside each fold (`{args.selector}`), variance prefilter {prefilter}.",
        report.df_to_md(agg, floatfmt=".3f"),
        report.figure_md(fig, "Balanced accuracy vs panel size (mean ± sd over folds)"),
        f"Smallest k with mean balanced accuracy ≥ 0.99: **{smallest_perfect}** (point estimate — see phase 06 for the certified k).",
        compare_md,
        "Selection consistency across folds:\n\n" + report.df_to_md(cons_summary, floatfmt=".2f"),
        "\n\n".join(pairs_md),
        f"Stability selection (k={args.stability_k}, {n_boot} animal-bootstraps): features selected in ≥ 80% of resamples:\n\n"
        + (report.df_to_md(panel, floatfmt=".2f") if len(panel) else "_none — the panel is not stable at this k_"),
    ]
    report.add_section(f"05 · Compact panels ({args.assay}/{source})", "\n\n".join(b for b in body if b),
                       params={"assay": args.assay, "model": args.model, "grid": args.grid, "quick": args.quick})
    print(agg.to_string())
    print(f"wrote {out}")


if __name__ == "__main__":
    main()
