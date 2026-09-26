#!/usr/bin/env python
"""Phase 05b — annotate a candidate panel. Annotation only: nothing here filters or re-ranks it.

Per panel feature: gene symbol; the tissue whose one-vs-rest score is highest (the class the
round-robin selector takes it for), that score, its rank inside that tissue's list and the runner-up
tissue; effect size = mean in the marker tissue − highest mean among the other tissues (log2 CPM for
counts, normalized units otherwise); whether the feature is training-regulated in the marker tissue
or anywhere (`training_regulated_features.csv`); its within-marker-tissue Pearson correlation with a
library QC metric (`--qc-col`, default pct_mrna from data/raw/meta/<ASSAY>.csv, when present) and the
largest |r| over all tissues; and two risk flags: training-regulated in the marker tissue (a risk for
the controls→trained shift, phase 08) and |r| ≥ --qc-risk with the QC metric (a risk for external
validation under a different library chemistry).

The one-vs-rest scores are computed on the full stacked matrix for annotation only (no evaluation
uses them). Inputs: results/05_panels/<assay>/candidate_panel.csv (or --panel-file / --stability-file
with --min-freq). Outputs: candidate_panel_annotated.csv next to the input; one REPORT.md section.
"""
from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.impute import SimpleImputer

from motrpac import cli, config as C, io, models, report


def main() -> None:
    ap = cli.common_parser("Annotate a candidate panel")
    ap.add_argument("--label", default="tissue")
    ap.add_argument("--panel-file", default=None, help="default: results/05_panels/<assay>/candidate_panel.csv")
    ap.add_argument("--stability-file", default=None, help="annotate a stability_k*.csv instead (with --min-freq)")
    ap.add_argument("--min-freq", type=float, default=0.0, help="keep features with selection_frequency >= this")
    ap.add_argument("--qc-col", default="pct_mrna")
    ap.add_argument("--qc-risk", type=float, default=0.5, help="|r| with the QC metric that raises the flag")
    args = ap.parse_args()
    cli.banner("05_annotate_panel", args)
    source = cli.resolve_source(args.assay, args.source)
    out = cli.outdir(f"05_panels/{args.assay}", args.out)
    panel_path = Path(args.panel_file or args.stability_file or (out / "candidate_panel.csv"))
    panel = pd.read_csv(panel_path, dtype={"feature_ID": str})
    if args.min_freq > 0:
        panel = panel[panel["selection_frequency"] >= args.min_freq]
    panel = panel.drop_duplicates("feature_ID").reset_index(drop=True)
    pheno = io.load_pheno()
    om = io.stack_tissues(args.assay, tissues=cli.parse_tissues(args.tissues), source=source, join="inner", pheno=pheno,
                           complete=args.complete_features,
                          drop_incomplete_samples=args.drop_incomplete_samples)
    y = om.meta[args.label].astype(str).to_numpy()
    classes = sorted(np.unique(y))
    ids = [f for f in panel["feature_ID"] if f in om.X.columns]
    missing = [f for f in panel["feature_ID"] if f not in om.X.columns]
    if missing:
        print(f"  {len(missing)} panel features not in the stacked matrix: {missing[:5]}")
    if not ids:
        raise SystemExit("no panel feature found in the stacked matrix")

    # one-vs-rest scores on the full matrix (annotation only)
    full = SimpleImputer(strategy="median").fit_transform(om.X.to_numpy(dtype=float))
    sel = models.RoundRobinSelector(k=len(ids)).fit(full, y)
    scores = sel.scores_
    col = {f: j for j, f in enumerate(om.X.columns)}
    ranks = {i: np.argsort(-scores[i]) for i in range(len(classes))}
    rank_pos = {i: {j: r + 1 for r, j in enumerate(ranks[i])} for i in range(len(classes))}
    top_full = set(np.asarray(om.X.columns)[sel.idx_])

    mu = om.X[ids].groupby(y).mean()  # tissue × feature, NaN-skipping
    reg = pd.read_csv(C.RAW_DIR / "training_regulated_features.csv", dtype=str, usecols=["feature_ID", "assay", "tissue"])
    reg = reg[reg["assay"] == args.assay]
    reg_map = reg.groupby("feature_ID")["tissue"].agg(lambda s: sorted(set(s))).to_dict()
    qc = None
    meta_path = C.META_DIR / f"{args.assay}.csv"
    if meta_path.exists():
        m = pd.read_csv(meta_path, dtype=str, low_memory=False)
        if "viallabel" in m.columns and args.qc_col in m.columns:
            qc = pd.to_numeric(m.drop_duplicates("viallabel").set_index("viallabel")[args.qc_col], errors="coerce").reindex(om.X.index)
    try:
        symbols = io.map_to_gene_symbols(ids) if args.assay == "TRNSCRPT" else pd.Series(ids, index=ids)
    except Exception as e:  # annotation is optional
        print(f"  (gene symbols unavailable: {e})")
        symbols = pd.Series(ids, index=ids)

    rows = []
    for f in ids:
        j = col[f]
        sc = scores[:, j]
        order = np.argsort(-sc)
        i1, i2 = int(order[0]), int(order[1]) if len(order) > 1 else int(order[0])
        t1 = classes[i1]
        means = mu[f]
        others = means.drop(t1)
        x = om.X[f]
        r_in, r_max, t_max = np.nan, np.nan, None
        if qc is not None:
            for t in classes:
                mk = (y == t) & x.notna().to_numpy() & qc.notna().to_numpy()
                if mk.sum() >= 8 and x[mk].std() > 0 and qc[mk].std() > 0:
                    r = float(np.corrcoef(x[mk], qc[mk])[0, 1])
                    if t == t1:
                        r_in = r
                    if np.isnan(r_max) or abs(r) > abs(r_max):
                        r_max, t_max = r, t
        reg_t = reg_map.get(f, [])
        rows.append({
            "feature_ID": f, "gene_symbol": symbols.get(f) if hasattr(symbols, "get") else None,
            "selection_frequency": panel.set_index("feature_ID").loc[f, "selection_frequency"] if "selection_frequency" in panel.columns else np.nan,
            "marker_tissue": t1, "ovr_score": float(sc[i1]), "rank_in_tissue": int(rank_pos[i1][j]),
            "runner_up_tissue": classes[i2], "runner_up_score": float(sc[i2]),
            "in_full_data_roundrobin_top_k": f in top_full,
            "mean_in_marker_tissue": float(means[t1]), "next_highest_tissue": str(others.idxmax()),
            "next_highest_mean": float(others.max()), "effect_size": float(means[t1] - others.max()),
            "regulated_in_marker_tissue": t1 in reg_t, "regulated_tissues": ";".join(reg_t),
            f"r_{args.qc_col}_in_marker_tissue": r_in, f"max_abs_r_{args.qc_col}": r_max, "tissue_of_max_r": t_max,
            "risk_T7_regulated": t1 in reg_t,
            "risk_qc_correlated": bool(not np.isnan(r_in) and abs(r_in) >= args.qc_risk),
        })
    ann = pd.DataFrame(rows).sort_values(["marker_tissue", "ovr_score"], ascending=[True, False])
    ann.to_csv(out / (panel_path.stem + "_annotated.csv"), index=False)
    uncovered = [c for c in classes if c not in set(ann["marker_tissue"])]
    n_reg_any = int((ann["regulated_tissues"] != "").sum())
    body = [
        f"{om.notes[0]}; {len(ann)} panel features from `{panel_path.name}`"
        + (f" ({len(missing)} not in the stacked matrix)" if missing else "") + ". One-vs-rest scores and means are "
        "computed on the full stacked matrix for annotation only; nothing here changes the panel.",
        report.df_to_md(ann[["feature_ID", "gene_symbol", "selection_frequency", "marker_tissue", "ovr_score", "rank_in_tissue",
                             "runner_up_tissue", "effect_size", "next_highest_tissue", "regulated_in_marker_tissue",
                             "regulated_tissues", f"r_{args.qc_col}_in_marker_tissue", f"max_abs_r_{args.qc_col}",
                             "tissue_of_max_r", "risk_T7_regulated", "risk_qc_correlated"]], floatfmt=".2f"),
        f"Tissues with no panel marker of their own: **{', '.join(uncovered) or 'none'}**. "
        f"Training-regulated in the marker tissue: **{int(ann['risk_T7_regulated'].sum())}** of {len(ann)} "
        f"(regulated in any tissue: {n_reg_any}) — a risk for the controls→trained shift (phase 08) because the "
        f"marker's level moves with training. |r| ≥ {args.qc_risk} with `{args.qc_col}` in the marker tissue: "
        f"**{int(ann['risk_qc_correlated'].sum())}** — a risk for external validation under a different library "
        "chemistry (docs/EXTERNAL_VALIDATION.md). Flags are annotation, not exclusion.",
    ]
    report.add_section(f"05 · Candidate panel annotation ({args.assay}/{source})", "\n\n".join(body),
                       params={"assay": args.assay, "panel": str(panel_path), "qc_col": args.qc_col})
    print(ann.to_string(index=False, max_colwidth=28))
    print(f"wrote {out / (panel_path.stem + '_annotated.csv')}")


if __name__ == "__main__":
    main()
