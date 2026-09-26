#!/usr/bin/env python
"""Phase 06/08 read-out: certificate distribution, coverage and the shift table pulled from
results/06_conformal/* and results/08_shift/*, with the interpretation for the hackathon.
Writes results/08_shift/readout_*.csv and appends one REPORT.md section.
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
    ap = cli.common_parser("Conformal + shift read-out")
    args = ap.parse_args()
    cli.banner("08_readout", args)
    out = cli.outdir("08_shift", args.out)
    body = []

    # ---- certificates --------------------------------------------------------------------
    certs = []
    for label, d in (("primary: one vial per animal, α=δ=0.10, 22 cal animals", "TRNSCRPT"),
                     ("secondary: pooled vials, α=δ=0.05, 30% cal (12 animals ≈ 216 vials)", "TRNSCRPT_pooled_a05")):
        v = read(R / "06_conformal" / d / "certificate_validity.csv")
        dist = read(R / "06_conformal" / d / "certificate_distribution.csv")
        fin = read(R / "06_conformal" / d / "certificate_final_repeats.csv")
        if v is None:
            continue
        ks = v["certified_k"].dropna()
        k_half = dist.loc[dist["frac_certified"] >= 0.5, "k"].min() if (dist["frac_certified"] >= 0.5).any() else None
        certs.append({"certificate": label, "splits": len(v), "certified": int(ks.size), "none": int(v["certified_k"].isna().sum()),
                      "k_min": int(ks.min()) if ks.size else None, "k_median": ks.median() if ks.size else None,
                      "k_max": int(ks.max()) if ks.size else None, "k_certified_in_half": k_half,
                      "validity_test_err_le_alpha": float(v["test_err_le_alpha"].dropna().astype(bool).mean()),
                      "deployment_k": ", ".join(f"{int(k)}×{n}" for k, n in fin["certified_k"].value_counts(dropna=False).sort_index().items() if not pd.isna(k))
                      + (f", none×{int(fin['certified_k'].isna().sum())}" if fin["certified_k"].isna().any() else ""),
                      "animal_any_wrong_certified": int(v["certified_k_animal_any_wrong"].notna().sum())})
    certs = pd.DataFrame(certs)
    certs.to_csv(out / "readout_certificates.csv", index=False)
    dist_p = read(R / "06_conformal" / "TRNSCRPT" / "certificate_distribution.csv")
    dist_s = read(R / "06_conformal" / "TRNSCRPT_pooled_a05" / "certificate_distribution.csv")
    dist = dist_p[["k", "frac_certified"]].rename(columns={"frac_certified": "primary_frac_certified"}).merge(
        dist_s[["k", "frac_certified"]].rename(columns={"frac_certified": "secondary_frac_certified"}), on="k", how="outer").sort_values("k")
    sizing = read(R / "06_conformal" / "TRNSCRPT" / "sizing_table.csv")
    cov = read(R / "06_conformal" / "TRNSCRPT" / "coverage.csv")
    cs = cov.groupby(["calibration", "method", "alpha"]).agg(coverage=("coverage", "mean"), sd=("coverage", "std"),
                                                             set_size=("avg_set_size", "mean"), frac_empty=("frac_empty", "mean"),
                                                             n_cal=("n_cal", "mean")).reset_index()
    pc = read(R / "06_conformal" / "TRNSCRPT" / "per_class_coverage.csv")
    p10 = pc[(pc["method"] == "lac") & (pc["alpha"] == 0.1)].groupby(["y_true", "calibration"])["coverage"].mean().unstack().reset_index()
    under = p10[p10[["pooled", "one_per_animal"]].min(axis=1) < 0.85].sort_values("pooled")
    body.append("### Certificates (T6)\n\n" + report.df_to_md(certs) + "\n\nFraction of fold × repeat splits certifying each k:\n\n"
                + report.df_to_md(dist, floatfmt=".2f") + "\n\n"
                "Read: with one vial per animal and 22 calibration animals the certificate at (0.10, 0.10) is a zero-error "
                "certificate — one wrong calibration animal at a panel size ends it there — so the certified k is a draw "
                "from a distribution, not a number: it lands at 15 or 20 in about half the splits and at nothing at all in 30% "
                "of them, while the test error at every certified k stayed below α (validity 1.00). The per-animal "
                "'any tissue wrong' loss never certifies: with ~18 vials per animal, roughly a third of animals have at "
                "least one wrong vial at any k. The pooled-vial certificate at (0.05, 0.05) uses ~216 vials as if they were "
                "independent (they are not: they come from 12 animals), lands at 15–20 in most splits and is the optimistic "
                "bracket. Zero-error sizing:\n\n" + report.df_to_md(sizing))
    body.append("### Coverage in distribution (T6, TRNSCRPT, held-out test animals)\n\n" + report.df_to_md(cs, floatfmt=".3f")
                + "\n\nPer-tissue LAC coverage at α = 0.1 below 0.85 under either calibration:\n\n" + report.df_to_md(under, floatfmt=".2f")
                + "\n\nRead: pooled LAC coverage sits on 1 − α (table above); one-vial-per-animal calibration with 22 points "
                "moves in steps of 1/23, and at α = 0.05 its threshold is the largest calibration score. APS never produces "
                "an empty set and over-covers. The under-covered tissues are the confusable ones (table above): the two "
                "muscles, the brain regions and the vena cava — the marginal guarantee holds while the muscles sit far below it.")

    # ---- shift ---------------------------------------------------------------------------
    for label, d in (("TRNSCRPT, k = 20, pooled source calibration", "TRNSCRPT"),
                     ("METAB core 9 (193 complete features), k = 10, pooled source calibration", "METAB_core9")):
        t = read(R / "08_shift" / d / "shift_table.csv")
        if t is None:
            continue
        cols = ["split", "arm", "n_test_animals", "n_cal", "unseen_classes", "accuracy_all", "bal_acc_seen", "coverage_source_id",
                "coverage_target_seen", "coverage_drop", "set_size_change", "cov_target_recal_N3", "cov_target_recal_N5"]
        cols = [c for c in cols if c in t.columns]
        pcs = read(R / "08_shift" / d / "shift_per_class.csv")
        low = pcs[(pcs["coverage"] < 0.7) & (pcs["n"] >= 10)].sort_values("coverage")
        ctt = pcs[(pcs["split"] == "train_control_test_trained") & (pcs["y_true"].isin(["SMLINT", "BAT", "SKM-GN", "SKM-VL", "VENACV"]))]
        ctt_w = ctt.pivot(index="y_true", columns="arm", values="coverage").reset_index() if len(ctt) else None
        body.append(f"### Shift tests (T7) — {label}\n\n" + report.df_to_md(t[cols], floatfmt=".3f")
                    + "\n\nPer-class coverage below 0.70 (n ≥ 10): "
                    + ("; ".join(f"{r['split']}/{r['arm']}: {r['y_true']} {r['coverage']:.2f}" for _, r in low.head(12).iterrows()) or "none")
                    + (("\n\nControls→trained, tissues with training-regulated panel markers and the confusable pair:\n\n" + report.df_to_md(ctt_w, floatfmt=".2f")) if ctt_w is not None else ""))
    # 0b: the confusable tissues under held-out time points next to their in-distribution coverage
    pcs_t = read(R / "08_shift" / "TRNSCRPT" / "shift_per_class.csv")
    pc06 = read(R / "06_conformal" / "TRNSCRPT" / "per_class_coverage.csv")
    if pcs_t is not None and pc06 is not None:
        watch = ["SKM-GN", "SKM-VL", "VENACV", "HIPPOC", "HYPOTH"]
        indist = pc06[(pc06["calibration"] == "pooled") & (pc06["method"] == "lac") & (pc06["alpha"] == 0.1)
                      & (pc06.get("conformal", "marginal") == "marginal")].groupby("y_true")["coverage"].mean()
        rows = []
        for t in watch:
            row = {"tissue": t, "in_distribution_full": float(indist.get(t, np.nan))}
            for split in ("holdout_group_8w", "holdout_group_1w", "train_control_test_trained"):
                for arm in ("full", "panel_k20"):
                    v = pcs_t[(pcs_t["split"] == split) & (pcs_t["arm"] == arm) & (pcs_t["y_true"] == t)]["coverage"]
                    row[f"{split.replace('holdout_group_', '').replace('train_control_test_trained', 'ctrl→trained')}_{arm}"] = float(v.iloc[0]) if len(v) else np.nan
            rows.append(row)
        w = pd.DataFrame(rows)
        w.to_csv(out / "readout_confusable_under_shift.csv", index=False)
        body.append("Confusable tissues: in-distribution coverage (T6, pooled LAC α = 0.1) next to held-out 8w, held-out 1w and "
                    "controls→trained:\n\n" + report.df_to_md(w, floatfmt=".2f")
                    + "\n\nRead: if a tissue is under-covered under the 8w/1w holdouts as well (same 40/10 geometry as a CV fold), "
                    "the controls→trained collapse is the confusable pair, not training; if it is covered there and collapses only "
                    "under controls→trained, training biology (or the 10-animal fit set) is the cause.")
    t_all = read(R / "08_shift" / "TRNSCRPT" / "shift_table.csv")
    if t_all is not None and "lac_frac_empty_target" in t_all.columns:
        e = t_all[["split", "arm", "coverage_target_seen", "lac_frac_empty_source", "lac_frac_empty_target",
                   "aps_coverage_target_seen", "aps_set_size_source", "aps_set_size_target"]]
        body.append("Empty-set rate (LAC) and APS under the same shifts (TRNSCRPT, k = 20):\n\n" + report.df_to_md(e, floatfmt=".3f")
                    + "\n\nRead: LAC loses coverage under shift by returning empty sets, APS never does and pays in set size.")
    body.append("Read (tables above; numbers deliberately not repeated here, since they moved with the 2026-09-25 quantile "
                "fix, see results/QUANTILE_FIX_CHANGES.md): the guarantee breaks before the classifier does. Transcripts at "
                "k = 20: accuracy on seen classes holds under every shift, but coverage on the target falls below 1 − α when a "
                "model trained on one sex meets the other, and LAC sets go empty, because the source-calibrated threshold "
                "assumes the source's confidence. Held-out 8-week animals behave like an in-distribution fold (the same 40/10 "
                "geometry); held-out 1-week and controls→trained lose a little more. What breaks first is the confusable "
                "classes (the muscle pair, vena cava). Metabolomics is the sharper example: accuracy stays perfect under "
                "male→female while coverage falls far below 1 − α. Recalibrating on 3 target animals (pooled vials) brings "
                "every row back to about 1 − α, and 5 animals do no better. An animal-level certificate under the sex shift "
                "would need 22 zero-error calibration animals (the zero-error sizing at α = δ = 0.1).")

    report.add_section("06/08 · Certificate and shift read-out (real data)", "\n\n".join(body),
                       params={"inputs": "results/06_conformal/*, results/08_shift/*"})
    print("\n".join(b[:400] for b in body))
    print(f"wrote {out}/readout_*.csv and the REPORT.md section")


if __name__ == "__main__":
    main()
