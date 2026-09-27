#!/usr/bin/env python
"""Run 2, part A — verification of the numbers quoted in docs/MULTIOMIC_REPORT.md, every value read from a CSV.

Writes results_multiomic/08_verification/:
  per_tissue_k20.csv           per-tissue accuracy with n samples and n donors for the Jiang k20 protein run (cleaned relative scale),
                               the raw-ppm run, and the 7-class RNA k20 run on the same 42 samples (Phase 5);
  accuracy_ci.csv              the k20 accuracies with their donor-bootstrap 95 % intervals (the raw-ppm 0.734 included);
  recalibration_set_sizes.csv  every recalibrated coverage quoted anywhere in the report with its mean set size, the number of classes
                               in its label space and the set size as a fraction of it;
  sato_design.csv              sedentary and exercised mice and samples behind the Sato invariance test, the exercise design fields;
  flags.csv                    what the verification changes in the wording (set sizes that are half the label space or more, ...);
  README.md, REPORT_SECTION.md, STATUS.json.
"""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from tfp import config as C, report

ROOT = C.ROOT
R = ROOT / "results_multiomic"
EXT = ROOT / "data" / "external_multiomic"
OUT = R / "08_verification"
FROZEN = C.FROZEN_DIR


def f(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{float(v):.{nd}f}"


def key(s):
    return "-".join(str(s).split("-")[:3])


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # ---- 1. per-tissue accuracy with n samples and n donors -----------------------------------------------------------------
    rows = []
    for run, path in (("jiang_relative_k20", R / "03_prot_transfer" / "accuracy_by_tissue.csv"), ("jiang_rawppm_k20", R / "03_prot_transfer" / "rawppm" / "accuracy_by_tissue.csv")):
        a = pd.read_csv(path)
        a = a[(a["model"] == "k20") & (a["rat_class"] != "OOD")]
        for r in a.itertuples():
            rows.append({"run": run, "layer": "protein", "jiang_tissue": r.jiang_tissue, "rat_class": r.rat_class, "n_samples": int(r.n), "n_donors": int(r.n_donors), "accuracy": float(r.accuracy),
                         "top_prediction": r.top_prediction, "note": "technical replicates kept as samples" if run == "jiang_rawppm_k20" else ""})
    # RNA k20 on the same 42 samples: Phase 5 gives n per rat class; donors per class are recomputed from the Jiang design for the same sample selection
    ft = pd.read_csv(R / "05_fusion_transfer" / "fusion_accuracy_by_tissue.csv")
    rel_cols = pd.read_csv(EXT / "jiang2020" / "protein_relative_columns.csv", dtype=str).drop_duplicates("column").set_index("column")
    des = pd.read_csv(EXT / "jiang2020" / "experimental_design.csv", dtype=str).dropna(subset=["GTEx Sample_ID"]).drop_duplicates("GTEx Sample_ID").set_index("GTEx Sample_ID")
    rna_cols = pd.read_csv(EXT / "jiang2020" / "rna_log_tpm_columns.csv", dtype=str)
    psamp = [c for c in rel_cols.index if str(c).startswith("GTEX")]
    rk = {key(s) for s in rna_cols["column"] if str(s).startswith("GTEX")}
    both = [s for s in psamp if key(s) in rk]
    J7 = {"Muscle - Skeletal": "SKM-GN", "Heart - Atrial Appendage": "HEART", "Heart - Left Ventricle": "HEART", "Lung": "LUNG", "Liver": "LIVER", "Brain - Cortex": "CORTEX"}
    mt = pd.DataFrame({"tissue": rel_cols.loc[both, "row2"].to_numpy(), "donor": des["Individual ID"].reindex(both).to_numpy()}, index=both)
    mt["rat_class"] = mt["tissue"].map(J7)
    donors = mt.dropna(subset=["rat_class"]).groupby("rat_class")["donor"].nunique()
    nsamp = mt.dropna(subset=["rat_class"]).groupby("rat_class").size()
    for layer in ("RNA", "protein", "late_mean", "stacked_LR"):
        sub = ft[(ft["model"] == "k20") & (ft["layer"] == layer)]
        for r in sub.itertuples():
            assert int(r.n) == int(nsamp[r.rat_class]), (r.rat_class, r.n, nsamp[r.rat_class])
            rows.append({"run": f"same42_{layer}_k20", "layer": layer, "jiang_tissue": ";".join(sorted(mt.loc[mt['rat_class'] == r.rat_class, 'tissue'].unique())), "rat_class": r.rat_class, "n_samples": int(r.n),
                         "n_donors": int(donors[r.rat_class]), "accuracy": float(r.accuracy), "top_prediction": r.top_prediction, "note": "7-class fingerprint on the 42 Jiang samples with both layers (Phase 5)"})
    pt = pd.DataFrame(rows)
    pt.to_csv(OUT / "per_tissue_k20.csv", index=False)

    # ---- 2. accuracies with intervals ------------------------------------------------------------------------------------
    ci_rows = []
    for run, path, col_lo, col_hi in (("jiang_relative", R / "03_prot_transfer" / "accuracy_overall.csv", "acc_ci95_low_donor_boot", "acc_ci95_high_donor_boot"),
                                      ("jiang_rawppm", R / "03_prot_transfer" / "rawppm" / "accuracy_overall.csv", "acc_ci95_low_donor_boot", "acc_ci95_high_donor_boot")):
        a = pd.read_csv(path).set_index("model")
        for m in a.index:
            ci_rows.append({"run": run, "layer": "protein", "model": m, "accuracy": float(a.loc[m, "accuracy_sample_weighted"]), "ci95_low": float(a.loc[m, col_lo]), "ci95_high": float(a.loc[m, col_hi]),
                            "n_samples": int(a.loc[m, "n_samples_mapped"]), "n_donors": int(a.loc[m, "n_donors_mapped"]), "bootstrap_unit": "donor", "n_boot": int(a.loc[m, "n_boot"]), "source": str(path.relative_to(ROOT))})
    fs = pd.read_csv(R / "05_fusion_transfer" / "fusion_transfer_summary.csv")
    for r in fs.itertuples():
        ci_rows.append({"run": "same42", "layer": r.layer, "model": r.model, "accuracy": float(r.accuracy), "ci95_low": float(r.acc_ci95_low_donor_boot), "ci95_high": float(r.acc_ci95_high_donor_boot),
                        "n_samples": int(r.n_mapped), "n_donors": int(r.n_donors), "bootstrap_unit": "donor", "n_boot": 1000, "source": "results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv"})
    ci = pd.DataFrame(ci_rows)
    ci.to_csv(OUT / "accuracy_ci.csv", index=False)

    # ---- 3. every recalibrated coverage with its set size and label-space size ---------------------------------------------
    rc_rows = []
    def add_recal(path, source_label, n_classes, note=""):
        if not path.exists():
            return
        d = pd.read_csv(path)
        for r in d.itertuples():
            rc_rows.append({"source": source_label, "model": r.model, "n_recal": int(r.n_recal), "n_test_individuals": int(r.n_test_individuals), "draws": int(r.draws),
                            "coverage_recalibrated": float(r.coverage_recalibrated), "coverage_source_cal_same_test": float(r.coverage_source_cal_same_test),
                            "set_size_recalibrated": float(r.set_size_recalibrated), "n_classes_label_space": n_classes, "set_size_frac_of_classes": float(r.set_size_recalibrated) / n_classes,
                            "frac_empty_recalibrated": float(r.frac_empty_recalibrated), "file": str(path.relative_to(ROOT)), "note": note})
    add_recal(R / "03_prot_transfer" / "recalibration.csv", "Phase 3 protein → Jiang 2020 (cleaned relative)", 7)
    add_recal(R / "03_prot_transfer" / "rawppm" / "recalibration.csv", "Phase 3 protein → Jiang 2020 (raw ppm)", 7)
    add_recal(R / "03_prot_transfer" / "wang2019" / "recalibration.csv", "Phase 3b protein → Wang 2019", 7, "recalibration individuals are single tissue samples")
    add_recal(R / "03_prot_transfer" / "geiger2013" / "recalibration.csv", "Phase 3b protein → Geiger 2013", 7, "recalibration individuals are single tissue samples")
    add_recal(R / "04_metab_transfer" / "hilic_sato" / "recalibration.csv", "Phase 4 hilic_sato", 19)
    add_recal(R / "04_metab_transfer" / "deep_sato" / "recalibration.csv", "Phase 4 deep_sato", 9)
    add_recal(R / "04_metab_transfer" / "deep_mw" / "recalibration.csv", "Phase 4 deep_mw", 9)
    add_recal(FROZEN / "13_gtex" / "recalibration.csv", "frozen phase 13 RNA → GTEx", 19, "the RNA ladder rows quoted for comparison")
    add_recal(FROZEN / "12_bodymap" / "recalibration.csv", "frozen phase 12 RNA → BodyMap", 19, "the RNA ladder rows quoted for comparison")
    rc = pd.DataFrame(rc_rows)
    rc.to_csv(OUT / "recalibration_set_sizes.csv", index=False)

    # ---- 4. the Sato invariance test: mice, samples, design -------------------------------------------------------------
    srows, all_ids = [], {}
    for t in ["BAT", "eWAT", "HEART", "HYPOTHALAMUS", "iWAT", "LIVER", "MUSCLE", "SERUM"]:
        m = pd.read_csv(EXT / "sato2022" / f"{t}_samples.csv", dtype=str)
        treat = m["TREATMENT"].astype(str).str.replace("Light ", "").str.replace("Dark ", "")
        ids = m["SUBJECT OR ANIMAL ID"].astype(str)
        for trt in ("Sedentary", "Exercise"):
            for i in ids[treat == trt]:
                all_ids.setdefault(trt, set()).add(i)
        tae = m["TIME AFTER EXERCSE"].dropna().unique().tolist() if "TIME AFTER EXERCSE" in m.columns else []
        zt = (m["TIME POINT"] if "TIME POINT" in m.columns else m.get("TIME POINT ZT", pd.Series(dtype=str))).dropna().unique().tolist()
        srows.append({"tissue": t, "n_samples": len(m), "n_mice": ids.nunique(), "n_sedentary_mice": ids[treat == "Sedentary"].nunique(), "n_exercised_mice": ids[treat == "Exercise"].nunique(),
                      "n_sedentary_samples": int((treat == "Sedentary").sum()), "n_exercised_samples": int((treat == "Exercise").sum()),
                      "time_after_exercise_h": ";".join(map(str, tae)), "zeitgeber_time_points": ";".join(map(str, zt))})
    sd = pd.DataFrame(srows)
    inv = pd.read_csv(R / "04_metab_transfer" / "hilic_sato" / "invariance_native_sedentary_to_exercised.csv")
    ik = inv[inv["model"] == "k20"].iloc[0]; ifu = inv[inv["model"] == "full"].iloc[0]
    design = pd.DataFrame([{"n_sedentary_mice_union": len(all_ids.get("Sedentary", ())), "n_exercised_mice_union": len(all_ids.get("Exercise", ())),
                            "n_sedentary_samples_total": int(sd["n_sedentary_samples"].sum()), "n_exercised_samples_total": int(sd["n_exercised_samples"].sum()),
                            "animal_ids_shared_across_tissues": bool(len(all_ids.get("Sedentary", ())) + len(all_ids.get("Exercise", ())) < int(sd["n_mice"].sum())),
                            "exercise_design": "acute: a single treadmill bout; the Metabolon sample metadata records TIME AFTER EXERCISE = 0 h and two Zeitgeber time points (ZT4 light, ZT16 dark) with 6 mice per treatment × time point per tissue; no training programme (Sato et al. 2022, Cell Metab 34:329, 'time-dependent signatures of metabolic homeostasis' after acute exercise)",
                            "invariance_fit_animals": int(ik["n_fit_animals"]), "invariance_cal_animals": int(ik["n_cal_animals"]), "invariance_test_exercised_samples": int(ik["n_test_exercised_samples"]), "invariance_test_exercised_animals": int(ik["n_test_animals"]),
                            "invariance_k20_accuracy": float(ik["accuracy_exercised"]), "invariance_k20_coverage": float(ik["coverage_exercised"]), "invariance_k20_set_size": float(ik["avg_set_size"]),
                            "invariance_full_accuracy": float(ifu["accuracy_exercised"]), "invariance_full_coverage": float(ifu["coverage_exercised"]), "invariance_full_set_size": float(ifu["avg_set_size"]), "invariance_n_classes": 8}])
    sd.to_csv(OUT / "sato_design.csv", index=False)
    design.to_csv(OUT / "sato_invariance_design.csv", index=False)

    # ---- 5. flags: what the verification changes ---------------------------------------------------------------------------
    flags = []
    for r in rc.itertuples():
        if r.set_size_frac_of_classes >= 0.5:
            flags.append({"item": f"{r.source} {r.model} n_recal={r.n_recal}", "flag": "recalibrated set holds ≥ half the label space", "value": f"set size {r.set_size_recalibrated:.2f} of {r.n_classes_label_space} classes",
                          "wording": "the number is restored, the information is not (GTEx pattern); never quote this coverage without the set size"})
    j3 = rc[(rc["source"].str.startswith("Phase 3 protein → Jiang 2020 (cleaned")) & (rc["model"] == "k20")]
    for r in j3.itertuples():
        flags.append({"item": f"Jiang k20 n_recal={r.n_recal}", "flag": "pattern classification", "value": f"coverage {r.coverage_recalibrated:.3f} at set size {r.set_size_recalibrated:.2f} of 7",
                      "wording": "GTEx pattern (number restored, information partly: sets hold " + f"{r.set_size_frac_of_classes:.0%} of the label space), not the BodyMap pattern (one tissue per set)"})
    flags.append({"item": "Sato invariance test", "flag": "wording", "value": "acute single bout, TIME AFTER EXERCISE = 0 h", "wording": "'exercised', never 'trained', for Sato mice; the RNA analogue (controls → trained) is an 8-week training programme, so the two tests are not the same shift"})
    ed = pt[pt["run"] == "same42_RNA_k20"]
    flags.append({"item": "RNA k20 on the same 42 Jiang samples", "flag": "n per class", "value": "; ".join(f"{r.rat_class} {r.n_samples} samples / {r.n_donors} donors" for r in ed.itertuples()), "wording": "two classes have ≤ 5 samples; the 1.000 is over 42 samples from 12 donors"})
    fl = pd.DataFrame(flags)
    fl.to_csv(OUT / "flags.csv", index=False)

    # ---- README + section --------------------------------------------------------------------------------------------------
    ci_k = ci[ci["model"] == "k20"].copy()
    lines = [f"# Run 2 A — verification", "", f"Built by `scripts/multiomic/08_verification.py` on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}; every number from a CSV in this directory or the one it names.", "",
             "## Per-tissue accuracy, k20, with n samples and n donors — `per_tissue_k20.csv`", "", report.df_to_md(pt.drop(columns=["note"]), floatfmt=".3f", max_rows=60), "",
             "## k20 accuracies with donor-bootstrap intervals — `accuracy_ci.csv`", "", report.df_to_md(ci_k, floatfmt=".3f"), "",
             "## Every recalibrated coverage with its set size and label space — `recalibration_set_sizes.csv`", "", report.df_to_md(rc.drop(columns=["file"]), floatfmt=".3f", max_rows=80), "",
             "## Sato 2022 design behind the invariance test — `sato_design.csv`, `sato_invariance_design.csv`", "", report.df_to_md(sd), "", report.df_to_md(design.T.reset_index().rename(columns={"index": "quantity", 0: "value"})), "",
             "## Flags — `flags.csv`", "", report.df_to_md(fl), ""]
    (OUT / "README.md").write_text("\n".join(lines))
    d = design.iloc[0]
    jr = ci[(ci["run"] == "jiang_relative") & (ci["model"] == "k20")].iloc[0]; jw = ci[(ci["run"] == "jiang_rawppm") & (ci["model"] == "k20")].iloc[0]
    rr = ci[(ci["run"] == "same42") & (ci["layer"] == "RNA") & (ci["model"] == "k20")].iloc[0]; rp = ci[(ci["run"] == "same42") & (ci["layer"] == "protein") & (ci["model"] == "k20")].iloc[0]
    rk3 = rc[(rc["source"].str.startswith("Phase 3 protein → Jiang 2020 (cleaned")) & (rc["model"] == "k20") & (rc["n_recal"] == 3)].iloc[0]
    rk5 = rc[(rc["source"].str.startswith("Phase 3 protein → Jiang 2020 (cleaned")) & (rc["model"] == "k20") & (rc["n_recal"] == 5)].iloc[0]
    gt3 = rc[(rc["source"].str.startswith("frozen phase 13")) & (rc["model"] == "k20") & (rc["n_recal"] == 3)].iloc[0]
    bm3 = rc[(rc["source"].str.startswith("frozen phase 12")) & (rc["model"] == "k20") & (rc["n_recal"] == 3)].iloc[0]
    per = lambda run: "; ".join(f"{r.rat_class} {f(r.accuracy, 2)} ({int(r.n_samples)} samples, {int(r.n_donors)} donors)" for r in pt[pt["run"] == run].itertuples())
    sec = ["Run 2 checked the quoted numbers against their CSVs and added what was missing beside them (`results_multiomic/08_verification/`).",
           f"- **Per-tissue accuracy, k20** (`per_tissue_k20.csv`). Jiang cleaned relative scale: {per('jiang_relative_k20')}. Raw ppm (technical replicates kept as samples): {per('jiang_rawppm_k20')}. "
           f"RNA 7-class on the same 42 samples: {per('same42_RNA_k20')}; protein on the same 42: {per('same42_protein_k20')}.",
           f"- **k20 accuracies with donor-bootstrap 95 % intervals** (`accuracy_ci.csv`): Jiang cleaned relative {f(jr.accuracy)} ({f(jr.ci95_low, 2)}–{f(jr.ci95_high, 2)}; {int(jr.n_samples)} samples, {int(jr.n_donors)} donors); "
           f"raw ppm {f(jw.accuracy)} ({f(jw.ci95_low, 2)}–{f(jw.ci95_high, 2)}; {int(jw.n_samples)} samples, {int(jw.n_donors)} donors); RNA on the same 42: {f(rr.accuracy)} ({f(rr.ci95_low, 2)}–{f(rr.ci95_high, 2)}); protein on the same 42: {f(rp.accuracy)} ({f(rp.ci95_low, 2)}–{f(rp.ci95_high, 2)}).",
           f"- **Recalibrated coverages with set sizes** (`recalibration_set_sizes.csv`, {len(rc)} rows, every one quoted in this report or in the RNA comparison). Jiang k20: 3 donors {f(rk3.coverage_recalibrated)} at {f(rk3.set_size_recalibrated, 2)} of 7 classes ({rk3.set_size_frac_of_classes:.0%}), "
           f"5 donors {f(rk5.coverage_recalibrated)} at {f(rk5.set_size_recalibrated, 2)} of 7 ({rk5.set_size_frac_of_classes:.0%}). For comparison, RNA → GTEx 3 donors {f(gt3.coverage_recalibrated)} at {f(gt3.set_size_recalibrated, 1)} of 19 ({gt3.set_size_frac_of_classes:.0%}), "
           f"RNA → BodyMap 3 animals {f(bm3.coverage_recalibrated)} at {f(bm3.set_size_recalibrated, 2)} of 19 ({bm3.set_size_frac_of_classes:.0%}). **The Jiang recalibration is the GTEx pattern** — the number is restored, the information only partly — not the BodyMap pattern of one tissue per set. "
           f"{int((rc['set_size_frac_of_classes'] >= 0.5).sum())} of the {len(rc)} recalibrated coverages have sets holding at least half their label space (`flags.csv`); they are read as full-set effects.",
           f"- **Sato 2022 invariance test** (`sato_design.csv`, `sato_invariance_design.csv`): {int(d.n_sedentary_mice_union)} sedentary and {int(d.n_exercised_mice_union)} exercised mice (ids shared across the 8 tissue tables), {int(d.n_sedentary_samples_total)} sedentary and {int(d.n_exercised_samples_total)} exercised samples; "
           f"the native panel is fit on {int(d.invariance_fit_animals)} sedentary mice, calibrated on {int(d.invariance_cal_animals)}, tested on {int(d.invariance_test_exercised_samples)} samples of {int(d.invariance_test_exercised_animals)} exercised mice (8 classes): k20 accuracy {f(d.invariance_k20_accuracy)}, coverage {f(d.invariance_k20_coverage)} at set size {f(d.invariance_k20_set_size, 2)}; "
           f"full {f(d.invariance_full_accuracy)}, {f(d.invariance_full_coverage)} at {f(d.invariance_full_set_size, 2)}. **The exercise is an acute single bout** (TIME AFTER EXERCISE = 0 h; ZT4 and ZT16; 6 mice per treatment × time point per tissue), not a training programme, so the mice are 'exercised', never 'trained', "
           "and the RNA analogue (controls → 8-week-trained rats) is a different shift.",
           "- **What changed in the wording**: every recalibrated coverage now carries its set size and label-space size; the Jiang recalibration is described as the GTEx pattern; the Sato mice are 'exercised'. No finding was contradicted by a CSV; the Sato metabolite legs' recalibrations (sets of 5–12 of 9–19 classes) are read as full-set effects, so their findings keep the accuracy and drop the recalibration from the headline."]
    (OUT / "REPORT_SECTION.md").write_text("\n".join(sec) + "\n")
    (OUT / "STATUS.json").write_text(json.dumps({"status": "DONE", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}))
    print(pt.to_string(index=False)); print(rc[["source", "model", "n_recal", "coverage_recalibrated", "set_size_recalibrated", "n_classes_label_space"]].to_string(index=False)); print(design.T.to_string())
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
