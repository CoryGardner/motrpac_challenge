#!/usr/bin/env python
"""Multiomic Phase 6 — identifiability of the external designs (metadata only).

The audit (phase 16) showed that in MoTrPAC every processing batch holds one tissue (TMT plex, RNA extraction plate, library
batch …), so the tissue axis is not separable from processing inside the study. Here the same framing (`tfp.batch.nesting_table`,
`tfp.batch.estimable_pairs`) is applied to the sample metadata of every external dataset obtained in Phase 2:
  Jiang 2020      TMT run (56 plexes of 10 channels) × TMT tag × donor × tissue (Table S1 A / mmc2);
  Wang 2019       MS experiment id per tissue (Table EV1 A);
  Geiger 2013     one run per tissue, no per-sample batch variable shipped;
  Sato 2022       Metabolon ROUND / RUN DAY per sample (per-tissue tables);
  MW ST003188     Batch per sample (mwtab).
For each batch variable: levels per tissue, tissues per level, Cramér's V with tissue, and the number of tissue pairs sharing a
level (estimable pairs need a shared level of every variable listed). MoTrPAC's own numbers are read from
results_frozen/16_identifiability/ for the comparison table.
Outputs: results_multiomic/06_external_identifiability/ (nesting_<dataset>.csv, estimable_pairs.csv, design_comparison.csv,
README.md, REPORT_SECTION.md, FINDINGS.json, STATUS.json).
"""
from __future__ import annotations

import json
import time

import numpy as np
import pandas as pd

from tfp import batch, config as C, report

ROOT = C.ROOT
EXT = ROOT / "data" / "external_multiomic"
OUT = ROOT / "results_multiomic" / "06_external_identifiability"
FROZEN = C.FROZEN_DIR / "16_identifiability"


def f(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{float(v):.{nd}f}"


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    designs = {}
    # Jiang 2020: study samples (reference channels excluded), tissue = Tissue Name
    des = pd.read_csv(EXT / "jiang2020" / "experimental_design.csv", dtype=str)
    j = des[des["Tissue Name"].astype(str).str.lower() != "reference"].dropna(subset=["Tissue Name"]).copy()
    j = j.rename(columns={"Run": "tmt_run", "TMT Tag": "tmt_tag", "Individual ID": "donor", "Tissue Name": "tissue"})
    designs["Jiang2020"] = (j, ["tmt_run", "tmt_tag", "donor"], ["tmt_run"], "TMT run = plex (10 channels + 1 reference); the batch variable of a TMT design")
    # Wang 2019: one MS experiment per tissue sample
    w = pd.read_csv(EXT / "wang2019" / "sample_information.csv", dtype=str).rename(columns={"Tissue name": "tissue", "MS Experiment ID": "ms_experiment", "Gender of donor": "sex"})
    w = w.dropna(subset=["tissue", "ms_experiment"])
    designs["Wang2019"] = (w, ["ms_experiment"], ["ms_experiment"], "one MS experiment (label-free run) per tissue sample")
    # Sato 2022: ROUND / RUN DAY per sample, tissues in separate Metabolon tables
    srows = []
    for t in ["BAT", "eWAT", "HEART", "HYPOTHALAMUS", "iWAT", "LIVER", "MUSCLE", "SERUM"]:
        m = pd.read_csv(EXT / "sato2022" / f"{t}_samples.csv", dtype=str)
        srows.append(pd.DataFrame({"tissue": t, "round": m.get("ROUND"), "run_day": m.get("RUN DAY"), "animal": m.get("SUBJECT OR ANIMAL ID"), "treatment": m.get("TREATMENT")}))
    s = pd.concat(srows, ignore_index=True)
    s["round_tissue_table"] = s["tissue"] + ":" + s["round"].astype(str)   # a Metabolon round is reported per tissue table
    designs["Sato2022"] = (s, ["round", "run_day", "animal"], ["round"], "Metabolon ROUND (per tissue table) and RUN DAY where given; ROUND labels repeat across tissue tables the way MoTrPAC's S1–S6 do, so a shared label is not evidence of a shared run")
    # MW ST003188
    mw = pd.read_csv(EXT / "mw_ST003188" / "samples.csv", dtype=str)
    mw = mw[mw["mouse"].notna() & ~mw["source"].isin(["Pool", "Blank"])].rename(columns={"organ": "tissue"})
    designs["MW_ST003188"] = (mw, ["batch", "mouse"], ["batch"], "Batch from the mwtab SUBJECT_SAMPLE_FACTORS")
    rows_est, rows_n = [], []
    for name, (m, variables, est_vars, note) in designs.items():
        variables = [v for v in variables if v in m.columns and m[v].notna().any()]
        est_vars = [v for v in est_vars if v in variables]
        tab = batch.nesting_table(m, "tissue", variables)
        tab.insert(0, "dataset", name)
        tab.to_csv(OUT / f"nesting_{name}.csv", index=False)
        rows_n.append(tab)
        pairs, total = batch.estimable_pairs(m, "tissue", est_vars) if est_vars else ([], 0)
        n_t = int(m["tissue"].nunique())
        rows_est.append({"dataset": name, "n_samples": int(len(m)), "n_tissues": n_t, "n_pairs_total": n_t * (n_t - 1) // 2, "n_pairs_estimable": len(pairs), "variables_used": ";".join(est_vars),
                         "frac_pairs_estimable": len(pairs) / (n_t * (n_t - 1) / 2) if n_t > 1 else np.nan, "note": note,
                         "estimable_pairs": ";".join(f"{a}|{b}" for a, b in pairs)})
        print(f"  {name}: {len(m)} samples, {n_t} tissues; estimable pairs {len(pairs)}/{n_t * (n_t - 1) // 2} on {est_vars}")
    est = pd.DataFrame(rows_est)
    est.to_csv(OUT / "estimable_pairs.csv", index=False)
    nest = pd.concat(rows_n, ignore_index=True)
    nest.to_csv(OUT / "nesting_all.csv", index=False)
    # Geiger: recorded as unavailable
    unavailable = pd.DataFrame([{"dataset": "Geiger2013", "reason": "one pooled sample per tissue, no per-sample run/batch variable in the supplementary tables"}])
    unavailable.to_csv(OUT / "unavailable.csv", index=False)
    # comparison with MoTrPAC (frozen phase 16)
    comp = []
    fz = pd.read_csv(FROZEN / "estimable_pairs.csv")
    for _, r in fz.iterrows():
        nz = pd.read_csv(FROZEN / f"nesting_{r['assay']}.csv") if (FROZEN / f"nesting_{r['assay']}.csv").exists() else None
        main_var = r["variables_used"].split(";")[0]
        vrow = nz[nz["variable"] == main_var].iloc[0] if nz is not None and (nz["variable"] == main_var).any() else None
        comp.append({"dataset": f"MoTrPAC {r['assay']}", "batch_variable": main_var, "n_samples": int(r["n_samples"]), "n_tissues": int(r["n_tissues"]), "n_levels": int(vrow["n_levels"]) if vrow is not None else np.nan,
                     "max_tissues_per_level": int(vrow["max_tissues_per_level"]) if vrow is not None else np.nan, "cramers_v": float(vrow["cramers_v"]) if vrow is not None else np.nan,
                     "n_pairs_estimable": int(r["n_pairs_estimable"]), "n_pairs_total": int(r["n_pairs_total"]), "source": "results_frozen/16_identifiability"})
    for _, r in est.iterrows():
        v = r["variables_used"].split(";")[0] if r["variables_used"] else ""
        vrow = nest[(nest["dataset"] == r["dataset"]) & (nest["variable"] == v)]
        vrow = vrow.iloc[0] if len(vrow) else None
        comp.append({"dataset": r["dataset"], "batch_variable": v, "n_samples": r["n_samples"], "n_tissues": r["n_tissues"], "n_levels": int(vrow["n_levels"]) if vrow is not None else np.nan,
                     "max_tissues_per_level": int(vrow["max_tissues_per_level"]) if vrow is not None else np.nan, "cramers_v": float(vrow["cramers_v"]) if vrow is not None else np.nan,
                     "n_pairs_estimable": r["n_pairs_estimable"], "n_pairs_total": r["n_pairs_total"], "source": "results_multiomic/06_external_identifiability"})
    comp = pd.DataFrame(comp)
    comp["frac_pairs_estimable"] = comp["n_pairs_estimable"] / comp["n_pairs_total"].replace(0, np.nan)
    comp.to_csv(OUT / "design_comparison.csv", index=False)
    print(comp.to_string(index=False))
    jr = comp[comp["dataset"] == "Jiang2020"].iloc[0]
    jn = nest[(nest["dataset"] == "Jiang2020") & (nest["variable"] == "tmt_run")].iloc[0]
    jd = nest[(nest["dataset"] == "Jiang2020") & (nest["variable"] == "donor")].iloc[0]
    crossed = comp[(comp["cramers_v"] < 0.999) & (comp["n_pairs_estimable"] > 0)]
    lines = [f"# Phase 6 — identifiability of the external designs (metadata only)", "", f"Built by `scripts/multiomic/06_external_identifiability.py` on {time.strftime('%Y-%m-%d %H:%M UTC', time.gmtime())}; numbers from the CSVs here and from `results_frozen/16_identifiability/`.", "",
             "## Design comparison — `design_comparison.csv`", "", report.df_to_md(comp, floatfmt=".3f"), "",
             "## Per variable — `nesting_all.csv`", "", report.df_to_md(nest, floatfmt=".3f"), "",
             "## Estimable pairs — `estimable_pairs.csv`", "", report.df_to_md(est.drop(columns=["estimable_pairs"]), floatfmt=".3f"), "",
             "Not assessable: " + "; ".join(f"{r.dataset}: {r.reason}" for r in unavailable.itertuples()) + " (`unavailable.csv`).", "",
             "## Reading", "",
             f"- **Jiang 2020 is the positive counterexample.** Its TMT runs each hold up to {int(jn['max_tissues_per_level'])} different tissues (Cramér's V between run and tissue {f(jn['cramers_v'])}, {int(jn['n_levels'])} runs, "
             f"{int(jr['n_pairs_estimable'])} of {int(jr['n_pairs_total'])} tissue pairs share a run); a donor also spans up to {int(jd['max_levels_per_tissue']) if 'max_levels_per_tissue' in jd else 'n/a'} runs. Tissue is therefore separable from plex inside that study — the same TMT chemistry as MoTrPAC, a different allocation of samples to plexes.",
             "- MoTrPAC proteomics (plex_id, V = 1.000, 0 of 21 pairs), MoTrPAC RNA-seq (1 of 171 pairs), Wang 2019 (one run per tissue), MW ST003188 (one batch per organ) are the nested designs.",
             "- Nesting is a choice of design, not a property of the assay: TMT can hold ten tissues in one plex, and when it does the audit's confound disappears.", ""]
    (OUT / "README.md").write_text("\n".join(lines))
    sec = ["- question · does any external dataset have tissue crossed with its processing batch (pre-registration f), and what does that say about the audit's conclusion?",
           f"- data · sample metadata of Jiang 2020 ({int(jr['n_samples'])} samples), Wang 2019, Sato 2022, MW ST003188; MoTrPAC layers from `results_frozen/16_identifiability/` (`results_multiomic/06_external_identifiability/design_comparison.csv`, `nesting_all.csv`, `estimable_pairs.csv`).",
           "- design · `tfp.batch.nesting_table` and `estimable_pairs` unchanged: levels per tissue, tissues per level, Cramér's V with tissue, tissue pairs sharing a level of every batch variable.",
           f"- result · **(f) {'PASS' if (crossed['dataset'] == 'Jiang2020').any() else 'FAIL'}**: " + "; ".join(f"{r.dataset} — {r.batch_variable}: {int(r.n_levels)} levels, up to {int(r.max_tissues_per_level)} tissues per level, Cramér's V {f(r.cramers_v)}, {int(r.n_pairs_estimable)} of {int(r.n_pairs_total)} pairs estimable" for r in crossed.itertuples())
           + f". Nested designs: " + "; ".join(f"{r.dataset} ({r.batch_variable}: V {f(r.cramers_v)}, {int(r.n_pairs_estimable)}/{int(r.n_pairs_total)})" for r in comp[~comp.index.isin(crossed.index)].itertuples())
           + ". MoTrPAC IMMUNO (Luminex plates holding several tissues) was already the audit's one crossed layer; Sato's ROUND is a per-tissue-table label and is listed, not claimed.",
           "- what it does not show · that MoTrPAC's fingerprint is or is not batch; only that the confound is a property of MoTrPAC's sample-to-plex allocation, not of TMT proteomics. Sato's ROUND labels repeat across tissue tables like MoTrPAC's S1–S6 and are not evidence of shared runs; Geiger has no per-sample batch variable.",
           "- files · `results_multiomic/06_external_identifiability/README.md`."]
    (OUT / "REPORT_SECTION.md").write_text("\n".join(sec) + "\n")
    finds = []
    if (crossed["dataset"] == "Jiang2020").any():
        r = crossed[crossed["dataset"] == "Jiang2020"].iloc[0]
        finds.append({"rank": 2, "text": f"**An external TMT design with tissue crossed with plex exists.** In Jiang 2020 each TMT run holds up to {int(r['max_tissues_per_level'])} tissues from several donors ({int(r['n_levels'])} runs, {int(r['n_samples'])} samples, "
                                          f"Cramér's V run × tissue {f(r['cramers_v'])}, {int(r['n_pairs_estimable'])} of {int(r['n_pairs_total'])} tissue pairs estimable within a run) versus MoTrPAC PROT (V = 1.000, 0 of 21). "
                                          f"Nesting is a choice of design, not a property of the assay. — `results_multiomic/06_external_identifiability/design_comparison.csv`"})
    (OUT / "FINDINGS.json").write_text(json.dumps(finds, indent=1))
    (OUT / "STATUS.json").write_text(json.dumps({"status": "DONE", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
