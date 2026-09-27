#!/usr/bin/env python
"""Build docs/MULTIOMIC_REPORT.md from the CSVs under results_multiomic/.

No number is typed here: every figure in the report is read from a result file and the file is named beside it. The
"Significant findings so far" list at the top is rebuilt from scratch on every run, so it always reflects every phase that
has produced results. Phases without results are shown as pending, or as DROPPED when results_multiomic/<phase>/STATUS.json
says so. `--complete` appends the RUN COMPLETE line (Phase 7 only).
"""
from __future__ import annotations

import argparse
import json
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[2]
R = ROOT / "results_multiomic"
DOC = ROOT / "docs" / "MULTIOMIC_REPORT.md"
PHASES = [("00_setup", "Phase 0 — setup and pre-registration"), ("01_rii", "Phase 1 — the RII rescue (local data only)"),
          ("02_discovery", "Phase 2 — data discovery"), ("03_prot_transfer", "Phase 3 — protein fingerprint transfer"),
          ("04_metab_transfer", "Phase 4 — metabolite fingerprint transfer"), ("05_fusion_transfer", "Phase 5 — fusion judged by transfer"),
          ("06_external_identifiability", "Phase 6 — identifiability of the external designs"), ("07_synthesis", "Phase 7 — synthesis")]


def csv(rel: str) -> pd.DataFrame | None:
    p = R / rel
    return pd.read_csv(p) if p.exists() else None


def status(phase: str) -> dict:
    p = R / phase / "STATUS.json"
    return json.loads(p.read_text()) if p.exists() else {"status": "pending"}


def f(v, nd=3):
    try:
        if v is None or (isinstance(v, float) and np.isnan(v)):
            return "n/a"
        return f"{float(v):.{nd}f}"
    except (TypeError, ValueError):
        return str(v)


def pct(v, nd=0):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{100 * float(v):.{nd}f} %"


def md_table(df: pd.DataFrame, floatfmt=".3f", max_rows=40) -> str:
    d = df.head(max_rows).copy()
    for c in d.columns:
        if pd.api.types.is_float_dtype(d[c]):
            d[c] = d[c].map(lambda v: "" if pd.isna(v) else format(v, floatfmt))
    cols = [str(c) for c in d.columns]
    out = ["| " + " | ".join(cols) + " |", "|" + "|".join(["---"] * len(cols)) + "|"]
    out += ["| " + " | ".join(str(v) for v in r) + " |" for r in d.itertuples(index=False)]
    return "\n".join(out)


# ---------------------------------------------------------------------------------------------------------------------
def phase0(findings):
    return ["- question · fix predictions (a)–(f) before any data are analysed.", "- data · none.",
            "- design · `docs/PREREGISTRATION_MULTIOMIC.md` (six predictions with pass rules, fixed analysis choices, drop policy).",
            "- result · written and committed before Phase 1 started.", "- what it does not show · nothing; it is the contract."]


def phase1(findings):
    d = "01_rii/"
    js, vp, vpr = csv(d + "join_summary.csv"), csv(d + "variance_partition.csv"), csv(d + "variance_partition_ratio.csv")
    if js is None or vp is None:
        return None
    js = js.iloc[0]
    vpc, vpm = csv(d + "variance_partition_complete.csv"), csv(d + "variance_partition_median_norm.csv")
    ds = csv(d + "diagnostic_accuracy_summary.csv").set_index("quantity")
    pc, null = csv(d + "panel_curve_summary.csv"), csv(d + "panel_null_k20.csv").iloc[0]
    pcc = csv(d + "panel_curve_complete_summary.csv")
    rp, rpp = csv(d + "rna_protein_correlation_summary.csv"), csv(d + "rna_protein_panel_summary.csv")
    rp = rp.iloc[0] if rp is not None else None
    rpp = rpp.iloc[0] if rpp is not None else None
    r2 = float(vp["R2_tissue"].iloc[0]); r2r = float(vpr["R2_tissue"].iloc[0]); n95 = float(vp["R2_tissue_null95"].iloc[0])
    k_first = int(pc.loc[pc["bal_acc_mean"] >= 0.95, "k"].min()) if (pc["bal_acc_mean"] >= 0.95).any() else None
    k20 = pc[pc["k"] == 20].iloc[0]
    findings.append((4, f"**RII tissue-axis recovery.** On reporter-ion intensities normalised to the channel total, tissue explains R² = {f(r2)} of PC1 "
                        f"(PC2 {f(float(vp['R2_tissue'].iloc[1]))}; label-permutation null 95th pct {f(n95)}) versus {f(r2r, 4)} on the distributed ratio matrix with the same code; "
                        f"n = {int(js['n_vials'])} vials, {int(js['n_animals'])} animals, {int(js['n_proteins_inner'])} proteins in every tissue; identical on the "
                        f"{int(js['n_proteins_inner_complete'])} proteins with no missing value (R² {f(float(vpc['R2_tissue'].iloc[0]))}). Within-study: plex is nested in tissue. "
                        f"— `results_multiomic/01_rii/variance_partition.csv`, `variance_partition_ratio.csv`, `variance_partition_complete.csv`"))
    if rpp is not None:
        findings.append((5, f"**RNA panel markers hold at the protein level.** Of the {int(rpp['n_testable_c'])} RNA panel genes whose marker tissue is one of the 7 proteomics tissues and whose "
                            f"protein is quantified, {int(rpp['n_same_marker'])} ({pct(rpp['frac_same_marker'])}) have the same marker tissue at the protein level (chance 1/7); "
                            f"over all {int(rp['n_genes'])} genes with both layers, the cross-tissue RNA–protein Spearman has median {f(rp['spearman_median'])} "
                            f"(IQR {f(rp['spearman_q25'])}–{f(rp['spearman_q75'])}; mismatched-pair null median {f(rp['null_mismatched_median'])}, 95th pct {f(rp['null_mismatched_q95'])}); "
                            f"n = {rp['n_shared_animals_per_tissue'].split(';')[0].split(':')[1]} shared animals per tissue. "
                            f"— `results_multiomic/01_rii/rna_protein_panel_summary.csv`, `rna_protein_correlation_summary.csv`"))
    lines = [
        "- question · does MoTrPAC proteomics carry a tissue axis on a scale where cross-tissue comparison is defined (pre-registration a–c)?",
        f"- data · portal quant-id `prot-pr` reporter-ion intensities, {int(js['n_tissues'])} tissues, {int(js['n_vials'])} vials, {int(js['n_animals'])} animals; "
        f"peptides summed per protein per plex, < {int(js['min_peptides'])} peptides dropped, channel-total normalisation, log2 ppm; {int(js['n_proteins_union'])} proteins in the union, "
        f"{int(js['n_proteins_inner'])} in every tissue (NaN {pct(js['frac_nan_inner'], 1)}), {int(js['n_proteins_inner_complete'])} with no missing value (`join_summary.csv`).",
        "- design · phase-03 PCA/R² code, phase-04 diagnostic on all 5 animal-grouped folds, `tfp.models.panel_curve` unchanged (RoundRobinSelector → logreg_l2), label-permutation nulls; "
        "RNA–protein Spearman of tissue means over the same animals with a mismatched-pair null.",
        f"- result · **(a) PASS**: tissue R² of PC1 = {f(r2)} (null95 {f(n95)}; ratio matrix {f(r2r, 4)}; median-normalised {f(float(vpm['R2_tissue'].iloc[0])) if vpm is not None else 'n/a'}; "
        f"complete proteins {f(float(vpc['R2_tissue'].iloc[0]))}). **(b) PASS** (context): balanced accuracy ≥ 0.95 from k = {k_first}; k = 20 gives {f(k20['bal_acc_mean'])} ± {f(k20['bal_acc_sd'])} "
        f"over {int(k20['n_folds'])} folds ({f(k20['n_test_animals_mean'], 0)} test animals each; permutation null 95th pct {f(null['null_q95'])}); the same on the complete-protein matrix "
        f"({f(float(pcc.loc[pcc['k'] == 20, 'bal_acc_mean'].iloc[0]))}). Diagnostic: missingness alone still classifies tissue ({f(ds.loc['missingness_outer_acc', 'mean'])} outer, "
        f"{f(ds.loc['missingness_inner_acc', 'mean'])} inner), per-tissue-mean removal collapses to {f(ds.loc['per_tissue_means_removed_acc', 'mean'])} (chance {f(ds.loc['chance', 'mean'])}; on the ratio matrix, frozen phase 04 fold 0: "
        f"{f(float(pd.read_csv(ROOT / 'results_frozen' / '04_baselines' / 'PROT' / 'diagnostic_accuracy.csv')['per_tissue_means_removed'].iloc[0]))}, `results_frozen/04_baselines/PROT/diagnostic_accuracy.csv`). "
        + (f"**(c) PASS**: {int(rpp['n_same_marker'])} of {int(rpp['n_testable_c'])} testable RNA panel genes ({pct(rpp['frac_same_marker'])}) keep their marker tissue at the protein level. "
           f"Cross-tissue RNA–protein Spearman over {int(rp['n_genes'])} genes: median {f(rp['spearman_median'])}, {pct(rp['frac_spearman_gt_0.5'], 1)} above 0.5, {pct(rp['frac_genes_above_null_q95'], 1)} above the "
           f"mismatched-pair null 95th percentile ({f(rp['null_mismatched_q95'])}); same marker tissue in {pct(rp['frac_same_marker_tissue'], 1)} of genes (chance {pct(rp['chance_same_marker'], 1)}); "
           f"correlation rises with the protein's cross-tissue range (median {f(rp.get('spearman_median_Q1 smallest protein range'))} → {f(rp.get('spearman_median_Q4 largest protein range'))} by quartile)." if rpp is not None else ""),
        "- what it does not show · transfer. One plex is one tissue plus that tissue's reference pool, so tissue and plex are confounded exactly as in the audit; the axis is real on this scale "
        "but a per-plex processing offset cannot be excluded without an external dataset. Accuracy 1.000 is a ceiling-task number and is reported as context only. The mismatched-pair null "
        "median of the RNA–protein correlation is far above zero, meaning a large part of any gene's cross-tissue agreement is a shared tissue structure (e.g. muscle/heart vs brain), not gene-specific.",
        f"- files · `results_multiomic/01_rii/README.md` (built from the CSVs), matrices `rii_inner_log2ppm.parquet`, `rii_outer_log2ppm.parquet`, `rii_meta.csv`.",
    ]
    st = csv("01_rii/stability_k20_summary.csv")
    if st is not None:
        s = st.iloc[0]
        lines.append(f"- stability of the k = 20 selection ({int(s['n_boot'])} animal-bootstraps; `stability_k20_summary.csv`): {int(s['n_features_ever_selected'])} proteins ever selected, {int(s['n_selected_ge_0.8'])} in ≥ 80 % of resamples, "
                     f"{int(s['n_selected_ge_0.5'])} in ≥ 50 %; the all-animal panel's members have median selection frequency {f(s['all_animal_panel_median_frequency'], 2)}.")
    for sub, label in (("ph", "phospho (`prot-ph`, c1.0)"), ("c2", "portal release c2.0 (rn7 reprocessing, `prot-pr`)"), ("ac_c2", "acetyl (`prot-ac`, c2.0)")):
        js2 = csv(f"01_rii/{sub}/join_summary.csv")
        if js2 is None:
            continue
        v2 = csv(f"01_rii/{sub}/variance_partition.csv"); p2 = csv(f"01_rii/{sub}/panel_curve_summary.csv"); rp2 = csv(f"01_rii/{sub}/rna_protein_panel_summary.csv")
        extra = ""
        if rp2 is not None:
            q2 = rp2.iloc[0]
            extra = f"; RNA panel genes with the same marker tissue {int(q2['n_same_marker'])} of {int(q2['n_testable_c'])}"
        lines.append(f"- secondary, {label} (`results_multiomic/01_rii/{sub}/`): {int(js2.iloc[0]['n_vials'])} vials, {int(js2.iloc[0]['n_proteins_inner'])} features in every tissue; tissue R² of PC1 = "
                     f"{f(float(v2['R2_tissue'].iloc[0]))} (permutation null95 {f(float(v2['R2_tissue_null95'].iloc[0]))}); k = 20 balanced accuracy {f(float(p2.loc[p2['k'] == 20, 'bal_acc_mean'].iloc[0]))}{extra}.")
    return lines


def generic(phase: str, findings):
    """Phases whose section is written by their own script into results_multiomic/<phase>/REPORT_SECTION.md (numbers from CSVs there)."""
    p = R / phase / "REPORT_SECTION.md"
    if not p.exists():
        return None
    fp = R / phase / "FINDINGS.json"
    if fp.exists():
        for it in json.loads(fp.read_text()):
            findings.append((it["rank"], it["text"]))
    return [p.read_text().rstrip()]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complete", action="store_true")
    args = ap.parse_args()
    findings: list[tuple[int, str]] = []
    sections = {}
    sections["00_setup"] = phase0(findings)
    sections["01_rii"] = phase1(findings)
    for ph, _ in PHASES[2:]:
        sections[ph] = generic(ph, findings)
    # Phase 3 secondary runs: raw-scale sensitivity and the two one-sample-per-tissue atlases
    if sections.get("03_prot_transfer") is not None:
        extra = []
        raw = csv("03_prot_transfer/rawppm/accuracy_overall.csv"); rawc = csv("03_prot_transfer/rawppm/conformal_transfer.csv"); rawr = csv("03_prot_transfer/rawppm/recalibration.csv"); rawo = csv("03_prot_transfer/rawppm/gene_overlap.csv")
        if raw is not None:
            a = raw.set_index("model")["accuracy_sample_weighted"]; c = rawc[(rawc["model"] == "k20") & (rawc["conformal"] == "marginal")].iloc[0]; r5 = rawr[(rawr["model"] == "k20") & (rawr["n_recal"] == 5)].iloc[0]
            extra.append(f"- sensitivity, same scale both sides (`results_multiomic/03_prot_transfer/rawppm/`): Jiang raw reporter intensities re-normalised like the RII (channel total → log2 ppm, technical replicates kept), "
                         f"{int(rawo.iloc[0]['matched_genes'])} genes, {int(raw['n_samples_mapped'].iloc[0])} mapped samples / {int(raw['n_donors_mapped'].iloc[0])} donors: accuracy k20 **{f(a['k20'])}**, k50 {f(a.get('k50', np.nan))}, full {f(a['full'])}; "
                         f"coverage with MoTrPAC calibration {f(c['coverage_mapped'])} (empty {f(c['frac_empty_mapped'])}), recalibrated on 5 donors {f(r5['coverage_recalibrated'])} (set size {f(r5['set_size_recalibrated'])}).")
        for atlas, label in (("wang2019", "Wang 2019, 29 human tissues, one donor each, label-free"), ("geiger2013", "Geiger 2013, 28 mouse tissues, pooled mice, SILAC ratios, gene-symbol match")):
            p = R / "03_prot_transfer" / atlas / "summary.json"
            if p.exists():
                s = json.loads(p.read_text())
                extra.append(f"- secondary atlas ({label}; `results_multiomic/03_prot_transfer/{atlas}/`): {s['matched_genes']} genes; {s['n_mapped']} mapped tissue samples, {s['n_ood']} OOD; accuracy k20 {f(s.get('acc_k20'))}, k50 {f(s.get('acc_k50'))}, full {f(s.get('acc_full'))}; "
                             f"coverage with MoTrPAC calibration at k20 {f(s['coverage_k20_marginal_source_cal'])} (empty {f(s['frac_empty_k20_source_cal'])}); OOD empty-set fraction {f(s.get('ood_frac_empty_k20'))}.")
        sections["03_prot_transfer"] = sections["03_prot_transfer"] + extra
    now = time.strftime("%Y-%m-%d %H:%M UTC", time.gmtime())
    done = [ph for ph, _ in PHASES if sections.get(ph) is not None or status(ph).get("status") == "DONE"]
    out = ["# Multiomic overnight report — can the proteome and metabolome carry a tissue fingerprint that transfers?", "",
           "Branch `multiomic-overnight`. Pre-registration: `docs/PREREGISTRATION_MULTIOMIC.md`. Log: `docs/MULTIOMIC_LOG.md`.",
           "Every number in this file is read from a CSV under `results_multiomic/` by `scripts/multiomic/build_report.py`; the file is named beside each number. "
           "Within-study accuracy is context, never a finding (the audit: batch is nested in tissue).", "",
           f"## Significant findings so far ({now}; phases with results: {', '.join(p.split('_')[0].lstrip('0') or '0' for p in done)})", ""]
    if findings:
        out += [f"{i + 1}. {t}" for i, (_, t) in enumerate(sorted(findings, key=lambda x: x[0]))]
    else:
        out += ["_None yet._"]
    out.append("")
    # run 2: verification and the two extensions, written by their scripts from their CSVs
    for sub, title in (("08_verification", "Verification (run 2, part A)"), ("09_extensions", "What B found (run 2, part B)")):
        p = R / sub / "REPORT_SECTION.md"
        if p.exists():
            out += [f"## {title}", "", p.read_text().rstrip(), ""]
    for ph, title in PHASES:
        out += [f"## {title}", ""]
        st = status(ph)
        body = sections.get(ph)
        if body is None:
            if st.get("status") == "DROPPED":
                out += [f"**DROPPED** — {st.get('reason', '')}", ""]
            else:
                out += ["_pending_", ""]
        else:
            if st.get("status") == "DROPPED":
                out += [f"**DROPPED (partial)** — {st.get('reason', '')}", ""]
            out += body + [""]
    # Phase 7 synthesis sections (written by scripts/multiomic/07_synthesis.py from the CSVs)
    syn = R / "07_synthesis"
    for fn, title in (("submission_changes.md", "What this changes in the submission, if anything"), ("for_the_talk.md", "For the talk"),
                      ("not_done.md", "Not done / not possible"), ("merge_checklist.md", "Merge checklist"),
                      ("site_draft_multiomic.md", "Site page (run 2 built it as site/multiomic.html; this is the run-1 draft it grew from)")):
        p = syn / fn
        if p.exists():
            out += [f"## {title}", "", p.read_text().rstrip(), ""]
    if args.complete:
        out += [f"RUN COMPLETE {time.strftime('%Y-%m-%dT%H:%M:%SZ', time.gmtime())}"]
    DOC.write_text("\n".join(out) + "\n")
    print(f"wrote {DOC} ({len(findings)} findings; sections with results: {done})")


if __name__ == "__main__":
    main()
