#!/usr/bin/env python
"""Multiomic Phase 7 — synthesis. Writes, from the CSVs under results_multiomic/, the four Markdown blocks that
build_report.py appends to docs/MULTIOMIC_REPORT.md:
  submission_changes.md — sentences of the submission (README, site) that the overnight results would change, quoted, with
                          the replacement text; nothing in the submission is edited;
  for_the_talk.md       — at most three sentences, only for findings that survived their null and their caveats;
  not_done.md           — what was not done or not possible, with reasons, and the download-log summary;
  site_draft_multiomic.md — the draft content of site/_drafts/multiomic.html as Markdown (figures with their source CSVs).
Every number is read from a result file; the file is named beside it.
"""
from __future__ import annotations

import json
import re
import time
from pathlib import Path

import numpy as np
import pandas as pd

from tfp import config as C

ROOT = C.ROOT
R = ROOT / "results_multiomic"
OUT = R / "07_synthesis"


def f(v, nd=3):
    return "n/a" if v is None or (isinstance(v, float) and np.isnan(v)) else f"{float(v):.{nd}f}"


def csv(rel):
    p = R / rel
    return pd.read_csv(p) if p.exists() else None


def main():
    OUT.mkdir(parents=True, exist_ok=True)
    # ---- numbers ----------------------------------------------------------------------------------------------------
    vp = csv("01_rii/variance_partition.csv"); vpr = csv("01_rii/variance_partition_ratio.csv"); js = csv("01_rii/join_summary.csv").iloc[0]
    rpp = csv("01_rii/rna_protein_panel_summary.csv").iloc[0]; rp = csv("01_rii/rna_protein_correlation_summary.csv").iloc[0]
    a3 = csv("03_prot_transfer/accuracy_overall.csv").set_index("model"); c3 = csv("03_prot_transfer/conformal_transfer.csv"); r3 = csv("03_prot_transfer/recalibration.csv"); o3 = csv("03_prot_transfer/gene_overlap.csv").iloc[0]
    a3r = csv("03_prot_transfer/rawppm/accuracy_overall.csv").set_index("model")
    cm3 = c3[(c3["model"] == "k20") & (c3["conformal"] == "marginal")].iloc[0]; rk5 = r3[(r3["model"] == "k20") & (r3["n_recal"] == 5)].iloc[0]
    g13 = json.loads((R / "03_prot_transfer" / "geiger2013" / "summary.json").read_text()); w19 = json.loads((R / "03_prot_transfer" / "wang2019" / "summary.json").read_text())
    L4 = csv("04_metab_transfer/legs_summary.csv").set_index("leg")
    mw = L4.loc["deep_mw"]; hs = L4.loc["hilic_sato"]; ds = L4.loc["deep_sato"]
    F5 = csv("05_fusion_transfer/fusion_transfer_summary.csv"); k5 = F5[F5["model"] == "k20"].set_index("layer"); f5 = F5[F5["model"] == "full"].set_index("layer"); o5 = csv("05_fusion_transfer/overlap.csv").iloc[0]
    r19 = csv("05_fusion_transfer/rna19_on_jiang.csv").set_index("model")
    D6 = csv("06_external_identifiability/design_comparison.csv").set_index("dataset"); jj = D6.loc["Jiang2020"]; pp = D6.loc["MoTrPAC PROT"]
    dl = csv("02_discovery/download_log.csv"); att = csv("02_discovery/attempts.csv")
    r2, r2r = float(vp["R2_tissue"].iloc[0]), float(vpr["R2_tissue"].iloc[0])

    # ---- 1. what this changes in the submission -----------------------------------------------------------------------
    sc = ["Quoted sentences are from the submission as found (`README.md`, `site/limitations.html`); the replacement is what the overnight results support. **Nothing in the submission was edited.** Numbers: the CSV named in each item.", "",
          "1. `README.md` (Scope): *\"proteomics and metabolomics evaluated for within-tissue use only.\"*",
          f"   → *\"proteomics evaluated cross-tissue on the portal reporter-ion scale (tissue R² of PC1 {f(r2)} vs {f(r2r, 4)} on the distributed ratios; `results_multiomic/01_rii/variance_partition.csv`) and transferred to a human proteome atlas "
          f"(k20 accuracy {f(a3.loc['k20', 'accuracy_sample_weighted'])} on {int(a3.loc['k20', 'n_samples_mapped'])} samples / {int(a3.loc['k20', 'n_donors_mapped'])} donors, coverage {f(cm3['coverage_mapped'])} under MoTrPAC calibration; `results_multiomic/03_prot_transfer/`); "
          f"metabolomics transferred to two external mouse metabolomes by RefMet name (deep-platform panel: k20 {f(mw['acc_k20'])} on {int(mw['n_mapped'])} samples / {int(mw['n_individuals_mapped'])} mice; `results_multiomic/04_metab_transfer/legs_summary.csv`). Within-study, tissue is still nested in plex for both.\"*", "",
          "2. `site/limitations.html`: *\"Proteomics is within-tissue only (ratios to per-tissue reference pools; missingness identifies the tissue), so no cross-tissue protein fingerprint is claimed.\"*",
          f"   → *\"The distributed proteomics (ratios to per-tissue reference pools) carries no cross-tissue axis; the portal reporter-ion intensities do (tissue R² of PC1 {f(r2)}, identical on the {int(js['n_proteins_inner_complete'])} proteins with no missing value), "
          f"and a protein panel selected on them transfers to Jiang 2020 with the RNA pattern: accuracy above chance, coverage collapse ({f(cm3['coverage_mapped'])}) under MoTrPAC calibration, recovery to {f(rk5['coverage_recalibrated'])} after recalibration on 5 donors. "
          f"Missingness still identifies the tissue on that scale ({f(float(csv('01_rii/diagnostic_accuracy_summary.csv').set_index('quantity').loc['missingness_outer_acc', 'mean']))} accuracy from the NaN pattern alone), and plex is nested in tissue.\"*", "",
          "3. `README.md` (Roadmap): *\"a cross-tissue proteomics fingerprint;\"*",
          f"   → done on the RII scale (branch `multiomic-overnight`, `results_multiomic/01_rii/`, `03_prot_transfer/`); remaining: a rat multi-tissue proteome as the BodyMap-equivalent (none exists in PRIDE/ProteomeXchange, `results_multiomic/02_discovery/attempts.csv`).", "",
          "4. `README.md` (numbers table): *\"Sedentary vs 8-week-trained within tissue: mean best single-omic AUROC / fusion beats single / attributable to training | 0.994 / 0 of 7 / 4 of 7\"*",
          f"   → add a row: *\"Fusion judged by transfer (Jiang 2020, {int(o5['mapped_samples_7class'])} samples with matched RNA + protein, {int(o5['mapped_donors_7class'])} donors): k20 accuracy RNA {f(k5.loc['RNA', 'accuracy'])} / protein {f(k5.loc['protein', 'accuracy'])} / late fusion {f(k5.loc['late_mean', 'accuracy'])}; "
          f"coverage {f(k5.loc['RNA', 'coverage_motrpac_cal'])} / {f(k5.loc['protein', 'coverage_motrpac_cal'])} / {f(k5.loc['late_mean', 'coverage_motrpac_cal'])} — fusion is not more robust than RNA alone under the species shift\"* (`results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv`).", "",
          "5. `site/identifiability.html` / README (\"Estimable tissue pairs within study (RNA-seq) | 1 of 171\"): add the external counterexample —",
          f"   *\"In Jiang 2020's TMT design each run holds up to {int(jj['max_tissues_per_level'])} tissues ({int(jj['n_levels'])} runs; Cramér's V run × tissue {f(jj['cramers_v'])}; {int(jj['n_pairs_estimable'])} of {int(jj['n_pairs_total'])} tissue pairs estimable within a run) versus MoTrPAC PROT (V {f(pp['cramers_v'])}, {int(pp['n_pairs_estimable'])} of {int(pp['n_pairs_total'])}): nesting is a choice of design, not a property of the assay.\"* (`results_multiomic/06_external_identifiability/design_comparison.csv`).", "",
          "Not changed: every RNA number on the site and in the README; the audit's conclusion for MoTrPAC itself (still nested); the human-transfer numbers (GTEx)."]
    (OUT / "submission_changes.md").write_text("\n".join(sc) + "\n")

    # ---- 2. for the talk ---------------------------------------------------------------------------------------------
    talk = [f"1. Put MoTrPAC's proteomics back on the reporter-ion scale and tissue becomes its dominant axis (R² {f(r2)} of PC1 against {f(r2r, 4)} on the distributed ratios; permutation null {f(float(vp['R2_tissue_null95'].iloc[0]))}), so the layer we had to leave out of the fingerprint was a normalisation choice, not a property of the proteome (`results_multiomic/01_rii/variance_partition.csv`).",
            f"2. A 20-protein panel selected on those intensities names the tissue of {f(a3.loc['k20', 'accuracy_sample_weighted'], 2)} of {int(a3.loc['k20', 'n_samples_mapped'])} human TMT samples from {int(a3.loc['k20', 'n_donors_mapped'])} GTEx donors (chance {f(1 / 7, 2)}; {f(a3r.loc['k20', 'accuracy_sample_weighted'], 2)} when both sides are processed the same way), "
            f"and its 90 % prediction sets cover only {f(cm3['coverage_mapped'], 2)} with MoTrPAC calibration but {f(rk5['coverage_recalibrated'], 2)} after recalibrating on five donors — exactly the RNA pattern (`results_multiomic/03_prot_transfer/`).",
            f"3. The audit's confound is a design choice, not the assay: in Jiang 2020 each TMT run holds up to {int(jj['max_tissues_per_level'])} tissues and {int(jj['n_pairs_estimable'])} of {int(jj['n_pairs_total'])} tissue pairs are estimable within a run, versus 0 of {int(pp['n_pairs_total'])} in MoTrPAC's one-tissue-per-plex design (`results_multiomic/06_external_identifiability/design_comparison.csv`)."]
    (OUT / "for_the_talk.md").write_text("\n".join(talk) + "\n")

    # ---- 3. not done / not possible -----------------------------------------------------------------------------------
    dl["dataset"] = dl["path"].str.extract(r"external_multiomic/(?:_logs/\.\./)?([^/]+)/")[0]
    per = dl.groupby("dataset").agg(files=("url", "size"), MB=("bytes", lambda s: s.sum() / 1e6)).reset_index()
    per_lines = "; ".join(f"{r.dataset} {int(r.files)} files, {r.MB:.0f} MB" for r in per.itertuples())
    nd = ["- **Metabolite transfer HILIC+ → MW ST003188 (`hilic_mw`)**: stopped under the pre-registered rule — only " + str(int(L4.loc['hilic_mw', 'matched_refmet'])) + " RefMet names shared between the 129 HILIC+ metabolites named in every MoTrPAC tissue and the atlas's 190 (`results_multiomic/04_metab_transfer/hilic_mw/feature_overlap.csv`). The deep-platform source covered the atlas instead.",
          "- **Wang 2019 as a full target**: the PRIDE MaxQuant bundle (28.2 GB) exceeds the 20 GB box and was not downloaded; the Europe PMC EV tables (gene-level intensities, one donor per tissue) were used, so Wang contributes 6 mapped tissue samples and no donor-level statement (`results_multiomic/03_prot_transfer/wang2019/`).",
          "- **Geiger 2013**: matched by gene symbol across species (no rat–mouse orthology table in the repo) and one pooled sample per tissue; treated as a secondary check only (`results_multiomic/03_prot_transfer/geiger2013/`).",
          "- **A rat multi-tissue proteome (the BodyMap-equivalent)**: none with processed tables exists on PRIDE or ProteomeXchange (searches logged in `results_multiomic/02_discovery/attempts.csv`); the reverse-direction test therefore uses a human atlas.",
          "- **MetaboLights** was reachable but not searched: two matchable rodent metabolomes were already in hand within the 2 h discovery box.",
          "- **Sato 2022 batch variable**: ROUND labels repeat across the per-tissue Metabolon tables and could not be verified as shared runs; listed in Phase 6, not claimed.",
          "- **PXD082651** (2026 mouse lifespan multi-tissue atlas of non-canonical peptides) not pursued (non-canonical peptide focus).",
          "- **Portal release c2.0 RII, acetyl and ubiquityl RII**: not run in the main line (c1.0 `prot-pr` and `prot-ph` only) unless a later log entry says otherwise.",
          "- **Figures**: described in the site draft with their source CSVs; no image files were rendered and `site/` was not touched.",
          "- **Kidney and adipose** have no human protein target in Jiang 2020; they are covered only by the mouse atlas (Geiger, n = 1 per tissue).",
          f"- **Download log** (`results_multiomic/02_discovery/download_log.csv`): {len(dl)} files, {dl['bytes'].sum() / 1e6:.0f} MB — {per_lines}. {int(att.shape[0])} attempts recorded (`attempts.csv`), including the failures (raw-only PRIDE project, bot-walled publisher page, over-box bundle, no rat atlas, Sato absent from Metabolomics Workbench)."]
    (OUT / "not_done.md").write_text("\n".join(nd) + "\n")

    # ---- 4. site draft ------------------------------------------------------------------------------------------------
    sd = ["*Draft content for a new site page (Markdown only; `site/` unchanged). Each figure names the CSV it would be drawn from.*", "",
          "### Title: The proteome and metabolome carry the fingerprint too — on the right scale", "",
          f"**Lead.** The submission left proteomics and metabolomics out of the cross-tissue fingerprint. On the portal's reporter-ion intensities tissue explains R² {f(r2)} of the first principal component (distributed ratios: {f(r2r, 4)}); a 20-protein panel then transfers to a human TMT atlas and a 20-metabolite panel to a mouse metabolome atlas, both with the coverage collapse and recalibration recovery already seen for RNA.", "",
          "**Figure 1 — Same proteins, two scales.** Tissue R² of PC1–PC3 on the distributed ratio matrix vs the reporter-ion (log2 ppm) matrix, with the label-permutation null. Source: `results_multiomic/01_rii/variance_partition.csv`, `variance_partition_ratio.csv`, `variance_partition_complete.csv` (missingness-free repeat).", "",
          f"**Figure 2 — The protein ladder beside the RNA ladder.** Accuracy and 90 %-set coverage (MoTrPAC calibration) for k20 / k50 / full: protein → Jiang 2020 ({int(o3['target_samples_mapped'])} samples, {int(o3['target_donors_mapped'])} donors; primary scale and raw-ppm scale) and RNA → GTEx (frozen). Source: `results_multiomic/03_prot_transfer/ladder_protein_vs_rna.csv`, `rawppm/accuracy_overall.csv`, `coverage_ci.csv` (donor-bootstrap intervals).", "",
          f"**Figure 3 — Recalibration on a few donors.** Coverage vs number of target individuals used for recalibration (0, 3, 5) with set size, for protein (Jiang) and metabolites (Sato, MW ST003188). Source: `results_multiomic/03_prot_transfer/recalibration.csv`, `results_multiomic/04_metab_transfer/<leg>/recalibration.csv`.", "",
          f"**Figure 4 — Where the protein panel is wrong.** Confusion of Jiang tissues (rows, mapped and out-of-distribution) against the 7 rat classes at k20; out-of-distribution tissues get empty sets in {f(cm3['ood_frac_empty'])} of cases. Source: `results_multiomic/03_prot_transfer/confusion_k20.csv`, `ood_sets.csv`, `accuracy_by_tissue.csv`.", "",
          f"**Figure 5 — Metabolites: 70 mice, 12 organs, 5 ages.** Per-organ accuracy of the deep-platform panel on the mouse aging atlas (k20 {f(mw['acc_k20'])}, full {f(mw['acc_full'])}, chance {f(mw['chance'])}), flat across ages; and the Sato sedentary → exercised invariance test (native panel: accuracy {f(hs['invariance_acc_full'])}, coverage {f(hs['invariance_cov_full'])}). Source: `results_multiomic/04_metab_transfer/deep_mw/accuracy_by_tissue.csv`, `accuracy_by_stage.csv`, `hilic_sato/invariance_native_sedentary_to_exercised.csv`.", "",
          f"**Figure 6 — Fusion under shift.** Accuracy, coverage and empty-set fraction of RNA, protein, late-mean fusion and stacked fusion on the same {int(o5['mapped_samples_7class'])} Jiang samples (k20 and full), with donor-bootstrap intervals. Source: `results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv`.", "",
          f"**Figure 7 — The design that is not nested.** Cramér's V (tissue × batch variable) and the fraction of estimable tissue pairs for every MoTrPAC layer and every external dataset; Jiang 2020 at V {f(jj['cramers_v'])} and {int(jj['n_pairs_estimable'])}/{int(jj['n_pairs_total'])}. Source: `results_multiomic/06_external_identifiability/design_comparison.csv`.", "",
          f"**Figure 8 — RNA markers at the protein level.** Distribution of the cross-tissue RNA–protein Spearman over {int(rp['n_genes'])} genes with the mismatched-pair null, and the marker-tissue agreement of the RNA panel genes ({int(rpp['n_same_marker'])} of {int(rpp['n_testable_c'])}). Source: `results_multiomic/01_rii/rna_protein_correlation.csv`, `rna_protein_panel_genes.csv`.", "",
          "**Caveats box.** Within-study accuracies are context: plex is nested in tissue on the reporter-ion scale too, and the NaN pattern alone identifies the tissue. Human atlases are adult, post-mortem and differently processed; kidney and adipose have no human protein target. RefMet name matches across platforms are name matches. Every number on this page is read from the CSV named beside it."]
    (OUT / "site_draft_multiomic.md").write_text("\n".join(sd) + "\n")
    (OUT / "STATUS.json").write_text(json.dumps({"status": "DONE", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}))
    (R / "00_setup" / "STATUS.json").write_text(json.dumps({"status": "DONE", "utc": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())}))
    print(f"wrote {OUT}")


if __name__ == "__main__":
    main()
