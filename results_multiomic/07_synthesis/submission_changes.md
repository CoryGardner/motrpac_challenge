Quoted sentences are from the submission as found (`README.md`, `site/limitations.html`); the replacement is what the overnight results support. **Nothing in the submission was edited.** Numbers: the CSV named in each item.

1. `README.md` (Scope): *"proteomics and metabolomics evaluated for within-tissue use only."*
   → *"proteomics evaluated cross-tissue on the portal reporter-ion scale (tissue R² of PC1 0.991 vs 0.0009 on the distributed ratios; `results_multiomic/01_rii/variance_partition.csv`) and transferred to a human proteome atlas (k20 accuracy 0.455 on 44 samples / 13 donors, coverage 0.068 under MoTrPAC calibration; `results_multiomic/03_prot_transfer/`); metabolomics transferred to two external mouse metabolomes by RefMet name (deep-platform panel: k20 0.637 on 490 samples / 70 mice; `results_multiomic/04_metab_transfer/legs_summary.csv`). Within-study, tissue is still nested in plex for both."*

2. `site/limitations.html`: *"Proteomics is within-tissue only (ratios to per-tissue reference pools; missingness identifies the tissue), so no cross-tissue protein fingerprint is claimed."*
   → *"The distributed proteomics (ratios to per-tissue reference pools) carries no cross-tissue axis; the portal reporter-ion intensities do (tissue R² of PC1 0.991, identical on the 2393 proteins with no missing value), and a protein panel selected on them transfers to Jiang 2020 with the RNA pattern: accuracy above chance, coverage collapse (0.068) under MoTrPAC calibration, recovery to 0.917 after recalibration on 5 donors at 2.95 of 7 classes per set (the GTEx pattern, not the BodyMap one). Missingness still identifies the tissue on that scale (1.000 accuracy from the NaN pattern alone), and plex is nested in tissue."*

3. `README.md` (Roadmap): *"a cross-tissue proteomics fingerprint;"*
   → done on the RII scale (branch `multiomic-overnight`, `results_multiomic/01_rii/`, `03_prot_transfer/`); remaining: a rat multi-tissue proteome as the BodyMap-equivalent (none exists in PRIDE/ProteomeXchange, `results_multiomic/02_discovery/attempts.csv`).

4. `README.md` (numbers table): *"Sedentary vs 8-week-trained within tissue: mean best single-omic AUROC / fusion beats single / attributable to training | 0.994 / 0 of 7 / 4 of 7"*
   → add a row: *"Fusion judged by transfer (Jiang 2020, 42 samples with matched RNA + protein, 12 donors): k20 accuracy RNA 1.000 / protein 0.476 / late fusion 0.738; coverage 0.452 / 0.071 / 0.048 — fusion is not more robust than RNA alone under the species shift"* (`results_multiomic/05_fusion_transfer/fusion_transfer_summary.csv`).

5. `site/identifiability.html` / README ("Estimable tissue pairs within study (RNA-seq) | 1 of 171"): add the external counterexample —
   *"In Jiang 2020's TMT design each run holds up to 8 tissues (56 runs; Cramér's V run × tissue 0.324; 424 of 528 tissue pairs estimable within a run) versus MoTrPAC PROT (V 1.000, 0 of 21): nesting is a choice of design, not a property of the assay."* (`results_multiomic/06_external_identifiability/design_comparison.csv`).

Not changed: every RNA number on the site and in the README; the audit's conclusion for MoTrPAC itself (still nested); the human-transfer numbers (GTEx).
