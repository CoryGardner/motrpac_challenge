# Figures for slides

Rendered by `scripts/multiomic/figures_png.py` from `site/data/multiomic.json` (the same numbers as `site/multiomic.html`; provenance in `site/data/provenance.json`, ids `mo_*`).

- `fig0_pca_two_scales.png`: the proteomics problem in one picture: PC1 vs PC2 per vial, as distributed (each sample ÷ a same-tissue reference pool) and rebuilt from the reporter-ion intensities (each sample ÷ its own total signal); from `results_multiomic/01_rii/pca_scores_ratio.csv` and `pca_scores_rii.csv` (`scripts/multiomic/01b_pca_scores.py`; render: `python scripts/multiomic/figures_png.py --only fig0`).
- `fig1_same_proteins_two_scales.png`
- `fig2_protein_ladder_beside_rna.png`
- `fig3_recalibration_coverage_and_set_size.png`
- `fig4_fusion_under_shift.png`
- `fig5_design_comparison.png`
- `fig6_rna_markers_at_protein_level.png`
