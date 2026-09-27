# Phase 40 — product validation

Written by `scripts/40_product_validation.py`; no model is fit. The within_all BodyMap probabilities reproduce `results_frozen/31_site_regen/12_bodymap/scores_target_probs.csv` (max |Δp| 6.7e-16).

- `scaling.csv`: BodyMap under three scalings (whole set z-scored together / each organ alone / MoTrPAC reference z-scoring).
- `flag_rates.csv`: false-Mismatch rate on correct labels and detection of 1,000 simulated label swaps (seed 20260927), α ∈ {0.05, 0.10, 0.20}; in-study (held-out k20 scores, phase-06 design) and BodyMap (the product model).
- `venacv_cases.csv`: held-out vena cava vials called brown fat, their sets and flag status under the claim 'vena cava'.
- `pca_*.csv`: the reference map (MoTrPAC z-space of the model, 2 components).
