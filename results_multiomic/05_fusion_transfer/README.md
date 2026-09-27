# Phase 5 — fusion judged by transfer (Jiang 2020, matched RNA + protein)

Built by `scripts/multiomic/05_fusion_transfer.py` on 2026-09-27 07:41 UTC; numbers from the CSVs here.

- 182 Jiang samples (13 donors, 32 tissues) have both layers; **42 from 12 donors map to the 7-class fingerprint** (CORTEX, HEART, LIVER, LUNG, SKM-GN present; KIDNEY, WAT-SC absent). RNA: 13800 matched genes; protein: 2731 (`overlap.csv`).
- MoTrPAC 7-tissue source: RNA 350 vials / 50 animals; RII 420 vials / 60 animals; calibration on 15 animals held out of both layers (105 animal × tissue rows).

## Accuracy, coverage and empty sets on the same mapped Jiang samples — `fusion_transfer_summary.csv`

| model | layer | n_mapped | n_donors | accuracy | acc_ci95_low_donor_boot | acc_ci95_high_donor_boot | coverage_motrpac_cal | cov_ci95_low_donor_boot | cov_ci95_high_donor_boot | frac_empty | avg_set_size | ood_frac_empty | n_cal_rows | cal_accuracy |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| k20 | RNA | 42 | 12 | 1.000 | 1.000 | 1.000 | 0.452 | 0.348 | 0.581 | 0.548 | 0.452 | 0.979 | 105 | 1.000 |
| k20 | protein | 42 | 12 | 0.476 | 0.386 | 0.571 | 0.071 | 0.000 | 0.150 | 0.929 | 0.071 | 0.964 | 105 | 1.000 |
| k20 | late_mean | 42 | 12 | 0.738 | 0.611 | 0.857 | 0.048 | 0.000 | 0.119 | 0.952 | 0.048 | 0.979 | 105 | 1.000 |
| k20 | stacked_LR | 42 | 12 | 0.738 | 0.611 | 0.857 | 0.048 | 0.000 | 0.119 | 0.952 | 0.048 | 0.979 | 105 | 1.000 |
| k50 | RNA | 42 | 12 | 1.000 | 1.000 | 1.000 | 0.238 | 0.133 | 0.342 | 0.762 | 0.238 | 0.979 | 105 | 1.000 |
| k50 | protein | 42 | 12 | 0.476 | 0.381 | 0.568 | 0.048 | 0.000 | 0.119 | 0.952 | 0.048 | 0.979 | 105 | 1.000 |
| k50 | late_mean | 42 | 12 | 0.833 | 0.703 | 0.933 | 0.048 | 0.000 | 0.119 | 0.952 | 0.048 | 0.979 | 105 | 1.000 |
| k50 | stacked_LR | 42 | 12 | 0.833 | 0.703 | 0.933 | 0.048 | 0.000 | 0.119 | 0.952 | 0.048 | 0.979 | 105 | 1.000 |
| full | RNA | 42 | 12 | 1.000 | 1.000 | 1.000 | 0.024 | 0.000 | 0.075 | 0.976 | 0.024 | 1.000 | 105 | 1.000 |
| full | protein | 42 | 12 | 0.524 | 0.395 | 0.636 | 0.048 | 0.000 | 0.119 | 0.952 | 0.048 | 1.000 | 105 | 1.000 |
| full | late_mean | 42 | 12 | 1.000 | 1.000 | 1.000 | 0.048 | 0.000 | 0.119 | 0.952 | 0.048 | 1.000 | 105 | 1.000 |
| full | stacked_LR | 42 | 12 | 1.000 | 1.000 | 1.000 | 0.048 | 0.000 | 0.119 | 0.952 | 0.048 | 1.000 | 105 | 1.000 |

Rule (pre-registration): a fusion is 'more robust' only if accuracy AND coverage ≥ the better single layer (RNA at k20) and empty-set fraction ≤ it. At k20: late_mean NO, stacked_LR NO.

Per tissue — `fusion_accuracy_by_tissue.csv`:

| model | rat_class | n | RNA | late_mean | protein | stacked_LR |
|---|---|---|---|---|---|---|
| full | CORTEX | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| full | HEART | 17 | 1.00 | 1.00 | 0.00 | 1.00 |
| full | LIVER | 5 | 1.00 | 1.00 | 1.00 | 1.00 |
| full | LUNG | 8 | 1.00 | 1.00 | 1.00 | 1.00 |
| full | SKM-GN | 10 | 1.00 | 1.00 | 0.70 | 1.00 |
| k20 | CORTEX | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| k20 | HEART | 17 | 1.00 | 0.35 | 0.00 | 0.35 |
| k20 | LIVER | 5 | 1.00 | 1.00 | 0.80 | 1.00 |
| k20 | LUNG | 8 | 1.00 | 1.00 | 0.88 | 1.00 |
| k20 | SKM-GN | 10 | 1.00 | 1.00 | 0.70 | 1.00 |
| k50 | CORTEX | 2 | 1.00 | 1.00 | 1.00 | 1.00 |
| k50 | HEART | 17 | 1.00 | 0.59 | 0.00 | 0.59 |
| k50 | LIVER | 5 | 1.00 | 1.00 | 0.20 | 1.00 |
| k50 | LUNG | 8 | 1.00 | 1.00 | 1.00 | 1.00 |
| k50 | SKM-GN | 10 | 1.00 | 1.00 | 0.90 | 1.00 |

## Confusion structure across layers — `confusion_structure.csv`

| model | where | n_offdiag_cells | pearson_hard_offdiag | pearson_soft_offdiag | per_tissue_accuracy_pearson | rna_oof_accuracy | protein_oof_accuracy |
|---|---|---|---|---|---|---|---|
| k20 | Jiang 2020 (transfer) | 30 |  | 0.259 |  |  |  |
| full | Jiang 2020 (transfer) | 30 |  | 0.065 |  |  |  |
| k20 | MoTrPAC 7 tissues (animal-grouped out-of-fold) | 42 |  | 0.618 |  | 1.000 | 1.000 |
| full | MoTrPAC 7 tissues (animal-grouped out-of-fold) | 42 |  | 0.806 |  | 1.000 | 1.000 |

Off-diagonal correlation of the row-normalised confusion matrices (hard = argmax calls; soft = mean class probability per true tissue). On Jiang the true classes are the 5 mapped tissues; locally all 7.

## The 19-class RNA fingerprint on the same Jiang RNA — `rna19_on_jiang.csv`

| model | n_mapped_19class | n_donors | accuracy | n_tissues_mapped | rna_genes |
|---|---|---|---|---|---|
| k20 | 76 | 13 | 0.526 | 13 | 14492 |
| k50 | 76 | 13 | 0.763 | 13 | 14492 |
| full | 76 | 13 | 0.816 | 13 | 14492 |

## What this does not show

- 42 samples from 12 donors, five classes, two of them tiny (cortex, liver): the fusion comparison is directional, not a precise estimate.
- The protein layer here is the authors' cleaned relative abundance; Phase 3's raw-scale sensitivity gave higher protein accuracy, so the protein arm is conservative.

_Run time 0.2 min._