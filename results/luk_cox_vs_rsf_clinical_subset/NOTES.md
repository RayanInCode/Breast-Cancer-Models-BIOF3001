# GBSG2: Cox and random survival forest

Script owner: Abigail Luk.

The script is `scripts/luk_cox_vs_rsf_clinical_subset.R`. It compares a Cox model with a random survival forest on the German Breast Cancer Study Group cohort (GBSG2, 686 patients, 299 events). Each model is fitted to all eight predictors, to a clinical subset, and to a baseline subset. The same 80/20 split is used for all six fits.

The test patients come from a new stratified draw of `TH.data::GBSG2`, because that table has no patient id. The row numbers are in `thdata_row_split.csv` in this folder. That file is not the group's `patient_split.csv`.

## What the script does

1. Loads GBSG2. `cens = 1` is recurrence or death. `cens = 0` is censored. Follow-up (`time`) is in days.
2. Draws an 80/20 split, stratified on the event, with seed 42. That is 548 training patients (239 events) and 138 test patients (60 events). The row ids are saved in `thdata_row_split.csv`.
3. Fits six models on the training patients only.
4. Scores those models on the test patients at 3 years (day 1095).

Full predictor set: age, menopausal status, tumour size, grade, number of positive nodes, progesterone receptor, oestrogen receptor, hormonal therapy.

Clinical subset: age, tumour size, positive nodes, grade, hormonal therapy. The two receptor columns are left out because, in this 1984-89 cohort, they were a quantitative cytosol assay rather than part of the routine pathology report. Menopausal status is not in this subset. Hormonal therapy is kept because it is not an assay. It is the treatment recorded in the trial, so it is not a fact known before that treatment was chosen.

Baseline subset: age, menopausal status, tumour size, grade, positive nodes. The receptor assays are left out for the same reason as in the clinical subset. Menopausal status is put back because it is recorded in the history and does not need an assay. Hormonal therapy is left out because it is the treatment assigned in the trial. The clinical subset answers whether the record still ranks patients after the assays are removed. The baseline subset answers a different question: what is known from the history and the surgical specimen before that treatment assignment is known.

The forest uses 500 trees and the default settings in `randomForestSRC`. Cox is an unpenalized model with Efron ties.

## Test-set results

| Model | C-index | AUROC at 3 years | Brier score at 3 years | Calibration slope |
|---|---|---|---|---|
| Cox, full | 0.674 | 0.704 | 0.202 | 0.774 |
| Cox, clinical subset | 0.663 | 0.706 | 0.199 | 0.852 |
| Cox, baseline subset | 0.663 | 0.701 | 0.204 | 0.821 |
| Forest, full | 0.709 | 0.719 | 0.200 | 0.543 |
| Forest, clinical subset | 0.701 | 0.724 | 0.198 | 0.709 |
| Forest, baseline subset | 0.695 | 0.714 | 0.204 | 0.673 |

A C-index or AUROC of 0.5 is chance and 1 is perfect. Both smaller sets stay close to the full record on the ranking measures. The clinical subset is the closer of the two: its 3-year AUROC is 0.706 for Cox and 0.724 for the forest, against 0.701 and 0.714 for the baseline subset. Taking out the assigned treatment and putting menopausal status back does not improve on the clinical subset. The forest still ranks the test patients a little higher than Cox within each variable set. The Brier score is the squared error of the predicted 3-year event risk. All six models sit together, near 0.20. With 60 test events, gaps of about 0.01 are not a clear separation.

The calibration slope should be 1. The Cox slopes are 0.77, 0.85, and 0.82. The forest slopes are lower, 0.54 on the full record, 0.71 on the clinical subset, and 0.67 on the baseline subset, so the forest risk differences are too large on the test patients. Standard errors are in `test_set_metrics.csv`.

## Figures

`dca_cox_clinical_subset.png` is the decision curve at 3 years for the clinical-subset Cox model. `dca_cox_baseline_subset.png` is the same curve for the baseline-subset Cox model. On each plot the other two lines are treat-all and treat-none. The plotted range of thresholds is 5% to 50%. The forest is not on these curves.

`shap_rsf_full_record.png` shows SHAP values for the full-record forest. Forty of the 138 test patients are explained. A positive value pushes that patient's risk score up. The largest mean absolute values are positive nodes and progesterone receptor, then tumour size, hormonal therapy, and oestrogen receptor (`shap_mean_abs.csv`).

## How to run

From the repository root, with R 4.5:

```bash
Rscript scripts/luk_cox_vs_rsf_clinical_subset.R
```

Packages: `TH.data`, `survival`, `randomForestSRC`, `rsample`, `timeROC`, `pec`, `dcurves`, `shapr`, `ggplot2`. The SHAP step takes about three minutes. The script writes the table and figures back into this folder.
