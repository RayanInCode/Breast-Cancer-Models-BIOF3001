# GBSG2: Cox and random survival forest

Script owner: Abigail Luk.

In this repository the script is `scripts/gbsg_cox_rsf.R`. The table and figures in this folder are the saved test-set results. From the repository root, run `Rscript scripts/gbsg_cox_rsf.R`. A fresh run writes a new `outputs/` folder in the working directory.

The script to mark is `gbsg_cox_rsf.R`. It compares a Cox model with a random survival forest on the German Breast Cancer Study Group cohort (GBSG2, 686 patients, 299 events). Each model is fitted to all eight predictors and to a smaller clinical set. The same 80/20 split is used for all four fits.

`Breast-Cancer-Models-BIOF3001/` is a local copy of the group repository. It is not the script for this comparison. The test patients in `gbsg_cox_rsf.R` come from a new stratified draw of `TH.data::GBSG2`, because that table has no patient id.

## What the script does

1. Loads GBSG2. `cens = 1` is recurrence or death. `cens = 0` is censored. Follow-up (`time`) is in days.
2. Draws an 80/20 split, stratified on the event, with seed 42. That is 548 training patients (239 events) and 138 test patients (60 events). The row ids are saved in `outputs/split.csv`.
3. Fits four models on the training patients only.
4. Scores those models on the test patients at 3 years (day 1095).

Full predictor set: age, menopausal status, tumour size, grade, number of positive nodes, progesterone receptor, oestrogen receptor, hormonal therapy.

Clinical subset: age, tumour size, positive nodes, grade, hormonal therapy. The two receptor columns are left out because, in this 1984-89 cohort, they were a quantitative cytosol assay rather than part of the routine pathology report. Menopausal status is not in this subset. Hormonal therapy is kept because it is not an assay. It is the treatment recorded in the trial, so it is not a fact known before that treatment was chosen.

The forest uses 500 trees and the default settings in `randomForestSRC`. Cox is an unpenalized model with Efron ties.

## Test-set results

| Model | C-index | AUROC at 3 years | Brier score at 3 years | Calibration slope |
|---|---|---|---|---|
| Cox, full | 0.674 | 0.704 | 0.202 | 0.774 |
| Cox, clinical subset | 0.663 | 0.706 | 0.199 | 0.852 |
| Forest, full | 0.709 | 0.719 | 0.200 | 0.543 |
| Forest, clinical subset | 0.701 | 0.724 | 0.198 | 0.709 |

A C-index or AUROC of 0.5 is chance and 1 is perfect. The clinical subset is close to the full record on both ranking measures. The forest ranks the test patients a little higher than Cox. The Brier score is the squared error of the predicted 3-year event risk. The four models sit together, near 0.20.

The calibration slope should be 1. The Cox slopes are 0.77 and 0.85. The forest slopes are lower, 0.54 on the full record and 0.71 on the clinical subset, so the forest risk differences are too large on the test patients. Standard errors are in `outputs/metrics.csv`.

## Figures

`outputs/dca_cox_subset.png` is the decision curve at 3 years for the clinical-subset Cox model. The other two lines are treat-all and treat-none. The plotted range of thresholds is 5% to 50%.

`outputs/shap_rsf_full.png` shows SHAP values for the full-record forest. Forty of the 138 test patients are explained. A positive value pushes that patient's risk score up. The largest mean absolute values are positive nodes and progesterone receptor, then tumour size, hormonal therapy, and oestrogen receptor (`outputs/shap_importance.csv`).

## How to run

From this folder, with R 4.5:

```bash
Rscript gbsg_cox_rsf.R
```

Packages: `TH.data`, `survival`, `randomForestSRC`, `rsample`, `timeROC`, `pec`, `dcurves`, `shapr`, `ggplot2`. The SHAP step takes about three minutes. The metrics table is written again to `outputs/metrics.csv`.
