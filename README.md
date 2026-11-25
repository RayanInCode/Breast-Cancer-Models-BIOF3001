# Breast-Cancer-Models-BIOF3001
Evaluation of multiple breast cancer prognostic models

## GBSG survival models and results (RayanInCode)

Compares penalised Cox, random survival forest, XGBoost AFT and survival SVM
on GBSG (686 patients).

- [Notebook](scripts/gbsg_survival_models.ipynb) · [HTML results](results/survival_comparison/gbsg_survival_models.html)
- [All outputs](results/survival_comparison/): figures, tables, predictions and fitted models.
- [Reliability summary](results/survival_comparison/reliability/reliability_summary.md): bootstrap, calibration and SHAP analyses.

To rerun with Python 3.12:

```bash
python3.12 -m venv .venv-survival
.venv-survival/bin/python -m pip install -r requirements_survival.txt
.venv-survival/bin/python -m ipykernel install --user --name biof3001-survival
```

Select the `biof3001-survival` Jupyter kernel and run the notebook top to bottom.
Methods, references and limitations are documented in the notebook.

## Cox versus random survival forest (Abigail Luk)

`scripts/luk_cox_vs_rsf_clinical_subset.R` compares an unpenalized Cox model with a random survival forest on GBSG2, on all eight predictors, on the clinical subset (age, tumour size, nodes, grade, hormonal therapy), and on the baseline subset (age, menopause, tumour size, nodes, grade). Notes and saved results are in `results/luk_cox_vs_rsf_clinical_subset/`. From the repository root:

```bash
Rscript scripts/luk_cox_vs_rsf_clinical_subset.R
```

`thdata_row_split.csv` in that folder is this script's own 80/20 draw of `TH.data::GBSG2`, stored by row number. It is separate from `results/survival_comparison/patient_split.csv`, which is the group holdout stored by patient id.
