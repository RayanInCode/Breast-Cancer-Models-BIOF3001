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

Compares an unpenalised Cox model with a random survival forest on GBSG2
(686 patients): the full record, a clinical subset, and a baseline subset.

- [Script](scripts/luk_cox_vs_rsf_clinical_subset.R) · [Notes](results/luk_cox_vs_rsf_clinical_subset/NOTES.md)
- [All outputs](results/luk_cox_vs_rsf_clinical_subset/): test-set metrics, the row-number split, decision curves and SHAP.
- Code added for this comparison is that script only. The notebooks in `scripts/` and everything under `results/survival_comparison/` were already in the repository.

To rerun from the repository root, with R 4.5:

```bash
Rscript scripts/luk_cox_vs_rsf_clinical_subset.R
```

The script writes back into `results/luk_cox_vs_rsf_clinical_subset/`.
`thdata_row_split.csv` is this script's own 80/20 draw of `TH.data::GBSG2`.
It is separate from `results/survival_comparison/patient_split.csv`.
