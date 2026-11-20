# Breast-Cancer-Models-BIOF3001
Evaluation of multiple breast cancer prognostic models

The executed [survival-model notebook](scripts/gbsg_survival_models.ipynb) follows
the supplied 13-step prediction workflow. It compares a survival random forest,
XGBoost AFT and a linear survival SVM with a penalised Cox baseline, using package
implementations. The supplied dataset is `data/gbsg.csv` (686 patients).

Open the [HTML version](results/survival_comparison/gbsg_survival_models.html) to
read the results without running Python. The existing exploratory notebook is
preserved in `scripts/gbsg_exploratory.ipynb`.

To rerun, use Python 3.12, install `requirements_survival.txt`, and select that
Python environment as the Jupyter kernel. An environment is already available at
`.venv-survival/bin/python` on this computer. `requirements_survival_lock.txt`
records the installed dependencies used for execution.

```bash
python3.12 -m venv .venv-survival
.venv-survival/bin/python -m pip install -r requirements_survival.txt
.venv-survival/bin/python -m ipykernel install --user --name biof3001-survival
```

Run the notebook from top to bottom. It creates
`results/survival_comparison/`, containing figures, validation tables, holdout
predictions and fitted models with preprocessing. It reserves 138 patients for
the holdout and performs nested cross-validation on the other 548. All four saved
models are fitted to the development set; the holdout is retained for evaluation.

The outcome is time to recurrence **or death**, with right censoring. SVM provides
risk rankings, so its probability-based metrics and decision curves are omitted.
Paper and package references, model assumptions and limitations are documented
in the notebook. These are new GBSG fits informed by the papers, rather than
reproductions of the authors' fitted models.

The additional reliability analysis includes 2,000 paired holdout bootstrap
comparisons with Cox, calibration intervals, Brier skill against a development-set
Kaplan–Meier reference, and 50 bootstrap development refits with fixed settings.
SHAP explains the same 50 holdout patients for all four models, using 30 development
patients as the background. Probability models explain three-year event risk;
SVM explains its ranking score. Results and plots are saved in
`results/survival_comparison/reliability/`. The analyses preserve the original
model choice and are exploratory extensions of the already inspected holdout.
Read the [reliability conclusions](results/survival_comparison/reliability/reliability_summary.md)
for the results of these checks.

## Cox versus random survival forest (Abigail Luk)

`scripts/luk_cox_vs_rsf_clinical_subset.R` compares an unpenalized Cox model with a random survival forest on GBSG2, once with all eight predictors and once with the clinical subset (age, tumour size, nodes, grade, hormonal therapy). Notes and saved results are in `results/luk_cox_vs_rsf_clinical_subset/`. From the repository root:

```bash
Rscript scripts/luk_cox_vs_rsf_clinical_subset.R
```

`thdata_row_split.csv` in that folder is this script's own 80/20 draw of `TH.data::GBSG2`, stored by row number. It is separate from `results/survival_comparison/patient_split.csv`, which is the group holdout stored by patient id.
