# Objective
# Compare Cox proportional hazards with a random survival forest on GBSG2.
# Each model is fitted twice: once with all eight predictors, and once with
# the clinical subset. All four fits use the same training patients and are
# scored on the same test patients.
#
# Reported on the test set: Harrell's C-index, time-dependent AUROC at 3 years
# (1095 days), Brier score at 3 years, calibration slope, a decision curve at
# 3 years for the clinical-subset Cox model, and SHAP values for the
# full-record forest.
#
# Script owner: add your name.
# Coding assistance was used to draft this file. Paste the prompts in README.md
# before submission.

# 3 years, in days. GBSG2 records follow-up in days.
HORIZON <- 1095
SEED <- 42
NTREE <- 500

suppressPackageStartupMessages({
  library(survival)
  library(randomForestSRC)
  library(pec)
  library(rsample)
  library(timeROC)
  library(dcurves)
  library(shapr)
  library(ggplot2)
  library(TH.data)
})

dir.create("outputs", showWarnings = FALSE)

# GBSG2 from TH.data: 686 patients, 299 events.
# cens = 1 is the event (recurrence or death). cens = 0 is censored.
# time is the follow-up in days.
data(GBSG2, package = "TH.data")
gbsg <- GBSG2
stopifnot(nrow(gbsg) == 686, sum(gbsg$cens) == 299)

# horTh and menostat are already factors. tgrade arrives as an ordered factor,
# which would enter the Cox model as a polynomial contrast. It is reset to an
# unordered factor so grade II and grade III are separate indicators, with
# grade I as the reference.
gbsg$horTh <- factor(gbsg$horTh, levels = c("no", "yes"))
gbsg$menostat <- factor(gbsg$menostat, levels = c("Pre", "Post"))
gbsg$tgrade <- factor(gbsg$tgrade, levels = c("I", "II", "III"), ordered = FALSE)
gbsg$event <- factor(gbsg$cens, levels = c(0, 1))
gbsg$horizon_days <- HORIZON
gbsg$row_id <- seq_len(nrow(gbsg))

# Full set: every predictor in the file.
# Clinical subset: age, tumour size, positive nodes, grade, and hormonal
# therapy. Progesterone receptor (progrec) and oestrogen receptor (estrec)
# are left out. In this 1984-89 cohort those values were a quantitative
# cytosol assay, not a line on a routine pathology report. Menopausal status
# is also left out of this subset. horTh is the treatment recorded in the
# trial. It is included here because it is not an assay, not because it is
# known before the treatment decision.
full_vars <- c("age", "menostat", "tsize", "tgrade", "pnodes", "progrec", "estrec", "horTh")
clinical_vars <- c("age", "tsize", "pnodes", "tgrade", "horTh")

formula_of <- function(vars) {
  as.formula(paste("Surv(time, cens) ~", paste(vars, collapse = " + ")))
}

# One 80/20 split, stratified on the event. prop = 0.8 is the training share.
# The same row ids are used for every model below. Seed 42 makes the draw
# repeatable. GBSG2 has no patient id, so these 138 test patients are the
# patients from this draw, not a split that was stored by id.
set.seed(SEED)
split <- initial_split(gbsg, prop = 0.8, strata = event)
train <- training(split)
test <- testing(split)
stopifnot(length(intersect(train$row_id, test$row_id)) == 0)

message(sprintf(
  "Split: %d train (%d events), %d test (%d events).",
  nrow(train), sum(train$cens), nrow(test), sum(test$cens)
))

write.csv(
  data.frame(
    row_id = c(train$row_id, test$row_id),
    sample = rep(c("train", "test"), c(nrow(train), nrow(test)))
  ),
  "outputs/split.csv",
  row.names = FALSE
)

# Fitted on the training patients only. pec needs the Cox model frame, so
# x = TRUE and y = TRUE are stored with the fit. The forest uses 500 trees
# and the package defaults for the other settings. Its seed is fixed so the
# trees do not change between runs.
cox_full <- coxph(formula_of(full_vars), data = train, ties = "efron", x = TRUE, y = TRUE)
cox_subset <- coxph(formula_of(clinical_vars), data = train, ties = "efron", x = TRUE, y = TRUE)

rsf_full <- rfsrc(
  formula_of(full_vars), data = train,
  ntree = NTREE, seed = SEED, importance = FALSE
)
rsf_subset <- rfsrc(
  formula_of(clinical_vars), data = train,
  ntree = NTREE, seed = SEED, importance = FALSE
)

# Survival at day 1095 from the step function: the last recorded time at or
# before the horizon. For a Cox model, surv is one column per patient.
# For the forest, survival is one row per patient.
surv_at_cox <- function(fit, newdata, horizon) {
  sf <- survfit(fit, newdata = newdata)
  idx <- findInterval(horizon, sf$time)
  if (idx < 1) {
    return(rep(1, nrow(newdata)))
  }
  if (is.matrix(sf$surv)) sf$surv[idx, ] else rep(sf$surv[idx], nrow(newdata))
}

surv_at_rsf <- function(fit, newdata, horizon) {
  pred <- predict(fit, newdata = newdata)
  idx <- findInterval(horizon, pred$time.interest)
  if (idx < 1) {
    return(rep(1, nrow(newdata)))
  }
  pred$survival[, idx]
}

# Event risk, not survival. Decision curves and the Brier score treat a
# larger number as a higher chance of the event.
risk_of <- function(surv) pmin(pmax(1 - surv, 0), 1)

# Ranking scores. A larger value means higher risk.
# Cox: linear predictor. Forest: mortality from randomForestSRC.
# concordance() otherwise reads a larger value as longer survival, so
# reverse = TRUE is set when the C-index is calculated.
rank_score <- list(
  "Cox, full" = predict(cox_full, newdata = test, type = "lp"),
  "Cox, clinical subset" = predict(cox_subset, newdata = test, type = "lp"),
  "RSF, full" = predict(rsf_full, newdata = test)$predicted,
  "RSF, clinical subset" = predict(rsf_subset, newdata = test)$predicted
)

risk_3y <- list(
  "Cox, full" = risk_of(surv_at_cox(cox_full, test, HORIZON)),
  "Cox, clinical subset" = risk_of(surv_at_cox(cox_subset, test, HORIZON)),
  "RSF, full" = risk_of(surv_at_rsf(rsf_full, test, HORIZON)),
  "RSF, clinical subset" = risk_of(surv_at_rsf(rsf_subset, test, HORIZON))
)

# Harrell's C-index on the test patients.
c_index_of <- function(risk) {
  concordance(Surv(time, cens) ~ risk, data = test, reverse = TRUE)$concordance
}

# Time-dependent AUROC at day 1095. The marker is the risk score.
# Weighting is the marginal Kaplan-Meier weights from timeROC.
auroc_of <- function(risk) {
  roc <- timeROC(
    T = test$time, delta = test$cens, marker = risk,
    cause = 1, weighting = "marginal", times = HORIZON, iid = FALSE
  )
  as.numeric(roc$AUC[length(roc$AUC)])
}

# Brier score at day 1095 on the test patients. splitMethod = "none" means
# pec does not refit the models; it only scores the fits made above.
# The censoring distribution is the marginal Kaplan-Meier.
brier_fit <- pec(
  object = list(
    "Cox, full" = cox_full,
    "Cox, clinical subset" = cox_subset,
    "RSF, full" = rsf_full,
    "RSF, clinical subset" = rsf_subset
  ),
  formula = Surv(time, cens) ~ 1,
  data = test,
  times = HORIZON,
  cens.model = "marginal",
  exact = FALSE,
  splitMethod = "none"
)

brier_of <- function(name) {
  err <- brier_fit$AppErr[[name]]
  as.numeric(err[length(err)])
}

# Calibration slope on the test patients. Target value is 1.
# Cox: coefficient of the linear predictor in a Cox model fitted to the
# test outcome. Under proportional hazards this slope is not tied to one day.
# Forest: the score is the cumulative hazard at day 1095. It is multiplied by
# the coefficient from the same regression on the training patients, which
# forces the training slope to 1. The test coefficient is then the slope.
# A slope below 1 means the risk differences are too large on the test set.
calibration_of <- function(score_test, score_train = NULL) {
  scored <- test
  scored$lp <- score_test
  if (!is.null(score_train)) {
    trained <- train
    trained$lp <- score_train
    train_slope <- coef(coxph(Surv(time, cens) ~ lp, data = trained))[["lp"]]
    scored$lp <- score_test * as.numeric(train_slope)
  }
  fit <- coxph(Surv(time, cens) ~ lp, data = scored)
  c(slope = unname(coef(fit)[["lp"]]), se = sqrt(vcov(fit)[["lp", "lp"]]))
}

rsf_chf <- function(fit, newdata) {
  pred <- predict(fit, newdata = newdata)
  idx <- max(1, findInterval(HORIZON, pred$time.interest))
  pred$chf[, idx]
}

cal <- list(
  "Cox, full" = calibration_of(rank_score[["Cox, full"]]),
  "Cox, clinical subset" = calibration_of(rank_score[["Cox, clinical subset"]]),
  "RSF, full" = calibration_of(rsf_chf(rsf_full, test), rsf_chf(rsf_full, train)),
  "RSF, clinical subset" = calibration_of(
    rsf_chf(rsf_subset, test), rsf_chf(rsf_subset, train)
  )
)

metrics <- data.frame(
  model = names(rank_score),
  variable_set = c("full", "clinical subset", "full", "clinical subset"),
  n_train = nrow(train),
  n_test = nrow(test),
  events_test = sum(test$cens),
  c_index = vapply(rank_score, c_index_of, numeric(1)),
  auroc_3y = vapply(rank_score, auroc_of, numeric(1)),
  brier_3y = vapply(names(rank_score), brier_of, numeric(1)),
  calibration_slope = vapply(cal, function(x) x["slope"], numeric(1)),
  calibration_slope_se = vapply(cal, function(x) x["se"], numeric(1)),
  row.names = NULL
)

write.csv(metrics, "outputs/metrics.csv", row.names = FALSE)
message("\nTest-set metrics at ", HORIZON, " days:")
print(
  metrics[, c("model", "c_index", "auroc_3y", "brier_3y", "calibration_slope")],
  digits = 3, row.names = FALSE
)

# Decision curve at 3 years for the clinical-subset Cox model.
# cox_subset_risk_3y is already a probability of the event by day 1095.
# It is not passed through as_probability, because that option would refit
# a Cox model and would not use this probability.
# "Treat all" and "Treat none" are added by dca().
dca_data <- test
dca_data$cox_subset_risk_3y <- risk_3y[["Cox, clinical subset"]]
dca_fit <- dca(
  Surv(time, cens) ~ cox_subset_risk_3y,
  data = dca_data,
  time = HORIZON,
  thresholds = seq(0.05, 0.50, by = 0.01),
  label = list(cox_subset_risk_3y = "Cox, clinical subset")
)
dca_plot <- plot(dca_fit) +
  ggplot2::labs(
    title = "Decision curve at 3 years, test set",
    subtitle = "Cox model using age, size, nodes, grade, and tamoxifen"
  )
ggplot2::ggsave("outputs/dca_cox_subset.png", dca_plot, width = 7, height = 5, dpi = 150)

# SHAP for the full-record forest.
# randomForestSRC is not a built-in shapr model, so predict_model returns the
# mortality score. fastshap does not install on R 4.5, which is why shapr is used.
# Forty test patients are explained, drawn with seed 42. All 138 would require
# a separate prediction for every feature coalition. phi0 is the mean mortality
# on the training patients, the prediction when no feature is known.
predict_rsf <- function(model, newdata) {
  as.numeric(predict(model, newdata = as.data.frame(newdata))$predicted)
}

x_full <- full_vars
set.seed(SEED)
explain_rows <- sample(seq_len(nrow(test)), size = min(40, nrow(test)))
phi0 <- mean(predict_rsf(rsf_full, train[, x_full]))

shap_fit <- explain(
  model = rsf_full,
  x_explain = test[explain_rows, x_full],
  x_train = train[, x_full],
  approach = "independence",
  phi0 = phi0,
  n_MC_samples = 25,
  seed = SEED,
  verbose = "basic",
  predict_model = predict_rsf
)

shap_values <- as.data.frame(shap_fit$shapley_values_est)
feature_cols <- setdiff(names(shap_values), c("explain_id", "none"))
mean_abs <- sort(
  vapply(shap_values[feature_cols], function(z) mean(abs(z)), numeric(1)),
  decreasing = TRUE
)
top5 <- names(mean_abs)[seq_len(min(5, length(mean_abs)))]

# Colour is the raw value of that variable for the same patient.
# Factors are coded as their level number so they can share one colour scale.
shap_long <- do.call(rbind, lapply(top5, function(var) {
  data.frame(
    variable = var,
    shap = shap_values[[var]],
    value_num = as.numeric(test[[var]][explain_rows]),
    stringsAsFactors = FALSE
  )
}))
shap_long$variable <- factor(shap_long$variable, levels = rev(top5))

shap_plot <- ggplot(shap_long, aes(x = shap, y = variable, colour = value_num)) +
  geom_vline(xintercept = 0, colour = "grey70") +
  geom_jitter(height = 0.15, width = 0, size = 1.6, alpha = 0.85) +
  scale_colour_viridis_c(option = "C", name = "Value") +
  labs(
    title = "SHAP values for the full-record random survival forest",
    subtitle = sprintf(
      "Top 5 variables, %d test patients. A positive value raises the risk score.",
      length(explain_rows)
    ),
    x = "SHAP value (forest mortality)",
    y = NULL
  ) +
  theme_minimal(base_size = 12)

ggplot2::ggsave("outputs/shap_rsf_full.png", shap_plot, width = 7, height = 4.5, dpi = 150)
write.csv(
  data.frame(variable = names(mean_abs), mean_abs_shap = unname(mean_abs)),
  "outputs/shap_importance.csv",
  row.names = FALSE
)

message("\nWrote outputs/metrics.csv")
message("Wrote outputs/dca_cox_subset.png")
message("Wrote outputs/shap_rsf_full.png")
