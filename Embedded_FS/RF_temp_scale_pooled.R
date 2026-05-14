# RF_temp_scale_pooled.R
# -----------------------------------------------------------------------------
# Curiosity-only sanity check for Fig. 11's calibration pattern.
#
# Pools out-of-fold (OOF) calibration predictions across ALL (seed x dAb)
# combinations into one large calibration set (~3000 predictions), fits a
# single global per-class T_k, then applies it to each (seed x dAb) held-out
# final_test prediction set and reports pooled ECE / NLL before vs after.
#
# This is NOT a calibration recipe to deploy -- pooling across dAbs averages
# over structurally different failure modes (e.g. dAb16's Delta collapse vs
# dAb5's near-perfect curves). Its only purpose is to test whether the
# WT-under-confident / Beta-over-confident pattern in Fig. 11 survives at large
# sample size, or whether it was a small-N artifact.
#
# Outputs:
#   Embedded_FS/rf_temp_scale_pooled_T.csv   -- the four global T_k values
#   Console: pooled ECE/NLL/accuracy before vs after, plus T_k SDs from the
#   per-(seed,dAb) fits for comparison.
# -----------------------------------------------------------------------------

.libPaths("C:/Users/Zach/Documents/R/win-library/4.1")
setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")

library(tidyverse)
library(tidymodels)
library(ranger)

source("Embedded_FS/rf_tidymodels.R")
source("Visualization_and_preprocessing/S_protein_screen_preprocessing.R")
# Utility functions are redefined inline below to avoid re-running the
# 40-fit sweep that lives at the top level of RF_temp_scale.R.

VARIANT_LEVELS <- c("WT", "Beta", "Delta", "Omicron")
TOP_DABS <- c("dAb5", "dAb6", "dAb8", "dAb13",
              "dAb16", "dAb17", "dAb18", "dAb19")
SEEDS    <- c(1969, 1986, 2024, 2025, 2026)

fit_per_class_T <- function(probs_df, truth, classes) {
  p_mat     <- as.matrix(probs_df[, paste0(".pred_", classes), drop = FALSE])
  truth_idx <- match(as.character(truth), classes)
  logits    <- log(pmax(p_mat, 1e-10))
  nll <- function(log_T) {
    T_k    <- exp(log_T)
    scaled <- sweep(logits, 2, log(T_k), `-`)
    scaled <- scaled - apply(scaled, 1, max)
    p      <- exp(scaled); p <- p / rowSums(p)
    -mean(log(pmax(p[cbind(seq_along(truth_idx), truth_idx)], 1e-10)))
  }
  fit <- optim(rep(0, length(classes)), nll, method = "BFGS",
               control = list(maxit = 200))
  setNames(exp(fit$par), classes)
}

apply_per_class_T <- function(probs_df, T_k, classes) {
  p_mat  <- as.matrix(probs_df[, paste0(".pred_", classes), drop = FALSE])
  logits <- log(pmax(p_mat, 1e-10))
  scaled <- sweep(logits, 2, log(T_k), `-`)
  scaled <- scaled - apply(scaled, 1, max)
  p_cal  <- exp(scaled); p_cal <- p_cal / rowSums(p_cal)
  colnames(p_cal) <- paste0(".pred_", classes)
  out <- probs_df
  out[, paste0(".pred_", classes)] <- as_tibble(p_cal)
  out$.pred_class <- factor(classes[max.col(p_cal)], levels = classes)
  out
}

expected_calibration_error <- function(probs_df, truth, classes, n_bins = 10) {
  p_mat     <- as.matrix(probs_df[, paste0(".pred_", classes), drop = FALSE])
  truth_idx <- match(as.character(truth), classes)
  pred_idx  <- max.col(p_mat)
  conf      <- p_mat[cbind(seq_len(nrow(p_mat)), pred_idx)]
  correct   <- as.integer(pred_idx == truth_idx)
  bins <- cut(conf, breaks = seq(0, 1, length.out = n_bins + 1),
              include.lowest = TRUE, labels = FALSE)
  N <- length(conf); ece <- 0
  for (b in seq_len(n_bins)) {
    in_bin <- which(bins == b); if (length(in_bin) == 0) next
    ece <- ece + (length(in_bin) / N) *
                 abs(mean(conf[in_bin]) - mean(correct[in_bin]))
  }
  ece
}

negloglik_one <- function(probs_df, truth, classes) {
  p_mat     <- as.matrix(probs_df[, paste0(".pred_", classes), drop = FALSE])
  truth_idx <- match(as.character(truth), classes)
  -mean(log(pmax(p_mat[cbind(seq_along(truth_idx), truth_idx)], 1e-10)))
}

# -----------------------------------------------------------------------------
# Re-run the pipeline once but only collect (1) pooled OOF preds and
# (2) per-fit final_test predictions (already in rf_temp_scale_results.rds,
# but we'd rather not depend on a separate file). Same ~2 min runtime.
# -----------------------------------------------------------------------------

dab_list <- main_preprocessor()

oof_preds_all  <- list()
test_preds_all <- list()

for (seed in SEEDS) {
  for (dab in TOP_DABS) {
    cat(sprintf("[seed=%d  dab=%s]\n", seed, dab))

    res <- run_rf_fs(dab_list[[dab]], seed = seed)
    target        <- "target_vector"
    training_data <- training(res$final_split)
    final_test    <- testing(res$final_split)
    best_features <- res$best_features
    formula       <- as.formula(
                       paste(target, "~", paste(best_features, collapse = " + ")))

    set.seed(seed + 1L)
    folds <- vfold_cv(training_data, v = 5, strata = target)

    rf_spec <- rand_forest(mtry  = min(res$best_mtry, length(best_features)),
                           min_n = res$best_min_n,
                           trees = 500) |>
      set_engine("ranger", seed = seed) |>
      set_mode("classification")

    rf_wf <- workflow() |>
      add_recipe(recipe(formula, data = training_data) |> step_zv(all_predictors())) |>
      add_model(rf_spec)

    cv_fit <- fit_resamples(rf_wf, folds,
                            control = control_resamples(save_pred = TRUE),
                            metrics = metric_set(accuracy, roc_auc))
    cv_preds <- collect_predictions(cv_fit)

    oof_preds_all[[length(oof_preds_all) + 1]] <- cv_preds |>
      select(target_vector, starts_with(".pred_")) |>
      mutate(seed = seed, dab = dab, .before = 1)

    final_probs <- augment(res$final_fit, new_data = final_test) |>
      select(target_vector, starts_with(".pred_")) |>
      mutate(seed = seed, dab = dab, .before = 1)
    test_preds_all[[length(test_preds_all) + 1]] <- final_probs
  }
}

oof_pool  <- bind_rows(oof_preds_all)
test_pool <- bind_rows(test_preds_all)

cat(sprintf("\nPooled OOF calibration set: %d predictions\n",  nrow(oof_pool)))
cat(sprintf("Pooled test set:           %d predictions\n", nrow(test_pool)))

# -----------------------------------------------------------------------------
# Fit one global T_k on the pooled OOF predictions
# -----------------------------------------------------------------------------

T_k_global <- fit_per_class_T(oof_pool, oof_pool$target_vector, VARIANT_LEVELS)

cat("\n=== Global per-class T_k (fit on", nrow(oof_pool), "pooled OOF preds) ===\n")
print(round(T_k_global, 3))

write.csv(tibble(variant = names(T_k_global), T_global = T_k_global),
          "Embedded_FS/rf_temp_scale_pooled_T.csv", row.names = FALSE)

# -----------------------------------------------------------------------------
# Apply global T_k to each (seed, dAb) final_test, compute pooled metrics
# -----------------------------------------------------------------------------

test_pool_cal <- apply_per_class_T(test_pool, T_k_global, VARIANT_LEVELS)

raw_ece <- expected_calibration_error(test_pool,     test_pool$target_vector,     VARIANT_LEVELS)
cal_ece <- expected_calibration_error(test_pool_cal, test_pool_cal$target_vector, VARIANT_LEVELS)
raw_nll <- negloglik_one(test_pool,     test_pool$target_vector,     VARIANT_LEVELS)
cal_nll <- negloglik_one(test_pool_cal, test_pool_cal$target_vector, VARIANT_LEVELS)
raw_acc <- mean(as.character(test_pool$.pred_class)     == as.character(test_pool$target_vector))
cal_acc <- mean(as.character(test_pool_cal$.pred_class) == as.character(test_pool_cal$target_vector))

cat("\n=== Pooled-test metrics with global T_k ===\n")
print(tibble(
  raw_ece = raw_ece, cal_ece = cal_ece, ece_delta = cal_ece - raw_ece,
  raw_nll = raw_nll, cal_nll = cal_nll, nll_delta = cal_nll - raw_nll,
  raw_acc = raw_acc, cal_acc = cal_acc
), width = Inf)

# -----------------------------------------------------------------------------
# Sanity check: if the pattern in Fig. 11 was real, expect T_WT > 1
# (under-confident -> sharpen -> wait, actually under-confident wants T < 1
# to push extreme probs further). Note direction:
#   under-confident => raw probs too uniform => want to sharpen => T < 1
#   over-confident  => raw probs too peaked  => want to flatten => T > 1
# Fig. 11 read: WT under-confident -> expect T_WT < 1; Beta over-confident in
# 0.2-0.4 -> expect T_Beta > 1.
# -----------------------------------------------------------------------------
cat("\nDirection check (T<1 = sharpen, T>1 = flatten):\n")
for (k in VARIANT_LEVELS) {
  direction <- ifelse(T_k_global[k] < 1, "sharpen", "flatten")
  cat(sprintf("  %-7s T = %.3f  -> %s\n", k, T_k_global[k], direction))
}
