# RF_temp_scale.R
# -----------------------------------------------------------------------------
# Per-class temperature scaling for the RF classifier on top dAbs.
#
# Pipeline per (seed x dAb):
#   1. Run the standard run_rf_fs to get best_features and a (training, final_test) split.
#   2. Refit RF on training_data using only best_features via 5-fold CV to obtain
#      out-of-fold (OOF) probability predictions. These OOF preds are the calibration
#      set -- they are not used to fit the final classifier, so no leakage.
#   3. Fit four scalar temperatures T_k (one per variant) by minimising NLL on the
#      OOF predictions: scaled logit = log(p) - log(T_k), then re-softmax.
#   4. Apply T_k to the final_fit's probability predictions on the held-out final_test.
#   5. Compare ECE / NLL / argmax accuracy before vs after.
#
# Repeat across 5 seeds (1969, 1986, 2024, 2025, 2026) and 8 top dAbs to report
# the variability of T_k and of the calibration improvement.
#
# Outputs:
#   Embedded_FS/rf_temp_scale_results.rds   -- per-(seed,dAb) results tibble with raw/cal probs
#   Embedded_FS/rf_temp_scale_summary.csv   -- per-dAb means and SDs
#   Thesis_figures/figs/11b_temp_scale_T_k.png         -- T_k variability boxplots
#   Thesis_figures/figs/11c_temp_scale_ece_before_after.png -- ECE before/after
# -----------------------------------------------------------------------------

.libPaths("C:/Users/Zach/Documents/R/win-library/4.1")
setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")

library(tidyverse)
library(tidymodels)
library(ranger)

source("Embedded_FS/rf_tidymodels.R")
source("Visualization_and_preprocessing/S_protein_screen_preprocessing.R")

VARIANT_LEVELS <- c("WT", "Beta", "Delta", "Omicron")
TOP_DABS <- c("dAb5", "dAb6", "dAb8", "dAb13",
              "dAb16", "dAb17", "dAb18", "dAb19")
SEEDS    <- c(1969, 1986, 2024, 2025, 2026)
N_BINS_ECE <- 10

# -----------------------------------------------------------------------------
# Calibration utilities
# -----------------------------------------------------------------------------

# Fit one T_k per class by NLL minimisation. Parameterise as exp(log_T) to keep
# T_k > 0 without bounded optimisation.
fit_per_class_T <- function(probs_df, truth, classes) {
  p_mat     <- as.matrix(probs_df[, paste0(".pred_", classes), drop = FALSE])
  truth_idx <- match(as.character(truth), classes)
  if (any(is.na(truth_idx))) stop("Some truth labels missing from `classes`.")

  logits <- log(pmax(p_mat, 1e-10))

  nll <- function(log_T) {
    T_k    <- exp(log_T)
    scaled <- sweep(logits, 2, log(T_k), `-`)
    scaled <- scaled - apply(scaled, 1, max)          # numerically stable softmax
    p      <- exp(scaled); p <- p / rowSums(p)
    -mean(log(pmax(p[cbind(seq_along(truth_idx), truth_idx)], 1e-10)))
  }

  fit <- optim(rep(0, length(classes)), nll, method = "BFGS",
               control = list(maxit = 200))
  setNames(exp(fit$par), classes)
}

# Apply learned T_k to a fresh probability tibble; returns a new tibble with the
# same .pred_<class> columns plus an updated .pred_class.
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

# Expected Calibration Error using max-prob bins (Guo et al. 2017, Eq. 3).
expected_calibration_error <- function(probs_df, truth, classes, n_bins = 10) {
  p_mat     <- as.matrix(probs_df[, paste0(".pred_", classes), drop = FALSE])
  truth_idx <- match(as.character(truth), classes)

  pred_idx <- max.col(p_mat)
  conf     <- p_mat[cbind(seq_len(nrow(p_mat)), pred_idx)]
  correct  <- as.integer(pred_idx == truth_idx)

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
# Main loop
# -----------------------------------------------------------------------------

dab_list <- main_preprocessor()
results  <- list()

for (seed in SEEDS) {
  for (dab in TOP_DABS) {
    cat(sprintf("[seed=%d  dab=%s]\n", seed, dab))

    # 1. Standard FS pipeline -------------------------------------------------
    res <- run_rf_fs(dab_list[[dab]], seed = seed)

    target        <- "target_vector"
    training_data <- training(res$final_split)
    final_test    <- testing(res$final_split)
    best_features <- res$best_features
    formula       <- as.formula(
                       paste(target, "~", paste(best_features, collapse = " + ")))

    # 2. 5-fold CV on training_data with best_features -> OOF probability preds
    set.seed(seed + 1L)
    folds <- vfold_cv(training_data, v = 5, strata = target)

    # NB: parsnip's classification mode + ranger engine returns class
    # probabilities via augment automatically. Do NOT pass probability = TRUE
    # explicitly -- that breaks the parsnip wrapper's prediction code path.
    rf_spec <- rand_forest(mtry  = min(res$best_mtry, length(best_features)),
                           min_n = res$best_min_n,
                           trees = 500) |>
      set_engine("ranger", seed = seed) |>
      set_mode("classification")

    rf_wf <- workflow() |>
      add_recipe(recipe(formula, data = training_data) |> step_zv(all_predictors())) |>
      add_model(rf_spec)

    # roc_auc in the metric set forces probability predictions to be stored
    # alongside the class predictions when save_pred = TRUE.
    cv_fit <- fit_resamples(rf_wf, folds,
                            control = control_resamples(save_pred = TRUE),
                            metrics = metric_set(accuracy, roc_auc))
    cv_preds <- collect_predictions(cv_fit)

    # 3. Fit per-class T_k on OOF predictions ---------------------------------
    T_k <- fit_per_class_T(cv_preds, cv_preds[[target]], VARIANT_LEVELS)

    # 4. Apply T_k to final_fit's predictions on the held-out final_test ------
    final_probs     <- augment(res$final_fit, new_data = final_test)
    final_probs_cal <- apply_per_class_T(final_probs, T_k, VARIANT_LEVELS)

    truth   <- final_test[[target]]
    raw_ece <- expected_calibration_error(final_probs,     truth, VARIANT_LEVELS, N_BINS_ECE)
    cal_ece <- expected_calibration_error(final_probs_cal, truth, VARIANT_LEVELS, N_BINS_ECE)
    raw_nll <- negloglik_one(final_probs,     truth, VARIANT_LEVELS)
    cal_nll <- negloglik_one(final_probs_cal, truth, VARIANT_LEVELS)
    raw_acc <- mean(as.character(final_probs$.pred_class)     == as.character(truth))
    cal_acc <- mean(as.character(final_probs_cal$.pred_class) == as.character(truth))

    results[[length(results) + 1]] <- tibble(
      seed = seed, dab = dab,
      n_train = nrow(training_data), n_test = nrow(final_test),
      n_features = length(best_features),
      T_WT = T_k["WT"], T_Beta = T_k["Beta"],
      T_Delta = T_k["Delta"], T_Omicron = T_k["Omicron"],
      raw_ece = raw_ece, cal_ece = cal_ece,
      raw_nll = raw_nll, cal_nll = cal_nll,
      raw_acc = raw_acc, cal_acc = cal_acc,
      raw_probs = list(final_probs |> select(target_vector, starts_with(".pred_"))),
      cal_probs = list(final_probs_cal |> select(target_vector, starts_with(".pred_")))
    )
  }
}

results_df <- bind_rows(results)
saveRDS(results_df, "Embedded_FS/rf_temp_scale_results.rds")

# -----------------------------------------------------------------------------
# Summaries
# -----------------------------------------------------------------------------

per_dab <- results_df |>
  group_by(dab) |>
  summarise(across(c(T_WT, T_Beta, T_Delta, T_Omicron,
                     raw_ece, cal_ece, raw_nll, cal_nll, raw_acc, cal_acc),
                   list(mean = ~mean(.x), sd = ~sd(.x))),
            .groups = "drop")
write.csv(per_dab, "Embedded_FS/rf_temp_scale_summary.csv", row.names = FALSE)

cat("\n=== Per-class T_k (mean +/- SD across all 5 seeds x 8 dAbs = 40 fits) ===\n")
results_df |>
  summarise(across(c(T_WT, T_Beta, T_Delta, T_Omicron),
                   list(mean = ~mean(.x), sd = ~sd(.x), median = ~median(.x)))) |>
  print(width = Inf)

cat("\n=== Calibration improvement (mean +/- SD across 40 fits) ===\n")
results_df |>
  summarise(raw_ece_mean = mean(raw_ece), raw_ece_sd = sd(raw_ece),
            cal_ece_mean = mean(cal_ece), cal_ece_sd = sd(cal_ece),
            raw_nll_mean = mean(raw_nll), cal_nll_mean = mean(cal_nll),
            raw_acc_mean = mean(raw_acc), cal_acc_mean = mean(cal_acc)) |>
  print(width = Inf)

# -----------------------------------------------------------------------------
# Plots
# -----------------------------------------------------------------------------

dir.create("Thesis_figures/figs", showWarnings = FALSE, recursive = TRUE)

T_long <- results_df |>
  pivot_longer(c(T_WT, T_Beta, T_Delta, T_Omicron),
               names_to = "variant", values_to = "T") |>
  mutate(variant = factor(str_remove(variant, "^T_"), levels = VARIANT_LEVELS))

p_T <- ggplot(T_long, aes(variant, T, fill = variant)) +
  geom_hline(yintercept = 1, linetype = 2, color = "grey50") +
  geom_boxplot(alpha = 0.5, outlier.shape = NA, width = 0.5) +
  geom_jitter(width = 0.15, alpha = 0.7, size = 1.6) +
  scale_y_log10() +
  scale_fill_manual(values = c(WT = "#4C72B0", Beta = "#DD8452",
                               Delta = "#55A868", Omicron = "#C44E52"),
                    guide = "none") +
  labs(x = NULL, y = "Learned temperature T_k (log scale)",
       title = "Per-class temperature variability across 5 seeds x 8 top dAbs",
       subtitle = "Each point = one (seed, dAb) fit. T = 1 is no scaling; T > 1 flattens, T < 1 sharpens.") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave("Thesis_figures/figs/11b_temp_scale_T_k.png", p_T,
       width = 9, height = 5, dpi = 300)

ece_long <- results_df |>
  select(seed, dab, raw_ece, cal_ece) |>
  pivot_longer(c(raw_ece, cal_ece), names_to = "stage", values_to = "ece") |>
  mutate(stage = factor(stage, levels = c("raw_ece", "cal_ece"),
                        labels = c("Uncalibrated", "Per-class T")))

p_ece <- ggplot(ece_long,
                aes(stage, ece, group = interaction(seed, dab),
                    color = dab)) +
  geom_line(alpha = 0.6) +
  geom_point(size = 2, alpha = 0.85) +
  scale_color_brewer(palette = "Dark2", name = "dAb") +
  labs(x = NULL, y = "Expected Calibration Error (ECE, lower is better)",
       title = "ECE before vs after per-class temperature scaling",
       subtitle = "One line per (seed, dAb); 5 seeds x 8 dAbs = 40 lines") +
  theme_bw() +
  theme(panel.grid.minor = element_blank())
ggsave("Thesis_figures/figs/11c_temp_scale_ece_before_after.png", p_ece,
       width = 9, height = 5, dpi = 300)

cat("\nDone. Saved:\n",
    "  Embedded_FS/rf_temp_scale_results.rds\n",
    "  Embedded_FS/rf_temp_scale_summary.csv\n",
    "  Thesis_figures/figs/11b_temp_scale_T_k.png\n",
    "  Thesis_figures/figs/11c_temp_scale_ece_before_after.png\n", sep = "")
