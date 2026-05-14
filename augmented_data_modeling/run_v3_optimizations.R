suppressPackageStartupMessages({
  library(tidymodels)
  library(themis)
  library(ranger)
  library(ggplot2)
  library(tidyr)
  library(dplyr)
  library(purrr)
})

setwd("augmented_data_modeling")

out_dir = "v3_optimizations"
dir.create(out_dir, showWarnings = FALSE)

SEEDS = c(1969, 1986, 2024, 2025, 2026)

# ============================================================================
# Factorial sweep over outlier-removal x synthetic-augmentation methods, with
# RF held fixed as the workhorse (winner of the model-family comparison in
# run_v3_model_comparison_multiseed.R). Every cell is retuned per seed on a
# small mtry x min_n grid so augmentation can shift the optimum.
#
# Outlier removal is fit on TRAIN ONLY (never the test set). Augmentation is
# a recipe step with skip = TRUE so it does not leak into evaluation.
# ============================================================================

# ----- Data (seed-independent) -----
raw = read.csv("combined_augment_with_batch.csv") |>
  as_tibble() |>
  mutate(target_vector = factor(target_vector),
         batch         = factor(batch)) |>
  select(-cAb13)

conc_threshold = 0.9

clean = raw |>
  filter(Conc. > conc_threshold) |>
  select(-Replicate) |>
  mutate(strata = interaction(target_vector, batch, drop = TRUE))

cab_cols = grep("^cAb", names(clean), value = TRUE)

cat("rows after threshold:", nrow(clean), "\n")
cat("seeds:", paste(SEEDS, collapse = ", "), "\n\n")

# ============================================================================
# Outlier removal helpers (fit on training set only)
# ============================================================================

# Per-class, per-feature 1.5xIQR rule. Flag a row if any cab feature falls
# outside [Q1 - 1.5*IQR, Q3 + 1.5*IQR] within its class.
flag_iqr = function(df, feat_cols, k = 1.5) {
  df_idx = df |> mutate(.row = row_number())
  flagged = df_idx |>
    pivot_longer(all_of(feat_cols), names_to = "feat", values_to = "val") |>
    group_by(target_vector, feat) |>
    mutate(q1  = quantile(val, 0.25, na.rm = TRUE),
           q3  = quantile(val, 0.75, na.rm = TRUE),
           iqr = q3 - q1,
           lo  = q1 - k * iqr,
           hi  = q3 + k * iqr,
           bad = val < lo | val > hi) |>
    ungroup() |>
    group_by(.row) |>
    summarise(any_bad = any(bad, na.rm = TRUE), .groups = "drop")
  flagged$.row[flagged$any_bad]
}

# Per-class Mahalanobis distance with chi-square 0.975 cutoff. Falls back to
# generalized inverse if the class covariance is singular.
flag_mahalanobis = function(df, feat_cols, alpha = 0.975) {
  X = as.matrix(df[, feat_cols])
  cls = df$target_vector
  cutoff = qchisq(alpha, df = length(feat_cols))
  bad = logical(nrow(df))
  for (lv in levels(cls)) {
    idx = which(cls == lv)
    if (length(idx) <= length(feat_cols) + 1) next
    Xc = X[idx, , drop = FALSE]
    mu = colMeans(Xc)
    S  = cov(Xc)
    Sinv = tryCatch(solve(S), error = function(e) MASS::ginv(S))
    diffs = sweep(Xc, 2, mu, "-")
    d2 = rowSums((diffs %*% Sinv) * diffs)
    bad[idx] = d2 > cutoff
  }
  which(bad)
}

apply_outlier = function(df, method, feat_cols) {
  drops = switch(method,
    none        = integer(0),
    iqr         = flag_iqr(df, feat_cols),
    mahalanobis = flag_mahalanobis(df, feat_cols),
    stop("unknown outlier method: ", method)
  )
  list(
    df      = if (length(drops)) df[-drops, ] else df,
    n_drop  = length(drops),
    n_keep  = nrow(df) - length(drops)
  )
}

# ============================================================================
# Recipe builders -- augmentation as the only varying step
# ============================================================================

base_recipe = function(d) {
  recipe(target_vector ~ ., data = d) |>
    step_rm(Conc., batch) |>
    step_zv(all_numeric_predictors()) |>
    step_impute_knn(all_numeric_predictors(), neighbors = 3)
}

make_recipe = function(d, aug) {
  rec = base_recipe(d)
  switch(aug,
    none     = rec,
    smote    = rec |> step_smote(target_vector, over_ratio = 1, skip = TRUE),
    bsmote   = rec |> step_bsmote(target_vector, over_ratio = 1, skip = TRUE),
    rose     = rec |> step_rose(target_vector, skip = TRUE),
    upsample = rec |> step_upsample(target_vector, over_ratio = 1, skip = TRUE),
    stop("unknown augmentation: ", aug)
  )
}

# ============================================================================
# Spec / grid / metrics
# ============================================================================

p_pred = length(cab_cols)
mtry_values = sort(unique(pmax(1L, c(2L, floor(sqrt(p_pred)), floor(p_pred / 3L)))))

rf_grid = expand.grid(mtry = mtry_values, min_n = c(2, 5, 10))

build_rf_spec = function(seed) {
  rand_forest(mtry = tune(), min_n = tune(), trees = 500) |>
    set_engine("ranger", importance = "permutation", seed = seed) |>
    set_mode("classification")
}

class_metrics = metric_set(accuracy, kap, f_meas, sens, spec)
prob_metrics  = metric_set(roc_auc, mn_log_loss)
cv_metrics    = metric_set(accuracy, kap, f_meas)

OUTLIER_METHODS = c("none", "iqr", "mahalanobis")
AUG_METHODS     = c("none", "smote", "bsmote", "rose", "upsample")

# ============================================================================
# Per-seed runner
# ============================================================================

run_one_seed = function(seed) {
  cat("\n##############################################\n")
  cat("# SEED ", seed, "\n", sep = "")
  cat("##############################################\n")
  set.seed(seed)

  split    = initial_split(clean, prop = 0.8, strata = strata)
  train_df_full = training(split) |> select(-strata)
  test_df       = testing(split)  |> select(-strata)

  per_seed = list()
  for (om in OUTLIER_METHODS) {
    or = apply_outlier(train_df_full, om, cab_cols)
    train_df = or$df
    cat("  outlier=", om, " | dropped ", or$n_drop,
        " of ", nrow(train_df_full), " train rows\n", sep = "")

    folds = vfold_cv(train_df, v = 5, strata = target_vector)
    spec  = build_rf_spec(seed)

    for (aug in AUG_METHODS) {
      cell_id = paste(om, aug, sep = "/")
      cat("    aug=", aug, "\n", sep = "")

      rec = make_recipe(train_df, aug)
      wf  = workflow() |> add_recipe(rec) |> add_model(spec)

      tune_res = tune_grid(
        wf, resamples = folds, grid = rf_grid,
        metrics = cv_metrics,
        control = control_grid(save_pred = FALSE, verbose = FALSE)
      )

      best        = select_best(tune_res, metric = "accuracy")
      final_wf    = finalize_workflow(wf, best)
      final_model = fit(final_wf, data = train_df)

      aug_pred = augment(final_model, new_data = test_df)
      pcols = setdiff(grep("^\\.pred_", names(aug_pred), value = TRUE), ".pred_class")

      cls_res = aug_pred |> class_metrics(truth = target_vector, estimate = .pred_class)
      prb_res = aug_pred |> prob_metrics(truth = target_vector, !!!syms(pcols))
      metrics = bind_rows(cls_res, prb_res) |>
        mutate(outlier = om, aug = aug, seed = seed, cell = cell_id,
               n_train_kept = or$n_keep, n_train_drop = or$n_drop)

      preds = aug_pred |>
        select(target_vector, .pred_class, all_of(pcols), Conc., batch) |>
        mutate(outlier = om, aug = aug, seed = seed, cell = cell_id)

      best_long = best |>
        pivot_longer(-.config, names_to = "param", values_to = "value") |>
        mutate(value = as.character(value),
               outlier = om, aug = aug, seed = seed, cell = cell_id)

      drop_row = tibble(outlier = om, aug = aug, seed = seed, cell = cell_id,
                        n_train_kept = or$n_keep, n_train_drop = or$n_drop,
                        n_train_full = nrow(train_df_full))

      per_seed[[cell_id]] = list(metrics = metrics, preds = preds,
                                 best = best_long, drops = drop_row)
    }
  }
  per_seed
}

# ============================================================================
# Run all seeds
# ============================================================================
all_runs = lapply(SEEDS, run_one_seed)

flat = function(field) {
  bind_rows(unlist(lapply(all_runs, function(seed_runs) {
    lapply(seed_runs, `[[`, field)
  }), recursive = FALSE))
}

metrics_all = flat("metrics")
preds_all   = flat("preds")
best_all    = flat("best")
drops_all   = flat("drops") |> distinct()

write.csv(metrics_all, file.path(out_dir, "metrics_per_seed.csv"), row.names = FALSE)
write.csv(best_all,    file.path(out_dir, "best_hyperparams_per_seed.csv"), row.names = FALSE)
write.csv(drops_all,   file.path(out_dir, "outlier_drop_counts.csv"), row.names = FALSE)

# ----- Aggregate metrics: mean and SE across seeds, per (outlier, aug) cell ----
metrics_summary = metrics_all |>
  group_by(outlier, aug, .metric) |>
  summarise(
    n_seeds = n(),
    mean    = mean(.estimate),
    se      = sd(.estimate) / sqrt(n()),
    .groups = "drop"
  )

cat("\n\n===== Aggregated metrics across", length(SEEDS), "seeds =====\n")
print(metrics_summary |>
        filter(.metric %in% c("accuracy", "kap", "f_meas", "roc_auc", "mn_log_loss")) |>
        arrange(.metric, desc(mean)),
      n = Inf)

write.csv(metrics_summary, file.path(out_dir, "metrics_summary.csv"), row.names = FALSE)

# Wide form table (one row per (outlier,aug), one column per metric)
metrics_wide = metrics_summary |>
  filter(.metric %in% c("accuracy", "kap", "f_meas", "roc_auc", "mn_log_loss")) |>
  mutate(label = paste0(round(mean, 3), " ± ", round(se, 3))) |>
  select(outlier, aug, .metric, label) |>
  pivot_wider(names_from = .metric, values_from = label)
write.csv(metrics_wide, file.path(out_dir, "metrics_summary_wide.csv"), row.names = FALSE)
cat("\n--- Wide table (mean ± SE) ---\n")
print(metrics_wide, n = Inf)

# ============================================================================
# Plots
# ============================================================================

outlier_levels = OUTLIER_METHODS
aug_levels     = AUG_METHODS

# ----- Heatmap of mean accuracy / kap / roc_auc / log_loss -----
heat_df = metrics_summary |>
  filter(.metric %in% c("accuracy", "kap", "roc_auc", "mn_log_loss")) |>
  mutate(outlier = factor(outlier, levels = outlier_levels),
         aug     = factor(aug,     levels = aug_levels))

p_heat = ggplot(heat_df, aes(aug, outlier, fill = mean)) +
  geom_tile(color = "white") +
  geom_text(aes(label = sprintf("%.3f", mean)), size = 3) +
  facet_wrap(~ .metric, scales = "free") +
  scale_fill_viridis_c(option = "D") +
  labs(title    = "Outlier x Augmentation factorial (RF, 5 seeds)",
       subtitle = paste0("Mean test metric across seeds | Conc. > ",
                         conc_threshold, " ng/ml"),
       x = "Augmentation", y = "Outlier removal") +
  theme_minimal()
ggsave(file.path(out_dir, "metrics_heatmap.png"), p_heat,
       width = 11, height = 8, dpi = 150)

# ----- Per-seed dotplot -----
p_dot = metrics_all |>
  filter(.metric %in% c("accuracy", "kap", "roc_auc", "mn_log_loss")) |>
  mutate(outlier = factor(outlier, levels = outlier_levels),
         aug     = factor(aug,     levels = aug_levels),
         cell    = paste(outlier, aug, sep = "/")) |>
  ggplot(aes(.estimate, cell, color = outlier, shape = aug)) +
  geom_point(size = 2.6, alpha = 0.85) +
  facet_wrap(~ .metric, scales = "free_x") +
  labs(title = "Per-seed metric scatter across outlier x aug cells",
       x = "Test metric", y = NULL) +
  theme_minimal() +
  theme(legend.position = "right")
ggsave(file.path(out_dir, "per_seed_scatter.png"), p_dot,
       width = 12, height = 8, dpi = 150)

# ----- Mean accuracy bar chart with SE error bars (cell-level) -----
keep_metrics = c("accuracy", "kap", "f_meas", "roc_auc", "mn_log_loss")
p_bars = metrics_summary |>
  filter(.metric %in% keep_metrics) |>
  mutate(outlier = factor(outlier, levels = outlier_levels),
         aug     = factor(aug,     levels = aug_levels)) |>
  ggplot(aes(aug, mean, fill = outlier)) +
  geom_col(position = position_dodge(width = 0.85)) +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se),
                position = position_dodge(width = 0.85),
                width = 0.3) +
  facet_wrap(~ .metric, scales = "free_y") +
  labs(title    = "RF performance by outlier x augmentation",
       subtitle = paste0("Mean ± SE across ", length(SEEDS),
                         " seeds | Conc. > ", conc_threshold, " ng/ml"),
       x = "Augmentation", y = NULL, fill = "Outlier") +
  theme_minimal()
ggsave(file.path(out_dir, "metrics_bars.png"), p_bars,
       width = 12, height = 8, dpi = 150)

# ----- Outlier drop counts (sanity panel) -----
p_drops = drops_all |>
  filter(outlier != "none") |>
  ggplot(aes(factor(seed), n_train_drop, fill = outlier)) +
  geom_col(position = position_dodge(width = 0.8)) +
  geom_text(aes(label = n_train_drop),
            position = position_dodge(width = 0.8),
            vjust = -0.3, size = 3) +
  labs(title = "Training rows flagged as outliers per seed",
       x = "Seed", y = "Rows dropped", fill = "Method") +
  theme_minimal()
ggsave(file.path(out_dir, "outlier_drop_counts.png"), p_drops,
       width = 9, height = 5, dpi = 150)

cat("\nWrote outputs to ", out_dir, "/\n", sep = "")
