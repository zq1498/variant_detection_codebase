suppressPackageStartupMessages({
  library(tidymodels)
  library(themis)
  library(ranger)
  library(xgboost)
  library(glmnet)
  library(kernlab)
  library(ggplot2)
  library(tidyr)
  library(purrr)
})

setwd("augmented_data_modeling")

out_dir = "v3_outputs_model_comparison_multiseed"
dir.create(out_dir, showWarnings = FALSE)

SEEDS = c(1969, 1986, 2024, 2025, 2026)

# ============================================================================
# Multi-seed version of run_v3_model_comparison.R. Each seed gets a fresh
# initial_split + vfold_cv so we capture train/test variance, not just CV
# variance. Tuning grid, recipes, and specs are identical across seeds.
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
classes  = levels(clean$target_vector)

cat("rows after threshold:", nrow(clean), "\n")
cat("seeds:", paste(SEEDS, collapse = ", "), "\n\n")

# ----- Recipes -----
make_rec_tree = function(d) {
  recipe(target_vector ~ ., data = d) |>
    step_rm(Conc., batch) |>
    step_zv(all_numeric_predictors()) |>
    step_impute_knn(all_numeric_predictors(), neighbors = 3) |>
    step_smote(target_vector, over_ratio = 1, skip = TRUE)
}
make_rec_scaled = function(d) {
  recipe(target_vector ~ ., data = d) |>
    step_rm(Conc., batch) |>
    step_zv(all_numeric_predictors()) |>
    step_impute_knn(all_numeric_predictors(), neighbors = 3) |>
    step_normalize(all_numeric_predictors()) |>
    step_smote(target_vector, over_ratio = 1, skip = TRUE)
}

# ----- Grids -----
p_pred = ncol(clean) - 4L  # subtract target_vector, Conc., batch, strata
mtry_values = sort(unique(pmax(1L, c(2L, floor(sqrt(p_pred)),
                                     floor(p_pred / 3L), floor(p_pred / 2L)))))

build_models = function(seed) {
  list(
    RF = list(
      spec = rand_forest(mtry = tune(), min_n = tune(), trees = 500) |>
        set_engine("ranger", importance = "permutation", seed = seed) |>
        set_mode("classification"),
      grid = expand.grid(mtry = mtry_values, min_n = c(2, 5, 10)),
      rec  = "tree"
    ),
    XGBoost = list(
      spec = boost_tree(trees = 500, tree_depth = tune(), learn_rate = tune(),
                        mtry = tune(), min_n = tune()) |>
        set_engine("xgboost", verbosity = 0) |>
        set_mode("classification"),
      grid = expand.grid(tree_depth = c(2, 4, 6),
                         learn_rate = c(0.05, 0.1),
                         mtry       = mtry_values,
                         min_n      = c(5, 10)),
      rec  = "tree"
    ),
    ElasNet = list(
      spec = multinom_reg(penalty = tune(), mixture = tune()) |>
        set_engine("glmnet") |>
        set_mode("classification"),
      grid = expand.grid(penalty = 10 ^ seq(-4, 0, length.out = 8),
                         mixture = c(0, 0.25, 0.5, 0.75, 1)),
      rec  = "scaled"
    ),
    SVM_RBF = list(
      spec = svm_rbf(cost = tune(), rbf_sigma = tune()) |>
        set_engine("kernlab") |>
        set_mode("classification"),
      grid = expand.grid(cost = c(0.5, 1, 2, 4, 8),
                         rbf_sigma = c(0.01, 0.05, 0.1, 0.25)),
      rec  = "scaled"
    )
  )
}

class_metrics = metric_set(accuracy, kap, f_meas, sens, spec)
prob_metrics  = metric_set(roc_auc, mn_log_loss)
cv_metrics    = metric_set(accuracy, kap, f_meas)

# ----- Per-seed runner -----
run_one_seed = function(seed) {
  cat("\n##############################################\n")
  cat("# SEED ", seed, "\n", sep = "")
  cat("##############################################\n")
  set.seed(seed)

  split    = initial_split(clean, prop = 0.8, strata = strata)
  train_df = training(split) |> select(-strata)
  test_df  = testing(split)  |> select(-strata)

  rec_tree   = make_rec_tree(train_df)
  rec_scaled = make_rec_scaled(train_df)
  rec_lookup = list(tree = rec_tree, scaled = rec_scaled)

  folds = vfold_cv(train_df, v = 5, strata = target_vector)
  models = build_models(seed)

  per_seed = list()
  for (name in names(models)) {
    cfg = models[[name]]
    cat("  -- ", name, " (", nrow(cfg$grid), " configs)\n", sep = "")

    wf = workflow() |>
      add_recipe(rec_lookup[[cfg$rec]]) |>
      add_model(cfg$spec)

    tune_res = tune_grid(
      wf, resamples = folds, grid = cfg$grid,
      metrics = cv_metrics,
      control = control_grid(save_pred = FALSE, verbose = FALSE)
    )

    best        = select_best(tune_res, metric = "accuracy")
    final_wf    = finalize_workflow(wf, best)
    final_model = fit(final_wf, data = train_df)

    aug = augment(final_model, new_data = test_df)
    pcols = setdiff(grep("^\\.pred_", names(aug), value = TRUE), ".pred_class")

    cls_res = aug |> class_metrics(truth = target_vector, estimate = .pred_class)
    prb_res = aug |> prob_metrics(truth = target_vector, !!!syms(pcols))
    metrics = bind_rows(cls_res, prb_res) |>
      mutate(model = name, seed = seed)

    preds = aug |>
      select(target_vector, .pred_class, all_of(pcols), Conc., batch) |>
      mutate(model = name, seed = seed)

    roc = preds |>
      roc_curve(truth = target_vector, !!!syms(pcols)) |>
      mutate(model = name, seed = seed)

    best_long = best |>
      pivot_longer(-.config, names_to = "param", values_to = "value") |>
      mutate(value = as.character(value), model = name, seed = seed)

    per_seed[[name]] = list(metrics = metrics, preds = preds,
                            roc = roc, best = best_long)
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
roc_all     = flat("roc")
best_all    = flat("best")

write.csv(metrics_all, file.path(out_dir, "metrics_per_seed.csv"), row.names = FALSE)
write.csv(best_all,    file.path(out_dir, "best_hyperparams_per_seed.csv"), row.names = FALSE)

# ----- Aggregate metrics: mean and SE across seeds -----
metrics_summary = metrics_all |>
  group_by(model, .metric) |>
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

# Wide form for the README-style table
metrics_wide = metrics_summary |>
  filter(.metric %in% c("accuracy", "kap", "f_meas", "roc_auc", "mn_log_loss")) |>
  mutate(label = paste0(round(mean, 3), " ± ", round(se, 3))) |>
  select(model, .metric, label) |>
  pivot_wider(names_from = .metric, values_from = label)
write.csv(metrics_wide, file.path(out_dir, "metrics_summary_wide.csv"), row.names = FALSE)
cat("\n--- Wide table (mean ± SE) ---\n")
print(metrics_wide)

# ----- Plot palette -----
model_colors = c(RF = "#266B6E", XGBoost = "#7BA7AB",
                 ElasNet = "#B5403F", SVM_RBF = "#1F4E50")
model_order = names(model_colors)

# ----- Aggregate metric bar chart with SE error bars -----
keep_metrics = c("accuracy", "kap", "f_meas", "roc_auc", "mn_log_loss")
p_metrics = metrics_summary |>
  filter(.metric %in% keep_metrics) |>
  mutate(model = factor(model, levels = model_order),
         .metric = factor(.metric, levels = keep_metrics)) |>
  ggplot(aes(model, mean, fill = model)) +
  geom_col() +
  geom_errorbar(aes(ymin = mean - se, ymax = mean + se), width = 0.3) +
  geom_text(aes(label = round(mean, 3)), vjust = -1.6, size = 3) +
  facet_wrap(~ .metric, scales = "free_y") +
  scale_fill_manual(values = model_colors) +
  labs(title    = "Model comparison — combined augmented data",
       subtitle = paste0("Mean ± SE across ", length(SEEDS),
                         " seeds | Conc. > ", conc_threshold, " ng/ml"),
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "metrics.png"), p_metrics,
       width = 11, height = 7, dpi = 150)

# ----- ROC curves: one panel per class, color by model, faint lines per seed,
#       plus a per-model average step-curve over a common FPR grid.
fpr_grid = seq(0, 1, by = 0.01)
roc_avg = roc_all |>
  group_by(model, .level, seed) |>
  arrange(specificity) |>
  group_modify(~ tibble(
    fpr = fpr_grid,
    tpr = approx(x    = 1 - .x$specificity,
                 y    = .x$sensitivity,
                 xout = fpr_grid,
                 ties = "ordered",
                 rule = 2)$y
  )) |>
  ungroup() |>
  group_by(model, .level, fpr) |>
  summarise(tpr = mean(tpr), .groups = "drop")

p_roc = ggplot() +
  geom_path(data = roc_all |> mutate(model = factor(model, levels = model_order)),
            aes(1 - specificity, sensitivity,
                color = model, group = interaction(model, seed)),
            linewidth = 0.4, alpha = 0.35) +
  geom_path(data = roc_avg |> mutate(model = factor(model, levels = model_order)),
            aes(fpr, tpr, color = model),
            linewidth = 1) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  facet_wrap(~ .level) +
  coord_equal() +
  scale_color_manual(values = model_colors) +
  labs(title    = "ROC curves by model and variant (5 seeds)",
       subtitle = paste0("Faint lines = per-seed; bold = mean | Conc. > ",
                         conc_threshold, " ng/ml"),
       x = "1 - Specificity", y = "Sensitivity",
       color = NULL) +
  theme_minimal()
ggsave(file.path(out_dir, "roc_curves.png"), p_roc,
       width = 10, height = 8, dpi = 150)

# ----- Per-concentration accuracy across seeds -----
per_conc_seed = preds_all |>
  mutate(Conc. = round(Conc., 4),
         correct = as.integer(.pred_class == target_vector)) |>
  group_by(model, seed, Conc.) |>
  summarise(num_correct     = sum(correct),
            total           = n(),
            percent_correct = num_correct / total * 100,
            .groups         = "drop")

per_conc_summary = per_conc_seed |>
  group_by(model, Conc.) |>
  summarise(
    n_seeds   = n(),
    mean_pct  = mean(percent_correct),
    se_pct    = sd(percent_correct) / sqrt(n()),
    .groups = "drop"
  )
write.csv(per_conc_summary, file.path(out_dir, "per_conc_summary.csv"), row.names = FALSE)
write.csv(per_conc_seed,    file.path(out_dir, "per_conc_per_seed.csv"), row.names = FALSE)

p_per_conc = per_conc_summary |>
  mutate(model = factor(model, levels = model_order),
         Conc. = factor(Conc.)) |>
  ggplot(aes(Conc., mean_pct, fill = model)) +
  geom_col() +
  geom_errorbar(aes(ymin = pmax(mean_pct - se_pct, 0),
                    ymax = pmin(mean_pct + se_pct, 100)),
                width = 0.3) +
  geom_text(aes(label = round(mean_pct, 1)), vjust = 1.5, color = "white", size = 2.7) +
  facet_wrap(~ model, ncol = 1) +
  scale_fill_manual(values = model_colors) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.1),
                     limits = c(0, 110)) +
  labs(title    = "Accuracy by concentration — model comparison (5 seeds)",
       subtitle = paste0("Mean ± SE across seeds | Conc. > ",
                         conc_threshold, " ng/ml"),
       x        = "Concentration (ng/ml)",
       y        = "Percent Correct") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "percent_correct_by_conc.png"), p_per_conc,
       width = 9, height = 11, dpi = 150)

# ----- Per-seed dotplot of accuracy / kap (rank stability check) -----
p_dot = metrics_all |>
  filter(.metric %in% c("accuracy", "kap", "roc_auc", "mn_log_loss")) |>
  mutate(model = factor(model, levels = model_order)) |>
  ggplot(aes(.estimate, model, color = model)) +
  geom_point(size = 3, alpha = 0.8) +
  facet_wrap(~ .metric, scales = "free_x") +
  scale_color_manual(values = model_colors) +
  labs(title = "Per-seed metric scatter (rank stability)",
       x = "Test metric", y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "per_seed_scatter.png"), p_dot,
       width = 10, height = 6, dpi = 150)

cat("\nWrote outputs to ", out_dir, "/\n", sep = "")
