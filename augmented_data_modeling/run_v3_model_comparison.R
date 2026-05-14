suppressPackageStartupMessages({
  library(tidymodels)
  library(themis)
  library(ranger)
  library(xgboost)
  library(glmnet)
  library(kernlab)
  library(ggplot2)
  library(tidyr)
})

set.seed(2024)

setwd("augmented_data_modeling")

out_dir = "v3_outputs_model_comparison"
dir.create(out_dir, showWarnings = FALSE)

# ============================================================================
# Compares RF, XGBoost, multinomial elastic net, and SVM-RBF on the combined
# augmented data. All four models share:
#   - same Conc. threshold (0.9 ng/ml, applied pre-split)
#   - same train/test split, same CV folds
#   - same downstream evaluation (test-set metrics, ROC, per-conc accuracy)
# Recipes differ only where the model genuinely needs it (RF doesn't normalize;
# the others do).
# ============================================================================

# ----- Data + split -----
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

cat("rows after threshold:", nrow(clean), "\n")

split    = initial_split(clean, prop = 0.8, strata = strata)
train_df = training(split) |> select(-strata)
test_df  = testing(split)  |> select(-strata)

cab_cols = grep("^cAb", names(train_df), value = TRUE)
classes  = levels(train_df$target_vector)
prob_col_names = paste0(".pred_", classes)

# ----- Two recipes: tree-based vs scale-sensitive -----
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

rec_tree   = make_rec_tree(train_df)
rec_scaled = make_rec_scaled(train_df)

# ----- Folds (shared across models) -----
folds = vfold_cv(train_df, v = 5, strata = target_vector)

p_pred = ncol(train_df) - 3L  # subtract target_vector, Conc., batch
mtry_values = sort(unique(pmax(1L, c(2L, floor(sqrt(p_pred)),
                                     floor(p_pred / 3L), floor(p_pred / 2L)))))

cv_metrics = metric_set(accuracy, kap, f_meas)

# ============================================================================
# Model specs + grids
# ============================================================================
spec_rf = rand_forest(mtry = tune(), min_n = tune(), trees = 500) |>
  set_engine("ranger", importance = "permutation", seed = 2024) |>
  set_mode("classification")
grid_rf = expand.grid(mtry = mtry_values, min_n = c(2, 5, 10))

spec_xgb = boost_tree(
  trees      = 500,
  tree_depth = tune(),
  learn_rate = tune(),
  mtry       = tune(),
  min_n      = tune()
) |>
  set_engine("xgboost", verbosity = 0) |>
  set_mode("classification")
grid_xgb = expand.grid(
  tree_depth = c(2, 4, 6),
  learn_rate = c(0.05, 0.1),
  mtry       = mtry_values,
  min_n      = c(5, 10)
)

spec_glm = multinom_reg(penalty = tune(), mixture = tune()) |>
  set_engine("glmnet") |>
  set_mode("classification")
grid_glm = expand.grid(
  penalty = 10 ^ seq(-4, 0, length.out = 8),
  mixture = c(0, 0.25, 0.5, 0.75, 1)   # 0 = ridge, 1 = LASSO
)

spec_svm = svm_rbf(cost = tune(), rbf_sigma = tune()) |>
  set_engine("kernlab") |>
  set_mode("classification")
grid_svm = expand.grid(
  cost      = c(0.5, 1, 2, 4, 8),
  rbf_sigma = c(0.01, 0.05, 0.1, 0.25)
)

models = list(
  RF      = list(spec = spec_rf,  rec = rec_tree,   grid = grid_rf,  metric_tune = "accuracy"),
  XGBoost = list(spec = spec_xgb, rec = rec_tree,   grid = grid_xgb, metric_tune = "accuracy"),
  ElasNet = list(spec = spec_glm, rec = rec_scaled, grid = grid_glm, metric_tune = "accuracy"),
  SVM_RBF = list(spec = spec_svm, rec = rec_scaled, grid = grid_svm, metric_tune = "accuracy")
)

# ----- Helpers -----
class_metrics = metric_set(accuracy, kap, f_meas, sens, spec)
prob_metrics  = metric_set(roc_auc, mn_log_loss)

eval_test = function(final_model, test_df, model_name) {
  aug = augment(final_model, new_data = test_df)
  pcols = grep("^\\.pred_", names(aug), value = TRUE)
  pcols = setdiff(pcols, ".pred_class")

  cls_res = aug |> class_metrics(truth = target_vector, estimate = .pred_class)
  prb_res = aug |> prob_metrics(truth = target_vector, !!!syms(pcols))
  metrics = bind_rows(cls_res, prb_res) |> mutate(model = model_name)

  preds = aug |>
    select(target_vector, .pred_class, all_of(pcols), Conc., batch) |>
    mutate(model = model_name)

  list(metrics = metrics, preds = preds, fit = final_model)
}

# ============================================================================
# Tune + finalize each model
# ============================================================================
results = list()

for (name in names(models)) {
  cfg = models[[name]]
  cat("\n==========================================================\n")
  cat("Tuning ", name, "  (", nrow(cfg$grid), " configs x 5 folds)\n", sep = "")
  cat("==========================================================\n")

  wf = workflow() |>
    add_recipe(cfg$rec) |>
    add_model(cfg$spec)

  t0 = Sys.time()
  tune_res = tune_grid(
    wf,
    resamples = folds,
    grid      = cfg$grid,
    metrics   = cv_metrics,
    control   = control_grid(save_pred = FALSE, verbose = FALSE)
  )
  cat("  tune time:", round(as.numeric(difftime(Sys.time(), t0, units = "secs")), 1), "s\n")

  best     = select_best(tune_res, metric = cfg$metric_tune)
  cat("  best:\n"); print(best)

  final_wf    = finalize_workflow(wf, best)
  final_model = fit(final_wf, data = train_df)

  res = eval_test(final_model, test_df, name)
  res$best     = best
  res$tune_res = tune_res
  results[[name]] = res

  cat("\n  test metrics:\n")
  print(res$metrics)
}

# ============================================================================
# Combined tables + plots
# ============================================================================
all_metrics = bind_rows(lapply(results, `[[`, "metrics"))
all_preds   = bind_rows(lapply(results, `[[`, "preds"))

cat("\n\n===== Combined test metrics =====\n")
print(all_metrics |> pivot_wider(names_from = .metric, values_from = .estimate))

write.csv(all_metrics, file.path(out_dir, "metrics_long.csv"), row.names = FALSE)
write.csv(
  all_metrics |> select(model, .metric, .estimate) |>
    pivot_wider(names_from = .metric, values_from = .estimate),
  file.path(out_dir, "metrics_wide.csv"), row.names = FALSE
)

# Color palette consistent across plots
model_colors = c(
  RF      = "#266B6E",
  XGBoost = "#7BA7AB",
  ElasNet = "#B5403F",
  SVM_RBF = "#1F4E50"
)
model_order = names(model_colors)

# ----- Metric bar chart (faceted by metric) -----
keep_metrics = c("accuracy", "kap", "f_meas", "roc_auc", "mn_log_loss")
p_metrics = all_metrics |>
  filter(.metric %in% keep_metrics) |>
  mutate(model = factor(model, levels = model_order),
         .metric = factor(.metric, levels = keep_metrics)) |>
  ggplot(aes(model, .estimate, fill = model)) +
  geom_col() +
  geom_text(aes(label = round(.estimate, 3)), vjust = -0.3, size = 3) +
  facet_wrap(~ .metric, scales = "free_y") +
  scale_fill_manual(values = model_colors) +
  labs(title    = "Model comparison — combined augmented data",
       subtitle = paste0("Held-out test set | Conc. > ", conc_threshold, " ng/ml"),
       x = NULL, y = NULL) +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "metrics.png"), p_metrics,
       width = 11, height = 7, dpi = 150)

# ----- ROC curves: facet by class, one line per model -----
roc_df = bind_rows(lapply(names(results), function(name) {
  preds = results[[name]]$preds
  pcols = grep("^\\.pred_", names(preds), value = TRUE)
  pcols = setdiff(pcols, ".pred_class")
  preds |>
    roc_curve(truth = target_vector, !!!syms(pcols)) |>
    mutate(model = name)
}))

p_roc = roc_df |>
  mutate(model = factor(model, levels = model_order)) |>
  ggplot(aes(1 - specificity, sensitivity, color = model)) +
  geom_path(linewidth = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray60") +
  facet_wrap(~ .level) +
  coord_equal() +
  scale_color_manual(values = model_colors) +
  labs(title    = "ROC curves by model and variant",
       subtitle = paste0("Held-out test set | Conc. > ", conc_threshold, " ng/ml"),
       x = "1 - Specificity", y = "Sensitivity",
       color = NULL) +
  theme_minimal()
ggsave(file.path(out_dir, "roc_curves.png"), p_roc,
       width = 10, height = 8, dpi = 150)

# ----- Per-concentration accuracy by model -----
per_conc = all_preds |>
  mutate(Conc. = round(Conc., 4),
         correct = as.integer(.pred_class == target_vector)) |>
  group_by(model, Conc.) |>
  summarise(num_correct     = sum(correct),
            total           = n(),
            percent_correct = round(num_correct / total * 100, 2),
            .groups         = "drop") |>
  arrange(model, Conc.)

write.csv(per_conc, file.path(out_dir, "per_conc.csv"), row.names = FALSE)

p_per_conc = per_conc |>
  mutate(model = factor(model, levels = model_order),
         Conc. = factor(Conc.)) |>
  ggplot(aes(Conc., percent_correct, fill = model)) +
  geom_col() +
  geom_text(aes(label = percent_correct), vjust = 1.5, color = "white", size = 2.7) +
  facet_wrap(~ model, ncol = 1) +
  scale_fill_manual(values = model_colors) +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001),
                     limits = c(0, 110)) +
  labs(title    = "Accuracy by concentration — model comparison",
       subtitle = paste0("Held-out test set | Conc. > ", conc_threshold, " ng/ml"),
       x        = "Concentration (ng/ml)",
       y        = "Percent Correct") +
  theme_classic() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "percent_correct_by_conc.png"), p_per_conc,
       width = 9, height = 11, dpi = 150)

# ----- Confusion matrices (4-panel) -----
cm_df = bind_rows(lapply(names(results), function(name) {
  preds = results[[name]]$preds
  preds |>
    count(target_vector, .pred_class, name = "n") |>
    complete(target_vector, .pred_class, fill = list(n = 0)) |>
    mutate(model = name)
}))

p_cm = cm_df |>
  mutate(model = factor(model, levels = model_order)) |>
  ggplot(aes(.pred_class, target_vector, fill = n)) +
  geom_tile() +
  geom_text(aes(label = n), color = "white", size = 4) +
  facet_wrap(~ model, ncol = 2) +
  scale_fill_gradient(low = "#7BA7AB", high = "#0E2E2F") +
  labs(title    = "Confusion matrices — held-out test set",
       subtitle = paste0("Conc. > ", conc_threshold, " ng/ml"),
       x        = "Predicted", y = "Truth") +
  theme_minimal() +
  theme(legend.position = "none")
ggsave(file.path(out_dir, "confusion_matrices.png"), p_cm,
       width = 9, height = 8, dpi = 150)

# ----- Best hyperparameters per model -----
best_df = bind_rows(lapply(names(results), function(name) {
  results[[name]]$best |> mutate(model = name) |>
    pivot_longer(-c(.config, model), names_to = "param", values_to = "value") |>
    mutate(value = as.character(value))
}))
write.csv(best_df, file.path(out_dir, "best_hyperparams.csv"), row.names = FALSE)

cat("\nWrote outputs to ", out_dir, "/\n", sep = "")
