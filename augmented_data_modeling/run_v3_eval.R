suppressPackageStartupMessages({
  library(tidymodels)
  library(themis)
  library(ranger)
  library(ggplot2)
})

set.seed(2024)

setwd("augmented_data_modeling")

raw = read.csv("combined_augment_with_batch.csv") |>
  as_tibble() |>
  mutate(
    target_vector = factor(target_vector),
    batch         = factor(batch)
  ) |>
  select(-cAb13)

cat("rows:", nrow(raw), "  cols:", ncol(raw), "\n")
cat("class counts:\n"); print(table(raw$target_vector))
cat("batch counts:\n"); print(table(raw$batch))

conc_threshold = 0.9

clean = raw |>
  filter(Conc. > conc_threshold) |>
  select(-Replicate) |>
  mutate(strata = interaction(target_vector, batch, drop = TRUE))

cat("after threshold:", nrow(clean), "rows\n")

split    = initial_split(clean, prop = 0.8, strata = strata)
train_df = training(split) |> select(-strata)
test_df  = testing(split)  |> select(-strata)

rf_recipe = recipe(target_vector ~ ., data = train_df) |>
  step_rm(Conc., batch) |>
  step_zv(all_numeric_predictors()) |>
  step_impute_knn(all_numeric_predictors(), neighbors = 3) |>
  step_smote(target_vector, over_ratio = 1, skip = TRUE)

rf_spec = rand_forest(mtry = tune(), min_n = tune(), trees = 500) |>
  set_engine("ranger", importance = "permutation", seed = 2024) |>
  set_mode("classification")

rf_wf = workflow() |>
  add_recipe(rf_recipe) |>
  add_model(rf_spec)

folds = vfold_cv(train_df, v = 5, strata = target_vector)

p_pred = ncol(train_df) - 3L
mtry_values = sort(unique(pmax(1L, c(2L, floor(sqrt(p_pred)),
                                     floor(p_pred / 3L), floor(p_pred / 2L)))))

rf_grid = expand.grid(
  mtry  = mtry_values,
  min_n = c(2, 5, 10)
)

cv_metrics = metric_set(accuracy, kap, f_meas)

cat("\nTuning RF over grid (", nrow(rf_grid), " combos x 5 folds)...\n", sep = "")
tune_res = tune_grid(
  rf_wf,
  resamples = folds,
  grid      = rf_grid,
  metrics   = cv_metrics,
  control   = control_grid(save_pred = FALSE, verbose = FALSE)
)

cat("\nTop CV results:\n")
print(tune_res |> show_best(metric = "accuracy", n = 5))

best     = select_best(tune_res, metric = "accuracy")
final_wf = finalize_workflow(rf_wf, best)
cat("\nbest:\n"); print(best)

final_model = fit(final_wf, data = train_df)
test_aug    = augment(final_model, new_data = test_df)

prob_cols = grep("^\\.pred_", names(test_aug), value = TRUE)
prob_cols = setdiff(prob_cols, ".pred_class")

class_metrics = metric_set(accuracy, kap, f_meas, sens, spec)
prob_metrics  = metric_set(roc_auc, mn_log_loss)

class_results = test_aug |>
  class_metrics(truth = target_vector, estimate = .pred_class)

prob_results = test_aug |>
  prob_metrics(truth = target_vector, !!!syms(prob_cols))

perf = bind_rows(class_results, prob_results)
cat("\n===== Performance on held-out test set =====\n")
print(perf)

cat("\n===== Confusion matrix (test set) =====\n")
print(conf_mat(test_aug, truth = target_vector, estimate = .pred_class))

cat("\n===== Per-concentration accuracy =====\n")
percent_correct_by_conc = test_aug |>
  mutate(Conc. = round(Conc., 4),
         correct = as.integer(.pred_class == target_vector)) |>
  group_by(Conc.) |>
  summarise(num_correct     = sum(correct),
            total           = n(),
            percent_correct = round(num_correct / total * 100, 2),
            .groups         = "drop") |>
  arrange(Conc.)
print(percent_correct_by_conc, n = Inf)

# ---- Plots ----
dir.create("v3_outputs", showWarnings = FALSE)

# ROC curves
roc_data = test_aug |>
  roc_curve(truth = target_vector, !!!syms(prob_cols))

p_roc = autoplot(roc_data) +
  labs(title = "ROC curves — one-vs-all by variant",
       subtitle = paste0("Conc. Threshold = ", conc_threshold, " ng/ml")) +
  theme_minimal()
ggsave("v3_outputs/v3_roc_curves.png", p_roc, width = 7, height = 6, dpi = 150)

# Confusion matrix heatmap
p_cm = conf_mat(test_aug, truth = target_vector, estimate = .pred_class) |>
  autoplot(type = "heatmap") +
  labs(title = "Confusion matrix — test set",
       subtitle = paste0("Conc. Threshold = ", conc_threshold, " ng/ml"))
ggsave("v3_outputs/v3_conf_mat.png", p_cm, width = 6, height = 5, dpi = 150)

# Per-concentration percent correct
p_pc = percent_correct_by_conc |>
  mutate(Conc. = factor(Conc.)) |>
  ggplot(aes(Conc., percent_correct)) +
  geom_col(fill = "#266B6E") +
  geom_text(aes(label = percent_correct),
            vjust = 1.5, color = "white") +
  scale_y_continuous(labels = scales::number_format(accuracy = 0.001)) +
  labs(title    = "Combined set of augmented data",
       subtitle = paste0("Conc. Threshold = ", conc_threshold, " ng/ml"),
       x        = "Concentration",
       y        = "Percent Correct") +
  theme_classic()
ggsave("v3_outputs/v3_percent_correct_by_conc.png", p_pc, width = 9, height = 6, dpi = 150)

# Feature importance
extracted = extract_fit_engine(final_model)
imp_df = tibble(feature    = names(extracted$variable.importance),
                importance = extracted$variable.importance) |>
  arrange(desc(importance))

p_imp = imp_df |>
  ggplot(aes(reorder(feature, importance), importance)) +
  geom_col(fill = "#266B6E") +
  coord_flip() +
  theme_minimal() +
  labs(x = NULL, title = "Permutation importance — final RF")
ggsave("v3_outputs/v3_importance.png", p_imp, width = 6, height = 5, dpi = 150)

cat("\nWrote v3_outputs/v3_roc_curves.png, v3_conf_mat.png, v3_percent_correct_by_conc.png, v3_importance.png\n")
