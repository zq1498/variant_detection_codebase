library(tidymodels)
library(ranger)
library(mRMRe)

# Ranks antibody pairs via minimum Redundancy Maximum Relevance (mRMR).
# Greedy forward selection maximizes MI with the target while penalizing
# redundancy with already-selected features. Model-free: no classifier is fit.
# Target factor is integer-encoded for MI estimation; ordering is arbitrary and
# does not affect the joint distribution used in MI computation.
rank_mrmr_features = function(train_data, target) {
  feature_cols = setdiff(colnames(train_data), target)

  # Zero-variance features have zero MI with any target — drop them first
  feature_sds  = sapply(train_data[feature_cols], sd)
  feature_cols = feature_cols[feature_sds > 0]

  # mRMRe requires numeric, ordered factor, or Surv columns.
  # Build with direct indexing: features (double) + target (double integer codes).
  # Target goes last so solutions() indices map directly to feature column positions.
  mrmr_df = as.data.frame(train_data[, feature_cols, drop = FALSE])
  mrmr_df[[target]] = as.numeric(as.integer(as.factor(train_data[[target]])))

  target_col_idx = ncol(mrmr_df)
  dd             = mRMR.data(data = mrmr_df)
  result         = mRMR.classic(data           = dd,
                                target_indices = as.integer(target_col_idx),
                                feature_count  = as.integer(length(feature_cols)))

  # solutions() returns 1-indexed column positions in mrmr_df;
  # since features are columns 1:p and target is column p+1, these are feature indices
  selected_idx = solutions(result)[[1]]

  tibble(
    feature    = colnames(mrmr_df)[selected_idx],
    mrmr_rank  = seq_along(selected_idx),
    mrmr_score = scores(result)[[1]]
  )
}


make_mrmr_eval_recipe = function(formula, data) {
  recipe(formula, data = data) |>
    step_zv(all_numeric_predictors())
}


#' run_mrmr_fs
#'
#' Three-level split pipeline for mRMR-based feature selection, mirroring the
#' structure of run_rf_fs and run_lasso_fs.
#'
#'   FULL DATA
#'     ├── FINAL_TEST              (held out until the very end)
#'     └── TRAINING
#'           ├── FEATURE_VAL_TEST  (held out for feature count selection)
#'           └── FEATURE_VAL_TRAIN
#'                 → mRMR ranking on full FEATURE_VAL_TRAIN (no tuning — filter method)
#'                 → evaluate feature subsets on FEATURE_VAL_TEST via RF
#'         Retrain on all TRAINING with best N
#'         Final evaluation on FINAL_TEST
#'
#' mRMR is a filter method: ranking is model-independent. Feature-count evaluation
#' uses RF with sqrt(p) mtry (standard classification default, no tuning). This
#' mirrors the common "mRMR pre-filter → classifier evaluation" pattern in the
#' literature (Peng et al. 2005; De Jay et al. 2013).
#'
#' @param dab_df              Single dAb dataframe from main_preprocessor()
#' @param target              Name of the outcome column
#' @param threshold           Rows at or below this concentration are excluded
#' @param prop_final_test     Fraction of full data held out as the final test set
#' @param prop_feature_test   Fraction of TRAINING held out for feature count selection
#' @param n_top_features      Feature counts to evaluate
#' @param n_trees             Number of trees for the RF evaluation model
#' @param min_n               Node size for the RF evaluation model
#' @param seed                RNG seed for reproducibility
#'
#' @return Named list mirroring run_rf_fs output. feature_ranking carries mrmr_rank
#'   and mrmr_score instead of importance. No tune_results (filter method, no tuning).
run_mrmr_fs = function(dab_df,
                       target            = "target_vector",
                       threshold         = 0,
                       prop_final_test   = 0.20,
                       prop_feature_test = 0.20,
                       n_top_features    = c(3, 4, 5, 6, 7, 8, 10),
                       n_trees           = 500,
                       min_n             = 5,
                       seed              = 2024) {

  set.seed(seed)

  clean_df = dab_df |>
    filter(Conc. > threshold) |>
    select(-c(Conc., Replicates)) |>
    mutate(across(all_of(target), as.factor))

  # ---- Level 1: Final holdout test set ------------------------------------
  final_split   = initial_split(clean_df, prop = 1 - prop_final_test, strata = target)
  training_data = training(final_split)
  final_test    = testing(final_split)

  # ---- Level 2: Feature selection validation split -----------------------
  feature_split     = initial_split(training_data, prop = 1 - prop_feature_test,
                                    strata = target)
  feature_val_train = training(feature_split)
  feature_val_test  = testing(feature_split)

  # ---- Feature ranking via mRMR (no hyperparameters to tune) -------------
  feature_ranking = rank_mrmr_features(feature_val_train, target)

  n_available   = nrow(feature_ranking)
  effective_max = min(n_available, max(n_top_features))
  valid_n       = sort(unique(c(n_top_features[n_top_features <= effective_max],
                                effective_max)))

  # ---- Feature count selection: train on feature_val_train, score on
  #      feature_val_test using RF with sqrt(p) mtry ------------------------
  fv_metrics = metric_set(accuracy, kap, f_meas)

  feature_val_results = map_dfr(valid_n, function(n) {
    top_features = feature_ranking$feature[seq_len(n)]
    formula      = as.formula(paste(target, "~", paste(top_features, collapse = " + ")))

    fit = workflow() |>
      add_recipe(make_mrmr_eval_recipe(formula, feature_val_train)) |>
      add_model(
        rand_forest(mtry = max(1L, floor(sqrt(n))), min_n = min_n, trees = n_trees) |>
          set_engine("ranger", importance = "permutation", seed = seed) |>
          set_mode("classification")
      ) |>
      fit(data = feature_val_train)

    augment(fit, new_data = feature_val_test) |>
      fv_metrics(truth = !!sym(target), estimate = .pred_class) |>
      mutate(n_features = n)
  })

  # ---- Auto-select best N by feature_val_test accuracy -------------------
  best_n = feature_val_results |>
    filter(.metric == "accuracy") |>
    slice_max(.estimate, n = 1, with_ties = FALSE) |>
    pull(n_features)

  best_features = feature_ranking$feature[seq_len(best_n)]

  # ---- Final model: retrain on all TRAINING data -------------------------
  final_formula = as.formula(paste(target, "~", paste(best_features, collapse = " + ")))

  final_fit = workflow() |>
    add_recipe(make_mrmr_eval_recipe(final_formula, training_data)) |>
    add_model(
      rand_forest(mtry  = max(1L, floor(sqrt(best_n))),
                  min_n = min_n,
                  trees = n_trees) |>
        set_engine("ranger", seed = seed) |>
        set_mode("classification")
    ) |>
    fit(data = training_data)

  augmented_final_test = augment(final_fit, new_data = final_test)

  # ---- Final-test curve: retrain on training_data at each valid_n ---------
  final_test_results = map_dfr(valid_n, function(n) {
    top_features = feature_ranking$feature[seq_len(n)]
    formula      = as.formula(paste(target, "~",
                                    paste(top_features, collapse = " + ")))

    fit = workflow() |>
      add_recipe(make_mrmr_eval_recipe(formula, training_data)) |>
      add_model(
        rand_forest(mtry  = max(1L, floor(sqrt(n))),
                    min_n = min_n,
                    trees = n_trees) |>
          set_engine("ranger", seed = seed) |>
          set_mode("classification")
      ) |>
      fit(data = training_data)

    augment(fit, new_data = final_test) |>
      fv_metrics(truth = !!sym(target), estimate = .pred_class) |>
      mutate(n_features = n)
  })

  # ---- Build output --------------------------------------------------------
  run_summary = tibble(
    seed                 = seed,
    best_n_features      = best_n,
    best_features        = list(best_features),
    n_features_available = n_available
  )

  list(
    run_summary         = run_summary,
    final_split         = final_split,
    feature_split       = feature_split,
    feature_ranking     = mutate(feature_ranking,     seed = seed),
    feature_val_results = mutate(feature_val_results, seed = seed),
    final_test_results  = mutate(final_test_results,  seed = seed),
    best_n_features     = best_n,
    best_features       = best_features,
    final_fit           = final_fit,
    test_metrics        = mutate(
                            metrics(augmented_final_test,
                                    truth    = !!sym(target),
                                    estimate = .pred_class),
                            seed = seed
                          ),
    test_conf_mat       = conf_mat(augmented_final_test,
                                   truth    = !!sym(target),
                                   estimate = .pred_class)
  )
}


#' run_mrmr_fs_all_dabs
#'
#' Applies run_mrmr_fs across a named list of dAb dataframes. Appends a `dab`
#' column to all tidy tibbles so results from multiple calls (different seeds)
#' can be stacked directly.
#'
#' @param dab_df_list Named list of dAb dataframes (e.g. from main_preprocessor())
#' @param ...         Passed to run_mrmr_fs
run_mrmr_fs_all_dabs = function(dab_df_list, ...) {
  imap(dab_df_list, function(df, name) {
    message("mRMR FS: ", name)
    result = run_mrmr_fs(df, ...)

    result$run_summary         = mutate(result$run_summary,         dab = name)
    result$feature_ranking     = mutate(result$feature_ranking,     dab = name)
    result$feature_val_results = mutate(result$feature_val_results, dab = name)
    result$final_test_results  = mutate(result$final_test_results,  dab = name)
    result$test_metrics        = mutate(result$test_metrics,        dab = name)

    result
  })
}


#' plot_perf_by_features
#'
#' Identical to the version in rf_tidymodels.R and lasso_tidymodels.R.
#' Works with mRMR output since feature_val_results has the same schema.
plot_perf_by_features = function(feature_val_results, title = NULL,
                                 facet_by_dab = FALSE) {
  plot_df = feature_val_results |>
    filter(.metric %in% c("accuracy", "kap", "f_meas")) |>
    mutate(metric_label = case_when(
      .metric == "accuracy" ~ "Mean Accuracy",
      .metric == "kap"      ~ "Mean Kappa",
      .metric == "f_meas"   ~ "Mean Average F1"
    ))

  group_vars = if (facet_by_dab) c("dab", "n_features", "metric_label") else
                                 c("n_features", "metric_label")

  summary_df = plot_df |>
    group_by(across(all_of(group_vars))) |>
    summarise(mean_val = mean(.estimate),
              se_val   = sd(.estimate) / sqrt(n()),
              .groups  = "drop")

  p = ggplot(summary_df,
             aes(x = n_features, y = mean_val,
                 color = metric_label, group = metric_label)) +
    geom_line() +
    geom_point() +
    geom_errorbar(aes(ymin = mean_val - se_val,
                      ymax = mean_val + se_val),
                  width = 0.3) +
    scale_x_reverse(breaks = sort(unique(summary_df$n_features),
                                  decreasing = TRUE),
                    limits = c(max(summary_df$n_features),
                               min(summary_df$n_features))) +
    scale_y_continuous(limits = c(0, 1)) +
    scale_color_manual(
      values = c("Mean Accuracy"   = "#6BAED6",
                 "Mean Average F1" = "#2171B5",
                 "Mean Kappa"      = "#08306B"),
      name   = "Performance Metric"
    ) +
    labs(title = title,
         x     = "Number of antibody pairs in the model",
         y     = "Mean Performance Metric Value") +
    theme_bw() +
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank())

  if (facet_by_dab) p = p + facet_wrap(~dab, ncol = 4)

  p
}
