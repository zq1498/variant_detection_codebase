library(tidymodels)
library(glmnet)

# Preprocessing: BoxCox then normalize. BoxCox requires strictly positive inputs —
# if any fluorescence readings are <= 0, swap step_BoxCox for step_YeoJohnson.
make_lasso_recipe = function(formula, data) {
  recipe(formula, data = data) |>
    step_BoxCox(all_numeric_predictors()) |>
    step_normalize(all_numeric_predictors())
}

# Ranks antibody pairs by mean absolute LASSO coefficient across variant classes.
# Run on feature_val_train — downstream of both holdout splits.
rank_lasso_features = function(train_data, target, best_penalty, mixture = 1) {
  formula = as.formula(paste(target, "~ ."))

  fit_obj = workflow() |>
    add_recipe(make_lasso_recipe(formula, train_data)) |>
    add_model(
      multinom_reg(penalty = best_penalty, mixture = mixture) |> set_engine("glmnet", standardize = FALSE)
    ) |>
    fit(data = train_data)

  coef_list = coef(extract_fit_engine(fit_obj), s = best_penalty)

  map_dfr(coef_list,
    \(mat) tibble(feature = rownames(mat), coef = as.vector(as.matrix(mat))),
    .id = "class"
  ) |>
    filter(feature != "(Intercept)") |>
    group_by(feature) |>
    summarise(importance = mean(abs(coef)), .groups = "drop") |>
    arrange(desc(importance))
}


#' run_lasso_fs
#'
#' Three-level split pipeline to select antibody pairs via LASSO:
#'
#'   FULL DATA
#'     ├── FINAL_TEST              (held out until the very end)
#'     └── TRAINING
#'           ├── FEATURE_VAL_TEST  (held out for feature count selection)
#'           └── FEATURE_VAL_TRAIN
#'                 ├── n_param_opt_split > 0: vfold_cv → tune_grid → best lambda
#'                 └── n_param_opt_split = 0: cv.glmnet internal CV → lambda.1se
#'                 → feature ranking on full FEATURE_VAL_TRAIN
#'                 → evaluate feature subsets on FEATURE_VAL_TEST
#'         Retrain on all TRAINING with best N + best lambda
#'         Final evaluation on FINAL_TEST
#'
#' @param dab_df              Single dAb dataframe from main_preprocessor()
#' @param target              Name of the outcome column
#' @param threshold           Rows at or below this concentration are excluded
#' @param prop_final_test     Fraction of full data held out as the final test set
#' @param prop_feature_test   Fraction of TRAINING held out for feature count selection
#' @param n_param_opt_split   CV folds for lambda tuning on FEATURE_VAL_TRAIN.
#'                            0 = delegate to cv.glmnet's internal CV (uses all of
#'                            FEATURE_VAL_TRAIN; preferred when n is very small).
#'                            Note: when 0, BoxCox params are estimated on all of
#'                            FEATURE_VAL_TRAIN before glmnet's internal folds —
#'                            a minor approximation acceptable for fluorescence data.
#' @param n_top_features      Feature counts to evaluate; values exceeding the
#'                            number of nonzero LASSO features are dropped
#' @param lambda_levels       Grid size for penalty tuning (n_param_opt_split > 0 only)
#' @param lambda_range        log10 bounds for the penalty grid (n_param_opt_split > 0 only)
#' @param mixture             Elastic net mixing parameter: 1 = pure LASSO (default,
#'                            maximally sparse), 0 = ridge (no zeros). Values between 0
#'                            and 1 give elastic net — the L1 component still forces zeros,
#'                            but the L2 term stabilizes selection among correlated features.
#'                            Use mixture = 1 for feature selection; other values for
#'                            robustness comparisons.
#' @param seed                RNG seed for reproducibility
#'
#' @return Named list. All tidy tibbles carry a `seed` column so results from
#'   multiple runs can be stacked with bind_rows(). run_lasso_fs_all_dabs()
#'   additionally appends a `dab` column. Key fields:
#'     run_summary          — one-row tibble; primary unit for cross-run stacking
#'     feature_ranking      — importance scores per pair
#'     feature_val_results  — accuracy/kap on FEATURE_VAL_TEST per feature count
#'     test_metrics         — final accuracy/kap on FINAL_TEST
#'     test_conf_mat        — confusion matrix on FINAL_TEST
#'     final_fit, tune_results (NULL when n_param_opt_split = 0),
#'     final_split, feature_split
run_lasso_fs = function(dab_df,
                        target             = "target_vector",
                        threshold          = 0,
                        prop_final_test    = 0.20,
                        prop_feature_test  = 0.20,
                        n_param_opt_split  = 4,
                        n_top_features     = c(3, 4, 5, 6, 7, 8, 10),
                        lambda_levels      = 30,
                        lambda_range       = c(-4, 0),
                        mixture            = 1,
                        seed               = 2024) {

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

  # ---- Lambda tuning -------------------------------------------------------
  if (n_param_opt_split == 0) {
    # Preprocess once on all of feature_val_train, then let cv.glmnet handle
    # its own internal fold structure. This keeps the full feature_val_train
    # available rather than carving out another explicit split.
    prep_obj    = prep(make_lasso_recipe(as.formula(paste(target, "~ .")),
                                         feature_val_train),
                       training = feature_val_train)
    baked_train = bake(prep_obj, new_data = feature_val_train)

    x_train = as.matrix(select(baked_train, -all_of(target)))
    y_train = baked_train[[target]]

    cv_fit       = cv.glmnet(x_train, y_train, family = "multinomial", alpha = mixture,
                             type.multinomial = "grouped", type.measure = "class",
                             standardize = FALSE)
    best_penalty = cv_fit$lambda.1se
    tune_results = NULL

  } else {
    lambda_folds = vfold_cv(feature_val_train, v = n_param_opt_split, strata = target)

    tune_wf = workflow() |>
      add_recipe(make_lasso_recipe(as.formula(paste(target, "~ .")), feature_val_train)) |>
      add_model(multinom_reg(penalty = tune(), mixture = mixture) |> set_engine("glmnet", standardize = FALSE))

    lambda_grid = grid_regular(penalty(range = lambda_range), levels = lambda_levels)

    tune_results = tune_grid(
      tune_wf,
      resamples = lambda_folds,
      grid      = lambda_grid,
      metrics   = metric_set(accuracy, kap)
    )

    # 1-SE rule: sparsest model within one SE of peak CV accuracy
    best_penalty = select_by_one_std_err(
      tune_results, metric = "accuracy", desc(penalty)
    )$penalty
  }

  # ---- Feature ranking on full feature_val_train at best lambda -----------
  feature_ranking = rank_lasso_features(feature_val_train, target, best_penalty, mixture)

  n_nonzero = sum(feature_ranking$importance > 0)
  if (n_nonzero == 0) stop(
    "LASSO selected no nonzero features at best_penalty = ", best_penalty,
    ". Try widening lambda_range or reducing n_param_opt_split."
  )
  # effective_max: evaluate up to 15 (or max of n_top_features grid), capped at
  # the number of nonzero features available. Always included so every dAb has
  # a datapoint at its own ceiling even when that ceiling is below 15.
  effective_max = min(n_nonzero, max(n_top_features))
  valid_n = sort(unique(c(n_top_features[n_top_features <= effective_max],
                          effective_max)))

  # ---- Feature count selection: train on feature_val_train, score on
  #      feature_val_test — a clean holdout that lambda tuning never touched --
  fv_metrics = metric_set(accuracy, kap, f_meas)

  feature_val_results = map_dfr(valid_n, function(n) {
    top_features = feature_ranking$feature[seq_len(n)]
    formula = as.formula(
      paste(target, "~", paste(top_features, collapse = " + "))
    )

    fit = workflow() |>
      add_recipe(make_lasso_recipe(formula, feature_val_train)) |>
      add_model(
        multinom_reg(penalty = best_penalty, mixture = mixture) |> set_engine("glmnet", standardize = FALSE)
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
  # feature_val_test rows are now included; lambda and features are already
  # fixed so this does not constitute leakage.
  final_formula = as.formula(
    paste(target, "~", paste(best_features, collapse = " + "))
  )

  final_fit = workflow() |>
    add_recipe(make_lasso_recipe(final_formula, training_data)) |>
    add_model(
      multinom_reg(penalty = best_penalty, mixture = mixture) |> set_engine("glmnet", standardize = FALSE)
    ) |>
    fit(data = training_data)

  augmented_final_test = augment(final_fit, new_data = final_test)

  # ---- Final-test curve: retrain on training_data at each valid_n ---------
  # Evaluates every feature count on the true holdout so performance-by-feature
  # plots are on the same data as test_metrics (the table). best_n is the
  # selected point on this curve; feature_val_results was the selection tool.
  final_test_results = map_dfr(valid_n, function(n) {
    top_features = feature_ranking$feature[seq_len(n)]
    formula      = as.formula(paste(target, "~",
                                    paste(top_features, collapse = " + ")))

    fit = workflow() |>
      add_recipe(make_lasso_recipe(formula, training_data)) |>
      add_model(
        multinom_reg(penalty = best_penalty, mixture = mixture) |>
          set_engine("glmnet", standardize = FALSE)
      ) |>
      fit(data = training_data)

    augment(fit, new_data = final_test) |>
      fv_metrics(truth = !!sym(target), estimate = .pred_class) |>
      mutate(n_features = n)
  })

  # ---- Build output --------------------------------------------------------
  # run_summary is the primary unit for stacking across seeds and dAbs.
  # best_features is stored as a list-column so summaries can be row-bound
  # directly and then unnest_longer("best_features") for frequency analysis.
  run_summary = tibble(
    seed               = seed,
    best_penalty       = best_penalty,
    mixture            = mixture,
    best_n_features    = best_n,
    best_features      = list(best_features),
    n_nonzero_features = n_nonzero
  )

  list(
    run_summary         = run_summary,
    final_split         = final_split,
    feature_split       = feature_split,
    tune_results        = tune_results,
    best_penalty        = best_penalty,
    feature_ranking     = mutate(feature_ranking,     seed = seed),
    feature_val_results  = mutate(feature_val_results,  seed = seed),
    final_test_results   = mutate(final_test_results,   seed = seed),
    best_n_features      = best_n,
    best_features        = best_features,
    final_fit            = final_fit,
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


#' run_lasso_fs_all_dabs
#'
#' Applies run_lasso_fs across a named list of dAb dataframes. Appends a `dab`
#' column to all tidy tibbles so results from multiple calls (different seeds)
#' can be stacked directly.
#'
#' To compile uncertainty estimates across seeds:
#'   runs = map(c(2024, 2025, 2026), \(s) run_lasso_fs_all_dabs(dab_list, seed = s))
#'   all_summaries  = map(runs, \(r) map_dfr(r, \(x) x$run_summary))  |> bind_rows()
#'   all_metrics    = map(runs, \(r) map_dfr(r, \(x) x$test_metrics)) |> bind_rows()
#'
#' @param dab_df_list Named list of dAb dataframes (e.g. from main_preprocessor())
#' @param ...         Passed to run_lasso_fs
run_lasso_fs_all_dabs = function(dab_df_list, ...) {
  imap(dab_df_list, function(df, name) {
    message("LASSO FS: ", name)
    result = run_lasso_fs(df, ...)

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
#' Line plot of mean accuracy, kappa, and macro-F1 vs. feature count, with
#' ±1 SE error bars aggregated across all (seed × dAb) combinations. Mirrors
#' the engineering-constraint view: fewer pairs on the right, more on the left.
#'
#' Identical function exists in rf_tidymodels.R — works with either output.
#'
#' @param feature_val_results Stacked feature_val_results tibble from
#'   run_lasso_fs_all_dabs() calls across multiple seeds. Must contain columns
#'   n_features, .metric, .estimate. A `dab` column is required when
#'   facet_by_dab = TRUE.
#' @param title    Optional plot title string
#' @param facet_by_dab If TRUE, produces a per-dAb faceted plot instead of
#'   one aggregate panel; error bars then reflect SE across seeds only
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
