library(tidymodels)
library(Boruta)
library(ranger)

# Boruta uses RF internally — RF is invariant to monotone transforms,
# so no BoxCox or normalization needed. step_zv removes zero-variance predictors.
make_boruta_recipe = function(formula, data) {
  recipe(formula, data = data) |>
    step_zv(all_numeric_predictors())
}


# Runs Boruta on train_data and returns all features ranked by mean importance.
# Tentative features are resolved via TentativeRoughFix (rough majority vote).
# Returns list: $ranking (tibble with feature, importance, decision) and
# $boruta_out (the resolved Boruta object for inspection / plotting).
# Boruta uses the randomForest backend (ntree arg); the eval RF uses ranger.
rank_boruta_features = function(train_data, target, n_trees, max_runs,
                                p_value, seed) {
  set.seed(seed)

  x = as.data.frame(select(train_data, -all_of(target)))
  y = train_data[[target]]

  boruta_fit = Boruta(x       = x,
                      y       = y,
                      maxRuns = max_runs,
                      pValue  = p_value,
                      doTrace = 0,
                      ntree   = n_trees)

  boruta_fit = TentativeRoughFix(boruta_fit)

  stats = attStats(boruta_fit)

  ranking = tibble(
    feature    = rownames(stats),
    importance = stats$meanImp,
    decision   = as.character(stats$decision)
  ) |>
    mutate(decision_rank = if_else(decision == "Confirmed", 1L, 2L)) |>
    arrange(decision_rank, desc(importance)) |>
    select(-decision_rank)

  list(ranking = ranking, boruta_out = boruta_fit)
}


#' run_boruta_fs
#'
#' Three-level split pipeline to select antibody pairs via Boruta,
#' mirroring the structure of run_rf_fs and run_lasso_fs:
#'
#'   FULL DATA
#'     ├── FINAL_TEST              (held out until the very end)
#'     └── TRAINING
#'           ├── FEATURE_VAL_TEST  (held out for feature count selection)
#'           └── FEATURE_VAL_TRAIN
#'                 ├── Boruta (max_runs shadow-feature RF iterations)
#'                 │     → confirmed feature ranking
#'                 ├── n_param_opt_split > 0: vfold_cv → tune_grid (eval RF mtry, min_n)
#'                 └── n_param_opt_split = 0: OOB error grid → best eval_mtry
#'                 → evaluate confirmed-feature subsets on FEATURE_VAL_TEST
#'         Retrain on all TRAINING with best N + best eval params
#'         Final evaluation on FINAL_TEST
#'
#' @param dab_df              Single dAb dataframe from main_preprocessor()
#' @param target              Name of the outcome column
#' @param threshold           Rows at or below this concentration are excluded
#' @param prop_final_test     Fraction of full data held out as the final test set
#' @param prop_feature_test   Fraction of TRAINING held out for feature count selection
#' @param n_param_opt_split   CV folds for eval RF mtry tuning on FEATURE_VAL_TRAIN.
#'                            0 = OOB error grid search over eval_mtry_grid (preferred
#'                            when n is very small).
#' @param n_top_features      Feature counts to evaluate; values exceeding the number
#'                            of Boruta-confirmed features are dropped
#' @param max_runs            Boruta maxRuns: maximum shadow-feature RF iterations
#' @param p_value             Boruta pValue: two-sided significance threshold for
#'                            feature accept/reject decisions
#' @param n_trees             Trees for both Boruta's internal RF and the eval RF
#' @param eval_mtry_grid      mtry values for the eval RF when n_param_opt_split = 0.
#'                            NULL = auto-compute from p: c(2, floor(sqrt(p)),
#'                            floor(p/3), floor(p/2))
#' @param eval_min_n          Node size for the eval RF; fixed for OOB path, searched
#'                            alongside mtry in the tune_grid path
#' @param n_mtry_levels       Grid resolution for eval mtry when n_param_opt_split > 0
#' @param min_n_candidates    min_n values searched when n_param_opt_split > 0
#' @param seed                RNG seed for reproducibility
#'
#' @return Named list. All tidy tibbles carry a `seed` column so results from
#'   multiple runs can be stacked with bind_rows(). run_boruta_fs_all_dabs()
#'   additionally appends a `dab` column. Key fields:
#'     run_summary          - one-row tibble; primary unit for cross-run stacking
#'     boruta_out           - resolved Boruta object (after TentativeRoughFix)
#'     feature_ranking      - importance and decision per pair (Confirmed first)
#'     feature_val_results  - accuracy/kap on FEATURE_VAL_TEST per feature count
#'     test_metrics         - final accuracy/kap on FINAL_TEST
#'     test_conf_mat        - confusion matrix on FINAL_TEST
#'     final_fit, tune_results (tibble of OOB errors or tune object),
#'     final_split, feature_split
run_boruta_fs = function(dab_df,
                         target             = "target_vector",
                         threshold          = 0,
                         prop_final_test    = 0.20,
                         prop_feature_test  = 0.20,
                         n_param_opt_split  = 0,
                         n_top_features     = c(3, 4, 5, 6, 7, 8, 10),
                         max_runs           = 100,
                         p_value            = 0.01,
                         n_trees            = 500,
                         eval_mtry_grid     = NULL,
                         eval_min_n         = 5,
                         n_mtry_levels      = 5,
                         min_n_candidates   = c(2, 5, 10, 20),
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

  # ---- Boruta feature ranking on full feature_val_train ------------------
  boruta_result   = rank_boruta_features(feature_val_train, target,
                                         n_trees, max_runs, p_value, seed)
  feature_ranking = boruta_result$ranking
  boruta_out      = boruta_result$boruta_out

  n_confirmed = sum(feature_ranking$decision == "Confirmed")
  if (n_confirmed == 0) stop(
    "Boruta confirmed no features at p_value = ", p_value,
    ". Try increasing p_value or max_runs."
  )

  # effective_max: evaluate up to max(n_top_features), capped at confirmed count.
  # Always included so every dAb has a datapoint at its own ceiling.
  effective_max = min(n_confirmed, max(n_top_features))
  valid_n = sort(unique(c(n_top_features[n_top_features <= effective_max],
                          effective_max)))

  # ---- Eval RF hyperparameter tuning on full feature_val_train ------------
  if (n_param_opt_split == 0) {
    prep_obj    = prep(make_boruta_recipe(as.formula(paste(target, "~ .")),
                                          feature_val_train),
                       training = feature_val_train)
    baked_train = bake(prep_obj, new_data = feature_val_train)

    x_train = as.data.frame(select(baked_train, -all_of(target)))
    y_train = baked_train[[target]]

    p = ncol(x_train)
    if (is.null(eval_mtry_grid)) {
      eval_mtry_grid = sort(unique(pmax(1L, c(2L, floor(sqrt(p)),
                                              floor(p / 3L), floor(p / 2L)))))
    }

    tune_results = map_dfr(eval_mtry_grid, \(m) {
      fit = ranger(x              = x_train,
                   y              = y_train,
                   num.trees      = n_trees,
                   mtry           = m,
                   min.node.size  = eval_min_n,
                   classification = TRUE,
                   seed           = seed)
      tibble(mtry = m, oob_error = fit$prediction.error)
    })

    best_mtry  = tune_results$mtry[which.min(tune_results$oob_error)]
    best_min_n = eval_min_n

  } else {
    param_folds = vfold_cv(feature_val_train, v = n_param_opt_split, strata = target)

    tune_wf = workflow() |>
      add_recipe(make_boruta_recipe(as.formula(paste(target, "~ .")),
                                    feature_val_train)) |>
      add_model(
        rand_forest(mtry = tune(), min_n = tune(), trees = n_trees) |>
          set_engine("ranger") |>
          set_mode("classification")
      )

    p = ncol(feature_val_train) - 1L
    param_grid = expand.grid(
      mtry  = floor(seq(2, max(2, floor(p / 2)), length.out = n_mtry_levels)),
      min_n = min_n_candidates
    ) |> distinct()

    tune_results = tune_grid(
      tune_wf,
      resamples = param_folds,
      grid      = param_grid,
      metrics   = metric_set(accuracy, kap)
    )

    best_params = select_best(tune_results, metric = "accuracy")
    best_mtry   = best_params$mtry
    best_min_n  = best_params$min_n
  }

  # ---- Feature count selection: train on feature_val_train, score on
  #      feature_val_test - a clean holdout that Boruta and mtry tuning never touched
  fv_metrics = metric_set(accuracy, kap, f_meas)

  feature_val_results = map_dfr(valid_n, function(n) {
    top_features   = feature_ranking$feature[seq_len(n)]
    formula        = as.formula(
      paste(target, "~", paste(top_features, collapse = " + "))
    )
    effective_mtry = min(best_mtry, n)

    fit = workflow() |>
      add_recipe(make_boruta_recipe(formula, feature_val_train)) |>
      add_model(
        rand_forest(mtry = effective_mtry, min_n = best_min_n, trees = n_trees) |>
          set_engine("ranger", seed = seed) |>
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
  # feature_val_test rows are now included; Boruta, mtry, and features are
  # already fixed so this does not constitute leakage.
  final_formula = as.formula(
    paste(target, "~", paste(best_features, collapse = " + "))
  )

  final_fit = workflow() |>
    add_recipe(make_boruta_recipe(final_formula, training_data)) |>
    add_model(
      rand_forest(mtry  = min(best_mtry, best_n),
                  min_n = best_min_n,
                  trees = n_trees) |>
        set_engine("ranger", seed = seed) |>
        set_mode("classification")
    ) |>
    fit(data = training_data)

  augmented_final_test = augment(final_fit, new_data = final_test)

  # ---- Final-test curve: retrain on training_data at each valid_n ---------
  # Evaluates every feature count on the true holdout so performance-by-feature
  # plots are on the same data as test_metrics (the table). best_n is the
  # selected point on this curve; feature_val_results was the selection tool.
  final_test_results = map_dfr(valid_n, function(n) {
    top_features   = feature_ranking$feature[seq_len(n)]
    formula        = as.formula(paste(target, "~",
                                      paste(top_features, collapse = " + ")))
    effective_mtry = min(best_mtry, n)

    fit = workflow() |>
      add_recipe(make_boruta_recipe(formula, training_data)) |>
      add_model(
        rand_forest(mtry  = effective_mtry,
                    min_n = best_min_n,
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
  # run_summary is the primary unit for stacking across seeds and dAbs.
  # best_features stored as a list-column for unnest_longer() frequency analysis.
  run_summary = tibble(
    seed            = seed,
    best_eval_mtry  = best_mtry,
    best_eval_min_n = best_min_n,
    best_n_features = best_n,
    best_features   = list(best_features),
    n_confirmed     = n_confirmed
  )

  list(
    run_summary          = run_summary,
    final_split          = final_split,
    feature_split        = feature_split,
    boruta_out           = boruta_out,
    tune_results         = tune_results,
    best_eval_mtry       = best_mtry,
    best_eval_min_n      = best_min_n,
    feature_ranking      = mutate(feature_ranking,     seed = seed),
    feature_val_results  = mutate(feature_val_results, seed = seed),
    final_test_results   = mutate(final_test_results,  seed = seed),
    best_n_features      = best_n,
    best_features        = best_features,
    final_fit            = final_fit,
    test_metrics         = mutate(
                             metrics(augmented_final_test,
                                     truth    = !!sym(target),
                                     estimate = .pred_class),
                             seed = seed
                           ),
    test_conf_mat        = conf_mat(augmented_final_test,
                                    truth    = !!sym(target),
                                    estimate = .pred_class)
  )
}


#' run_boruta_fs_all_dabs
#'
#' Applies run_boruta_fs across a named list of dAb dataframes. Appends a `dab`
#' column to all tidy tibbles so results from multiple calls (different seeds)
#' can be stacked directly.
#'
#' To compile uncertainty estimates across seeds:
#'   runs = map(c(1969, 1986, 2024), \(s) run_boruta_fs_all_dabs(dab_list, seed = s))
#'   all_summaries  = map(runs, \(r) map_dfr(r, \(x) x$run_summary))  |> bind_rows()
#'   all_metrics    = map(runs, \(r) map_dfr(r, \(x) x$test_metrics)) |> bind_rows()
#'
#' @param dab_df_list Named list of dAb dataframes (e.g. from main_preprocessor())
#' @param ...         Passed to run_boruta_fs
run_boruta_fs_all_dabs = function(dab_df_list, ...) {
  imap(dab_df_list, function(df, name) {
    message("Boruta FS: ", name)
    result = run_boruta_fs(df, ...)

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
#' +/-1 SE error bars aggregated across all (seed x dAb) combinations. Mirrors
#' the engineering-constraint view: fewer pairs on the right, more on the left.
#'
#' Identical function exists in rf_tidymodels.R and lasso_tidymodels.R -- works
#' with any of their outputs.
#'
#' @param feature_val_results Stacked feature_val_results tibble from
#'   run_boruta_fs_all_dabs() calls across multiple seeds. Must contain columns
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
