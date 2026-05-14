library(tidymodels)
library(discrim)   # discrim_flexible()
library(earth)     # evimp()

# earth loads MASS, which masks dplyr::select. Pin it here so all select()
# calls in this file resolve to dplyr regardless of load order.
select = dplyr::select

# scale + center matches the original caret preprocessing for MARS.
make_mars_recipe = function(formula, data) {
  recipe(formula, data = data) |>
    step_normalize(all_numeric_predictors())
}

# Ranks antibody pairs by GCV-based variable importance from the earth submodel
# inside the mda::fda fit. Run on feature_val_train — downstream of both holdout
# splits. Features earth never selected receive importance = 0 and sort last.
rank_mars_features = function(train_data, target, best_num_terms, best_prod_degree) {
  all_features = setdiff(names(train_data), target)
  formula      = as.formula(paste(target, "~ ."))

  fit = workflow() |>
    add_recipe(make_mars_recipe(formula, train_data)) |>
    add_model(
      discrim_flexible(num_terms = best_num_terms, prod_degree = best_prod_degree) |>
        set_engine("earth") |>
        set_mode("classification")
    ) |>
    fit(data = train_data)

  # discrim_flexible wraps mda::fda; the earth submodel lives in $fit
  earth_fit = extract_fit_engine(fit)$fit
  imp_mat   = earth::evimp(earth_fit, trim = FALSE)

  used = tibble(feature    = rownames(imp_mat),
                importance = imp_mat[, "gcv"])

  tibble(feature = all_features) |>
    left_join(used, by = "feature") |>
    mutate(importance = replace_na(importance, 0)) |>
    arrange(desc(importance))
}


#' run_mars_fs
#'
#' Three-level split pipeline to select antibody pairs via MARS (Flexible
#' Discriminant Analysis via earth), mirroring run_rf_fs / run_lasso_fs:
#'
#'   FULL DATA
#'     ├── FINAL_TEST              (held out until the very end)
#'     └── TRAINING
#'           ├── FEATURE_VAL_TEST  (held out for feature count selection)
#'           └── FEATURE_VAL_TRAIN
#'                 ├── n_param_opt_split > 0: vfold_cv → tune_grid (num_terms,
#'                 │   prod_degree)
#'                 └── n_param_opt_split = 0: conservative defaults, no CV
#'                 → feature ranking (GCV importance) on full FEATURE_VAL_TRAIN
#'                 → evaluate feature subsets on FEATURE_VAL_TEST
#'         Retrain on all TRAINING with best N + best params
#'         Final evaluation on FINAL_TEST
#'
#' Uses discrim::discrim_flexible(engine = "earth") so that the 4-class outcome
#' is handled natively. parsnip::mars() is binary-only; discrim_flexible is the
#' tidymodels-recommended multiclass MARS path.
#'
#' @param dab_df               Single dAb dataframe from main_preprocessor()
#' @param target               Name of the outcome column
#' @param threshold            Rows at or below this concentration are excluded
#' @param prop_final_test      Fraction of full data held out as the final test set
#' @param prop_feature_test    Fraction of TRAINING held out for feature count selection
#' @param n_param_opt_split    CV folds for hyperparameter tuning on FEATURE_VAL_TRAIN.
#'                             0 = skip tuning and use (num_terms=10, prod_degree=1).
#' @param n_top_features       Feature counts to evaluate; values exceeding the
#'                             number of features available are dropped
#' @param num_terms_grid       num_terms (nprune) values to search
#' @param prod_degree_candidates  Interaction degree values to search (1 = main
#'                             effects only, 2 = pairwise)
#' @param seed                 RNG seed for reproducibility
#'
#' @return Named list with the same structure as run_rf_fs / run_lasso_fs. Key fields:
#'     run_summary          — one-row tibble; primary unit for cross-run stacking
#'     feature_ranking      — GCV importance per pair
#'     feature_val_results  — accuracy/kap/f_meas on FEATURE_VAL_TEST per feature count
#'     test_metrics         — final accuracy/kap on FINAL_TEST
#'     test_conf_mat        — confusion matrix on FINAL_TEST
#'     final_fit, tune_results (NULL when n_param_opt_split = 0),
#'     final_split, feature_split
run_mars_fs = function(dab_df,
                       target                 = "target_vector",
                       threshold              = 0,
                       prop_final_test        = 0.20,
                       prop_feature_test      = 0.20,
                       n_param_opt_split      = 5,
                       n_top_features         = c(3, 4, 5, 6, 7, 8, 10),
                       num_terms_grid         = c(5, 10, 15, 20),
                       prod_degree_candidates = c(1, 2),
                       seed                   = 2024) {

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

  # ---- Hyperparameter tuning -----------------------------------------------
  if (n_param_opt_split == 0) {
    best_num_terms   = 10
    best_prod_degree = 1
    tune_results     = NULL

  } else {
    param_folds = vfold_cv(feature_val_train, v = n_param_opt_split, strata = target)

    tune_wf = workflow() |>
      add_recipe(make_mars_recipe(as.formula(paste(target, "~ .")),
                                  feature_val_train)) |>
      add_model(
        discrim_flexible(num_terms = tune(), prod_degree = tune()) |>
          set_engine("earth") |>
          set_mode("classification")
      )

    param_grid = expand.grid(
      num_terms   = num_terms_grid,
      prod_degree = prod_degree_candidates
    )

    tune_results = tune_grid(
      tune_wf,
      resamples = param_folds,
      grid      = param_grid,
      metrics   = metric_set(accuracy, kap)
    )

    best_params      = select_best(tune_results, metric = "accuracy")
    best_num_terms   = best_params$num_terms
    best_prod_degree = best_params$prod_degree
  }

  # ---- Feature ranking on full feature_val_train at best params -----------
  feature_ranking = rank_mars_features(feature_val_train, target,
                                       best_num_terms, best_prod_degree)

  n_available   = nrow(feature_ranking)
  effective_max = min(n_available, max(n_top_features))
  valid_n       = sort(unique(c(n_top_features[n_top_features <= effective_max],
                                effective_max)))

  # ---- Feature count selection: train on feature_val_train, score on
  #      feature_val_test — a clean holdout that tuning never touched ----------
  fv_metrics = metric_set(accuracy, kap, f_meas)

  feature_val_results = map_dfr(valid_n, function(n) {
    top_features = feature_ranking$feature[seq_len(n)]
    formula      = as.formula(paste(target, "~",
                                    paste(top_features, collapse = " + ")))

    fit = workflow() |>
      add_recipe(make_mars_recipe(formula, feature_val_train)) |>
      add_model(
        discrim_flexible(num_terms   = best_num_terms,
                         prod_degree = best_prod_degree) |>
          set_engine("earth") |>
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
  final_formula = as.formula(
    paste(target, "~", paste(best_features, collapse = " + "))
  )

  final_fit = workflow() |>
    add_recipe(make_mars_recipe(final_formula, training_data)) |>
    add_model(
      discrim_flexible(num_terms   = best_num_terms,
                       prod_degree = best_prod_degree) |>
        set_engine("earth") |>
        set_mode("classification")
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
      add_recipe(make_mars_recipe(formula, training_data)) |>
      add_model(
        discrim_flexible(num_terms   = best_num_terms,
                         prod_degree = best_prod_degree) |>
          set_engine("earth") |>
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
    best_num_terms       = best_num_terms,
    best_prod_degree     = best_prod_degree,
    best_n_features      = best_n,
    best_features        = list(best_features),
    n_features_available = n_available
  )

  list(
    run_summary         = run_summary,
    final_split         = final_split,
    feature_split       = feature_split,
    tune_results        = tune_results,
    best_num_terms      = best_num_terms,
    best_prod_degree    = best_prod_degree,
    feature_ranking     = mutate(feature_ranking,     seed = seed),
    feature_val_results  = mutate(feature_val_results,  seed = seed),
    final_test_results   = mutate(final_test_results,   seed = seed),
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


#' run_mars_fs_all_dabs
#'
#' Applies run_mars_fs across a named list of dAb dataframes. Appends a `dab`
#' column to all tidy tibbles so results from multiple calls (different seeds)
#' can be stacked directly.
#'
#' To compile uncertainty estimates across seeds:
#'   runs = map(c(1969, 1986, 2024), \(s) run_mars_fs_all_dabs(dab_list, seed = s))
#'   all_summaries  = map(runs, \(r) map_dfr(r, \(x) x$run_summary))  |> bind_rows()
#'   all_metrics    = map(runs, \(r) map_dfr(r, \(x) x$test_metrics)) |> bind_rows()
#'
#' @param dab_df_list Named list of dAb dataframes (e.g. from main_preprocessor())
#' @param ...         Passed to run_mars_fs
run_mars_fs_all_dabs = function(dab_df_list, ...) {
  imap(dab_df_list, function(df, name) {
    message("MARS FS: ", name)
    result = run_mars_fs(df, ...)

    result$run_summary          = mutate(result$run_summary,          dab = name)
    result$feature_ranking      = mutate(result$feature_ranking,      dab = name)
    result$feature_val_results  = mutate(result$feature_val_results,  dab = name)
    result$final_test_results   = mutate(result$final_test_results,   dab = name)
    result$test_metrics         = mutate(result$test_metrics,         dab = name)

    result
  })
}


#' plot_perf_by_features
#'
#' Line plot of mean accuracy, kappa, and macro-F1 vs. feature count, with
#' ±1 SE error bars aggregated across all (seed × dAb) combinations.
#'
#' Identical function exists in rf_tidymodels.R and lasso_tidymodels.R.
#'
#' @param feature_val_results Stacked feature_val_results tibble from
#'   run_mars_fs_all_dabs() calls across multiple seeds. Must contain columns
#'   n_features, .metric, .estimate. A `dab` column is required when
#'   facet_by_dab = TRUE.
#' @param title        Optional plot title string
#' @param facet_by_dab If TRUE, produces a per-dAb faceted plot; error bars
#'   then reflect SE across seeds only
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
