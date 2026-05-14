library(tidyverse)

# TODO: SES integration
# SES (SES_FS.Rmd) does not yet produce a feature_ranking tibble in the standard
# format. When ready, the most natural importance score is -log10(ses_pvalues),
# which is available for all predictors (not just selected ones), giving a
# continuous importance score compatible with the other methods. Unselected
# features will naturally fall to the bottom of the ranking via their high
# p-values. Draft entry for the method_rankings list:
#
#   ses_rankings = map_dfr(all_ses_results, \(x) {
#     tibble(
#       dab       = x$dab,          # needs to be added in run_ses_fs_all_dabs
#       seed      = x$seed,
#       feature   = names(x$ses_pvalues),
#       importance = -log10(x$ses_pvalues + 1e-300)  # guard against p = 0
#     )
#   })

# Normalize one method's feature_ranking tibble to (dab, feature, importance, method)
# where importance is always higher-is-better.
# If a `seed` column is present, averages across seeds before ranking.
# Boruta: Rejected features are excluded (majority-vote across seeds).
# mRMR: mrmr_rank (ascending) is negated to match the higher-is-better convention.
standardize_ranking = function(ranking_df, method_name) {
  df           = ranking_df
  has_seed     = "seed"      %in% colnames(df)
  has_decision = "decision"  %in% colnames(df)
  has_mrmr     = "mrmr_rank" %in% colnames(df)

  if (has_seed) {
    if (has_decision) {
      df = df |>
        group_by(dab, feature) |>
        summarise(
          importance = mean(importance),
          decision   = names(which.max(table(decision))),
          .groups    = "drop"
        )
    } else if (has_mrmr) {
      df = df |>
        group_by(dab, feature) |>
        summarise(mrmr_rank = mean(mrmr_rank), .groups = "drop")
    } else {
      df = df |>
        group_by(dab, feature) |>
        summarise(importance = mean(importance), .groups = "drop")
    }
  }

  if (has_decision) {
    df = df |> filter(decision != "Rejected") |> select(-decision)
  }

  if (has_mrmr) {
    df = df |>
      mutate(importance = -mrmr_rank) |>
      select(-mrmr_rank, -any_of("mrmr_score"))
  }

  df |> select(dab, feature, importance) |> mutate(method = method_name)
}


# Assign Borda rank within each (dab, method).
# Rank 1 = highest importance. Ties share the average rank.
assign_borda_rank = function(df) {
  df |>
    group_by(dab, method) |>
    mutate(borda_rank = rank(-importance, ties.method = "average")) |>
    ungroup()
}


#' borda_aggregate
#'
#' Aggregate multiple feature selection rankings into a single Borda count
#' composite score. Within each dAb, each method's features are ranked by
#' importance (rank 1 = most important). Ranks are summed across methods so
#' a lower Borda score means consistently higher priority across methods.
#' Features absent from a method receive a penalty rank equal to that
#' method's feature count for that dAb plus 1.
#'
#' @param method_rankings Named list of feature_ranking tibbles, one per
#'   method. Names are used as method labels in the output.
#'   Each tibble must contain: dab, feature, and either
#'     - `importance`  (higher = better; LASSO, RF, MARS, PLS, Boruta)
#'     - `mrmr_rank`   (lower = better; mRMR)
#'   A `seed` column is allowed and averaged over before ranking.
#'   A `decision` column (Boruta) causes Rejected features to be excluded.
#'
#' @param min_methods Minimum number of methods in which a pair must appear
#'   to be included in the output. Pairs absent from a method still receive
#'   a penalty rank but count towards this threshold only if they were
#'   observed (non-absent) in at least min_methods methods. Default 1.
#'
#' @return Tibble with columns: dab, feature, borda_score (sum of per-method
#'   ranks, lower = higher priority), n_methods (number of methods where the
#'   pair was observed), one rank_<method> column per method, and
#'   composite_rank (within-dab ordinal rank sorted by borda_score).
borda_aggregate = function(method_rankings, min_methods = 1) {
  standardized = imap_dfr(method_rankings, standardize_ranking)
  ranked       = assign_borda_rank(standardized)

  all_pairs   = ranked |> distinct(dab, feature)
  all_methods = ranked |> distinct(method)

  method_feature_counts = ranked |>
    count(dab, method, name = "n_features")

  full_grid = tidyr::crossing(all_pairs, all_methods) |>
    left_join(
      ranked |> select(dab, feature, method, borda_rank),
      by = c("dab", "feature", "method")
    ) |>
    left_join(method_feature_counts, by = c("dab", "method")) |>
    mutate(
      borda_rank = if_else(is.na(borda_rank), as.double(n_features) + 1, borda_rank)
    ) |>
    select(-n_features)

  rank_wide = full_grid |>
    pivot_wider(names_from = method, values_from = borda_rank, names_prefix = "rank_")

  n_observed = ranked |> count(dab, feature, name = "n_methods")

  rank_wide |>
    left_join(n_observed, by = c("dab", "feature")) |>
    filter(n_methods >= min_methods) |>
    mutate(borda_score = rowSums(across(starts_with("rank_")), na.rm = TRUE)) |>
    group_by(dab) |>
    mutate(composite_rank = rank(borda_score, ties.method = "average")) |>
    ungroup() |>
    arrange(dab, borda_score)
}


#' plot_borda_heatmap
#'
#' Tile heatmap of scaled Borda scores (importance) across dAbs and pairs.
#' Within each dAb, scores are rescaled 0–1 so colour reflects relative
#' priority rather than absolute score.
#'
#' @param composite Output tibble from borda_aggregate().
#' @param top_n     Show only the top N pairs per dAb. NULL shows all.
plot_borda_heatmap = function(composite, top_n = NULL) {
  df = composite

  if (!is.null(top_n)) {
    df = df |>
      group_by(dab) |>
      slice_min(borda_score, n = top_n, with_ties = FALSE) |>
      ungroup()
  }

  df = df |>
    group_by(dab) |>
    mutate(
      scaled_importance = 1 - (borda_score - min(borda_score)) /
        pmax(max(borda_score) - min(borda_score), 1e-9)
    ) |>
    ungroup() |>
    # Strip dAb prefix from feature name so shared cAbs appear as one row
    mutate(cab = sub("^dAb[^_]+_", "", feature))

  # Order cAb rows by mean scaled importance across all dAbs
  cab_order = df |>
    group_by(cab) |>
    summarise(mean_imp = mean(scaled_importance), .groups = "drop") |>
    arrange(mean_imp) |>
    pull(cab)

  df = df |> mutate(cab = factor(cab, levels = cab_order))

  ggplot(df, aes(x = dab, y = cab, fill = scaled_importance)) +
    geom_tile(color = "white") +
    scale_fill_gradient(low = "#FFFFFF", high = "#012169", name = "Scaled Rank") +
    labs(
      title    = "Antibody Importance Across Feature Selection Methods",
      subtitle = "Scaled Borda rank of capture and detection pairs",
      x = "Detection Antibody",
      y = "Capture Antibody"
    ) +
    theme_classic() +
    theme(
      axis.text.x  = element_text(angle = 45, vjust = 1, hjust = 1, size = 12),
      axis.text.y  = element_text(size = 10),
      axis.title   = element_text(size = 14),
      plot.title   = element_text(size = 14, face = "bold"),
      legend.title = element_text(size = 11)
    )
}


#' plot_borda_bars
#'
#' Bar charts of top N pairs by Borda composite score, faceted by dAb.
#' Lower scores are better; bars are plotted in ascending score order within
#' each facet so the best pair appears at the top.
#'
#' @param composite Output from borda_aggregate().
#' @param top_n     Number of top pairs to show per dAb. Default 10.
plot_borda_bars = function(composite, top_n = 10) {
  df = composite |>
    group_by(dab) |>
    slice_min(borda_score, n = top_n, with_ties = FALSE) |>
    arrange(dab, borda_score) |>
    mutate(
      feat_key = paste0(dab, "___", feature),
      feat_key = factor(feat_key, levels = feat_key)
    ) |>
    ungroup()

  ggplot(df, aes(x = feat_key, y = borda_score)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    scale_x_discrete(labels = \(x) sub(".*___", "", x)) +
    facet_wrap(~dab, scales = "free", ncol = 4) +
    labs(
      title = paste0("Top ", top_n, " pairs by Borda composite rank, per dAb"),
      x     = NULL,
      y     = "Borda score (lower = higher priority)"
    ) +
    theme_bw() +
    theme(
      strip.text = element_text(size = 11),
      axis.text  = element_text(size = 9)
    )
}


# =============================================================================
# Usage example
# =============================================================================
#
# source("Embedded_FS/lasso_tidymodels.R")
# source("Embedded_FS/rf_tidymodels.R")
# source("Embedded_FS/mars_tidymodels.R")
# source("Embedded_FS/pls_tidymodels.R")
# source("Embedded_FS/mrmr_tidymodels.R")
# source("Embedded_FS/boruta_tidymodels.R")
# source("Embedded_FS/composite_score.R")
#
# Extract feature rankings from each method's all-dabs run:
#   lasso_rankings  = map_dfr(all_lasso_results,  \(x) x$feature_ranking)
#   rf_rankings     = map_dfr(all_rf_results,     \(x) x$feature_ranking)
#   mars_rankings   = map_dfr(all_mars_results,   \(x) x$feature_ranking)
#   pls_rankings    = map_dfr(all_pls_results,    \(x) x$feature_ranking)
#   mrmr_rankings   = map_dfr(all_mrmr_results,   \(x) x$feature_ranking)
#   boruta_rankings = map_dfr(all_boruta_results, \(x) x$feature_ranking)
#
# Aggregate — require a pair to appear in at least 2 methods:
#   composite = borda_aggregate(
#     list(
#       lasso  = lasso_rankings,
#       rf     = rf_rankings,
#       mars   = mars_rankings,
#       pls    = pls_rankings,
#       mrmr   = mrmr_rankings,
#       boruta = boruta_rankings
#     ),
#     min_methods = 2
#   )
#
# Visualise:
#   plot_borda_heatmap(composite, top_n = 15)
#   plot_borda_bars(composite, top_n = 10)
#
# Save:
#   write.xlsx(composite, "pair_subset_results/borda_composite.xlsx")
