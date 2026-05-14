setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")

source("Visualization_and_preprocessing/S_protein_screen_preprocessing.R")
source("Embedded_FS/rf_tidymodels.R")
source("Embedded_FS/lasso_tidymodels.R")

split_dab_df_list = main_preprocessor()
split_dab_df_list = split_dab_df_list[-c(1, 4, 6)]

seeds = c(1969, 1986, 2024)

# ---- RF ---------------------------------------------------------------------
cat("====== RF: 3 seeds ======\n")
rf_runs = lapply(seeds, function(s) {
  cat("  seed", s, "\n")
  run_rf_fs_all_dabs(split_dab_df_list, seed = s)
})
saveRDS(rf_runs, "Embedded_FS/all_rf_runs_multiseed.rds")

rf_metrics = lapply(rf_runs, function(run) {
  map_dfr(run, \(x) x$test_metrics)
}) |> bind_rows()

rf_summary = rf_metrics |>
  filter(.metric %in% c("accuracy", "kap")) |>
  group_by(dab, .metric) |>
  summarise(mean = mean(.estimate),
            se   = sd(.estimate) / sqrt(n()),
            .groups = "drop") |>
  pivot_wider(names_from  = .metric,
              values_from = c(mean, se),
              names_glue  = "{.metric}_{.value}") |>
  arrange(desc(accuracy_mean))

cat("\nRF mean ± SE across 3 seeds:\n")
print(as.data.frame(rf_summary))

# ---- LASSO ------------------------------------------------------------------
cat("\n====== LASSO: 3 seeds ======\n")
lasso_runs = lapply(seeds, function(s) {
  cat("  seed", s, "\n")
  run_lasso_fs_all_dabs(split_dab_df_list, seed = s)
})
saveRDS(lasso_runs, "Embedded_FS/all_lasso_runs_multiseed.rds")

lasso_metrics = lapply(lasso_runs, function(run) {
  map_dfr(run, \(x) x$test_metrics)
}) |> bind_rows()

lasso_summary = lasso_metrics |>
  filter(.metric %in% c("accuracy", "kap")) |>
  group_by(dab, .metric) |>
  summarise(mean = mean(.estimate),
            se   = sd(.estimate) / sqrt(n()),
            .groups = "drop") |>
  pivot_wider(names_from  = .metric,
              values_from = c(mean, se),
              names_glue  = "{.metric}_{.value}") |>
  arrange(desc(accuracy_mean))

cat("\nLASSO mean ± SE across 3 seeds:\n")
print(as.data.frame(lasso_summary))

# ---- Plots ------------------------------------------------------------------
library(ggplot2)

rf_fvr    = map_dfr(rf_runs,    \(run) map_dfr(run, \(x) x$final_test_results))
lasso_fvr = map_dfr(lasso_runs, \(run) map_dfr(run, \(x) x$final_test_results))

# RF plots
p_rf_agg = plot_perf_by_features(
  rf_fvr,
  title        = "RF — Performance by Feature Subset (all dAbs, 3 seeds)",
  facet_by_dab = FALSE
)
ggsave("Embedded_FS/rf_perf_aggregate.png", p_rf_agg,
       width = 8, height = 5, dpi = 300)

p_rf_dab = plot_perf_by_features(
  rf_fvr,
  title        = "RF — Performance by Feature Subset (per dAb, 3 seeds)",
  facet_by_dab = TRUE
)
ggsave("Embedded_FS/rf_perf_by_dab.png", p_rf_dab,
       width = 14, height = 8, dpi = 300)

# LASSO plots
p_lasso_agg = plot_perf_by_features(
  lasso_fvr,
  title        = "LASSO — Performance by Feature Subset (all dAbs, 3 seeds)",
  facet_by_dab = FALSE
)
ggsave("Embedded_FS/lasso_perf_aggregate.png", p_lasso_agg,
       width = 8, height = 5, dpi = 300)

p_lasso_dab = plot_perf_by_features(
  lasso_fvr,
  title        = "LASSO — Performance by Feature Subset (per dAb, 3 seeds)",
  facet_by_dab = TRUE
)
ggsave("Embedded_FS/lasso_perf_by_dab.png", p_lasso_dab,
       width = 14, height = 8, dpi = 300)

cat("\nPlots saved:\n")
cat("  Embedded_FS/rf_perf_aggregate.png\n")
cat("  Embedded_FS/rf_perf_by_dab.png\n")
cat("  Embedded_FS/lasso_perf_aggregate.png\n")
cat("  Embedded_FS/lasso_perf_by_dab.png\n")
