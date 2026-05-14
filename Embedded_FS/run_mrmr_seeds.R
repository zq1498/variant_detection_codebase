setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")

source("Visualization_and_preprocessing/S_protein_screen_preprocessing.R")
source("Embedded_FS/mrmr_tidymodels.R")

split_dab_df_list = main_preprocessor()
split_dab_df_list = split_dab_df_list[-c(1, 4, 6)]

seeds = c(1969, 1986, 2024)

# ---- mRMR -------------------------------------------------------------------
cat("====== mRMR: 3 seeds ======\n")
mrmr_runs = lapply(seeds, function(s) {
  cat("  seed", s, "\n")
  run_mrmr_fs_all_dabs(split_dab_df_list, seed = s)
})
saveRDS(mrmr_runs, "Embedded_FS/all_mrmr_runs_multiseed.rds")

mrmr_metrics = lapply(mrmr_runs, function(run) {
  map_dfr(run, \(x) x$test_metrics)
}) |> bind_rows()

mrmr_summary = mrmr_metrics |>
  filter(.metric %in% c("accuracy", "kap")) |>
  group_by(dab, .metric) |>
  summarise(mean = mean(.estimate),
            se   = sd(.estimate) / sqrt(n()),
            .groups = "drop") |>
  pivot_wider(names_from  = .metric,
              values_from = c(mean, se),
              names_glue  = "{.metric}_{.value}") |>
  arrange(desc(accuracy_mean))

cat("\nmRMR mean ± SE across 3 seeds:\n")
print(as.data.frame(mrmr_summary))

# ---- Plots ------------------------------------------------------------------
library(ggplot2)

mrmr_fvr = map_dfr(mrmr_runs, \(run) map_dfr(run, \(x) x$final_test_results))

p_mrmr_agg = plot_perf_by_features(
  mrmr_fvr,
  title        = "mRMR — Performance by Feature Subset (all dAbs, 3 seeds)",
  facet_by_dab = FALSE
)
ggsave("Embedded_FS/mrmr_perf_aggregate.png", p_mrmr_agg,
       width = 8, height = 5, dpi = 300)

p_mrmr_dab = plot_perf_by_features(
  mrmr_fvr,
  title        = "mRMR — Performance by Feature Subset (per dAb, 3 seeds)",
  facet_by_dab = TRUE
)
ggsave("Embedded_FS/mrmr_perf_by_dab.png", p_mrmr_dab,
       width = 14, height = 8, dpi = 300)

cat("\nPlots saved:\n")
cat("  Embedded_FS/mrmr_perf_aggregate.png\n")
cat("  Embedded_FS/mrmr_perf_by_dab.png\n")
