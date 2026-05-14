setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")

source("Visualization_and_preprocessing/S_protein_screen_preprocessing.R")
source("Embedded_FS/mars_tidymodels.R")
source("Embedded_FS/pls_tidymodels.R")

split_dab_df_list = main_preprocessor()
split_dab_df_list = split_dab_df_list[-c(1, 4, 6)]

seeds = c(1969, 1986, 2024)

# ---- MARS -------------------------------------------------------------------
cat("====== MARS: 3 seeds ======\n")
mars_runs = lapply(seeds, function(s) {
  cat("  seed", s, "\n")
  run_mars_fs_all_dabs(split_dab_df_list, seed = s)
})
saveRDS(mars_runs, "Embedded_FS/all_mars_runs_multiseed.rds")

mars_metrics = lapply(mars_runs, function(run) {
  map_dfr(run, \(x) x$test_metrics)
}) |> bind_rows()

mars_summary = mars_metrics |>
  filter(.metric %in% c("accuracy", "kap")) |>
  group_by(dab, .metric) |>
  summarise(mean = mean(.estimate),
            se   = sd(.estimate) / sqrt(n()),
            .groups = "drop") |>
  pivot_wider(names_from = .metric,
              values_from = c(mean, se),
              names_glue  = "{.metric}_{.value}") |>
  arrange(desc(accuracy_mean))

cat("\nMARS mean ± SE across 3 seeds:\n")
print(as.data.frame(mars_summary))

# ---- PLS-DA -----------------------------------------------------------------
cat("\n====== PLS-DA: 3 seeds ======\n")
pls_runs = lapply(seeds, function(s) {
  cat("  seed", s, "\n")
  run_pls_fs_all_dabs(split_dab_df_list, seed = s)
})
saveRDS(pls_runs, "Embedded_FS/all_pls_runs_multiseed.rds")

pls_metrics = lapply(pls_runs, function(run) {
  map_dfr(run, \(x) x$test_metrics)
}) |> bind_rows()

pls_summary = pls_metrics |>
  filter(.metric %in% c("accuracy", "kap")) |>
  group_by(dab, .metric) |>
  summarise(mean = mean(.estimate),
            se   = sd(.estimate) / sqrt(n()),
            .groups = "drop") |>
  pivot_wider(names_from = .metric,
              values_from = c(mean, se),
              names_glue  = "{.metric}_{.value}") |>
  arrange(desc(accuracy_mean))

cat("\nPLS-DA mean ± SE across 3 seeds:\n")
print(as.data.frame(pls_summary))
