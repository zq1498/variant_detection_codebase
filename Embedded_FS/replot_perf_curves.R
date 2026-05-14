setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")

library(tidyverse)
library(tidymodels)
library(ggplot2)
source("Embedded_FS/rf_tidymodels.R")
source("Embedded_FS/lasso_tidymodels.R")
source("Embedded_FS/mars_tidymodels.R")
source("Embedded_FS/pls_tidymodels.R")

rf_runs    = readRDS("Embedded_FS/all_rf_runs_multiseed.rds")
lasso_runs = readRDS("Embedded_FS/all_lasso_runs_multiseed.rds")
mars_runs  = readRDS("Embedded_FS/all_mars_runs_multiseed.rds")
pls_runs   = readRDS("Embedded_FS/all_pls_runs_multiseed.rds")

rf_fvr    = map_dfr(rf_runs,    \(run) map_dfr(run, \(x) x$final_test_results))
lasso_fvr = map_dfr(lasso_runs, \(run) map_dfr(run, \(x) x$final_test_results))
mars_fvr  = map_dfr(mars_runs,  \(run) map_dfr(run, \(x) x$final_test_results))
pls_fvr   = map_dfr(pls_runs,   \(run) map_dfr(run, \(x) x$final_test_results))

save_plot = function(fvr, prefix, model_label) {
  p_agg = plot_perf_by_features(fvr,
    title        = paste(model_label, "— Performance by Feature Subset (all dAbs, 3 seeds)"),
    facet_by_dab = FALSE)
  ggsave(paste0("Embedded_FS/", prefix, "_perf_aggregate.png"), p_agg,
         width = 8, height = 5, dpi = 300)

  p_dab = plot_perf_by_features(fvr,
    title        = paste(model_label, "— Performance by Feature Subset (per dAb, 3 seeds)"),
    facet_by_dab = TRUE)
  ggsave(paste0("Embedded_FS/", prefix, "_perf_by_dab.png"), p_dab,
         width = 14, height = 8, dpi = 300)

  cat("Saved:", prefix, "aggregate and per-dAb plots\n")
}

save_plot(rf_fvr,    "rf",    "RF")
save_plot(lasso_fvr, "lasso", "LASSO")
save_plot(mars_fvr,  "mars",  "MARS")
save_plot(pls_fvr,   "pls",   "PLS-DA")
