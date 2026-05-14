setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")

library(tidyverse)
library(tidymodels)
source("Embedded_FS/mars_tidymodels.R")
source("Embedded_FS/pls_tidymodels.R")

mars_runs = readRDS("Embedded_FS/all_mars_runs_multiseed.rds")
pls_runs  = readRDS("Embedded_FS/all_pls_runs_multiseed.rds")

# Stack final_test_results across all seeds and dAbs.
# These are retrained on all TRAINING data and scored on FINAL_TEST at each
# feature count, so the curves are consistent with the summary table (which
# shows the same evaluation at the selected best N).
mars_fvr = map_dfr(mars_runs, \(run) map_dfr(run, \(x) x$final_test_results))
pls_fvr  = map_dfr(pls_runs,  \(run) map_dfr(run, \(x) x$final_test_results))

# ---- MARS plots -------------------------------------------------------------
p_mars_agg = plot_perf_by_features(
  mars_fvr,
  title        = "MARS — Performance by Feature Subset (all dAbs, 3 seeds)",
  facet_by_dab = FALSE
)
ggsave("Embedded_FS/mars_perf_aggregate.png", p_mars_agg,
       width = 8, height = 5, dpi = 300)

p_mars_dab = plot_perf_by_features(
  mars_fvr,
  title        = "MARS — Performance by Feature Subset (per dAb, 3 seeds)",
  facet_by_dab = TRUE
)
ggsave("Embedded_FS/mars_perf_by_dab.png", p_mars_dab,
       width = 14, height = 8, dpi = 300)

# ---- PLS-DA plots -----------------------------------------------------------
p_pls_agg = plot_perf_by_features(
  pls_fvr,
  title        = "PLS-DA — Performance by Feature Subset (all dAbs, 3 seeds)",
  facet_by_dab = FALSE
)
ggsave("Embedded_FS/pls_perf_aggregate.png", p_pls_agg,
       width = 8, height = 5, dpi = 300)

p_pls_dab = plot_perf_by_features(
  pls_fvr,
  title        = "PLS-DA — Performance by Feature Subset (per dAb, 3 seeds)",
  facet_by_dab = TRUE
)
ggsave("Embedded_FS/pls_perf_by_dab.png", p_pls_dab,
       width = 14, height = 8, dpi = 300)

cat("Plots saved:\n")
cat("  Embedded_FS/mars_perf_aggregate.png\n")
cat("  Embedded_FS/mars_perf_by_dab.png\n")
cat("  Embedded_FS/pls_perf_aggregate.png\n")
cat("  Embedded_FS/pls_perf_by_dab.png\n")
