setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")

library(tidyverse)
library(tidymodels)
library(ranger)
library(glmnet)
library(earth)
library(mRMRe)
library(Boruta)

source("Embedded_FS/lasso_tidymodels.R")
source("Embedded_FS/rf_tidymodels.R")
source("Embedded_FS/mars_tidymodels.R")
source("Embedded_FS/pls_tidymodels.R")
source("Embedded_FS/mrmr_tidymodels.R")
source("Embedded_FS/boruta_tidymodels.R")
source("Embedded_FS/composite_score.R")

cat("Loading multiseed results...\n")
lasso_runs  = readRDS("Embedded_FS/all_lasso_runs_multiseed.rds")
rf_runs     = readRDS("Embedded_FS/all_rf_runs_multiseed.rds")
mars_runs   = readRDS("Embedded_FS/all_mars_runs_multiseed.rds")
pls_runs    = readRDS("Embedded_FS/all_pls_runs_multiseed.rds")
mrmr_runs   = readRDS("Embedded_FS/all_mrmr_runs_multiseed.rds")
boruta_runs = readRDS("Embedded_FS/all_boruta_runs_multiseed.rds")

cat("Extracting feature rankings...\n")
extract_rankings = function(runs) map_dfr(runs, \(run) map_dfr(run, \(x) x$feature_ranking))

lasso_rankings  = extract_rankings(lasso_runs)
rf_rankings     = extract_rankings(rf_runs)
mars_rankings   = extract_rankings(mars_runs)
pls_rankings    = extract_rankings(pls_runs)
mrmr_rankings   = extract_rankings(mrmr_runs)
boruta_rankings = extract_rankings(boruta_runs)

cat("Running Borda aggregation (min_methods = 3)...\n")
composite = borda_aggregate(
  list(
    lasso  = lasso_rankings,
    rf     = rf_rankings,
    mars   = mars_rankings,
    pls    = pls_rankings,
    mrmr   = mrmr_rankings,
    boruta = boruta_rankings
  ),
  min_methods = 3
)

cat("\nTop 5 pairs per dAb:\n")
composite |>
  group_by(dab) |>
  slice_min(borda_score, n = 5, with_ties = FALSE) |>
  select(dab, feature, borda_score, composite_rank, n_methods) |>
  print(n = Inf)

cat("\nSaving plots...\n")
p_heat = plot_borda_heatmap(composite, top_n = 15)
ggsave("Embedded_FS/outputs/v2024/borda_composite_heatmap.png",
       p_heat, width = 12, height = 9, dpi = 300)
cat("  Saved: borda_composite_heatmap.png\n")

p_bars = plot_borda_bars(composite, top_n = 8)
ggsave("Embedded_FS/outputs/v2024/borda_composite_bars.png",
       p_bars, width = 14, height = 9, dpi = 300)
cat("  Saved: borda_composite_bars.png\n")

cat("\nDone.\n")
