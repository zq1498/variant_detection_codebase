setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")

source("Visualization_and_preprocessing/S_protein_screen_preprocessing.R")
source("Embedded_FS/mars_tidymodels.R")
source("Embedded_FS/pls_tidymodels.R")

split_dab_df_list = main_preprocessor()
split_dab_df_list = split_dab_df_list[-c(1, 4, 6)]

# ---- MARS all dAbs ----------------------------------------------------------
cat("\n====== MARS: all dAbs ======\n")
all_mars_results = run_mars_fs_all_dabs(split_dab_df_list, seed = 2024)
saveRDS(all_mars_results, "Embedded_FS/all_mars_results.rds")

mars_summary = imap_dfr(all_mars_results, function(res, name) {
  bind_cols(tibble(dab = name,
                   best_num_terms   = res$best_num_terms,
                   best_prod_degree = res$best_prod_degree,
                   best_n_features  = res$best_n_features,
                   best_features    = paste(res$best_features, collapse = ", ")),
            res$test_metrics |> select(.metric, .estimate) |>
              pivot_wider(names_from = .metric, values_from = .estimate))
})
cat("\nMARS results summary:\n")
print(as.data.frame(mars_summary))

# ---- PLS all dAbs -----------------------------------------------------------
cat("\n====== PLS-DA: all dAbs ======\n")
all_pls_results = run_pls_fs_all_dabs(split_dab_df_list, seed = 2024)
saveRDS(all_pls_results, "Embedded_FS/all_pls_results.rds")

pls_summary = imap_dfr(all_pls_results, function(res, name) {
  bind_cols(tibble(dab           = name,
                   best_num_comp = res$best_num_comp,
                   best_n_features = res$best_n_features,
                   best_features   = paste(res$best_features, collapse = ", ")),
            res$test_metrics |> select(.metric, .estimate) |>
              pivot_wider(names_from = .metric, values_from = .estimate))
})
cat("\nPLS-DA results summary:\n")
print(as.data.frame(pls_summary))
