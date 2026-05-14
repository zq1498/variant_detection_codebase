setwd("C:/Users/Zach/Documents/Duke/covspot_analysis/screen_analysis/variant_detection_codebase")
source("Visualization_and_preprocessing/S_protein_screen_preprocessing.R")

mars = readRDS("Embedded_FS/all_mars_results.rds")
pls  = readRDS("Embedded_FS/all_pls_results.rds")

for (dab in c("dAb18", "dAb19")) {
  cat("\n========", dab, "========\n")

  cat("\n-- MARS feature_val_results --\n")
  print(mars[[dab]]$feature_val_results |>
    tidyr::pivot_wider(names_from = .metric, values_from = .estimate) |>
    dplyr::arrange(n_features))

  cat("\n-- MARS test split size --\n")
  cat("train n:", nrow(rsample::training(mars[[dab]]$final_split)),
      "  test n:", nrow(rsample::testing(mars[[dab]]$final_split)), "\n")

  cat("\n-- PLS feature_val_results --\n")
  print(pls[[dab]]$feature_val_results |>
    tidyr::pivot_wider(names_from = .metric, values_from = .estimate) |>
    dplyr::arrange(n_features))

  cat("\n-- PLS test split size --\n")
  cat("train n:", nrow(rsample::training(pls[[dab]]$final_split)),
      "  test n:", nrow(rsample::testing(pls[[dab]]$final_split)), "\n")
}
