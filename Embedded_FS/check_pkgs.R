pkgs <- c("psych", "sda", "tidymodels", "discrim", "earth", "plsmod", "mixOmics",
          "mda",      # discrim_flexible backend
          "hardhat",  # tune() sentinel
          "ranger",   # used in rf_tidymodels
          "glmnet"    # used in lasso_tidymodels
          )
for (p in pkgs) {
  status <- system.file(package = p)
  cat(p, if (nchar(status) > 0) "OK" else "MISSING", "\n")
}
