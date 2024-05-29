library(tidyverse)
library(caret)
library(ggplot2)
library(stats)
library(psych)
library(glmnet)
library(gmodels)
library(gridExtra)
library(doParallel)
library(roxygen2)
library(RColorBrewer)
library(ggpubr)


setwd("..")
source("Visualization_and_preprocessing/S_protein_screen_preprocessing.R")
source("Embedded_FS/log_regression_accessory_functions.r")
source("Embedded_FS/lasso_resampling_functions.R")
split_dab_df_list = main_preprocessor()
split_dab_df_list = split_dab_df_list[-c(1, 4, 6)]


# Careful with this function name, dangerously close to glm.control from stats package
# In this case n_top_features is actually an array of features
glm_control_multi_fs = function(df,
                                target,
                                k1,
                                k2,
                                num_repeats,
                                threshold,
                                n_top_features = c()
) {
  cl = makeCluster(6)
  registerDoParallel(cl)
  
  external_perf_df = NULL
  global_perf_df = NULL
  
  external_conf_table_list = list()
  global_conf_table_list = list()
  
  external_opt_model_list = list()
  global_opt_model_list = list()
  
  average_coef_plot_list = list()
  
  threshold_df = df %>% filter(Conc. > threshold) %>%
    select(-c(Conc., Replicates))
  
  # Adding a dummy value at index 1 that gets replaced later when I know what 
  # the max number of features is
  num_feature_array = c(1, n_top_features)
  
  for (rep in 1:num_repeats) {
    global_kfold = createFolds(threshold_df[[target]], k = k1)
    
    for (fold1 in 1:k1) {
      global_train = threshold_df[-global_kfold[[fold1]], ]
      global_test = threshold_df[global_kfold[[fold1]], ]
      
      external_resample = resampling_function(
        target = target,
        dataset = global_train,
        n_folds = k2,
        perf_df = external_perf_df,
        model_list = external_opt_model_list,
        conf_table_list = external_conf_table_list,
        transformation = c("scale", "center", "BoxCox")
      )
      external_conf_table_list = external_resample$conf_table_list
      external_opt_model_list = external_resample$opt_model_list
      external_perf_df = external_resample$perf_df
      
      # SORTING COEFS AND DOWNSELECTING GLOBAL DF BASED ON EXTERNAL RESAMPLING
      #   Easier just to call this regardless since it's really fast      
      coefs_call = get_cv_coefs(external_opt_model_list, plot = T)
      average_coef_plot_list[[length(average_coef_plot_list) + 1]] = coefs_call$coef_plot
      average_external_coefs = coefs_call$coef_df
      
      average_abs_external_coefs = average_external_coefs %>% mutate_if(is.numeric, abs)
      
      sorted_coefs = average_abs_external_coefs %>% select_if(is.numeric) %>%
        mutate(sum = rowSums(.)) %>%
        mutate('Pair' = average_abs_external_coefs$Pair) %>%
        arrange(desc(sum))
      
      # Append the max number of features to the n_top_features list
      max_features = length(sorted_coefs$Pair) # the only other column is the target vector
      num_feature_array[1] = max_features
      # Remove any values from the array that are below the max feature number
      # If this occurs, then at the end of the function, a warning is printed out
      num_feature_array = num_feature_array[which(num_feature_array <= max_features)]
      
      global_resample = multi_feature_global_resampling(
        target = target,
        train = global_train,
        test = global_test,
        sorted_coefs = sorted_coefs,
        num_feature_list = num_feature_array,
        perf_df = global_perf_df,
        model_list = global_opt_model_list,
        conf_table_list = global_conf_table_list,
        transformation = c("scale", "center", "BoxCox")
      )
      global_conf_table_list = global_resample$conf_table_list
      global_opt_model_list = global_resample$opt_model_list
      global_perf_df = global_resample$perf_df
    }
  }
  
  # Averaging the performance based on the recorded metrics of the external
  # k-fold cv
  external_iterations = append(c(1:nrow(external_perf_df)),
                               c("Average", "Stdev"))
  global_iterations = append(c(1:nrow(global_perf_df)),
                             c("Average", "Stdev"))
  
  external_perf_df = bind_rows(
    external_perf_df,
    colMeans(external_perf_df, na.rm = T),
    sapply(external_perf_df, sd, na.rm = T)
  )
  external_perf_df = cbind(external_perf_df, external_iterations)
  
  global_perf_df = bind_rows(
    global_perf_df,
    colMeans(global_perf_df, na.rm = T),
    sapply(global_perf_df, sd, na.rm = T)
  )
  global_perf_df = cbind(global_perf_df, global_iterations)
  
  stopCluster(cl)
  registerDoSEQ()
  
  if(length(num_feature_array) < length(n_top_features)){
    cat("WARNING: \n The max number of features was below some of the feature numbers your specified")
  }
  
  return(list("external_conf_tables" = external_conf_table_list,
              "external_opt_models" = external_opt_model_list,
              "external_performance" = external_perf_df,
              "global_conf_tables" = global_conf_table_list,
              "global_opt_models" = global_opt_model_list,
              "global_performance" = global_perf_df,
              "sorted_coefs" = sorted_coefs$Pair
  ))
}

all_models = model_call_loop(df_list = split_dab_df_list,
                             model_type = 'lasso_fs',
                             k1 = 5,
                             k2 = 5,
                             num_repeats = 3,
                             threshold = 0,
                             n_top_features = c(10, 8, 7, 6, 5, 4, 3)
)


all_results = sapply(all_models, function(x) x['global_performance'])

sapply(names(all_results),
       function(x) write.table(all_results[[x]], file = paste(x,"txt",sep=".")) )
