library(tidyverse)
library(caret)
library(recipes)


# Provided a dataset, this function can identify which values are outliers based
# on the distributions of values with the same concentration and variant. I then
# provides a data frame that contains relevant information about identified
# outliers that can be used to find and replace that outlier later on
check_outliers = function(dataset, alpha = 0.05) {
  flag_df = NULL
  current_variant = unique(dataset$target_vector)
  current_conc = unique(dataset$Conc.)
  
  cols_to_check = dataset %>% select(starts_with("cAb"))
  cab_names = colnames(cols_to_check)
  col_num = 1
  
  # Adapting function to replace grubbs.test. Basically, I need to set a p-value
  # to 0.01 if a value within the array is higher or lower than 3D*mean, then
  # return 'higher' or 'lower' depending on what value (probably need to choose the
  # value with the greatest distance from the mean so I don't get two returns)
  # for (array in cols_to_check){
  #   arr_mean = mean(array, na.rm=T)
  #   arr_sd = sd(array, na.rm = T)
  #   upper_bound = arr_mean + 2*arr_sd
  #   lower_bound = arr_mean - 2*arr_sd
  #   
  #   if (lower_bound < 0){
  #     lower_bound = 0
  #   }
  #   max_val = max(array)
  #   min_val = min(array)
  #   
  #   max_test = ifelse(max_val > upper_bound, 1, 0)
  #   min_test = ifelse(min_val < lower_bound, 1, 0)
  #   
  #   if(max_val < 0){
  #     max_test = 0
  #     min_test = 0
  #   }
  #   # print(paste('max_test: ',max_test, "\n", "min_test: ", min_test))
  #   # The original function requires a p_value as a check for whether to add to 
  #   # the flag df but we can check if either of the values are true
  #   if(max_test == 1 | min_test == 1){
  #     max_dist = max_val - arr_mean
  #     min_dist = arr_mean - min_val
  #     higher_lower = ifelse(abs(max_dist) > abs(min_dist), 'higher', 'lower')
  #     
  #     flagged_cAb = cab_names[[col_num]]
  #     flag_report = data.frame(
  #       variant = current_variant,
  #       concentration = current_conc,
  #       cAb = flagged_cAb,
  #       higher_lower = higher_lower)
  #       
  #     if (is.null(flag_df)) {
  #       flag_df = flag_report
  #     }
  #     else{
  #       flag_df = bind_rows(flag_df, flag_report)
  #     }
  #   }
  #   col_num = col_num + 1
  # }
  
  for (array in cols_to_check) {
    test = grubbs.test(array, type = 10)
    test_stat = test$statistic
    test_alpha = test$p.value
    
    if (test_alpha <= alpha & sum(array) != 0) {
      alternative_hypo = test$alternative
      higher_or_lower = ifelse(grepl('highest', alternative_hypo), 'higher', 'lower')
      
      flagged_cAb = cab_names[[col_num]]
      flag_report = data.frame(
        variant = current_variant,
        concentration = current_conc,
        cAb = flagged_cAb,
        higher_lower = higher_or_lower
      )
      
      if (is.null(flag_df)) {
        flag_df = flag_report
      }
      else{
        flag_df = bind_rows(flag_df, flag_report)
      }
    }
    col_num = col_num + 1
  }
  return(flag_df)
}

# Provided with a data frame containing information about individual outliers (called flag_df here), it will find those outlier values and add additional information to the original flag_df, mainly the specific outlier value and the mean of the group without the outlier included
enhance_flag_df = function(flag_df, dataset){
  # Loop over all rows of flag df, then build an array row by row that contains
  # the max and mean (na.rm = T), then add that to the original data frame and 
  # return it
  num_rows = nrow(flag_df)
  enhanced_df = NULL
  for (iter in 1:num_rows){
    current_outlier_row = flag_df[iter,]
    
    current_outlier_array = dataset %>% filter(target_vector == current_outlier_row$variant,
                                               Conc. == current_outlier_row$concentration) %>% 
      select(current_outlier_row$cAb)
    
    outlier_cab = current_outlier_row$cAb
    
    outlier_col = current_outlier_array[current_outlier_row$cAb]
    
    if (current_outlier_row$higher_lower == 'higher'){
      # Excluding the maximum value, obtain the mean of the column
      outlier_value = max(outlier_col)
      temp_col = outlier_col[[outlier_cab]]
      arr_without_max = temp_col[temp_col <  outlier_value]
      mean_val = mean(arr_without_max, na.rm = T)
    }
    
    else{
      # Excluding the minimum value, obtain the mean of the column
      outlier_value = min(outlier_col)
      temp_col = outlier_col[[outlier_cab]]
      arr_without_min = temp_col[temp_col > outlier_value]
      mean_val = mean(arr_without_min, na.rm = T)
    }
    
    extra_info=data.frame(outlier_value = outlier_value,
                          mean_value = mean_val)
    extra_info$mean_value = ifelse(is.nan(extra_info$mean_value), 0, extra_info$mean_value)
    
    extra_info_df = bind_cols(current_outlier_row, extra_info)
    
    if (is.null(enhanced_df)){
      enhanced_df = extra_info_df
    }
    else{
      enhanced_df = bind_rows(enhanced_df, extra_info_df)
    }
  }
  return(enhanced_df)
}

# Given an enhanced flag_df and the original dataset, this function will find
# all outlier values and replace them with the mean value or NA if you'd like to
# use an alternative imputation method. It returns a dataset of the same
# dimensions as the original, with the outlier values changed.
replace_outliers = function(dataset, enhanced_flag_df, impute_method='mean'){
  row_numbers = 1:nrow(dataset)
  dataset = bind_cols("row_number" = row_numbers, dataset)
  matching_row_arr = list()
  
  for(i in 1:nrow(enhanced_flag_df)) {
    current_conc = enhanced_flag_df$concentration[i]
    current_variant = enhanced_flag_df$variant[i]
    outlier_col = enhanced_flag_df$cAb[i]
    outlier_val = enhanced_flag_df$outlier_value[i]
    # Replacement value can either be the mean value contained in enhanced_flag_df
    # or it can just replace with NA
    replacement_value = ifelse(impute_method=='mean', enhanced_flag_df$mean_value[i], NA)
    
    # Locate the subset of rows that contain the outlier value
    matching_subset = dataset %>% 
      filter(Conc. == current_conc, target_vector == current_variant)
    
    # Find the single row within the matching subset where the outlier col has
    # the outlier val, then take the first column (the row number) and use that
    # to set the value to the mean in the original dataset
    matching_row = as.tibble(matching_subset[matching_subset[outlier_col]==outlier_val])
    matching_row_num = as.integer(matching_subset[matching_subset[outlier_col]==outlier_val][1])
    
    matching_row_arr = c(matching_row_arr, matching_row_num)

    # Replace outlier value with replacement value using the row number found above
    dataset[matching_row_num, outlier_col] = replacement_value
  }
  return(list('dataset' = dataset,
              'row_nums_arr' = matching_row_arr))
}


#------------------------------------------------------------------------------#
# Outlier control function
main_outliers_detection_function = function(dataset, threshold, impute = 'mean') {
  outlier_test_df_list = dataset %>% filter(Conc. > threshold) %>% group_by(Conc., target_vector) %>% group_split()
  flag_df = NULL
  
  for (df in outlier_test_df_list) {
    outlier_check = check_outliers(df)
    
    if (is.null(flag_df)) {
      flag_df = outlier_check
    }
    
    else{
      flag_df = bind_rows(flag_df, outlier_check)
    }
  }
  # Add specific outlier value and mean value of outlier group without the outlier
  # value to the flag_df
  print(flag_df)
  enhanced_flag_df = enhance_flag_df(flag_df, dataset)
  
  # Now that we've returned a more robust version of the flag_df, we can use that
  # information to iterate over the entire data set and replace outlier values 
  # with either the mean value of their group or with NA (for MICE)
  
  outlier_removed_df = replace_outliers(dataset, enhanced_flag_df, impute_method=impute)
  outlier_removed_df = outlier_removed_df$dataset %>% select(-row_number)
  return(outlier_removed_df)
}

imputation_function = function(test_set, training_set, impute_method = 'knn'){
  # Default imputation method is knn but can also do bag and linear imputation
  impute_rec = recipe(
    paste("target_vector","~",".") %>% as.formula(),
    data=training_set
  ) %>% step_rm(removals='Conc.')
  
  if (impute_method == 'knn'){
    imputed_train = impute_rec %>% 
      step_impute_knn(all_predictors(), neighbors = 3) %>% 
      prep(training=training_set) %>% 
      bake(training_set)
    
    imputed_test = impute_rec %>% 
      step_impute_knn(all_predictors(), neightbors = 3) %>% 
      prep(training=training_set) %>% 
      bake(test_set)
  }
  
  if (impute_method == 'bag')
    impute_train = impute_rec %>% 
      step_impute_bag(all_predictors(), trees = 25) %>%
      prep(training=training_set) %>% 
      bake(train)
    
    impute_test = impute_rec %>% 
      step_impute_bag(all_predictors(), trees = 25) %>% 
      prep(training=training_set) %>% 
      bake(test_set)
    
    if (impute_method == 'linear')
      impute_train = impute_rec %>% 
      step_impute_linear(all_predictors()) %>% 
      prep(training=training_set) %>% 
      bake(train)
    
      impute_test = impute_rec %>% 
        step_impute_linear(all_predictors()) %>% 
        prep(training=train) %>% 
        bake(test_set)
      
    return(list('impute_train' = impute_train,
                'impute_test' = impute_test))
}