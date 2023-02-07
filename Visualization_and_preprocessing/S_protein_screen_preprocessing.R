library(tidyverse)
library(psych)
library(sda)


# NOTES: DOSE RESPONSE IS DETERMINED BY SELECTING VARIANT CURVES AND TAKING ONLY THOSE WHERE THE MAX IS GREATER THAN 10,000 AND THE MIN IS LESS THAN 1000
# SIGNIFICANCE IS DETERMINED BY CAT SCORES, WHICH ARE CONVERTED TO P-VALUES USING D.O.F OF 2*NUM_SAMPLES-4 BASED ON THE D.O.F CALCULATION FOR DEPENDENT PAIRS
# ANY PAIRS WITH SIGNIFICANTLY ASSOCIATED WITH ANY VARIANT (P<0.1) ARE RETAINED.

# Operating under the condition that a good dose response generally has a low value of log fluorescence less than 3
# and a high value of log fluorescence above 4. Very rarely dose a good dose response have a base noise higher than 1000
# or a high value less than 10,000.

# To identify this, this function simply jerry rigs a "by column" tapply without using a for loop. Basically we call sapply
# so that each column vector is then fed into a tapply function which in turn calls the actual effector function. I use 
# tapply because its a convenient way of checking vectors by some external categorical variable (in our case the variant) Importantly,
# the tapply function does do some work in that it determines if any of the variants return true when passed to the dr check function.

dose_response_check = function(dataframe, cols_to_check, target_vector){
  # Accessory functions to dose_response_check
  dr_check2 = function(vec){
    bool = ifelse(max(vec) > 4 & min(vec) < 3, TRUE, FALSE)
    return(bool)
  }
  
  dr_check_tapply = function(col){
    tapply_index = dataframe[target_vector]
    dif_bool = tapply(col, tapply_index, dr_check2)
    if(any(dif_bool == TRUE)){
      return(T)
    }
    else{
      return(F)
    }
  }
  
  
  dr_cols = sapply(dataframe[cols_to_check], dr_check_tapply)
  dr_pairs = names(dr_cols[dr_cols == T])
  dr_df = dataframe[dr_pairs]
  ancillary_cols = dataframe %>% dplyr::select(-cols_to_check)
  dr_df = cbind(ancillary_cols, dr_df)
  return(dr_df)
}



# This function is completely deprecated. It was inappropriate to use CAT scores
# on my unsplit dataset, since CAT scores apply an LDA-like model to identify
# significant correlations.

cat_score_downselect = function(dataframe, cols_to_exclude, sig_level, target_vector){
  numeric_df = dataframe %>% dplyr::select(-cols_to_exclude)
  
  cat_matrix = catscore(Xtrain = as.matrix(numeric_df), L = dataframe[[target_vector]])
  
  deg_freedom = nrow(dataframe) * 2 - 4
  
  p_score_matrix = as.data.frame(2*pt(cat_matrix[,(1:4)], deg_freedom, lower.tail = F))
  
  sig_predictors = apply(p_score_matrix, 1, function(vec) ifelse(min(vec) <= sig_level, TRUE, FALSE))
  
  sig_pairs = names(sig_predictors[sig_predictors==T])
  significant_df = dataframe[sig_pairs]
  
  ancillary_cols = dataframe %>% dplyr::select(cols_to_exclude)
  sig_downselected_df = cbind(ancillary_cols, significant_df)
  
  return(sig_downselected_df)
}





# Ideally, we would like to use a single detection antibody in any final test
# Additionally, our screen experiments are constructed to test the performance
# of a detection antibody against every other capture antibody, so we can only explicitly
# claim efficacy in that context

# Therefore, we should split the dataset into N subsets where N=number of dAbs

split_dAb_dfs = function(dataframe_to_split, cols_to_exclude){
  # :param dataframe_to_split: Each feature name must follow the pattern of "dab_cab",
  # since the regex will fail otherwise: dataframe
  # :param cols_to_exclude: a list of column names that you don't want to include in the regex,
  # each will eventually be rebound to the split dataframes: list of strings
  # :return: a list where each element is a dataframe containing only the information
  # from a single dAb
  numeric_screen = dataframe_to_split %>% dplyr::select(-cols_to_exclude)
  pair_names = colnames(numeric_screen)
  excluded_cols = dataframe_to_split[cols_to_exclude]
  
  # For future me: regex patter just says find the underscore, and extract everything on the 
  # right side of it, then replace with an empty string
  dAb_names = unique(sub("_.*","",pair_names))

  dAb_dataframe_list = list()
  dataframe_list_ind = 1
  for (dAb in dAb_names){
    dAb_cols = numeric_screen[str_detect(pair_names, dAb)]
    split_dab_df = cbind(excluded_cols, dAb_cols)
    dAb_dataframe_list[[dataframe_list_ind]] = split_dab_df
    
    dataframe_list_ind = dataframe_list_ind + 1
  }
  
  return(dAb_dataframe_list)

}


main_preprocessor = function(){
  path = "data/covid_screen_data_cleaned.csv"
  screen_df = read.csv(path)
  screen_df$target_vector = as.factor(screen_df$target_vector)
  colnames(screen_df)[1] = "Conc."
  print(screen_df)
  
  cols_to_check = colnames(screen_df %>% dplyr::select(-Conc.) %>% select_if(is.numeric))
  
  preselected_screen_df = dose_response_check(screen_df, cols_to_check, "target_vector")
  
  cols_to_exclude = c("Conc.","Replicates","target_vector")
  
  # Removed because CAT scores are inappropriate
  # significant_df = cat_score_downselect(preselected_screen_df, cols_to_exclude, sig_level = 0.1, "target_vector")
  
  split_dabs = split_dAb_dfs(preselected_screen_df, cols_to_exclude)
  
  return(split_dabs)
}




