# Imports -- straight copied from log regression parent script, should trim unnecessary imports
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


# Contains all functions necessary for handling the output of the log regression
# control functions

# While the control functions actually train and cross-validate the models,
# these accessory functions are necessary to actually parse the outputs and
# obtain informative results



#' @title get_cv_coefs
#' @usage: Takes a list of CV.glmnet objects, obtains the coefficients from every 
#' single one, then averages them. This is used when you want to take a set of 
#' cross-validation runs and visualize the average feature importance (by beta
#' coefs) obtained over that cross-validation
#'
#' @param model_list: list of glmnet models of type multinom  
#' @param condense: (default T) determines whether or not the zero beta rows 
#' should be removed from the final matrix 
#' @param plot: (default F) Determines whether you want Rstudio to plot the output 
#'
#' @return: a list with two elements, one is a data frame of the averaged coefs
#' and the other is a list of plots of these coefs used for visualization 
get_cv_coefs = function(model_list,
                        condense = T,
                        plot = F) {
  
  coef_df = NULL
  for (model in model_list) {
    betas = coef(model, s = 'lambda.1se')
    coefs = data.frame(Pair = rownames(betas[[1]]))
    for (num in 1:4) {
      coefs = cbind(coefs, betas[[num]][, 1])
      colnames(coefs)[num + 1] = rownames(as.matrix(betas[num]))
    }
    if (condense == T) {
      coefs = filter(coefs, coefs[2] != 0)
    }
    if (is.null(coef_df)) {
      coef_df = coefs
    }
    else{
      coef_df = rbind(coef_df, coefs)
    }
  }
  coef_df = coef_df %>% group_by(Pair) %>% summarize_all(mean)
  
  if (plot == T) {
    coef_plot = coef_plot(coef_df)
  }
  return(list("coef_df" = coef_df,
              "coef_plot" = coef_plot))
}


#' @title coef plot
#' @usage Takes a data frame of coefficients (such as the one created by 
#' get_cv_coefs) and returns the plots of those coefficients. Coupled to 
#' get_cv_coefs in that it is called to create the coef plots that get_cv_coefs
#' returns
#' @seealso get_cv_coefs
#' 
#' @param coef_df: a data frame of beta coefficients
#' 
#' @return: a list of plots that visualize the beta coefficients
coef_plot = function(coef_df) {
  
  coef_df = coef_df %>% mutate_if(is.numeric, abs)
  coef_plot = coef_df %>% pivot_longer(col = -Pair,
                                       names_to = 'variant') %>%
    ggplot(aes(
      x = reorder(Pair,-value),
      y = value,
      fill = variant
    )) +
    geom_bar(position = 'dodge', stat = 'identity') +
    labs(# title = 'Average (absolute value) lasso coefficients by variant',
         x = "Pairs",
         y = 'Coefficient') +
    theme_classic() + 
    theme(axis.text.x = element_text(
      angle = 45,
      vjust = 1,
      hjust = 1
    )) +
    scale_fill_brewer(palette = 'Dark2')
  return(coef_plot)
}


#' @title extract_performance_metrics
#' @usage Used by the multi_performance_viz function to obtain performance metrics
#' from each element within the list. Currently is only compatible with specific
#' template of performance metric output (rows: iteration, accuracy, kappa, delta_F1, 
#' omi_F1, wt_F1, beta_F1, average_F1)
#' @seealso multi_performance_viz
#'   
#' @param list_of_results: A list of model results as produced by any of the main
#' glm control functions
#' @param ordered_list_of_dabs: A list of strings wherein each element is a dab
#' that I want I'm going to obtain performance metrics for
#' @param object_to_obtain: (default = 'performance_metrics) used to index the 
#' list of results (in which each element is actually a list of results in itself)
#' and obtain specifically obtain the intended object
#'
#' @return: A data frame of performance metrics for each of the dabs in the
#' provided ordered list of dabs
#' @note Really this function is a just a subfunction called by some of the other
#' accessory function. To avoid further abstraction, DO NOT DECOMPOSE THIS FUNCTION.
#' 
#' @note I desperately want to remove the ordered list of dabs and replace with
#' a simpler input, but I really don't have a plan for that
extract_performance_metrics = function(list_of_results,
                                       ordered_list_of_dabs,
                                       object_to_obtain = 'performance_metrics') {
  # Getting averages
  # Getting Kappa and accuracy
  performance_metric_list = sapply(list_of_results, function (x) x[object_to_obtain])
  averages = sapply(performance_metric_list, function (x) x %>% filter(iteration == 'Averages'))
  
  Accuracy = unlist(averages[row.names(averages) == 'accuracy'])
  Kappa = unlist(averages[row.names(averages) == 'kappa'])
  
  # Produces a list (# dabs * 4 in length) that contains the values in the following arrangement:
  #     dab1_delta, dab1_omi, dab1_wt, dab1_beta, dab2_delta, dab2_omi, dab2_wt dab2_beta, ...
  #     the order is dependent on the rownames ordering in the matrix, NOT the ordering in the list below
  f1_metrics = c("Delta_F1", "WT_F1", "Omi_F1", "Beta_F1")
  f1s = unlist(averages[rownames(averages) %in% f1_metrics, ])
  
  # This will be the final dataframe that contains all the values
  dabwise_performance = data.frame("dAb" = ordered_list_of_dabs)
  
  # Eventually will cbind this dataframe to the dabwise performance
  temp_f1_df = data.frame(matrix(nrow = 1, ncol = 4))
  colnames(temp_f1_df) = rownames(averages)[rownames(averages) %>% str_detect("F1")]
  
  # Obtain the first 4 values, these are the associated vals from the matrix corresponding to the first dAb
  for (val in 1:length(ordered_list_of_dabs)) {
    start_ind = 1 + length(f1_metrics) * (val - 1)
    temp_metrics = f1s[start_ind:(start_ind + 4)]
    temp_f1_df = rbind(temp_f1_df, temp_metrics)
  }
  temp_f1_df = na.omit(temp_f1_df)
  temp_perf_metrics = cbind(Accuracy, Kappa, temp_f1_df) %>% signif(digits =
                                                                      3)
  
  dabwise_performance = cbind(dabwise_performance, temp_perf_metrics)
  
  # Might add standard deviations later but for now Im leaving it out because it's a huge pain in the ass
  return(dabwise_performance)
}


# MULTI_PERFORMANCE_VIZ -> USED TO OBTAIN PERFORMANCE METRICS FROM A LIST OF RESULTS
# THE LIST OF RESULTS SHOULD FOLLOW THE SAME FORM AS THOSE PRODUCED IN GLM_CV FUNCTIONS
# HANDLES PLOTTING BUT CALLS THE EXTRACT_PERFORMANCE_METRICS FUNCTION TO ACTUALLY OBTAIN VALUES TO PLOT
#' @title multi_performance_viz
#' @usage Used to obtain performance metrics from a list of results. The list of
#' results should follow the same form as those produced in the glm_cv functions.
#' Also handles plotting but calls the extract_performance_metrics function to
#' actually obtain values to plot
#' @seealso extract_performance_metrics
#' 
#' @param list_of_results A list of model results as produced by any of the main
#' glm control functions
#' @param ordered_list_of_dabs A list of strings wherein each element is a dab
#' that I want I'm going to obtain performance metrics for
#' @param object_to_obtain (default = 'performance_metrics) used to index the 
#' list of results (in which each element is actually a list of results in itself)
#' and obtain specifically obtain the intended object
#'
#' @return A plot that contains a subplot of bar charts describing the performance
#' of each detection antibody
multi_performance_viz = function(list_of_results,
                                 ordered_list_of_dabs,
                                 object_to_obtain = 'performance_metrics') {
  # Obtaining performance metrics
  perf_metrics = extract_performance_metrics(list_of_results, ordered_list_of_dabs, object_to_obtain)
  long_performance_metrics = perf_metrics %>% pivot_longer(cols = -dAb, names_to = 'metric')
  
  # Creating a bar graphs with grouped cols for each dab and all dabs on a single 
  # graph organized by greatest kappa
  col_plot = ggplot(long_performance_metrics, aes(x = dAb, y = value, fill = metric)) +
    geom_bar(stat = "identity", position = "dodge") +
    scale_fill_brewer(palette = 'Dark2')
  print(col_plot)
  
  col_plot_2 = ggplot(long_performance_metrics,
                      aes(x = metric, y = value, fill = metric)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = value, y = value, vjust = 2)) +
    scale_fill_brewer(palette = 'Dark2', type = "qual") +
    facet_wrap( ~ dAb)
  
  return(col_plot_2)
  # I want the top N cAbs from each dab
}


# GLM PERFROMANCE FUNCTION
# Used in GLM_FS function to obtain performance metrics based on model predictions
#' @title get_glm_performance
#' @usage called in all main glm function to obtain desired performance metrics
#' and confusion matrices of a model based on the predicted and ground truth y 
#' values
#' @param predictions predicted class labels from a model
#' @param y_test ground truth class labels
#'
#' @return a data frame of performance metrics (currently accuracy, kappa, and F1)
get_glm_performance = function(predictions, y_test){
  
  if(length(levels(predictions)) == length(levels(y_test))){
    conf_matrix = confusionMatrix(y_test, predictions, mode = 'everything')
    
    Delta_F1 = conf_matrix$byClass[,'F1']['Class: Delta'][[1]]
    Omi_F1 = conf_matrix$byClass[,'F1']['Class: Omicron'][[1]]
    WT_F1 = conf_matrix$byClass[,'F1']['Class: WT'][[1]]
    Beta_F1 = conf_matrix$byClass[,'F1']['Class: Beta'][[1]]
    
    perf_metrics = data.frame("accuracy" = conf_matrix$overall[[1]],
                              "kappa" = conf_matrix$overall[[2]],
                              "Delta_F1" = Delta_F1,
                              "Omi_F1" = Omi_F1,
                              "WT_F1" = WT_F1,
                              "Beta_F1" = Beta_F1,
                              "Average_f1" = mean(c(Delta_F1, Omi_F1, WT_F1, Beta_F1), na.rm = T))
    return(list("perf_metrics" = perf_metrics,
                "conf_matrix" = conf_matrix))
  }
  else{
    reference = y_test
    conf_matrix = CrossTable(predictions,
                             y_test,
                             prop.chisq = F)
    table = conf_matrix$t
    
    # Convert table to df to make it easier to work with
    N = length(reference)
    tab_df = as.data.frame(table) %>% mutate_if(is.factor, as.character)
    
    # Accuracy is the sum of when the x and y are the same divided by the total
    accuracy = sum(tab_df$Freq[tab_df$x==tab_df$y]) / N
    
    # kappa is 
    # po is the accuracy of the model
    po = sum(tab_df$Freq[tab_df$x==tab_df$y])/sum(tab_df$Freq)
    
    # pe is the agreement between the probabilities of each class based on the reference and the predicted
    true_beta_percent = sum(tab_df$Freq[tab_df$y=='Beta']) / N
    pred_beta_percent = sum(tab_df$Freq[tab_df$x=='Beta']) / N 
    
    true_omicron_percent = sum(tab_df$Freq[tab_df$y=='Omicron']) / N 
    pred_omicron_percent = sum(tab_df$Freq[tab_df$x=='Omicron']) / N 
    
    true_delta_percent = sum(tab_df$Freq[tab_df$y=='Delta']) / N 
    pred_delta_percent = sum(tab_df$Freq[tab_df$x=='Delta']) / N 
    
    true_wt_percent = sum(tab_df$Freq[tab_df$y=='WT']) / N 
    pred_wt_percent = sum(tab_df$Freq[tab_df$x=='WT']) / N 
    
    pe = ((true_beta_percent * pred_beta_percent) / N^2) +
      ((true_omicron_percent * pred_omicron_percent) / N^2) + 
      ((true_delta_percent * pred_delta_percent) / N^2) +
      ((true_wt_percent * pred_wt_percent) / N^2)
    
    kap = (po - pe)/(1 - pe)
    
    perf_metrics = data.frame("accuracy" = accuracy,
                              "kappa" = kap,
                              "Delta_F1" = NA,
                              "Omi_F1" = NA,
                              "WT_F1" = NA,
                              "Beta_F1" = NA,
                              "Average_f1" = NA)
    
    return(list("perf_metrics" = perf_metrics,
                "conf_matrix" = table))
  }
}



multi_fs_parser = function(mixed_performance_df, dab){
  
  std.err = function(x, na.rm = T) sd(x, na.rm)/sqrt(length(x))
  
  fs_data = mixed_performance_df %>% filter(!global_iterations %in% c("Average","Stdev"))
  
  fs_data_avg = fs_data %>%
    select(accuracy, kappa, Average_f1, num_features) %>% 
    group_by(num_features) %>%
    summarise(across(everything(), .f = list(mean = mean, se = std.err), na.rm = T))
  
  colnames(fs_data_avg) = c(
    "num_features",
    "Mean Accuracy",
    "accuracy_se",
    "Mean Kappa",
    "kappa_se",
    "Mean Average F1",
    "Average_f1_se"
  )
  # Braindead way of doing this becuase gathers been deprecated for forever,
  # pivot_longer would be better but I can't be bothered right now
  fs_avg_long = fs_data_avg %>%
    gather(mean_type, mean_val, "Mean Accuracy", "Mean Kappa", "Mean Average F1") %>%
    gather(se_type, se_val, accuracy_se, kappa_se, Average_f1_se) %>% 
    distinct(mean_val, .keep_all = T)
  

  viz = ggplot(fs_avg_long, aes(num_features, y = mean_val, colour=mean_type)) + 
    geom_line() + 
    geom_point() + 
    geom_errorbar(aes(ymin=mean_val-se_val, ymax=mean_val+se_val), width = 0.2)+
    scale_x_reverse(breaks = seq(max(fs_data$num_features), min(fs_data$num_features),  by = -1))+
    ylim(0.2, 1.0) +
    labs(title = paste('Average Performance by number of pairs used in the model - ', dab),
         x = "Number of antibody pairs in model",
         y = "Mean performance metric value",
         color = "Performance Metric") +
    theme_classic()+
    scale_color_brewer(palette = 'Dark2')
  
  return(list("avg_df" = fs_avg_long,
              "plot" = viz))
}


model_call_loop = function(df_list,
                           model_type,
                           k1 = 5,
                           k2 = 5,
                           num_repeats = 3,
                           threshold = 0,
                           n_top_features = c()
                           ){
  result_list = list()
  if (model_type == "lasso_fs"){
    for (df in df_list){
      result_list[[length(result_list) + 1]] = glm_control_multi_fs(df = df, 
                                                                  target = "target_vector", 
                                                                  k1 = 5, 
                                                                  k2 = 5, 
                                                                  num_repeats = 3, 
                                                                  threshold = 0,
                                                                  n_top_features = c(10, 8, 7, 6, 5, 4, 3))
    }
  }
  return(result_list)
}

  















