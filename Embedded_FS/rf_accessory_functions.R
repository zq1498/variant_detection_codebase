library(tidyverse)
library(caret)
library(stats)
library(psych)
library(gmodels)
library(gridExtra)
library(doParallel)
library(roxygen2)
library(RColorBrewer)


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
get_rf_performance = function(predictions, y_test){
  
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


rf_coef_visualization = function(coef_df){
  # Calculate the stderrs and averages, bind them to the original data set, then select them
  # (I only bind to the data set so that it remains in the same data structure). Pivot them
  # longer, then wider again so that I have a pairs, average, and se.
  std.err = function(x, na.rm = T) sd(x, na.rm)/sqrt(length(x))
  summarized_coef_df = coef_df %>% 
    group_by(n_features) %>% 
    summarise(across(everything(), .f = list(mean = mean, se = std.err), na.rm = T))
  
  # From the summarized dataframe, select the means and add the number of features.
  # Pivot the data frame to have pairs as a column, add a measure value, then pivot back.
  # Remove the _mean from the end of each pair so that the pairs can be easily joined later.
  coef_df_means = summarized_coef_df %>% 
    select(ends_with("_mean")) %>% 
    bind_cols("n_features" = as.factor(summarized_coef_df$n_features)) %>% 
    pivot_longer(cols = -n_features,
                 names_to = "Pairs",
                 values_to = "Average") %>% 
    bind_cols("Measure" = 'Average') %>% 
    pivot_wider(names_from = Measure,
                values_from = "Average") %>% 
    mutate(across("Pairs", ~ gsub("_mean", "", .)))
  
  coef_df_se = summarized_coef_df %>% select(ends_with("_se")) %>% 
    bind_cols("n_features" = as.factor(summarized_coef_df$n_features)) %>% 
    pivot_longer(cols = -n_features,
                 names_to = "Pairs",
                 values_to = "Standard_errors") %>% 
    bind_cols("Measure" = 'se') %>% 
    pivot_wider(names_from = Measure,
                values_from = "Standard_errors") %>% 
    mutate(across("Pairs", ~ gsub("_se", "", .))) 
  
  modded_df = inner_join(coef_df_se, coef_df_means, by = c("Pairs", "n_features"))
  
  # Grabbing dAb names from modded_df and stripping them for nice plotting
  dab_name = regexpr(".+(?=_)",modded_df$Pairs, perl = T)
  dab_name = unique(regmatches(modded_df$Pairs, dab_name))
  
  strip_cAbs = gsub("^[^c]*","",modded_df$Pairs, perl = T) #!(c.*)
  modded_df$Pairs = strip_cAbs
  
  coef_plot = ggplot(modded_df, aes(x=reorder(Pairs, -Average), y=Average))+
    geom_bar(stat='identity', fill = "#0577B1")+
    geom_errorbar(aes(ymin=Average-se, ymax=Average+se), width=0.2)+
    labs(x="Pair", y="Average Importance", 
         title = paste("Average feature importance:", dab_name))+
    facet_wrap(~fct_rev(n_features), scales = "free_x")+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))

  
  return(coef_plot)
}



rf_performance_plotter = function(perf_df){
  
  blue_palette = c(
    `lighter_blue` = "#06a8f9",
      `ligh_blue` = "#0577B1",
      `blue` = "#005587",
      `dark_blue` = "#012169"
  ) 
  
  # dab_check = regexpr(".+(?=_)",modded_df$Pairs, perl = T)
  # dab_name = unique(regmatches(modded_df$Pairs, dab_check))

  std.err = function(x, na.rm = T) sd(x, na.rm)/sqrt(length(x))

  fs_data = perf_df %>% 
    select(accuracy, kappa, Average_f1, n_features) %>% 
    group_by(n_features) %>% 
    summarise(across(everything(), .f = list(mean = mean, se = std.err), na.rm = T))

  colnames(fs_data) = c(
    "num_features",
    "Mean Accuracy",
    "accuracy_se",
    "Mean Kappa",
    "kappa_se",
    "Mean Average F1",
    "Average_f1_se"
  )
  fs_avg_long = fs_data %>% 
    gather(mean_type, mean_val, "Mean Accuracy", "Mean Kappa", "Mean Average F1") %>% 
    gather(se_type, se_val, accuracy_se, kappa_se, Average_f1_se) %>% 
    distinct(mean_val, .keep_all = T)

  viz = ggplot(fs_avg_long, aes(num_features, y = mean_val, colour = mean_type)) +
    geom_line()+
    geom_point()+
    geom_errorbar(aes(ymin=mean_val-se_val, ymax=mean_val+se_val), width = 0.2)+
    scale_x_reverse(breaks = seq(max(fs_data$num_features), min(fs_data$num_features), b = -1))+
    ylim(0.2,1.0)+
    labs(tite = paste('Average Performace by number of features'),
         x = 'Number of antibody pairs in the model',
         y = "Mean Performance Metric Value",
         color = "Performance Metric")+
    theme_classic()+
    scale_color_manual(values = c(blue_palette[[1]], blue_palette[[3]], blue_palette[[4]]))
  
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
  if (model_type == "rf_fs"){
    for (df in df_list){
      result_list[[length(result_list) + 1]] = rf_control_function(df = df, 
                                                                    target = "target_vector", 
                                                                    k1 = 5, 
                                                                    k2 = 5, 
                                                                    num_repeats = 3, 
                                                                    threshold = 0,
                                                                    n_top_features = c(10, 8, 7, 6, 5, 4, 3)
                                                                   )
    }
  }
  else if (model_type == 'mars'){
    for (df in df_list){
      result_list[[length(result_list) + 1]] = mars_control_function(df = df, 
                                                                   target = "target_vector", 
                                                                   k1 = 5, 
                                                                   k2 = 5, 
                                                                   num_repeats = 3, 
                                                                   threshold = 0,
                                                                   n_top_features = c(10, 8, 7, 6, 5, 4, 3)
                                                                   )
    }
  }
  return(result_list)
}

