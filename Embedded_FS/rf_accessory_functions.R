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