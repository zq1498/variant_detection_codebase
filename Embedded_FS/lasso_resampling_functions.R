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


#' resampling_function
#'
#' @param target the name of the target vector: str
#' @param dataset data that you want to split and use in the K-fold CV loop
#' @param n_folds how many folds you want to use
#' @param perf_df a current data frame of performance metrics, must be stored
#' after the function returns and passed in again on the next loop iteration. The
#' data added by this function is obtained from the get_glm_performance function
#' @param model_list a current list in which each element is a model object 
#' returned by glmnet
#' @param conf_table_list a current list of confusion matrices in which each
#' element is a confusion matrix returned by either the caret function or manually
#' constructed by the get_glm_performance function 
#' @param transformation a vector containing strings that denote the set of 
#' transformations you wish to apply to your dataset. Passed to the preProcess
#' Caret function
#'
#' @return three objects, (1) an updated model_list, (2) and updates confusion
#' matrix list, and (3) and updated perfromance metric list
#' c("scale", "center", "BoxCox")
#' @example external_resample = resampling_function(
#'  target = target,
#'  dataset = global_train,
#'  n_folds = k2,
#'  perf_df = external_perf_df,
#'  model_list = external_opt_model_list,
#'  conf_table_list = external_conf_table_list,
#'  transformation = c("scale", "center", "BoxCox")
#' )
resampling_function = function(target,
                               dataset,
                               n_folds,
                               perf_df,
                               model_list,
                               conf_table_list,
                               transformation
                               ) {
  # Create model matrix -- .-1 term removes the intercept
  model_matrix_formula = paste(target, "~", ".-1") %>% as.formula()
  
  folds = createFolds(dataset[[target]], k = n_folds)
  
  for (fold in 1:n_folds) {
    train = dataset[-folds[[fold]],]
    test = dataset[folds[[fold]],]
    
    preprocess = preProcess(train, method = transformation)
    train = predict(preprocess, train)
    test = predict(preprocess, test)
    
    x_train_matrix = model.matrix(model_matrix_formula, train)
    y_train = train[[target]]
    
    x_test_matrix = model.matrix(model_matrix_formula, test)
    y_test = test[[target]]
    
    # Create initial models This creates a set of models trained with
    # varying lambdas. From here, I can obtain the optimal lambda and its
    # associated coefficients. I train this internal loop without an
    # intercept because I don't actually care what it is, nor does it matter
    # if the model is made worse by it.
    
    # Hardcoding alpha = 1 since I always want lasso, would need to be changed
    # for elastic net
    model = cv.glmnet(
      x_train_matrix,
      y_train,
      family = 'multinomial',
      alpha = 1,
      type.multinomial = 'grouped',
      nfolds = 5,
      type.measure = 'class',
      intercept = F,
      parallel = T
    )
    opt_predict = predict(
      model,
      s = model$lambda.1se,
      newx = x_test_matrix,
      type = 'class'
    )
    
    opt_predict = as.factor(opt_predict)
    perf_call = get_glm_performance(opt_predict, y_test)
    
    conf_matrix = perf_call$conf_matrix
    conf_table_list[[length(conf_table_list) + 1]] = conf_matrix
    model_list[[length(model_list) + 1]] = model
    
    if (is.null(perf_df)) {
      perf_df = perf_call$perf_metrics
    }
    else{
      perf_df = rbind(perf_df, perf_call$perf_metrics)
    }
  }
  
  return(
    list(
      "conf_table_list" = conf_table_list,
      "opt_model_list" = model_list,
      "perf_df" = perf_df
    )
  )
}


#' global_resampling_function
#'
#' @param target the name of the target vector: str
#' @param train training data created by the global folds within the parent function
#' @param test testing data created by the global folds within the parent function
#' @param best_features a vector if pair names derived from the sorted_coefs
#' object created during the parent function run. The parent function controls
#' whether or not this item is a subset of the original dataset (if fs == T) or
#' if it is the entire dataset (if fs == F)
#' @param perf_df a current data frame of performance metrics, must be stored
#' after the function returns and passed in again on the next loop iteration. The
#' data added by this function is obtained from the get_glm_performance function
#' @param model_list a current list in which each element is a model object 
#' returned by glmnet
#' @param conf_table_list a current list of confusion matrices in which each
#' element is a confusion matrix returned by either the caret function or manually
#' constructed by the get_glm_performance function 
#' @param transformation a vector containing strings that denote the set of 
#' transformations you wish to apply to your dataset. Passed to the preProcess
#' Caret function
#'
#' @return
#' @note I think there is a way to restructure this file such that I use the same
#' call for both the internal and external resample, but the key problem is that
#' I'm not sure how to keep track of the the global and external folds simultaneously
#' without getting way overcomplicated. This is fine for now though.
#'
#' @examples global_resample = global_resampling_function(
#' target = "target_vector",
#' train = global_train,
#' test = global_test,
#' best_features = sorted_coefs$Pair[1:n_top_features],
#' perf_df = global_perf_df,
#' model_list = global_opt_model_list,
#' conf_table_list = global_conf_table_list,
#' transformation = c("scale", "center")
global_resampling_function = function(target,
                                      train,
                                      test,
                                      sorted_coefs,
                                      perf_df,
                                      model_list,
                                      conf_table_list,
                                      transformation) {
  # The only difference between this and the normal resampling function is
  # that this function doesn't use k-fold CV, it just takes a train and test df
  # directly. It also downselects the df based on best features
  
  downselect_train = cbind(train[[target]], train[best_features])
  colnames(downselect_train)[1] = target
  
  downselect_test = cbind(test[[target]], test[best_features])
  colnames(downselect_test)[1] = target
  
  preprocess = preProcess(downselect_train, method = transformation)
  downselect_train = predict(preprocess, downselect_train)
  downselect_test = predict(preprocess, downselect_test)
  
  model_matrix_formula = paste(target, "~", ".-1") %>% as.formula()
  
  x_train_matrix = model.matrix(model_matrix_formula, downselect_train)
  y_train = downselect_train[[target]]
  
  x_test_matrix = model.matrix(model_matrix_formula, downselect_test)
  y_test = downselect_test[[target]]
  
  model = cv.glmnet(
    x_train_matrix,
    y_train,
    family = 'multinomial',
    type.multinomial = 'grouped',
    type.measure = 'class',
    nfold = 5,
    parallel = T,
    alpha = 1,
    intercept = F
  )
  
  opt_predict = predict(model,
                        newx = x_test_matrix,
                        type = 'class',
                        s = model$lambda.1se)
  opt_predict = as.factor(opt_predict)
  
  perf_call = get_glm_performance(opt_predict, y_test)
  
  if (is.null(perf_df)) {
    perf_df = perf_call$perf_metrics
  }
  else{
    perf_df = rbind(perf_df, perf_call$perf_metrics)
  }
  
  conf_table_list[[length(conf_table_list) + 1]] = perf_call$conf_matrix
  model_list[[length(model_list) + 1]] = model
  
  return(
    list(
      "conf_table_list" = conf_table_list,
      "opt_model_list" = model_list,
      "perf_df" = perf_df
    )
  )
  
}


global_resampling_temp = function(target,
                                  train,
                                  test,
                                  sorted_coefs,
                                  num_feature_list,
                                  perf_df,
                                  model_list,
                                  conf_table_list,
                                  transformation) {
  
  # Prototype version of a resampling function that loops over a set of n_features
  # A new training set needs to be created each time I think?
  
  # Only need to preprocess the data once
  preprocess = preProcess(train, method = transformation)
  train = predict(preprocess, train)
  test = predict(preprocess, test)
  
  model_matrix_formula = paste(target, "~", ".-1") %>% as.formula()
  tester = 0
  # Starting for loop to iterate over the n_feature_list
  for (n_features in num_feature_list) {
    features = sorted_coefs$Pair[1:n_features]
    downselect_train = train %>% select(target, features)
    downselect_test = test %>% select(target, features)
    
    x_train_matrix = model.matrix(model_matrix_formula, downselect_train)
    y_train = downselect_train[[target]]
    
    x_test_matrix = model.matrix(model_matrix_formula, downselect_test)
    y_test = downselect_test[[target]]
    
    model = cv.glmnet(
      x_train_matrix,
      y_train,
      family = 'multinomial',
      type.multinomial = 'grouped',
      type.measure = 'class',
      nfold = 5,
      parallel = T,
      alpha = 1,
      intercept = F
    )
    
    opt_predict = predict(model,
                          newx = x_test_matrix,
                          type = 'class',
                          s = model$lambda.1se)
    
    opt_predict = as.factor(opt_predict)
    
    perf_call = get_glm_performance(opt_predict, y_test)
    
    dat = perf_call$perf_metrics %>% bind_cols(num_features = n_features)
    
    if (is.null(perf_df)) {
      perf_df = dat
    }
    else{
      perf_df = rbind(perf_df, dat)
    }
    conf_table_list[[length(conf_table_list) + 1]] = perf_call$conf_matrix
    model_list[[length(model_list) + 1]] = model
  }
  return(
    list(
      "conf_table_list" = conf_table_list,
      "opt_model_list" = model_list,
      "perf_df" = perf_df
    )
  )
  
}


