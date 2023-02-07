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


external_resample = function(target,
                             global_train,
                             n_folds,
                             external_perf_df,
                             external_model_list,
                             external_conf_table_list) {
  # Create model matrix -- .-1 term removes the intercept
  model_matrix_formula = paste(target, "~", ".-1") %>% as.formula()
  
  external_folds = createFolds(global_train[[target]], k = n_folds)
  
  for (fold in 1:n_folds) {
    external_train = global_train[-external_folds[[fold]],]
    external_test = global_train[external_folds[[fold]],]
    
    x_external_train_matrix = model.matrix(model_matrix_formula, external_train)
    y_external_train = external_train[[target]]
    
    x_external_test_matrix = model.matrix(model_matrix_formula, external_test)
    y_external_test = external_test[[target]]
    
    # Create initial models This creates a set of models trained with
    # varying lambdas. From here, I can obtain the optimal lambda and its
    # associated coefficients. I train this internal loop without an
    # intercept because I don't actually care what it is, nor does it matter
    # if the model is made worse by it.
    
    # Hardcoding alpha = 1 since I always want lasso, would need to be changes
    # for elastic net
    external_model = cv.glmnet(
      x_external_train_matrix,
      y_external_train,
      family = 'multinomial',
      alpha = 1,
      type.multinomial = 'grouped',
      nfolds = 5,
      type.measure = 'class',
      intercept = F,
      parallel = T
    )
    print(external_model)
    opt_external_predict = predict(
      external_model,
      s = external_model$lambda.1se,
      newx = x_external_test_matrix,
      type = 'class'
    )
    
    opt_external_predict = as.factor(opt_external_predict)
    perf_call = get_glm_performance(opt_external_predict, y_external_test)
    
    external_conf_matrix = perf_call$conf_matrix
    external_conf_table_list[[length(external_conf_table_list) + 1]] = external_conf_matrix
    external_model_list[[length(external_model_list) + 1]] = external_model
    
    if (is.null(external_perf_df)) {
      external_perf_df = perf_call$perf_metrics
    }
    else{
      external_perf_df = rbind(external_perf_df, perf_call$perf_metrics)
    }
  }
  
  return(
    list(
      "external_conf_table_list" = external_conf_table_list,
      "external_opt_model_list" = external_model_list,
      "external_perf_df" = external_perf_df
    )
  )

}

resampling_function = function(target,
                               dataset,
                               n_folds,
                               perf_df,
                               model_list,
                               conf_table_list
                               ) {
  # Create model matrix -- .-1 term removes the intercept
  model_matrix_formula = paste(target, "~", ".-1") %>% as.formula()
  
  folds = createFolds(dataset[[target]], k = n_folds)
  
  for (fold in 1:n_folds) {
    train = dataset[-folds[[fold]],]
    test = dataset[folds[[fold]],]
    
    x_train_matrix = model.matrix(model_matrix_formula, train)
    y_train = train[[target]]
    
    x_test_matrix = model.matrix(model_matrix_formula, test)
    y_test = test[[target]]
    
    # Create initial models This creates a set of models trained with
    # varying lambdas. From here, I can obtain the optimal lambda and its
    # associated coefficients. I train this internal loop without an
    # intercept because I don't actually care what it is, nor does it matter
    # if the model is made worse by it.
    
    # Hardcoding alpha = 1 since I always want lasso, would need to be changes
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


global_resampling_function = function(target,
                                      train,
                                      test,
                                      best_features,
                                      perf_df,
                                      model_list,
                                      conf_table_list) {
  # The only difference between this and the normal resampling function is
  # that this function doesn't use k-fold CV, it just takes a train and test df
  # directly. It also downselects the df based on best features
  downselect_train = cbind(train[[target]], train[best_features])
  colnames(downselect_train)[1] = target
  
  downselect_test = cbind(test[[target]], test[best_features])
  colnames(downselect_test)[1] = target
  
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


