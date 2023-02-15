library(tidyverse)
library(caret)
library(ggplot2)
library(stats)
library(psych)
library(gmodels)
library(gridExtra)
library(doParallel)
library(roxygen2)

rf_resampling_function = function(target,
                                  dataset,
                                  n_folds,
                                  perf_df,
                                  coef_df,
                                  conf_table_list,
                                  transformation){
  # Function is almost identical to the lasso resample function, with the 
  # exception that instead of a model list we pass a data frame of coefficients.
  # The model list in the lasso resampling function is essentially only used to
  # generate coefficients, so this cuts out the middle man
  model_formula = paste(target, "~", ".") %>% as.formula()
  
  folds = createFolds(dataset[[target]], k = n_folds)
  
  train_contrl = trainControl(method = "repeatedcv",
                          number = 5,
                          repeats = 3)
  
  for (fold in 1:n_folds){
    train = dataset[-folds[[fold]], ]
    test = dataset[folds[[fold]], ]
    
    preprocess = preProcess(train, method = transformation)
    train = predict(preprocess, train)
    test = predict(preprocess, test)
    
    x_train = train %>% select(-target)
    y_train = train[[target]]
    
    x_test = test %>% select(-target)
    y_test = test[[target]]
    
    rf_train = train(model_formula,
                     data = train,
                     method = "ranger",
                     importance = "impurity",
                     trControl = train_contrl,
                     verbose = T)
    
    opt_predict = predict(rf_train$finalModel, x_test,
                          type = 'response')
    opt_predict = opt_predict$predictions %>% as.factor()
    
    perf_call = get_rf_performance(opt_predict, y_test)
    
    coef_imp = varImp(rf_train)
    # Transposing the matrix is perfect, since it makes the row names the colnames
    # which can then be passed to bind_rows. Bind_rows is ideal because it matches
    # by column name.
    coef_imp = as.data.frame(t(coef_imp$importance))
    
    if (is.null(perf_df)){
      perf_df = perf_call$perf_metrics
    }
    else{
      perf_df = bind_rows(perf_df, perf_call$perf_metrics)
    }
    
    if (is.null(coef_df)){
      coef_df = coef_imp
    }
    else{
      coef_df = coef_df %>% bind_rows(coef_imp)
    }
    conf_table_list[[length(conf_table_list)+1]] = perf_call$conf_matrix
    
  }
  return(list("perf_df" = perf_df,
              "coef_df" = coef_df,
              "conf_table" = conf_table_list))
}


multi_feature_global_rf_resampling = function(target,
                                              train,
                                              test,
                                              sorted_coefs,
                                              num_feature_list,
                                              perf_df,
                                              coef_df,
                                              conf_table_list,
                                              transformation) {

  
  train_contrl = trainControl(method = "repeatedcv",
                              number = 5,
                              repeats = 3)
  
  model_formula = paste(target, "~", ".") %>% as.formula()
  
  for (n_features in num_feature_list){
    preprocess = preProcess(train, method = transformation)
    train = predict(preprocess, train)
    test = predict(preprocess, test)
    
    features = sorted_coefs$Pairs[1:n_features]
    downselect_train = train %>% select(target, features)
    downselect_test = test %>% select(target, features)
    
    x_test = downselect_test %>% select(-target)
    y_test = downselect_test[[target]]
    
    rf_train = train(model_formula,
                     data = downselect_train,
                     method = "ranger",
                     importance = "impurity",
                     trControl = train_contrl,
                     verbose = T)

    opt_predict = predict(rf_train$finalModel, x_test,
                          type = 'response')
    opt_predict = opt_predict$predictions %>% as.factor()

    perf_call = get_rf_performance(opt_predict, y_test)
    
    coef_imp = varImp(rf_train)
    # Transposing the matrix is perfect, since it makes the row names the colnames
    # which can then be passed to bind_rows. Bind_rows is ideal because it matches
    # by column name.
    coef_imp = as.data.frame(t(coef_imp$importance))
    
    if (is.null(perf_df)){
      perf_df = perf_call$perf_metrics
    }
    else{
      perf_df = bind_rows(perf_df, perf_call$perf_metrics)
    }
    
    if (is.null(coef_df)){
      coef_df = coef_imp
    }
    else{
      coef_df = coef_df %>% bind_rows(coef_imp)
    }
    conf_table_list[[length(conf_table_list)+1]] = perf_call$conf_matrix
  }
  return(list("perf_df" = perf_df,
              "coef_df" = coef_df,
              "conf_table" = conf_table_list))

}


