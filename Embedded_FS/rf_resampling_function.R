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
                                  transformation){
  model_formula = paste(target, "~", ".") %>% as.formula()
  
  folds = createFolds(dataset[[target]], k = n_folds)
  
  train_contrl = trainControl(method = "repeatedcv",
                          number = 5,
                          repeats = 3)
  
  for (fold in 1:n_folds){
    train = dataset[-folds[[fold]], ]
    test = dataset[folds[[fold]], ]
    
    preprocess = preProcess(train, method = tramsformation)
    train = predict(preprocess, train)
    test = predict(preprocess, test)
    
    rf_train = train(model_formula,
                     data = train,
                     method = "ranger",
                     importance = "impurity",
                     trControl = train_contrl,
                     verbose = T)
    
    opt_predict = predict()
    
    
  }
}


#shit ass rf_train function of yore for reference
```{r, ranger_cv_function}
# Outer cross-validation loop is the training data (recombinant) and testing data is 
# the clinical data

# The secondary cross-validation loop k-fold cross-validation of the recombinant data. This 
# includes validation of the model trained in the tertiary loop, and subsequent feature
# selection

# The tertiary cross-validation loop is k-fold CV via caret or glmnet

recombinant_cv = function(clinical_df, recombinant_df, target, class_dist = F,
                          model_type, tune_grid, train_control){
  # I can pretty easily modify this function to be a generalized method by adding
  # a tune grid and a model method argument.
  # Class distance indicates whether the training and testing sets should be 
  # modified to include the Mahalanobis class distances as additional predictors
  confusion_list = list()
  model_list = list()
  mod_index = 1
  
  best_tunes = NULL
  internal_pred_df = NULL
  
  if (class_dist == T){
    class_distances = getCentroids(recombinant_df, clinical_df, target)
    clinical_df = cbind(clinical_df, class_distances$test_distances)
    recombinant_df = cbind(recombinant_df, class_distances$train_distances)
  }
  
  x_clinical_test = clinical_df %>% dplyr::select(-target)
  y_clinical_test = clinical_df$Variant
  
  global_train_folds = createFolds(recombinant_df[[target]], k=5)
  
  for (fold in 1:length(global_train_folds)){
    global_train = recombinant_df[-global_train_folds[[fold]],]
    global_test = recombinant_df[global_train_folds[[fold]],]
    
    x_global_train = global_train %>% dplyr::select(-target)
    y_global_train = global_train[[target]]
    
    x_global_test = global_test %>% dplyr::select(-target)
    y_global_test = global_test[[target]]
    
    internal_mod = train(x = x_global_train,
                         y = y_global_train,
                         # preProc = c("center", "scale"),
                         method = model_type,
                         trControl = train_control,
                         metric = 'Kappa',
                         tuneGrid = tune_grid)
    
    # Using the final model obtained in the cross validation, predict on the left out data
    internal_preds = predict(internal_mod$finalModel, x_global_test, 
                             type = 'response')
    
    accuracy = sum(ifelse(internal_preds$predictions == y_global_test,
                          1, 0))/
      length(y_global_test)
    
    tune_acc = cbind(internal_mod$bestTune, accuracy)
    
    # Adds best tuning parameters onto the tuning parameters dataframe
    if (is.null(best_tunes)){
      best_tunes = tune_acc
    }
    else {
      best_tunes = rbind(best_tunes, tune_acc)
    }
    
    confusionTable = confusionMatrix(data = internal_preds$prediction, 
                                     reference = y_global_test)
    model_list[[mod_index]] = c(model_list, internal_mod)
    confusion_list[[mod_index]] = c(confusion_list, confusionTable)
    mod_index = mod_index + 1
  }
  
  x_recombinant_full = recombinant_df %>% dplyr::select(-target)
  y_recombinant_full = recombinant_df[[target]]
  
  clinical_model = train(x = x_recombinant_full,
                         y = y_recombinant_full,
                         trControl = train_control,
                         method = model_type,
                         metric = 'Kappa',
                         tuneGrid = tune_grid)
  
  clinical_preds = predict(clinical_model$finalModel, x_clinical_test,
                           type = 'response')
  
  clinical_confusion = confusionMatrix(data = clinical_preds$predictions,
                                       reference = y_clinical_test)
  
  return(list("model_list" = model_list,
              "tuner" = best_tunes,
              "confusionMats" = confusion_list,
              "final_model" = c(clinical_model),
              "testConfusionMat" = clinical_confusion))
}
subsetted_clinical_df = known_var_clinical_data %>% filter(`1B2` > 3.5)

clinical_model_data = subsetted_clinical_df %>% dplyr::select(-`1B2`)

tuneGrid = expand.grid(mtry = c(1:3),
                       min.node.size = c(1:5),
                       splitrule = 'gini')

trainControl = trainControl(method = 'repeatedcv',
                            number = 5,
                            repeats = 3)
# adaptive = list(min = 5, alpha = 0.05,
#                 method = 'BT', complete = T),
# search = 'random')

rf_train = recombinant_cv(clinical_model_data, full_recombinant_data, "Variant",
                          model_type = 'ranger', tune_grid = tuneGrid, 
                          train_control = trainControl)


whole_set_performance = rf_train$testConfusionMat
whole_set_model = rf_train$final_model
whole_set_performance
whole_set_model
```