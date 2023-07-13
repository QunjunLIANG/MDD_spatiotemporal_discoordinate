# leave one site off
RandomLOSOCVpipeline <- function(data, site){
  library(caret)
  library(randomForest)
  library(tidyverse)
  
  # split the data
  data_train <- data %>% filter(center != site)
  data_test <- data %>% filter(center == site)

  # clean the site variable
  data_train <- data_train %>% select(-center)
  data_test <- data_test %>% select(-center)

  # train random forest
  fit_control <- trainControl(method = 'cv', p = .7)
  model_train <- train(group_ind ~.,
                       data = data_train,
                       method = "rf",
                       tuneGrid = expand.grid(mtry = c(5:28)),
                       trControl = fit_control)
  # record the number of select features
  n_feature <- model_train[["finalModel"]][["mtry"]]

  # the ACC of training model
  acc_all <- model_train[["results"]][["Accuracy"]]
  acc_opt <- max(acc_all) %>% round(digits = 3)
  kapa_all <- model_train[["results"]][["Kappa"]]
  kapa_opt <- max(kapa_all) %>% round(digits = 3)

  # record variable importantce
  var_imp <- varImp(model_train)
  var_imp <- var_imp$importance %>%
    mutate(Overall = rank(Overall)) %>%
    t(.) %>% as.data.frame()

  # model efficiency
  predict_result <- predict(model_train, data_test[,-which(names(data_test)=='group_ind')])
  model_perf <- confusionMatrix(predict_result, data_test$group_ind)
  model_perf_overall <- model_perf$overall %>% as.data.frame() %>%
    t()%>% as.data.frame()
  model_perf_byClass <- model_perf$byClass %>% as.data.frame()%>%
    t()%>% as.data.frame()

  outData <- data.frame(n_feature = n_feature,
                        train_acc = acc_opt,
                        train_kapa = kapa_opt) %>%
    cbind(var_imp,model_perf_overall,model_perf_byClass)

  return(outData)
}

# random 2 fold split
Random2FCVpipeline <- function(data, randomSeed, n_dev){
  library(caret)
  library(randomForest)
  library(tidyverse)
  
  set.seed(randomSeed)
  # 2 fold split
  split_ind <- createFolds(data[,n_dev], k=2)
  data_train <- data[split_ind$Fold1,]
  data_test <- data[split_ind$Fold2, ]
  # train random forest
  fit_control <- trainControl(method = 'cv', p = .7)
  model_train <- train(group_ind ~.,
                       data = data_train,
                       method = "rf",
                       tuneGrid = expand.grid(mtry = c(5:28)),
                       trControl = fit_control)
  # record the number of select features
  n_feature <- model_train[["finalModel"]][["mtry"]]

  # the ACC of training model
  acc_all <- model_train[["results"]][["Accuracy"]]
  acc_opt <- max(acc_all) %>% round(digits = 3)
  kapa_all <- model_train[["results"]][["Kappa"]]
  kapa_opt <- max(kapa_all) %>% round(digits = 3)

  # record variable importantce
  var_imp <- varImp(model_train)
  var_imp <- var_imp$importance %>%
    mutate(Overall = rank(Overall)) %>%
    t(.) %>% as.data.frame()

  # model efficiency
  predict_result <- predict(model_train, data_test[,-n_dev])
  model_perf <- confusionMatrix(predict_result, data_test$group_ind)
  model_perf_overall <- model_perf$overall %>% as.data.frame() %>%
    t()%>% as.data.frame()
  model_perf_byClass <- model_perf$byClass %>% as.data.frame()%>%
    t()%>% as.data.frame()

  outData <- data.frame(n_feature = n_feature,
                        train_acc = acc_opt,
                        train_kapa = kapa_opt) %>%
    cbind(var_imp,model_perf_overall,model_perf_byClass)

  return(outData)
}

# site-based 2 fold split
RandomSite2FCVpipeline <- function(data, randomSeed){
  library(caret)
  library(randomForest)
  library(tidyverse)
  
  set.seed(randomSeed)
  # 2 fold split
  split_ind <- groupKFold(data$center, k=2)
  data_train <- data[split_ind$Fold1,]
  data_test <- data[split_ind$Fold2, ] %>% select(-center)
  # train random forest - site based 10 fold
  train_split <- groupKFold(data_train$center)
  data_train <- data_train %>% select(-center)
  fit_control <- trainControl(method = 'cv', 
                              index = train_split, number = 2)
  model_train <- train(group_ind ~.,
                       data = data_train,
                       method = "rf",
                       tuneGrid = expand.grid(mtry = c(5:28)),
                       trControl = fit_control)
  # record the number of select features
  n_feature <- model_train[["finalModel"]][["mtry"]]
  
  # the ACC of training model
  acc_all <- model_train[["results"]][["Accuracy"]]
  acc_opt <- max(acc_all) %>% round(digits = 3)
  kapa_all <- model_train[["results"]][["Kappa"]]
  kapa_opt <- max(kapa_all) %>% round(digits = 3)
  
  # record variable importantce
  var_imp <- varImp(model_train)
  var_imp <- var_imp$importance %>%
    mutate(Overall = rank(Overall)) %>%
    t(.) %>% as.data.frame()
  
  # model efficiency
  predict_result <- predict(model_train, data_test[,-which(names(data_test)=='group_ind')])
  model_perf <- confusionMatrix(predict_result, data_test$group_ind)
  model_perf_overall <- model_perf$overall %>% as.data.frame() %>%
    t()%>% as.data.frame()
  model_perf_byClass <- model_perf$byClass %>% as.data.frame()%>%
    t()%>% as.data.frame()
  
  outData <- data.frame(n_feature = n_feature,
                        train_acc = acc_opt,
                        train_kapa = kapa_opt) %>%
    cbind(var_imp, model_perf_overall, model_perf_byClass)
  
  return(outData)
}