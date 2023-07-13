###################
#
# This script contains the pipeline for machine learning
# approaches trained by the
# nested LOOCV and nested two-fold cross-validation
#
# Liang Qunjun 

library(tidyverse)
library(bruceR)
library(caret)
library(ggstatsplot)
library(viridis)
library(pROC)
library(RColorBrewer)
library(ggeasy)
library(ROCR)
library(purrr)
library(furrr)
library(randomForest)

############################################################
#
# Data preprocessing
#
#############################################################

dat_info <- rio::import('inputs/MDD_HC_match_data.xlsx') %>%
  select(participant_id, gender, age, QC_bold_fd_mean) %>%
  mutate(gender_ind = factor(ifelse(gender == 'male',0,1))) %>%
  select(-gender)
 
dat_ML <- rio::import('output/Topology_Distance_machineLearning.csv')%>%
  left_join(dat_info) %>% 
  mutate(group_ind = factor(ifelse(group == 'MDD', 1, 0))) %>% 
  select(-participant_id, -group)
 
## data preprocessing - scale the features
preProc <- dat_ML %>%
  preProcess(method = c('center', 'scale'))
dat_predict <- predict(preProc, dat_ML)

############################################################
#
# Randomforest - LOOCV
#
#############################################################

## define function
MachineLeanringPip_LOOCV <- function(data, n_sbj){
  dat_train <- data[-n_sbj,]
  dat_test <- data[n_sbj,]
  # model training also with LOOCV
  fit_control <- trainControl(method = 'LOOCV')
  model_train <- train(group_ind ~.,
                       data = dat_train,
                       method = "rf",
                       tuneGrid = expand.grid(mtry = 5:26),
                       trControl = fit_control)
  # record the number of select features
  n_feature <- model_train[["finalModel"]][["mtry"]]

  # the ACC of training model
  acc_all <- model_train[["results"]][["Accuracy"]]
  acc_opt <- acc_all[n_feature-4]
  kapa_all <- model_train[["results"]][["Kappa"]]
  kapa_opt <- kapa_all[n_feature-4]

  # record variable importantce
  var_imp <- varImp(model_train)
  var_imp <- var_imp$importance %>%
    mutate(Overall = rank(Overall)) %>%
    t(.) %>% as.data.frame()

  # colllect the results
  predict_result <- predict(model_train, dat_test)
  predict_result_prob <- predict(model_train, dat_test, type = 'prob') %>%
    as.data.frame()
  colnames(predict_result_prob) <- c('pred_zero','pred_one')
  outData <- data.frame(sbj = n_sbj,
                        pred_group = predict_result,
                        n_feature = n_feature,
                        train_acc = acc_opt,
                        train_kapa = kapa_opt) %>%
    cbind(var_imp,predict_result_prob)

  return(outData)
}

# calling furrr for parallel calculaton
print("Start LOOCV pipeline")
date()
plan(multisession)
ML_result <- future_map_dfr(c(1:length(dat_predict$group_ind)),
                           ~ MachineLeanringPip_LOOCV(data = dat_predict, .x),
                           .progres = T)
ML_result['actual'] <- dat_predict$group_ind
write_csv(ML_result, file = 'output/MachineLearning_LOOCV_result.csv')
print("Finished LOOCV pipeline")
date()

############################################################
#
# Randomforest - Random 2-fold CV
#
#############################################################
Random2FCVpipeline <- function(data, randomSeed, n_dev = 31){
  set.seed(randomSeed)
  # 2 fold split
  split_ind <- createFolds(data[,n_dev], k=2)
  data_train <- data[split_ind$Fold1,]
  data_test <- data[split_ind$Fold2, ]
  # train random forest
  fit_control <- trainControl(method = 'cv', 
                              p = .7)
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

# choose 500 random seeds 
seeds <- sample(1:1000, 500)
print("Start Randomly assign machine learning pipeline")
date()
#plan(multisession)
result_random <- map_dfr(seeds, 
                         ~Random2FCVpipeline(data = dat_predict, .x,
                                             n_dev = ncol(dat_predict)),
                                            .progress = T)
write_csv(result_random, 'output/MachineLearning_2FoldCV_result.csv')
print("All finished!")
date()  


