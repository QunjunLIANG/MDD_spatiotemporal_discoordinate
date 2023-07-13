BrainAgeRegression_LOOCV <- function(data, out_sample, tuneTable, 
                                     n_dev, model){
  
  library(tidyverse)
  library(caret)
  library(e1071)
  
  svm_model <- train(
    age_delta ~ ., data = data[-out_sample,],
    method = model,
    tuneGrid = tuneTable,
    trControl = trainControl(method = "LOOCV")
  )
  ## predictions
  predicted_data <- predict(svm_model, data[out_sample,-n_dev])
  
  outTable <- data.frame(
    out_sample = out_sample,
    pred_value = predicted_data,
    true_value = data[out_sample, n_dev]
  )
  
  return(outTable)
}

BrainAgeRegression_randomSplit <- function(data, randomSeed, model,
                                           tuneTable, n_dev){
  library(tidyverse)
  library(caret)
  library(e1071)
  
  set.seed(randomSeed)
  
  # split into train & test dataset
  # 2 fold split
  split_ind <- createFolds(data[,n_dev], k=2)
  data_train <- data[split_ind$Fold1,]
  data_test <- data[split_ind$Fold2, ]
  
  svm_model <- train(
    age_delta ~ ., data = data_train,
    method = model,
    tuneGrid = tuneTable,
    trControl = trainControl(method = "cv", p = 0.7)
  )
  
  ## predictions
  predicted_data <- predict(svm_model, data_test[,-n_dev])
  
  ## obtain model performance
  cor_test <- cor.test(predicted_data, data_test[,n_dev])
  cor_value <- as.numeric(cor_test$estimate)
  cor_p <- cor_test$p.value
  
  ## calculate prediction error and RMSE
  error <- predicted_data - data_test[,n_dev]
  mse <- mean(error^2)
  rmse <- sqrt(mse)
  
  ## calculate efficiency
  mean_value <- mean(data_test[,n_dev])
  baseline_predicted <- rep(mean_value, length(data_test[,n_dev]))
  baseline_error <- baseline_predicted - data_test[,n_dev]
  
  outTable <- data.frame(
    seed = randomSeed,
    cor = cor_value,
    cor_p = cor_p,
    mse = mse,
    rmse = rmse,
    baseline_error = sqrt(mean(baseline_error^2))
  )
  
  return(outTable)
}
