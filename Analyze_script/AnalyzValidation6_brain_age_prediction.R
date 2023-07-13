#####################################
# Brain age prediction analysis
#
# This script used to run the machine learning 
# pipeline in predicting the predicted brain age difference (PAD)
#
# Random forest is used as the classier and a feature selection 
# process will be conducted before model training.
#
# Liang Qunjun 2023/04/23

library(tidyverse)
library(rio)
library(ggstatsplot)
library(patchwork)
library(ggstatsplot)
library(ggeasy)
library(caret)
library(e1071)
library(performance)
library(glmnet)
library(purrr)

# load the data
source('Analyze_script/function_BrainAgeRegression.R')

## load demographic infomation
info_dat <- import('inputs/MDD_validation_dataset.xlsx')

## load the predicted brain ages
pad_file <- list.files(path = 'inputs/brainageR/', full.names = T)
pad_raw <- read_csv(pad_file)
pad_raw <- pad_raw %>% mutate(participant_id = str_extract(File, pattern = "sub-[0-9]{3}"))

## merge the demographic info to the predited brain age
pad_use <- pad_raw %>% filter(participant_id %in% info_dat$participant_id)
pad_use <- pad_use %>% left_join(info_dat)

# visualize the cor mat 

# plot use for publish -- scatter plot
ggplot(pad_use, aes(x = age, y = brain.predicted_age)) +
  geom_point(size = 3, alpha = 0.7) +
  geom_smooth(method = 'lm') +
  xlab('Chronological age') + ylab('Predicted age') +
  theme_classic() +
  easy_text_size(15)

ggsave(filename = 'output/Brainage_cor_point.png',
       height = 4, width = 5)

mean(pad_use$age)
mean(pad_use$brain.predicted_age)

# data preparetion before perform Machine learning  ----------------------------

dat_ml <- rio::import('output/validation_TDtopoDistance_machineLearning.csv')

## data preprocessing - scale the features
preProc <- dat_ml %>%
  preProcess(method = c('center', 'scale'))
dat_ml<- predict(preProc, dat_ml)

## age delta
pad_use <- pad_use %>% mutate(age_delta = age - brain.predicted_age)

## merge predicted brain age
dat_ml_use <- dat_ml %>% left_join(pad_use) %>% 
  dplyr::select(1:29, age, gender, QC_bold_fd_mean, age_delta) %>%
  mutate(gender_ind = factor(ifelse(gender=='male','0','1'))) %>%
  dplyr::select(-gender)

### scale the FD mean and age
dat_ml_use$QC_bold_fd_mean <- scale(dat_ml_use$QC_bold_fd_mean)[1:nrow(dat_ml_use)]
dat_ml_use$age <- scale(dat_ml_use$age)[1:nrow(dat_ml_use)]

# feature selection by correlation ---------------------------------------------
ggcorrmat(dat_ml_use[,c(2:29,32)],
          p.adjust.method = 'none')

feature_select_tmp <- ggcorrmat(dat_ml_use[,c(2:29,32)],
                                p.adjust.method = 'none', 
                                output = "dataframe") %>%
  filter(parameter2 == 'age_delta', (abs(estimate) >= 0.1)) %>% 
  .$parameter1

dat_ml_fs <- dat_ml_use %>% dplyr::select(feature_select_tmp, QC_bold_fd_mean, age_delta)


####################################
#
# Predict with random forest
#
####################################
source('Analyze_script/function_BrainAgeRegression.R')

# global variable
id_dependent <- 31 # column number for the dependent variable
tuneValue <- expand.grid(mtry = c(5:28)) # hyperparamter optimize
predict_model <- 'rf'
dat_ml_rf <- dat_ml_use[,-1]

# use LOOCV pipeline
result_loocv <- map_dfr(1:nrow(dat_ml_fs),
                        ~BrainAgeRegression_LOOCV(dat_ml_rf, n_dev = id_dependent,
                                                  out_sample = .x,
                                                  model = predict_model,
                                                  tuneTable = tuneValue),
                        .progress=T)
# export the LOOCV result
write_csv(result_loocv, file = bruceR::Glue('output/Brainage_prediction_LOOCV_{predict_model}.csv'))

# use random split pipeline
seed_select <- sample(1:1000, size = 500)
result_ranSplit <- map_dfr(seed_select,
                           ~BrainAgeRegression_randomSplit(dat_ml_rf, n_dev = id_dependent,
                                                           randomSeed = .x,
                                                           model = predict_model,
                                                           tuneTable = tuneValue),
                           .progress=T)
write_csv(result_ranSplit, file = bruceR::Glue('output/Brainage_prediction_randomSplit_{predict_model}.csv'))

####################################
#
# Result visualization
#
####################################

# load the brain age prediction
pred_age_loocv <- rio::import('output/Brainage_prediction_LOOCV_rf.csv')
pred_age_random <- rio::import('output/Brainage_prediction_randomSplit_rf.csv')

cor.test(pred_age_loocv$pred_value, pred_age_loocv$true_value)
ggplot(data = pred_age_loocv, aes(x = pred_value, y = true_value)) +
  geom_point(size =2.5, alpha = .7) +
  geom_abline(slope = 1, linetype = 'dashed', color = 'steelblue') +
  xlab('Predicted brain-PAD') + ylab('Actual brain-PAD') +
  theme_classic() + easy_text_size(15)
ggsave('output/Brainage_prediction_LOOCV_randomForest.png', width = 6, height = 5)

ggplot(data = pred_age_random) +
  geom_histogram(aes(x = cor), bins = 20, 
                 fill = 'steelblue', color = 'white') +
  ylab('Frequency') + xlab('Correlation') + ggtitle('Correlation distribution') +
  theme_classic() + easy_text_size(15) + easy_center_title()
ggsave('output/Brainage_prediction_randomSplit_correlation.png', width = 6, height = 5)

ggplot(data = pred_age_random) +
  geom_histogram(aes(x = rmse), bins = 20, 
                 fill = 'steelblue', color = 'white') +
  ylab('Frequency') + xlab('Root Mean Squared Error') + ggtitle('RMSE distribution') +
  theme_classic() + easy_text_size(15) + easy_center_title()
ggsave('output/Brainage_prediction_randomSplit_RMSE.png', width = 6, height = 5)
