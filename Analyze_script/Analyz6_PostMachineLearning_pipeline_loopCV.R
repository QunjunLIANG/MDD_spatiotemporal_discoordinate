###################
#
# This script is used to 
# collect the results from prvious 
# machine learning pipeline and 
# visualize for publish
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
library(purrr)
library(patchwork)
library(wordcloud)

############################################################
#
# Model performance - 2-fold CV
#
#############################################################

dat_result <- rio::import('output/MachineLearning_2FoldCV_result.csv')

psych::describe(dat_result)

# plot the ACC, sensitivity and specificity
p_acc <- ggplot(dat_result) +
  geom_histogram(bins = 15, 
                 aes(x = Accuracy), 
                 fill = 'steelblue',
                 color="white", alpha = 0.8) +
  ylab('Frequency') + xlab('Value') + ggtitle('Accuracy distribution') +
  theme_classic(base_size = 15) +
  easy_center_title()
p_acc

p_sens <- ggplot(dat_result) +
  geom_histogram(bins = 10, 
                 aes(x = Sensitivity), 
                 fill = 'steelblue',
                 color="white", alpha = 0.8) +
  ylab('Frequency') + xlab('Value') + ggtitle('Sensitivity distribution') +
  theme_classic(base_size = 15) +
  easy_center_title()
p_sens  

p_spc <- ggplot(dat_result) +
  geom_histogram(bins = 10, 
                 aes(x = Specificity), 
                 fill = 'steelblue',
                 color="white", alpha = 0.8) +
  ylab('Frequency') + xlab('Value') + ggtitle('Specificity distribution') +
  theme_classic(base_size = 15) +
  easy_center_title()
p_spc

p_performance <- p_acc + p_sens + p_spc
ggsave(filename = 'output/MachineLearning_performance.png', width = 15, height = 5)

# feature importance 

feaImp <- dat_result[,c(4:31)] %>%
  summarise_each(funs(mean)) %>%
  melt()
wordcloud(feaImp$variable, feaImp$value, 
          colors = hcl.colors(30, palette = "Set2"),
          scale = c(3,.2), rot.per = .4)

############################################################
#
# Model performance - LOOVC
#
#############################################################

dat_result <- read_csv('output/MachineLearning_LOOCV_result.csv')
model_perf <- confusionMatrix(factor(dat_result$pred_group),
                              factor(dat_result$actual))
## reshape to visual
v_perf <- model_perf$overall
v_perf <- data.frame(variable = names(v_perf),
                     value = v_perf)
c_perf <- model_perf$byClass
c_perf <- data.frame(variable = names(c_perf),
                     value = c_perf)

m_perf <- rbind(v_perf, c_perf)[c(1,2,8,9),]
m_perf$value <- round(m_perf$value, digits = 2)
p_loocv <- ggplot(m_perf, aes(x = variable, y = value)) +
  geom_bar(stat = 'identity', width = .5,
           fill = '#e27283', alpha =.8) +
  geom_text(aes(label = value),vjust = -0.5) +
  ggtitle('Model performance in LOOCV') +
  ylim(0, .8) +
  theme_classic(base_size = 15) + 
  easy_remove_x_axis(what = 'title')+
  easy_center_title() +
  easy_rotate_x_labels(angle = -30)

p_loocv + p_acc + p_sens + p_spc
ggsave(filename = 'output/MachineLearning_performance_LOOCV.png',
       height = 6, width = 9)

# Feature importance

feaImp <- dat_result[,c(6:33)] %>%
  summarise_each(funs(mean)) %>%
  melt()

wordcloud(feaImp$variable, feaImp$value, 
          colors = hcl.colors(30, palette = "Set2"),
          scale = c(3,0.5), rot.per = .4)
