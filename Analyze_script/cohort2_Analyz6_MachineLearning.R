############################################
#
# Machine learning pipeline for
# REST.
#
# The pipeline is the identified pipeline used
# in the main dataset.
#
# Liang Qunjun 2023/04/06

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
library(patchwork)
library(wordcloud)

# load the data for machine learning

dat_ml_raw <- rio::import('output/REST_spttopo_machineLearning.csv')

# Testing the difference acorss the SPT features

dat_ml_long <- dat_ml_raw %>% reshape2::melt(id.vars = 1:7)

MANOVA(data = dat_ml_long,
       subID = 'ID', dv = 'value', between = 'group',
       within = 'variable', covariate = c('Age','gender','fd_mean'),
       file = 'output/REST_repeatedMANOVA.doc') %>%
  EMMEANS(effect = 'group', by = 'variable')


## visualization with matrix
dat_vis <- dat_ml_raw %>% dplyr::select(group, ID, 8:ncol(dat_ml_raw)) %>% data.table::melt(id.vars = 1:2) %>% group_by(group, variable) %>%
  summarise(mean_dist = mean(value)) %>% mutate(scale_dist = (mean_dist-mean(mean_dist))/sd(mean_dist))

Net_id <- dat_vis$variable %>% as.character()
Net_id <- strsplit(Net_id, split = '_') %>% unlist()
fromNet <- Net_id[seq(1,length(Net_id)-1,2)]
toNet <-  Net_id[seq(2,length(Net_id),2)]

dat_vis['fromNet'] <- fromNet
dat_vis['toNet'] <- toNet

dat_vis$fromNet <- factor(dat_vis$fromNet, levels = c("somMot","VIS",'DAN',"VAN","FPC","Limbic","DMN"))
dat_vis$toNet <- factor(dat_vis$toNet, levels = c("somMot","VIS",'DAN',"VAN","FPC","Limbic","DMN"))

# building HC distance matrix

NetDistancePlot <- function(network_ind){
  library(ggcorrplot)
  # make an empty matrix
  mat <- matrix(nrow = length(unique(network_ind$fromNet)), 
                ncol = length(unique(network_ind$fromNet)))
  colnames(mat) <- unique(network_ind$fromNet)
  rownames(mat) <- unique(network_ind$fromNet)
  
  for (i in 1:nrow(network_ind)) {
    net_tmp <- network_ind[i,]
    from_tmp <- net_tmp$fromNet
    to_tmp <- net_tmp$toNet
    mat[from_tmp, to_tmp] <- net_tmp$scale_dist
    mat[ to_tmp,from_tmp] <- net_tmp$scale_dist
  }
  
  p <- ggcorrplot(mat, type = 'lower', show.diag = T, legend.title = "Distance",
                  ggtheme = theme_minimal(base_size = 20)) +
    scale_fill_gradient2(
      low = 'steelblue',
      mid = 'white',
      high = "darkred") +
    easy_text_size(20) + easy_add_legend_title('Distance')
  
  return(list(distMat=mat, distPlot=p))
}

## HC matrix
network_hc <- dat_vis %>% filter(group == 'HC')
distMat_list <- NetDistancePlot(network_ind = network_hc)
disMat_hc <- distMat_list$distMat
disMat_hc_plot <- distMat_list$distPlot
## MDD matrix
network_mdd <- dat_vis %>% filter(group == 'MDD')
distMat_list <- NetDistancePlot(network_ind = network_mdd)
disMat_mdd <- distMat_list$distMat
disMat_mdd_plot <- distMat_list$distPlot

## difference in distance matrix 
disMat_diff <- disMat_mdd - disMat_hc
p_groupDiff <- ggcorrplot(disMat_diff, type = 'lower', show.diag = T, legend.title = "Distance",
           ggtheme = theme_minimal(base_size = 20)) +
  scale_fill_gradient2(
    low = 'steelblue',
    mid = 'white',
    high = "darkred") +
  easy_text_size(20) + easy_add_legend_title('Distance')

## export the PNGs
Cairo::Cairo(file = 'output/REST_DistanceMat_spttopo_HC.png')
print(disMat_hc_plot)
dev.off()

Cairo::Cairo(file = 'output/REST_DistanceMat_spttopo_MDD.png')
print(disMat_mdd_plot)
dev.off()

Cairo::Cairo(file = 'output/REST_DistanceMat_spttopo_groupDiff.png')
print(p_groupDiff)
dev.off()

## data preprocessing - scale the features -------------------------------------
preProc <- dat_ml_raw %>%
  preProcess(method = c('center', 'scale'))
dat_predict <- predict(preProc, dat_ml_raw) %>%
  mutate(group_ind = factor(ifelse(group=='MDD',1,0)),
         gender_ind = factor(ifelse(gender=='male',1,0))) %>%
  select(-ID, -gender, -group, -education)
dat_predict$center <- factor(dat_predict$center)

# plot mosaic
dat_ml_ft <- dat_ml_raw %>% mutate(gender = ifelse(gender == 'male','M','F'))
p_mosaic <- vcd::mosaic(~ center + group +gender , data = dat_ml_ft,
            main = "Sampling in Sites", shade = TRUE, legend = TRUE)

Cairo::Cairo(height = 820, file = 'output/REST_demographic_Chisquare_mosaic.png')
vcd::mosaic(~ center + group +gender , data = dat_ml_ft,
            main = "Sampling in Sites", shade = TRUE, legend = TRUE)
dev.off()

## load the pipelines
source('Analyze_script/function_RESTMachineLearningPipeline.R')

############################################################
#
# Randomforest - leave one site off
#
#############################################################

print('Running LOSOCV pipeline')

use_site <- unique(dat_ml_raw$center)

#plan(multisession)
result_random <- map_dfr(use_site,
                         ~RandomLOSOCVpipeline(data = dat_predict, .x),
                         .progress = T)

write_csv(result_random, 'output/REST_Machinelearning_LOSOCVResult.csv')
date()

############################################################
#
# Randomforest - random 2-fold CV
#
#############################################################
print('Running randomly assign pipeline')
dat_predict <- dat_predict %>% select(-center)

# choose 500 random seeds
seeds <- sample(1:1000, 500)
print("Start Randomly assign machine learning pipeline")
date()
#plan(multisession)
result_random <- map_dfr(seeds,
                                ~Random2FCVpipeline(data = dat_predict, .x), .progress = T)

write_csv(result_random, 'output/REST_Machinelearning_randomlyResult.csv')
print("All finished!")
date()

############################################################
#
# Result visualization
#
############################################################

result_table <- rio::import('output/REST_Machinelearning_randomlyResult.csv')

mean(result_table$Accuracy)
psych::describe(result_table$Accuracy)
psych::describe(result_table$Sensitivity)
psych::describe(result_table$Specificity)

p_acc <- ggplot(result_table) +
  geom_histogram(aes(x = Accuracy), fill = 'steelblue', color = 'white', bins = 20) +
  ylab('Frequency') + xlab('value') + ggtitle('Accuracy distribution') +
  theme_classic() + easy_text_size(15) + easy_center_title()

p_spe <- ggplot(result_table) +
  geom_histogram(aes(x = Specificity), fill = 'steelblue', color = 'white', bins = 20) +
  ylab('Frequency') + xlab('value') + ggtitle('Sensitivity distribution') +
  theme_classic() + easy_text_size(15)  + easy_center_title()

p_sen <- ggplot(result_table) +
  geom_histogram(aes(x = Sensitivity), fill = 'steelblue', color = 'white', bins = 20) +
  ylab('Frequency') + xlab('value') + ggtitle('Specificity distribution') +
  theme_classic() + easy_text_size(15) + easy_center_title()

result_table <- rio::import('output/REST_Machinelearning_LOSOCVResult.csv')
result_table <- result_table %>% select(Accuracy, Sensitivity, Specificity, Kappa) %>%
  summarise_each(funs(mean)) %>%
  data.table::melt() %>% mutate(value = round(value, digits = 3)) 
result_table$variable <- factor(result_table$variable, levels = c('Accuracy', 'Kappa', 'Sensitivity','Specificity'))

p_loocv <- ggplot(result_table, aes(x = variable, y = value)) +
  geom_bar(stat = 'identity', width = .5,
           fill = '#e27283', alpha =.8) +
  geom_text(aes(label = value),vjust = -0.5) +
  ggtitle('Model performance in LOSCV') +
  ylim(0, .7) +
  theme_classic() + easy_text_size(15) +
  easy_remove_x_axis(what = 'title')+
  easy_center_title() +
  easy_rotate_x_labels(angle = -30)


p_loocv + p_acc + p_sen + p_spe

ggsave(filename = 'output/REST_MachineLearning_result.png', height = 6, width = 9)

# Feature importance
result_table <- rio::import('output/REST_Machinelearning_randomlyResult.csv')
feaImp <- result_table[,c(6:33)] %>%
  summarise_each(funs(mean)) %>%
  melt()

wordcloud(feaImp$variable, feaImp$value, 
          colors = hcl.colors(30, palette = "Set2"),
          scale = c(3,0.5), rot.per = .4)

result_table <- rio::import('output/REST_Machinelearning_LOSOCVResult.csv')
feaImp <- result_table[,c(6:33)] %>%
  summarise_each(funs(mean)) %>%
  melt()

wordcloud(feaImp$variable, feaImp$value, 
          colors = hcl.colors(30, palette = "Set2"),
          scale = c(3,0.5), rot.per = .4)
