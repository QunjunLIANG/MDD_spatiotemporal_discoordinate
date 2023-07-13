####################
#
# This script is used to plot the TD topology
# Heatmap for distance and shape the 
# data for further machine learning pipeline
#
# Liang Qunjun 2023/03/19


library(tidyverse)
library(bruceR)
library(RColorBrewer)
library(ggsci)
library(ggeasy)
library(patchwork)
library(aplot)
library(pheatmap)
library(viridis)
library(ggstatsplot)
library(ggridges)
library(lmerTest)
library(rio)


# indicate the subject name -------------------------------------------------

dat_collect <- import('output/validation_TD_topo_Gradient_collect_long.csv')
graddist_collect <- import('output/validation_TDtopology_distance_collection.csv')

## load the distance estimation functions
source('Analyze_script/function_HeatMapforDistance.R')


######################################################################
# Heat map of the average distance between MDD and HC
#
########################################################################
# MDD STP distance plot matrix 

cor_mat_mdd <- graddist_collect %>% filter(group == 'MDD') %>% group_by(from, to) %>%
  summarize(stpdist = mean(stpdist)) 

HeatMapforDistance(cor_mat = cor_mat_mdd, dat_collect = dat_collect, subject = 'sub-007', 
                   colorbar = colorpanel(35,"steelblue","white","darkred"),
                   breaks = seq(0,3.5,0.1),
                   title = "TD-based gradient distance in MDD sub-dataset", 
                   outFileName = "output/validation_Heatmap_TDtopo_distance_MDD.png")

## obtain the within-network value
model1_dat <- graddist_collect %>% filter(from_net == to_net) %>% 
  group_by(from_net, participant_id) %>%
  summarize_each(funs = funs(mean))

# preparing for machine learning --------------------------------------
## reshape the data
dat_within <- model1_dat %>% select(-from,-to,-to_net) %>%
  spread(key = from_net, value = stpdist)
colnames(dat_within)[3:9] <- paste0('intra_',colnames(dat_within)[3:9])

## slice the STP distance to fit mutual network distance
dat_interDMN <- graddist_collect %>% filter(from_net == 'DMN', to_net != 'DMN') %>%
  group_by(group, participant_id, to_net) %>% summarize_each(funs = funs(mean)) %>% 
  mutate(from_net = 'DMN')
dat_interFPN <- graddist_collect %>% filter(from_net == 'control',
                                            to_net %in% c("visual","somMot","dorsalAttn","salience","limbic")) %>%
  group_by(group, participant_id, to_net) %>% summarize_each(funs = funs(mean)) %>% 
  mutate(from_net = 'control')
dat_interLimb <- graddist_collect %>% filter(from_net == 'limbic',
                                             to_net %in% c("visual","somMot","dorsalAttn","salience")) %>%
  group_by(group, participant_id, to_net) %>% summarize_each(funs = funs(mean)) %>% 
  mutate(from_net = 'limbic')
dat_interSal <- graddist_collect %>% filter(from_net == 'salience',
                                            to_net %in% c("visual","somMot","dorsalAttn")) %>%
  group_by(group, participant_id, to_net) %>% summarize_each(funs = funs(mean)) %>% 
  mutate(from_net = 'salience')
dat_interDosAtt <- graddist_collect %>% filter(from_net == 'dorsalAttn',
                                               to_net %in% c("visual","somMot")) %>%
  group_by(group, participant_id, to_net) %>% summarize_each(funs = funs(mean)) %>% 
  mutate(from_net = 'dorsalAttn')
dat_interSoM <- graddist_collect %>% filter(from_net == 'somMot',
                                            to_net %in% c("visual")) %>%
  group_by(group, participant_id, to_net) %>% summarize_each(funs = funs(mean)) %>% 
  mutate(from_net = 'somMot')

dat_mutNet <- rbind(dat_interDMN, dat_interFPN, dat_interLimb, 
                    dat_interDosAtt, dat_interSal, dat_interSoM)
dat_mutNet <- dat_mutNet %>% select(-from, -to) %>% 
  pivot_wider(names_from = c(from_net, to_net), values_from = stpdist)


dat_logsitic <- dat_within[,-2] %>% left_join(dat_mutNet[,-1])

rio::export(dat_logsitic, file = 'output/validation_TDtopoDistance_machineLearning.csv')
