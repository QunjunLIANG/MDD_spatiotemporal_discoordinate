####################################
#
# Checking the reproducibility of the 
# spatiotemporal topology
#

library(tidyverse)
library(rio)
library(bruceR)
library(ggeasy)
library(lattice)
library(rasterVis)
library(raster)

# indicate the path
tdTopo_path <- '/home/lqj/MDD_patient/NIfTI_convert/derivatives/REST/SPT_topology'
sbj_file <- 'inputs/REST_match_data.xlsx'
hc_file <- 'inputs/REST_match_data_HC.tsv'
mdd_file <- 'inputs/REST_match_data_MDD.tsv'

# load the subject list
sbj_list <- rio::import(sbj_file)
hc_list <- rio::import(hc_file)
mdd_list <- rio::import(mdd_file)

# sample 100 subjects from HC and MDD
source('Analyze_script/function_RESTSPTtopoogy_reliability.R')

#########################################
#
# Compare with Main dataset
#
##########################################

# load SPT in REST
spt_rest_list <- list.files(tdTopo_path, pattern = '_spt_distance.csv', full.names = T)
spt_rest_dat <- read_csv(spt_rest_list, col_types = 'f')
spt_rest_dat_mean <- spt_rest_dat %>% group_by(parcel) %>%summarise_all(mean) %>% 
  dplyr::select(paste0('parcel',1:214)) %>% as.matrix()
# load SPT in Main dataset
spt_main_list <- list.files('output/SPT_topology_Power264/', pattern = '_spt_distance.csv', full.names = T)
spt_main_dat <- read_csv(spt_main_list, col_types = 'f')
spt_main_dat_mean <- spt_main_dat %>% group_by(parcel) %>%summarise_all(mean) %>% 
  dplyr::select(paste0('parcel',1:214)) %>% as.matrix()

rlist <- list(raster(spt_rest_dat_mean), raster(spt_main_dat_mean))
names(rlist) <- c('REST dataset', 'Main dataset')

p <- levelplot(stack(rlist), layout = c(2,1),
               at = seq(0, 3, .05))
p

Cairo::Cairo(file = 'output/REST_Main_Power264_matrix.png')
print(p)
dev.off()
# correlation between two dataset 
cor.test(spt_main_dat_mean, spt_rest_dat_mean)
# Pearson's product-moment correlation
# 
# data:  spt_main_dat_mean and spt_rest_dat_mean
# t = 316.16, df = 45794, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.8252362 0.8309917
# sample estimates:
#       cor 
# 0.8281358 

# visualization 
spt_rest_mean_low <- spt_rest_dat_mean[lower.tri(spt_rest_dat_mean)]
spt_main_mean_low <- spt_main_dat_mean[lower.tri(spt_main_dat_mean)]
spt_main_dataFrame <- data.frame('con' = factor(1:length(spt_main_mean_low)),
                                 'Main' = spt_main_mean_low)
spt_rest_dataFrame <- data.frame('con' = factor(1:length(spt_main_mean_low)),
                                 'REST' = spt_rest_mean_low)
spt_vis <- left_join(spt_main_dataFrame, spt_rest_dataFrame)

p_dotplot <- ggplot(spt_vis, aes(x = Main, y = REST)) +
  geom_point(size = 3, alpha = .15, 
             shape = 4, fill = 'black') +
  geom_smooth(method = 'lm') +
  xlab('Main dataset') + ylab('Validation dataset') +
  theme_classic(base_size = 20) +
  easy_center_title()
p_dotplot

Cairo::Cairo(file = 'output/REST_Main_Power264_dot.png')
print(p_dotplot)
dev.off()

ggsave('output/REST_Main_Power264_dot.png', height = 4,width = 5)

#########################################
#
# Sample size method
#
##########################################

hc_select <- sample(size = length(hc_list$ID), replace = F, hc_list$ID)
hc_compare <- hc_select[(length(hc_select)/2+1):length(hc_select)]
## the first half of hc data
hc_compare_list <- as.data.frame(list(hc_compare), col.names = 'ID') %>% 
  mutate(dir = paste0(tdTopo_path,'/',ID,'_spt_distance.csv'))
hc_comp_dat <- hc_compare_list$dir %>%
  map_dfr(~read_csv(.x, show_col_types = F), .progress=T) 
hc_comp_dat$parcel <- factor(hc_comp_dat$parcel, levels = paste0('parcel',1:214))
### group mean for the first half
hc_comp_dat_mean <- hc_comp_dat %>% group_by(parcel) %>%summarise_all(mean) %>% 
  select(paste0('parcel',1:214)) %>% as.matrix()
lattice::levelplot(hc_comp_dat_mean)

result_relia <- map_dfr(c(5:91),
                        ~ STPtopologyReliability_upsample(hc_select = hc_select,
                                                                   nsubject = .x, 
                                                                   hc_comp_dat_mean = hc_comp_dat_mean,
                                                                   tdTopo_path = tdTopo_path),
                        .progress =T)


result_relia2 <- result_relia %>% mutate(step = 1:nrow(result_relia))

p_relia <- ggplot(result_relia2, aes(x = step, y = cor, group = 1)) +
  geom_point(color = 'black', size = 1.5, alpha = 0.7) +
  stat_smooth() +
  geom_hline(yintercept = .9, linetype = 2, color = 'grey30')+
  geom_vline(xintercept = 21, linetype = 2, color = 'grey30')+
  scale_x_continuous(breaks = seq(5,nrow(result_relia),5)) +
  xlab('number of subjects') + ylab('correlation') +
  theme_classic() +
  easy_text_size(20)
ggsave(filename = 'output/SPTtopology_reliablity.png',
       width = 7.89, height = 4.95)

Cairo::Cairo(file = 'output/SPTtopology_reliablity.png')
print(p_relia)
dev.off()
