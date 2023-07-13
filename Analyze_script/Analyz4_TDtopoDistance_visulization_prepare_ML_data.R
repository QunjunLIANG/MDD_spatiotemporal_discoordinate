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
library(ggcorrplot)
library(rasterVis)
library(raster)

# indicate the subject name -------------------------------------------------

dat_collect <- import('output/TD_topo_Gradient_collect_long.csv')
graddist_collect <- import('output/TDtopology_distance_collection.csv')

## load the distance estimation functions
source('Analyze_script/function_HeatMapforDistance.R')


######################################################################
# Heat map of the average distance between MDD and HC
#
########################################################################
# HC STP distance plot matrix 

cor_mat_hc <- graddist_collect %>% filter(group == 'HC') %>% group_by(from, to) %>%
  summarize(stpdist = mean(stpdist)) 





HeatMapforDistance(cor_mat = cor_mat_hc, dat_collect = dat_collect, colorbar = colorpanel(35,"steelblue","white","darkred"),
                   breaks = seq(0,3.5,0.1),
                   title = "Spatiotemporal distance in HC", outFileName = "output/Heatmap_TD_distance_HC.png")

# MDD STP distance plot matrix 

cor_mat_mdd <- graddist_collect %>% filter(group == 'MDD') %>% group_by(from, to) %>%
  summarize(stpdist = mean(stpdist)) 

HeatMapforDistance(cor_mat = cor_mat_mdd, dat_collect = dat_collect, colorbar = colorpanel(35,"steelblue","white","darkred"),
                   breaks = seq(0,3.5,0.1),
                   title = "Spatiotemporal distance distance in MDD", outFileName = "output/Heatmap_TD_distance_MDD.png")

# MDD - HC 
cor_mat_diff <- cor_mat_mdd
cor_mat_diff$stpdist <- cor_mat_mdd$stpdist - cor_mat_hc$stpdist

HeatMapforDistance(cor_mat = cor_mat_diff, dat_collect = dat_collect, colorbar = colorpanel(35,"steelblue","white","darkred"),
                   title = "Spatiotemporal distance distance for MDD > HC", outFileName = "output/Heatmap_TD_distance_difference.png")


# plot the matrix with Raster ----------------------------------------
mat_dir <- list.files('output/TDtopology_distance_individual/', full.names = T)

## load the MDD matrices
mat_mdd <- map(mat_dir[1:91], ~ read_csv(.x, col_names = F, skip = 1) %>% dplyr::select(-1) )
mat_mdd_mean <- mat_mdd %>% reduce(`+`) / length(mat_mdd) 

## load the HC matrices
mat_hc <- map(mat_dir[92:181], ~ read_csv(.x, col_names = F, skip = 1) %>% dplyr::select(-1) )
mat_hc_mean <- mat_hc %>% reduce(`+`) / length(mat_hc) 

##  plot the matrices with raster
myCol <- gplots::colorpanel(140, "steelblue", "white" ,"darkred")
rlist <- list(raster(as.matrix(mat_mdd_mean)), raster(as.matrix(mat_hc_mean)))
names(rlist) <- c('MDD SPT matrix', 'HC SPT matrix')
p <- levelplot(stack(rlist), layout = c(2,1), scales = list(cex = 1.5),
               at = seq(0, 3.5, .1), col.regions = myCol)
p

Cairo::CairoPNG(filename = 'output/MainData_SPTmatrices.png', width = 800, height = 600)
print(p)
dev.off()

myCol <- gplots::colorpanel(200, "blue", "white" ,"red")
rlist <- list(raster(as.matrix(mat_mdd_mean) - as.matrix(mat_hc_mean)))
p <- levelplot(stack(rlist), scales = list(cex = 2), 
               margin = FALSE,
               at = seq(-1, 1, .01), 
               col.regions = myCol, colorkey = list(space = 'right', labels = list(cex = 1.5)))
p
Cairo::CairoPNG(filename = 'output/MainData_SPTmatrices_GAP.png', width = 640, height = 480)
print(p)
dev.off()

# preparing for machine learning --------------------------------------
## reshape the data
dat_within <- model1_dat %>% select(-from,-to,-to_net) %>% spread(key = from_net, value = stpdist)
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
dat_mutNet <- dat_mutNet %>% select(-from, -to) %>% pivot_wider(names_from = c(from_net, to_net), values_from = stpdist)


dat_logsitic <- dat_within %>% left_join(dat_mutNet)

rio::export(dat_logsitic, file = 'output/Topology_Distance_machineLearning.csv')

#######################
#
# visualization
#
#####################

dat_ML <- rio::import('output/Topology_Distance_machineLearning.csv')

## Anova 
dat_ML %>% data.table::melt(id.vars = 1:2) %>%
  MANOVA(data = ., dv = 'value', between = c('group', 'variable'), file = 'output/ANOVA_result_spt_topology.doc') %>%
  EMMEANS(effect = 'group', by = 'variable')

## visualization with matrix
dat_vis <- dat_ML %>% data.table::melt(id.vars = 1:2) %>% group_by(group, variable) %>%
  summarise(mean_dist = mean(value)) %>% mutate(scale_dist = (mean_dist-mean(mean_dist))/sd(mean_dist))

Net_id <- dat_vis$variable %>% as.character()
Net_id <- strsplit(Net_id, split = '_') %>% unlist()
fromNet <- Net_id[seq(1,length(Net_id)-1,2)]
toNet <-  Net_id[seq(2,length(Net_id),2)]

dat_vis['fromNet_raw'] <- fromNet
dat_vis['toNet'] <- toNet

dat_vis <- dat_vis %>% mutate(fromNet = ifelse(fromNet_raw=='intra',toNet,fromNet_raw))
dat_vis$fromNet <- factor(dat_vis$fromNet, levels = c("somMot","visual",'dorsalAttn',"salience","control","limbic","DMN"))
dat_vis$toNet <- factor(dat_vis$toNet, levels = c("somMot","visual",'dorsalAttn',"salience","control","limbic","DMN"))

# building HC distance matrix

NetDistancePlot <- function(network_ind){
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
    easy_text_size(15) + easy_add_legend_title('Distance')
  
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

ggsave(filename = 'output/DistanceMat_spttopo_HC.png', plot = disMat_hc_plot,
       height = 4, width = 5)
ggsave(filename = 'output/DistanceMat_spttopo_MDD.png', plot = disMat_mdd_plot,
       height = 4, width = 5)

## difference in distance matrix 
disMat_diff <- disMat_mdd - disMat_hc
ggcorrplot(disMat_diff, type = 'lower', show.diag = T, legend.title = "Distance",
          ggtheme = theme_minimal(base_size = 20)) +
  scale_fill_gradient2(
    low = 'steelblue',
    mid = 'white',
    high = "darkred") +
  easy_text_size(15) + easy_add_legend_title('Distance')
ggsave(filename = 'output/DistanceMat_spttopo_MDD_HC.png',
       height = 4, width = 5)
