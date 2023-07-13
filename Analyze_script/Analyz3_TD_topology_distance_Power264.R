###############################
#
# This script used to calculate the 
# spatiotemporal topology based on TD
# gradient

library(tidyverse)
library(rio)
library(bruceR)

#########################################################
# Load the subject data
#
#########################################################

info_dat <- rio::import('inputs/MDD_HC_match_data.xlsx')

dir_tdg <- '/home/lqj/MDD_patient/NIfTI_convert/derivatives/timeseries_for_TD_Power264/time_lag_estimation/TDEGradient_alignment'
dir_fg <- '/home/lqj/MDD_patient/NIfTI_convert/derivatives/timeseries_for_gradient_Power264/individual_gradient'
dir_out <- 'output/SPT_topology_Power264'
## load the distance estimation functions
source('Analyze_script/function_EstimateDistance.R')

# create the output directory
if (!dir.exists(dir_out)) {
  print('create output directory')
  dir.create(dir_out, recursive = T)
}

# estimate the spatiotemporal topology and
# obtain the group-mean of the gradients in HC and MDD, respectively

## initial data frame for group average 
grp_mean_hc <- matrix(0, nrow = 214, ncol = 3) %>% as.data.frame()
grp_mean_mdd <- matrix(0, nrow = 214, ncol = 3) %>% as.data.frame()

colnames(grp_mean_hc) <- c('gradient1','gradient2','td_topo')
colnames(grp_mean_mdd) <- c('gradient1','gradient2','td_topo')

## main work
for (i in 1:nrow(info_dat)) {
  id <- info_dat$participant_id[i]
  group_tmp <- info_dat$group[i]
  Print("Estimate SPT topology in {id}")
  dat_tdg <- rio::import(file.path(dir_tdg,Glue('{id}_gradient_kernel-normalized_angle_method-dm_align-procrustes.csv')))
  dat_fg <- rio::import(file.path(dir_fg,Glue('{id}_gradient_export.csv')))
  this_data <- data.frame(gradient1 = dat_fg$gradient1, 
                          gradient2 = dat_fg$gradient2,
                          td_topo = dat_tdg$gradient1)
  dat_tdtopo <- Estimate_TDtopologyDistance(this_data, n_parcel = 214)
  if (group_tmp=='HC') {
    grp_mean_hc <- (grp_mean_hc + this_data)/2
  }else{
    grp_mean_mdd <- (grp_mean_mdd + this_data)/2
  }
  rio::export(dat_tdtopo, file = file.path(dir_out, Glue('{id}_spt_distance.csv')))
}
## export the group mean of STP topology
rio::export(grp_mean_hc, file = file.path(dir_out,'HC_group_mean.csv'))
rio::export(grp_mean_mdd, file = file.path(dir_out,'MDD_group_mean.csv'))

# visualization
## obtain the color representations of two group of subjects
isfahan <- MetBrewer::met.brewer("Isfahan1")
length(isfahan)
isfahan[1]

hc_color <- colorspace::lighten(isfahan[2], 0.35)
mdd_color <-  colorspace::lighten(isfahan[6], 0.35)

## load network identity
net_id <- import('output/Power2011ROI_clean_filtered.xlsx') %>% 
  filter(! network %in% c('Subcortical','Cerebellar','Memory','Uncertain'))
net_id$network <- factor(net_id$network, levels = c('VIS', 'somMot', 
                                                    'VAN', 'DAN', 'FPC', 
                                                    'Limbic', 'DMN'))

## network color 
net_color <- c("#a153a2","#6fabd2","#2c8b4b","#b77fb4","#e7edca","#edaf5e","#e27283")

## Plot in 3D -----------------------------------------------------------------
library(plotly)
library(ggstatsplot)
ggscatterstats(data = grp_mean_hc,
               x = gradient2, y = td_topo)

dat_tmp <- rbind(grp_mean_hc, grp_mean_mdd) %>% mutate(group = rep(c('HC','MDD'), each = 214))
plot_ly(data = dat_tmp,
        x = ~ gradient1, 
        y = ~ gradient2,
        z = ~ td_topo,
        color = ~group,
        colors = c(hc_color, mdd_color),
        marker = list(size = 6, line = list(color = 'grey',
                                            width = 2)) 
) %>%
  layout(scene = list(xaxis = list(title = 'functional gradient 1'),
                      yaxis = list(title = 'functional gradient 2'),
                      zaxis = list(title = 'TD-based gradient'))
  )


dat_hc <- grp_mean_hc %>% mutate(network = net_id$network)
plot_ly(data = dat_hc,
        x = ~ gradient1, 
        y = ~ gradient2,
        z = ~ td_topo,
        color = ~network,
        colors = net_color,
        marker = list(size = 6, line = list(color = 'grey',
                                            width = 2)) 
) %>%
  layout(scene = list(xaxis = list(title = 'functional gradient 1'),
                      yaxis = list(title = 'functional gradient 2'),
                      zaxis = list(title = 'TD-based gradient'),
                      title = 'HC')
  )

dat_mdd <- grp_mean_mdd %>% mutate(network = net_id$network)
plot_ly(data = dat_mdd,
        x = ~ gradient1, 
        y = ~ gradient2,
        z = ~ td_topo,
        color = ~network,
        colors = net_color,
        marker = list(size = 6, line = list(color = 'grey',
                                            width = 2)) 
) %>%
  layout(scene = list(xaxis = list(title = 'functional gradient 1'),
                      yaxis = list(title = 'functional gradient 2'),
                      zaxis = list(title = 'TD-based gradient'))
  )

# plot in heatmap -------------------------------------------------------
hc_grp <- file.path(dir_out,'HC_group_mean.csv') %>% import()
mdd_grp <- file.path(dir_out,'MDD_group_mean.csv') %>% import()

# estimte the group-mean distance
hc_grp_distance <- Estimate_TDtopologyDistance(hc_grp, n_parcel = 214)
mdd_grp_distance <- Estimate_TDtopologyDistance(mdd_grp, n_parcel = 214)

# plot the distance in HC-----------------------------------------------
cor_mat_hc <- hc_grp_distance[,c(-1)]
cor_mat_hc <- as.matrix(cor_mat_hc)

rownames(cor_mat_hc) <- colnames(cor_mat_hc)
annot_row <- data.frame('Network'=net_id$network)
rownames(annot_row) <- colnames(cor_mat_hc)

annot_row$Network <- factor(annot_row$Network, levels = c('somMot', 'VIS', 
                                                          'VAN', 'DAN', 'FPC', 
                                                          'Limbic', 'DMN'))
net_color <- list(Network = c(DMN="#e27283", FPC="#edaf5e", Limbic="#e7edca",DAN="#2c8b4b",
                              VAN="#b77fb4", VIS="#a153a2", somMot = "#6fabd2"))
library(pheatmap)
pheatmap(cor_mat_hc, 
         annotation_row = annot_row,
         annotation_col = annot_row,
         annotation_colors = net_color,
         border_color = NA,
         breaks = seq(0,6.5,0.05),
         color = gplots::colorpanel(120,"steelblue","white","darkred"),
         show_rownames = F, show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         main = 'REST HC SPT topology',
         filename = 'output/Main_Power264_hc_spttopo_distance.png',
         width = 12, height = 10)

# plot the distance in MDD -----------------------------------------------
cor_mat_mdd <- mdd_grp_distance[,c(-1)]
cor_mat_mdd <- as.matrix(cor_mat_mdd)

rownames(cor_mat_mdd) <- colnames(cor_mat_mdd)

pheatmap(cor_mat_mdd, 
         annotation_row = annot_row,
         annotation_col = annot_row,
         annotation_colors = net_color,
         border_color = NA,
         breaks = seq(0,6.5,0.05),
         color = gplots::colorpanel(120,"steelblue","white","darkred"),
         show_rownames = F, show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         main = 'REST MDD SPT topology',
         filename = 'output/Main_Power264_mdd_spttopo_distance.png',
         width = 12, height = 10)

# plot the distance for MDD>HC -----------------------------------------------
cor_mat_gap <- cor_mat_mdd - cor_mat_hc

pheatmap(cor_mat_gap, 
         annotation_row = annot_row,
         annotation_col = annot_row,
         annotation_colors = net_color,
         border_color = NA,
         # breaks = breaks
         color = gplots::colorpanel(35,"steelblue","white","darkred"),
         show_rownames = F, show_colnames = F,
         cluster_rows = F, cluster_cols = F,
         main = 'REST SPT topology (MDD>HC)',
         filename = 'output/Main_Power264_gap_spttopo_distance.png',
         width = 12, height = 10)
