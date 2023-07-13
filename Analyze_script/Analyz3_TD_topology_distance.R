####################

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

#########################################################
# Load the subject data
#
#########################################################

info_dat <- import('MDD_HC_match_data.xlsx')

## if the file of dat_collect did not exists, please run Analyze1 script
dat_collect <- import('output/TD_topo_Gradient_collect_long.csv')

## load the distance estimation functions
source('Analyze_script/function_EstimateDistance.R')


######################################################################
# calculate the raw Gradient distance for each subject
#
########################################################################

# graddist_collect <- import('output/TDtopology_distance_collection.csv')

if (!exists('graddist_collect')) {
  dist_outDir <- 'output/TDtopology_distance_individual'
  if (!dir.exists(dist_outDir)) {
    dir.create(dist_outDir)
  }
  
  graddist_collect <- data.frame()
  for (i in info_dat$participant_id) {
    
    stpdist_tmp <- dat_collect %>% filter(participant_id == i) 
    net_from <- data.frame('from' = paste0('parcel',c(1:nrow(stpdist_tmp))), 'from_net'=stpdist_tmp$network)
    net_to <- data.frame('to' = paste0('parcel',c(1:nrow(stpdist_tmp))), 'to_net'=stpdist_tmp$network)
    
    stpdist_tmp <- Estimate_TDtopologyDistance(stpdist_tmp)
    export(stpdist_tmp, file = file.path(dist_outDir, paste0(i,'_TDEFC_distance.csv')))
    
    stpdist_tmp <- melt(stpdist_tmp)
    
    colnames(stpdist_tmp) <- c('from','to','stpdist')
    
    stpdist_tmp <- stpdist_tmp %>% left_join(net_from) %>% left_join(net_to)
    
    stpdist_tmp['participant_id'] <- i
    
    graddist_collect <- rbind(graddist_collect, stpdist_tmp)
  }
  
  graddist_collect <- graddist_collect %>% left_join(net_from) %>% left_join(net_to)
  graddist_collect <- graddist_collect %>% left_join(info_dat[,c(1,71)])
  
  export(graddist_collect, file = 'output/TDtopology_distance_collection.csv')
}

#######################################################
# visualization the TDE and gradient topology
#
########################################################
# color identification --------------------------------------------------------
#
# obtain the color representations of two group of subjects
isfahan <- MetBrewer::met.brewer("Isfahan1")
length(isfahan)
isfahan[1]

hc_color <- colorspace::lighten(isfahan[2], 0.35)
mdd_color <-  colorspace::lighten(isfahan[6], 0.35)

## network color 
net_color <- c("#a153a2","#6fabd2","#2c8b4b","#b77fb4","#e7edca","#edaf5e","#e27283")

## 3D plotly ------------------------------------------------------------------
dat_collect$network <- factor(dat_collect$network,  levels = c("visual", "somMot",'dorsalAttn',"salience","control","limbic","DMN"))
dat_tmp <- dat_collect %>% group_by(group, parcel) %>%
  summarise_each(funs(mean))

library(plotly)
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

dat_hc <- dat_collect %>% filter(group == 'HC') %>% group_by(parcel, network) %>%
  summarise_each(funs(mean))
dat_mdd <- dat_collect %>% filter(group == 'MDD') %>% group_by(parcel, network) %>%
  summarise_each(funs(mean))

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
                      zaxis = list(title = 'TD-based gradient'))
  )

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
