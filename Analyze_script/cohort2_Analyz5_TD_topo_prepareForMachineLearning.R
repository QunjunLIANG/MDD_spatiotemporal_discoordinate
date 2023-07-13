library(tidyverse)
library(rio)

# load data
info_dat <- import('REST_match_data.xlsx')
dir_topo <- '/home/lqj/MDD_patient/NIfTI_convert/derivatives/REST/SPT_topology'
## load network identity
net_id <- import('output/Power2011ROI_clean_filtered.xlsx') %>% 
  filter(! network %in% c('Subcortical','Cerebellar','Memory','Uncertain'))
net_id$network <- factor(net_id$network, levels = c('somMot', 'VIS', 
                                                    'VAN', 'DAN', 'FPC', 
                                                    'Limbic', 'DMN'))
fromNet <- data.frame('fromParcel'=paste0('parcel',1:nrow(net_id)),
                      'fromNet'=net_id$network)
toNet <- data.frame('toParcel'=paste0('parcel',1:nrow(net_id)),
                      'toNet'=net_id$network)

# obtain the network-based distance 
source('Analyze_script/function_RESTobtainNetDistance.R')
dat_netDist <- purrr::map_dfr(info_dat$ID, 
                              ~ ConstructNetDistance(dir_topo = dir_topo, .x),
                              .progress = T)
dat_netDist_demo <- info_dat %>% 
  select(ID, Age, gender, education, group, center, fd_mean) %>%
  left_join(dat_netDist)

rio::export(dat_netDist_demo, file = 'output/REST_spttopo_machineLearning.csv')

# statistics ------------------------------------------------------------------

dat_netDist_demo <- rio::import('output/REST_spttopo_machineLearning.csv')

library(bruceR)

## longer the data
dat_netDist_long <- dat_netDist_demo %>% 
  melt(id.vars = 1:7, 
       variable.name = 'link', value.name='distance')

## between group mean-difference analysis
bet_group <- ggstatsplot::grouped_ggbetweenstats(
  data = dat_netDist_long,
  x = group,
  y = distance,
  grouping.var = link,
  bf.message = F
)
ggsave(filename = 'output/REST_between_diff.png',
       plot = bet_group, height = 17, width = 25)


## divide the intra-network and inter-network link
intr_net_ind <- c('DAN_DAN', 'DMN_DMN', 'FPC_FPC', 
                  'Limbic_Limbic', 'somMot_somMot', 
                  'VAN_VAN', 'VIS_VIS')
dat_netDist_intra <- dat_netDist_long %>% filter(link %in% intr_net_ind)

bruceR::MANOVA(data = dat_netDist_intra,
               subID = 'ID', dv = 'distance', 
               between = 'group', within = 'link', covariate = 'gender',
               file = 'output/REST_spttopo_MANOVA.doc')
