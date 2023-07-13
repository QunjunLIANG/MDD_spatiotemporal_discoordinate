
ConstructNetDistance <- function(dir_topo, sbjID){
  library(tidyverse)
  # load data
  dat_raw <- paste0(sbjID,'_spt_distance.csv') %>%
    file.path(dir_topo, .) %>% rio::import()
  
  dat_tmp <- melt(dat_raw, variable.name = 'toParcel', value.name = 'distance')
  colnames(dat_tmp)[1] <- 'fromParcel'
  dat_tmp$fromParcel <- factor(dat_tmp$fromParcel, levels = paste0('parcel',1:nrow(net_id)))
  dat_tmp$toParcel <- factor(dat_tmp$toParcel, levels = paste0('parcel',1:nrow(net_id)))
  dat_tmp <- dat_tmp %>% left_join(fromNet) %>% left_join(toNet)
  
  # summarize the distance based on network identification
  dat_tmp_net <- dat_tmp %>% filter(fromParcel != toParcel) %>%
    group_by(fromNet, toNet) %>%
    summarise(mean_dis=mean(distance))
  dat_tmp_net <- dat_tmp_net %>% rowwise() %>%
    mutate(link = paste(sort(c(fromNet,toNet)),collapse = '_'))
  result <- aggregate(mean_dis ~ link, data = dat_tmp_net, mean)
  result <- result %>% pivot_wider(names_from = link, values_from = mean_dis) %>%
    mutate(ID=sbjID)
  
  return(result)
}