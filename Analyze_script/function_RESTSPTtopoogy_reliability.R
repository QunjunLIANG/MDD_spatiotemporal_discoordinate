
STPtopologyReliability_upsample <- function(hc_select, nsubject, hc_comp_dat_mean, tdTopo_path){
  
    print(Glue("this is the {nsubject} subject step"))
  
    hc_test <- hc_select[1:nsubject]
    hc_test_list <- as.data.frame(list(hc_test), col.names = 'ID') %>% 
      mutate(dir = paste0(tdTopo_path,'/',ID,'_spt_distance.csv'))
    
    hc_test_dat <- hc_test_list$dir %>%
      map_dfr(~read_csv(.x, show_col_types = F)) 
    hc_test_dat$parcel <- factor(hc_test_dat$parcel, 
                                 levels = paste0('parcel',1:214))
    
    hc_test_dat_mean <- hc_test_dat %>% group_by(parcel) %>% 
      summarise_all(mean) %>% 
      select(paste0('parcel',1:214)) %>% as.matrix()
    
    cor_tmp <- cor.test(hc_comp_dat_mean, hc_test_dat_mean, method = 'pearson')
    cor_tmp_nonpa <- cor.test(hc_comp_dat_mean, hc_test_dat_mean, method = 'spearman')
    
    dat_tmp <- data.frame(nsbj = nsubject,
                          cor = cor_tmp$estimate,
                          low_conf = cor_tmp$conf.int[1],
                          up_conf = cor_tmp$conf.int[2],
                          rho = cor_tmp_nonpa$estimate)
    return(dat_tmp)
}


SPTtopologyReliability_resample <- function(hc_list_select, mdd_list_select, tdTopo_path){
  
  hc_list <- as.data.frame(list(hc_list_select), col.names = 'ID')
  mdd_list <- as.data.frame(list(mdd_list_select), col.names = 'ID')
  
  ## obtain the hc data path
  hc_file_select <- hc_list %>% 
    mutate(dir = paste0(tdTopo_path,'/',ID,'_spt_distance.csv'))
  ## the first half of hc data
  hc_d1 <- hc_file_select$dir[1:91] %>%
    map_dfr(~read_csv(.x, show_col_types = F)) 
  hc_d1$parcel <- factor(hc_d1$parcel, levels = paste0('parcel',1:214))
  ### group mean for the first half
  hc_d1_mean <- hc_d1 %>% group_by(parcel) %>%summarise_all(mean) %>% 
    select(paste0('parcel',1:214)) %>% as.matrix()
  
  ## the second half of hc data
  hc_d2 <- hc_file_select$dir[92:182] %>%
    map_dfr(~read_csv(.x, show_col_types = F), .progress = T) 
  hc_d2$parcel <- factor(hc_d1$parcel, levels = paste0('parcel',1:214))
  ### group mean for the second half
  hc_d2_mean <- hc_d2 %>% group_by(parcel) %>%summarise_all(mean) %>% 
    select(paste0('parcel',1:214)) %>% as.matrix()
  # similarity in hc dataset
  cor_hc <- cor.test(hc_d1_mean, hc_d2_mean)
  
  ## mdd comparison 
  mdd_file_select <- mdd_list %>% 
    mutate(dir = paste0(tdTopo_path,'/',ID,'_spt_distance.csv'))
  mdd_d1 <- mdd_file_select$dir %>%
    map_dfr(~read_csv(.x, show_col_types = F)) 
  mdd_d1$parcel <- factor(mdd_d1$parcel, levels = paste0('parcel',1:214))
  mdd_d1_mean <- mdd_d1 %>% group_by(parcel) %>%summarise_all(mean) %>% 
    select(paste0('parcel',1:214)) %>% as.matrix()
  # similarity in hc dataset
  cor_hc_mdd <- cor.test(hc_d1_mean, mdd_d1_mean)
  
  return(data.frame(hc_hc_cor = cor_hc$estimate,
                    hc_low_conf = cor_hc$conf.int[1],
                    hc_up_conf = cor_hc$conf.int[2],
                    hc_mdd_cor = cor_hc_mdd$estimate,
                    mdd_low_conf = cor_hc_mdd$conf.int[1],
                    mdd_up_conf = cor_hc_mdd$conf.int[2]))
}