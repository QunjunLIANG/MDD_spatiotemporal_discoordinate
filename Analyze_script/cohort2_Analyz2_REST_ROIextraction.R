#######################################
#
# This script is used to extract the 
# timeseries from the selection subjects

library(tidyverse)
library(rio)
library(R.matlab)
library(bruceR)
library(purrr)

###############################################
#
# Prepare the Power264 brain network identity
#
###############################################

dat_power <- rio::import('output/Power2011ROI_Module-2hv15xc.xls',
                    col_names = c('parcel_id','net_id','net_color','net_name'))

dat_power <- dat_power %>% mutate(network = ifelse(net_name=='Default mode', 'DMN',
                                                   ifelse(net_name=='Sensory/somatomotor Hand', 'somMot',
                                                          ifelse(net_name=='Sensory/somatomotor Mouth', 'somMot',
                                                                 ifelse(net_name=='Cingulo-opercular Task Control', 'VAN',
                                                                        ifelse(net_name=='Fronto-parietal Task Control', 'FPC',
                                                                               ifelse(net_name=='Dorsal attention', 'DAN', 
                                                                                      ifelse(net_name=="Ventral attention", 'Limbic',
                                                                                             ifelse(net_name=='Auditory', 'VAN', 
                                                                                                    ifelse(net_name=='Salience', 'VAN',
                                                                                                           ifelse(net_name=='Visual', 'VIS', net_name)))))))))))
export(dat_power, 'output/Power2011ROI_clean.csv')
export(dat_power, 'output/Power2011ROI_clean.tsv')
export(dat_power, 'output/Power2011ROI_clean.xlsx')

## filter out the Cerebellar and Subcortical regions
dat_power_filtered <- dat_power %>% 
  filter(! network %in% c('Subcortical','Cerebellar','Memory','Uncertain'))

export(dat_power_filtered, 'output/Power2011ROI_clean_filtered.csv')
export(dat_power_filtered, 'output/Power2011ROI_clean_filtered.tsv')
export(dat_power_filtered, 'output/Power2011ROI_clean_filtered.xlsx')

# indicate which columns to use in REST signal matrix
power_filtered_index <- 1569 + dat_power_filtered$parcel_id

## export the net identity for time delay matrix
dat_power_filtered <- dat_power_filtered %>% mutate(net_id_new=ifelse(network=='somMot',1,
                                                                      ifelse(network=='VIS',2,
                                                                             ifelse(network=='VAN',3,
                                                                                    ifelse(network=='DAN',4,
                                                                                           ifelse(network=='FPC',5,
                                                                                                  ifelse(network=='Limbic',6,7)))))))

write_csv(dat_power_filtered$net_id_new %>% as.data.frame(), col_names = F,
          file ='output/Power264_net_identity.csv')

###############################################
#
# Prepare the Dosenbach160 brain network identity
#
#############################################

dat_dosen <- rio::import('output/Dosenbach160_net_identity.csv',
                         header = T)

# fix the parcelID
dat_dosen$V1 <- dat_dosen$V1 +1
colnames(dat_dosen)[1] <- 'parcel_id'

# 

#################################################
#
# Extract power264 signal with filtered indeces
#
#################################################

sbj_list <- import('REST_match_data.xlsx')
sbj_list <- sbj_list$ID

# extract the timeseries 
source('Analyze_script/function_RESTsignalROIoperate.R')
map(sbj_list, ~ SliceParcelREST(sbj_name = .x,
                              matDir = 'output/REST/ROISignals_FunImgARCWF/',
                              outDir = '/home/lqj/MDD_patient/NIfTI_convert/derivatives/REST/power264_timeseries',
                              col_range = power_filtered_index,
                              outname = 'Power264'),
    .progress = T)

map(sbj_list, ~ SliceParcelREST_csv(sbj_name = .x,
                                matDir = 'output/REST/ROISignals_FunImgARCWF/',
                                outDir = '/home/lqj/MDD_patient/NIfTI_convert/derivatives/REST/power264_timeseries',
                                col_range = power_filtered_index,
                                outname = 'Power264'),
    .progress = T)
