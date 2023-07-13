##################################################
#
# Data selection for REST
#
# Liang Qunjun  2023/03/30

library(tidyverse)
library(rio)
library(VIM)
library(mice)
library(data.table)
library(table1)
library(tableone)
library(MatchIt)
library(purrr)

###############################################
#
# shape the demographic information
#
###############################################

# load the raw data of HC
dat_hc_raw <- import('REST-meta-MDD-PhenotypicData_WithHAMDSubItem_V4.xlsx',
                     sheet = 'Controls')
VIM::matrixplot(dat_hc_raw)
colnames(dat_hc_raw)

## shape the HC raw data 
dat_hc_mod <- dat_hc_raw %>% mutate(group = 'HC')
dat_hc_mod['center'] <- str_extract(dat_hc_mod$ID, pattern = 'S[0-9]{1,2}')
dat_hc_mod['subject'] <- paste('sub', str_extract(dat_hc_mod$ID, pattern = '[0-9]{4}$'), sep = '-')
dat_hc_mod <- dat_hc_mod %>% mutate(education = ifelse(`Education (years)` == 0, 1,
                                                       ifelse(`Education (years)` > 0 & `Education (years)` <= 6, 2,
                                                              ifelse(`Education (years)` > 6 & `Education (years)` <= 9, 3,
                                                                     ifelse(`Education (years)` > 9 & `Education (years)` <= 12, 4,
                                                                            ifelse(`Education (years)` > 12 & `Education (years)` <= 16, 5, 6))))))

## check the distribution of Age, education and gender
hist(dat_hc_mod$Age, breaks = 40)
table(dat_hc_mod$education)
table(dat_hc_mod$Sex)

# load the raw data of MDD
dat_mdd_raw <- import('REST-meta-MDD-PhenotypicData_WithHAMDSubItem_V4.xlsx',
                  sheet = 'MDD', na=c('',-9999,'[]'))
VIM::matrixplot(dat_mdd_raw)
colnames(dat_mdd_raw)

## shape the MDD raw data
dat_mdd_mod <- dat_mdd_raw %>% select(1:8)
colnames(dat_mdd_mod)


dat_mdd_mod['center'] <- str_extract(dat_mdd_mod$ID, pattern = 'S[0-9]{1,2}')
dat_mdd_mod['subject'] <- paste('sub', str_extract(dat_mdd_mod$ID, pattern = '[0-9]{4}$'), sep = '-')

dat_mdd_mod <- dat_mdd_mod %>% mutate(group = 'MDD') %>%
  mutate(education = ifelse(`Education (years)` == 0, 1,
                            ifelse(`Education (years)` > 0 & `Education (years)` <= 6, 2,
                                   ifelse(`Education (years)` > 6 & `Education (years)` <= 9, 3,
                                          ifelse(`Education (years)` > 9 & `Education (years)` <= 12, 4,
                                                 ifelse(`Education (years)` > 12 & `Education (years)` <= 16, 5, 6))))))
VIM::matrixplot(dat_mdd_mod)

# pre-filter the demographic data
dat_hc_mod_prefilter <- dat_hc_mod %>% mutate(prefilter = ifelse(Age >= 18, 1, 0))
dat_mdd_mod_prefilter <- dat_mdd_mod %>% mutate(prefilter = ifelse(Age >=18 & HAMD >= 17, 1, 0))

## combine the pre-filter data
dat_demog_prefilter <- data.table::rbindlist(list(dat_hc_mod_prefilter, dat_mdd_mod_prefilter),
                                             fill = T)

#######################################################################
#
# Combine the image quality control index to demegraphic
#
######################################################################

# load the image quality checklist
dat_img_quality <- import('output/REST/RawFunImgQC.tsv')

# load the FD and average within subject
fd_list <- list.files('output/REST/FD/', full.names = T)
fd_sbj_name <- str_extract(fd_list, pattern = 'S[0-9]{1,2}-[0-9]-[0-9]{4}')
## define a function to load fd
LoadFD <- function(path_to_fd){
  sbj_name <- str_extract(path_to_fd, pattern = 'S[0-9]{1,2}-[0-9]-[0-9]{4}')
  dat_fd_tmp <- read_csv(file = path_to_fd, col_names = 'fd', show_col_types = F) %>%
    summarise(fd_mean = mean(fd), timepoint = n(), fd_precent = mean(fd>0.2)) %>%
    mutate(ID = sbj_name)
  return(dat_fd_tmp)
}

dat_fd_raw <- map_dfr(fd_list, ~LoadFD(path_to_fd = .x), .progress = T)

# combine the image quality and fd to the main demographic data
colnames(dat_img_quality)[1] <- 'ID'
colnames(dat_fd_raw)[4] <- 'ID'
dat_demo_img_fd <- dat_demog_prefilter %>% left_join(dat_img_quality) %>% left_join(dat_fd_raw)
export(dat_demo_img_fd, file = 'output/REST_demog_imag_fd.xlsx')

#######################################################################
#
# Merge the ROI signal matrix dimension to main data, checking 
# the completely of the parcellations
#
######################################################################

# load the essential functions
source('Analyze_script/function_RESTsignalROIoperate.R')

sbj_list <- dat_demo_img_fd$ID
roiDir <- 'output/REST/ROISignals_FunImgARCWF/'

### please NOTE!!!! power_filtered_ind was obtain from Analyz2
roi_signal_check <- map_dfr(sbj_list, ~ RecordRoIdimensionREST(sbj_name = .x,
                                                               matDir = roiDir,
                                                               col_range = power_filtered_index),
                            .progress=T)

# merge the data into the main data 
dat_demo_img_fd_roicheck <- dat_demo_img_fd %>% left_join(roi_signal_check)
export(dat_demo_img_fd_roicheck, file = 'output/REST_demog_qc_fd_roicheck.xlsx')

######################################################################
#
# Filtering the data with the following principal:
#   1. age >= 18;
#   2. MDD patients with HAMD score > 17;
#   3. QC score > 2 (poor);
#   4. mean FD <= 0.2;
#   5. FD spike <= 30%
#   6. timepoints >= 190
#   7. ROI matrix = 1833 columns
#   8. No zero signal in Power264 parcellation
#
######################################################################

# filter the data to clean
dat_final_use <- dat_demo_img_fd_roicheck %>% filter(prefilter == 1, 
                                            `QC Score` > 2,
                                            fd_mean <= 0.2, 
                                            fd_precent <= 0.3,
                                            timepoints >= 190,
                                            Ncol == 1833,
                                            parcel_complete == 0)
export(dat_final_use, file = 'output/REST_filtered.xlsx')

# summary the data 
dat_final_use <- dat_final_use %>% mutate(gender = ifelse(Sex == 1, 'male', 'female'))
dat_final_use$education <- factor(dat_final_use$education, levels = c(1:6))
dat_final_use$center <- factor(dat_final_use$center, 
                               levels = unique(dat_final_use$center))
table1::table1(data = dat_final_use, 
               ~ gender + Age + education + center | group )

######################################################################
#
# Match the group according to gender, age and education
# using propensity score matching
#
######################################################################

dat_final_use$education <- as.numeric(dat_final_use$education)
confounds <- c('gender', 'Age', 'center')
table_ori <- tableone::CreateTableOne(vars = confounds,
                            strata = 'group',
                            data = dat_final_use, test = T)
print(table_ori)

dat_final_use <- dat_final_use %>% mutate(group_encode = ifelse(group == 'HC', 0, 1))
psm_matchit <- matchit(group_encode ~ gender + Age + center,
                       method = 'nearest', distance = 'glm',
                       ratio = 1, replace =F, caliper = .2,data = dat_final_use)

summary(psm_matchit, un = F)

plot(psm_matchit, type = 'jitter')
# output the match data
psm_matchit_data <- match.data(psm_matchit)
psm_matchit_data$education <- factor(psm_matchit_data$education)
table1::table1(data = psm_matchit_data, 
               ~ gender + Age + education + center | group )

# add participant_id
psm_matchit_data <- psm_matchit_data %>% 
  mutate(participant_id = paste0('sub-', 
                                 formatC(1:nrow(psm_matchit_data), 
                                               width = 4,
                                               flag = "0")))

export(psm_matchit_data, file = 'REST_match_data.xlsx')
psm_matchit_data %>% 
  select(ID) %>%
  write_tsv(file = 'REST_match_data_all.csv')

psm_matchit_data %>% filter(group == 'HC') %>% 
  select(ID) %>%
  write_tsv(file = 'REST_match_data_HC.tsv')

psm_matchit_data %>% filter(group == 'MDD') %>%
  select(ID) %>%
  write_tsv(file = 'REST_match_data_MDD.tsv')
