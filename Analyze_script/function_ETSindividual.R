################################
#
# This function is used to estimte 
# the ETS based on Faskowitz, 2020, nature neuroscience
#
# Liang Qunjun 2023/06/05


EST_individual <- function(file_name, standard=FALSE){
  #############
  #
  # Computation for co-fluctuation 
  # netwrok
  
  library(tidyverse)
  library(caret)
  
  # obtain the subject name 
  sbjname <- str_extract(file_name, pattern = 'sub-[0-9]*')
  print(paste0("start working in ", sbjname))
  
  # load the raw signal extracted over Schafer-400
  bold_ts <- read_csv(file_name, col_names = F, show_col_types = FALSE) # delet the first 5 volumns
  
  # z-score for each region
  bold_ts_z <- bold_ts
  if (standard) {
    z_method <- bold_ts %>% preProcess(method = c('center','scale'))
    bold_ts_z <- predict(z_method, bold_ts)
  }
  bold_ts_z <- t(as.matrix(bold_ts_z)) # transform the dataframe to matrix
  
  # calculate element-wise product
  element_prod <- matrix(data = NA, nrow = 79800, ncol = 235)
  
  pos <- 1
  for (i in 1:nrow(bold_ts_z)) {
    for (j in i:nrow(bold_ts_z)) {
      if (i == j) {
        next
      }
      ets_tmp <- bold_ts_z[i,] * bold_ts_z[j,]
      element_prod[pos,] <- ets_tmp
      pos <- pos+1
    }
  }
  
  dim(element_prod)
  
  # calculate whole-brain co-fluctuation
  rss <- apply(element_prod, 2, function(x) {
    sqrt(sum(x^2))
  })
  
  # find the peak and trough
  t_peak <- c()
  t_trough <- c()
  
  for (i in 2:(length(rss)-1)) {
    if (rss[i] > rss[i-1] & rss[i] > rss[i+1]) {
      t_peak <- c(t_peak, i)
    }
    if (rss[i] < rss[i-1] & rss[i] < rss[i+1]) {
      t_trough <- c(t_trough, i)
    }
  }
  
  # calculate the duration of peak-trough and their coefficient of variation
  peak_amplitude <- mean(rss[t_peak])
  trough_duration <- mean(diff(t_trough))
  cv_amplitude <- sd(rss[t_peak])/mean(rss[t_peak])
  cv_duration <- sd(diff(t_trough))/mean(diff(t_trough))
  
  # export the result
  rss_df <- data.frame(participant_id = sbjname, 
                       peak_amplitude = peak_amplitude,
                       trough_duration = trough_duration,
                       cv_amplitude = cv_amplitude,
                       cv_duration = cv_duration,
                       as.data.frame(t(rss)))
  
  return(rss_df)
}