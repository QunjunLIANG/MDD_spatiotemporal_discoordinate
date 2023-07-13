# function to estimate STP 3D distance 
Estimate_TDtopologyDistance <- function(data, n_parcel=400) {
  
  library(tidyverse)
  data_scaled <- scale(data[,c("gradient1","gradient2","td_topo")])
  
  pariwise_dist <- as.matrix(dist(data_scaled))
  
  result_table <- data.frame(parcel = paste0('parcel', seq(nrow(data))))
  
  for (i in seq(nrow(data))) {
    result_table[,i+1] <- pariwise_dist[i,]
  }
  
  colnames(result_table)[2:ncol(result_table)] <- paste0('parcel',1:n_parcel)
  return(result_table)
}


# function to estimate STP 3D distance 
Estimate_Gradientdistance <- function(data, n_parcel=400) {
  
  library(tidyverse)
  data_scaled <- scale(data[,c("gradient1","gradient2","gradient3")])
  
  pariwise_dist <- as.matrix(dist(data_scaled))
  
  result_table <- data.frame(parcel = paste0('parcel', seq(nrow(data))))
  
  for (i in seq(nrow(data))) {
    result_table[,i+1] <- pariwise_dist[i,]
  }
  
  colnames(result_table)[2:ncol(result_table)] <- paste0('parcel',1:n_parcel)
  return(result_table)
}