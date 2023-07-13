#' Function to plot the heat map for distance.
#'
#' This function takes a number as input and returns its square.
#'
#' @param cor_mat Distance dataframe estimated in Analyze3.
#' @param dat_collect A dataframe contains subject name and network identification.
#' @param colorbar A vector of color to use for heatmap's colorbar.
#' @param title Title of the heatmap.
#' @param outFileName Where to generate the heatmap.
#' @param breaks speficy the data scale in heatmap.


HeatMapforDistance <- function(cor_mat, dat_collect, colorbar, title, outFileName, breaks = NULL,
                               subject = 'sub-004') {
  library(tidyverse)
  library(pheatmap)
  library(viridis)
  
  cor_mat['to_ind'] <- str_extract(cor_mat$to, pattern = '[0-9]{1,3}') %>% as.numeric()
  cor_mat['from_ind'] <- str_extract(cor_mat$from, pattern = '[0-9]{1,3}') %>% as.numeric()
  cor_mat <- cor_mat %>% arrange(from_ind, to_ind) %>%
    mutate(from = paste0('parcel_', formatC(from_ind, width = 3, flag = "0"))) %>%
    mutate(to = paste0('parcel_', formatC(to_ind, width = 3, flag = "0")))
  # arrange by network
  annot_row <- data.frame('Network'=dat_collect %>% filter(participant_id == subject) %>% .$network)
  
  test <- cor_mat[,c(-4,-5)] %>% spread(to, value = stpdist) 
  cor_mat <- test[,c(-1)]
  cor_mat <- as.matrix(cor_mat)
  
  rownames(cor_mat) <- colnames(cor_mat)
  rownames(annot_row) <- colnames(cor_mat)
  
  annot_col <- data.frame('Network'=dat_collect %>% filter(participant_id == subject) %>% .$network)
  rownames(annot_col) <- colnames(cor_mat)
  
  net_color <- list(Network = c(DMN="#e27283", control="#edaf5e", limbic="#e7edca",dorsalAttn="#2c8b4b",
                                salience="#b77fb4", visual="#a153a2", somMot = "#6fabd2"))
  
  if (!is.null(breaks)) {
    pheatmap(cor_mat, 
             annotation_row = annot_row,
             annotation_col = annot_col,
             annotation_colors = net_color,
             breaks = breaks,
             color = colorbar,
             show_rownames = F, show_colnames = F,
             cluster_rows = F, cluster_cols = F,
             main = title,
             filename = outFileName,
             width = 12, height = 10)
  }else{
    pheatmap(cor_mat, 
             annotation_row = annot_row,
             annotation_col = annot_col,
             annotation_colors = net_color,
             color = colorbar,
             show_rownames = F, show_colnames = F,
             cluster_rows = F, cluster_cols = F,
             main = title,
             filename = outFileName,
             width = 12, height = 10)
  }
  
}

