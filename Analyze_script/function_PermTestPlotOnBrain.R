#' Statistical test for mean differnce of the input feature.
#'
#' This function takes a number as input and returns its square.
#'
#' @param atlas A Schaefer 400 7 network.
#' @param dat_feature A dataframe contains features.
#' @param feature One of: "TDE", "gradient", "TDtopology".
#' @param hemi_position One of: "horizontal", "vertical", "montage".
#' @param plot_title A string title for the plot.
#' @param lengend_title A string for plot's legend.
#' @param p_adj One of "holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none"
#' @return A plot of brain surface.
#' @examples
#' NULL
 
PermTestPlotOnBrain <- function(atlas, dat_feature, hemi_position, plot_title, line_size=0,
                                lengend_title="Permutation t value", n_permu = 10000, p_adj = "fdr"){
  
  ## You need the following package(s) to generate the plot:
  library(ggseg)
  library(ggsegSchaefer)
  library(tidyverse)
  library(ggeasy)
  library(viridis)
  library(exactRankTests) 
  
  ## T-test for each parcel in gradient 1 -----------------------
  t_collect <- rep(NA, times = 400)
  p_collect <- rep(NA, times = 400)
  for (i in seq_along(1:400)) {
    dat_tmp <- dat_feature %>% filter(gradient == 'gradient1') %>% select(all_of(i), group)
    colnames(dat_tmp) <- c('value','group')
    result_tmp <- perm.test.formula(value ~ group, data = dat_tmp,  nperm = n_permu)
    t_collect[i] <- result_tmp$statistic
    p_collect[i] <- result_tmp$p.value
  }
  
  ## collect the result
  p_collect_adj <- p.adjust(p_collect, method = p_adj)
  which(p_collect_adj <= 0.05)
  which(p_collect <= 0.05)
  
  t_collect_p <- t_collect
  t_collect_p[which(p_collect > 0.05)] = 0
  
  t_collect_p_adj <- t_collect
  t_collect_p_adj[which(p_collect_adj > 0.05)] = 0
  
  stat_result <- data.frame(parcel = paste0('parcel', 1:length(p_collect)),
                            t = t_collect,
                            t_p = t_collect_p,
                            t_p_adj = t_collect_p_adj,
                            p = p_collect,
                            p_fdr = p_collect_adj
  )
  
  stat_result_uncor <- data.frame(parcel = paste0('parcel', 1:length(p_collect)),
                                  t = t_collect_p,
                                  p = p_collect
  )
  
  stat_result_cor <- data.frame(parcel = paste0('parcel', 1:length(p_collect)),
                                t = t_collect_p_adj,
                                p = p_collect_adj
  )
  
  ## plot the t-test result -----------------------------
  data_uncor <-  atlas
  data_uncor <- data_uncor %>% left_join(stat_result_uncor)
  data_uncor['group'] = 'uncorrectP'
  
  data_cor <-  atlas
  data_cor <- data_cor %>% left_join(stat_result_cor)
  data_cor['group'] = 'correctP'
  
  atlas_data = rbind(data_uncor, data_cor)
  
  ## plot the t value on brain surface
  brain_pos <- position_brain("horizontal")
  ## Plot atlas:
  p_this <- ggplot() + geom_brain(
    atlas       = atlas_data,
    mapping     = aes(fill=t),
    position    = brain_pos,
    color       ='black',
    size        = line_size,
    show.legend = T) +
    ggtitle(plot_title) +
    facet_wrap(~group, drop = F) +
    scale_fill_viridis_c(option= "H") +
    theme_void() + easy_legend_at(to = 'bottom') +
    easy_center_title() + easy_add_legend_title(lengend_title) +
    easy_text_size(which = c("plot.title","legend.text", "legend.title", "strip.text"),
                   size = 15)
  
  if (hemi_position == "vertical") {
    brain_pos <- position_brain(hemi + side ~ . )
    p_this <- ggplot() + geom_brain(
      atlas       = atlas_data,
      mapping     = aes(fill=t),
      position    = brain_pos,
      color       ='black',
      size        = line_size,
      show.legend = T) +
      ggtitle(plot_title) +
      facet_grid(.~group) +
      scale_fill_viridis_c(option= "H") +
      theme_void() + easy_legend_at(to = 'bottom') +
      easy_center_title() + easy_add_legend_title(lengend_title) +
      easy_text_size(which = c("plot.title","legend.text", "legend.title", "strip.text"),
                     size = 15)
  }else if(hemi_position == "montage"){
    brain_pos <- position_brain(side ~ hemi )
    p_this <- ggplot() + geom_brain(
      atlas       = atlas_data,
      mapping     = aes(fill=t),
      position    = brain_pos,
      color       ='black',
      size        = line_size,
      show.legend = T) +
      ggtitle(plot_title) +
      facet_grid(.~group) +
      scale_fill_viridis_c(option= "H") +
      theme_void() + easy_legend_at(to = 'bottom') +
      easy_center_title() + easy_add_legend_title(lengend_title) +
      easy_text_size(which = c("plot.title","legend.text", "legend.title", "strip.text"),
                     size = 15)
    }
  
  return(list(statistic = stat_result, plot = p_this))
}

