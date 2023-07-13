#' Calculate the square of a number.
#'
#' This function takes a number as input and returns its square.
#'
#' @param atlas A Schaefer 400 7 network.
#' @param dat_feature A dataframe contains features.
#' @param feature One of: "TDE", "gradient", "TDtopology".
#' @param hemi_position One of: "horizontal", "vertical", "montage".
#' @param plot_title A string title for the plot.
#' @param lengend_title A string for plot's legend.
#' @return A plot of brain surface.
#' @examples
#' PlotFeatureOnSurface_group(atlas_data, dat_tdtopo, "TDtopology", "TD gradient in groups", "Gradient score")

PlotFeatureOnSurface_group <- function(atlas, dat_feature, feature, 
                                       hemi_position, plot_title, lengend_title, 
                                       gradient=NULL,
                                       line_size=0,
                                       show_legend = T) {
  
  library(ggseg)
  library(ggsegSchaefer)
  library(tidyverse)
  library(ggeasy)
  library(viridis)
  
  ## reshape the atlas
  atlas_data = atlas
  
  if (feature == "TDE") {
    feature_mean_mdd <- dat_feature %>% filter(group == 'MDD') %>% 
      summarise_each(funs(mean)) %>% select(-group, -participant_id) %>%
      t() %>% as.data.frame()
    feature_mean_hc <- dat_feature %>% filter(group == 'HC') %>% 
      summarise_each(funs(mean)) %>% select(-group, -participant_id) %>%
      t() %>% as.data.frame()
  }else if(feature == "gradient" || feature == "TDtopology"){
    feature_mean_mdd <- dat_feature %>% filter(group == 'MDD', gradient == 'gradient1') %>% 
      summarise_each(funs(mean)) %>%
      select(-group, -participant_id, -gradient) %>%
      t() %>% as.data.frame()
    feature_mean_hc <- dat_feature %>% filter(group == 'HC', gradient == 'gradient1') %>% 
      summarise_each(funs(mean)) %>% 
      select(-group, -participant_id, -gradient) %>%
      t() %>% as.data.frame()

  if(gradient=="grdient2"){
      feature_mean_mdd <- dat_feature %>% filter(group == 'MDD', gradient == 'gradient2') %>% 
        summarise_each(funs(mean)) %>%
        select(-group, -participant_id, -gradient) %>%
        t() %>% as.data.frame()
      feature_mean_hc <- dat_feature %>% filter(group == 'HC', gradient == 'gradient2') %>% 
        summarise_each(funs(mean)) %>% 
        select(-group, -participant_id, -gradient) %>%
        t() %>% as.data.frame()
    }

  }
  
  ## rename the colname
  feature_mean_mdd['parcel'] <- rownames(feature_mean_mdd)
  colnames(feature_mean_mdd)[1] <- 'value'
  feature_mean_hc['parcel'] <- rownames(feature_mean_hc)
  colnames(feature_mean_hc)[1] <- 'value'
  
  ## combine the feature to the atlas data
  data_mdd = atlas_data %>% left_join(feature_mean_mdd)
  data_mdd['group'] = 'MDD'
  data_hc = atlas_data %>% left_join(feature_mean_hc)
  data_hc['group'] = 'HC'
  
  ## combine mdd and hc data
  atlas_data = rbind(data_mdd, data_hc)
  
  brain_pos <- position_brain("horizontal")
  
  colorbar_width = .1
  ## Plot atlas:
  p_this <- ggplot() + geom_brain(
    atlas       = atlas_data,
    mapping     = aes(fill=value),
    position    = brain_pos,
    color       ='black',
    size        = line_size,
    show.legend = show_legend) +
    ggtitle(plot_title) +
    facet_grid(group~.) +
    scale_fill_viridis_c(option= "H") +
    theme_void() +
    theme(legend.key.width = unit(colorbar_width, "cm")) +
    easy_move_legend(to = 'bottom') + easy_adjust_legend("center") +
    easy_center_title() + easy_add_legend_title(lengend_title) +
    easy_text_size(which = c("plot.title", "legend.title", "strip.text"),
                   size = 15)
  
  if (hemi_position == "vertical") {
    brain_pos <- position_brain(hemi + side ~ . )
    ## Plot atlas:
    p_this <- ggplot() + geom_brain(
      atlas       = atlas_data,
      mapping     = aes(fill=value),
      position    = brain_pos,
      color       ='black',
      size        = line_size,
      show.legend = show_legend) +
      ggtitle(plot_title) +
      facet_grid(.~group) +
      scale_fill_viridis_c(option= "H") +
      theme_void() +
      theme(legend.key.width = unit(colorbar_width, "cm")) +
      easy_move_legend(to = 'bottom') + easy_adjust_legend("center") +
      easy_center_title() + easy_add_legend_title(lengend_title) +
      easy_text_size(which = c("plot.title", "legend.title", "strip.text"),
                     size = 15)
  }else if(hemi_position == "montage"){
    brain_pos <- position_brain(side ~ hemi)
    p_this <- ggplot() + geom_brain(
      atlas       = atlas_data,
      mapping     = aes(fill=value),
      position    = brain_pos,
      color       ='black',
      size        = line_size,
      show.legend = show_legend) +
      ggtitle(plot_title) +
      facet_grid(.~group) +
      scale_fill_viridis_c(option= "H") +
      theme_void() +
      theme(legend.key.width = unit(colorbar_width, "cm")) +
      easy_move_legend(to = 'bottom') + easy_adjust_legend("center") +
      easy_center_title() + easy_add_legend_title(lengend_title) +
      easy_text_size(which = c("plot.title", "legend.title", "strip.text"),
                     size = 15)
  }
  
  return(p_this)
}

