library(rio)
library(tidyverse)
library(ggridges)
library(plotly)
library(RColorBrewer)
library(ggsci)
library(ggeasy)
library(patchwork)
library(bruceR)
library(ggstatsplot)

dat_all_long <- rio::import('output/TD_topo_Gradient_collect_long.csv')
dat_all_long$network <- factor(dat_all_long$network, levels = c("somMot","visual",'dorsalAttn',"salience","control","limbic","DMN"))

# color identification --------------------------------------------------------
#
# obtain the color representations of two group of subjects
isfahan <- MetBrewer::met.brewer("Isfahan1")
length(isfahan)
isfahan[1]

hc_color <- colorspace::lighten(isfahan[2], 0.35)
mdd_color <-  colorspace::lighten(isfahan[6], 0.35)

## network color 
net_color <- c("#a153a2","#6fabd2","#2c8b4b","#b77fb4","#e7edca","#edaf5e","#e27283")

# TD gradient distribution ------------------------------------------------------
## group difference
dat_hc <- dat_all_long %>% filter(group == 'HC') %>% group_by(parcel) %>%
  summarize_each(funs = funs(mean))
dat_mdd <- dat_all_long %>% filter(group == 'MDD') %>% group_by(parcel) %>%
  summarize_each(funs = funs(mean))
dat_tmp <- dat_all_long %>% group_by(group,parcel) %>% summarize_each(funs = funs(mean))

p_td_topo_distribution <- ggplot(data = dat_tmp, aes(x = td_topo, y =..density.., fill = group)) + 
  geom_histogram(data = dat_hc, 
                 bins = 50, 
                 aes(x = td_topo, y =..density..), 
                 fill = hc_color,
                 color="white", alpha = 0.8)+
  geom_histogram(data = dat_mdd, 
                 bins = 50, aes(x = td_topo, y = ..density..),
                 fill = mdd_color,
                 color="white", alpha = 0.8) +
  geom_density(aes(x = td_topo, y =-..density.., fill = group),
               color = 'white', alpha = 0.6, size = 1) +
  xlab('gradient score') + ylab('Density') + 
  ggtitle('Temporal gradient') +
  scale_fill_manual(values = c(hc_color, mdd_color)) +
  theme_classic() + theme(legend.position = c(0.85, 0.85)) +
  easy_text_size(size = 15) + easy_center_title()
p_td_topo_distribution

# gradient1 distribution ---------------------------------------------------------

## group difference
p_grad_distribution <- ggplot(data = dat_tmp, aes(x = gradient1, y =..density.., fill = group)) + 
  geom_histogram(data = dat_hc, 
                 bins = 50, 
                 aes(x = gradient1, y =..density..), 
                 fill = hc_color,
                 color="white", alpha = 0.8)+
  geom_histogram(data = dat_mdd, 
                 bins = 50, aes(x = gradient1, y = ..density..),
                 fill = mdd_color,
                 color="white", alpha = 0.8) +
  geom_density(aes(x = gradient1, y =-..density.., fill = group),
               color = 'white', alpha = 0.6, size = 1) +
  xlab('gradient score') + ylab('Density') + 
  ggtitle('Functional gradient 1') +
  scale_fill_manual(values = c(hc_color, mdd_color)) +
  theme_classic() + theme(legend.position = c(0.85, 0.85)) +
  easy_text_size(size = 15) + easy_center_title()
p_grad_distribution

# gradient2 distribution ---------------------------------------------------------

## group difference
p_grad2_distribution <- ggplot(data = dat_tmp, aes(x = gradient2, y =..density.., fill = group)) + 
  geom_histogram(data = dat_hc, 
                 bins = 50, 
                 aes(x = gradient2, y =..density..), 
                 fill = hc_color,
                 color="white", alpha = 0.8)+
  geom_histogram(data = dat_mdd, 
                 bins = 50, aes(x = gradient2, y = ..density..),
                 fill = mdd_color,
                 color="white", alpha = 0.8) +
  geom_density(aes(x = gradient2, y =-..density.., fill = group),
               color = 'white', alpha = 0.6, size = 1) +
  xlab('gradient score') + ylab('Density') + 
  ggtitle('Functional gradient 2') +
  scale_fill_manual(values = c(hc_color, mdd_color)) +
  theme_classic() + theme(legend.position = c(0.85, 0.85)) +
  easy_text_size(size = 15) + easy_center_title()
p_grad2_distribution

# Network specific distributions ------------------------------------------

dat_tmp <- dat_all_long[,c(-8,-9)] %>% 
  melt(id.vars = c('participant_id', "group", 'parcel', 'network'), 
       variable.name = "measure", value.name = 'value') %>%
  group_by(group,parcel,network, measure) %>%
  summarize_each(funs = funs(mean))

facet_labels <- c(gradient1 = "Functional gradient1",
                  gradient2 = "Functional gradient2",
                  td_topo="Temporal gradient")
p_net_distribution <- dat_tmp %>% filter(measure!='gradient3') %>%
  ggplot(data = ., aes(y=network, x=value, fill = group)) +
  geom_density_ridges(scale = .8, alpha = 0.8,
                      size = 0.3, 
                      color = 'gray') +
  scale_fill_manual(values = c(hc_color, mdd_color)) +
  facet_grid(~measure, labeller = labeller(measure = facet_labels),
             scales = "free",
             space = "fixed") +
  xlab('Gradient score') + ylab('Density') + 
  theme_classic(base_size = 15, base_line_size = .5) + 
  theme(strip.background = element_rect(colour = 'white'),
        strip.text.x = element_text(size = 18),
        axis.text.y.left = element_text(size = 10)
  )

p_net_distribution

# brain surface plot ------------------------------------------------------
library(ggseg)
library(ggsegSchaefer)
## load the funtion from outer R scritp
source('Analyze_script/function_PlotOnBrainSurface_group.R')
source('Analyze_script/function_PermTestPlotOnBrain.R')

net_id <- read_csv('Yeo_net_identity_7.csv', col_names = 'network_ind') 
net_parcel <- read_csv('Yeo_7net_roiAnnotation.csv', col_names = 'region')
net_ana <- cbind(net_id, net_parcel)
net_ana$network_ind <- factor(net_ana$network_ind)
net_ana['parcel'] <- paste0('parcel', 1:400)
net_ana <- net_ana %>% mutate(network = ifelse(str_detect(region, pattern = 'Vis'),'visual',
                                               ifelse(str_detect(region, pattern = 'SomMot'),'somMot',
                                                      ifelse(str_detect(region, pattern = 'DorsAttn'),'dorsalAttn',
                                                             ifelse(str_detect(region, pattern = 'SalVentAttn'),'salience',
                                                                    ifelse(str_detect(region, pattern = 'Limbic'),'limbic',
                                                                           ifelse(str_detect(region, pattern = 'Cont'),'control',
                                                                                  'DMN')))))))




## Load in atlas data provided by ggseg package
atlas      = as_tibble(schaefer7_400)

## Select atlas region names and hemisphere so that we can add the values
## we want to plot:
region     = atlas$region
hemi       = atlas$hemi
data       = distinct(na.omit(data.frame(region,hemi))) #remove NA and duplicate regions

#####################################
#
# plot the Yeo 7 network template
#
######################################

data_Yeo7 <- data %>% left_join(net_ana)
atlas_temp = left_join(atlas,data_Yeo7)

atlas_data = left_join(atlas,data_Yeo7)
atlas_data$network[which(is.na(atlas_data$network))] <- "Medial Wall"
atlas_data$network <- factor(atlas_data$network, levels = c("visual","somMot",'dorsalAttn',"salience","limbic","control","DMN",'Medial Wall'))

p_yeo7 <- ggplot() + geom_brain(
  atlas       = atlas_data,
  mapping     = aes(fill=network),
  position    = position_brain(hemi~ side),
  color       ='black',
  size        = 0,
  show.legend = T) +
  ggtitle('Yeo 7 network template') +
  scale_fill_manual(values = c("#a153a2","#6fabd2","#2c8b4b","#b77fb4","#e7edca","#edaf5e","#e27283",'darkgrey')) +
  easy_all_text_size(13) +
  theme_void() +# easy_legend_at('bottom') +
  easy_center_title() + easy_add_legend_title('Network')
p_yeo7

#####################################
#
# plot TD topology in brain
#
#####################################

dat_file_grad <- list.files('inputs/gradient_inputs', pattern = 'sub*', full.names = T)

dat_file_tde_topo <- list.files('inputs/TDE_Gradient', pattern = 'sub.*align-procrustes.csv', full.names = T)


## load the data 
if (!exists('dat_tdtopo')) {
  dat_td_topo <- data.frame()
  for (i in dat_file_tde_topo) {
    sbjID_tmp <- str_extract(i, pattern = 'sub-[0-9]{3}')
    dat_tmp <- read_csv(i,show_col_types = FALSE) %>% t() %>%
      as.data.frame()
    colnames(dat_tmp) <- paste0('parcel',1:400)
    dat_tmp['participant_id'] <- sbjID_tmp
    dat_tmp['gradient'] <- c("gradient1","gradient2",'gradient3')
    dat_td_topo <- rbind(dat_td_topo, dat_tmp)
  }
  
  dat_td_topo <- dat_td_topo %>% left_join(info_dat[,c(1,71)])
}

p_TD_topo <- PlotFeatureOnSurface_group(atlas_temp, dat_td_topo, "TDtopology", "montage", 
                                        "Temporal gradient in groups", "Gradient score")
p_TD_topo

#####################################
#
# plot gradient 1&2 in brain
#
#####################################

## load the data 
if (!exists('dat_grad')) {
  dat_grad <- data.frame()
  for (i in dat_file_grad) {
    sbjID_tmp <- str_extract(i, pattern = 'sub-[0-9]{3}')
    dat_tmp <- read_csv(i,show_col_types = FALSE) %>% t() %>%
      as.data.frame()
    colnames(dat_tmp) <- paste0('parcel',1:400)
    dat_tmp['participant_id'] <- sbjID_tmp
    dat_tmp['gradient'] <- c("gradient1","gradient2",'gradient3')
    dat_grad <- rbind(dat_grad, dat_tmp)
  }
  
  dat_grad <- dat_grad %>% left_join(info_dat[,c(1,71)])
}

p_grad_1 <- PlotFeatureOnSurface_group(atlas_temp, dat_grad, "gradient", "montage", 
                                     "Fcuntional gradient 1 in groups", "Gradient score",
                                     gradient = 'gradient1')
p_grad_1

p_grad_2 <- PlotFeatureOnSurface_group(atlas_temp, dat_grad, "gradient", "montage", 
                                       "Fcuntional gradient 2 in groups", "Gradient score",
                                       gradient = 'grdient2')
p_grad_2


ggsave(plot = p_td_topo_distribution, 
       filename = 'output/plot_td_gradient.png',
       width = 4.54, height = 3.9)
ggsave(plot = p_grad2_distribution, 
       filename = 'output/plot_grad2_gradient.png',
       width = 4.54, height = 3.9)
ggsave(plot = p_grad_distribution, 
       filename = 'output/plot_grad1_gradient.png',
       width = 4.54, height = 3.9)
ggsave(plot = p_net_distribution, 
       filename = 'output/plot_net_distribution.png',
       width = 10, height = 4)
ggsave(plot = p_yeo7, 
       filename = 'output/plot_brain_yeo.png',
       width = 5, height = 4)
ggsave(plot = p_grad_1, 
       filename = 'output/plot_brain_gradient1.png',
       width = 5, height = 4)
ggsave(plot = p_TD_topo, 
       filename = 'output/plot_brain_TDgradient.png',
       width = 5, height = 4)
ggsave(plot = p_grad_2, 
       filename = 'output/plot_brain_gradient2.png',
       width = 5, height = 4)

#####################################
#
# plot TD & functional gradient correlation
#
#####################################

# calculate the mean gradients in HC and MDD, respectively
dat_hc_mean <- dat_all_long %>% filter(group=='HC') %>% group_by(parcel) %>% 
  summarise(mean_FG1 = mean(gradient1), mean_FG2 = mean(gradient2),
            mean_FG3 = mean(gradient3), mean_TD = mean(td_topo)) %>%
  mutate(group = 'HC')
dat_mdd_mean <- dat_all_long %>% filter(group=='MDD') %>% group_by(parcel) %>% 
  summarise(mean_FG1 = mean(gradient1), mean_FG2 = mean(gradient2),
            mean_FG3 = mean(gradient3), mean_TD = mean(td_topo)) %>%
  mutate(group = 'MDD')

dat_all_mean <- rbind(dat_hc_mean, dat_mdd_mean)

ggplot(dat_all_mean, aes(x = mean_FG1, y = mean_TD, color = group)) +
  geom_point(size = 2, shape = 4, alpha = .7) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = c(hc_color, mdd_color)) +
  xlab('Functional gradient score') + ylab('Temporal gradinet score') +
  theme_classic() +
  easy_text_size(15)
ggsave('output/Main_cor_TD_FG1.png', width = 7, height = 5)

ggplot(dat_all_mean, aes(x = mean_FG2, y = mean_TD, color = group)) +
  geom_point(size = 2, shape = 4, alpha = .7) +
  geom_smooth(method = 'lm') +
  scale_color_manual(values = c(hc_color, mdd_color)) +
  xlab('Functional gradient score') + ylab('Temporal gradinet score') +
  theme_classic() +
  easy_text_size(15)
ggsave('output/Main_cor_TD_FG2.png', width = 7, height = 5)


# calculate the correlation between gradients in HC and MDD, respectively
dat_hc_mean <- dat_all_long %>% filter(group=='HC') %>% group_by(participant_id) %>% 
  summarise(cor_FG1 = cor(gradient1, td_topo), 
            cor_FG2 = cor(gradient2, td_topo)) %>%
  mutate(group = 'HC')
dat_mdd_mean <- dat_all_long %>% filter(group=='MDD') %>% group_by(participant_id) %>% 
  summarise(cor_FG1 = cor(gradient1, td_topo), 
            cor_FG2 = cor(gradient2, td_topo)) %>%
  mutate(group = 'MDD')
dat_all_mean <- rbind(dat_hc_mean, dat_mdd_mean)

t.test(cor_FG1 ~ group, data = dat_all_mean)
t.test(cor_FG2 ~ group, data = dat_all_mean)
ggstatsplot::ggbetweenstats(
  data = dat_all_mean,
  x = group,
  y = cor_FG1
)
