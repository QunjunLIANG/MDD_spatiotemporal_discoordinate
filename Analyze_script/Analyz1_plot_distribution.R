###################
#
# This pipeline is used to visualize the 
# distribution of TD gradient and functional gradient
#
# Qunjun Liang 2023/02/23

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

# indicate the subject name -------------------------------------------------
info_dat <- import('MDD_HC_match_data.xlsx')

dat_file_grad <- list.files('inputs/gradient_inputs', pattern = 'sub*', full.names = T)

dat_file_tde_topo <- list.files('inputs/TDE_Gradient', pattern = 'sub.*align-procrustes.csv', full.names = T)

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

net_ana$network <- factor(net_ana$network, levels = c("somMot","visual",'dorsalAttn',"salience","control","limbic","DMN"))

outDir <- 'output'
if (!dir.exists(outDir)) {
  dir.create(outDir)
}

# load subject data -------------------------------------------------------------
## load the TD topology data
if (!exists('dat_td_topo')) {
  dat_td_topo <- data.frame()
  for (i in dat_file_tde_topo) {
    sbjID_tmp <- str_extract(i, pattern = 'sub-[0-9]{3}')
    dat_tmp <- read_csv(i) %>% t() %>% as.data.frame()
    colnames(dat_tmp) <- paste0('parcel',1:400)
    dat_tmp['participant_id'] <- sbjID_tmp
    dat_tmp['gradient'] <- c("gradient1","gradient2",'gradient3')
    dat_td_topo <- rbind(dat_td_topo, dat_tmp)
  }
  
  dat_td_topo <- dat_td_topo %>% left_join(info_dat[,c(1,71)])
}

dat_td_topo_long <- dat_td_topo %>% filter(gradient == 'gradient1') %>%
  select(- gradient) %>%
  melt(id.vars = c('participant_id', "group"), 
       variable.name = "parcel", value.name = 'td_topo')

## load the functional gradient data 
if (!exists('dat_grad')) {
  dat_grad <- data.frame()
  for (i in dat_file_grad) {
    sbjID_tmp <- str_extract(i, pattern = 'sub-[0-9]{3}')
    dat_tmp <- read_csv(i) %>% t() %>% as.data.frame()
    colnames(dat_tmp) <- paste0('parcel',1:400)
    dat_tmp['participant_id'] <- sbjID_tmp
    dat_tmp['gradient'] <- c("gradient1","gradient2",'gradient3')
    dat_grad <- rbind(dat_grad, dat_tmp)
  }
  
  dat_grad <- dat_grad %>% left_join(info_dat[,c(1,71)])
}
dat_grad_long <- dat_grad %>% filter(gradient == 'gradient1') %>%
  select(- gradient) %>%
  melt(id.vars = c('participant_id', "group"), 
       variable.name = "parcel", value.name = 'gradient')

dat_all_long <- dat_grad_long %>% left_join(dat_td_topo_long)
dat_all_long <- dat_all_long %>% left_join(net_ana)

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
  xlab('TD-based gradient') + ylab('Density') + 
  scale_fill_manual(values = c(hc_color, mdd_color)) +
  theme_classic() + theme(legend.position = c(0.85, 0.85)) +
  easy_text_size(size = 15)
p_td_topo_distribution

# gradient1 distribution ---------------------------------------------------------

## group difference
p_grad_distribution <- ggplot(data = dat_tmp, aes(x = gradient, y =..density.., fill = group)) + 
  geom_histogram(data = dat_hc, 
                 bins = 50, 
                 aes(x = gradient, y =..density..), 
                 fill = hc_color,
                 color="white", alpha = 0.8)+
  geom_histogram(data = dat_mdd, 
                 bins = 50, aes(x = gradient, y = ..density..),
                 fill = mdd_color,
                 color="white", alpha = 0.8) +
  geom_density(aes(x = gradient, y =-..density.., fill = group),
               color = 'white', alpha = 0.6, size = 1) +
  xlab('Functional gradient score') + ylab('Density') + 
  scale_fill_manual(values = c(hc_color, mdd_color)) +
  theme_classic() + theme(legend.position = c(0.85, 0.85)) +
  easy_text_size(size = 15)
p_grad_distribution

# Network specific distributions ------------------------------------------

dat_tmp <- dat_all_long[,c(-6,-7)] %>% 
  melt(id.vars = c('participant_id', "group", 'parcel', 'network'), 
       variable.name = "measure", value.name = 'value') %>%
  group_by(group,parcel,network, measure) %>%
  summarize_each(funs = funs(mean))

facet_labels <- c(gradient = "Functional gradient", td_topo="TD-based gradient")
p_net_distribution <- ggplot(data = dat_tmp, aes(y=network, x=value, fill = group)) +
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

# correlation between TD and gradient -------------------------------------

dat_tmp <- dat_all_long %>% group_by(group,participant_id) %>%
  summarise(across(c("gradient1", "gradient2", "td_topo"), mean)) %>%
  ungroup()

grouped_ggscatterstats(dat_tmp,
               x = gradient1,
               y = td_topo,
               grouping.var = group,
               bf.message = F,
               xlab = "Functional gradient",
               ylab = "TD-based gradient",
               marginal = F,
               ggtheme = theme_classic(base_size = 16)
               )

p_corplot <- ggplot(dat_tmp, aes(x = gradient, y = td_topo)) +
  geom_point(size = 2.4, alpha = 0.5) +
  geom_smooth(method = 'lm', color ='red') +
  scale_color_manual(values = c(hc_color, mdd_color)) +
  facet_grid(~group, scale = 'free') +
  ylab("TD-based gradient") + xlab("Functional gradient") +
  theme_classic(base_size = 15, base_line_size = .7) +
  theme(strip.background = element_rect(colour = 'white'),
        strip.text.x = element_text(size = 18, face = 'bold')) +
  easy_remove_legend()
p_corplot

# combine the plots -------------------------------------------------------
p_distribution <- (p_td_topo_distribution + p_grad_distribution)/p_net_distribution/p_corplot
p_distribution + plot_annotation(tag_levels = c('a', '1')) & 
  theme(plot.tag = element_text(size = 22, face = 'bold'))
ggsave('output/Distribution_plots.png', width = 9, height = 10)

# export the collected data -------------------------------------------------

dat_grad_long <- dat_grad %>% 
  melt(id.vars = c('participant_id', "group", "gradient"), 
       variable.name = "parcel", value.name = 'value') %>%
  pivot_wider(names_from = gradient, values_from = value)

dat_all_long <- dat_grad_long %>% left_join(dat_td_topo_long)
dat_all_long <- dat_all_long %>% left_join(net_ana)

export(dat_all_long, file = file.path(outDir, 'TD_topo_Gradient_collect_long.csv'))
