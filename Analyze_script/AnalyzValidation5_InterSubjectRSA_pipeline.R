#################################
#
# IS-RSA pipeline 
#
# This script is used to explore the relationship
# between the gradient topology of the brain and the 
# HAMD score.
#
# The pipeline is composed of 5 parts:
#   1. plot the correlation matrix between gradient topology and HAMD score.
#   2. plot the example neruo activity similarity matrix.
#   3. plot the behavior similarity matrix in NN model and AK model, respectively.
#   4. run IS-RSA corresponding to NN and AK model.
#   5. Permutation Test.

library(tidyverse)
library(bruceR)
library(ggstatsplot)
library(ggeasy)
library(ggridges)
library(plotly)
library(ggeasy)
library(Cairo)
library(gplots)
library(webshot) # use to save wordcloud

dat_beha_file  <- 'inputs/MDD_validation_dataset.xlsx' # data file contains the HAMD score for each patients
dat_STdist     <- 'output/validation_TDtopoDistance_machineLearning.csv' # data file contains the TD gradient distance features

cor_method     <-  "spearman" # correlation method to use, with spearman correlation.
feature_thresh <-  0.1 # threshold for feature selection using correlation analysis.
n_permutations <- 5000

source('Analyze_script/function_EstimateBehaviorSimilarity.R')
source('Analyze_script/function_SearchingISRSA.R')

## plot the correlation matrix in STP distance #################################

## load the features
dat_ML <- import(dat_STdist) 

## load the demographic/clinical information
dat_info <- import(dat_beha_file) %>%
  select(participant_id, wave1_scale)

## combine the predict and
dat_ML <- dat_ML %>% left_join(dat_info)

## data preprocessing - scale the features
dat_ML <- dat_ML[,(2:ncol(dat_ML))]

feature_select_tmp <- ggcorrmat(dat_ML,
                                sig.level = feature_thresh,
                                type = 'nonparameteric',
                                p.adjust.method = 'none', 
                                output = "dataframe") %>%
  filter_if(parameter2 == 'wave1_scale',(abs(estimate) >= feature_thresh)) %>% 
  .$parameter1

dat_predict <- dat_ML %>% select(feature_select_tmp, wave1_scale)

## rank the HAMD scores
dat_beha <- dat_predict$wave1_scale %>% rank()

# plot the brain activity similarity matrix ######################################
## brain activity similarity
dat_brain <- dat_predict[,-ncol(dat_predict)]
dat_brain_cor <- cor(dat_brain %>% t())

## sort the matrix for plotting
dat_brain_cor_sort <- sort_square_mtx(dat_brain_cor, dat_beha)
Cairo(20,20, units="cm", file=bruceR::Glue("output/ISRSA_brain_similarity_matrtix.png"), 
      type="png", bg="white", dpi=300)
fields::image.plot(dat_brain_cor_sort,main = "Neuro activity similarity",
                   legend.width=1, 
                   horizontal=F, legend.shrink=1, 
                   col = colorpanel(100,"steelblue","white","darkred"),
                   legend.lab=c("Similarity"),
                   legend.cex = 1.5,
                   axes = F,
                   axis.args=list(at=c(0.9, -0.9), labels=c("High","Low"), tick=F, cex.axis=1.3),
                   legend.line = -1)
dev.off()

# plot the behavior similarity matrix ##########################################

## plot the Euclidean model 
dat_beha_dist <- EstimateBehaviorSimilarity(dat_beha, type = 'Euclidean') 
dat_beha_dist_sort <- sort_square_mtx(dat_beha_dist, dat_beha)
Cairo(20,20, units="cm", file=bruceR::Glue("output/ISRSA_behavior_similarity_matrtix_Euclidean.png"), 
      type="png", bg="white", dpi=300)
fields::image.plot(dat_beha_dist_sort,main = "Nearest-neighbor model",
                   legend.width=1, 
                   horizontal=F, legend.shrink=1, 
                   col = colorpanel(100,"steelblue","white","darkred"),
                   legend.lab=c("Similarity"),
                   legend.cex = 1.5,
                   axes = F,
                   axis.args=list(at=c(-0.1, -4), labels=c("High","Low"), tick=F, cex.axis=1.3),
                   legend.line = -1)
dev.off()

### calculate the ISRSA
dat_beha_rsa <- dat_beha_dist[upper.tri(dat_beha_dist)]
dat_brain_rsa <- dat_brain_cor[upper.tri(dat_brain_cor)]
rsa_Eud <- cor.test(dat_beha_dist, dat_brain_cor, method = cor_method)

## plot the A-K model
dat_beha_dist <- EstimateBehaviorSimilarity(dat_beha, type = 'AK') 
dat_beha_dist_sort <- sort_square_mtx(dat_beha_dist, dat_beha)
Cairo(20,20, units="cm", file=bruceR::Glue("output/ISRSA_behavior_similarity_matrtix_AK.png"), 
      type="png", bg="white", dpi=300)
fields::image.plot(dat_beha_dist_sort,main = "Anna-Karenina model",
                   legend.width=1, 
                   horizontal=F, legend.shrink=1, 
                   col = colorpanel(100,"steelblue","white","darkred"),
                   legend.lab=c("Similarity"),
                   legend.cex = 1.5,
                   axes = F,
                   axis.args=list(at=c(-0.1, -4.6), labels=c("High","Low"), tick=FALSE, cex.axis=1.3),
                   legend.line = -1)
dev.off()

### calculate the ISRSA
dat_beha_rsa <- dat_beha_dist[upper.tri(dat_beha_dist)]
dat_brain_rsa <- dat_brain_cor[upper.tri(dat_brain_cor)]
rsa_AK <- cor.test(dat_beha_dist, dat_brain_cor, method = cor_method)


# permutation test for Euclidean distance
res_perm <- replicate(n_permutations, {
  var_perm <- sample(dat_beha_rsa, replace = FALSE)
  cor(dat_brain_rsa, var_perm, method = cor_method)
})

p_value <- sum(res_perm >= rsa_Eud$estimate) / n_permutations

data.frame(cor = res_perm) %>%
  ggplot() +
  geom_histogram(aes(x = cor), fill = 'steelblue',
                 color = 'white', bins = 15) +
  geom_vline(xintercept = rsa_Eud$estimate, color = 'darkred') +
  ylab('Frequency') + xlab('Similarity') + ggtitle('Similarity distribution in 5000 permutations') +
  theme_classic(base_size = 15) +
  easy_center_title()
ggsave(filename = 'output/ISRSA_Euclidean_permutationTest.png')