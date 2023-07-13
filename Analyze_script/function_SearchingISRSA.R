
SearchingISRSA <- function(data, threshold, var.y="wave1_scale", randomSeed=130,
                           feature_thresh=0.1, dist_type='Euclidean',
                           cor_method='spearman', n_permutations=5000) {
  library(tidyverse)
  library(ggstatsplot)
  source('Analyze_script/function_EstimateBehaviorSimilarity.R')
  
  ## separate the dataset into sub-dataset
  set.seed(randomSeed)
  
  dat_train <- NA
  dat_test <- NA
  
  if (threshold <= 1) {
    rsa_ind_tmp <- caret::createDataPartition(data[var.y]%>% unlist(), p = threshold)
    dat_train <- data[rsa_ind_tmp$Resample1, ]
    dat_test <- data[-rsa_ind_tmp$Resample1, ]
  }else{
    rowN <- 1:nrow(data)
    Resample1 <- sample(rowN, size = threshold)
    dat_train <- data[Resample1, ]
    dat_test <- data[-Resample1, ]
  }
  
  
  ## data preprocessing - scale the features
  dat_train <- dat_train[,(2:ncol(dat_train))]
  
  ## feature selection with correlation
  # feature_select_tmp <- ggcorrmat(dat_train,
  #                                 sig.level = feature_thresh,
  #                                 type = 'nonparameteric',
  #                                 p.adjust.method = 'none', 
  #                                 output = "dataframe") %>%
  #   filter(parameter2 == var.y, p.value <= feature_thresh) %>% 
  #   .$parameter1
  feature_select_tmp <- ggcorrmat(dat_train,
                                  sig.level = feature_thresh,
                                  type = 'nonparameteric',
                                  p.adjust.method = 'none',
                                  output = "dataframe") %>%
    filter(parameter2 == var.y, (abs(estimate) >= feature_thresh)) %>%
    .$parameter1
  
  if (is_empty(feature_select_tmp) | length(feature_select_tmp) < 3) {
    result_dat <- data.frame(propotion=threshold, 
                             n_train = nrow(dat_train),
                             n_test = nrow(dat_test),
                             feature_thresh = feature_thresh,
                             cor_method = cor_method,
                             n_feature = 0,
                             dist_method = dist_type,
                             rsa=0, 
                             p_value=1)
    return(result_dat)
  }else{
    
    dat_test <- dat_test %>% select(feature_select_tmp, var.y)
    
    ## shape brain and behavior data
    dat_brain_act <- dat_test[,-ncol(dat_test)]
    dat_brain_corMat <- cor(dat_brain_act %>% t())
    dat_brain_dist <- dat_brain_corMat[lower.tri(dat_brain_corMat)]
    
    dat_scale <- dat_test[var.y] %>% rank(max(dat_test[var.y])-.)
    dat_scale_dist <- EstimateBehaviorSimilarity(dat_scale, type = dist_type) # for one of "AK" or "Euclidean"
    dat_scale_dist <- dat_scale_dist[lower.tri(dat_scale_dist)]
    
    ## calculate the RSA
    rsa_tmp <- cor.test(dat_brain_dist, dat_scale_dist, method = cor_method)
    
    ## permutation test for 5000 permutes
    
    # res_perm <- numeric(n_permutations)
    # for (i in 1:n_permutations) {
    #   var_perm <- sample(dat_brain_dist, replace = FALSE)
    #   res_perm[i] <- cor(dat_scale_dist, var_perm, method = cor_method)
    # }
    res_perm <- replicate(n_permutations, {
      var_perm <- sample(dat_brain_dist, replace = FALSE)
      cor(dat_scale_dist, var_perm, method = cor_method)
    })
    
    p_value <- sum(abs(res_perm) >= abs(rsa_tmp$estimate)) / n_permutations
    
    result_dat <- data.frame(propotion=threshold, 
                             n_train = nrow(dat_train),
                             n_test = nrow(dat_test),
                             feature_thresh = feature_thresh,
                             cor_method = cor_method,
                             n_feature = length(feature_select_tmp),
                             dist_method = dist_type,
                             rsa=rsa_tmp$estimate, 
                             p_value=p_value)
    
    
    return(result_dat)
  }
}

sort_square_mtx <- function(mtx, vct){
  # Sorts rows/columns of a matrix according to a separate vector.
  inds <- order(vct)
  mtx_sorted <- mtx[inds, inds]
  
  return(mtx_sorted)
}

ObtainImportFeature <- function(data, threshold, feature_thresh=0.1, var.y="wave1_scale"){
  ## separate the dataset into sub-dataset
  set.seed(130)
  rsa_ind_tmp <- caret::createDataPartition(data[var.y]%>% unlist(), p = threshold)
  dat_train <- data[rsa_ind_tmp$Resample1, ]
  
  ## data preprocessing - scale the features
  dat_train <- dat_train[,(2:ncol(dat_train))]
  
  ## feature selection with correlation
  # feature_select_tmp <- ggcorrmat(dat_train,
  #                                 sig.level = feature_thresh,
  #                                 type = 'nonparameteric',
  #                                 p.adjust.method = 'none', 
  #                                 output = "dataframe") %>%
  #   filter(parameter2 == var.y, p.value <= feature_thresh) %>% 
  #   .$parameter1
  
  feature_select_tmp <- ggcorrmat(dat_train,
                                  sig.level = feature_thresh,
                                  type = 'nonparameteric',
                                  p.adjust.method = 'none',
                                  output = "dataframe") %>%
    filter(parameter2 == var.y, (estimate >= feature_thresh | estimate <= - feature_thresh)) %>%
    .$parameter1
  
  return(feature_select_tmp)
}
