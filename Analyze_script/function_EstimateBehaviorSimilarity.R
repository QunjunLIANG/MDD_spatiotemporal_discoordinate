#' Function to estimate similarity matrix in behavior.
#'
#' Function is adjusted from https://naturalistic-data.org/content/Intersubject_RSA.html#hypothesis-tests
#'
#' @param behav_rank A vector contains the RANKs of subjects' behavioral outcome or measures.
#' @param type A string to indicate which distance should be calculate, alternative in "Eculidean" or "AK"


EstimateBehaviorSimilarity <- function(behav_rank, type){
  
  behav_sim_annak <- matrix(nrow = length(behav_rank), ncol = length(behav_rank))
  
  if (type == 'Euclidean') { # use Eclidean distance
    behav_sim_annak <- dist(behav_rank) %>% as.matrix()
    ## turn the distance matrix to similarity matrix
    behav_sim_annak <-  -1 * behav_sim_annak / sd(behav_sim_annak)
  }else{
    ## calculate the A-K distance
    ## AK distance did not need to reverse to similarity matrix because
    ## the low value gets the high rank in R. 
    for (i in 1:length(behav_rank)) {
      for (j in 1:length(behav_rank)) {
        if (i < j) {
          sim_ij <- mean(c(behav_rank[i], behav_rank[j]))/length(behav_rank)
          behav_sim_annak[i,j] <- sim_ij
          behav_sim_annak[j,i] <- sim_ij
        } else if (i == j) {
          behav_sim_annak[i,j] <- 1
        }
      }
    }
  }
  
  return(behav_sim_annak)
}