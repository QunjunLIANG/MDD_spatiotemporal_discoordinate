#############################
#
# Functions used to check the ROI matrix
# and extract the specific parcellation

RecordRoIdimensionREST <- function(sbj_name, matDir,
                                   col_range=NULL){
  library(bruceR)
  library(R.matlab)
  library(purrr)
  input_tmp <- file.path(matDir,Glue('ROISignals_{sbj_name}.mat')) %>%
    readMat()
  input_tmp <- input_tmp$ROISignals
  if (is_empty(col_range)) {
    output_tmp <- data.frame(ID = sbj_name, 
                             timepoints = nrow(input_tmp),
                             Ncol = ncol(input_tmp))
  }else{
    if (max(col_range)>ncol(input_tmp)) {
      output_tmp <- data.frame(ID = sbj_name, 
                               timepoints = nrow(input_tmp),
                               Ncol = ncol(input_tmp),
                               parcel_complete = NA)
    }else{
      parcel_mat <- input_tmp[,col_range]
      zero_col <- apply(parcel_mat, 2, function(x) all(x==0))
      zero_ind <- mean(zero_col)
      output_tmp <- data.frame(ID = sbj_name, 
                               timepoints = nrow(input_tmp),
                               Ncol = ncol(input_tmp),
                               parcel_complete = zero_ind)
    }
  }
  
  return(output_tmp)
}

SliceParcelREST <- function(sbj_name, matDir, outDir,
                                 col_range,
                                 outname){
  library(bruceR)
  library(R.matlab)
  # create the directory if not exists
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = T)
  }
  
  input_tmp <- file.path(matDir,Glue('ROISignals_{sbj_name}.mat')) %>%
    readMat()
  input_tmp <- input_tmp$ROISignals
  output_tmp <- input_tmp[,col_range]
  writeMat(data = output_tmp, 
           file.path(outDir, Glue('{sbj_name}_{outname}.mat')))
}

SliceParcelREST_csv <- function(sbj_name, matDir, outDir,
                            col_range,
                            outname){
  library(readr)
  library(bruceR)
  library(R.matlab)
  # create the directory if not exists
  if (!dir.exists(outDir)) {
    dir.create(outDir, recursive = T)
  }
  
  input_tmp <- file.path(matDir,Glue('ROISignals_{sbj_name}.mat')) %>%
    readMat()
  input_tmp <- input_tmp$ROISignals
  output_tmp <- input_tmp[,col_range] %>%as.data.frame()
  write_csv(output_tmp, col_names = F,
           file.path(outDir, Glue('{sbj_name}_{outname}.csv')))
}
