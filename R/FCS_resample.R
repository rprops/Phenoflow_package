#' Resampling function for flowSet objects
#'
#' This function resamples each sample to an equal number of cells. Options are to resample 
#' with or withour replacement, and to resample on a per sample basis or across the entire dataset
#' (normalization to equal sample size).
#' @param x flowSet object.
#' @param sample Desired sample size. Defaults to minimum sample size.
#' @param replace Do you want to resample with or without replacement? Defaults to FALSE, which is without replacement.
#' @param rarefy Do you want each sample resampled without adjusting sample size? Default to FALSE.
#' @param progress Should progress be reported? Defaults to yes.
# #' @param keep Should samples with lower sample sizes than the specified `sample` be kept? Defaults to `FALSE`.
#' @keywords resampling, fcm
#' @examples 
#' # Load raw data (imported using flowCore)
#' data(flowData)
#' flowData <- FCS_resample(flowData, replace=TRUE)
#' @export

FCS_resample <- function(x, sample = 0, replace = FALSE, rarefy = FALSE, 
                         progress = TRUE) {
  min <- min(flowCore::fsApply(x = x, FUN = function(x) nrow(x), 
                              use.exprs = TRUE))
  max <- max(flowCore::fsApply(x = x, FUN = function(x) nrow(x), 
                               use.exprs = TRUE))
  if (sample == 0) 
    sample <- min
  
  ## Remove all .fcs files with less observations than the specified
  ## sample
    x <- x[flowCore::fsApply(x = x, FUN = function(x) nrow(x), use.exprs = TRUE) >= 
                sample]
  
  if (rarefy == FALSE) {
    for (i in 1:length(x)) {
      flowCore::exprs(x[[i]]) <- flowCore::exprs(x[[i]])[base::sample(1:nrow(flowCore::exprs(x[[i]])), 
                                                                      size = sample, replace = replace), ]
    }
    if (progress == TRUE) {
      cat(paste0("Your samples range between ", min, 
                 " and ", max, " cells\n"))
      cat(paste0("Your samples were randomly subsampled to ", sample, 
                 " cells\n"))
    }
  } else {
    for (i in 1:length(x)) {
      flowCore::exprs(x[[i]]) <- flowCore::exprs(x[[i]])[base::sample(1:nrow(flowCore::exprs(x[[i]])), 
                                                                      size = nrow(flowCore::exprs(x[[i]])), replace = replace), 
                                                         ]
    }
    if (progress == TRUE) 
      cat(paste0("Your samples were randomly subsampled to their respective sample size\n"))
  }
  return(x)
}