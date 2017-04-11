#' FCS_clean function for FCM data
#'
#' This function denoises a flowSet object using the flowClean package. Observations that
#' do not need meet the criteria will be removed. Only works for samples with more than
#' 30,000 observations.
#' @param x flowSet containing the samples to analyse. Cannot have been subsetted to 
#' a lower number of parameters.
#' @param cleanparam Indices of parameters to be used for removing errant observations. 
#' Defaults to 9 and 11 for FL1-H and FL3-H on the BD C6 Accuri FCM.
#' @keywords denoise, fcm, alpha
#' @examples
#' ## Full data processing example
#' data(flowData)
#' cleanparam <- c(9,11)
#' flowData_cleaned <- flowCore::fsApply(x = flowData, FUN = function(x) FCS_clean(x, cleanparam))
#' @export

FCS_clean <- function(x, cleanparam = c(9,11)){
  if(nrow(exprs(x)) > 30000){
    GoodCells <- flowClean::clean(x, vectMarkers = cleanparam, nCellCutoff = 500,
                                  binSize = 0.01, returnVector = TRUE,
                                  filePrefixWithDir = paste(x@description$FILENAME, "_cleaned.fcs", 
                                                            sep = ""), ext = "")
    flowCore::exprs(x) <- flowCore::exprs(x)[GoodCells < 10000, ]
    
  }
  return(x)
}
