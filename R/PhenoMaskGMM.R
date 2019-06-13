#' Apply PhenoGMM made model to new data.
#'
#' @param fcs_x flowSet object on which to apply the GMM mask
#' @param gmm model mask used to assign cluster allocations
#' @param fcs_scale Should data be scaled/normalized by row and column before running GMM? Defaults to FALSE. 
#' Be aware that if input modelled used rescaled data as input/option this should be set to TRUE for predicting on new data.
#' @importFrom BiocGenerics unique colnames
#' @keywords fingerprint
#' @importFrom mclust predict.Mclust
#' @importFrom magrittr %>%
#' @examples
#' data(flowData_transformed)
#' testGMM <- PhenoGMM(flowData_transformed, downsample = 1e3, nG = 128, param = c("FL1-H", "FL3-H"))
#' testPred <- PhenoMaskGMM(flowData_transformed, gmm = testGMM)
#' @export

PhenoMaskGMM <- function(fcs_x, gmm, fcs_scale = FALSE){
  # profvis({ # for profiling performance
  # Subset to parameters available in model
  fcs_x <- fcs_x[, attributes(gmm[[1]])$param]
  
  # Extract colmeans and colsd used in model creation
  # to normalize data
  fcs_x_colM <- attributes(gmm[[1]])$parameter_mean
  fcs_x_sd <- attributes(gmm[[1]])$parameter_sd
  
  if(fcs_scale){
    # Normalize input data and assign cluster allocations
    fcs_x_t <- flowCore::fsApply(fcs_x, FUN = function(x) {
      tmp_nm <- base::sweep(x, 2, fcs_x_colM)/c(fcs_x_sd)
      data.frame(table(mclust::predict.Mclust(gmm[[2]], tmp_nm)$classification))
    }, use.exprs = TRUE, simplify = FALSE)
  } else {
    # Do not normalize input data and assign cluster allocations
    fcs_x_t <- flowCore::fsApply(fcs_x, FUN = function(x) {
      data.frame(table(mclust::predict.Mclust(gmm[[2]], x)$classification))
    }, use.exprs = TRUE, simplify = FALSE)
  }
  
  # Make mixture contigency table
  gmm_pred <- suppressWarnings(dplyr::bind_rows(fcs_x_t, .id = "id") %>%
    tidyr::spread(Var1, Freq))
  colnames(gmm_pred)[1] <- "Sample_names"
  
  # Replace NAs with zeros
  gmm_pred[is.na(gmm_pred)] <- 0L
  
  # Use model to predict cluster abundances of all data
  fp_return <- list()
  fp_return[[1]] <- gmm_pred # model applied to input sample data
  fp_return[[2]] <- gmm[[2]] # model for future predictions

  
  # Add GMM parameters as attributes to data.table object
  attr(fp_return[[1]], "nG") <- attributes(gmm[[1]])$nG
  attr(fp_return[[1]], "param") <- attributes(gmm[[1]])$param
  attr(fp_return[[1]], "parameter_mean") <- fcs_x_colM
  attr(fp_return[[1]], "parameter_sd") <- fcs_x_sd
  
  # }) # enable for profviz
  
  return(fp_return)
}
