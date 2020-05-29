#' Train GMM-fitted model to FCS data. 
#'
#' @param fcs_x flowSet object with input data on which the model should be built 
#' @param downsample Indicate to which sample size individual samples should be downsampled. 
#' By default no downsampling is performed
#' @param nG Number of mixtures to use. Defaults to 128.
#' @param auto_nG TRUE/FALSE. Option to choose best number of mixtures from 1:nG based on BIC. 
#' Defaults to FALSE which forces nG clusters.
#' @param nG_interval if auto_nG = TRUE, specify the intervals from nG_interval:nG 
#' to calculate BIC for. Defaults to 4.
#' @param param parameters to be used in the mixture modeling.
#' @param fcs_scale Should data be scaled/normalized by row and column before running GMM? Defaults to FALSE.
#' @param diagnostic_plot Specify whether a diagnostic plot should be made showing the cluster allocation
#' of each cell in the specified parameter space.
#' @importFrom BiocGenerics unique colnames
#' @importFrom mclust Mclust predict.Mclust
#' @importFrom magrittr %>%
#' @importFrom flowCore exprs
#' @keywords fingerprint
#' @examples
#' data(flowData_transformed)
#' testGMM <- PhenoGMM(flowData_transformed, downsample = 1e3, 
#' nG = 30,
#' auto_nG = TRUE,
#' nG_interval = 10,
#' param = c("FL1-H", "FL3-H"))
#' testPred <- PhenoMaskGMM(flowData_transformed, gmm = testGMM)
#' @export

PhenoGMM <-
  function(fcs_x,
    param,
    downsample = 0,
    nG = 128,
    auto_nG = FALSE,
    nG_interval = 4,
    fcs_scale = FALSE,
    diagnostic_plot = FALSE) {
    
  # profvis({ # for profiling performance
  
  # Select parameters of interest
  fcs_x <- fcs_x[, param]
  
  # Downsample if necessary
  if (downsample != 0)
    fcs_x_sb <-
      Phenoflow::FCS_resample(fcs_x, sample = downsample, replace = TRUE)
  else
    fcs_x_sb <- fcs_x
  
  # Merge all samples 
  fcs_m <- Phenoflow::FCS_pool(fcs_x_sb, stub = "*")
  fcs_m <- suppressWarnings(fcs_m[, param])
  
  # Register colmeans and colsd for applying the standardization on
  # new data
  fcs_m_colM <- flowCore::fsApply(fcs_m, FUN = function(x) colMeans(x), 
    use.exprs = TRUE)
  fcs_m_sd <- flowCore::fsApply(fcs_m, FUN = function(x) apply(x, 2, sd), 
    use.exprs = TRUE)
  
  if(fcs_scale){
    # Standardize data (center & scale)
    fcs_m <- flowCore::fsApply(fcs_m, FUN = function(x) scale(x), use.exprs = TRUE)
    # Start performing MClust for GMM estimation
    if (auto_nG) {
      BIC <- mclustBIC(fcs_m, G = seq(from = nG_interval, to = nG, by = nG_interval))
      gmm_clust <- Mclust(data = fcs_m, x = BIC)
    } else
      gmm_clust <- Mclust(data = fcs_m, G = nG)
    # Normalize input data and assign cluster allocations
    fcs_x_t <- flowCore::fsApply(fcs_x, FUN = function(x) {
      tmp_nm <- base::sweep(x, 2, fcs_m_colM)/c(fcs_m_sd)
      data.frame(table(predict.Mclust(gmm_clust, tmp_nm)$classification))
    }, use.exprs = TRUE, simplify = FALSE)
  } else {
    if (auto_nG) {
      # Start performing MClust for GMM estimation
      BIC <- mclustBIC(exprs(fcs_m[[1]]), G = seq(from = nG_interval, to = nG, by = nG_interval))
      gmm_clust <- Mclust(data = exprs(fcs_m[[1]]), x = BIC)
    } else
      gmm_clust <- Mclust(data = exprs(fcs_m[[1]]), G = nG)
    # Normalize input data and assign cluster allocations
    fcs_x_t <- flowCore::fsApply(fcs_x, FUN = function(x) {
      data.frame(table(predict.Mclust(gmm_clust, x)$classification))
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
  fp_return[[2]] <- gmm_clust # model for future predictions
  
  # Add GMM parameters as attributes to data.table object
  attr(fp_return[[1]], "nG") <- nG
  attr(fp_return[[1]], "auto_nG") <- auto_nG
  attr(fp_return[[1]], "nG_interval") <- nG_interval
  attr(fp_return[[1]], "param") <- param
  attr(fp_return[[1]], "parameter_mean") <- fcs_m_colM
  attr(fp_return[[1]], "parameter_sd") <- fcs_m_sd
  
  # }) # enable for profviz
  
  if(diagnostic_plot){
    gmm_pred
  }
  # Return fingerprint object
  return(fp_return)
}
