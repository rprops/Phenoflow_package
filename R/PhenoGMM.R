#' Train GMM-fitted model to FCS data. 
#'
#' @param fcs_x flowSet object with input data on which the model should be built 
#' @param downsample Indicate to which sample size individual samples should be downsampled. 
#' By default no downsampling is performed
#' @param nG Number of mixtures to use. Defaults to 128.
#' @param param parameters to be used in the mixture modeling.
#' @importFrom BiocGenerics unique colnames
#' @importFrom mclust Mclust predict.Mclust
#' @importFrom tidyr "%>%"
#' @keywords fingerprint
#' @examples
#' data(flowData_transformed)
#' testGMM <- PhenoGMM(flowData_transformed, downsample = 1e3, nG = 128, param = c("FL1-H", "FL3-H"))
#' testPred <- PhenoMaskGMM(flowData_transformed, gmm = testGMM)
#' @export

PhenoGMM <- function(fcs_x, param, downsample = 0, nG = 128){
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
  
  # Standardize data (center & scale)
  fcs_m <- flowCore::fsApply(fcs_m, FUN = function(x) scale(x), use.exprs = TRUE)

  # Start performing MClust for GMM estimation
  gmm_clust <- Mclust(data = fcs_m, G = nG)
  
  # Normalize input data and assign cluster allocations
  fcs_x_t <- flowCore::fsApply(fcs_x, FUN = function(x) {
    tmp_nm <- base::sweep(x, 2, fcs_m_colM)/c(fcs_m_sd)
    data.frame(table(predict.Mclust(gmm_clust, tmp_nm)$classification))
    }, use.exprs = TRUE, simplify = FALSE)
  
  # Make mixture contigency table
  gmm_pred <- suppressWarnings(dplyr::bind_rows(fcs_x_t, .id = "id") %>%
    tidyr::spread(Var1, Freq))
  colnames(gmm_pred)[1] <- "sample_names"
  
  # Replace NAs with zeros
  gmm_pred[is.na(gmm_pred)] <- 0L
  
  # Use model to predict cluster abundances of all data
  fp_return <- list()
  fp_return[[1]] <- gmm_pred # model applied to input sample data
  fp_return[[2]] <- gmm_clust # model for future predictions
  
  # Add GMM parameters as attributes to data.table object
  attr(fp_return[[1]], "nG") <- nG
  attr(fp_return[[1]], "param") <- param
  attr(fp_return[[1]], "parameter_mean") <- fcs_m_colM
  attr(fp_return[[1]], "parameter_sd") <- fcs_m_sd
  
  # }) # enable for profviz
  
  # Return fingerprint object
  return(fp_return)
}
