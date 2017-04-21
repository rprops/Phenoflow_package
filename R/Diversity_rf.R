#' Diversity_rf function for FCM data
#'
#' This function calculates Hill diversity metrics from FCM data. This function differs from
#' the Diversity() function in that it resamples (with replacement) all individual samples and 
#' averages out the diversity over all subsamples. This function is recommended in case there there are
#' differences in sample size (nr. of cells). Analysis time is approximately 10s/resample run (R) on default
#' settings for the flowData example.
#' @param x flowSet containing the samples to analyse.
#' @param d Rounding factor for density values. Defaults to 4.
#' @param R Number of resampling runs to conduct on individual samples. Defaults to 100
#' @param R.b Number of bootstraps to conduct on the fingerprint (requires less). Defaults to 100
#' @param bw Bandwidth used in the kernel density estimation. Defaults to 0.01 which is ideal for normalized
#' FCM data (i.e., all parameter values are within [0,1]).
#' @param nbin Resolution of the binning grid. Defaults to 128 bins which corresponds to a 128x128 binning grid.
#' @param param Parameter vector indicating on which parameters the diversity should be estimated. An example input would be:
#' c("FL1-H", "FL3-H", "SSC-H", "FSC-H").
#' @param parallel Should the calculation be parallelized? Defaults to FALSE
#' @param ncores How many cores should be used in case of parallel computation?
#' Defaults to 2.
#' @param cleanFCS Indicate whether outlier removal should be conducted prior to diversity assessment (flowClean package). 
#' Defaults to FALSE. Requires not subsetted flowSet object (i.e. no parameters should be removed or otherwise flowClean functions fail).
#' Only samples with > 30,000 cells will be cleaned.
#' @param cleanparam Indices of parameters to be used for removing errant observations. 
#' Defaults to 9 and 11 for FL1-H and FL3-H on the BD C6 Accuri FCM.
#' @keywords diversity, fcm, alpha
#' @examples
#' ## Full data processing example
#' 
#' # Load raw data (imported using flowCore)
#' data(flowData)
#' 
#' # Take subsample
#' flowData <- flowData[1:5]
#' 
#' # Asinh transform and select parameters of interest (cells were stained with Sybr Green I).
#' flowData_transformed <- flowCore::transform(flowData,`FL1-H`=asinh(`FL1-H`),
#'        `SSC-H`=asinh(`SSC-H`),
#'        `FL3-H`=asinh(`FL3-H`),
#'        `FSC-H`=asinh(`FSC-H`))
#' param=c('FL1-H', 'FL3-H','SSC-H','FSC-H', 'Time')
#'
#' # Create a PolygonGate for denoising the dataset
#' # Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
#' sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3),ncol=2, nrow=4)
#' colnames(sqrcut1) <- c('FL1-H','FL3-H')
#' polyGate1 <- flowCore::polygonGate(.gate=sqrcut1, filterId = 'Total Cells')
#'
#' # Gating quality check
#' flowViz::xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=polyGate1,
#'          scales=list(y=list(limits=c(0,14)),
#'          x=list(limits=c(6,16))),
#'          axis = lattice::axis.default, nbin=125,
#'          par.strip.text=list(col='white', font=2, cex=2), smooth=FALSE)
#'
#' # Isolate only the cellular information based on the polyGate1
#' flowData_transformed <- flowCore::Subset(flowData_transformed, polyGate1)
#'
#' # Save one object with all parameters and one with the subset of parameters
#' flowData_transformed_all = flowData_transformed
#' flowData_transformed_subset = flowData_transformed[,param]
#'
#' # Normalize parameter values to [0,1] interval based on max. value across parameters
#' summary <- flowCore::fsApply(x=flowData_transformed_all, FUN=function(x) apply(x, 2, max), use.exprs=TRUE)
#' max = max(summary[,9])
#' mytrans <- function(x) x/max
#' flowData_transformed_all <- flowCore::transform(flowData_transformed_all,`FL1-H`=mytrans(`FL1-H`),
#'          `FL3-H`=mytrans(`FL3-H`),
#'          `SSC-H`=mytrans(`SSC-H`),
#'          `FSC-H`=mytrans(`FSC-H`))
#'
#' summary <- flowCore::fsApply(x=flowData_transformed_subset, FUN=function(x) apply(x,2,max), use.exprs=TRUE)
#' max = max(summary[,1])
#' mytrans <- function(x) x/max
#' flowData_transformed_subset <- flowCore::transform(flowData_transformed_subset,`FL1-H`=mytrans(`FL1-H`),
#'          `FL3-H`=mytrans(`FL3-H`),
#'          `SSC-H`=mytrans(`SSC-H`),
#'          `FSC-H`=mytrans(`FSC-H`))
#'
#' # Calculate diversity for first 5 samples without cleaning data
#' Diversity_rf(flowData_transformed_subset[1:5], param = param, R = 3, R.b = 3)
#'
#' # Calculate diversity for first 5 samples with cleaning data if needed
#' # Requires flowSet with all parameters included. Only samples with more than
#' # 30,000 cells will be evaluated.
#' Diversity_rf(flowData_transformed_all[1:5], param = param, R = 3, R.b = 3,
#' cleanFCS = TRUE, cleanparam = c(9,11))

#' @export

Diversity_rf_test <- function(x, d = 4, R = 100, R.b = 100, bw = 0.01, nbin = 128, 
                              param, cleanFCS = FALSE, cleanparam = c(9,11), ncores,
                              parallel = FALSE) {
  
  if(cleanFCS == TRUE){
    cat(date(), paste0("--- Using the following parameters for removing errant collection events\n in samples with > 30,000 cells: ", colnames(x)[cleanparam[1]], " ", 
                       colnames(x)[cleanparam[2]], "\n"))
    x <- flowCore::fsApply(x = x, FUN = function(x) FCS_clean(x, cleanparam))
    # Save parameters used to filter
    paramfilter <- colnames(x)[cleanparam]
    x <- x[,param]
    cat(date(), paste0("--- Done with cleaning data\n"))
  }
  if(parallel == FALSE){
    for (i in 1:R) {
      cat(date(), paste0("--- Starting resample run ", i, "\n"))
      tmp <- FCS_resample(x, rarefy = TRUE, replace = TRUE, progress = FALSE)
      tmp.basis <- flowFDA::flowBasis(tmp, param = param, nbin = nbin, bw = bw, 
                                      normalize = function(x) x)
      tmp.diversity <- Diversity(tmp.basis, plot = FALSE, d = d, R = R.b, 
                                 progress = FALSE)
      rm(tmp, tmp.basis)
      if (i == 1) 
        results <- cbind(tmp.diversity) else {
          results <- rbind(results, tmp.diversity)
        }
      rm(tmp.diversity)
    }
    
  } else {
    cl <- parallel::makeCluster(ncores)
    doParallel::registerDoParallel(cl)
    cat(date(), "--- Using", ncores, "cores for calculations\n")
    # log.socket <- make.socket(port = 4000)
    
    results <- foreach::foreach(i = 1:R, .combine = rbind, .packages = c("flowCore", "Phenoflow")) %dopar% {
      # cat(date(), paste0("--- Starting resample run ", i, "\n"))
      tmp <- FCS_resample(x, rarefy = TRUE, replace = TRUE, progress = FALSE)
      tmp.basis <- flowFDA::flowBasis(tmp, param = param, nbin = nbin, bw = bw, 
                                      normalize = function(x) x)
      tmp.diversity <- Phenoflow::Diversity(tmp.basis, plot = FALSE, d = d, R = R.b, 
                                            progress = FALSE)
    }
  }
  if(parallel == TRUE){
    cat(date(), "--- Closing workers\n")
    parallel::stopCluster(cl)
  }
  results.sd <- by(results[, c(2, 3, 5)], INDICES = factor(results$Sample_name), 
                   FUN = function(x) apply(x, 2, sd))
  results.sd <- do.call(rbind, results.sd)
  colnames(results.sd) <- c("sd.D0", "sd.D1", "sd.D2")
  results.m <- by(results[, c(2, 3, 5)], INDICES = factor(results$Sample_name), 
                  FUN = colMeans)
  results.m <- do.call(rbind, results.m)
  results <- data.frame(Sample_names = flowCore::sampleNames(x), results.m, 
                        results.sd)
  # Add parameters as attributes to dataframe
  attr(results, "R") <- R
  attr(results, "R.b") <- R.b
  attr(results, "bw") <- bw
  attr(results, "nbin") <- nbin
  attr(results, "d") <- d
  attr(results, "cleanFCS") <- cleanFCS
  if(cleanFCS == TRUE) attr(results, "cleanparam") <- paramfilter
  
  cat(date(), paste0("--- Alpha diversity metrics (D0,D1,D2) have been computed after ", 
                     R, " bootstraps\n"))
  return(results)
}