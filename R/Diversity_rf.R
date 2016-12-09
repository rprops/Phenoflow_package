#' Diversity_rf function for FCM data
#'
#' This function calculates Hill diversity metrics from FCM data. This function differs from
#' the Diversity() function in that it resamples (with replacement) all individual samples and 
#' averages out the diversity over all subsamples. This function is recommended in case there there are
#' differences in sample size (nr. of cells). Analysis time is approximately 1 min/resample run.
#' @param x flowSet containing the samples to analyse.
#' @param d Rounding factor for density values. Defaults to 4.
#' @param R Number of resampling runs to conduct on individual samples. Defaults to 100
#' @param R.b Number of bootstraps to conduct on the fingerprint (requires less). Defaults to 100
#' @param bw Bandwidth used in the kernel density estimation. Defaults to 0.01 which is ideal for normalized
#' FCM data (i.e., all parameter values are within [0,1]).
#' @param nbin Resolution of the binning grid. Defaults to 128 bins which corresponds to a 128x128 binning grid.
#' @param param Parameter vector indicating on which parameters the diversity should be estimated. An example input would be:
#' c("FL1-H", "FL3-H", "SSC-H", "FSC-H").
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
#' param=c('FL1-H', 'FL3-H','SSC-H','FSC-H')
#' flowData_transformed = flowData_transformed[,param]
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
#'  # Isolate only the cellular information based on the polyGate1
#'  flowData_transformed <- flowCore::Subset(flowData_transformed, polyGate1)
#'  
#'  # Normalize parameter values to [0,1] interval based on max. value across parameters
#'  summary <- flowCore::fsApply(x=flowData_transformed,FUN=function(x) apply(x,2,max),use.exprs=TRUE)
#'  max = max(summary[,1])
#'  mytrans <- function(x) x/max
#'  flowData_transformed <- flowCore::transform(flowData_transformed,`FL1-H`=mytrans(`FL1-H`),
#'          `FL3-H`=mytrans(`FL3-H`), 
#'          `SSC-H`=mytrans(`SSC-H`),
#'          `FSC-H`=mytrans(`FSC-H`))
#'  
#'  # Calculate diversity
#'  Diversity_rf(flowData_transformed, param = param, R = 3, R.b = 3)
#' @export

Diversity_rf <- function(x, d = 4, R = 100, R.b = 100, bw = 0.01, nbin = 128, 
                         param) {
  for (i in 1:R) {
    cat(date(), paste0("---- Starting resample run ", i, "\n"))
    tmp <- FCS_resample(x, rarefy = TRUE, replace = TRUE, progress = FALSE)
    tmp.basis <- flowBasis(tmp, param = param, nbin = nbin, bw = bw, 
                           normalize = function(x) x)
    tmp.diversity <- Diversity(tmp.basis, plot = FALSE, d = d, R = R.b, 
                               progress = FALSE)
    if (i == 1) 
      results <- cbind(tmp.diversity) else {
        results <- rbind(results, tmp.diversity)
      }
    rm(tmp, tmp.basis, tmp.diversity)
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
  cat(date(), paste0("---- Alpha diversity metrics (D0,D1,D2) have been computed after ", 
                     R, " bootstraps\n"))
  return(results)
}