#' Diversity_rf function for FCM data
#'
#' This function calculates Hill diversity metrics from FCM data. This function 
#' differs from the Diversity() function in that it resamples (with replacement) 
#' all individual samples and averages out the diversity over all subsamples. 
#' This function is recommended in case there are differences in sample 
#' size (nr. of cells). Analysis time is approximately 10s/resample run (R) on
#' default settings for the flowData example.
#' @param x flowSet containing the samples to analyse.
#' @param d Rounding factor for density values. Defaults to 4.
#' @param R Number of resampling runs to conduct on individual samples. Defaults 
#'  to 100
#' @param R.b Number of bootstraps to conduct on the fingerprint (requires less). 
#'  Defaults to 100
#' @param bw Bandwidth used in the kernel density estimation. Defaults to 0.01 
#'  which is ideal for normalized FCM data (i.e., all parameter values are 
#'  within [0,1]).
#' @param nbin Resolution of the binning grid. Defaults to 128 bins which 
#'  corresponds to a 128x128 binning grid.
#' @param param Parameter vector indicating on which parameters the diversity 
#'  should be estimated. An example input would be: c("FL1-H", "FL3-H", "SSC-H", 
#'  "FSC-H"). In addition the first parameter in this vector will be used to 
#'  normalize the data. Please make sure this is the primary fluorescence 
#'  channel if available. For example FL1-H for SYBR Green stained bacterial 
#'  cells measured on an Accuri C6.
#' @param parallel Should the calculation be parallelized? Defaults to FALSE
#' @param ncores How many cores should be used in case of parallel computation?
#' Defaults to 1.
#' @param cleanFCS Indicate whether outlier removal should be conducted prior 
#'  to diversity assessment (flowAI package). Defaults to FALSE. I would 
#'  recommend to make sure samples have > 500 cells. Will denoise based on 
#'  the parameters specified in `param`.
#' @param timesplit Fraction of timestep used in flowAI for denoising. Please 
#'  consult the `flowAI::flow_auto_qc` function for more information.
#' @param TimeChannel Name of time channel in the FCS files. This can differ 
#'  between flow cytometers. Defaults to "Time". You can check this by: 
#'  colnames(flowSet).
#' @keywords diversity, fcm, alpha
#' @importFrom foreach %dopar%
#' @importFrom flowCore fsApply
#' @importFrom BiocGenerics unique colnames
#' @importFrom flowAI flow_auto_qc
#' @importFrom flowFDA flowBasis
#' @examples
#' # Full data processing example
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
#' 
#' # Create a PolygonGate for denoising the dataset
#' # Define coordinates for gate in sqrcut1 in format: c(x,x,x,x,y,y,y,y)
#' sqrcut1 <- matrix(c(8.75,8.75,14,14,3,7.5,14,3),ncol=2, nrow=4)
#' colnames(sqrcut1) <- c('FL1-H','FL3-H')
#' polyGate1 <- flowCore::polygonGate(.gate=sqrcut1, filterId = 'Total Cells')
#' 
#' # Gating quality check
#' flowViz::xyplot(`FL3-H` ~ `FL1-H`, data=flowData_transformed[1], filter=polyGate1,
#'        scales=list(y=list(limits=c(0,14)),
#'        x=list(limits=c(6,16))),
#'        axis = lattice::axis.default, nbin=125,
#'        par.strip.text=list(col='white', font=2, cex=2), smooth=FALSE)
#' 
#' # Isolate only the cellular information based on the polyGate1
#' flowData_transformed <- flowCore::Subset(flowData_transformed, polyGate1)
#' 
#' # Calculate diversity for first 5 samples without cleaning data
#' Diversity_rf(flowData_transformed[1:5], param = param, R = 3, R.b = 3,
#'              cleanFCS = FALSE)
#' @export

Diversity_rf <- function(x, d = 4, R = 100, R.b = 100, bw = 0.01, nbin = 128, 
                          param, cleanFCS = FALSE, ncores=1,
                          parallel = FALSE, 
                          timesplit = 0.1,
                          TimeChannel = "Time") {
  
  ### Normalizing ##############################################################
  summary_x <- flowCore::fsApply(x = x, FUN = function(x) apply(x, 2, max), 
                                 use.exprs = TRUE)
  max <- base::max(summary_x[, param[1]])
  if (round(max,0) > 1) { 
    cat(paste0("-------------------------------------------------------------------------------------------------", "\n"))
    cat(date(), paste0("--- Normalizing your FCS data based on maximum ", 
                       param[1]," value\n"))
    for (i in 1:length(param)) cat(paste0("--- Maximum ", param[i],
                                         " before normalizing: ", 
                                         round(base::max(summary_x[, param[i]]),
                                               2),"\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------", "\n"))
    for (pm in param) {
      if (pm != TimeChannel) {
        myTrans <- flowCore::transformList(pm, function(x) x/max)
        x <- flowCore::transform(x, myTrans)
      }
    }
    summary_x <- flowCore::fsApply(x = x, FUN = function(x) apply(x, 2, max), 
                                   use.exprs = TRUE)
    for (i in 1:length(param)) cat(paste0("--- Maximum ", param[i],
                                         " after normalizing: ", 
                                         round(base::max(summary_x[, param[i]]),
                                               2),"\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------", "\n \n"))
    ##############################################################################
  } else cat(paste0("--- parameters are already normalized at: ", 
                    base::max(summary_x[, param[1]]),"\n"))
  
  if (cleanFCS == TRUE) {
    cat(paste0("-------------------------------------------------------------------------------------------------", "\n"))
    cat(date(), paste0("--- Using the following parameters for removing errant collection events\n in samples:\n \n"))
    cat(paste0(param),"\n \n")
    cat(date(), paste0("--- Scatter parameters will be automatically excluded", "\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------"))
    cat("\n", paste0("Please cite:", "\n"))
    cat("\n", paste0("Monaco et al., flowAI: automatic and interactive anomaly discerning tools for flow cytometry data,\n Bioinformatics, Volume 32, Issue 16, 15 August 2016, Pages 2473-2480, \n https://doi.org/10.1093/bioinformatics/btw191", "\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------", "\n \n"))
    
    # Extract parameters not to base denoising on
    param_f <- BiocGenerics::unique(gsub(param, 
                                         pattern = "-H|-A|-W", 
                                         replacement = ""))
    filter_param <- BiocGenerics::colnames(x)
    filter_param <- BiocGenerics::unique(gsub(filter_param, 
                                              pattern = "-H|-A|-W", 
                                              replacement = ""))
    filter_param <- filter_param[!filter_param %in% param_f & 
                                   filter_param != TimeChannel]
    filter_param <- c(filter_param, "FSC", "SSC")
    
    # Exclude all scatter information from denoising
    add_measuredparam <- base::unique(gsub("^.*-([A-Z])$","\\1",param))[1]
    
    # Denoise with flowAI
    x <- flowAI::flow_auto_qc(x, alphaFR = 0.01,
                              folder_results = "QC_flowA",
                              fcs_highQ = "HighQ",
                              output = 1,
                              timeCh = TimeChannel,
                              ChExcludeFM = paste0(param_f[param_f %in% c("FSC",
                                                                    "SSC")],
                                            "-", add_measuredparam),
                              ChExcludeFS = filter_param,
                              second_fractionFR = timesplit
    )
    
    # Subset data to relevant data
    x <- x[, param]
    
    # Change characters in parameter description from character back to numeric
    # Otherwise nothing that follows will work
    for (i in 1:length(x)) {
      x[[i]]@parameters@data[,3] <- as.numeric(x[[i]]@parameters@data[,3])
      x[[i]]@parameters@data[,4] <- as.numeric(x[[i]]@parameters@data[,4])
      x[[i]]@parameters@data[,5] <- as.numeric(x[[i]]@parameters@data[,5])
      
    }
    cat("\n", paste0("-----------------------------------------------------------------------------------------------------"), sep = "")
    cat("\n", date(), paste0(" --- Done with cleaning data\n"), sep = "")
  }
  
  # Subset data
  x <- x[, param]
  
  # Calculate diversity
  if (parallel == FALSE) {
    for (i in 1:R) {
      cat(date(), paste0("--- Starting resample run ", i, "\n"))
      tmp <- Phenoflow::FCS_resample(x, rarefy = TRUE, replace = TRUE, 
                                     progress = FALSE)
      tmp.basis <- flowFDA::flowBasis(tmp, param = param[param != TimeChannel], 
                                      nbin = nbin, bw = bw, 
                                      normalize = function(x) x)
      tmp.diversity <- Phenoflow::Diversity(tmp.basis, plot = FALSE, 
                                            d = d, R = R.b, 
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

    results <- foreach::foreach(i = 1:R, .combine = rbind,
                                .packages = c("flowCore", 
                                              "Phenoflow")) %dopar% {
      # cat(date(), paste0("--- Starting resample run ", i, "\n"))
      tmp <- Phenoflow::FCS_resample(x, rarefy = TRUE, replace = TRUE, 
                                     progress = FALSE)
      tmp.basis <- flowFDA::flowBasis(tmp, param = param[param != TimeChannel], 
                                      nbin = nbin, bw = bw, 
                                      normalize = function(x) x)
      tmp.diversity <- Phenoflow::Diversity(tmp.basis, plot = FALSE, 
                                            d = d, R = R.b, 
                                            progress = FALSE)
    }
  }
  if (parallel == TRUE) {
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
  results <- data.frame(Sample_names = rownames(results.m), results.m, 
                        results.sd)
  # Add parameters as attributes to dataframe
  attr(results, "R") <- R
  attr(results, "R.b") <- R.b
  attr(results, "bw") <- bw
  attr(results, "nbin") <- nbin
  attr(results, "d") <- d
  attr(results, "cleanFCS") <- cleanFCS
  if (cleanFCS) attr(results, "cleanparam") <- param
  
  cat(date(), paste0("--- Alpha diversity metrics (D0,D1,D2) have been computed after ", 
                     R, " bootstraps\n"))
  cat(paste0("-----------------------------------------------------------------------------------------------------", "\n \n"))
  return(results)
}
