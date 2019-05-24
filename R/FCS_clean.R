#' FCS_clean function for FCM data
#'
#' This function denoises a flowSet object using the flowAI package. 
#' Observations that #' do not need meet the criteria will be removed. 
#' Only works for samples with more than 1000 observations.
#' @param x flowSet containing the samples to analyse. 
#'  Cannot have been subsetted to  a lower number of parameters.
#' @param param Parameters to base removing of errant observartions on.
#' @param timesplit Fraction of timestep used in flowAI for denoising. 
#'  Please consult the `flowAI::flow_auto_qc` function for more information.
#' @param TimeChannel Name of time channel in the FCS files. This can differ 
#'  between flow cytometers. Defaults to "Time". You can check this by: 
#'  colnames(flowSet).
#' @importFrom flowAI flow_auto_qc
#' @importFrom BiocGenerics unique colnames
#' @keywords denoise, fcm, alpha
#' @examples
#' ## Full data processing example
#' data(flowData)
#' flowData_cleaned <- FCS_clean(flowData[1:3])
#' @export

FCS_clean <- function(x, 
                      param = c("FL1-H", "FL3-H", "FSC-H", "SSC-H"),
                      timesplit = 0.1,
                      TimeChannel = "Time")
  {
  # Extract sample names
  sam_names <- flowCore::sampleNames(x) 
    # Extract parameters not to base denoising on
  param_f <- BiocGenerics::unique(gsub(param, pattern = "-H|-A|-W", 
                                       replacement = ""))
  filter_param <- BiocGenerics::colnames(x)
  filter_param <- BiocGenerics::unique(gsub(filter_param, 
                                            pattern = "-H|-A|-W", 
                                            replacement = ""))
  filter_param <- filter_param[!filter_param %in% param_f & 
                                 filter_param!= TimeChannel]
  filter_param <- c(filter_param, "FSC", "SSC")
  # Exclude all scatter information from denoising
  add_measuredparam <- base::unique(gsub(".*-([A-Z])$","\\1",param))[1]
    # Denoise with flowAI
  x <- flowAI::flow_auto_qc(x, alphaFR = 0.01,
                            folder_results = "QC_flowAI",
                            fcs_highQ = "HighQ",
                            output = 1,
                            ChExcludeFM = paste0(param_f[!param_f %in% 
                                                    c("FSC","SSC")],"-",
                                          add_measuredparam),
                            timeCh=TimeChannel,
                            ChExcludeFS = filter_param,
                            second_fractionFR = timesplit
  )
  return(x)
}
