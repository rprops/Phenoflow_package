#' Predict cell labels based on a model constructed using RandomF_FCS function.
#'
#' @param x Random forest model outputted from the RandomF_FCS function.
#' @param new_data flowSet containing the data to be predicted.
#' @param cleanFCS Indicate whether outlier removal should be conducted prior to model prediction.
#' Defaults to FALSE. I would recommend to make sure samples have > 500 cells. Will denoise based on the parameters specified in `param`.
#' @param param Parameters required to denoise the new_data
#' @param TimeChannel Name of time channel in the FCS files. This can differ between flow cytometers. Defaults to "Time". You can check this by: colnames(flowSet).
#' @param timesplit Fraction of timestep used in flowAI for denoising. Please consult the `flowAI::flow_auto_qc` function for more information. 
#' @importFrom flowAI flow_auto_qc
#' @keywords prediction, random forest, fcm
#' @examples 
#' # Load raw data (imported using flowCore)
#' data(flowData)
#' 
#' # Format necessary metadata
#' metadata <- data.frame(names = flowCore::sampleNames(flowData), 
#' do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData),"_"), rbind)))
#' colnames(metadata) <- c("name", "Cycle_nr", "Location", "day", 
#' "timepoint", "Staining", "Reactor_phase", "replicate")
#' 
#' # Run Random Forest classifier to predict the Reactor phase based on the
#' # single-cell FCM data
#' model_rf <- RandomF_FCS(flowData, sample_info = metadata, target_label = "Reactor_phase",
#' downsample = 10)
#' 
#' # Make a model prediction on new data and report contigency table of predictions
#' model_pred <- RandomF_predict(x = model_rf[[1]], new_data =  flowData[1], cleanFCS = FALSE)
#' @export

RandomF_predict <- function(x, new_data, cleanFCS = FALSE,
                            param = c("FL1-H", "FL3-H", "FSC-H", "SSC-H"),timesplit = 0.1,
                            TimeChannel = "Time") {
  if(cleanFCS == TRUE){
    cat(paste0("-------------------------------------------------------------------------------------------------", "\n"))
    cat(date(), paste0("--- Using the following parameters for removing errant collection events\n in samples:\n \n"))
    cat(paste0(param),"\n \n")
    cat(date(), paste0("--- Scatter parameters will be automatically excluded", "\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------"))
    cat("\n", paste0("Please cite:", "\n"))
    cat("\n", paste0("Monaco et al., flowAI: automatic and interactive anomaly discerning tools for flow cytometry data,\n Bioinformatics, Volume 32, Issue 16, 15 August 2016, Pages 2473-2480, \n https://doi.org/10.1093/bioinformatics/btw191", "\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------", "\n \n"))
    
    # Extract sample names
    sam_names <- flowCore::sampleNames(new_data) 
    
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
    new_data <- flowAI::flow_auto_qc(new_data, alphaFR = 0.01,
                              folder_results = "QC_flowAI",
                              fcs_highQ = "HighQ",
                              output = 1,
                              ChFM = paste0(param_f[!param_f %in% 
                                                      c("FSC","SSC")],
                                            "-", add_measuredparam),
                              timeCh=TimeChannel,
                              ChRemoveFS = filter_param,
                              second_fractionFR = timesplit
    )
    
    # Subset data to relevant data
    new_data <- new_data[, param]
    
    # Add sample names again since the QC function removes these
    new_data@phenoData@data$name  <- sam_names
    flowCore::sampleNames(new_data) <- sam_names
    
    # Change characters in parameter description from character back to numeric
    # Otherwise nothing that follows will work
    for(i in 1:length(new_data)){
      new_data[[i]]@parameters@data[,3] <- as.numeric(new_data[[i]]@parameters@data[,3])
      new_data[[i]]@parameters@data[,4] <- as.numeric(new_data[[i]]@parameters@data[,4])
      new_data[[i]]@parameters@data[,5] <- as.numeric(new_data[[i]]@parameters@data[,5])
      
    }
  }
  
  # Extract parameters used in the model construction
  param2 <- base::colnames(x$trainingData)[-1]
  new_data <- new_data[, param2]
  
  # Step 2: merge all cells in a single data.table (faster than dataframe)
  for(i in 1:length(new_data)){
    tmp <- data.table::data.table(flowCore::exprs(new_data[[i]]))
    tmp_labels <- base::rep(flowCore::sampleNames(new_data)[i], 
                            nrow(flowCore::exprs(new_data[[i]])))
    if(i == 1){
      full_data <- tmp
      full_samples <- tmp_labels
    } else {
      full_data <- base::rbind(full_data, tmp)
      full_samples <- base::c(full_samples, tmp_labels)
    }
  }
  full_data <- base::droplevels(full_data)
  
  # Step 3 run predictions
  new_data_pred <- stats::predict(x, newdata = full_data)
  res_data_pred <- data.table::data.table(Sample = full_samples, Predicted_label = new_data_pred) 
  
  # Step 4 output table with predictions
  return_df <- base::data.frame(base::table(res_data_pred))
  colnames(return_df)[3] <- "Counts"
  
  return(return_df)
}


