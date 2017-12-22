#' Random Forest classifier for supervised demarcation of groups using flow cytometry data.
#'
#' @param x flowSet object where the necessary metadata for classification is 
#' included in the phenoData ac
#' @param sample_info Sample information necessary for the classification, has to 
#' contain a column named "name"
#' which matches the samplenames of the FCS files stored in the flowSet.
#' @param target_label column name of the sample_info dataframe that should be 
#' predicted based on the flow cytometry data.
#' @param downsample Indicate to which sample size should be downsampled. 
#' By default samples are downsampled to the sample size of the sample with the 
#' lowest number of cells.
#' Defaults to sample level.
#' @param param Parameters to base classification on.
#' @param p_train Percentage of the data set that should be used for training the model.
#' @param seed Set random seed to be used during the analysis. Put at 777 by default.
#' @param cleanFCS Indicate whether outlier removal should be conducted prior to diversity assessment (flowAI package). 
#' Defaults to FALSE. I would recommend to make sure samples have > 500 cells. Will denoise based on the parameters specified in `param`.
#' @param timesplit Fraction of timestep used in flowAI for denoising. Please consult the `flowAI::flow_auto_qc` function for more information.
#' @param TimeChannel Name of time channel in the FCS files. This can differ between flow cytometers. Defaults to "Time". You can check this by: colnames(flowSet).
#' @keywords resampling, fcm
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
#' RandomF_FCS(flowData, sample_info = metadata, target_label = "Reactor_phase",
#' downsample = 10)
#' @export

RandomF_FCS <- function(x, sample_info, target_label, downsample = 0, 
                       classification_type = "sample",
                       param = c("FL1-H", "FL3-H", "FSC-H", "SSC-H"),
                       p_train = 0.75, seed = 777,
                       cleanFCS = FALSE,
                       timesplit = 0.1,
                       TimeChannel = "Time") {
  # Set seed
  set.seed(seed)
  
  if(cleanFCS == TRUE){
    cat(paste0("-------------------------------------------------------------------------------------------------", "\n"))
    cat(date(), paste0("--- Using the following parameters for removing errant collection events\n in samples:\n \n"))
    cat(paste0(param),"\n \n")
    cat(date(), paste0("--- Scatter parameters will be automatically excluded", "\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------"))
    cat("\n", paste0("Please cite:", "\n"))
    cat("\n", paste0("Monaco et al., flowAI: automatic and interactive anomaly discerning tools for flow cytometry data,\n Bioinformatics, Volume 32, Issue 16, 15 August 2016, Pages 2473–2480, \n https://doi.org/10.1093/bioinformatics/btw191", "\n"))
    cat(paste0("-------------------------------------------------------------------------------------------------", "\n \n"))
    
    # Extract sample names
    sam_names <- flowCore::sampleNames(x) 
    
    # Extract parameters not to base denoising on
    param_f <- BiocGenerics::unique(gsub(param, pattern = "-H|-A", replacement = ""))
    filter_param <- BiocGenerics::colnames(x)
    filter_param <- BiocGenerics::unique(gsub(filter_param, pattern = "-H|-A", replacement = ""))
    filter_param <- filter_param[!filter_param %in% param_f & filter_param!= TimeChannel]
    filter_param <- c(filter_param, "FSC", "SSC")# Exclude all scatter information from denoising
    
    # Denoise with flowAI
    x <- flowAI::flow_auto_qc(x, alphaFR = 0.01,
                              folder_results = "QC_flowA",
                              fcs_highQ = "HighQ",
                              output = 1,
                              ChRemoveFS = filter_param,
                              second_fractionFR = timesplit
    )
    
    # Subset data to relevant data
    x <- x[, param]
    
    # Add sample names again since the QC function removes these
    x@phenoData@data$name  <- sam_names
    flowCore::sampleNames(x) <- sam_names
    
    # Change characters in parameter description from character back to numeric
    # Otherwise nothing that follows will work
    for(i in 1:length(x)){
      x[[i]]@parameters@data[,3] <- as.numeric(x[[i]]@parameters@data[,3])
      x[[i]]@parameters@data[,4] <- as.numeric(x[[i]]@parameters@data[,4])
      x[[i]]@parameters@data[,5] <- as.numeric(x[[i]]@parameters@data[,5])
      
    }
  }
  # Step 0: Format metadata
  Biobase::pData(x) <- base::cbind( Biobase::pData(x), 
                    sample_info[base::order(base::match(as.character(sample_info[, "name"]),
                                            as.character(Biobase::pData(x)[, "name"]))), 
                                target_label])
  base::colnames(Biobase::pData(x))[2] <- target_label
  
  # Step 1: Subset FCS files on parameters of interest
  x <- x[, param]
  
  # Step 1b: downsample
  cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  x <- Phenoflow::FCS_resample(x, sample = downsample)
  cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  
  # Step 2: add group labels present in pData to each cell and
  # merge all cells in a single data.table (faster than dataframe)
  for(i in 1:length(x)){
    tmp <- data.table::data.table(label = Biobase::pData(x[i])[, target_label], 
                                  flowCore::exprs(x[[i]]))
    if(i == 1){
      full_data <- tmp
    } else full_data <- base::rbind(full_data, tmp)
  }
  full_data <- base::droplevels(full_data)
  
  # if(classification_type == "single-cell"){
  #   # Step 2: add sample labels to each cell
  #   flowCore::fsApply(x, use.exprs = TRUE)
  #   Biobase::pData(x)  
  # }
  
  # Step 3: Set model parameters
  fitControl <- caret::trainControl( ## 10-fold CV
    method = "repeatedcv",
    number = 10,
    ## repeated ten times
    repeats = 3)
  
  # Step 4: Create data partitions
  train_data <- full_data[caret::createDataPartition(full_data$label,
                                                     p = p_train)[[1]], ]
  test_data <- full_data[caret::createDataPartition(full_data$label, 
                                                    p = (1 - p_train))[[1]], ]
  
  # Step 5: Train Random Forest classifier on training set
  metric <- "Accuracy"
  mtry <- round(base::sqrt(ncol(train_data)), 0)
  tunegrid <- base::expand.grid(.mtry = mtry)
  cat(date(), paste0("--- Training Random Forest classifier on ",
                     100 * p_train,
                     "% training set with options: \n",
                     "\t-Performance metric: ", metric, "\n",
                     "\t-Number of trees: ", 500, "\n",
                     "\t-mtry: ", round(mtry, 2), "\n",
                     "\t-method: ", fitControl$number, "x ", fitControl$method, "\n",
                     "\t-repeats: ", fitControl$repeats, "\n"))
  cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  RF_train <- caret::train(label~., data = train_data, method = "rf", 
                      metric = metric, tuneGrid = tunegrid, 
                      trControl = fitControl, ntree = 500)
  print(RF_train)
  
  # Step 6: Accuracy on test set
  cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  performance_metrics <- data.frame(metric = 1, 
                                    n_cells = c(table(test_data$label)), 
                                    label = levels(test_data$label))
  for(n_label in 1:length(unique(test_data$label))){
    tmp <- test_data[test_data$label == unique(test_data$label)[n_label], ]
    tmp_pred <- stats::predict(RF_train, newdata = tmp)
    performance_metrics$metric[n_label] <- round(100*base::sum(tmp_pred == tmp$label)/nrow(test_data),2)
  }
  colnames(performance_metrics)[1] <- metric
  
  cat(date(), paste0("--- Accuracy on ",
                     100 * (1 - p_train),
                     "% test set: \n"))
  print(performance_metrics)
  cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  
  # Final step: Return model for further applications
  return(RF_train)
}

