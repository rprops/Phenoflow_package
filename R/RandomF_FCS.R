#' Random Forest classifier for supervised demarcation of groups using flow cytometry data.
#'
#' @param x flowSet object where the necessary metadata for classification is included in the phenoData ac
#' @param sample_info Sample information necessary for the classification, has to contain a column named "name"
#' which matches the samplenames of the FCS files stored in the flowSet.
#' @param target_label Label that you should be predicted based on the flow cytometry data.
#' @param downsample Indicate to which sample size should be downsampled. 
#' By default samples are downsampled to the sample size of the sample with the lowest number of cells.
#' @param classification_type Whether classification is desired on "sample" level or on "single-cell" level. 
#' Defaults to sample level.
#' @param param Parameters to base classification on.
#' @param p_train Percentage of the data set that should be used for training the model.
#' @param seed Set random seed to be used during the analysis. Put at 777 by default.
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
#' downsample = 500)
#' @export

RandomF_FCS <- function(x, sample_info, target_label, downsample = 0, 
                       classification_type = "sample",
                       param = c("FL1-H", "FL3-H", "FSC-H", "SSC-H"),
                       p_train = 0.75, seed = 777) {
  # Set seed
  set.seed(seed)
  
  # Step 0: Format metadata
  pData(x) <- base::cbind(pData(x), 
                    sample_info[base::order(base::match(as.character(sample_info[, "name"]),
                                            as.character(pData(x)[, "name"]))), 
                                target_label])
  base::colnames(pData(x))[2] <- target_label
  
  # Step 1: Subset FCS files on parameters of interest
  x <- x[, param]
  
  # Step 1b: downsample
  x <- Phenoflow::FCS_resample(x, sample = downsample)
  
  # Step 2: add group labels present in pData to each cell and
  # merge all cells in a single data.table (faster than dataframe)
  if(classification_type == "sample"){
    for(i in 1:length(x)){
      tmp <- data.table::data.table(label = flowCore::pData(x[i])[, target_label], 
                                    flowCore::exprs(x[[i]]))
      if(i == 1){
        full_data <- tmp
      } else full_data <- base::rbind(full_data, tmp)
    }
  }
  full_data <- base::droplevels(full_data)
  
  # if(classification_type == "single-cell"){
  #   # Step 2: add sample labels to each cell
  #   flowCore::fsApply(x, use.exprs = TRUE)
  #   flowCore::pData(x)  
  # }
  
  # Step 3: Set model parameters
  fitControl <- caret::trainControl( ## 10-fold CV
    method = "repeatedcv",
    number = 10,
    ## repeated ten times
    repeats = 3,
    search = "random")
  
  # Step 4: Create data partitions
  train_data <- full_data[caret::createDataPartition(full_data$label,
                                                     p = p_train)[[1]], ]
  test_data <- full_data[caret::createDataPartition(full_data$label, 
                                                    p = (1-p_train))[[1]], ]
  
  # Step 5: Train Random Forest classifier on training set
  metric <- "Accuracy"
  mtry <- base::sqrt(ncol(train_data))
  tunegrid <- base::expand.grid(.mtry = mtry)
  RF_train <- caret::train(label~., data = train_data, method = "rf", 
                      metric = metric, tuneLength=15, 
                      trControl = fitControl)
  print(RF_train)
  plot(RF_train)
  
  # Step 6: Accuracy on test set
  cat(date(), paste0("--- Accuracy on ",
                     100*(1-p_train),
                     "% test set: ",  
                     100*round(base::sum(stats::predict(RF_train, newdata = test_data) == test_data$label)/
                       nrow(test_data),2), "%\n"))

  
  # Final step: Return model for further applications
  return(RF_train)
}

