#' Random Forest classifier for supervised demarcation of groups using flow cytometry data.
#'
#' @param x flowSet object where the necessary metadata for classification is 
#' included in the phenoData.
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
#' @param cleanFCS Indicate whether outlier removal should be conducted prior to model estimation. 
#' Defaults to FALSE. I would recommend to make sure samples have > 500 cells. Will denoise based on the parameters specified in `param`.
#' @param timesplit Fraction of timestep used in flowAI for denoising. Please consult the `flowAI::flow_auto_qc` function for more information.
#' @param TimeChannel Name of time channel in the FCS files. This can differ between flow cytometers. Defaults to "Time". You can check this by: colnames(flowSet).
#' @param plot_fig Should the confusion matrix and the overall performance statistics on the test data partition be visualized?
#' Defaults to FALSE.
#' @importFrom BiocGenerics unique colnames
#' @importFrom flowAI flow_auto_qc 
#' @importFrom caret trainControl createDataPartition train confusionMatrix
#' @keywords random forest, fcm
#' @examples 
#' 
#' # 1. Example with environmental data:
#' 
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
#' print(model_pred)
#' 
#' # 2. Example with synthetic community data
#' # Load flow cytometry data of two strains with each 5,000 cells measured
#' data(flowData_ax)
#' 
#' # Quickly generate the necesary metadata
#' metadata_syn <- data.frame(name = flowCore::sampleNames(flowData_ax),
#'                        labels = flowCore::sampleNames(flowData_ax))
#' 
#' # Run Random forest model on 100 cells of each strain
#' model_rf_syn <- RandomF_FCS(flowData_ax, sample_info = metadata_syn, target_label = "labels",
#'                         downsample = 100, plot_fig = TRUE)
#'                         
#' # Make predictions on each of the samples or on new data of the mixed communities
#' model_pred_syn <- RandomF_predict(x = model_rf_syn[[1]], new_data =  flowData_ax, cleanFCS = FALSE)
#' print(model_pred_syn)

#' @export

RandomF_FCS <- function(x, sample_info, target_label, downsample = 0, 
                       classification_type = "sample",
                       param = c("FL1-H", "FL3-H", "FSC-H", "SSC-H"),
                       p_train = 0.75, seed = 777,
                       cleanFCS = FALSE,
                       timesplit = 0.1,
                       TimeChannel = "Time",
                       plot_fig = FALSE,
                       method = "rf") {
  # Set seed
  set.seed(seed)
  
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
                              ChFM = paste0(param_f[!param_f %in% 
                                                      c("FSC","SSC")],"-",
                                            add_measuredparam),
                              timeCh=TimeChannel,
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
  trainIndex <- caret::createDataPartition(full_data$label, p = p_train)
  
  train_data <- full_data[trainIndex$Resample1, ]
  test_data <- full_data[-trainIndex$Resample1, ]
  
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
  RF_train <- caret::train(label~., data = train_data, method = method, 
                           metric = metric, tuneGrid = tunegrid, 
                           trControl = fitControl, ntree = 500)
  print(RF_train)
  # cat(date(), paste0("--- Training Xgboost classifier on ",
  #                    100 * p_train,
  #                    "% training set with options: \n",
  #                    "\t-Performance metric: ", metric, "\n",
  #                    "\t-Number of trees: ", 500, "\n",
  #                    "\t-mtry: ", round(mtry, 2), "\n",
  #                    "\t-method: ", fitControl$number, "x ", fitControl$method, "\n",
  #                    "\t-repeats: ", fitControl$repeats, "\n"))
  # cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  # bstSparse <- xgboost(data = train_data$data, 
  #                      label = train_data$label, 
  #                      max.depth = 2, eta = 1, 
  #                      nthread = 2, nround = 2, 
  #                      objective = "binary:logistic")
  # print(bstSparse)

  
  # Step 6: Accuracy on test set
  cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  performance_metrics <- data.frame(metric = 1, 
                                    n_cells = c(table(test_data$label)), 
                                    label = levels(test_data$label))
  for(n_label in 1:length(unique(test_data$label))){
    tmp <- test_data[test_data$label == unique(test_data$label)[n_label], ]
    tmp_pred <- stats::predict(RF_train, newdata = tmp)
    index <- performance_metrics$label == unique(test_data$label)[n_label]
    performance_metrics$metric[index] <- 
      round(100*base::sum(tmp_pred == tmp$label)/performance_metrics$n_cells[index],2)
  }
  colnames(performance_metrics)[1] <- "Accuracy"
  
  cat(date(), paste0("--- Performance on ",
                     100 * (1 - p_train),
                     "% test set: \n"))
  print(performance_metrics)
  cat(paste0("-----------------------------------------------------------------------------------------------------\n"))
  
  # For confusion matrix:
  RF_pred <- stats::predict(RF_train, newdata = test_data)
  
  # Create dataframe containing descision boundary of final model
    ## Extract minimum and maximum values of FCS files
  # summary <- flowCore::fsApply(x = x, FUN = function(x) base::apply(x, 2, base::max), use.exprs = TRUE)
  # maxval <- base::apply(summary, 2, max)
  # minval <- base::apply(summary, 2, min)
  # desc_coord <- list()
  # for(i_param in 1:length(param[1:2])){
  #   desc_coord[[i_param]] <- seq(from = minval[i_param],
  #                                by = ((maxval[i_param] - minval[i_param])/64), 
  #                          to = maxval[i_param])
  # }
  #   ## Create vectors containing the parameter value bins for prediction
  # desc_pred <- base::expand.grid(desc_coord)
  # base::colnames(desc_pred) <- param
  # 
  #   ## Predict this dummy data
  # RF_pred_desc <- stats::predict(RF_train, newdata = desc_pred)
  # RF_pred_desc <- cbind(desc_pred, RF_pred_desc)
  # colnames(RF_pred_desc) <- c(param, "label")
  
  # Make list containing model, confusion matrix, summary statistics and descision boundary
  results_list <- list()
  results_list[[1]] <- RF_train
  results_list[[2]] <- caret::confusionMatrix(data = RF_pred, test_data$label)
  # results_list[[3]] <- RF_pred_desc

  # Return diagnostic plots (confusion matrix + descision boundary)
  if(plot_fig == TRUE){
    # Confusion matrix plot
    mytable <- base::round(data.frame(performance = results_list[[2]]$overall), 2)
    
    p_conf <- ggplot2::ggplot(data.frame(results_list[[2]]$table), 
                              ggplot2::aes(x = Prediction, y = Reference, fill = 100*Freq/sum(Freq)))+
      ggplot2::geom_raster()+
      ggplot2::geom_text(ggplot2::aes(label = round(100*Freq/sum(Freq), 0)), size = 6)+
      ggplot2::scale_fill_distiller( name = "% of total cells\n classified\n") +
      ggplot2::theme_bw()+
      ggplot2::scale_x_discrete(position = "top")+
      ggplot2::theme(axis.title=ggplot2::element_text(size=16), 
                     strip.text.x=ggplot2::element_text(size=14),
                     legend.title=ggplot2::element_text(size=14), 
                     legend.text=ggplot2::element_text(size=14),
                     axis.text.y =  ggplot2::element_text(size=13),
                     axis.text.x = ggplot2::element_text(size=13, angle = 55, hjust = 0),
                     title= ggplot2::element_text(size=20),
                     plot.margin = ggplot2::unit(c(1.1,1.1,1.1,1.1), "cm"),
                     panel.grid.major = ggplot2::element_blank(), 
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = "transparent",colour = NA),
                     plot.background =ggplot2::element_rect(fill = "transparent",colour = NA)
            
      )
      
    p_conf_table <- ggplot2::ggplot(data.frame(results_list[[2]]$table), 
                              ggplot2::aes(x = Prediction, y = Reference, fill = 100*Freq/sum(Freq)))+
      ggplot2::theme(axis.title=ggplot2::element_text(size=16), 
                     strip.text.x=ggplot2::element_text(size=14),
                     legend.title=ggplot2::element_text(size=14), 
                     legend.text=ggplot2::element_text(size=14),
                     axis.text.y =  ggplot2::element_text(size=13),
                     axis.text.x = ggplot2::element_text(size=13, angle = 55, hjust = 0),
                     title=ggplot2::element_text(size=20),
                     plot.margin = ggplot2::unit(c(1.1,1.1,1.1,1.1), "cm"),
                     panel.grid.major = ggplot2::element_blank(), 
                     panel.grid.minor = ggplot2::element_blank(),
                     panel.border = ggplot2::element_blank(),
                     panel.background = ggplot2::element_rect(fill = "transparent",colour = NA),
                     plot.background = ggplot2::element_rect(fill = "transparent",colour = NA)
      )+
      ggplot2::ylab("")+
      ggplot2::xlab("")+
      ggplot2::annotation_custom(gridExtra::tableGrob(mytable))
    
    print(cowplot::plot_grid(p_conf, p_conf_table, align = "h", ncol = 2, rel_widths = c(1/2, 1/5)))
    
    # Show positions of labelled cells for a stratified sample of the test data
    
    
    
    
    
    # Descision boundary plot
    # df_desc <- unique(results_list[[3]][, c(param[1:2], "label")]) # Retain only unique values on two primary parameters
    # p_desc <- ggplot2::ggplot(df_desc, ggplot2::aes(x = `FL1-H`, 
    #                                                  y = `FL3-H`,
    #                                                  fill = label))+
    #   geom_tile()+
    #   scale_fill_manual(values = c("red", "blue"))
    # Use grid.arrange to put them in one figure
    
  }
  
  # Final step: Return list with final model, confusion matrix, summary statistics
  # and descision boundary for further applications
  
  return(results_list)
}

