% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RandomF_predict.R
\name{RandomF_predict}
\alias{RandomF_predict}
\title{Predict cell labels based on a model constructed using RandomF_FCS function.}
\usage{
RandomF_predict(
  x,
  new_data,
  cleanFCS = FALSE,
  param = c("FL1-H", "FL3-H", "FSC-H", "SSC-H"),
  timesplit = 0.1,
  TimeChannel = "Time"
)
}
\arguments{
\item{x}{Random forest model outputted from the RandomF_FCS function.}

\item{new_data}{flowSet containing the data to be predicted.}

\item{cleanFCS}{Indicate whether outlier removal should be conducted prior to model prediction.
Defaults to FALSE. I would recommend to make sure samples have > 500 cells. Will denoise based on the parameters specified in `param`.}

\item{param}{Parameters required to denoise the new_data}

\item{timesplit}{Fraction of timestep used in flowAI for denoising. Please consult the `flowAI::flow_auto_qc` function for more information.}

\item{TimeChannel}{Name of time channel in the FCS files. This can differ between flow cytometers. Defaults to "Time". You can check this by: colnames(flowSet).}
}
\description{
Predict cell labels based on a model constructed using RandomF_FCS function.
}
\examples{
# Load raw data (imported using flowCore)
data(flowData)

# Format necessary metadata
metadata <- data.frame(names = flowCore::sampleNames(flowData), 
do.call(rbind, lapply(strsplit(flowCore::sampleNames(flowData),"_"), rbind)))
colnames(metadata) <- c("name", "Cycle_nr", "Location", "day", 
"timepoint", "Staining", "Reactor_phase", "replicate")

# Run Random Forest classifier to predict the Reactor phase based on the
# single-cell FCM data
model_rf <- RandomF_FCS(flowData, sample_info = metadata, target_label = "Reactor_phase",
downsample = 10)

# Make a model prediction on new data and report contigency table of predictions
model_pred <- RandomF_predict(x = model_rf[[1]], new_data =  flowData[1], cleanFCS = FALSE)
}
\keyword{fcm}
\keyword{forest,}
\keyword{prediction,}
\keyword{random}
