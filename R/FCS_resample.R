#' Resampling function for flowSet objects
#'
#' This function resamples each sample to an equal number of cells.
#' @param x flowSet object.
#' @param sample Desired sample size. Defaults to minimum sample size.
#' @param replace Do you want to resample with or without replacement? Defaults to FALSE, which is without replacement.
#' @keywords resampling, fcm
#' @export
#' @examples
#' FCS_resample()

FCS_resample <- function(x, sample=0, replace=FALSE){
  sample_distr <- data.frame(counts=fsApply(x,FUN=function(x) nrow(x),use.exprs=TRUE))
  p1 <- easyGgplot2::ggplot2.histogram(data=sample_distr , xName='counts',
                          fill="white", color="black",
                          linetype="longdash",addMeanLine=TRUE, meanLineColor="red",
                          meanLineType="dashed", meanLineSize=1)+
    theme_bw() + labs(y="Frequency", title="Original count distribution")
  if(sample==0) sample <- min(flowCore::fsApply(x=x,FUN=function(x) nrow(x),use.exprs=TRUE))
  ## Remove all .fcs files with less observations than the specified sample
  x <- x[fsApply(x=x,FUN=function(x) nrow(x),use.exprs=TRUE)>sample]
  for(i in 1:length(x)){
    flowCore::exprs(x[[i]]) <- flowCore::exprs(x[[i]])[base::sample(1:nrow(flowCore::exprs(x[[i]])), sample, replace=replace),]
  }
  print(p1)
  cat(paste0("Your samples were randomly subsampled to ",sample," cells"))
  return (x)
}