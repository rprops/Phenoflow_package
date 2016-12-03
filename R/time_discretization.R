#' Time discretization function for FCM data
#'
#' This function subsets the input frames based on analysis time. Developped for analysis of on-line time series data.
#' @param x flowSet of FCM data
#' @param create Do you want to make a new folder for the new flowframes? Defaults to FALSE.
#' @param analysis.length Dataframe with $time defining the total analysis time for each sample
#' @param start Vector of length n(x) that indicates for each sample at what time point it should start discretisizing.
#' (e.g., first 10 minutes are irrelevant, start = 10*60)
#' @param time.interval Bin size of each new FCS file. For example, time.interval = 10 will make new FCS files of 10 second intervals.
#' @param time.step discrete unit of time used in the FCS files (FCS files don't register in seconds). This value is used to transform
#' the discrete time units into seconds.
#' @param height Bottom and top of the time window utilized (use minimum and maximum parameter intensity). Default on c(0,200)
#' @param trigger Primary parameter of signal detection (e.g. FL1-H for SYBR Green). This is used to construct the binning window with the height parameter.
#' @keywords online, fcm, time series analysis FCM
#' @examples
#' # To be added in the near future
#' @export

time_discretization <- function(x, analysis.length, create=FALSE, start=0, time.interval, height = c(0,200), trigger = "FL1-H",
                                time.step=0.1){
  x <- transform(x,`Time`=(`Time`- min(`Time`))*time.step)
  for(j in 1:length(x)){
    number <- max(round((round(analysis.length/time.interval,0)+1)/10,0))
    old.wd <- getwd()
    if(create) {
      dir.create(paste(strsplit(rownames(analysis.length)[j],".fcs")[[1]], paste(time.interval), sep="_"))
      setwd(paste(strsplit(rownames(analysis.length)[j],".fcs")[[1]], paste(time.interval), sep="_"))
    }
    if(is.integer(analysis.length$time[j]/time.interval)) teller <- analysis.length$time[j]/time.interval
    else teller <- round(analysis.length$time[j]/time.interval,0)+1
    if(length(start)>0) teller = teller - start[j]/time.interval
    if(max(start) == 0) {
      start<-c()
      start[1:length(x)]<-0
    }
    res <- 0
    for(i in 1:teller){
      bottom <- (i-1)*time.interval + res + start[j]
      top <- i*time.interval + start[j]
      time.gate <- flowCore::rectangleGate(filterId = "Time discretization", "Time" = c(bottom, top), trigger = height)
      res <- 0.1
      flowData.temp <- flowCore::Subset(x[j],time.gate)
      flowData.temp[[1]]@description$`$VOL` <- as.numeric(as.numeric(x[[j]]@description$`$VOL`)*(time.interval)/(analysis.length$time[j]))
      print(as.numeric(as.numeric(x[[j]]@description$`$VOL`)*(time.interval)/(analysis.length$time[j])))
      flowCore::write.FCS(x=flowData.temp[[1]], filename=paste(i+(10*number),time.interval,paste(rownames(analysis.length)[j]), sep="_"), what="numeric")
    }
    setwd(old.wd)
  }
}
