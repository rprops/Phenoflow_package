#' Time discretization function for FCM data
#'
#' This function subsets the input frames based on analysis time. Developped for analysis of on-line time series data.
#' @param x flowSet of FCM data
#' @param create Do you want to make a new folder for the new flowframes? Defaults to FALSE.
#' @param analysis.length Dataframe with $time defining the total analysis time for each sample
#' @param start Vector of length n(x) that indicates for each sample at what time point it should start discretisizing.
#' (e.g., first 10 minutes are irrelevant, start = 10*60s)
#' @param time.interval Bin size of each new FCS file. For example, time.interval = 10 will make new FCS files of 10 second intervals.
#' @keywords online, fcm, time series analysis FCM
#' @export
#' @examples
#' time.discretization()

time.discretization <- function(x,analysis.length,create=FALSE,start=0,time.interval){
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
      time.gate <- rectangleGate(filterId = "Time discretization", "Time" = c(bottom, top), "FL1-H" = c(0, 1))
      res <- 0.1
      flowData.temp <- Subset(x[j],time.gate)
      flowData.temp[[1]]@description$`$VOL` <- as.numeric(as.numeric(x[[j]]@description$`$VOL`)*(time.interval)/(analysis.length$time[j]))
      print(as.numeric(as.numeric(x[[j]]@description$`$VOL`)*(time.interval)/(analysis.length$time[j])))
      write.FCS(x=flowData.temp[[1]], filename=paste(i+(10*number),time.interval,paste(rownames(analysis.length)[j]), sep="_"), what="numeric")
    }
    setwd(old.wd)
  }
}
