#' Small (inefficient) function for extracting sample names associated with replicates
#'
#' @param x Vector of Sample names of entire flowSet
#' @param n Number of replicates
#' @keywords resampling, fcm
#' trip()

trip_col <- function(x,n=3){
  y = c()
  j=1
  for(i in seq(1,length(x),by=n)){
    y[j] = x[(i+n-1)]
    j=j+1
  }
  return(y)
}
