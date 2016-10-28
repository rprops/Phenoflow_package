#' Small (inefficient) function for taking mean/sd
#'
#' @param x Input matrix
#' @param n Number of replicates
#' @keywords resampling, fcm
#' @export
#' @examples
#' trip()

trip <- function(x,n=3){
  y = c(); s = c()
  j=1
  for(i in seq(1,length(x),by=n)){
    y[j] = mean(x[i:(i+n-1)])
    s[j] = sd(x[i:(i+n-1)])
    j=j+1
  }
  y=cbind(y,s)
  return(y)
}