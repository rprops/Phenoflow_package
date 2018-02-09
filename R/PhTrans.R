#' Small function for normalizing channels based on maximum value of chosen
#' channel
#' 
#' @param channel FCM channel
#' @param max_v maximum value across chosen FCM channel
#' @keywords transformation, fcm
#' myTrans(`FL1-H`, max_v)
#' @export
#' 
PhTrans <- function(channel, max_v){
  channel_norm <- channel/max_v
  return(channel_norm)
}