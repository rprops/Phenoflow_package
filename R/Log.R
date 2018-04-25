#' Small function for printing progress of parallel processes
#'
#' @param text Input text
#' @param ... additional parameters passed on to sprintf
#' @keywords resampling, fcm
#' Log(x)

Log <- function(text, ...) {
  msg <- sprintf(paste0(as.character(Sys.time()), ": ", text, "\n"), ...)
  cat(msg)
  write.socket(log.socket, msg)
}