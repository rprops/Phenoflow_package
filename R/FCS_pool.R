#' Function to pool FCS files based on sample name patterns
#'
#' This function pools samples based on sample name patterns
#' @param x flowSet object containing the flowframes to merge
#' @param stub A vector of patterns by which to pool FCS files (e.g. replicates). 
#' For example: providing the pattern "CYCLUS1_BEKKEN_0_4h_SYBR_START" will merge all 
#' flowframes containing this pattern in their sample name.
#' @keywords pool, merge, FCS, preprocess
#' @export
#' 
FCS_pool <- function(x, stub){
  if(length(stub) == length(x)) cat("-- No samples to merge --")
  for(i in 1:length(stub)){
    index <- grep(stub[i], flowCore::sampleNames(x))
    temp <- flowCore::flowSet(as(x[index], "flowFrame"))
    flowCore::sampleNames(temp) <- as.character(stub[i])
    if(stub[i] == stub[1]){
      concat_x <- temp
    } else {
      concat_x <- flowCore::rbind2(concat_x, temp)
    }
  }
  return(concat_x)
}