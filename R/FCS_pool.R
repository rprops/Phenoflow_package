#' Function to pool FCS files based on sample name patterns
#'
#' This function pools samples based on sample name patterns
#' @param x flowSet object containing the flowframes to merge
#' @param stub A vector of patterns by which to pool FCS files (e.g. replicates). 
#' For example: providing the pattern "CYCLUS1_BEKKEN_0_4h_SYBR_START" will merge all 
#' flowframes containing this pattern in their sample name. 
#' **Notice**: You will have to provide a vector containing all the patterns to merge with
#' in case you want the full flowSet merged. For example: for stub = c("CYCLUS1_BEKKEN_0_4h_SYBR_START", 
#' "CYCLUS1_BEKKEN_0_10h_SYBR_START") `FCS_pool`` will **only** merge the samples containing those patterns, and it will
#' merge them per unique pattern.
#' @importFrom methods as
#' @importFrom flowCore sampleNames rbind2 flowSet fr_append_cols
#' @keywords pool, merge, FCS, preprocess
#' @examples 
#' 
#' @export
#' 
FCS_pool <- function(x, stub){
  if(length(stub) == length(x)) cat("-- No samples to merge --")
  for(i in 1:length(stub)){
    index <- grep(stub[i], flowCore::sampleNames(x))
    temp <- flowCore::flowSet(as(x[index], "flowFrame"))
    flowCore::sampleNames(temp) <- as.character(stub[i])
    if (length(index) == 1) {
      Original <- matrix(rep(1, dim(x[[index]])[1]))
      colnames(Original) <- "Original"
      temp <- fr_append_cols(x[[index]], cols = Original)
      temp <- as(temp, "flowSet")
      flowCore::sampleNames(temp) <- as.character(stub[i])
    }
    if(stub[i] == stub[1]){
      concat_x <- temp
    } else {
      concat_x <- flowCore::rbind2(concat_x, temp)
    }
  }
  return(concat_x)
}