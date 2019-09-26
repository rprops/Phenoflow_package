#' Unsupervised calculation of phenotypic diversity on gaussian mixture-based phenotypic fingerprints
#' 
#' This function calculates Hill diversity metrics from the mixtures calculated by PhenoGMM.
#' D0 is calculated as Observed richness. D1 (exponential of Shannon entropy)
#' and D2 (inverse Simpson index) are respectively Hill order 1 and 2.
#' Errors for D0, D1 and D2 are calculated by bootstrapping.
#' 
#' @param gmm list containing 1) contigency table and 2) the gmm model created by 
#' the PhenoGMM function
#' @param R Number of bootstraps to conduct. Defaults to 100
#' @param verbose Progress of function is reported. Defaults to FALSE
#' @keywords diversity, fcm, alpha, gmm, PhenoGMM
#' @importFrom stats sd
#' @importFrom phyloseq rarefy_even_depth otu_table
#' @examples 
#' data(flowData_transformed)
#' # Make model on training data
#' testGMM <- PhenoGMM(flowData_transformed[1:2], downsample = 1e3, nG = 128, param = c("FL1-H", "FL3-H"))
#' # Apply model to unseen/new data
#' testPred <- PhenoMaskGMM(flowData_transformed, gmm = testGMM)
#' # Calculate diversity for both contigency tables
#' Diversity_gmm(testGMM)
#' Diversity_gmm(testPred)
#' @export

Diversity_gmm <- function(gmm, R = 100, verbose = FALSE) {
  gmm <- gmm[[1]]
  n_samples <-  nrow(gmm)
  # Matrix for storing data
  DIV <- matrix(nrow = n_samples, ncol = 6)
  row.names(DIV) <- gmm$Sample_names
  
  # Diversity functions
  D0.boot <- function(x) sum(x != 0)
  D1.boot <- function(x) exp(-sum(x * log(x)))
  D2.boot <- function(x) 1/sum((x)^2)
  
  # Start resampling
  for (i in 1:n_samples) {
    temp.D0 <- c()
    temp.D1 <- c()
    temp.D2 <- c()
    temp.gmm <- gmm[i, , drop = FALSE]
    if (verbose)
      cat(
        paste0(
          "\tCalculating diversity for sample ",
          i,
          "/",
          n_samples,
          " --- ",
          gmm$sample_names[i],
          "\n"
        )
      )
    for (j in 1:R) {
      temp <- rarefy_even_depth(
          otu_table(temp.gmm[,-1], taxa_are_rows = FALSE),
          verbose = FALSE,
          replace = TRUE
        )
      # Calculate frequencies
      temp <- temp/sum(temp)
      # Calculate Diversities
      temp.D0 <- c(temp.D0, D0.boot(temp))
      temp.D1 <- c(temp.D1, D1.boot(temp))
      temp.D2 <- c(temp.D2, D2.boot(temp))
      # Store diversities at the end of resampling run
      if (j == R) {
        DIV[i, 1] <- mean(temp.D0)
        DIV[i, 2] <- stats::sd(temp.D0)
        DIV[i, 3] <- mean(temp.D1)
        DIV[i, 4] <- stats::sd(temp.D1)
        DIV[i, 5] <- mean(temp.D2)
        DIV[i, 6] <- stats::sd(temp.D2)
        remove(temp.D0, temp.D1, temp.D2)
      }
    }
  }
  DIV <- data.frame(Sample_names = row.names(DIV), DIV)
  colnames(DIV) = c("Sample_names", "D0", "sd.D0", "D1", "sd.D1", "D2", "sd.D2")
  cat(date(), "\tDone with all", n_samples, "samples\n")
  return(DIV)
}
