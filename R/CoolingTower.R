#' Flow cytometric fingerprint of cooling water microbial community
#'
#' Dataset containing the flow cytometric fingerprint of 41 samples 
#' of a secondary cooling water system, 128x128 binning grid, bw=0.01.
#' This is a subset of the total available dataset described in
#' doi: 10.1111/2041-210X.12607.
#'
#' @format fbasis object
#' \describe{
#'   \item{"@basis"}{Actual fingerprint - concatenated kernel density estimations}
#'   \item{"@param"}{Parameters used to generate the fingerprint}
#'   \item{"@pcom"}{Parameter combinations used to calculate kernel density 
#'   estimations}
#'   \item{"@nbin}{Number of bins used in kernel density estimation}
#'   \item{"@bw"}{Bandwidth used in kernel density estimation}
#' }
#' @source \url{http://dx.doi.org/10.5061/dryad.m1c04}
"CoolingTower"