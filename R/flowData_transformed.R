#' Flow cytometry test data of a cooling water microbial community, after transformation
#'
#' Dataset containing the multiparametric flow cytometry results of 41 samples.
#' This is a subset of the total available dataset described in
#' doi: 10.1111/2041-210X.12607.
#' Available parameters: FSC-A SSC-A FL1-A FL2-A FL3-A FL4-A FSC-H SSC-H 
#' FL1-H FL2-H FL3-H FL4-H Width Time
#' It has been transformed as outlined in  https://github.com/rprops/Phenoflow_package/wiki/1.-Phenotypic-diversity-analysis
#' In brief, the FSC-H, SSC-H, FL1-H and FL3-H parameters were asinh transformed
#' Subsequently, only cellular information was retained using a gating template
#' And finally the data was rescaled to the 0-1 range using the maximum value for the FL1-H fluoresecence.
#' @format A flowSet with 41 experiments
#' @source \url{http://dx.doi.org/10.5061/dryad.m1c04}
"flowData_transformed"