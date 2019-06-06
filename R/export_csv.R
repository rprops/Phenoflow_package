#' Exports the flow cytometric fingerprint or raw data to CSV.
#'
#' This function exports the fingerprint to a single CSV file or the fcs file 
#' associated with a flowFrame object to a CSV file .
#' @param x flowCore::flowFrame or flowFDA::flowBasis object
#' @param location absolute or relative path to where the resulting csv files 
#' are written to. If it does not exist, it will be created (provided the user
#' has write persmissions in that directory). Defaults to working directory.
#' @importFrom utils write.csv
#' @importFrom methods .hasSlot
#' @importFrom flowCore keyword
#' @keywords machine learning, fcm
#' @examples
#' 
#' data(CoolingTower)
#' data(flowData)
#' ### show class info to see if flowSet (collection of flowFrames) or flowBasis
#' class(flowData)
#' ###  Export csv files of the raw data
#' fsApply(flowData,export_csv,location="./raw_csv/")
#' ###  Export csv files of the fingerprint
#' export_csv(CoolingTower,location="./fingerprint/")
#' @export

export_csv <- function(x, location="./"){
  if(!dir.exists(location)){dir.create(path=location,recursive = TRUE)}
  if(class(x)=="flowFrame"){
    fframe <- x
    rawinfo <- as.data.frame(fframe@exprs)
    fname <- sub(".fcs","",sub(" ","_",keyword(fframe)$`$FIL`)) 
    #should we remove the "," here as well?
    write.csv(x = rawinfo,
              file = paste0(location,
                            fname,".csv"))
  }else{
    if(class(x)=="flowBasis"){
      if(.hasSlot(x,"fp"))
      {
        if(any(dim(counts(x@fp))==0)){
          fingerpr <- x@basis
          fname <- "fingerprint"
          write.csv(x = fingerpr,
                    file = paste0(location,
                                  fname,".csv"))
        }else{
          fingerpr <- counts(x@fp)
          fname <- "fingerprint_probability"
          write.csv(x = fingerpr,
                    file = paste0(location,
                                  fname,".csv"))
          }
        }else{
          fingerpr <- x@basis
          fname <- "fingerprint"
          write.csv(x = fingerpr,
                    file = paste0(location,
                                  fname,".csv"))
      }
    } else {
      print("x is not a flowBasis or flowFrame object")
    }
  }
}


