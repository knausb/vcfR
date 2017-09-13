#' @title Query the META section of VCF data
#' @name queryMETA
#' @rdname queryMETA
#' 
#' @description
#' Query the META section of VCF data for information about acronyms.
#'  
#' @param x an object of class vcfR or chromR.
#' @param element an acronym to search for in the META portion of the VCF data.
#' @param nice logical indicating whether to format the data in a 'nice' manner.
#' 
#' @details 
#' The META portion of VCF data defines acronyms that are used elsewhere in the data.
#' In order to better understand these acronyms they should be referenced.
#' This function facilitates looking up of acronyms to present their relevant information.
#' When 'element' is 'NULL' (the default), all acronyms from the META region are returned.
#' When 'element' is specified an attempt is made to return information about the provided element.
#' The function \code{grep} is used to perform this query.
#' If 'nice' is set to FALSE then the data is presented as it was in the file.
#' If 'nice' is set to TRUE the data is processed to make it appear more 'nice'.
#' 
#' @seealso \code{\link[base]{grep}}, \code{\link[base]{regex}}.
#' 
#' @examples
#' data(vcfR_test)
#' queryMETA(vcfR_test)
#' queryMETA(vcfR_test, element = "DP")
#' 
#' 
#' @export
#' 
queryMETA <- function(x, element = NULL, nice = TRUE){
  
  if( inherits(x, "chromR") ){
    x <- x@vcfR
  }  
  
  if( is.null(element) ){
    ID <- grep("=<ID=", x@meta, value = TRUE)
    ID <- grep("contig=<ID", ID, value = TRUE, invert = TRUE)
    
    if( nice ){
      ID <- nice(ID)
      ID <- lapply( ID, function(x){ x[1] } )
      ID <- unlist(ID)
    }
    myContigs <- grep("contig=<ID", x@meta)
    if( length(myContigs) > 0 ){
      ID <- c(ID, paste(length(myContigs), "contig=<IDs omitted from queryMETA"))
    }
    return(ID)
  }
  
  ID <- grep(element, x@meta, value = TRUE)
  if( nice ){
    ID <- nice(ID)
  }
  
  return(ID)
}


nice <- function(x){
  x <- sub("^##", "", x)
  x <- sub("<", "", x)
  x <- sub(">$", "", x)
  x <- sub("\"", "", x)
  x <- sub("\"$", "", x)
  x <- strsplit(x, split = ",")
  return(x)
}


