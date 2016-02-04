#' 
#' @title Query the gt slot
#' @name query.gt
#' @rdname query_gt
#' 
#' @description Query the 'gt' slot of objects of class vcfR
#' 
#' 
#' @aliases is.polymorphic
#' 
#' @param x an object of class vcfR
#' @param na.omit logical to omit missing data
#' 
#' @details 
#' The function \strong{is_polymorphic} returns a vector of logicals indicating whether a variant is polymorphic.
#' Only variable sites are reported in vcf files.
#' However, once someone manipulates a vcfR object, a site may become invariant.
#' For example, if a sample is removed it may result in a site becoming invariant.
#' This function queries the sites in a vcfR object and returns a vector of logicals (TRUE/FALSE) to indicate if they are actually variable.
#' 
#' @seealso
#' \code{\link{extract.gt}}
#' 
#' 
#' @export
is.polymorphic <- function(x, na.omit=FALSE){
  if(class(x) != "vcfR"){
    stop("Expected an object of class vcfR")
  }
  x <- extract.gt(x)
  
  test.poly <- function(x, na.omit=na.omit){
    if(na.omit == TRUE){
      x <- stats::na.omit(x)
    }
    sum(x[1] == x[-1]) < (length(x) - 1)
  }
  apply(x, MARGIN=1, test.poly, na.omit=na.omit)
}


#' @rdname query_gt
#' @aliases is_biallelic
#' 
#' @details 
#' The function \strong{is_bialleleic} returns a vector of logicals indicating whether a variant is biallelic.
#' Some analyses or downstream analyses only work with biallelic loci.
#' This function can help manage this.
#' 
#' @export
is.biallelic <- function(x){
  #  x <- as.character(x@fix$ALT)
  x <- as.character(x@fix[,'ALT'])
  x <- strsplit(x, split=",")
  lapply(x, length) == 1
}
