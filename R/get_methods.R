#' Get elements from the fixed region of a VCF file
#' 
#' Both chromR objects and vcfR objects contain a region with fixed variables.
#' These accessors allow you to isolate these variables from these objects.
#'   
#' @param x a vcfR or chromR object
#' @param getINFO logical specifying whether getFIX should return the INFO column
#'
#' @return a vector or data frame
#' @rdname getFIX
#' @export
#' @aliases getFIX,chromR-method getFIX,vcfR-method
#' @examples 
#' library("vcfR")
#' data("vcfR_example")
#' data("chromR_example")
# ' chrom <- create.chromR('sc50', seq=dna, vcf=vcf, ann=gff)
#' getFIX(vcf) %>% head
#' getFIX(chrom) %>% head
#' 
#' getCHROM(vcf) %>% head
#' getCHROM(chrom) %>% head
#' 
#' getPOS(vcf) %>% head
#' getPOS(chrom) %>% head
#' 
#' getID(vcf) %>% head
#' getID(chrom) %>% head
#' 
#' getREF(vcf) %>% head
#' getREF(chrom) %>% head
#' 
#' getALT(vcf) %>% head
#' getALT(chrom) %>% head
#' 
#' getQUAL(vcf) %>% head
#' getQUAL(chrom) %>% head
#' 
#' getFILTER(vcf) %>% head
#' getFILTER(chrom) %>% head
#' 
#' getINFO(vcf) %>% head
#' getINFO(chrom) %>% head
#' 
getFIX <- function(x, getINFO = FALSE) standardGeneric("getFIX")
#' @export
setGeneric("getFIX")

setMethod(
  f = "getFIX", 
  signature(x = "chromR"), 
  definition = function(x, getINFO = FALSE) {
    if(getINFO == TRUE){
      return(x@vcf@fix)
    } else {
      return(x@vcf@fix[,-8])
    }
})

setMethod(
  f = "getFIX", 
  signature(x = "vcfR"), 
  definition = function(x, getINFO = FALSE) {
    if(getINFO == TRUE){
      return(x@fix)
    } else {
      return(x@fix[,-8])
    }
})

#' @rdname getFIX
#' @export
#' @aliases getCHROM,chromR-method
#'    getCHROM,vcfR-method
getCHROM <- function(x) standardGeneric("getCHROM")
#' @export
setGeneric("getCHROM")

setMethod(
  f = "getCHROM", 
  signature(x = "chromR"), 
  definition = function(x) {
    x@vcf@fix[,"CHROM"]
})

setMethod(
  f = "getCHROM", 
  signature(x = "vcfR"), 
  definition = function(x) {
    x@fix[,"CHROM"]
})



#' @rdname getFIX
#' @export
#' @aliases getPOS,chromR-method
#'    getPOS,vcfR-method
getPOS <- function(x) standardGeneric("getPOS")
#' @export
setGeneric("getPOS")

setMethod(
  f = "getPOS", 
  signature(x = "chromR"), 
  definition = function(x) {
    as.integer(x@vcf@fix[,"POS"])
})

setMethod(
  f = "getPOS", 
  signature(x = "vcfR"), 
  definition = function(x) {
    as.integer(x@fix[,"POS"])
})

#' @rdname getFIX
#' @export
#' @aliases getQUAL,chromR-method
#'    getQUAL,vcfR-method
getQUAL <- function(x) standardGeneric("getQUAL")
#' @export
setGeneric("getQUAL")

setMethod(
  f = "getQUAL", 
  signature(x = "chromR"), 
  definition = function(x) {
    as.numeric(x@vcf@fix[,"QUAL"])
})

setMethod(
  f = "getQUAL", 
  signature(x = "vcfR"), 
  definition = function(x) {
    as.numeric(x@fix[,"QUAL"])
})

#' @rdname getFIX
#' @export
#' @aliases getALT,chromR-method
#'    getALT,vcfR-method
getALT <- function(x) standardGeneric("getALT")
#' @export
setGeneric("getALT")

setMethod(
  f = "getALT", 
  signature(x = "chromR"), 
  definition = function(x) {
    x@vcf@fix[,"ALT"]
})

setMethod(
  f = "getALT", 
  signature(x = "vcfR"), 
  definition = function(x) {
    x@fix[,"ALT"]
})

#' @rdname getFIX
#' @export
#' @aliases getREF,chromR-method
#'    getREF,vcfR-method
getREF <- function(x) standardGeneric("getREF")
#' @export
setGeneric("getREF")

setMethod(
  f = "getREF", 
  signature(x = "chromR"), 
  definition = function(x) {
    x@vcf@fix[,"REF"]
})

setMethod(
  f = "getREF", 
  signature(x = "vcfR"), 
  definition = function(x) {
    x@fix[,"REF"]
})

#' @rdname getFIX
#' @export
#' @aliases getID,chromR-method
#'    getID,vcfR-method
getID <- function(x) standardGeneric("getID")
#' @export
setGeneric("getID")

setMethod(
  f = "getID", 
  signature(x = "chromR"), 
  definition = function(x) {
    x@vcf@fix[,"ID"]
})

setMethod(
  f = "getID", 
  signature(x = "vcfR"), 
  definition = function(x) {
    x@fix[,"ID"]
})

#' @rdname getFIX
#' @export
#' @aliases getFILTER,chromR-method
#'    getFILTER,vcfR-method
getFILTER <- function(x) standardGeneric("getFILTER")
#' @export
setGeneric("getFILTER")

setMethod(
  f = "getFILTER", 
  signature(x = "chromR"), 
  definition = function(x) {
    x@vcf@fix[,"FILTER"]
})

setMethod(
  f = "getFILTER", 
  signature(x = "vcfR"), 
  definition = function(x) {
    x@fix[,"FILTER"]
})


#' @rdname getFIX
#' @export
#' @aliases getINFO,chromR-method
#'    getINFO,vcfR-method
getINFO <- function(x) standardGeneric("getINFO")
#' @export
setGeneric("getINFO")

setMethod(
  f = "getINFO", 
  signature(x = "chromR"), 
  definition = function(x) {
    x@vcf@fix[,"INFO"]
})

setMethod(
  f = "getINFO", 
  signature(x = "vcfR"), 
  definition = function(x) {
    x@fix[,"INFO"]
})

