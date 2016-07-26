#' Get elements from the fixed region of a VCF file
#' 
#' Both chromR objects and vcfR objects contain a region with fixed variables.
#' These accessors allow you to isolate these variables from these objects.
#' 
#' @param x a vcfR or chromR object
#' 
#' @return a vector or data frame
#' 
#' @rdname getFIX
#' @export
#' @aliases getFIX,chromR-method 
#'    getFIX,vcfR-method
#' @examples 
#' library("vcfR")
#' data("vcfR_example")
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
getFIX <- function(x) standardGeneric("getFIX")
#' @export
setGeneric("getFIX")

setMethod(
  f = "getFIX", 
  signature(x = "chromR"), 
  definition = function(x) {
    x@vcf@fix  
})

setMethod(
  f = "getFIX", 
  signature(x = "vcfR"), 
  definition = function(x) {
    x@fix  
})

#'
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



#'
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
    x@vcf@fix[,"POS"]
})

setMethod(
  f = "getPOS", 
  signature(x = "vcfR"), 
  definition = function(x) {
    x@fix[,"POS"]
})

#'
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
    x@vcf@fix[,"QUAL"]
})

setMethod(
  f = "getQUAL", 
  signature(x = "vcfR"), 
  definition = function(x) {
    x@fix[,"QUAL"]
})

#'
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

#'
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

#'
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

#'
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