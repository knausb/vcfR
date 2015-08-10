# vcf.R.

#' Variant call format files processed with vcfR.
#'
#' vcfR
#' 
#' The complete list of functions can be displayed with "library(help = vcfR)".
#'
#' Vignettes can be listed with: vignette(package='vcfR').
#'
#' Lubridate provides tools that make it easier to parse and 
#' manipulate dates. These tools are grouped below by common 
#' purpose. More information about each function can be found in 
#' its help documentation.
#'
#'
#' @references Brian J Knaus (2015). Variant call format files processed 
#' with vcfR. Journal TBA, NN(N),
#'   N-NN. \url{http://www.someurl.org/v40/i03/}.
#' @import ape
#' @docType package
#' @name vcfR
#' @rdname vcfR
#' @useDynLib vcfR
#' @importFrom Rcpp sourceCpp
NULL

#### Class definition. ####

#' @title vcfR class
#'
#' @name vcfR-class
#' @rdname vcfR-class
#'
#' @description
#' A class for storing vcf data.
#' 
#' @slot meta character vector for the meta (header) information
#' 
#' @slot fix  data.frame for the fixed information
#' @slot gt   data.frame for the genotype information 
# @slot fix  data.frame for the fixed information
# @slot gt   data.frame for the genotype information
#'
#' @details Defines a class for variant call format data.
#' A vcfR object contains three slots.  The first slot
#' is a character vector which holds the meta data.  The
#' second slot holds an eight column data.frame to hold the 
#' fixed data.  The third slot is a data.frame which holds
#' the genotype data.
#' @export
#' @import methods
setClass(
  Class="vcfR",
  representation=representation(
    meta="character",
    fix="matrix",
    gt="matrix"
#    fix="data.frame",
#    gt="data.frame"
  ),
  prototype=prototype(
#    fix = data.frame(matrix(ncol=8, nrow=0, 
#                            dimnames=list(c(),
#                                          c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'))),
#                     stringsAsFactors=FALSE)
    fix = matrix(ncol=8, nrow=0, 
                 dimnames=list(c(),
                               c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'))
                 )
  )
)

#### Generic methods. ####

#### Method show ####
setMethod(
  f="show",
  signature = "vcfR",
  definition=function(object){
    print("*** Class vcf, method Show ***")
    if(length(object@meta)>0){
      print("Meta")
      #      print(head(object@meta))
      head(object@meta)
#      print("\n\n")
      print("", quote=FALSE)
    }
    if(length(object@fix)>0){
      print(head(object@fix)[,1:7])
      print("Column 8 (info) omitted.")
#      print("\n")
      print("", quote=FALSE)
    }
    print("", quote=FALSE)
#    print("\n")
    print("******* End Show (vcf) *******")
  }
)

#### Method print ####
setMethod(
  f="print",
  signature="vcfR",
  #  definition=function (x,y,...){
  definition=function (x, ...){
    print("***** Object of class vcf *****")
    if(length(x@meta)>0){
      print("Meta")
      print(head(x@meta))
      print("", quote=FALSE)
#      print("\n\n")
    }
    if(length(x@fix)>0){
      print(head(x@fix)[,1:7])
      print("Column 8 (info) omitted.")
      print("", quote=FALSE)
#      print("\n")
    }
    if(length(x@gt)>0){
      print(head(x@gt))
      print("", quote=FALSE)
#      print("\n")
    }
    print("***** End print (vcf) *****")
  }
)

#### Method head ####
#' @rdname vcfR-methods
#' @export
#' @aliases head.vcfR
#' @title vcfR methods
#' 
#' @param x variant call format (vcf) file
#' @param n number of rows to print
#' 
setMethod(
  f="head",
  signature="vcfR",
  definition=function (x, n=6, ...){
    print("***** Object of class 'vcf' *****")
    print("***** Meta section *****")
    if(length(x@meta) > n){
      print(x@meta[1:n])
      print(paste("First", n, "rows."))
    } else {
      print(x@meta)
    }
    print("", quote=FALSE)
    #
    print("***** Fixed section *****")
    if(nrow(x@fix) >= n){
      print(x@fix[1:n,1:7])
    } else {
      print(x@fix[,1:7])
    }
    print("", quote=FALSE)
    #
    print("***** Genotype section *****")
    if(nrow(x@gt) >= n){
      if(ncol(x@gt)<6){
        print(x@gt[1:n,])
      } else {
        print(x@gt[1:n,1:6])
        print("First 6 columns only.")
      }
    } else {
      if(ncol(x@gt)<6){
        print(x@gt)
      } else {
        print(x@gt[,1:6])
      }
    }
    print("", quote=FALSE)
#    print("\n")
    print("Unique GT formats:")
    print(unique(as.character(x@gt[,1])))
    print("", quote=FALSE)
#    print("\n")
    #    print("***** Head not implemented *****\n")
  }
)

setGeneric("plot")
#### Method plot ####
#' @rdname vcfR-methods
#' @export
#' @aliases plot.vcfR
#' 
#' @param y not used
#' @param ... Arguments to be passed to methods
#' 
setMethod(
  f="plot",
  signature= "vcfR",
  definition=function(x, y, ...){
    #    print("***** Object of class 'vcf' *****\n")
    #    print("***** Plot not implemented *****\n")
    hist(x@fix$QUAL, col=5, main='Histogram of qualities', xlab='QUAL')
    rug(x@fix$QUAL)
  }
)

#### Method subset ####
#' @rdname vcfR-methods
#' @export
#' @aliases subset.vcf
#'
#' @param regex a regular expression to search for
#'
setMethod(
  f= "subset",
  signature= "vcfR",
  definition=function (x,regex,...){
    index <- grep(regex, x@fix[,1])
    x@gt <- x@gt[index,]
    x@fix <- x@fix[index,]
    x
  }
)


#### Method [] ####
#' @rdname vcfR-methods
#' @export
#' @aliases []
#'
#' @param i vector of rows (variants) to include
#' @param j vector of columns (samples) to include
#' @param drop delete the dimensions of an array which only has one level
#'
setMethod(
  f= "[",
  signature="vcfR",
  definition=function(x, i, j, drop){
    x@fix <- x@fix[i,]
    x@gt <- x@gt[i,j]
    return(x)
    #    if(i=="times"){return(x@times)}else {}
    #    if(i=="traj"){return(x@traj)}else {}
  }
)



#### EOF. ####



