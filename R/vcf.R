# vcf.R.
##### ##### ##### ##### #####
# Class definition.

#' @title vcf class
#'
#' @description
#' A class for storing vcf data.
#'
#' @details
#'
setClass(
  Class="vcf",
  representation=representation(
    meta="character",
    fix="data.frame",
    gt="data.frame"
  ),
  prototype=prototype(
    fix = data.frame(matrix(ncol=8, nrow=0, 
                            dimnames=list(c(),
                                          c('chrom','pos','id','ref','alt','qual','filter','info'))),
                       stringsAsFactors=FALSE)
  )
)

##### ##### ##### ##### #####
# Generic methods.

setMethod(
  f="show",
  signature = "vcf",
  definition=function(object){
    cat("*** Class vcf, method Show *** \n")
    if(length(object@meta)>0){
      cat("Meta\n")
      cat(head(object@meta))
      cat("\n\n")
    }
    if(length(object@fix)>0){
      print(head(object@fix)[,1:7])
      cat("Column 8 (info) omitted.")
      cat("\n")
    }
    cat("\n")
    cat("******* End Show (vcf) ******* \n")
  }
)

setMethod(
  f="print",
  signature="vcf",
  definition=function (x,y,...){
    cat("***** Object of class vcf *****\n")
    if(length(x@meta)>0){
      cat("Meta\n")
      cat(head(x@meta))
      cat("\n\n")
    }
    if(length(x@fix)>0){
      print(head(x@fix)[,1:7])
      cat("Column 8 (info) omitted.")
      cat("\n")
    }
    if(length(x@gt)>0){
      print(head(x@gt))
      cat("\n")
    }
    cat("***** End print (vcf) ***** \n")
  }
)

setMethod(
  f= "plot",
  signature= "vcf",
  definition=function (x,y,...){
    cat("***** Object of class 'vcf' *****\n")
    cat("***** Plot not implemented *****\n")
  }
)

##### ##### ##### ##### #####
# Data loading functions.

read.vcf<-function(x){
  vcf <- new(Class="vcf")
  i <- -1 # Line counter.
  j <- 0 # Success?
  tmp <- scan(x, what="character", sep="\n", skip=0, nlines=1, quiet=T, comment.char="")
  if(length(grep('^##', tmp)) >0 ){
    while(j == 0){
      i <- i+1
      tmp <- scan(x, what="character", sep="\n", skip=i, nlines=1, quiet=T, comment.char="")
      if(length(grep('^##',tmp)) == 0){j <- 1}
    }
    vcf@meta <- scan(x, what="character", sep="\n", skip=0, nlines=i, quiet=T, comment.char="")
    vcf@fix <- read.table(x,header=T,sep='\t',skip=i,comment.char='')
    vcf@gt <- vcf@fix[,9:ncol(vcf@fix)]
    vcf@fix <- vcf@fix[,1:8]
  } else if (length(grep('^#',tmp)) >0 ){
    # No meta region, but a header line.
    vcf@fix <- read.table(x,header=T,sep='\t',skip=0,comment.char='')
    vcf@gt <- vcf@fix[,9:ncol(vcf@fix)]
    vcf@fix <- vcf@fix[,1:8]
    colnames(vcf@fix) <- c('chrom','pos','id','ref','alt','qual','filter','info')
  } else {
    # No meta region or header line.
    vcf@fix <- read.table(x,header=F,sep='\t',skip=0,comment.char='')
    vcf@gt <- vcf@fix[,9:ncol(vcf@fix)]
    vcf@fix <- vcf@fix[,1:8]
    colnames(vcf@fix) <- c('chrom','pos','id','ref','alt','qual','filter','info')
  }
  return(vcf)
}

##### ##### ##### ##### #####
# EOF.
