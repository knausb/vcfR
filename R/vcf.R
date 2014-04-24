# vcf.R.

#' Variant call format files processed with vcfR.
#'
#' vcfR
#' 
#' The complete list of functions can be displayed with
#' "library(help = vcfR)".
#'
#' Lubridate provides tools that make it easier to parse and 
#' manipulate dates. These tools are grouped below by common 
#' purpose. More information about each function can be found in 
#' its help documentation.
#'
#' 
#' @references Brian Knaus (2014). Variant call format files processed 
#' with vcfR. Journal TBA, NN(N),
#'   N-NN. \url{http://www.someurl.org/v40/i03/}.
#' @import ape
#' @docType package
#' @name vcfR
#' @rdname vcfR
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
#' @slot fix  data.frame for the fixed information
#' @slot gt   data.frame for the genotype information
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

#### Generic methods. ####

#### Method show ####
setMethod(
  f="show",
  signature = "vcfR",
  definition=function(object){
    cat("*** Class vcf, method Show *** \n")
    if(length(object@meta)>0){
      cat("Meta")
#      cat(head(object@meta))
      head(object@meta)
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

#### Mehtod print ####
setMethod(
  f="print",
  signature="vcfR",
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

#### Method head ####
#' @rdname vcfR-methods
#' @export
#' @aliases head.vcfR
#' 
setMethod(
  f="head",
  signature="vcfR",
  definition=function (x, y, n=6, ...){
    cat("***** Object of class 'vcf' *****\n")
    cat("***** Meta section *****\n")
    if(length(x@meta) > n){
      print(x@meta[1:n])
      cat("First ", n, " rows.\n")
    } else {
      print(x@meta)
    }
    #
    cat("\n***** Fixed section *****\n")
    if(nrow(x@fix) >= n){
      print(x@fix[1:n,1:7])
    } else {
      print(x@fix[,1:7])
    }
    #
    cat("\n***** Genotype section *****\n")
    if(nrow(x@gt) >= n){
       if(ncol(x@gt)<6){
         print(x@gt[1:n,])
       } else {
         print(x@gt[1:n,1:6])
         cat("First 6 columns only.\n")
       }
    } else {
      if(ncol(x@gt)<6){
        print(x@gt)
      } else {
        print(x@gt[,1:6])
      }
    }
    cat("\n")
    cat("Unique GT formats:\n")
    print(unique(as.character(x@gt[,1])))
    cat("\n")
#    cat("***** Head not implemented *****\n")
  }
)

#### Method plot ####
#' @rdname vcfR-methods
#' @export
#' @aliases plot.vcfR
#' 
#' @param ... Arguments to be passed to methods
#' 
setMethod(
  f= "plot",
  signature= "vcfR",
  definition=function (x,y,...){
#    cat("***** Object of class 'vcf' *****\n")
#    cat("***** Plot not implemented *****\n")
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


#### Data loading functions. ####

#' @title vcfR methods
#' @rdname vcfR-methods
# @aliases read.vcf write.vcf
#' @aliases read.vcf
#' @export
#'
#' @description
#' Reads in a vcf file and stores it in a vcf class.
#'
#' @param x variant call format (vcf) file
#'
#' @details
#' Reads in a vcf file and stores it in a vcf class.  Once the number of lines the meta information contains the data is divided into three tables: meta data, fixed data and genotype data.
#'
#' @examples
#' library(vcfR)
#' data(vcfR_example)
#' head(pinf_vcf)
#' plot(pinf_vcf)
#' pinf_vcf[1:6,]
#' 
read.vcf<-function(x){
  vcf <- new(Class="vcfR")
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
    vcf@fix <- read.table(x, header=T, sep='\t', skip=i, comment.char='', colClasses = "character")
    vcf@gt <- vcf@fix[,9:ncol(vcf@fix)]
    vcf@fix <- vcf@fix[,1:8]
    vcf@fix$POS  <- as.integer(vcf@fix$POS)
    vcf@fix$QUAL <- as.integer(vcf@fix$QUAL)
  } else if (length(grep('^#',tmp)) >0 ){
    # No meta region, but a header line.
    vcf@fix <- read.table(x,header=T,sep='\t',skip=0,comment.char='', colClasses = "character")
    vcf@gt <- vcf@fix[,9:ncol(vcf@fix)]
    vcf@fix <- vcf@fix[,1:8]
#    colnames(vcf@fix) <- c('chrom','pos','id','ref','alt','qual','filter','info')
    vcf@fix$POS  <- as.integer(vcf@fix$POS)
    vcf@fix$QUAL <- as.integer(vcf@fix$QUAL)
  } else {
    # No meta region or header line.
    vcf@fix <- read.table(x,header=F,sep='\t',skip=0,comment.char='', colClasses = "character")
    vcf@gt <- vcf@fix[,9:ncol(vcf@fix)]
    vcf@fix <- vcf@fix[,1:8]
    colnames(vcf@fix) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
    vcf@fix$POS  <- as.integer(vcf@fix$POS)
    vcf@fix$QUAL <- as.integer(vcf@fix$QUAL)    
  }
  return(vcf)
}

#' @rdname vcfR-methods
#' @aliases write.vcf
#' 
# @usage write.vcf(xvcf, vfile)
#' 
#' @param xvcf a vcfR object
#' @param vfile an output filename
#' @param mask logical vector indicating rows to use
#' 
#' @export
#' 
write.vcf<-function(xvcf, vfile, mask=logical(0)){
  if(class(xvcf) == 'Chrom'){
    # Recast as a vcfR object.
    temp <- xvcf
    xvcf <- new(Class="vcfR")
    xvcf@meta <- temp@vcf.meta
    xvcf@fix <- temp@vcf.fix
    xvcf@gt <- temp@vcf.gt
    mask <- temp@var.info$mask
    rm(temp)
  }
  if(class(xvcf) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or Chrom.")
  }
  #
  if(length(mask) == 0){
#    mask <- 1:nrow(xvcf@fix)
    mask <- rep(TRUE, times=nrow(xvcf@fix))
  }
  #
  orig_scipen <- getOption("scipen")
  options(scipen=999)
  header <- c(names(xvcf@fix), names(xvcf@gt))
  header[1] <- paste("#",header[1],sep='')
  write.table(xvcf@meta, file = vfile, append = FALSE, quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".", row.names = FALSE,
              col.names = FALSE)
  write(header, file = vfile,
        ncolumns=length(header),
        append = TRUE,
        sep = "\t")
  write.table(cbind(xvcf@fix[mask,], xvcf@gt[mask,]), file = vfile, append = TRUE,
              quote = FALSE, sep = "\t",
              eol = "\n", na = "NA", dec = ".",
              row.names = FALSE,
              col.names = FALSE)
#              col.names = TRUE)
  options(scipen=orig_scipen)
}

#### EOF. ####