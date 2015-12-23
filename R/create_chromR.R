#' @title Create chromR object
#' @name create.chromR
#' @rdname create_chromR
#' @export
#' @aliases create.chromR
#'
#' @description
#' Creates and populates an object of class chromR.
#'
#' @param name a name for the chromosome (for plotting purposes)
#' @param seq a sequence as a DNAbin object
#' @param ann an annotation file (gff-like)
#' @param verbose should verbose output be printed to the console?
#' @param x an object of class chromR
#' @param vcf an object of class vcfR
#' @param gff a data.frame containing annotation data in the gff format
# @param ... arguments
#'
#' @details
#' Creates and names a chromR object from a name, a chromosome (an ape::DNAbin object), variant data (a vcfR object) and annotation data (gff-like).
#' The function \strong{create.chromR} is a wrapper which calls functions to populate the slots of the chromR object.
#' 
#' The function \strong{vcf2chromR} is called by create.chromR and transfers the data from the slots of a vcfR object to the slots of a chromR object.
#' It also tries to extract the 'DP' and 'MQ' fileds (when present) from the fix slot's INFO column.
#' It is not anticipated that a user would need to use this function directly, but its placed here in case they do.
#' 
#' The function \strong{seq2chromR} is currently defined as a generic function.
#' This may change in the future.
#' This function takes an object of class DNAbin and assigns it to the 'seq' slot of a chromR object.
#' 
#' The function \strong{ann2chromR} is called by create.chromR and transfers the information from a gff-like object to the 'ann' slot of a chromR object.
#' It is not anticipated that a user would need to use this function directly, but its placed here in case they do.
#' 
#' 
#' @seealso 
# \code{\link{seq2chromR}},
# \code{\link{vcf2chromR}},
#' \code{\link{chromR-class}},
#' \code{\link{vcfR-class}},
#' \code{\link[ape]{DNAbin}},
#' \href{http://www.1000genomes.org/wiki/analysis/variant\%20call\%20format/vcf-variant-call-format-version-41}{vcf format}, 
#' \href{http://www.sequenceontology.org/gff3.shtml}{gff3 format}
#' 
#' @examples
#' library(vcfR)
#' data(vcfR_example)
#' chrom <- create.chromR('sc50', seq=dna, vcf=vcf, ann=gff)
#' head(chrom)
#' chrom
# colnames(chrom)
#' plot(chrom)
#' 
#' chrom <- masker(chrom, min_QUAL = 1, min_DP = 300, max_DP = 700, min_MQ = 59, max_MQ = 61)
#' chrom <- proc.chromR(chrom, win.size=1000)
#'
#' plot(chrom)
# 
#' chromoqc(chrom)
# chromoqc(pinf_mt, xlim=c(25e+03, 3e+04), dot.alpha=99)
# 
# set.seed(10)
# x1 <- as.integer(runif(n=20, min=1, max=39000))
# y1 <- runif(n=length(x1), min=1, max=100)
# chromodot(pinf_mt, x1=x1, y1=y1)
#' 
#           1         2         3         4         5
#  12345678901234567890123456789012345678901234567890
# chromodot(pinf_mt, x1=x1, y1=y1, label1='My data',
#           x2=x1, y2=y1, label2='More of my data',
#           dot.alpha='ff')
#' 
# chromohwe(pinf_mt, dot.alpha='ff')
#' 
# chromopop(pinf_mt)
# gt <- extract.gt(pinf_mt)
# head(gt)
# tab <- variant.table(pinf_mt)
# head(tab)
# win <- window_table(pinf_mt)
# head(win)
# hist(tab$Ho - tab$He, col=5)
# # Note that this example is a mitochondrion, so this is a bit silly.
#' 
create.chromR <- function(name="CHROM1", vcf, seq=NULL, ann=NULL, verbose=TRUE){
  # Determine whether we received the expected classes.
  stopifnot(class(vcf) == "vcfR")

  # Initialize chromR object.  
  x <- new(Class="chromR")
  setName(x) <- name
  
  # Insert vcf into Chom.
  if(length(vcf)>0){
#    x <- vcf2chromR(x, vcf)
    x@vcf <- vcf
  }

  # Insert seq into chromR
  # Needs to handle lists and matrices of DNAbin
  # Matrices are better behaved.
  #
  if(is.null(seq)){
    POS <- getPOS(x)
    x@len <- POS[length(POS)]
#    x@len <- x@vcf.fix$POS[length(x@vcf.fix$POS)]
  } else if (class(seq)=="DNAbin"){
    x <- seq2chromR(x, seq)
  } else {
    stopifnot(class(seq)=="DNAbin")
  }

  # Annotations.
  if(!is.null(ann)){
#  if(nrow(ann) > 0){
    stopifnot(class(ann) == "data.frame")
    if(class(ann[,4]) == "factor"){ann[,4] <- as.character(ann[,4])}
    if(class(ann[,5]) == "factor"){ann[,5] <- as.character(ann[,5])}
    if(class(ann[,4]) == "character"){ann[,4] <- as.numeric(ann[,4])}
    if(class(ann[,5]) == "character"){ann[,5] <- as.numeric(ann[,5])}
    x@ann <- ann
  }

  # Report names of objects to user.
  if(verbose == TRUE){
    # Print names of elements to see if they match.
    message("Names in vcf:")
    chr_names <- unique(getCHROM(x))
    message(paste('  ', chr_names, sep=""))
#    message(paste('  ', unique(as.character(x@vcf.fix$CHROM)), sep=""))
    
    if(class(x@seq) == "DNAbin"){
      message("Names of sequences:")
      message(paste('  ', unique(labels(x@seq)), sep=""))

#      if(unique(as.character(x@vcf.fix$CHROM)) != unique(labels(x@seq))){
      if(chr_names != unique(labels(x@seq))){
        warning("
        Names in variant data and sequence data do not match perfectly.
        If you choose to proceed, we'll do our best to match the data.
        But prepare yourself for unexpected results.")
#        message("Names in variant file and sequence file do not match perfectly.")
#        message("If you choose to proceed, we'll do our best to match data.")
#        message("But prepare yourself for unexpected results.")
      }
    }

    if(nrow(x@ann) > 0){
      message("Names in annotation:")
      message(paste('  ', unique(as.character(x@ann[,1])), sep=""))
#      if(unique(as.character(x@vcf.fix$CHROM)) != unique(as.character(x@ann[,1]))){
      if(chr_names != unique(as.character(x@ann[,1]))){
        warning("
        Names in variant data and annotation data do not match perfectly.
        If you choose to proceed, we'll do our best to match the data.
        But prepare yourself for unexpected results.")
#        message("Names in variant file and annotation file do not match perfectly.\n")
#        message("If you choose to proceed, we'll do our best to match data.\n")
#        message("But prepare yourself for unexpected results.\n")
      }
    }
  }
  
  # Check to see if annotation positions exceed seq position.
  if( nrow(x@ann) > 0 ){
    if( max(as.integer(as.character(x@ann[,4]))) > x@len | max(as.integer(as.character(x@ann[,5]))) > x@len ){
      stop("Annotation positions exceed chromosome positions.  Is this the correct set of annotations?")
    }
  }
  
  if( verbose == TRUE ){
    message("Initializing var.info slot.")
  }
  x@var.info <- data.frame( CHROM = x@vcf@fix[,"CHROM"] , POS = as.integer(x@vcf@fix[,"POS"]) )
#  mq <- getINFO(x, element="MQ")
  mq <- extract.info(x, element = 'MQ', as.numeric = TRUE)
  if( length(mq) > 0 ){ x@var.info$MQ <- mq }
#  dp <- getDP(x)
  dp <- extract.info(x, element = 'DP', as.numeric = TRUE)
  if( length(dp) > 0 ){ x@var.info$DP <- dp }
  x@var.info$mask <- TRUE
  if( verbose == TRUE ){
    message("var.info slot initialized.")
  }
  return(x)
}



##### ##### ##### ##### #####
#
# chromR data loading functions
#
##### ##### ##### ##### #####

#' @rdname create_chromR
#' @export
#' @aliases chromR-methods vcf2chromR
#'
# @description
# Methods to work with objects of the chromR class
# Reads in a vcf file and stores it in a vcf class.
#'
# @param x an object of class chromR
#'
#'
vcfR2chromR <- function(x, vcf){
  x@vcf.fix <- as.data.frame(vcf@fix)
#  colnames(x@vcf.fix) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
#  x@vcf.fix[,2] <- as.numeric(x@vcf.fix[,2])
#  x@vcf.fix[,6] <- as.numeric(x@vcf.fix[,6])
  #
  for(i in 1:ncol(vcf@gt)){
    vcf@gt[,i] <- as.character(vcf@gt[,i])
  }
  x@vcf.gt <- vcf@gt
  #
  x@vcf.meta <- vcf@meta
  #
  # Initialize var.info slot
  x@var.info <- data.frame(matrix(ncol=5, nrow=nrow(vcf@fix)))
  names(x@var.info) <- c('CHROM', 'POS', 'mask', 'DP','MQ')
#  names(x@var.info) <- c('DP','MQ', 'mask')
  #
  x@var.info$CHROM <- x@vcf.fix$CHROM
  x@var.info$POS <- x@vcf.fix$POS
  x@var.info$mask <- rep(TRUE, times=nrow(x@vcf.fix))
  #
  if(length(grep("DP=", vcf@fix[,8])) > 0){
    x@var.info$DP <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(vcf@fix[,8]), ";"), function(x){grep("^DP=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
  }
  if(length(grep("MQ=", vcf@fix[,8])) > 0){
    x@var.info$MQ <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(vcf@fix[,8]), ";"), function(x){grep("^MQ=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
  }
  #
  # assign may be more efficient.
  return(x)
}


# Needs to handle lists and matrices of DNAbin.
# Matrices appear better behaved.
#
#' @rdname create_chromR
#' @export
#' @aliases seq2chromR
#' 
seq2chromR <- function(x, seq=NULL){
  # A DNAbin will store in a list when the fasta contains
  # multiple sequences, but as a matrix when the fasta
  # only contains one sequence.
  if(is.list(seq)){
    stopifnot(length(seq)==1)
    x@seq <- as.matrix(seq)
    x@len <- length(x@seq)
  } else if (is.matrix(seq)){
    stopifnot(nrow(seq)==1)
#    x@seq <- ape::as.DNAbin(as.character(seq)[1,])
#    dimnames(pinf_dna)[[1]][1]
    x@seq <- seq
    x@len <- length(x@seq)
  } else {
    stop("DNAbin is neither a list or matrix")
  }
  return(x)
}


#' @rdname create_chromR
#' @export
#' @aliases ann2chromR
#'
ann2chromR <- function(x, gff){
  x@ann <- as.data.frame(gff)
  colnames(x@ann) <- c('seqid','source','type','start','end','score','strand','phase','attributes')
  x@ann$start <- as.numeric(as.character(x@ann$start))
  x@ann$end   <- as.numeric(as.character(x@ann$end))
  return(x)
}



##### ##### ##### ##### #####
#
# Getters.
#
##### ##### ##### ##### #####

#' @rdname create_chromR
#' @export
#' @aliases getPOS
getFIX <- function(x){
  if(class(x) != "chromR"){stop("expecting object of class chromR")}
  x@vcf@fix
}

#' @rdname create_chromR
#' @export
#' @aliases getPOS
getCHROM <- function(x){
  if(class(x) != "chromR"){stop("expecting object of class chromR")}
  x@vcf@fix[,"CHROM"]
}


#' @rdname create_chromR
#' @export
#' @aliases getPOS
getPOS <- function(x){
  if(class(x) != "chromR"){stop("expecting object of class chromR")}
  as.integer(x@vcf@fix[,"POS"])
}


#' @rdname create_chromR
#' @export
#' @aliases getPOS
getQUAL <- function(x){
  if(class(x) != "chromR"){stop("expecting object of class chromR")}
  as.integer(x@vcf@fix[,"QUAL"])
}



##### ##### ##### ##### #####
# EOF.