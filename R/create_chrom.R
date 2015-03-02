


#' @title Create Chrom object
#' @name Create Chrom object
#' @rdname create_chrom
#' @export
#' @aliases create_chrom
#'
#' @description
#' Creates and populates an object of class Chrom.
#'
#' @param title a name for the chromosome (for plotting purposes)
#' @param seq a sequence as a DNAbin object
#' @param vcf a vcfR object
#' @param ann an annotation file (gff-like)
#' @param verbose should verbose output be printed to the console?
#' @param x an object of class Chrom
#'
#' @details
#' Creates and names a chrom object from a name, a chromosome (an ape::DNAbin object), variant data (a vcfR object) and annotation data (gff-like).
#' The function \strong{create_chrom} is a wrapper which calls functions to populate the slots of the Chrom object.
#' 
#' The function \strong{seq2chrom} is currently defined as a generic function.
#' This may change in the future.
#' This function takes an object of class DNAbin and assigns it to the 'seq' slot of a Chrom object.
#' 
#' @seealso 
# \code{\link{seq2chrome}},
# \code{\link{vcf2chrome}},
#' \code{\link{Chrom-class}},
#' \code{\link{vcfR-class}},
#' \code{\link[ape]{DNAbin}},
#' \href{http://www.1000genomes.org/wiki/analysis/variant\%20call\%20format/vcf-variant-call-format-version-41}{vcf format}, 
#' \href{http://www.sequenceontology.org/gff3.shtml}{gff3 format}
#' 
#' @examples
#' library(vcfR)
#' data(vcfR_example)
#' pinf_mt <- create_chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)
#' head(pinf_mt)
#' pinf_mt
#' names(pinf_mt)
#' plot(pinf_mt)
#' pinf_mt <- masker(pinf_mt)
#' pinf_mt <- proc_chrom(pinf_mt, win.size=1000)
# pinf_mt <- proc.chrom(pinf_mt, win.size=1000)
#'  
# str(pinf_mt)
#' plot(pinf_mt)
#' 
#' chromoqc(pinf_mt)
#' chromoqc(pinf_mt, xlim=c(25e+03, 3e+04), dot.alpha=99)
#' 
#' set.seed(10)
#' x1 <- as.integer(runif(n=20, min=1, max=39000))
#' y1 <- runif(n=length(x1), min=1, max=100)
#' chromodot(pinf_mt, x1=x1, y1=y1)
#' chromodot(pinf_mt, x1=x1, y1=y1, label1='My data', x2=x1, y2=y1, label2='More of my data', dot.alpha='ff')
#' 
#' chromohwe(pinf_mt, dot.alpha='ff')
#' 
#' chromopop(pinf_mt)
#' gt <- extract.gt(pinf_mt)
#' head(gt)
#' tab <- variant_table(pinf_mt)
#' head(tab)
#' win <- window_table(pinf_mt)
#' head(win)
# hist(tab$Ho - tab$He, col=5)
# # Note that this example is a mitochondrion, so this is a bit silly.
#' 
create_chrom <- function(title="CHROM1", vcf, seq=NULL, ann=NULL, verbose=TRUE){
  # Determine whether we received the expected classes.
  stopifnot(class(seq)=="DNAbin")
  if(!is.null(vcf)){stopifnot(class(vcf) == "vcfR")}
  #
  x <- new(Class="Chrom")
  setName(x) <- title
  
  # Insert vcf into Chom.
  if(length(vcf)>0){
    x <- vcf2chrom(x, vcf)
  }
  # Insert seq into Chrom
  if(class(seq)=="DNAbin" | is.null(seq)){
#    seq2chrom(x) <- seq
    x <- seq2chrom(x, seq)
  } else {
    # stop???
    message("** Error: seq is not of class DNAbin** \n")
    break
  }
  #  if(length(ann)>0){
  if(nrow(ann)>0 | is.null(ann)){
    x <- ann2chrom(x, ann)
  }
  
  # Report names of objects to user.
  if(verbose == TRUE){
    # Print names of elements to see if they match.
    message("Names of sequences:")
    message(paste('  ', unique(names(x@seq)), sep=""))
    message("Names in vcf:")
    message(paste('  ', unique(as.character(x@vcf.fix$CHROM)), sep=""))
#    message(unique(as.character(x@vcf.fix$CHROM)))
    message("Names in annotation:")
    message(paste('  ', unique(as.character(x@ann[,1])), sep=""))
    if(unique(names(x@seq)) != unique(as.character(x@vcf.fix$CHROM)) | unique(names(x@seq)) != unique(as.character(x@ann[,1]))){
      message("Names in sequence file, variant file or annotation file do not match perfectly.\n")
      message("If you choose to proceed, we'll do our best to match data.\n")
      message("But prepare yourself for unexpected results.\n")
    }
  }
  return(x)
}








#### Data loading functions. ####



#' @rdname create_chrom
#' @export
#' @aliases chrom-methods vcf2chrom
#'
# @description
# Methods to work with objects of the chrom class
# Reads in a vcf file and stores it in a vcf class.
#'
# @param x an object of class chrom
#' @param y a vcfR object
#' @param ... arguments
#'
#' @details
#' The function \strong{vcf2chrom} is called by create_chrome and transfers the data from the slots of a vcfR object to the slots of a Chrom object.
#' It also tries to move data from the 'DP' and 'MQ' portion of the fix region of the vcf to the var.info slot of the Chrom object.
#' It is not anticipated that a user would need to use this function directly, but its placed here in case they do.
#'
vcf2chrom <- function(x,y,...){
  x@vcf.fix <- as.data.frame(y@fix)
  #  colnames(x@vcf.fix) <- c('chrom','pos','id','ref','alt','qual','filter','info')
  colnames(x@vcf.fix) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
  x@vcf.fix[,2] <- as.numeric(x@vcf.fix[,2])
  x@vcf.fix[,6] <- as.numeric(x@vcf.fix[,6])
  #
  for(i in 1:ncol(y@gt)){
    y@gt[,i] <- as.character(y@gt[,i])
  }
  x@vcf.gt <- y@gt
  #
  x@vcf.meta <- y@meta
  #
  info <- matrix(ncol=2, nrow=nrow(y@fix))
  #  colnames(info) <- c('dp','mq')
  colnames(info) <- c('DP','MQ')
  if(length(grep("DP=", y@fix[,8])) > 0){
    info[,1] <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(y@fix[,8]), ";"), function(x){grep("^DP=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
    #    x@var.info$DP <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(y@fix[,8]), ";"), function(x){grep("^DP=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
  }
  if(length(grep("MQ=", y@fix[,8])) > 0){
    info[,2] <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(y@fix[,8]), ";"), function(x){grep("^MQ=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
    #    x@var.info$MQ <- unlist(lapply(strsplit(unlist(lapply(strsplit(as.character(y@fix[,8]), ";"), function(x){grep("^MQ=", x, value=TRUE)})),"="),function(x){as.numeric(x[2])}))
  }
  #  x@vcf.info <- as.data.frame(info)
  x@var.info <- as.data.frame(info)
  #
  #  x@mask <- rep(TRUE, times=nrow(x@vcf.fix))
  x@var.info$mask <- rep(TRUE, times=nrow(x@vcf.fix))
  # assign may be more efficient.
  return(x)
}



#' @rdname create_chrom
#' @export
#' @aliases seq2chrom
#' 
seq2chrom <- function(x, seq=NULL){
  if(is.null(seq)){
    x@len <- x@vcf.fix$POS[nrow(x@vcf.fix)]
  }
  
  # A DNAbin will store in a list when the fasta contains
  # multiple sequences, but as a matrix when the fasta
  # only contains one sequence.
  if(!is.list(class(as.character(seq)))){
    x@seq <- as.list(seq)
    x@len <-length(x@seq[[1]])
  } else {
    x@seq <-seq
    x@len <-length(x@seq[[1]])
  }

  return(x)
}






#' @rdname create_chrom
#' @export
#' @aliases ann2chrom
#'
# @description
# Reads in an annotation file and stores it in a Chrom class.
#'
# @param x an object of class chrom
# @param y an object
# @param ... arguments
#'
#' @details
#' The function \strong{ann2chrom} is called by create_chrome and transfers the information from a gff-like object to the 'ann' slot of a Chrom object.
#' It is not anticipated that a user would need to use this function directly, but its placed here in case they do.
#'
#'
ann2chrom <- function(x,y,...){
  x@ann <- as.data.frame(y)
  colnames(x@ann) <- c('seqid','source','type','start','end','score','strand','phase','attributes')
  x@ann$start <- as.numeric(as.character(x@ann$start))
  x@ann$end   <- as.numeric(as.character(x@ann$end))
  return(x)
}

