
#' @title Chrom methods
#' 
#' @param x an object of class Chrom


#' @rdname Chrom-methods
#' @export
#' @aliases names.chrom
#' 
setMethod(
  f="names",
  signature = "Chrom",
  definition=function(x){
    #    cat("**** Class Chrom, method names **** \n")
    #    cat("Sequence name: ", as.character(names(x@seq)), "\n")
    #    cat("First annotation name: ")
    #    print(as.character(x@ann[1,1]))
    #    cat("First variant name: ")
    #    print(as.character(x@vcf.fix[1,1]))
    #    cat("\n")
    #    cat("Sample names: \n")
    temp <- names(x@vcf.gt)[-1]
    temp
  }
)


#' @rdname Chrom-methods
#' @export
#' @aliases head.chrom
#' 
setMethod(
  f="head",
  signature = "Chrom",
  definition=function(x){
    message("*** Class Chrom, method head *** \n")
    message(paste("Name: ", x@name, "\n"))
    message(paste("Length: ", x@len, "\n"))
    message("\n")
    message("**** ** Sample names (Chrom) ** **** \n")
    print(names(x@vcf.gt)[-1])
    message("\n")
    message("**** ** Vcf fixed data (Chrom) ** **** \n")
    print(x@vcf.fix[1:6,1:7])
    message("\nFirst INFO record:\n")
    print(unlist(strsplit(as.character(x@vcf.fix$INFO[1]), split=";")))
    message("\n")
    message("**** ** Vcf genotype data (Chrom) ** **** \n")
    if(ncol(x@vcf.gt)>=6){
      message("**** **** * First 6 columns * **** **** \n")
      print(x@vcf.gt[1:6,1:6])
    } else {
      print(x@vcf.gt[1:6,])
    }
    message("\n")
    message("**** ** Var info (Chrom) ** **** \n")
    if(ncol(x@var.info)>=6){
      message("**** **** First 6 columns ***** **** \n")
      print(x@var.info[1:6,1:6])
    } else {
      print(x@var.info[1:6,])
    }
    message("\n")
    message("**** ** Vcf mask (Chrom) ** **** \n")
    message("Percent unmasked: ")
    message(100*(sum(x@var.info$mask)/length(x@var.info$mask)))
    message("\n")
    message("**** ** End head (Chrom) ** **** \n")
  }
)


#' @rdname Chrom-methods
#' @export
#' @aliases plot.chrom
#' 
setMethod(
  f= "plot",
  signature= "Chrom",
  definition=function (x,y,...){
    par(mfrow=c(2,2))
    if(sum(!is.na(x@var.info$DP[x@var.info$mask])) >= 1){
      hist(x@var.info$DP[x@var.info$mask], col=3, main="Depth (DP)", xlab="")
      rug(x@var.info$DP[x@var.info$mask])
    } else {
      plot(1:2,1:2, type='n')
      title(main="No depths found")
    }
    if(sum(!is.na(x@var.info$MQ[x@var.info$mask])) >= 1){
      hist(x@var.info$MQ[x@var.info$mask], col=4, main="Mapping quality (MQ)", xlab="")
      rug(x@var.info$MQ[x@var.info$mask])
    } else {
      plot(1:2,1:2, type='n')
      title(main="No mapping qualities found")
    }
    if(sum(!is.na(x@vcf.fix$QUAL[x@var.info$mask])) >= 1){
      hist(x@vcf.fix$QUAL[x@var.info$mask], col=5, main="Quality (QUAL)", xlab="")
      rug(x@vcf.fix$QUAL[x@var.info$mask])
    } else {
      plot(1:2,1:2, type='n')
      title(main="No qualities found")
    }
    if(length(x@win.info$variants)>0){
      hist(x@win.info$variants/x@win.info$length, col=6, main="Variant count (per window)", xlab="")
      rug(x@win.info$variants/x@win.info$length)
    } else {
      plot(1:2,1:2, type='n')
      title(main="No SNP densities found")
    }
    par(mfrow=c(1,1))
  }
)


##### ##### Accessors.  #####

#### Getter for "names" ####
setGeneric("getName",function(object){standardGeneric ("getName")})

setMethod("getName","Chrom",
          function(object){
            return(object@name)
          }
)

#### Setter for name. ####

setGeneric("setName<-",function(object,value){standardGeneric("setName<-")})

setReplaceMethod(
  f="setName",
  signature="Chrom",
  definition=function(object,value){
    object@name <-value
    return (object)
  }
)

#### Setter for seq. ####

setGeneric("seq2chrom<-",function(object,value){standardGeneric("seq2chrom<-")})

setReplaceMethod(
  f="seq2chrom",
  signature="Chrom",
  definition=function(object,value){
    # A DNAbin will store in a list when the fasta contains
    # multiple sequences, but as a matrix when the fasta
    # only contains one sequence.
    if(!is.list(class(as.character(value)))){
      object@seq <- as.list(value)
    } else {
      object@seq <-value      
    }
    object@len <-length(object@seq[[1]])
    return (object)
  }
)

##### ##### ##### ##### #####
#### Data loading functions. ####

#' @title Chrom methods
#' @rdname Chrom-methods
#' @export
#' @aliases chrom-methods vcf2chrom
#'
#' @description
#' Methods to work with objects of the chrom class
# Reads in a vcf file and stores it in a vcf class.
#'
#' @param x an object of class chrom
#' @param y an object
#' @param ... arguments
#'
#' @details
#' Reads in a vcf file and stores it in a Chrom class.
#' 
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


# @title Chrom methods
#' @rdname Chrom-methods
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
#' Reads in a vcf file and stores it in a Chrom class.
#' 
#'
ann2chrom <- function(x,y,...){
  x@ann <- as.data.frame(y)
  colnames(x@ann) <- c('seqid','source','type','start','end','score','strand','phase','attributes')
  x@ann$start <- as.numeric(as.character(x@ann$start))
  x@ann$end   <- as.numeric(as.character(x@ann$end))
  return(x)
}

# @title Chrom methods
#' @rdname Chrom-methods
#' @export
#' @aliases create.chrom
#'
# @description
# Creates and populates a Chrom class.
#'
#' @param name a name for the object
#' @param seq a sequence as a DNAbin object
#' @param vcf a vcfR object
#' @param ann an annotation file (gff-like)
#'
#' @details
#' Creates and names a chrom object from a name, a sequence (a DNA bin object),
#' vcf data (a vcfR object) and annotation data (gff-like).
#' 
#' @examples
#' library(vcfR)
#' data(vcfR_example)
#' pinf_mt <- create.chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)
#' head(pinf_mt)
#' pinf_mt
#' names(pinf_mt)
#' plot(pinf_mt)
#' pinf_mt <- masker(pinf_mt)
#' pinf_mt <- proc.chrom(pinf_mt, win.size=1000)
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
create.chrom <- function(name, seq, vcf=NULL, ann=NULL, verbose=TRUE){
  stopifnot(class(seq)=="DNAbin")
  if(!is.null(vcf)){stopifnot(class(vcf) == "vcfR")}
  #
  x <- new(Class="Chrom")
  setName(x) <- name
  if(class(seq)=="DNAbin"){
    seq2chrom(x) <- seq
  } else {
    message("** Error: seq is not of class DNAbin** \n")
    break
  }
  if(length(vcf)>0){
    x <- vcf2chrom(x, vcf)
  }
  #  if(length(ann)>0){
  if(nrow(ann)>0){
    x <- ann2chrom(x, ann)
  }
  if(verbose == TRUE){
    # Print names of elements to see if they match.
    message("Names of sequences:\n")
    print(unique(names(x@seq)))
    message("Names in vcf:\n")
    print(unique(as.character(x@vcf.fix$CHROM)))
    message("Names in annotation:\n")
    print(unique(as.character(x@ann[,1])))
    if(unique(names(x@seq)) != unique(as.character(x@vcf.fix$CHROM)) | unique(names(x@seq)) != unique(as.character(x@ann[,1]))){
      message("Names in sequence file, variant file or annotation file do not match perfectly.\n")
      message("If you choose to proceed, we'll do our best to match data.\n")
      message("But prepare yourself for unexpected results.\n")
    }
  }
  return(x)
}


