

#' @title Read and write vcf format files
#' @rdname io_vcfR
#' @export
#'
#' @description
#' Read and write files in the vcf format.
#' 
#' @param file A filename for a variant call format (vcf) file
#' @param x A vcfR object
# @param vfile an output filename
#' @param mask logical vector indicating rows to use
#' @param APPEND logical indicating whether to append to existing vcf file or write a new file
#'
#' @details
#' Reads in a vcf file and stores it in a vcf class.  Once the number of lines the meta information contains the data is divided into three tables: meta data, fixed data and genotype data.
#'
#' \code{read.vcf.devel} is a custom C++ implementation which may be faster than read.vcf.
#'
#' @seealso
#' \code{\link[PopGenome]{readVCF}}
# \code{\link[pegas]{read.vcf}}
#' \link[pegas]{read.vcf}
#' 
#' Bioconductor: \href{www.bioconductor.org/packages/3.0/bioc/html/VariantAnnotation.html}{VariantAnnotation::readVcf}
#'
#'
#'
#'
#' @examples
#' library(vcfR)
#' data(vcfR_example)
#' head(pinf_vcf)
#' plot(pinf_vcf)
#' pinf_vcf[1:6,]
#' 
read.vcf<-function(file){
  vcf <- new(Class="vcfR")
  i <- -1 # Line counter.
  j <- 0 # Success?
  tmp <- scan(file, what="character", sep="\n", skip=0, nlines=1, quiet=T, comment.char="")
  if(length(grep('^##', tmp)) >0 ){
    while(j == 0){
      i <- i+1
      tmp <- scan(file, what="character", sep="\n", skip=i, nlines=1, quiet=T, comment.char="")
      if(length(grep('^##',tmp)) == 0){j <- 1}
    }
    vcf@meta <- scan(file, what = "character", sep = "\n", skip = 0, nlines = i,
                     quiet = T, comment.char = "")
    vcf@fix <- read.table(file, header = TRUE, sep = '\t', skip = i, comment.char = '',
                          colClasses = "character", check.names = FALSE)
    vcf@gt <- vcf@fix[,9:ncol(vcf@fix)]
    vcf@fix <- vcf@fix[,1:8]
    vcf@fix$POS  <- as.integer(vcf@fix$POS)
    vcf@fix$QUAL <- as.integer(vcf@fix$QUAL)
  } else if (length(grep('^#',tmp)) >0 ){
    # No meta region, but a header line.
    vcf@fix <- read.table(file, header=T, sep='\t', skip=0,comment.char='', colClasses = "character")
    vcf@gt <- vcf@fix[,9:ncol(vcf@fix)]
    vcf@fix <- vcf@fix[,1:8]
    #    colnames(vcf@fix) <- c('chrom','pos','id','ref','alt','qual','filter','info')
    colnames(vcf@fix) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
    vcf@fix$POS  <- as.integer(vcf@fix$POS)
    vcf@fix$QUAL <- as.integer(vcf@fix$QUAL)
  } else {
    # No meta region or header line.
    vcf@fix <- read.table(file,header=F,sep='\t',skip=0,comment.char='', colClasses = "character")
    vcf@gt <- vcf@fix[,9:ncol(vcf@fix)]
    vcf@fix <- vcf@fix[,1:8]
    colnames(vcf@fix) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
    vcf@fix$POS  <- as.integer(vcf@fix$POS)
    vcf@fix$QUAL <- as.integer(vcf@fix$QUAL)    
  }
  if(names(vcf@fix)[1] != "CHROM"){
    names(vcf@fix)[1] <- "CHROM"
  }
  return(vcf)
}




#' @rdname io_vcfR
#' @aliases write.vcf
#' 
# @usage write.vcf(xvcf, vfile)
#' 
#' 
#' @export
#' 
#write.vcf<-function(xvcf, vfile, mask=logical(0), APPEND=FALSE){
write.vcf<-function(x, file="", mask=logical(0), APPEND=FALSE){
  if(class(x) == 'Chrom'){
    # Recast as a vcfR object.
    temp <- x
    x <- new(Class="vcfR")
    x@meta <- temp@vcf.meta
    x@fix <- temp@vcf.fix
    x@gt <- temp@vcf.gt
    mask <- temp@var.info$mask
    rm(temp)
  }
  if(class(x) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or Chrom.")
  }
  #
  if(length(mask) == 0){
    #    mask <- 1:nrow(xvcf@fix)
    mask <- rep(TRUE, times=nrow(x@fix))
  }
  #
  orig_scipen <- getOption("scipen")
  options(scipen=999)
  if(APPEND == FALSE){
    header <- c(names(x@fix), names(x@gt))
    header[1] <- paste("#",header[1],sep='')
    write.table(x@meta, file = file, append = FALSE, quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".", row.names = FALSE,
                col.names = FALSE)
    write(header, file = file,
          ncolumns=length(header),
          append = TRUE,
          sep = "\t")
    write.table(cbind(x@fix[mask,], x@gt[mask,]), file = file, append = TRUE,
                quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".",
                row.names = FALSE,
                col.names = FALSE)
    #              col.names = TRUE)
  } else if (APPEND == TRUE){
    write.table(cbind(x@fix[mask,], x@gt[mask,]), file = file, append = TRUE,
                quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".",
                row.names = FALSE,
                col.names = FALSE)    
  }
  options(scipen=orig_scipen)
}


#' @rdname io_vcfR
#' @aliases read.vcf.devel3
#' @export
#' 
read.vcf.devel <- function(file){
  vcf <- new(Class="vcfR")
  vcf@meta <- .Call('vcfR_readVcfHeader', PACKAGE = 'vcfR', file)
#  temp <- .Call('vcfR_readVcfBody', PACKAGE = 'vcfR', file)
#  vcf@fix <- as.data.frame(temp[-1,1:8])
#  names(vcf@fix) <- temp[1,1:8]
#  vcf@gt <- as.data.frame(temp[-1,-c(1:8)])
#  names(vcf@gt) <- temp[1,-c(1:8)]
  temp <- .Call('vcfR_readVcfBody2', PACKAGE = 'vcfR', file)
  vcf@fix <- temp[,1:8]
  vcf@gt <- temp[,9:ncol(temp)]
  return(vcf)
}


#' @rdname io_vcfR
#' @export
#' 
read.vcf.devel3 <- function(file, limit=1e9){
  vcf <- new(Class="vcfR")
  stats <- .Call('vcfR_vcf_stats', PACKAGE = 'vcfR', file)
#  element_size <- object.size(.Call('vcfR_ram_test', PACKAGE = 'vcfR'))
  element_size <- object.size(ram_test())
  
  ram_est <- stats['variants'] * stats['columns'] * element_size
  if(ram_est > limit){
    message(paste("The number of variants in your file is:", prettyNum(stats['variants'], big.mark=",")))
    message(paste("The number of samples in your file is:", prettyNum(stats['columns'] - 1, big.mark=",")))
    message(paste("This will result in an object of approximately:", round(ram_est/1e9, digits=3), "Gb in size"))
    stop("Object size limit exceeded")
  }
  
  vcf@meta <- .Call('vcfR_vcf_meta', PACKAGE = 'vcfR', file, stats)
  temp <- .Call('vcfR_vcf_body', PACKAGE = 'vcfR', file, stats)

  for(i in c(1, 3:5, 7:8)){
    temp[,i] <- as.character(temp[,i])
  }
  for(i in 9:ncol(temp)){
    temp[,i] <- as.character(temp[,i])
  }
  
  vcf@fix <- temp[,1:8]
#  vcf@fix$POS <- as.integer(as.character(vcf@fix$POS))
  vcf@gt <- temp[,9:ncol(temp)]
  
  return(vcf)
}

