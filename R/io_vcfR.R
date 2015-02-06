
#' Read vcf format file
#'
#' @rdname io_vcfR
#' 
read.vcf.devel <- function(x){
  vcf <- new(Class="vcfR")
  return(0)
}


#' @rdname io_vcfR
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
    vcf@meta <- scan(x, what = "character", sep = "\n", skip = 0, nlines = i,
                     quiet = T, comment.char = "")
    vcf@fix <- read.table(x, header = TRUE, sep = '\t', skip = i, comment.char = '',
                          colClasses = "character", check.names = FALSE)
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
    colnames(vcf@fix) <- c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO')
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




#' @rdname io_vcfR
#' @aliases write.vcf
#' 
# @usage write.vcf(xvcf, vfile)
#' 
#' @param xvcf a vcfR object
#' @param vfile an output filename
#' @param mask logical vector indicating rows to use
#' @param APPEND logical indicating whether to append to existing vcf file or write a new file
#' 
#' @export
#' 
#write.vcf<-function(xvcf, vfile, mask=logical(0)){
write.vcf<-function(xvcf, vfile, mask=logical(0), APPEND=FALSE){
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
  if(APPEND == FALSE){
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
  } else if (APPEND == TRUE){
    write.table(cbind(xvcf@fix[mask,], xvcf@gt[mask,]), file = vfile, append = TRUE,
                quote = FALSE, sep = "\t",
                eol = "\n", na = "NA", dec = ".",
                row.names = FALSE,
                col.names = FALSE)    
  }
  options(scipen=orig_scipen)
}



