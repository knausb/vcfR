

#' @title Read and write vcf format files
#' @rdname io_vcfR
#' @export
#'
#' @description
#' Read and write files in the vcf format.
#' 
#' @param file A filename for a variant call format (vcf) file
#' @param x An object of class vcfR or Chrom
# @param vfile an output filename
#' @param mask logical vector indicating rows to use
#' @param APPEND logical indicating whether to append to existing vcf file or write a new file
#' @param limit amount of memory not to exceed when reading in a file
#' @param verbose report verbose progress
#'
#' @details
#' Reads in a vcf file and stores it in a vcf class.  Once the number of lines the meta information contains the data is divided into three tables: meta data, fixed data and genotype data.
#'
# \code{read.vcf.devel} is a custom C++ implementation which may be faster than read.vcf.
#'
#' @seealso
# \code{\link[PopGenome]{readVCF}}
# \code{\link[pegas]{read.vcf}}
# \link[pegas]{read.vcf}
#'
#' CRAN:
#' \href{http://cran.r-project.org/web/packages/pegas/index.html}{pegas}::read.vcf,
#' \href{http://cran.r-project.org/web/packages/PopGenome/index.html}{PopGenome}::readVCF,
#' \href{http://cran.r-project.org/web/packages/data.table/index.html}{data.table}::fread
#' 
#' Bioconductor:
#' \href{http://www.bioconductor.org/packages/release/bioc/html/VariantAnnotation.html}{VariantAnnotation}::readVcf
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
read.vcf.R<-function(file){
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
#' @aliases write.vcf.R
#' 
# @usage write.vcf(xvcf, vfile)
#' 
#' 
#' @export
#' 
#write.vcf<-function(xvcf, vfile, mask=logical(0), APPEND=FALSE){
write.vcf.R<-function(x, file="", mask=logical(0), APPEND=FALSE){
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
#' @aliases read.vcf
#' @export
#' 
read.vcf <- function(file, limit=1e7, verbose = TRUE){
    
  if(file.access(file, mode = 0) != 0){
    stop(paste("File:", file, "does not appear to exist!"))
  }
  if(file.access(file, mode = 4) != 0){
    stop(paste("File:", file, "appears to exist but is not readable!"))
  }
  
  
  vcf <- new(Class="vcfR")
#  stats <- .Call('vcfR_vcf_stats', PACKAGE = 'vcfR', file)
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', file)
  
#  element_size <- object.size(.Call('vcfR_ram_test', PACKAGE = 'vcfR'))
#  element_size <- object.size(ram_test())
  
#  ram_est <- stats['variants'] * stats['columns'] * element_size
  ram_est <- stats['variants'] * stats['columns'] * 8 + 248

  if(ram_est > limit){
    message(paste("The number of variants in your file is:", prettyNum(stats['variants'], big.mark=",")))
    message(paste("The number of samples in your file is:", prettyNum(stats['columns'] - 1, big.mark=",")))
    message(paste("This will result in an object of approximately:", round(ram_est/1e9, digits=3), "Gb in size"))
    message("The function memory_plot() may help you estimate your memory needs.")
    stop("Object size limit exceeded")
  }
  
#  vcf@meta <- .Call('vcfR_vcf_meta', PACKAGE = 'vcfR', file, stats)
#  body <- .Call('vcfR_vcf_body', PACKAGE = 'vcfR', file, stats)

  vcf@meta <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', file, stats, as.numeric(verbose))
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', file, stats, as.numeric(verbose))

  for(i in c(1, 3:5, 7:8)){
    body[,i] <- as.character(body[,i])
  }
  for(i in 9:ncol(body)){
    body[,i] <- as.character(body[,i])
  }
  
  vcf@fix <- body[,1:8]
  vcf@fix$POS <- as.integer(as.character(vcf@fix$POS))
  vcf@fix$QUAL <- as.numeric(as.character(vcf@fix$QUAL))
  vcf@gt <- body[,9:ncol(body)]
  
  return(vcf)
}


#' @rdname io_vcfR
#' @export
#' @aliases write.vcf
#' 
write.vcf <- function(x, file = "", mask = FALSE, APPEND = FALSE){
  if(class(x) == "Chrom"){
    filter <- x@var.info$mask
    x <- chrom_to_vcfR(x)
    x@fix$FILTER[filter] <- "PASS"
  }
  if(class(x) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or Chrom.")
  }
  
  if(APPEND == FALSE){
    write(x@meta, file = file, append = FALSE)
    header <- c(names(x@fix), names(x@gt))
    header[1] <- "#CHROM"
    header <- paste(header, collapse="\t")
    write(header, file = file, append = TRUE)
  }
  
  if(mask == FALSE){
    test <- .Call('vcfR_write_vcf_body', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 0)
  } else if (mask == TRUE){
    test <- .Call('vcfR_write_vcf_body', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 1)
  }
}



#' @rdname io_vcfR
#' @export
#' @aliases write.vcf.gz
#' 
write.vcf.gz <- function(x, file = "", mask = FALSE, APPEND = FALSE){
  if(class(x) == "Chrom"){
    filter <- x@var.info$mask
    x <- chrom_to_vcfR(x)
    x@fix$FILTER[filter] <- "PASS"
  }
  if(class(x) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or Chrom.")
  }
  
  if(APPEND == FALSE){
    gz <- gzfile(file, "w")
    write(x@meta, gz)

    header <- c(names(x@fix), names(x@gt))
    header[1] <- "#CHROM"
    header <- paste(header, collapse="\t")
    write(header, gz)

    close(gz)
  }
  
  if(mask == FALSE){
    test <- .Call('vcfR_write_vcf_body_gz', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 0)
  } else if (mask == TRUE){
    test <- .Call('vcfR_write_vcf_body_gz', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 1)
  }
}






#' @rdname io_vcfR
#' @aliases memory_plot
#' 
#' @param exponent_range range of values to be used for object sizes (10^exponent_range)
#' 
#' @export
#' 
memory_plot <- function(exponent_range=2:6){
  msize <- 10^(exponent_range)
  osize <- vector(length=length(exponent_range))

  for(i in 1:length(exponent_range)){
#    osize[i] <- object.size(ram_test(nrow = msize[i], ncol = 1L))
    osize[i] <- object.size(.Call('vcfR_ram_test', PACKAGE = 'vcfR', nrow = msize[i], ncol = 1))
  }

#  plot(msize, osize, type='b', log='xy')
  plot(log10(msize), 
       log10(as.numeric(osize)/1e6), 
       type='b',
       #log='xy',
       xlab="log10(Matrix cells)",
       ylab="log10(object size(Mb))")

  x <- cbind(msize, osize)
  colnames(x) <- c("Matrix cells", "Memory (bytes)")
  return(x)
}



#' @rdname io_vcfR
#' @aliases write_var_info
#' 
#' @export
#' 
write_var_info <- function(x, file = "", mask = FALSE, APPEND = FALSE){
  if(class(x) == "vcfR"){
    stop("Unexpected class! Detected class vcfR. This class does not contain variant summaries.")
  }
  if(class(x) != "Chrom"){
    stop("Unexpected class! Expecting an object of class Chrom.")
  }
  
  if(mask == FALSE){
    write.table(x@var.info, file = file, append = APPEND, sep = ",", row.names = FALSE, col.names = !APPEND)
  } else if(mask == TRUE){
    write.table(x@var.info[x@var.info$mask,], file = file, append = APPEND, sep = ",", row.names = FALSE, col.names = !APPEND)
  }
}



#' @rdname io_vcfR
#' @aliases write_win_info
#' 
#' @export
#' 
write_win_info <- function(x, file = "", APPEND = FALSE){
  if(class(x) == "vcfR"){
    stop("Unexpected class! Detected class vcfR. This class does not contain window summaries.")
  }
  if(class(x) != "Chrom"){
    stop("Unexpected class! Expecting an object of class Chrom.")
  }
  
  write.table(x@win.info, file = file, append = APPEND, sep = ",", row.names = FALSE, col.names = !APPEND)
}


#' @rdname io_vcfR
#' @aliases write_fasta
#' 
#' @export
#' 
write_fasta <- function(x, file = "", gt_split = "|", rowlength=80, tolower=TRUE, verbose=TRUE, APPEND = FALSE){
  if(class(x) != "Chrom"){
    stop("Expected object of class Chrom")
  }
  if(APPEND == FALSE){
    if(file.exists(file)){
      file.remove(file)
    }
  }
  haps <- extract_haps(x, gt_split = gt_split)
  if(tolower == TRUE){
    haps <- apply(haps, MARGIN=2, tolower)
  }
  
  for(i in 1:ncol(haps)){
    seq <- as.character(x@seq)[1,]
    seq[x@vcf.fix$POS] <- haps[,i]
    invisible(.Call('vcfR_write_fasta', PACKAGE = 'vcfR', seq, colnames(haps)[i], file, rowlength, as.integer(verbose)))
  }
  
  #invisible(.Call('vcfR_write_fasta', PACKAGE = 'vcfR', seq, seqname, filename, rowlength, verbose))  
}
