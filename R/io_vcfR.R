

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
#' @param limit amount of memory (in bytes) not to exceed when reading in a file
#' @param verbose report verbose progress
#'
#' @details
#' The function \strong{read.vcf} reads in files in *.vcf (text) and *.vcf.gz (gzipped text) format and returns an object of class vcfR.
#' The parameter 'limit' is an attempt to keep the user from trying to read in a file which contains more data than there is memory to hold.
#' Based on the dimensions of the data matrix, an estimate of how much memory needed is made.
#' If this estimate exceeds the value of 'limit' an error is thrown and execution stops.
#' 
#' 
#' The function \strong{write.vcf} takes an object of either class vcfR or chromR and writes the vcf data to a vcf.gz file (gzipped text).
#' If the parameter 'mask' is set to FALSE, the entire object is written to file.
#' If the parameter 'mask' is set to TRUE and the object is of class chromR (which has a mask slot), this maske is used to subset the data.
#' If an index is supplied as 'mask', then this index is used, and recycled as necessary, to subset the data.
#' 
#' The function \strong{write_var_info} takes the variant information table from a chromR object and writes it as a comma delimited file. 
#' 
#' The function \strong{write_win_info} takes the window information table from a chromR object and writes it as a comma delimited file.
#' 
#' The function \strong{write_fasta} takes an object of class chromR and writes it to a fasta.gz (gzipped text) format file.
#' The sequence in the seq slot of the chromR object is used to fill in the invariant sites.
#' The parameter 'tolower', when set to TRUE, converts all the characters in teh sequence to lower case.
#' This is important because some software, such as ape::DNAbin, requires sequences to be in lower case.
#' 
#' 
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
# data(vcfR_example)
# head(pinf_vcf)
# plot(pinf_vcf)
# pinf_vcf[1:6,]
# 
#' @rdname io_vcfR
#' @aliases read.vcf
#' @export
#' 
read.vcf <- function(file, limit=1e7, verbose = TRUE){
#  require(memuse)
  
  if(file.access(file, mode = 0) != 0){
    stop(paste("File:", file, "does not appear to exist!"))
  }
  if(file.access(file, mode = 4) != 0){
    stop(paste("File:", file, "appears to exist but is not readable!"))
  }
  
  vcf <- new(Class="vcfR")
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', file)
  
#  ram_est <- stats['variants'] * stats['columns'] * 8 + 248
  ram_est <- memuse::howbig(stats['variants'], stats['columns'])
  
  if(ram_est@size > limit){
    message(paste("The number of variants in your file is:", prettyNum(stats['variants'], big.mark=",")))
    message(paste("The number of samples in your file is:", prettyNum(stats['columns'] - 1, big.mark=",")))
    message(paste("This will result in an object of approximately:", ram_est, "in size"))
    stop("Object size limit exceeded")
  }
  
  vcf@meta <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', file, stats, as.numeric(verbose))
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', file, stats, as.numeric(verbose))

  vcf@fix <- body[,1:8]
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
#    x <- chrom_to_vcfR(x)
    x <- x@vcf
#    x@fix$FILTER[filter] <- "PASS"
    x@fix[,'FILTER'] <- "PASS"
  }
  if(class(x) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or Chrom.")
  }
  
  if(APPEND == FALSE){
    gz <- gzfile(file, "w")
    write(x@meta, gz)

#    header <- c(names(x@fix), names(x@gt))
#    header[1] <- "#CHROM"
#    header <- paste(header, collapse="\t")
#    write(header, gz)

    close(gz)
  }
  
  if(mask == FALSE){
    test <- .Call('vcfR_write_vcf_body', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 0)
#    test <- .Call('vcfR_write_vcf_body_gz', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 0)
  } else if (mask == TRUE){
    test <- .Call('vcfR_write_vcf_body', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 1)
#    test <- .Call('vcfR_write_vcf_body_gz', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 1)
  }
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
#' @param gt_split character which delimits alleles in genotype
#' @param rowlength number of characters each row should not exceed
#' @param tolower convert all characters to lowercase (T/F)
#' 
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
#    seq[x@vcf.fix$POS] <- haps[,i]
    seq[x@var.info$POS] <- haps[,i]
    invisible(.Call('vcfR_write_fasta', PACKAGE = 'vcfR', seq, colnames(haps)[i], file, rowlength, as.integer(verbose)))
  }
  
  #invisible(.Call('vcfR_write_fasta', PACKAGE = 'vcfR', seq, seqname, filename, rowlength, verbose))  
}
