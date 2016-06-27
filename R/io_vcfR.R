#' @title Read and write vcf format files
#' @rdname io_vcfR
#' @name VCF input and output
#' 
#' @export
#'
#' @description
#' Read and files in the *.vcf structured text format, as well as the compressed *.vcf.gz format.
#' Write objects of class vcfR to *.vcf.gz.
#' 
#' @param file A filename for a variant call format (vcf) file.
#' @param limit amount of memory (in bytes) not to exceed when reading in a file.
#' @param nrows integer specifying the maximum number of rows (variants) to read in.
#' @param skip integer specifying the number of rows (variants) to skip before beginning to read data.
#' @param cols vector of column numbers to extract from file.
#' @param x An object of class vcfR or chromR.
# @param vfile an output filename.
#' @param mask logical vector indicating rows to use.
#' @param APPEND logical indicating whether to append to existing vcf file or write a new file.
#' 
#' @param verbose report verbose progress.
#'
#' @details
#' The function \strong{read.vcfR} reads in files in *.vcf (text) and *.vcf.gz (gzipped text) format and returns an object of class vcfR.
#' The parameter 'limit' is an attempt to keep the user from trying to read in a file which contains more data than there is memory to hold.
#' Based on the dimensions of the data matrix, an estimate of how much memory needed is made.
#' If this estimate exceeds the value of 'limit' an error is thrown and execution stops.
#' The user may increase this limit to any value, but is encourages to compare that value to the amout of available physical memory.
#' 
#' 
#' It is possible to input part of a VCF file by using the parameters nrows, skip and cols.
#' The first eight columns (the fix region) are part of the definition and will always be included.
#' Any columns beyond eight are optional (the gt region).
#' You can specify which of these columns you would like to input by setting the cols parameter.
#' If you want a usable vcfR object you will want to always include nine (the FORMAT column).
#' If you do not include column nine you may experience reduced functionality.
#' 
#' 
#' The function \strong{write.vcf} takes an object of either class vcfR or chromR and writes the vcf data to a vcf.gz file (gzipped text).
#' If the parameter 'mask' is set to FALSE, the entire object is written to file.
#' If the parameter 'mask' is set to TRUE and the object is of class chromR (which has a mask slot), this mask is used to subset the data.
#' If an index is supplied as 'mask', then this index is used, and recycled as necessary, to subset the data.
#' 
#' Because vcfR provides the opportunity to manipulate VCF data, it also provides the opportunity for the user to create invalid VCF files.
#' If there is a question regarding the validity of a file you have created one option is the \href{https://vcftools.github.io/perl_module.html#vcf-validator}{VCF validator} from VCF tools.
#' 
#' 
#' @return read.vcfR returns an object of class \code{\link{vcfR-class}}.
#' See the \strong{vignette:} \code{vignette('vcf_data')}.
#' The function write.vcf creates a gzipped VCF file.
#'
#' @seealso
# \code{\link[PopGenome]{readVCF}}
# \code{\link[pegas]{read.vcf}}
# \link[pegas]{read.vcf}
#'
#' CRAN:
#' \href{http://cran.r-project.org/package=pegas}{pegas}::read.vcf,
#' \href{http://cran.r-project.org/package=PopGenome}{PopGenome}::readVCF,
#' \href{http://cran.r-project.org/package=data.table}{data.table}::fread
#' 
#' Bioconductor:
#' \href{http://www.bioconductor.org/packages/release/bioc/html/VariantAnnotation.html}{VariantAnnotation}::readVcf
#'
#' Use: browseVignettes('vcfR') to find examples.
#'
#'
#' @examples 
#' data(vcfR_test)
#' vcfR_test
#' head(vcfR_test)
#' # CRAN requires developers to us a tempdir when writing to the filesystem.
#' # You may want to implement this example elsewhere.
#' orig_dir <- getwd()
#' temp_dir <- tempdir()
#' setwd( temp_dir )
#' write.vcf( vcfR_test, file = "vcfR_test.vcf.gz" )
#' vcf <- read.vcfR( file = "vcfR_test.vcf.gz", verbose = FALSE )
#' vcf
#' setwd( orig_dir )
#' 
#' 
#' @rdname io_vcfR
#' @aliases read.vcfR
#' @export
#' 
read.vcfR <- function(file, limit=1e7, nrows = -1, skip = 0, cols = NULL, verbose = TRUE){
#  require(memuse)
  
  # gzopen does not appear to deal well with tilde expansion.
  file <- path.expand(file)
  
  if(file.access(file, mode = 0) != 0){
    stop(paste("File:", file, "does not appear to exist!"))
  }
  if(file.access(file, mode = 4) != 0){
    stop(paste("File:", file, "appears to exist but is not readable!"))
  }
  
  vcf <- new(Class="vcfR")
  stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', file)
  
  if( is.null(cols) ){
    cols <- 1:stats['columns']
  }
  # Make sure we include the first nine columns.
  cols <- sort( unique( c(1:8, cols) ) )

  
#  ram_est <- stats['variants'] * stats['columns'] * 8 + 248
  ram_est <- memuse::howbig(stats['variants'], stats['columns'])
  
  if(ram_est@size > limit){
    message(paste("The number of variants in your file is:", prettyNum(stats['variants'], big.mark=",")))
    message(paste("The number of samples in your file is:", prettyNum(stats['columns'] - 1, big.mark=",")))
    message(paste("This will result in an object of approximately:", ram_est, "in size"))
    stop("Object size limit exceeded")
  }
  
  vcf@meta <- .Call('vcfR_read_meta_gz', PACKAGE = 'vcfR', file, stats, as.numeric(verbose))
  body <- .Call('vcfR_read_body_gz', PACKAGE = 'vcfR', file = file, stats = stats, 
                nrows = nrows, skip = skip, cols = cols, as.numeric(verbose))

  vcf@fix <- body[ ,1:8, drop=FALSE ]
  if( ncol(body) > 8 ){
    vcf@gt <- body[ , -c(1:8), drop=FALSE ]
  } else {
    vcf@gt <- matrix("a", nrow=0, ncol=0)
  }
  
  return(vcf)
}




#' @rdname io_vcfR
#' @export
#' @aliases write.vcf
#' 
write.vcf <- function(x, file = "", mask = FALSE, APPEND = FALSE){
  if(class(x) == "chromR"){
    if( mask == TRUE ){
      is.na( x@vcf@fix[,'FILTER'] ) <- TRUE
      x@vcf@fix[,'FILTER'][ x@var.info[,'mask'] ] <- 'PASS'
    }
    x <- x@vcf
  }
  if(class(x) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or chromR.")
  }
  
  if(APPEND == FALSE){
    gz <- gzfile(file, "w")
    
    if( length(x@meta) > 0 ){
      write(x@meta, gz)
    }
    
    header <- c(colnames(x@fix), colnames(x@gt))
    header[1] <- "#CHROM"
    header <- paste(header, collapse="\t")
    write(header, gz)
    close(gz)
  }
  
  if(mask == FALSE){
    test <- .Call('vcfR_write_vcf_body', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 0)
  } else if (mask == TRUE){
    test <- .Call('vcfR_write_vcf_body', PACKAGE = 'vcfR', fix = x@fix, gt = x@gt, filename = file, mask = 1)
  }
}


##### ##### ##### ##### #####
# EOF