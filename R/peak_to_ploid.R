
#' @title Convert allele balance peaks to ploidy
#' @name peak_to_ploid
#' @rdname peak_to_ploid
#' 
#' @description
#' Converts allele balance data produced by \code{freq_peak()} to a copy number by assinging the allele balance data (frequencies) to its closest expected ratio.
#'  
#' @param x an object produced by \code{freq_peak()}.
#' 
#' @details 
#' Converts allele balance data produced by \code{freq_peak()} to a copy number.
#' The expectation for five copies is that the allele balance data will be in fifths (1/5 or 4/5), four copies in quarters (1/4 or 3/4), three copies in thirds (1/3 or 2/3), and two copies in halves (1/2).
#' Based on this expectation we can assign a copy number to allele balance based on critical values that are half way between these expectation.
#' Allele balance values greater than or equal to 7/40 and less than 9/40 are called as five copies.
#' Allele balance values greater than or equal to 9/40 and less than 7/24 are called as four copies.
#' Allele balance values greater than or equal to 7/44 and less than 5/12 are called as three copies.
#' Allele balance values greater than or equal to 5/12 and less or equal to 7/12 are called as two copies.
#' Allele balance values greater than 7/12 and less than or equal to 17/24 are called as three copies.
#' Allele balance values greater than 17/24 and less than or equal to 31/40 are called as four copies.
#' Allele balance values greater than 31/40 and less than or equal to 33/40 are called as five copies.
#' Allele balance values greater than 33/40 or less than 7/40 are set as NA.
#' 
#' 
#' 
#' @examples
#' data(vcfR_example)
#' gt <- extract.gt(vcf)
#' # Censor non-heterozygous positions.
#' hets <- is_het(gt)
#' is.na(vcf@gt[,-1][!hets]) <- TRUE
#' # Extract allele depths.
#' ad <- extract.gt(vcf, element = "AD")
#' ad1 <- masplit(ad, record = 1)
#' ad2 <- masplit(ad, record = 2)
#' freq1 <- ad1/(ad1+ad2)
#' freq2 <- ad2/(ad1+ad2)
#' myPeaks1 <- freq_peak(freq1, getPOS(vcf))
#' # Censor windows with fewer than 20 heterozygous positions
#' is.na(myPeaks1$peaks[myPeaks1$counts < 20]) <- TRUE
#' # Convert peaks to ploidy call
#' peak_to_ploid(myPeaks1)
#' 
#' 
#' @export
peak_to_ploid <- function(x){
  
  # Validate our input
  if( class(x) != "list" | sum(names(x) == c("wins", "peaks", "counts")) != 3 ){
    msg <- "expecting a list with three elements named 'wins', 'peaks', and 'counts'"
    stop(msg)
  }
  
  # Initialize a result data structure.
#  gmat <- matrix(nrow=nrow(x$peaks), ncol=ncol(x$peaks))
#  colnames(gmat) <- colnames(x$peaks)
#  rownames(gmat) <- rownames(x$peaks)
  gmat <- x$peaks

  # Bin to ploidy
  #  critical <- 1/4 - (1/3-1/4)/2
  #  critical <- c(critical, 1/4 + (1/3-1/4)/2)
  #  critical <- c(critical, 1/2 - (2/3 - 1/2)/2)
  #  critical <- c(critical, 1/2 + (2/3 - 1/2)/2)
  #  critical <- c(critical, 3/4 - (1/3-1/4)/2)
  #  critical <- c(critical, 3/4 + (1/3-1/4)/2)

  critical <- c( 9/40, 7/24, 5/12, 7/12, 17/24, 31/40)
  
  gmat[ gmat <= 1 & gmat > critical[6] ] <- 5
  gmat[ gmat <= critical[6] & gmat > critical[5] ] <- 4
  gmat[ gmat <= critical[5] & gmat > critical[4] ] <- 3
  gmat[ gmat <= critical[4] & gmat >= critical[3] ] <- 2
  gmat[ gmat < critical[3] & gmat >= critical[2] ] <- 3
  gmat[ gmat < critical[2] & gmat >= critical[1] ] <- 4
  gmat[ gmat < critical[1] & gmat >= 0 ] <- 5
  
  is.na(gmat[x$peaks < 7/40 & !is.na(x$peaks)]) <- TRUE
  is.na(gmat[x$peaks > 33/40 & !is.na(x$peaks)]) <- TRUE

  return(gmat)
}

