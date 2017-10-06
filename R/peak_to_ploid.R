
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
#' @return A list consisting of two matrices containing the calls and the distance from expectation (i.e., confidence).
#' 
#' 
#' @seealso
#' freq_peak,
#' freq_peak_plot
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
#  if( class(x) != "list" | sum(names(x) == c("wins", "peaks", "counts")) != 3 ){
#    msg <- "expecting a list with three elements named 'wins', 'peaks', and 'counts'"
#    stop(msg)
#  }
  if( !inherits(x, "freq_peak") ){
    msg <- "expecting a freq_peak object."
  }
  
  
  # Initialize a result data structure.
#  gmat <- matrix(nrow=nrow(x$peaks), ncol=ncol(x$peaks))
#  colnames(gmat) <- colnames(x$peaks)
#  rownames(gmat) <- rownames(x$peaks)
  gmat <- x$peaks
  # Allele balance expectation
  abe <- matrix(ncol=ncol(gmat), nrow = nrow(gmat))

  # Bin to ploidy
  #  critical <- 1/4 - (1/3-1/4)/2
  #  critical <- c(critical, 1/4 + (1/3-1/4)/2)
  #  critical <- c(critical, 1/2 - (2/3 - 1/2)/2)
  #  critical <- c(critical, 1/2 + (2/3 - 1/2)/2)
  #  critical <- c(critical, 3/4 - (1/3-1/4)/2)
  #  critical <- c(critical, 3/4 + (1/3-1/4)/2)

  critical <- c( 9/40, 7/24, 5/12, 7/12, 17/24, 31/40)
  
  abe[ gmat <= 1 & gmat > critical[6] ] <- 4/5
  gmat[ gmat <= 1 & gmat > critical[6] ] <- 5
  abe[ gmat <= critical[6] & gmat > critical[5] ] <- 3/4
  gmat[ gmat <= critical[6] & gmat > critical[5] ] <- 4
  abe[ gmat <= critical[5] & gmat > critical[4] ] <- 2/3
  gmat[ gmat <= critical[5] & gmat > critical[4] ] <- 3
  abe[ gmat <= critical[4] & gmat >= critical[3] ] <- 1/2
  gmat[ gmat <= critical[4] & gmat >= critical[3] ] <- 2
  abe[ gmat < critical[3] & gmat >= critical[2] ] <- 1/3
  gmat[ gmat < critical[3] & gmat >= critical[2] ] <- 3
  abe[ gmat < critical[2] & gmat >= critical[1] ] <- 1/4
  gmat[ gmat < critical[2] & gmat >= critical[1] ] <- 4
  abe[ gmat < critical[1] & gmat >= 0 ] <- 1/5
  gmat[ gmat < critical[1] & gmat >= 0 ] <- 5
  
  is.na(gmat[x$peaks < 7/40  & !is.na(x$peaks)]) <- TRUE
  is.na(gmat[x$peaks > 33/40 & !is.na(x$peaks)]) <- TRUE
  is.na(abe[x$peaks  < 7/40  & !is.na(x$peaks)]) <- TRUE
  is.na(abe[x$peaks  > 33/40 & !is.na(x$peaks)]) <- TRUE

  # Distance from expectation
  dfe <- x$peaks - abe
  
  # Scale dfe by bin width
  dfe[ abe == 4/5 & !is.na(dfe) ] <- dfe[ abe == 4/5 & !is.na(abe) ] / (33/40 - 4/5)
  dfe[ abe == 3/4 & !is.na(dfe) & dfe > 0 ] <- dfe[ abe == 3/4 & !is.na(dfe) & dfe > 0 ] / (31/40 - 3/4)
  dfe[ abe == 3/4 & !is.na(dfe) & dfe < 0 ] <- dfe[ abe == 3/4 & !is.na(dfe) & dfe < 0 ] / (3/4 - 17/24)
  dfe[ abe == 2/3 & !is.na(dfe) & dfe > 0 ] <- dfe[ abe == 2/3 & !is.na(dfe) & dfe > 0 ] / (17/24 - 2/3)
  dfe[ abe == 2/3 & !is.na(dfe) & dfe < 0 ] <- dfe[ abe == 2/3 & !is.na(dfe) & dfe < 0 ] / (2/3 - 7/12)
  dfe[ abe == 1/2 & !is.na(dfe) ] <- dfe[ abe == 1/2 & !is.na(dfe) ] / (7/12 - 1/2)
  dfe[ abe == 1/3 & !is.na(dfe) & dfe > 0 ] <- dfe[ abe == 1/3 & !is.na(dfe) & dfe > 0 ] / (5/12 - 1/3)
  dfe[ abe == 1/3 & !is.na(dfe) & dfe < 0 ] <- dfe[ abe == 1/3 & !is.na(dfe) & dfe < 0 ] / (1/3 - 7/24)
  dfe[ abe == 1/4 & !is.na(dfe) & dfe > 0 ] <- dfe[ abe == 1/4 & !is.na(dfe) & dfe > 0 ] / (7/24 - 1/4)
  dfe[ abe == 1/4 & !is.na(dfe) & dfe < 0 ] <- dfe[ abe == 1/4 & !is.na(dfe) & dfe < 0 ] / (1/4 - 9/40)
  dfe[ abe == 1/5 & !is.na(dfe) ] <- dfe[ abe == 1/5 & !is.na(dfe) ] / (9/40 - 1/5)
  
  #return(gmat)
  list( calls = gmat, 
        #abe = abe,
        dfe = dfe)
}

