
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
#' Converts allele balance data produced by \code{freq_peak()} to copy number.
#' See the examples section for a graphical representation of the expectations and the bins around them.
#' Once a copy number has called a distance from expectation (dfe) is calculated as a form of confidence.
#' The bins around different copy numbers are of different width, so the dfe is scaled by its respective bin width.
#' This results in a dfe that is 0 when it is exactly at our expectation (high confidence) and at 1 when it is half way between two expectations (low confidence).
#' 
#' 
#' @seealso \code{freq_peak}, \code{freq_peak_plot}
#' 
#' 
#' @return A list consisting of two matrices containing the calls and the distance from expectation (i.e., confidence).
#' 
#' 
#' @examples
#' # Thresholds.
#' plot(c(0.0, 1), c(0,1), type = "n", xaxt = "n", xlab = "Expectation", ylab = "Allele balance")
#' myCalls <-  c(1/5, 1/4, 1/3, 1/2, 2/3, 3/4, 4/5)
#' axis(side = 1, at = myCalls, labels = c('1/5', '1/4', '1/3','1/2', '2/3', '3/4', '4/5'), las=2)
#' abline(v=myCalls)
#' abline(v=c(7/40, 9/40, 7/24, 5/12), lty=3, col ="#B22222")
#' abline(v=c(7/12, 17/24, 31/40, 33/40), lty=3, col ="#B22222")
#' text(x=7/40, y=0.1, labels = "7/40", srt = 90)
#' text(x=9/40, y=0.1, labels = "9/40", srt = 90)
#' text(x=7/24, y=0.1, labels = "7/24", srt = 90)
#' text(x=5/12, y=0.1, labels = "5/12", srt = 90)
#' text(x=7/12, y=0.1, labels = "7/12", srt = 90)
#' text(x=17/24, y=0.1, labels = "17/24", srt = 90)
#' text(x=31/40, y=0.1, labels = "31/40", srt = 90)
#' text(x=33/40, y=0.1, labels = "33/40", srt = 90)
#' 
#' # Prepare dta and visualize
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

