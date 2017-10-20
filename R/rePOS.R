
#' @title create non-overlapping positions
#' @name rePOS
#' @rdname rePOS
#' 
#' @description
#' Converts allele balance data produced by \code{freq_peak()} to a copy number by assinging the allele balance data (frequencies) to its closest expected ratio.
#'  
#' @param x a vcfR object
#' @param lens a data.frame describing the reference
#' @param buff an integer indicating buffer length
#' 
#' @details
#' Each chromosome in a genome typically begins with position one.
#' This creates a problem when plotting the data associated with each chromosome because the information will overlap.
#' This function uses the information in the data.frame \code{lens} to create a new coordinate system where chromosomes do not overlap.
#' 
#' The data.frame \strong{lens} should have a row for each chromosome and two columns.
#' The first column is the name of each chromosome as it appears in the vcfR object.
#' The second column is the length of each chromosome.
#' The parameter \strong{buff} indicates the length of a buffer to put in between each chromosome.
#' 
#' 
#' @return a vector of integers that represent the new coordinate system.
#' 
#' 
#' @examples
#' # Create some VCF data.
#' data(vcfR_example)
#' vcf1 <-vcf[1:500,]
#' vcf2 <-vcf[500:1500,]
#' vcf3 <- vcf[1500:2533]
#' vcf1@fix[,'CHROM'] <- 'chrom1'
#' vcf2@fix[,'CHROM'] <- 'chrom2'
#' vcf3@fix[,'CHROM'] <- 'chrom3'
#' vcf2@fix[,'POS'] <- as.character(getPOS(vcf2) - 21900)
#' vcf3@fix[,'POS'] <- as.character(getPOS(vcf3) - 67900)
#' vcf <- rbind2(vcf1, vcf2)
#' vcf <- rbind2(vcf, vcf3)
#' rm(vcf1, vcf2, vcf3)
#' 
#' # Create lens
#' lens <- data.frame(matrix(nrow=3, ncol=2))
#' lens[1,1] <- 'chrom1'
#' lens[2,1] <- 'chrom2'
#' lens[3,1] <- 'chrom3'
#' lens[1,2] <- 2200
#' lens[2,2] <- 4700
#' lens[3,2] <- 32089
#' 
#' # Illustrate the issue.
#' plot(getPOS(vcf), dp, col=as.factor(getCHROM(vcf)))
#' 
#' newPOS <- rePOS(vcf, lens)
#' dp <- extract.info(vcf, element="DP", as.numeric=TRUE)
#' plot(newPOS, dp, col=as.factor(getCHROM(vcf)))
#' 
#' 
#' @export
rePOS <- function(x, lens, buff = 0){
  if( class(x) == 'chromR' ){
    x <- x@vcfR
  }
  if( class(x) != 'vcfR' ){
    msg <- paste('expecting a chromR or vcfR object, received instead a', class(x))
    stop(msg)
  }
  
  # Check CHROM names.
  if( sum(lens[,1] %in% getCHROM(x)) != nrow(lens) ){
    msg <- "chromosome (CHROM) names in vcfR object and lens do not appear to match"
    stop(msg)
  }
  
  # Update lens with new starts.
  lens$new_start <- lens$X2
  
  
  
  oldPOS <- getPOS(x)
  newPOS <- oldPOS
  oldCHROM <- getCHROM(x)
  
  # rep(rownames(myM), times=myM[,1])
  
  return(newPOS)
}

