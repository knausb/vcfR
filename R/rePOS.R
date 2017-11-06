
#' @title Create non-overlapping positions (POS) for VCF data
#' @name rePOS
#' @rdname rePOS
#' 
#' @description
#' Converts allele balance data produced by \code{freq_peak()} to a copy number by assinging the allele balance data (frequencies) to its closest expected ratio.
#'  
#' @param x a vcfR object
#' @param lens a data.frame describing the reference
#' @param ret.lens logical specifying whether lens should be returned
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
#' 
#' The parameter \strong{buff} indicates the length of a buffer to put in between each chromosome.
#' This buffer may help distinguish chromosomes from one another.
#' 
#' In order to create the new coordinates the \code{lens} data.frame is updated with the new start positions.
#' The parameter \strong{}
#' 
#' 
#' @return Either a vector of integers that represent the new coordinate system or a list containing the vector of integers and the lens data.frame.
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
#' lens[1,2] <- 22000
#' lens[2,2] <- 47000
#' lens[3,2] <- 32089
#' 
#' # Illustrate the issue.
#' dp <- extract.info(vcf, element="DP", as.numeric=TRUE)
#' plot(getPOS(vcf), dp, col=as.factor(getCHROM(vcf)))
#' 
#' # Resolve the issue.
#' newPOS <- rePOS(vcf, lens)
#' dp <- extract.info(vcf, element="DP", as.numeric=TRUE)
#' plot(newPOS, dp, col=as.factor(getCHROM(vcf)))
#' 
#' # Illustrate the buffer
#' newPOS <- rePOS(vcf, lens, buff=10000)
#' dp <- extract.info(vcf, element="DP", as.numeric=TRUE)
#' plot(newPOS, dp, col=as.factor(getCHROM(vcf)))
#' 
#' 
#' @export
rePOS <- function(x, lens, ret.lens = FALSE, buff = 0){
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
  colnames(lens)[1:2] <- c('chrom', 'length')
  lens$new_start <- 0
  lens$new_start[1] <- 1
  lens$mids <- 0
  lens$mids[1] <- round(lens$length[1]/2)
  for(i in 2:nrow(lens)){
    lens$new_start[i] <- lens$new_start[i-1] + lens$length[i-1]
    # Apply buffer
    lens$new_start[i] <- lens$new_start[i] + buff
    # Midpoint
    lens$mids[i] <- lens$new_start[i] + lens$length[i]/2
  }
  
  # Apply new start to POS.
  oldPOS <- getPOS(x)
  oldCHROM <- getCHROM(x)
  # table converts our character vector to a factor.
  # This tends to sort things.
  # We want to retain the order so let's recast this ourselves.
  oldCHROM <- factor(oldCHROM, levels=unique(oldCHROM))
  myM <- as.matrix(table(oldCHROM))
  newPOS <- oldPOS + rep(lens$new_start, times=myM[,1]) - 1
  
  if(ret.lens == TRUE){
    return(list(newPOS=newPOS, lens=lens))
  } else {
    return(newPOS)
  }
}

