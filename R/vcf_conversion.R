#' @title Convert vcf data to other formats
#' @name Format conversion
#' @rdname vcf_conversion
#' @description
#' Convert vcfR object to other formats
#'  
#' @param x an object of class Chrom or vcfR



#' @rdname vcf_conversion
#' @aliases vcfR2genind
#' 
#' @param sep character (to be used in a regular expression) to delimit the alleles of genotypes
#' 
#' @export
vcfR2genind <- function(x, sep="[|/]") {
  x <- extract.gt(x)
  x <- adegenet::df2genind(t(x), sep=sep)
  x
}


#' @rdname vcf_conversion
#' @aliases vcfR2loci
#' 
#' @export
vcfR2loci <- function(x)
{
#  if(class(x) == "Chrom")
#  {
#    x <- x@vcf
#  }
  x <- extract.gt(x)
  # modified from pegas::as.loci.genind
  x <- as.data.frame(t(x))
  icol <- 1:ncol(x)
  for (i in icol) x[, i] <- factor(x[, i] )
  class(x) <- c("loci", "data.frame")
  attr(x, "locicol") <- icol
  x
}

