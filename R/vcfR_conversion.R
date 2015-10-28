#' @title Convert vcfR objects to other formats
#' @name Format conversion
#' @rdname vcfR_conversion
#' @description
#' Convert vcfR objects to objects supported by other R packages
#'  
#' @param x an object of class chromR or vcfR
#' 
#' @details 
#' After processing vcf data in vcfR, one will likely proceed to an analysis step.
#' Within R, three obvious choices are:
#' \href{http://cran.r-project.org/web/packages/pegas/index.html}{pegas},
#' \href{http://cran.r-project.org/web/packages/adegenet/index.html}{adegenet} 
#' and \href{http://cran.r-project.org/web/packages/poppr/index.html}{poppr}.
#' The package pegas uses objects of type loci.
#' The function \strong{vcfR2loci} calls extract.gt to create a matrix of genotypes which is then converted into an object of type loci.
#' 
#' The packages adegenet and poppr use the genind object.
#' The function \strong{vcfR2genind} uses extract.gt to create a matrix of genotypes and uses the adegenet function df2genind to create a genind object.
#' The package poppr additionally uses objects of class genclone.
#' A genind object can be converted to a genclone object with the function poppr::as.genclone.
#' 
#' 
#' 
#' 
#' 
#' @seealso
#' \code{\link{extract.gt}},
#' \code{\link[adegenet]{df2genind}},
#' \code{\link[adegenet]{genind}},
#' \href{http://cran.r-project.org/web/packages/pegas/index.html}{pegas},
#' \href{http://cran.r-project.org/web/packages/adegenet/index.html}{adegenet},
#' and 
#' \href{http://cran.r-project.org/web/packages/poppr/index.html}{poppr}.
#' To convert to objects of class \strong{DNAbin} see \code{\link{vcfR2DNAbin}}.
#'



#' @rdname vcfR_conversion
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


#' @rdname vcfR_conversion
#' @aliases vcfR2loci
#' 
#' @export
vcfR2loci <- function(x)
{
#  if(class(x) == "chromR")
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



