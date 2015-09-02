#' @title Convert vcfR objects to other formats
#' @name Format conversion
#' @rdname vcf_conversion
#' @description
#' Convert vcfR objects to objects supported by other packages
#'  
#' @param x an object of class Chrom or vcfR
#' 
#' @details 
#' After processing vcf data in vcfR, one will likely proceed to an analysis step.
#' Within R, three obvious choices are:
#' \href{http://cran.r-project.org/web/packages/pegas/index.html}{pegas},
#' \href{http://cran.r-project.org/web/packages/adegenet/index.html}{adegenet} 
#' and \href{http://cran.r-project.org/web/packages/poppr/index.html}{poppr}.
#' The package pegas uses objects of type loci.
#' The function vcfR2loci calls extract.gt to create a matrix of genotypes which is then converted into an object of type loci.
#' The packages adegenet and poppr use the genind object.
#' The function vcfR2genind uses extract.gt to create a matrix of genotypes and uses the adegenet function df2genind to create a genind object.
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
#'
#'



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


#' @rdname vcf_conversion
#' @aliases vcfR2DNAbin
#' 
#' @param extract_indels logical, at present, the only option is TRUE
#' @param consensus logical, at present, the only option is TRUE
#' @param extract_haps logical specifying whether to separate genotype into alleles based on a delimiting character
#' @param gt_split character to delimit alleles within genotypes
#' @param verbose logical specifying whether to produce verbose output
#' 
#' @details
#' The DNAbin object stores nucleotide sequence information.
#' This means that in order to convert vcf data to a nucleotide representation.
#' For haploid data, this is straight forward.
#' For diploid data there is the option of converting heterozygous sites to IUPAC ambiguity codes.
#' This results in one sequence per diploid individual.
#' Note that functions called downstream of this choice may handle IUPAC ambiguity codes in unexpected manners.
#' If you have phased data, an alternative is to extract the haplotypes from each genotype.
#' This should work for diploids and higher ploids.
#' 
#' 
#' @export
vcfR2DNAbin <- function( x, extract_indels = TRUE , consensus = TRUE, extract_haps = FALSE, gt_split="|", verbose = TRUE )
{
  if( class(x) == 'Chrom' )
  {
    x <- x@vcf
  }
  if( class(x) != 'vcfR' )
  {
    stop( "Expecting an object of class Chrom or vcfR" )
  }
  if( consensus == TRUE & extract_haps == TRUE){
    stop("consensus and extract_haps both set to true. These options are incompatible. A haplotype should not be ambiguous.")
  }
  
  x <- extract_indels(x)
  if( nrow(x@fix) < 1 ){
    return( NA )
  }
  
  if( extract_haps == FALSE){
    x <- extract.gt(x, return.alleles=TRUE)
  } else {
    x <- extract_haps(x, gt_split="|", verbose = verbose)
  }
  
  # If we have no variants (row) return NA.
  if( nrow(x) < 1 )
  {
    return( NA )
  }

  # Strategies to convert genotypes (with a delimiter) to a nucleotide.  
  ploid <- unlist( strsplit( x[!is.na(x)][1], split=gt_split ) )
  if( length(ploid) == 1 )
  {
    x[is.na(x)] <- 'n'
    x <- apply(x, MARGIN=2, tolower)
  } else if ( length(ploid) == 2 & consensus == TRUE )
  {
    x <- alleles_to_consensus( x, sep = gt_split, NA_to_n = TRUE )
    x <- apply(x, MARGIN=2, tolower)
  } else if (ploid == 2 & consensus == FALSE ){
    stop( "Ploidy is 2 but consensus is set to FALSE.\nGenotypes must be converted to a nucleotide.\nConsider setting consensus to TRUE." )
#    stop( "function only valid for diploid and haploid genotypes." )
  } else if ( ploid > 2 & extract_haps == FALSE ){
    stop( "Ploidy is greater than 2 but extract_haps is set to FALSE.\nGenotypes must be converted to a nucleotide.\nConsider setting extract_haps to TRUE." )
  }
  
  # DNAbin characters must be lower case.
  x <- ape::as.DNAbin(t(x))
  x
}


