#' @title Convert vcfR objects to other formats
#' @name Format conversion
#' @rdname vcf_conversion
#' @description
#' Convert vcfR objects to objects supported by other packages
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


#' @rdname vcf_conversion
#' @aliases vcfR2DNAbin
#' 
#' @param extract.indels logical, at present, the only option is TRUE
#' @param consensus logical, at present, the only option is TRUE
#' @param extract.haps logical specifying whether to separate genotype into alleles based on a delimiting character
#' @param gt.split character to delimit alleles within genotypes
#' @param ref.seq reference sequence for the region being converted
#' @param start.pos chromosomal position for the start of the ref.seq
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
#' A multiple sequence alignment for an entire chromosome contains a high degree of redundancy in invariant sites.
#' This is why vcf files contain only information on variants.
#' However, many analyses make use of invariant sites.
#' For example, likelihood inference of phylogeny (and therefore Bayesian methods as well) requires invariant sites.
#' By default, the function \strong{vcfR2DNAbin} by default will return an object of class DNAbin which only contains the variable positions.
#' In order to fill in the invariant sites, the parameter ref.seq should be provided.
#' When vcfR2DNAbin encounters missing data in the vcf data (NA) it is coded as an ambiguous nucleotide (n) in the DNAbin object.
#' Providing an entire chromosome may exceed available memory.
#' By providing a ref.seq for a region (of class ape::DNAbin), for example a gene, and the start.pos for that sequence (the end position will be determined from the sequence length), complete sequences can be created.
#' 
#' 
#' @export
vcfR2DNAbin <- function( x, extract.indels = TRUE , consensus = TRUE,
                         extract.haps = FALSE, gt.split="|",
                         ref.seq = NULL, start.pos = NULL,
                         verbose = TRUE )
{
  # Sanitize input.
  if( class(x) == 'chromR' )
  {
    x <- x@vcf
  }
  if( class(x) != 'vcfR' )
  {
    stop( "Expecting an object of class chromR or vcfR" )
  }
  if( consensus == TRUE & extract.haps == TRUE){
    stop("consensus and extract_haps both set to true. These options are incompatible. A haplotype should not be ambiguous.")
  }
  if( !is.null(start.pos) & class(start.pos) == "character" ){
    start.pos <- as.integer(start.pos)
  }
  
  # Check and sanitize ref.seq.
  if( class(ref.seq) != 'DNAbin' & !is.null(ref.seq) ){
    stop( paste("expecting ref.seq to be of class DNAbin but it is of class", class(ref.seq)) )
  }
  if( is.list(ref.seq) ){
    ref.seq <- as.matrix(ref.seq)
  }
  
  # Extract indels.
  if( extract.indels == TRUE ){
    x <- extract.indels(x)
  } else {
    stop("extract.indels == FALSE is not currently implemented.")
  }
  
  # If we removed all variants, return NA.
  if( nrow(x@fix) < 1 ){
    return( NA )
  }
  
  # Save POS in case we need it.
  pos <- as.numeric(x@fix[,'POS'])
  
  if( extract.haps == FALSE ){
    x <- extract.gt(x, return.alleles=TRUE)
  } else if( extract.haps == TRUE ){
#    x <- extract_haps(x, gt_split=gt.split, verbose = verbose)
    x <- extract.haps(x, gt.split =gt.split, verbose = verbose)
  } else {
    stop( "Invalid specification of extract.haps.\nShould be a logical." )
  }
  
  # If we have no variants (rows) return NA.
  if( nrow(x) < 1 )
  {
    return( NA )
  }

  # Strategies to convert genotypes (with a delimiter) to nucleotides.
  # If extract.haps was set to TRUE, then our data is effectively haploid now.
  ploid <- unlist( strsplit( x[!is.na(x)][1], split=gt.split ) )
  if( length(ploid) == 1 )
  {
    x[is.na(x)] <- 'n'
    if( nrow(x) > 1 ){
      x <- apply(x, MARGIN=2, tolower)
    } else if( nrow(x) == 1 ){
      x <- apply(x, MARGIN=2, tolower)
      x <- matrix(x, nrow=1, dimnames = list( NULL, names(x)))
    }
  } else if ( length(ploid) == 2 & consensus == TRUE )
  {
    x <- alleles2consensus( x, sep = gt.split, NA_to_n = TRUE )
    x <- apply(x, MARGIN=2, tolower)
  } else if ( ploid == 2 & consensus == FALSE ){
    stop( "Ploidy is 2 but consensus is set to FALSE.\nGenotypes must be converted to a nucleotide.\nConsider setting consensus to TRUE." )
#    stop( "function only valid for diploid and haploid genotypes." )
  } else if ( ploid > 2 & extract_haps == FALSE ){
    stop( "Ploidy is greater than 2 but extract_haps is set to FALSE.\nGenotypes must be converted to a nucleotide.\nConsider setting extract_haps to TRUE." )
  }
  
  # Fill with reference sequence.
  if( !is.null(ref.seq) ){
    if( class(ref.seq) == "list" & inherits(ref.seq, "DNAbin") ){
      ref.seq <- as.matrix(ref.seq)
    }
    
    # Subset the vcf data to our region of interest.
    end.pos <- start.pos + dim(ref.seq)[2] - 1
#    if( nrow(x) > 1 ){
      x <- x[ pos >= start.pos & pos <= end.pos, , drop = FALSE ]
#    } else if( nrow(x) == 1){
#      x <- x[ pos >= start.pos & pos <= end.pos,]
#      x <- matrix(x, nrow=1, dimnames = list( NULL, names(x)))
#    }
    
    # Subset pos to our region of interest
    # and set to new coordinate system.
    pos <- pos[ pos >= start.pos & pos <= end.pos ]
    pos <- pos - start.pos + 1
    
    # Create a matrix with the same number of columns as vcf data
    # and the number of rows to match the sequence length.
    out.seq <- matrix( as.character(ref.seq),
                       nrow = dim(ref.seq)[2],
                       ncol = ncol(x),
                       byrow = FALSE
    )
    colnames(out.seq) <- colnames(x)
    
#    pos <- pos - start.pos + 1
    out.seq[pos,] <- x
    x <- out.seq
  }

  # DNAbin characters must be lower case.
  x <- ape::as.DNAbin(t(x))
  x
}


