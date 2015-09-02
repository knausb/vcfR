#' @title Extract elements from the GT section of a vcf format object
#' @rdname extract_gt
#' 
#' @param x An object of class Chrom or vcfR 
#' @param element element to extract from vcf genotype data. Common options include "DP", "GT" and "GQ"
#' @param mask a logical indicating whether to apply the mask (TRUE) or return all variants (FALSE). Alternatively, a vector of logicals may be provided.
# @param as.matrix attempt to recast as a numeric matrix
#' @param verbose should verbose output be generated
#' 
#' @details
#' Note that when 'as.numeric' is set to 'TRUE' but the data are not actually numeric, unexpected results will likely occur.
#' 
#' The function \strong{extract.gt} isolates elements from the GT portion of vcf data.
#' Fields available for extraction are listed in the FORMAT column of the GT portion.
#' Because different vcf producing software produce different fields the options will vary by software.
#' 
#' 
# @export
#' 
#' @rdname extract_gt
#' 
#' @param as.numeric Logical, should the matrix be converted to numerics
#' @param return.alleles logical indicating whether to return the genotypes (0/1) or alleles (A/T)
#' @param allele.sep character which delimits the alleles in a genotype (/ or |)
#' 
#' @export
extract.gt <- function(x, element="GT", mask=FALSE, as.numeric=FALSE, return.alleles=FALSE, allele.sep="/" ){

  # Validate that we have an expected data structure
  if( class(x) != "Chrom" & class(x) != "vcfR" ){
    stop( "Expected an object of class Chrom or vcfR" )
  }
  
  # Catch unreasonable mask specification.
  if(class(x) == "vcfR"){
    if(length(mask) == 1 && mask == TRUE){
      # This condition does not appear to make 
      # sense and should be overridden.
      mask <- FALSE
    }
  }
  
  # If of class Chrom, extract the vcf
  if(class(x) == "Chrom"){
    tmpMask <- x@var.info$mask
    x <- x@vcf
  }

  # If a mask was specified in the call,
  # override the one from var.info
  if(length(mask) > 1){
    tmpMask <- mask
    mask <- TRUE
  }

  if(as.numeric == TRUE & return.alleles == TRUE ){
    stop("Invalid parameter choice, as.numeric and return.alleles can't both be true, alleles are characters!")
  }
  
  # If of class vcfR, call compiled code to extract field.
  if(class(x) == "vcfR"){
    if(colnames(x@gt)[1] != "FORMAT"){
      stop("First column is not named 'FORMAT', this is essential information.")
    }
#    .Call('vcfR_extract_haps', PACKAGE = 'vcfR', ref, alt, gt, gt_split, verbose)
#    outM <- .Call('vcfR_extract_GT_to_CM', PACKAGE = 'vcfR', x@gt, element)
    outM <- .Call('vcfR_extract_GT_to_CM2', PACKAGE = 'vcfR', x@fix, x@gt, element, allele.sep, return.alleles )
    
  }

  # If as.numeric is true, convert to a numeric matrix.
  if(as.numeric == TRUE){
    outM <- .Call('vcfR_CM_to_NM', PACKAGE = 'vcfR', outM)
  }

  # Apply mask.
  if(mask == TRUE){
    outM <- outM[tmpMask,]
  }

  return(outM)
}




#' @rdname extract_gt
#' @aliases extract_indels
#' @param return_indels logical indicating whether to return indels or not
#' 
#' @details 
#' The function \strong{extract_indels} is used to isolate indels from SNPs.
#' When the parameter return_indels is FALSE only SNPs will be returned.
#' When the parameter return_indels is TRUE only indels will be returned.
#'
#' 
#' @export
extract_indels <- function(x, return_indels=FALSE){
  if(class(x) == 'Chrom'){
    x <- x@vcf
  }
  if(class(x) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or Chrom.")
  }

  # Check reference for indels
  mask <- nchar(x@fix[,'REF']) > 1
  
  # Check alternate for indels
  mask[unlist(lapply(strsplit(x@fix[,'ALT'], split=","), function(x){max(nchar(x))})) > 1] <- TRUE

  if(return_indels == FALSE){
    x <- x[ !mask, , drop = FALSE ]
  } else {
    x <- x[ mask, , drop = FALSE ]
  }
  
  return(x)  
}



#' @rdname extract_gt
#' @aliases extract_info
#' 
#' @details 
#' The function \strong{extract_info} is used to isolate elements from the INFO column of vcf data.
#' 
#' @export
extract_info <- function(x, element, as.numeric=FALSE, mask=FALSE){
  values <- unlist(
    lapply(strsplit(unlist(
      lapply(strsplit(x@vcf.fix$INFO, split=";"),
             function(x){grep(paste("^", element, "=", sep=""), x, value=TRUE)})),
      split="="), function(x){x[2]})
    )

  if(as.numeric == TRUE){
    values <- as.numeric(values)
  }
  if(mask==TRUE){
    values <- values[x@var.info$mask]
  }
  values
}




#' @rdname extract_gt
#' @aliases extract_haps
#' @param gt_split character which delimits alleles in genotypes
#' 
#' @details 
#' The function \strong{extract_haps} uses extract.gt to isolate genotypes.
#' It then uses the information in the REF and ALT columns as well as an allele delimiter (gt_split) to split genotypes into their allelic state.
#' 
#' @export
extract_haps <- function(x, mask=FALSE, gt_split="|",verbose=TRUE){
  if(class(x) == "Chrom"){
    if(length(mask) == 1 && mask==TRUE){
      x <- chrom_to_vcfR(x, use.mask = TRUE)
    } else {
#      x <- chrom_to_vcfR(x)
      x <- x@vcf
    }
  }
  
  if(length(mask) > 1){
    x <- x[mask,]
  }
  
  gt <- extract.gt(x, element="GT")

#  Rcpp::StringMatrix extract_haps(Rcpp::StringVector ref,
#                                  Rcpp::StringVector alt,
#                                  Rcpp::StringMatrix gt,
#                                  char gt_split,
#                                  int vebosity) {
  
  
  haps <- .Call('vcfR_extract_haps', PACKAGE = 'vcfR', x@fix[,'REF'], x@fix[,'ALT'], gt, gt_split, 1)
  haps
}



#' @rdname extract_gt
#' @aliases is_polymorphic
#' @param na.omit logical to omit missing data
#' 
#' @details 
#' The function \strong{is_polymorphic} returns a vector of logicals indicating whether a variant is polymorphic.
#' Variants in vcf files are always polymorphic.
#' However, if the variants are censored somehow (set to NA) or samples are removed a variant may become invariant.
#' 
#' 
#' @export
is_polymorphic <- function(x, na.omit=FALSE){
  if(class(x) != "vcfR"){
    stop("Expected an object of class vcfR")
  }
  x <- extract.gt(x)
  
  test_poly <- function(x, na.omit=na.omit){
    if(na.omit == TRUE){
      x <- na.omit(x)
    }
    sum(x[1] == x[-1]) < (length(x) - 1)
  }
  apply(x, MARGIN=1, test_poly, na.omit=na.omit)
}


#' @rdname extract_gt
#' @aliases is_biallelic
#' 
#' @details 
#' The function \strong{is_bialleleic} returns a vector of logicals indicating whether a variant is biallelic.
#' Some analyses or downstream analyses only work with biallelic loci.
#' This function can help manage this.
#' 
#' @export
is_biallelic <- function(x){
#  x <- as.character(x@fix$ALT)
  x <- as.character(x@fix[,'ALT'])
  x <- strsplit(x, split=",")
  lapply(x, length) == 1
}



#' @rdname extract_gt
#' @aliases get.alleles
#' 
#' @param split character passed to strsplit to split the genotype into alleles
#' @param na.rm logical indicating whether to remove NAs
#' 
#' @details 
#' The function \strong{get.alleles} takes a vector of genotypes and returns the unique alleles.
#' 
#' @export
get.alleles <- function( x, split="/", na.rm = FALSE, as.numeric = FALSE ){
  x <- unlist(strsplit(x, split))
  if(na.rm == TRUE){
    x <- x[ x != "NA" ]
  }
  if(as.numeric == TRUE){
    x <- as.numeric(x)
  }
  x <- unique(x)
  x
}


#' @rdname extract_gt
#' @aliases alleles_to_consensus
#' 
#' @param x2 a matrix of alleles as genotypes (e.g., A/A, C/G, etc.)
#' @param sep a character which delimits the alleles in a genotype (/ or |)
#' @param NA_to_n logical indicating whether NAs should be scores as n
#' 
#' @details 
#' The function \strong{alleles_to_consensus} converts genotypes to a single consensus allele using IUPAC ambiguity codes for heterozygotes. 
#' 
#' 
#' @export
alleles_to_consensus <- function( x2, sep = "/", NA_to_n = TRUE ){
  lookup <- cbind(paste(c('A','C','G','T', 'A','T','C','G', 'A','C','G','T', 'A','G','C','T'),
                        c('A','C','G','T', 'T','A','G','C', 'C','A','T','G', 'G','A','T','C'),
                        sep=sep),
                  c('a','c','g','t', 'w','w','s','s', 'm','m','k','k', 'r','r','y','y'))
    
  for(i in 1:nrow( lookup ))
  {
    x2[ x2 == lookup[i,1] ] <- lookup[i,2]
  }
  if( NA_to_n == TRUE )
  {
    x2[ is.na(x2) ] <- 'n'
  }
  
  x2
}


