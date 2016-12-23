
#' @title Extract elements from vcfR objects
#' 
#'  
#' @rdname extract_gt
#' 
#' @description 
#' Extract elements from the 'gt' slot, convert extracted genotypes to their allelic state, extract indels from the data structure or extract elements from the INFO column of the 'fix' slot.
#' 
#' @param x An object of class chromR or vcfR 
#' @param element element to extract from vcf genotype data. Common options include "DP", "GT" and "GQ"
#' @param mask a logical indicating whether to apply the mask (TRUE) or return all variants (FALSE). Alternatively, a vector of logicals may be provided.
# @param as.matrix attempt to recast as a numeric matrix
#' @param verbose should verbose output be generated
#' @param as.numeric logical, should the matrix be converted to numerics
#' @param return.alleles logical indicating whether to return the genotypes (0/1) or alleles (A/T)
#' @param IDtoRowNames logical specifying whether to use the ID column from the FIX region as rownames
# @param allele.sep character which delimits the alleles in a genotype (/ or |), here this is not used for a regex (as it is in other functions)
#' @param extract logical indicating whether to return the extracted element or the remaining string
#' @param convertNA logical indicating whether to convert "." to NA.
#' 
#' @details
#' 
#' The function \strong{extract.gt} isolates elements from the 'gt' portion of vcf data.
#' Fields available for extraction are listed in the FORMAT column of the 'gt' slot.
#' Because different vcf producing software produce different fields the options will vary by software.
#' The mask parameter allows the mask to be implemented when using a chromR object.
#' The 'as.numeric' option will convert the results from a character to a numeric.
#' Note that if the data is not actually numeric, it will result in a numeric result which may not be interpretable.
#' The 'return.alleles' option allows the default behavior of numerically encoded genotypes (e.g., 0/1) to be converted to their nucleic acid representation (e.g., A/T).
# The allele.sep parameter allows the genotype delimiter to be specified.
#' Note that this is not used for a regular expression as similar parameters are used in other functions.
#' Extract allows the user to extract just the specified element (TRUE) or every element except the one specified.
#' 
#' Note that when 'as.numeric' is set to 'TRUE' but the data are not actually numeric, unexpected results will likely occur.
#' For example, the genotype field will typically be populated with values such as "0/1" or "1|0".
#' Although these may appear numeric, they contain a delimiter (the forward slash or the pipe) that is non-numeric.
#' This means that there is no straight forward conversion to a numeric and unexpected values should be expected.
#' 
#' 
#' @seealso
#' \code{\link{is.polymorphic}}
#' 
#' 
#' @examples 
#' data(vcfR_test)
#' gt <- extract.gt(vcfR_test)
#' gt <- extract.gt(vcfR_test, return.alleles = TRUE)
#' 
#' @export
extract.gt <- function(x, element="GT", 
                       mask=FALSE,
                       as.numeric=FALSE, 
                       return.alleles=FALSE,
                       IDtoRowNames = TRUE,
#                       allele.sep="/",
                       extract = TRUE,
                       convertNA = TRUE ){

  # Validate that we have an expected data structure
  if( class(x) != "chromR" & class(x) != "vcfR" ){
    stop( "Expected an object of class chromR or vcfR" )
  }
  
  # Catch unreasonable mask specification.
  if(class(x) == "vcfR"){
    if(length(mask) == 1 && mask == TRUE){
      # This condition does not appear to make 
      # sense and should be overridden.
      mask <- FALSE
    }
  }
  
  # If of class chromR, extract the vcf
  if(class(x) == "chromR"){
    tmpMask <- x@var.info$mask
    x <- x@vcf
  }

  # If a mask was specified in the call,
  # override the one from var.info
  if(length(mask) > 1){
    tmpMask <- mask
    mask <- TRUE
  }
  
  # Validate that the gt slot is a matrix
  if( class(x@gt) != "matrix" ){
    stop( paste("gt slot expected to be of class matrix. Instead found class", class(x@gt)) )
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
    outM <- .Call('vcfR_extract_GT_to_CM2', PACKAGE = 'vcfR',
                  x@fix,
                  x@gt,
                  element,
                  return.alleles, 
                  as.integer(extract), 
                  convertNA = as.numeric(convertNA) )
  }

  # If as.numeric is true, convert to a numeric matrix.
  if(as.numeric == TRUE){
    outM <- .Call('vcfR_CM_to_NM', PACKAGE = 'vcfR', outM)
  }
  
  # 
  if( IDtoRowNames == TRUE ){
    if( sum(is.na(x@fix[,'ID'])) > 0 ){
      x <- addID(x)
    }
    if( length(unique(x@fix[,'ID'])) != nrow(x@fix) ){
      stop('ID column contains non-unique names')
    }
    rownames(outM) <- x@fix[,'ID']
  }

  # Apply mask.
  if(mask == TRUE){
    outM <- outM[tmpMask,]
  }

  return(outM)
}



#' @rdname extract_gt
#' @aliases extract.haps
# @param gt.split character which delimits alleles in genotypes
#' @param unphased_as_NA logical specifying how to handle unphased genotypes
#' 
#' @details 
#' The function \strong{extract.haps} uses extract.gt to isolate genotypes.
#' It then uses the information in the REF and ALT columns as well as an allele delimiter (gt_split) to split genotypes into their allelic state.
#' Ploidy is determined by the first non-NA genotype in the first sample.
#' 
#' The VCF specification allows for genotypes to be delimited with a '|' when they are phased and a '/' when unphased.
#' This becomes important when dividing a genotype into two haplotypes.
#' When the alleels are phased this is straight forward.
#' When the alleles are unphased it presents a decision.
#' The default is to handle unphased data by converting them to NAs.
#' When unphased_as_NA is set to TRUE the alleles will be returned in the order they appear in the genotype.
#' This does not assign each allele to it's correct chromosome.
#' It becomes the user's responsibility to make informed decisions at this point.
#' 
#' 
#' @export
#extract.haps <- function(x, mask=FALSE, gt.split="|",verbose=TRUE){
extract.haps <- function(x, 
                         mask=FALSE, 
                         unphased_as_NA = TRUE, 
                         verbose=TRUE ){
  if(class(x) == "chromR"){
    if(length(mask) == 1 && mask==TRUE){
      x <- chromR2vcfR(x, use.mask = TRUE)
    } else {
      #      x <- chrom_to_vcfR(x)
      x <- x@vcf
    }
  }
  
  if(length(mask) > 1){
    x <- x[mask,]
  }

  # Determine ploidy  
  first.gt <- unlist(strsplit(x@gt[,-1][!is.na(x@gt[,-1])][1], ":"))[1]
#  ploidy <- length(unlist(strsplit(first.gt, split = gt.split, fixed = TRUE )))
  ploidy <- length(unlist(strsplit(first.gt, split = "[\\|/]" )))


  if( nrow( x@fix ) == 0 ){
    # No variants, return empty matrix.
    haps <- x@gt[ 0, -1 ]
  } else if ( ploidy == 1 ){
    haps <- extract.gt( x )
  } else if ( ploidy > 1 ) {
    gt <- extract.gt( x )
#    haps <- .Call('vcfR_extract_haps', PACKAGE = 'vcfR', 
#                  x@fix[,'REF'], x@fix[,'ALT'], 
#                  gt, gt.split, as.numeric(verbose))
    haps <- .Call('vcfR_extract_haps', PACKAGE = 'vcfR', 
                  x@fix[,'REF'], x@fix[,'ALT'],
                  gt, as.numeric(unphased_as_NA), as.numeric(verbose))
  } else {
    stop('Oops, we should never arrive here!')
  }

  haps
}




#' @rdname extract_gt
#' 
#' @aliases extract.indels
#' 
#' @param return.indels logical indicating whether to return indels or not
#' 
#' @details 
#' The function \strong{extract.indels} is used to remove indels from SNPs.
#' The function queries the 'REF' and 'ALT' columns of the 'fix' slot to see if any alleles are greater than one character in length.
#' When the parameter return_indels is FALSE only SNPs will be returned.
#' When the parameter return_indels is TRUE only indels will be returned.
#' 
#' 
#' @examples
#' data(vcfR_test)
#' getFIX(vcfR_test)
#' vcf <- extract.indels(vcfR_test)
#' getFIX(vcf)
#' vcf@fix[nrow(vcf@fix),'ALT'] <- ".,A"
#' vcf <- extract.indels(vcf)
#' getFIX(vcf)
#' 
#' data(vcfR_test)
#' extract.haps(vcfR_test, unphased_as_NA = FALSE)
#' extract.haps(vcfR_test)
#' 
#' 
#' @export
extract.indels <- function(x, return.indels=FALSE){
  if(class(x) == 'chromR'){
    x <- x@vcf
  }
  if(class(x) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or chromR.")
  }

  # Check reference for indels
  mask <- nchar(x@fix[,'REF']) > 1
  # Check reference for missing data.
#  mask[ grep(".", x@fix[,'REF'], fixed = TRUE) ] <- TRUE
  
  # Check alternate for indels
  mask[ unlist( lapply(
          strsplit(x@fix[,'ALT'], split=","), 
          function(x){ max(nchar(x)) > 1 }
  ) ) ] <- TRUE
  # Check alternate for missing data


  if(return.indels == FALSE){
    x <- x[ !mask, , drop = FALSE ]
  } else {
    x <- x[ mask, , drop = FALSE ]
  }
  
  return(x)  
}



#' @rdname extract_gt
#' @aliases extract.info
#' 
#' @details 
#' The function \strong{extract.info} is used to isolate elements from the INFO column of vcf data.
#' 
#' @export
extract.info <- function(x, element, as.numeric=FALSE, mask=FALSE){
  
  if( class(x) == 'chromR' ){
    mask <- x@var.info$mask
    x <- x@vcf
  }
  if( class(x) != 'vcfR' ){
    stop("Expecting an object of class vcfR or chromR.")
  }
  
#  values <- unlist(
#    lapply(strsplit(unlist(
#      lapply(strsplit(x@fix[,'INFO'], split=";"),             
#             function(x){grep(paste("^", element, "=", sep=""), x, value=TRUE)})),
#      split="="), function(x){x[2]})
#  )
  
  values <- strsplit(x@fix[,'INFO'], split=";")
  values <- lapply(values, function(x){grep(paste("^", element, "=", sep=""), x, value=TRUE)})
  values <- lapply(values, function(x){ unlist( strsplit(x, split="=") ) })
  values <- lapply(values, function(x){x[2]})
  values <- lapply(values, function(x){ if(is.null(x)){NA}else{x} })
  values <- unlist(values)
  

  if(as.numeric == TRUE){
    values <- as.numeric(values)
  }
#  if( mask != FALSE & !is.null(mask) ){
#    values <- values[x@var.info$mask]
#    values <- values[mask]
#  }
  values
}





##### ##### ##### ##### #####
# EOF.