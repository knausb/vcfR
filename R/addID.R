#' @title Populate the ID column of VCF data
#' @name addID
#' @rdname addID
#' 
#' @description
#' Populate the ID column of VCF data by concatenating the chromosome, position and optionally an index. 
#'  
#' @param x an object of class vcfR or chromR
#' 
#' @details 
#' Variant callers typically leave the ID column empty in VCF data.
#' However, the VCF data may potentially include variants with IDs as well as variants without.
#' This function populates the missing elements by concatenating the chromosome and position.
#' When this concatenation results in non-unique names, an index is added to force uniqueness.
#' 
#' @export
addID <- function(x){
  if( class(x) == 'chromR' ){
    ID <- x@vcf@fix[,'ID']
    CHROM <- x@vcf@fix[,'CHROM']
    POS <- x@vcf@fix[,'POS']
  } else if( class(x) == 'vcfR' ){
    ID <- x@fix[,'ID']
    CHROM <- x@fix[,'CHROM']
    POS <- x@fix[,'POS']
  } else {
    stop("expecting an object of class vcfR or chromR.")
  }
  
  if( sum(!is.na(ID)) < length(ID) ){
    ID[ is.na(ID) ] <- paste( CHROM[ is.na(ID) ], POS[ is.na(ID) ], sep="_" )
    if( length(unique(ID)) < length(ID) ){
      ID <- paste( ID, 1:length(ID), sep="_" )
    }
  }
  
  if( class(x) == 'chromR' ){
    x@vcf@fix[,'ID'] <- ID
  } else if( class(x) == 'vcfR' ){
    x@fix[,'ID'] <- ID
  } else {
    stop("expecting an object of class vcfR or chromR.")
  }
  
  return(x)
}



