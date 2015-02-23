
#' @title Convert chrom objects to vcfR objects
#' @rdname chrom_to_vcfR
#' @export
#'
#' @description
#' Convert chrom objects to vcfR objects.
#'
#' @param x Object of class chrom
#' @param use.mask Logical, determine if mask from chrom object should be used to subset vcf data
#'
#' @details
#' The chrom object is subset and recast as a vcfR object.  When use.mask is set 
#' to TRUE (the default), the object is subset to only the variants (rows) indicated
#' to include by the mask.  When use.mask is set to FALSE, all variants (rows) from 
#' the chrom object are included in the new vcfR object.
#' 
#' @return Returns an object of class vcfR. 
#' 
#'
chrom_to_vcfR <- function(x, use.mask=TRUE){
  if(class(x) != "chrom"){
    stop("Unexpected class! Expecting an object of class chrom.")
  }

  vcf <- new(Class = "vcfR")
  vcf@meta <- x@vcf.meta

  if(use.mask == TRUE){
    vcf@fix  <- x@vcf.fix[x@mask,]
    vcf@gt   <- x@vcf.gt[x@mask,]
  } else {
    vcf@fix  <- x@vcf.fix
    vcf@gt   <- x@vcf.gt
  }

  return(vcf)
}


