#' @title Extract elements from the GT section of a vcf format object
#' @rdname extract_gt
#' 
#' @param x An object of class chrom, vcfR or data.frame 
#' @param element element to extract from vcf genotype data. Common options include "DP", "GT" and "GQ"
#' @param mask a logical vector indicating which variants (rows) to include
#' @param as.matrix attempt to recast as a numeric matrix
#' 
#' @export
#' 
extract.gt <- function(x, element="GT", mask=logical(0), as.matrix=FALSE){
  if(class(x) == "Chrom"){
    vcf <- new(Class="vcfR")
    vcf@meta <- x@vcf.meta
    vcf@fix  <- x@vcf.fix
    vcf@gt   <- x@vcf.gt
    mask     <- x@var.info$mask
    x <- vcf
    rm(vcf)
  }
  if(class(x) != "vcfR"){stop("Expected an object of class vcfR or Chrom")}
  #
  # Make element regex more spercific.
  element <- paste("^[:]{0,1}", element, "[:]{0,1}$", sep="")
  #
  # Manage mask.
  if(length(mask) == 0){
    mask <- 1:nrow(x@gt)
  } else if (length(mask) > 0){
    # Use specified mask.
  } else if (sum(mask) > 0){
    # Use the mask in the Chom object.
    mask <- mask
  } else {
    stop("Unexpected mask.")
  }
  #
  # Create a function to get a single elements for a single variant (row).
  get.gt1 <- function(x, element="GT"){
    FORMAT <- unlist(strsplit(as.character(x[1]), ":"))
    x <- x[-1]
    pos <- grep(element, FORMAT)
    if(length(pos) == 0){
      x <- rep(NA, times=length(x))
    } else {
      x <- unlist(lapply(strsplit(as.character(x), ":"), function(x){x[pos]}))
    }
    is.na(x[x=="./."]) <- TRUE
    return(x)
  }
  #
  # Implement this function over the matrix.
  #  gt <- t(apply(x@gt[mask,], MARGIN=1, get.gt1, element=element))
  gt <- t(apply(x@gt, MARGIN=1, get.gt1, element=element))
  colnames(gt) <- names(x@gt)[-1]
  rownames(gt) <- rownames(x@gt)
  if(as.matrix==TRUE){
    tmp <- matrix(nrow=nrow(gt), ncol=ncol(gt))
    for(i in 1:ncol(gt)){
      tmp[,i] <- as.numeric(gt[,i]) 
    }
    colnames(tmp) <- colnames(gt)
    rownames(tmp) <- rownames(gt)
    gt <- tmp
  }
  return(gt)
}


#' @rdname extract_gt
#' 
#' @param as.numeric Logical, should the matrix be converted to numerics
#' @export
#extract.gt2 <- function(x, element="GT", mask=logical(0), as.matrix=FALSE){
extract.gt2 <- function(x, element="GT", mask=logical(0), as.numeric=FALSE){
  if(class(x) != "chrom" & class(x) != "vcfR" & class(x) != "data.frame"){
    stop("Expected an object of class chrom, vcfR or data.frame")
  }
  
  if(class(x) == "chrom"){
    x <- chrom_to_vcfR(x)
  }
  
  if(class(x) == "vcfR"){
#    outM <- .Call('vcfR_extractGT2NM', PACKAGE = 'vcfR', x@gt, element)
    outM <- .Call('vcfR_extract_GT_to_CM', PACKAGE = 'vcfR', x@gt, element)
  }
  
  if(class(x) == "data.frame"){
    outM <- .Call('vcfR_extractGT2NM', PACKAGE = 'vcfR', x, element)
    outM <- .Call('vcfR_extract_GT_to_CM', PACKAGE = 'vcfR', x, element)
  }

  if(as.numeric == TRUE){
    outM <- .Call('vcfR_CM_to_NM', PACKAGE = 'vcfR', outM)
  }

  return(outM)
}

