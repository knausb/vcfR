#' @title Extract elements from the GT section of a vcf format object
#' @rdname extract_gt
#' 
#' @param x An object of class Chrom, vcfR or data.frame 
#' @param element element to extract from vcf genotype data. Common options include "DP", "GT" and "GQ"
#' @param mask a logical indicating whether to apply the mask (TRUE) or return all variants (FALSE). Alternatively, a vector of logicals may be provided.
#' @param as.matrix attempt to recast as a numeric matrix
#' 
#' 
#' @details
#' Note that when 'as.numeric' is set to 'TRUE' but the data are not actually numeric, unexpected results will likely occur.
#' 
#' 
#' @export
#' 
extract.gt.allR <- function(x, element="GT", mask=logical(0), as.matrix=FALSE){
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
extract.gt <- function(x, element="GT", mask=FALSE, as.numeric=FALSE){
  if(class(x) != "Chrom" & class(x) != "vcfR" & class(x) != "data.frame"){
    stop("Expected an object of class Chrom, vcfR or data.frame")
  }
  
  if(class(x) == "Chrom"){
    tmpMask <- x@var.info$mask
    x <- chrom_to_vcfR(x)
  }
  
  if(length(mask) > 1){
    tmpMask <- mask
    mask <- TRUE
  }
  
  if(class(x) == "vcfR" | class(x) == "data.frame"){
    if(length(mask) == 1 & mask == TRUE){
      # This condition does not appear to make 
      # sense and should be overridden.
      mask <- FALSE
    }
  }

  if(class(x) == "vcfR"){
#    outM <- .Call('vcfR_extractGT2NM', PACKAGE = 'vcfR', x@gt, element)
    if(names(x@gt)[1] != "FORMAT"){
      stop("First column is not named 'FORMAT', this is essential information.")
    }
    outM <- .Call('vcfR_extract_GT_to_CM', PACKAGE = 'vcfR', x@gt, element)
  }
  
  if(class(x) == "data.frame"){
    if(names(x)[1] != "FORMAT"){
      stop("First column is not named 'FORMAT', this is essential information.")
    }
#    outM <- .Call('vcfR_extractGT2NM', PACKAGE = 'vcfR', x, element)
    outM <- .Call('vcfR_extract_GT_to_CM', PACKAGE = 'vcfR', x, element)
  }

  if(as.numeric == TRUE){
    outM <- .Call('vcfR_CM_to_NM', PACKAGE = 'vcfR', outM)
  }

  if(mask == TRUE){
    outM <- outM[tmpMask,]
  }

  return(outM)
}




#' @rdname extract_gt
#' @aliases extract_indels
#' @param return_indels logical indicating whether to return indels or not
#' 
#' @export
extract_indels <- function(x, return_indels=FALSE){
  if(class(x) == 'Chrom'){
    stop('extract_indels only works on vcfR objects!')
    # Recast as a vcfR object.
#    Chrom <- x
#    x <- new(Class="vcfR")
#    x@meta <- Chrom@vcf.meta
#    x@fix  <- Chrom@vcf.fix
#    x@gt   <- Chrom@vcf.gt
  }
  if(class(x) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or Chrom.")
  }
  
  mask <- nchar(x@fix$REF) > 1
  mask[unlist(lapply(strsplit(x@fix$ALT, split=","), function(x){max(nchar(x))})) > 1] <- TRUE
  
  if(return_indels == FALSE){
    x <- x[!mask,]
  } else {
    x <- x[mask,]
  }
  
#  if(length(grep("Chrom", ls())) > 0){
#    Chrom@vcf.fix <- x@fix
#    Chrom@vcf.gt <- x@gt
#    return(Chrom)
#  } else {
    return(x)  
#  }
}



#' @rdname extract_gt
#' @aliases extract_info
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

