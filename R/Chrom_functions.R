

#' @title Chrom_functions
#' @name Chrom functions
# @aliases Chrom functions
#' @rdname Chrom_functions
#' @description Functions which act on Chrom objects 



##### Set a mask #####

#' @rdname Chrom_functions
#' @export
#' @aliases masker
#' 
#' @param QUAL minimum variant quality
#' @param mindp minimum cumulative depth
#' @param maxdp maximum cumulative depth
#' @param minmq minimum mapping quality
#' @param maxmq maximum mapping quality
#' @param ... arguments to be passed to methods
#' 
#' @details
#' The function \strong{masker} creates a logical vector that determines which variants are masked.
#' By masking certain variants, instead of deleting them, it preserves the dimensions of the data structure until a change needs to be committed.
#' Variants can be masked based on the value of the QUAL column of the vcf object.
#' Experience seems to show that this value is either at its maximum (999) or a rather low value.
#' The maximum and minimum sequence depth can also be used (mindp and maxdp).
#' The default is to mask all variants with depths of less than the 0.25 quantile and greater than the 0.75 quantile (these are also known as the lower and upper quartile).
#' The minimum and maximum mapping qualities (minmq, maxmq) can also be used.
#' 
#' 
#' This vector is stored in the var.info$mask slot of a Chrom object.
#' 
masker <- function(x, QUAL=999, mindp=0.25, maxdp=0.75, minmq=20, maxmq=50, ...){
  quals  <- x@vcf.fix$QUAL
  info <- x@var.info[,grep("DP|MQ",names(x@var.info))]
  mask <- rep(TRUE, times=nrow(info))
  #
  if(sum(is.na(quals)) < length(quals)){
    mask[quals < QUAL] <- FALSE
  }
  #  if(sum(is.na(x@vcf.info$DP)) < length(x@vcf.info$DP)){
  #    mask[x@vcf.info$DP < quantile(x@vcf.info$DP, probs=c(mindp))] <- FALSE
  #    mask[x@vcf.info$DP > quantile(x@vcf.info$DP, probs=c(maxdp))] <- FALSE
  #  }
  #  if(sum(is.na(x@vcf.info$MQ)) < length(x@vcf.info$MQ)){
  #    mask[x@vcf.info$MQ < quantile(x@vcf.info$MQ, probs=c(minmq))] <- FALSE
  #    mask[x@vcf.info$MQ > quantile(x@vcf.info$MQ, probs=c(maxmq))] <- FALSE
  #  }
  if(sum(is.na(info$DP)) < length(info$DP)){
    mask[info$DP < quantile(info$DP, probs=c(mindp))] <- FALSE
    mask[info$DP > quantile(info$DP, probs=c(maxdp))] <- FALSE
  }
  if(sum(is.na(info$MQ)) < length(info$MQ)){
    mask[info$MQ < minmq] <- FALSE
    mask[info$MQ > maxmq] <- FALSE
    #    mask[info$MQ < quantile(info$MQ, probs=c(minmq))] <- FALSE
    #    mask[info$MQ > quantile(info$MQ, probs=c(maxmq))] <- FALSE
  }
  x@var.info$mask <- mask
  return(x)
}









#' @rdname Chrom_functions
#' 
#' @param x object of class Chrom
#' 
#' @export
#' @aliases variant_table
#' 
#' @details
#' The function \strong{variant_table} creates a data.frame containing information about variants.
#' 
variant_table <- function(x){
  tab <- x@var.info[x@var.info$mask,]
  tab <- cbind(rep(x@name, times=nrow(tab)), x@vcf.fix$QUAL[x@var.info$mask], tab)
  names(tab)[1] <- "chrom"
  names(tab)[2] <- "QUAL"
  tab
}




#' @rdname Chrom_functions
#' @export
#' @aliases window_table
#' @details
#' The funciton \strong{window_table}
#' 
window_table <- function(x){
  tab <- x@win.info
  tab <- cbind(rep(x@name, times=nrow(tab)), tab)
  names(tab)[1] <- "chrom"
  tab
}


#' @rdname Chrom_functions
#' @aliases window_table
#' 
#' @export
extract_indels <- function(x, return_indels=FALSE){
  if(class(x) == 'Chrom'){
    # Recast as a vcfR object.
    Chrom <- x
    x <- new(Class="vcfR")
    x@meta <- temp@vcf.meta
    x@fix <- temp@vcf.fix
    x@gt <- temp@vcf.gt
  }
  if(class(x) != "vcfR"){
    stop("Unexpected class! Expecting an object of class vcfR or Chrom.")
  }
  
  mask <- nchar(x@fix$REF) > 1
  mask[unlist(lapply(strsplit(x@fix$REF, split=","), function(x){max(nchar(x))})) > 1] <- TRUE

  if(return_indels == FALSE){
    x <- x[!mask,]
  } else {
    x <- x[mask,]
  }
  
  if(length(grep("Chrom", ls())) > 0){
    Chrom@vcf.fix <- x@fix
    Chrom@vcf.gt <- x@gt
    return(Chrom)
  } else {
    return(x)  
  }
}


# EOF.