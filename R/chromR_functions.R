

#' @title chromR_functions
#' @name chromR functions
# @aliases chromR functions
#' @rdname chromR_functions
#' @description Functions which act on chromR objects 



##### Set a mask #####

#' @rdname chromR_functions
#' @export
#' @aliases masker
#' 
#' @param min_QUAL minimum variant quality
#' @param min_DP minimum cumulative depth
#' @param max_DP maximum cumulative depth
#' @param min_MQ minimum mapping quality
#' @param max_MQ maximum mapping quality
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
#' This vector is stored in the var.info$mask slot of a chromR object.
#' 
#masker <- function(x, min_QUAL=999, min_DP=0.25, max_DP=0.75, minmq=20, maxmq=50, ...){
masker <- function(x, min_QUAL=1, min_DP=1, max_DP=1e4, min_MQ=20, max_MQ=100, ...){
  quals <- getQUAL(x)
#  quals  <- x@vcf.fix$QUAL
#  info <- x@var.info[,grep("DP|MQ",names(x@var.info)), drop=FALSE]
#  mask <- rep(TRUE, times=nrow(info))
  mask <- rep(TRUE, times=nrow(x@var.info))
  
  # Mask on QUAL
  if(sum(is.na(quals)) < length(quals)){
    mask[quals < min_QUAL] <- FALSE
  }

  # Mask on DP
  if( !is.null( x@var.info$DP ) ){
    if(sum(is.na(x@var.info$DP)) < length(x@var.info$DP)){
      mask[x@var.info$DP < min_DP] <- FALSE
      mask[x@var.info$DP > max_DP] <- FALSE
    }
  }
  
  if( !is.null( x@var.info$MQ ) ){
    if(sum(is.na(x@var.info$MQ)) < length(x@var.info$MQ)){
      mask[x@var.info$MQ < min_MQ] <- FALSE
      mask[x@var.info$MQ > max_MQ] <- FALSE
    }
  }
  x@var.info$mask <- mask
  return(x)
}









#' @rdname chromR_functions
#' 
#' @param x object of class chromR
#' 
#' @export
#' @aliases variant.table
#' 
#' @details
#' The function \strong{variant.table} creates a data.frame containing information about variants.
#' 
variant.table <- function(x){
  tab <- x@var.info[x@var.info$mask,]
#  tab <- cbind(rep(x@name, times=nrow(tab)), x@vcf.fix$QUAL[x@var.info$mask], tab)
#  names(tab)[1] <- "chrom"
#  names(tab)[2] <- "QUAL"
  tab
}




#' @rdname chromR_functions
#' @export
#' @aliases win.table
#' @details
#' The funciton \strong{win.table}
#' 
win.table <- function(x){
  tab <- x@win.info
  tab <- cbind(rep(x@name, times=nrow(tab)), tab)
  names(tab)[1] <- "chrom"
  tab
}




# EOF.
