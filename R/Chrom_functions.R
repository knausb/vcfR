

#' @title Chrom_functions
#' @aliases Chrom functions
#' @rdname Chrom_functions
#'  
#' 
#' 
#' 



##### Set a mask #####

#' @rdname Chrom_functions
#' @export
#' @aliases masker
#' 
#' @param QUAL variant quality
#' @param mindp minimum cumulative depth
#' @param maxdp maximum cumulative depth
#' @param minmq minimum mapping quality
#' @param maxmq maximum mapping quality
#' @param ... arguments to be passed to children functions
#' 
#masker <- function(x, QUAL=999, mindp=0.25, maxdp=0.75, minmq=0.25, maxmq=0.75, ...){
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
#' 
window_table <- function(x){
  tab <- x@win.info
  tab <- cbind(rep(x@name, times=nrow(tab)), tab)
  names(tab)[1] <- "chrom"
  tab
}


# EOF.