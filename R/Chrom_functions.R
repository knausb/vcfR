

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
#' This vector is stored in the var.info slot of a Chrom object.
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
#' @export
#' @aliases proc_chrom
#' 
#' @param verbose logical stating whether to produce verbose output
#'
#' 
#proc.chrom <- function(x, pop1=NA, pop2=NA, win.size=1000, max.win=10000, verbose=TRUE){
proc_chrom <- function(x, verbose=TRUE, ...){
  stopifnot(class(x) == "Chrom")
  #  x <- set.pop1(x, pop1)
  #  x <- set.pop2(x, pop2)
  #  ptime <- system.time(x@acgt.w <- regex.win(x))
  ptime <- system.time(x@seq.info$nuc.win <- regex.win(x))
  if(verbose==TRUE){
    message("Nucleotide regions complete.\n")
    print(ptime)
  }
  ptime <- system.time(x@seq.info$N.win <- regex.win(x, regex="[n]"))
  #  ptime <- system.time(x@n.w <- acgt.win(x, regex="[n]"))
  #  ptime <- system.time(x <- n.win(x))
  if(verbose==TRUE){
    message("N regions complete.\n")
    print(ptime)
  }
  ptime <- system.time(x@win.info <- var.win(x, ...))
  if(verbose==TRUE){
    message("Window analysis complete.\n")
    print(ptime)
  }
  #
  if(nrow(x@vcf.gt[x@var.info$mask,])>0){
    ptime <- system.time(x <- gt2popsum(x))
    if(verbose==TRUE){
      message("Population summary complete.\n")
      print(ptime)
    }
  }
  #  ptime <- system.time(x <- windowize(x, win.size=win.size, max.win=max.win))
  #  ptime <- system.time(x <- windowize(x))
  #  if(verbose==TRUE){
  #    cat("Sliding windows created.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- gc.win(x))
  #  if(verbose==TRUE){
  #    cat("Sliding GC windows complete.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- snp.win(x))
  #  if(verbose==TRUE){
  #    cat("Sliding SNP windows complete.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- vcf.fix2gt.m(x))
  #  if(verbose==TRUE){
  #    cat("Genotype matrix complete.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- gt.m2sfs(x))
  #  cat("gt.m2sfs is commented out\n")
  #  if(verbose==TRUE){
  #    cat("SFS complete.\n")
  #    print(ptime)
  #  }
  #  ptime <- system.time(x <- linkage(x))
  #  if(verbose==TRUE){
  #    cat("Linkage calculation complete.\n")
  #    print(ptime)
  #  }
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


# EOF.