




thetas <- function(x){
  #  print(x)
  rnum <- x[1]
  anum <- x[2]
  if(is.na(rnum)){return(c(NA,NA,NA))}
  n <- rnum + anum
  Si <- vector(mode="numeric", length=n)
  Si[anum] <- 1
  theta_w <- sum(1/1:(rnum+anum-1))^-1 * 1
  theta_pi <- (2*anum*rnum)/(n*(n-1))
  theta_h <- (2*1*anum^2)/(n*(n-1))
  return(c(theta_pi, theta_w, theta_h))
}




#' @rdname pop_gen_sum
#' @export
#' @aliases gt2popsum
#' 
#' @param deprecated logical specifying whether to run the function (FALSE) or present deprecation message (TRUE).
#' 
gt2popsum <- function(x, deprecated = TRUE){
  if(class(x) != "chromR"){stop("Object is not of class chromR")}
  #  stopifnot(class(x) == "chromR")
  #  gt <- extract.gt(x, element = "GT", mask = x@var.info$mask)
  #  stopifnot(length(grep("(1/1|0/0|0/1)", unique(as.vector(gt)))) == 3)
  #  gt <- x@gt.m
  #
  
  if(deprecated == TRUE){
    msg <- "This function has been deprecated since vcfR 1.8.0."
    msg <- paste(msg, "If you would like to advocate to have this function included in future versions of vcfR please contact the maintainer.")
    msg <- paste(msg, "Contact information for package maintainers can be found with maintainer('vcfR').")
    stop(msg)
  }
  
  
  hwe <- function(x){
    # Genotype counts
    n11 <- x[1]
    n1i <- x[2]
    nii <- x[3]
    n   <- sum(n11, n1i, nii)
    #
    # Allele count and frequency
    n1 <- 2*n11 + n1i
    p1 <- n1/(2*n)
    #
    # Probability
    num <- (factorial(n) * factorial(n1) * factorial(2*n - n1) * 2^(n1i))
    den <- (factorial((n1 - n1i)/2) * factorial(n1i) * factorial(n-(n1+n1i)/2) * factorial(2*n))
    prob <- num/den
    #
    # Disequilibrium
    Da <- n11/n - (n1/(2*n))^2
    # Chi-square
    chisq <- ((n*Da)^2)/(n*p1^2) + ((-2*n*Da)^2)/(2*n*p1*(1-p1)) + ((n*Da)^2)/(n*(1-p1)^2)
    p <- 1 - stats::pchisq(chisq, df=1)
    return(c(prob, Da, chisq, p))
  }
  #  tmp[gt == "0/0"] <- 0
  #  tmp[gt == "0/1"] <- 1
  #  tmp[gt == "1/0"] <- 1
  #  tmp[gt == "1/1"] <- 2
#  gt <- extract.gt(x, element = "GT", mask = rep(TRUE, times=nrow(x@var.info)))
  gt <- extract.gt(x, element = "GT")
  tmp <- matrix(ncol=ncol(gt), nrow=nrow(gt))
  tmp[gt == "0/0" | gt == "0|0"] <- 0
  tmp[gt == "0/1" | gt == "0|1"] <- 1
  tmp[gt == "1/0" | gt == "1|0"] <- 1
  tmp[gt == "1/1" | gt == "1|1"] <- 2
  #
  gt <- tmp
  rm(tmp)
  #
  mask <- x@var.info$mask
  summ <- matrix(ncol=19, nrow=nrow(gt), 
                 dimnames=list(c(),
                               c('n', 'RR','RA','AA','nAllele','nREF','nALT','Ho','He',
                                 'hwe.prob', 'hwe.Da', 'hwe.chisq', 'hwe.p',
                                 'Ne','theta_pi','theta_w','theta_h','tajimas_d', 'faywu_h'))
  )
  #
  # Homozygous for reference allele.
  summ[mask,'RR'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==0, na.rm=TRUE)}))
  # Heterozygote.
  summ[mask,'RA'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==1, na.rm=TRUE)}))
  # Homozygous for alternate allele.
  summ[mask,'AA'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==2, na.rm=TRUE)}))
  #
  summ[mask, 'n'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(!is.na(x))}))
  #
  summ[mask,'nREF'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                    function(x){sum(2*length(stats::na.omit(x))-sum(x), na.rm=TRUE)})
  )
  summ[mask,'nALT'] <- rowSums(gt[mask, , drop=FALSE], na.rm=TRUE)
  summ[,'nAllele'] <- summ[,'nREF']+summ[,'nALT']
  #
  # Observed heterozygosity
  summ[mask,'Ho'] <- unlist(apply(gt[mask, , drop=FALSE], MARGIN=1,
                                  function(x){sum(x==1, na.rm=TRUE)/length(stats::na.omit(x))}))
  #
  # Expected heterozygosity
  summ[,'He'] <- 1 - ((summ[,'nREF']/summ[,'nAllele'])^2 + (summ[,'nALT']/summ[,'nAllele'])^2)
  #
  summ[,'Ne'] <- 1/(1-summ[,'He'])
  #
  # Hardy-Weinberg Disequilibrium
  summ[mask,c('hwe.prob', 'hwe.Da', 'hwe.chisq', 'hwe.p')] <- t(apply(summ[mask,c('RR', 'RA', 'AA'), drop=FALSE], MARGIN=1, FUN=hwe))
  #
  # Thetas.
  summ[,c('theta_pi','theta_w','theta_h')] <- t(apply(summ[,c('nREF','nALT'), drop=FALSE], MARGIN=1,thetas))
  #summ[,7:9] <- t(apply(summ[,c('nREF','nALT'), drop=FALSE], MARGIN=1,thetas))
  #
  summ[,'tajimas_d'] <- summ[,'theta_pi'] - summ[,'theta_w']
  summ[,'faywu_h'] <- summ[,'theta_pi'] - summ[,'theta_h']
  #  summ[,10] <- summ[,7] - summ[,8]
  #  summ[,11] <- summ[,7] - summ[,9]
  #
  #  print(head(summ))
  #  x@vcf.stat <- as.data.frame(summ)
  x@var.info <- cbind(x@var.info, as.data.frame(summ))
  return(x)
}



#' @title Population genetics summaries
#' @name Population genetics summaries
#' @rdname pop_gen_sum
#' @aliases gt.to.popsum
#' 
#' @description Functions that make population genetics summaries 
#' 
#' @param x object of class chromR or vcfR
#' 
#' @details 
#' This function creates common population genetic summaries from either a chromR or vcfR object.
#' The default is to return a matrix containing allele counts, He, and Ne.
#' \strong{Allele_counts} is the a comma delimited string of counts.
#' The first position is the count of reference alleles, the second positions is the count of the first alternate alleles, the third is the count of second alternate alleles, and so on.
#' \strong{He} is the gene diversity, or heterozygosity, of the population.
#' This is \eqn{1 - \sum x^{2}_{i}}, or the probability that two alleles sampled from the population are different, following Nei (1973).
#' \strong{Ne} is the effective number of alleles in the population.
#' This is \eqn{1/\sum x^{2}_{i}} or one minus the homozygosity, from Nei (1987) equation 8.17. 
#' 
#' Nei, M., 1973. Analysis of gene diversity in subdivided populations. Proceedings of the National Academy of Sciences, 70(12), pp.3321-3323.
#' 
#' Nei, M., 1987. Molecular evolutionary genetics. Columbia University Press.
#' 
#' @examples
#' data(vcfR_test)
#' # Check the genotypes.
#' extract.gt(vcfR_test)
#' # Summarize the genotypes.
#' gt.to.popsum(vcfR_test)
#' 
#' @export
gt.to.popsum <- function(x){
#  if(class(x) != "chromR" | class(x) != "vcfR"){stop("Object is not of class chromR or vcfR")}
  if(!inherits(x, c('chromR', 'vcfR'))){
    stop("Object is not of class chromR or vcfR")
  }
  
  if(class(x) == "chromR"){
    var.info <- x@var.info
    # If summaries already exist, we'll remove them.
    x@var.info <- x@var.info[,grep("^n$|^Allele_counts$|^He$|^Ne$", colnames(x@var.info), invert = TRUE)]
  }
  if(class(x) == "vcfR"){
    # var.info <- matrix(nrow = nrow(x@fix), ncol = 5)
    # colnames(var.info) <- c('mask', "n", "Allele_counts", "He", "Ne")
    # var.info <- matrix(nrow = nrow(x@fix), ncol = 1)
    # colnames(var.info) <- c('mask')
    # var.info[,'mask', drop = FALSE] <- TRUE
    var.info <- matrix(TRUE, ncol=1, nrow=nrow(x@fix), dimnames = list(NULL, 'mask'))
  }
  
  # Extract genotypes from vcf.gt
  gt <- extract.gt(x, element="GT")
  
  var.info <- .gt_to_popsum(var_info=var.info, gt=gt)

  if(class(x) == 'chromR'){
    x@var.info <- var.info
    return(x)
  } else {
    return(var.info[,-1])
  }
}

