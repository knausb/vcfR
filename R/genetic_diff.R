

calc_jost <- function(x){
  # x is a list of subpopulations.
  # Each list element contains a data.frame with columns
  #  CHROM, POS, mask, n, Allele_counts, He, Ne.
  nPop <- length(x)
  nLoci <- nrow(x[[1]])
  
  # A matrix for heterozygosities.
  Hs <- matrix(nrow = nrow(x[[1]]), ncol = nPop)
  colnames(Hs) <- paste("Hs", names(x), sep = "_")
  
  # Find the maximum number of alleles.
  # We'll use this so we can store data in matrices.
  maxAlleles <- 0
  for(j in 1:nPop){
    tmp <- strsplit(as.character(x[[j]]$Allele_counts), split = ",")
    tmp <- max(unlist(lapply(tmp, function(x){length(x)})))
    if(tmp > maxAlleles){
      maxAlleles <- tmp
    }
  }

  # Hs is the heterozygosities for population j (created above).
  # subPop.l is a list of matricies that hold allele counts for each population.
  # Nj is the count or number of each allele in population j.
  subPop.l <- vector(mode = 'list', length = nPop)
  Nj <- matrix(nrow = nLoci, ncol = nPop)
  for(j in 1:nPop){
    subPop.l[[j]] <- matrix(0, nrow = nLoci, ncol = maxAlleles)
    ps <- strsplit(as.character(x[[j]]$Allele_counts), split = ",")
    lapply(as.list(1:nLoci), function(x){ subPop.l[[j]][x,1:length(ps[[x]])] <<- as.numeric(ps[[x]])})
    Nj[,j] <- rowSums(subPop.l[[j]])
     
     ps <- lapply(ps, function(x){as.numeric(x)/sum(as.numeric(x), na.rm = TRUE)})
     ps <- lapply(ps , function(x){1- sum(x^2)})
     Hs[,j] <- unlist(ps)
  #   Hs[,j] <- unlist( lapply(ps , function(x){1- sum(x^2)}) )
  }
   
  #
  Dg <- lapply(subPop.l, function(x){sweep(x, MARGIN = 1, STATS = rowSums(x, na.rm = TRUE), FUN = "/")})
  Dg <- Reduce('+', Dg)
#  Dg <- Reduce('+', lapply(subPop.l, function(x){x/sum(x)}))
  Dg <- Dg/nPop
  Dg <- Dg^2
  Dg <- 1/rowSums(Dg)
  
#  Hs2 <- Hs^2
#  Hs2[Hs2 < 0] <- 0
  Ha <- rowMeans(Hs, na.rm = TRUE)
  Da <- 1/(1-Ha)
  Db <- Dg/Da
  
  a <- matrix(0, nrow = nLoci, ncol = maxAlleles)
  b <- matrix(0, nrow = nLoci, ncol = maxAlleles)

  # Calculate a
  sum1 <- matrix(0, nrow = nLoci, ncol = maxAlleles)
  sum2 <- matrix(0, nrow = nLoci, ncol = maxAlleles)
  for(j in 1:nPop){
    tmp <- sweep(subPop.l[[j]], MARGIN = 1, STATS = Nj[,j], FUN = "/")
    sum1 <- sum1 + tmp
    sum2 <- sum2 + tmp^2
  }
  sum1 <- sum1^2
  a <- (sum1 - sum2)/(nPop-1)
  a <- rowSums(a)

  # Calculate b
  for(j in 1:nPop){
    tmp <- subPop.l[[j]] * (subPop.l[[j]] - 1)
    tmp[tmp<0] <- 0
    myDenom <- Nj[,j] * (Nj[,j] - 1)
    tmp <- sweep(tmp, MARGIN = 1, STATS = myDenom, FUN = "/")
    b <- b + tmp
  }
  b <- rowSums(b)
  
  Dest_Chao <- 1 - (a/b)

  myRet <- data.frame(Hs)
  myRet$a <- a
  myRet$b <- b
  myRet$Dest_Chao <- Dest_Chao
  myRet$Da <- Da
  myRet$Dg <- Dg
  myRet$Db <- Db
  return(myRet)
}


calc_nei <- function(x1, x2){
  # x1 is a data.frame for the total population.
  # x2 is a list of subpopulations.
  nPop <- length(x2)
  
  ps <- strsplit(as.character(x1$Allele_counts), split = ",")
  nAllele <- unlist(lapply(ps, function(x){ sum(as.numeric(x)) }))
  ps <- lapply(ps, function(x){as.numeric(x)/sum(as.numeric(x), na.rm = TRUE)})
  Ht <- unlist(lapply(ps , function(x){1- sum(x^2)}))

#  nAllele <- x1$n
  nAlleles <- matrix(nrow = length(nAllele), ncol = nPop)
  
  Hs <- matrix(nrow = nrow(x2[[1]]), ncol = nPop)
  colnames(Hs) <- paste("Hs", names(x2), sep = "_")
  
  Htmax <- vector("character", length = nrow(Hs))
  Hsize <- matrix(nrow=nrow(Hs), ncol=nPop)
  colnames(Hsize) <- paste("n", names(x2), sep = "_")
  
  for(i in 1:nPop){
    Htmax <- paste(Htmax, as.character(x2[[i]]$Allele_counts), sep = ",")
    ps <- strsplit(as.character(x2[[i]]$Allele_counts), split = ",")
    nAlleles[,i] <- unlist(lapply(ps, function(x){ sum(as.numeric(x)) }))
    Hsize[,i] <- unlist(lapply(ps, function(x){sum(as.numeric(x), na.rm = TRUE)}))
    ps <- lapply(ps, function(x){as.numeric(x)/sum(as.numeric(x), na.rm = TRUE)})
    ps <- lapply(ps , function(x){1- sum(x^2)})
    Hs[,i] <- unlist(ps)
#    nAlleles[,i] <- x2[[i]]$n
  }
  
  Htmax <- substring(Htmax, 2)
  ps <- strsplit(Htmax, split = ",")
  ps <- lapply(ps, function(x){as.numeric(x)/sum(as.numeric(x), na.rm = TRUE)})
  ps <- lapply(ps , function(x){1- sum(x^2)})
  Htmax <- unlist(ps)
  
#  Gst <- (Ht - rowMeans(Hs))/Ht
  Gst <- (Ht - rowSums(Hs * nAlleles)/nAllele)/Ht
  
#  Gstmax <- (Htmax - rowMeans(Hs))/ Htmax
  Gstmax <- (Htmax - rowSums(Hs * nAlleles)/nAllele)/ Htmax
  Gprimest <- Gst/Gstmax
  
  Hs <- cbind(Hs, Ht, Hsize, Gst, Htmax, Gstmax, Gprimest)
  return(Hs)
}



#' @title Genetic differentiation
#'
#' @name genetic_diff
#' @rdname genetic_diff
#' @aliases genetic_diff
#' @export
#' 
#' @description
#' Calculate measures of genetic differentiation.
#'
#' @param vcf a vcfR object
#' @param pops factor indicating populations
#' @param method the method to measure differentiation
#'
#' @details Measures of genetic differentiation, or fixation indicies, are commonly reported population genetic parameters.
#' This function reports genetic differentiation for all variants presented to it.
#' 
#' The method \strong{nei} returns Nei's Gst as well as Hedrick's G'st, a correction for high alleleism (Hedrick 2005).
#' Here it is calculated as in equation 2 from Hedrick (2005) with the exception that the heterozygosities are weighted by the number of alleles observed in each subpopulation.
#' This is similar to \code{hierfstat::pairwise.fst()} but by using the number of alleles instead of the number of individuals it avoids making an assumption about how many alleles are contributed by each individual.
#' G'st is calculated as in equation 4b from Hedrick (2005).
#' This method is based on heterozygosity where all of the alleles in a population are used to calculate allele frequecies.
#' This may make this a good choice when there is a mixture of ploidies in the sample.
#' 
#' The method \strong{jost} return's Jost's D as a measure of differentiation.
#' This is calculated as in equation 13 from Jost (2008).
#' Examples are available at Jost's website: \url{http://www.loujost.com}.
#' 
#' A nice review of Fst and some of its analogues can be found in Holsinger and Weir (2009).
#' 
#' @seealso poppr.amova in \href{https://cran.r-project.org/package=poppr}{poppr}, amova in \href{https://cran.r-project.org/package=ade4}{ade4}, amova in \href{https://cran.r-project.org/package=pegas}{pegas}, \href{https://cran.r-project.org/package=hierfstat}{hierfstat}, \href{https://cran.r-project.org/package=DEMEtics}{DEMEtics}, and, \href{https://cran.r-project.org/package=mmod}{mmod}.
#' 
#' 
#' @references 
#' Hedrick, Philip W. "A standardized genetic differentiation measure." Evolution 59.8 (2005): 1633-1638.
#' 
#' Holsinger, Kent E., and Bruce S. Weir. "Genetics in geographically structured populations: defining, estimating and interpreting FST." Nature Reviews Genetics 10.9 (2009): 639-650.
#' 
#' Jost, Lou. "GST and its relatives do not measure differentiation." Molecular ecology 17.18 (2008): 4015-4026.
#' 
#' Whitlock, Michael C. "G'ST and D do not replace FST." Molecular Ecology 20.6 (2011): 1083-1091.
#' 
#' 
#' @examples
#' data(vcfR_example)
#' myPops <- as.factor(rep(c('a','b'), each = 9))
#' myDiff <- genetic_diff(vcf, myPops, method = "nei")
#' colMeans(myDiff[,c(3:8,11)], na.rm = TRUE)
#' hist(myDiff$Gprimest, xlab = expression(italic("G'"["ST"])), 
#'      col='skyblue', breaks = seq(0, 1, by = 0.01))
#' 
#'
genetic_diff <- function(vcf, pops, method = "nei"){
  
  if( class(vcf) != "vcfR" ){
    stop( paste("Expecting an object of class vcfR, instead received:", class(vcf)) )
  }
  
  if( class(pops) != "factor" ){
    stop( paste("Expecting a factor, instead received:", class(pops)) )
  }
  
  method <- match.arg(method, choices = c('jost', 'nei'))
  
  nPop <- length(levels(pops))
  subpop.l <- vector('list', length = nPop)
  names(subpop.l) <- levels(pops)
  
  # Extract genotypes.
  gt <- extract.gt(vcf, element = "GT")

  # Assemble data for gt_to_popsum.
  var_info <- as.data.frame(vcf@fix[,1:2, drop = FALSE])
  if( is.null(var_info$mask) ){
    var_info$mask <- TRUE
  }
  
  # Get allele counts for total and subs.

  for(i in 1:nPop){
    subpop.l[[i]] <- .gt_to_popsum(var_info = var_info, 
                           gt = gt[,pops == levels(pops)[i], drop = FALSE]
                           )
  }
  
  if( method == "nei" ){
    tot <- .gt_to_popsum(var_info = var_info, gt = gt)
    
    gdiff <- calc_nei(tot, subpop.l)
  } else if( method == "jost" ){
    gdiff <- calc_jost(subpop.l)
#    warning('This methd is not currently implemented')
  }

  gdiff <- as.data.frame(gdiff)
  gdiff <- cbind(vcf@fix[,1:2, drop = FALSE], gdiff)
  return(gdiff)
}

