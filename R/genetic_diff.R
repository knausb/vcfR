


calc_jost <- function(x){
  # x is a list of subpopulations.
  nPop <- length(x)
  
  HprimeS <- matrix(nrow = nrow(x[[1]]), ncol = nPop + 1)
  colnames(HprimeS) <- c(paste("Hprimej", names(x), sep = "_"), "HprimeS")

  HprimeT <- vector(mode = "numeric", length = nrow(HprimeS))
  
  for(i in 1:length(x)){
    ps <- as.character(x[[i]]$Allele_counts)
    ps <- strsplit(ps, split = ",")
    ps <- lapply(ps, function(x){as.numeric(x)/sum(as.numeric(x), na.rm = TRUE)})
    Hprimej <- lapply(ps, function(x){ 1 - sum(x^2, na.rm = TRUE) } )
    Hprimej <- unlist(Hprimej)
    HprimeS[,i] <- Hprimej
    
    ps <- lapply(ps, function(x){ (1/(nPop) * sum(x, na.rm = TRUE))^2 })
    HprimeT <- HprimeT + unlist(ps)
  }

  HprimeS[, ncol(HprimeS)] <- 1/nPop * rowSums(HprimeS[,1:length(x) - 1])  
  
  HprimeT <- 1 - HprimeT

  Ntilde <- lapply(x, function(x){ ncol(x)^-1 })
  Ntilde <- (sum(unlist(Ntilde))/nPop)^-1

  # Here 2 implies diploidy
  Hs_est <- ((2 * Ntilde)/( 2 * Ntilde - 1)) * HprimeS[, 'HprimeS']
  Ht_est <- HprimeT + Hs_est/(2 * Ntilde * nPop)
  Dest <- (Ht_est - Hs_est)/(1 - Hs_est) * (nPop/(nPop - 1))
  
  HprimeS <- cbind(HprimeS, HprimeT, Hs_est, Ht_est, Dest)
  return(HprimeS)
}

calc_nei <- function(x1, x2){
  # x1 is a data.frame for the total population.
  # x2 is a list of subpopulations.
  nPop <- length(x2)
  
  ps <- strsplit(as.character(x1$Allele_counts), split = ",")
  ps <- lapply(ps, function(x){as.numeric(x)/sum(as.numeric(x), na.rm = TRUE)})
  Ht <- unlist(lapply(ps , function(x){1- sum(x^2)}))
  
  Hs <- matrix(nrow = nrow(x2[[1]]), ncol = nPop)
  colnames(Hs) <- paste("Hs", names(x2), sep = "_")
  
  Htmax <- vector("character", length = nrow(Hs))
  
  for(i in 1:nPop){
    Htmax <- paste(Htmax, as.character(x2[[i]]$Allele_counts), sep = ",")
    ps <- strsplit(as.character(x2[[i]]$Allele_counts), split = ",")
    ps <- lapply(ps, function(x){as.numeric(x)/sum(as.numeric(x), na.rm = TRUE)})
    ps <- lapply(ps , function(x){1- sum(x^2)})
    Hs[,i] <- unlist(ps)
  }
  
  Htmax <- substring(Htmax, 2)
  ps <- strsplit(Htmax, split = ",")
  ps <- lapply(ps, function(x){as.numeric(x)/sum(as.numeric(x), na.rm = TRUE)})
  ps <- lapply(ps , function(x){1- sum(x^2)})
  Htmax <- unlist(ps)
  
  Gst <- (Ht - rowMeans(Hs))/Ht
  Gstmax <- (Htmax - rowMeans(Hs))/ Htmax
  Gprimest <- Gst/Gstmax
  
  Hs <- cbind(Hs, Ht, Gst, Htmax, Gstmax, Gprimest)
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
#' Calculate measures of genetic differentiation..
#'
#' @param vcf a vcfR object
#' @param pops factor indicating populations
#' @param method the method to measure differentiation
#'
#' @details 
#' Measures of genetic differentiation, or fixation indicies, are commonly reported population genetic parameters.
#' This function reports genetic differentiation for all variants presented to it.
#' 
#' 
#' The method \strong{nei} returns Nei's Gst as well as Hedrick's correction for high alleleism (Hedrick 2005).
#' This method is based on heterozygosity where all of the alleles in a population are used to calculate allele frequecies.
#' This may make this a good choice when there is a mixture of ploidies in the sample.
#' 
#' 
#' @references 
#' Hedrick, Philip W. "A standardized genetic differentiation measure." Evolution 59.8 (2005): 1633-1638.
#' 
#' 
#' @seealso 
#' poppr.amova in \href{https://cran.r-project.org/package=poppr}{poppr}, amova in \href{https://cran.r-project.org/package=ade4}{ade4}, amova in \href{https://cran.r-project.org/package=pegas}{pegas}, \href{https://cran.r-project.org/package=hierfstat}{hierfstat}, \href{https://cran.r-project.org/package=poppr}{hierfstat}, \href{https://cran.r-project.org/package=DEMEtics}{DEMEtics}, and , \href{https://cran.r-project.org/package=mmod}{mmod}.
#' 
#' 
#' @examples
#' data(vcfR_example)
#' myPops <- as.factor(rep(c('a','b'), each = 9))
#' myDiff <- genetic_diff(vcf, myPops, method = "nei")
#' colMeans(myDiff[,c(3:6,9)], na.rm = TRUE)
#' 
#'
genetic_diff <- function(vcf, pops, method = "nei"){
  
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
#  popTot <- .Call("vcfR_gt_to_popsum", PACKAGE = "vcfR", var_info = var_info, gt = gt)
  
  for(i in 1:nPop){
    subpop.l[[i]] <- .Call("vcfR_gt_to_popsum", 
                           PACKAGE = "vcfR", 
                           var_info = var_info, 
                           gt = gt[,pops == levels(pops)[i], drop = FALSE]
                           )
  }
  
  
  if( method == "nei" ){
    tot <- .Call("vcfR_gt_to_popsum", 
                 PACKAGE = "vcfR", 
                 var_info = var_info, 
                 gt = gt
    )
    
    gdiff <- calc_nei(tot, subpop.l)
  } else if( method == "jost" ){
    #gdiff <- calc_jost(subpop.l)
    warning('This methd is not currently implemented')
  }

  gdiff <- as.data.frame(gdiff)
  gdiff <- cbind(vcf@fix[,1:2, drop = FALSE], gdiff)
  return(gdiff)
}

