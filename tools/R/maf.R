
#' @title Minor allele frequency
#' @name maf
#' @rdname maf
#' 
#' @description
#' Calculate the minor (or other) allele frequency.
#'  
#' @param x an object of class vcfR or chromR
#' @param element specify the allele number to return
#' 
#' @details 
#' The function maf() calculates the counts and frequency for an allele.
#' A variant may contain more than two alleles.
#' Rare alleles may be true rare alleles or the result of genotyping error.
#' In an attempt to address these competing issues we sort the alleles by their frequency and the report statistics based on their position.
#' For example, setting element=1 would return information about the major (most common) allele.
#' Setting element=2 returns information about the second allele.
#' 
#' @return 
#' a matrix of four columns.
#' The first column is the total number of alleles, the second is the number of NA genotypes, the third is the count and fourth the frequency.
#' 
#' @export
maf <- function(x, element=2){

  get_maf <- function(x, element=2){
    maf <- vector(mode='numeric', length=4)
    names(maf) <- c('nAllele', 'NA', 'Count', 'Frequency')
    x <- unlist(strsplit(x, split="[|/]"))
    maf['NA'] <- sum( is.na(x) )
    x <- table(x, useNA = "no")
    x <- sort(x, decreasing = TRUE, na.last = TRUE)
#    if( nrow(x) == 0 ){
    if( length(x) == 0 ){
      is.na(maf) <- TRUE
    } else {
      maf['nAllele'] <- sum(x)
      if( !is.na(x[element]) ){
        maf['Count'] <- x[element]
        maf['Frequency'] <-  x[element]/sum(x)
      }
    }
    return(maf)
  }

  gt <- extract.gt(x, element = "GT")
  maf <- t(apply(gt, MARGIN=1, get_maf, element=element))
  return(maf)
}

