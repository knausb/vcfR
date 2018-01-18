#' Example data from the Variant Effect Predictor (VEP).
#' 
#' Example data to use with unit tests.
#' 
#' \itemize{
#'   \item vep vcfR object
#' }
#' 
#' 
#' Output from the \href{https://uswest.ensembl.org/info/docs/tools/vep/index.html}{VEP} may include values with multiple equals signs.
#' This does not appear to conform with the VCF specification (at the time of writing this \href{http://samtools.github.io/hts-specs/}{VCF v4.3}).
#' But it appears fairly easy to accomodate.
#' This example data can be used to make unit tests to validate functionality.
#' 
#' 
#' @examples
#' data(vep)
#' vcfR2tidy(vep, info_only = TRUE)$fix
#' 
#' 
#' @docType data
#' @keywords datasets
#' @format A vcfR object
#' @name vep
#' @aliases vep
NULL