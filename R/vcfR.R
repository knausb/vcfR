# vcf.R.

#' Variant call format files processed with vcfR.
#'
#' vcfR provides a suite of tools for input and output of variant call format (VCF) files, manipulation of their content and visualization.
#' 
#' @details
#' 
#' 
#' \strong{File input and output} is facilitated with the functions \code{\link{read.vcfR}} and \code{\link{write.vcf}}.
#' Input of vcf format data results in an S4 \code{\link{vcfR-class}} object.
#' Objects of class vcfR can be manipulated with \code{\link{vcfR-method}}s and \code{\link{extract.gt}}.
#' Contents of the vcfR object can be visualized with the \code{\link{plot.vcfR}} method.
#' More complex visualizations can be created using a series of functions.
#' See \code{vignette(topic="sequence_coverage")} for an example.
#' Once manipulations are complete the object may be written to a *.vcf.gz format file using \code{\link{write.vcf}} or exported to objects supported by other R packages with \code{\link{vcfR2genind}} or \code{\link{vcfR2loci}}.
#' 
#' 
#' More complex visualization can be accomplished by converting a vcfR object to a \code{\link{chromR-class}} object.
#' An example exists on the \code{\link{create.chromR}} man page.
#' 
#' 
#' 
#' A \strong{complete list of functions} can be displayed with: library(help = vcfR).
#'
#' \strong{Vignettes} (documentation) can be listed with: \code{browseVignettes('vcfR')}.
#' 
#' 
#' Several example \strong{datasets} are included in vcfR.
#' \strong{vcfR_test} comes from the VCF specification and provides a vcfR object with a diversity of examples in a small dataset.
#' \strong{vcfR_example} is a subset of the pinfsc50 dataset that includes VCF, GFF and FASTA data for moderate sized testing.
#' The \href{https://cran.r-project.org/package=pinfsc50}{pinfsc50} dataset is available as a separate package and includes VCF, GFF and FASTA data for testing and benchmarking.
#'
# @references Brian J Knaus (2015). Variant call format files processed 
# with vcfR. Journal TBA, NN(N),
#   N-NN. \url{http://www.someurl.org/v40/i03/}.
#' 
#' 
#' @import pinfsc50
#' @import ape
#' @docType package
#' @name vcfR
#' @rdname vcfR
#' @useDynLib vcfR
#' @importFrom Rcpp sourceCpp
#' @importFrom stats setNames
#' 
#' 
NULL



#### #### #### #### ####
# EOF
