#' Test data for vcfR.
#' 
#' A test file containing a diversity of examples intended to test functionality.
#' 
#' \itemize{
#'   \item vcfR_test vcfR object
#' }
#' 
#' 
#' This data set began as the example (section 1.1) from The Variant Call Format Specification \href{http://samtools.github.io/hts-specs/}{VCFv4.3} .
#' This data consisted of 3 samples and 5 variants.
#' As I encounter examples that challenge the code in vcfR they can be added to this data set.
#' 
#' 
#' 
#' @examples
#' data(vcfR_test)
#' 
#' 
#' \dontrun{
#'   # When I add data it can be saved with this command.
#'   save(vcfR_test, file="data/vcfR_test.RData") 
#'   }
#' 
#' 
#' @docType data
#' @keywords datasets
#' @format A vcfR object
#' @name vcfR_test
#' @aliases vcf_test
NULL