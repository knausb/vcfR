#' @title Test data for vcfR
#'
#' @description A test file containing a diversity of examples intended to test functionality.
#'
#' \itemize{
#'   \item vcfR_test: vcfR object
#'   \item vcfR_test_snpEff: vcfR object from VCF annotated by \href{http://pcingola.github.io/SnpEff/}{snpEff} (v5.1d)
#'   \item vcfR_test_VEP: vcfR object from VCF annotated by \href{https://www.ensembl.org/info/docs/tools/vep/index.html}{VEP} (v106)
#' }
#'
#' This data set began as the example (section 1.1) from The Variant Call Format Specification \href{http://samtools.github.io/hts-specs/}{VCFv4.3} .
#' This data consisted of 3 samples and 5 variants.
#' As I encounter examples that challenge the code in vcfR they can be added to this data set.
#'
#' @examples
#' data(vcfR_test)
#' data(vcfR_test_snpEff)
#' data(vcfR_test_VEP)
#'
#' \dontrun{
#'   # When I add data it can be saved with this command.
#'   save(vcfR_test, file="data/vcfR_test.RData")
#' }
#'
#' @docType data
#' @keywords datasets
#' @format A vcfR object
#' @name vcfR_test
#' @rdname vcfR_test
#' @aliases vcf_test
NULL

#' @name vcfR_test_snpEff
#' @rdname vcfR_test
NULL

#' @name vcfR_test_VEP
#' @rdname vcfR_test
NULL
