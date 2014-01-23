#' Example data for vcfR.
#' 
#' An example dataset containing parts of the Phytophthora
#' infestans genome.
#' 
#' \itemize{
#'   \item pinf_dna   mitochondion IIa (GenBank: AY898627.1) as a DNAbin object (ape)
#'   \item pinf_gff   annotation data (gff-like) as a data.frame
#'   \item pinf_vcf   variant information as a vcfR object
#'
#' }
#' 
#' Note that it is encouraged to keep package contents small to facilitate easy
#' downloading and installation.  This is why a mitochondrion was chosen as an
#' example.  In practice I've used this package on supercontigs.  This package
#' was designed for much larger datasets in mind than in this example.
#' 
#' @examples
#' data(vcfR_example)
#' 
#' @docType data
#' @keywords datasets
#' @format A DNAbin object, a data.frame and a vcfR object
#' @name vcfR_example
#' @aliases pinf_dna pinf_vcf pinf_gff
NULL

