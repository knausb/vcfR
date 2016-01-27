#' Example data for vcfR.
#' 
#' An example dataset containing parts of the *Phytophthora infestans* genome.
#' 
#' \itemize{
#'   \item dna DNAbin object
#'   \item gff gff format data.frame
#'   \item vcf vcfR object
#   \item pinf_dna   mitochondion IIa (GenBank: AY898627.1) as a DNAbin object (ape)
#   \item pinf_gff   annotation data (gff-like) as a data.frame
#   \item pinf_vcf   variant information as a vcfR object
#' }
#' 
#' 
#' This data is a subset of the pinfsc50 dataset.
#' It has been subset to positions between 500 and 600 kbp.
#' The coordinate systems of the vcf and gff file have been altered by subtracting 500,000.
#' This results in a 100 kbp section of supercontig_1.50 that has positional data ranging from 1 to 100 kbp. 
#' 
#' Note that it is encouraged to keep package contents small to facilitate easy
#' downloading and installation.  This is why a mitochondrion was chosen as an
#' example.  In practice I've used this package on supercontigs.  This package
#' was designed for much larger datasets in mind than in this example.
#' 
#' @examples
#' data(vcfR_example)
#' 
#' 
#' 
#' 
#' @docType data
#' @keywords datasets
#' @format A DNAbin object, a data.frame and a vcfR object
#' @name vcfR_example
# @aliases pinf_dna pinf_vcf pinf_gff
#' @aliases dna gff vcf
NULL

