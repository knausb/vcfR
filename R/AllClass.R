

#### Class definition. ####

#' @title vcfR class
#'
#' @name vcfR-class
#' @rdname vcfR-class
#'
#' @description
#' An S4 class for storing vcf data.
#' 
#' 
#' @slot meta character vector for the meta information
#' @slot fix  matrix for the fixed information
#' @slot gt   matrix for the genotype information 
#'
#'
#' @details Defines a class for variant call format data.
#' A vcfR object contains three slots.
#' The first slot is a character vector which holds the meta data.
#' The second slot holds an eight column matrix to hold the fixed data.
#' The third slot is a matrix which holds the genotype data.
#' 
#' See \code{vignette('vcf_data')} for more information.
#' See the \href{http://samtools.github.io/hts-specs/}{vcf specification} for the file specification.
#' 
#' @export
#' @import methods
setClass(
  Class="vcfR",
  representation=representation(
    meta="character",
    fix="matrix",
    gt="matrix"
  ),
  prototype=prototype(
    fix = matrix(ncol=8, nrow=0, 
                 dimnames=list(c(),
                               c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'))
                 )
  )
)




##### ##### ##### ##### #####
#
# Class chromR
#
##### ##### ##### ##### #####



#### Class definition. ####

setOldClass("DNAbin")

#' @title chromR class
#'
#' @name chromR-class
#' @rdname chromR-class
#'
#' @description
#' A class for representing chromosomes (or contigs).
#'
#' @details Defines a class for chromosomal or contig data.
#' 
#' This object has a number of slots.
#' 
#' \itemize{
#'   \item \strong{name} name of the object (character)
#'   \item \strong{len} length of the sequence (integer)
#'   \item \strong{window_size} window size for windowing analyses (integer)
#'   
#'   \item \strong{seq} object of class ape::DNAbin
#'   \item \strong{vcf} object of class vcfR
#'   \item \strong{ann} annotation data in a gff-like data.frame
#'
#'   \item \strong{var.info} a data.frame containing information on variants
#'   \item \strong{win.info} a data.frame containing information on windows
#'   \item \strong{seq.info} a list containing information on the sequence
#'      
#   \item gt.m matrix of genotypes
#   
#   \item mask a logical vector to indicate masked variants
#' }
#' 
# More descriptions can be put here.
#' 
#' The \strong{seq} slot contains an object of class ape::DNAbin.
#' A DNAbin object is typically either a matrix or list of DNAbin objects.
#' The matrix form appears to be better behaved than the list form.
#' Because of this behavior this slot should be the matrix form.
#' When this slot is not populated it is of class "NULL" instead of "DNAbin".
#' Note that characters need to be lower case when inserted into an object of class DNAbin.
#' The function \code{\link[base]{tolower}} can facilitate this.
#' 
#' 
#' The \strong{vcf} slot is an object of class vcfR \code{\link{vcfR-class}}.
#' 
#' The \strong{ann} slot is a data.frame containing \href{http://www.sequenceontology.org/gff3.shtml}{gff format} data.
#' When this slot is not populated it has nrows equal to zero.
#' 
#' The \strong{var.info} slot contains a data.frame containing information about variants.
#' Every row of this data.frame is a variant.
#' Columns will typically contain the chromosome name, the position of the variant (POS), the mask as well as any other per variant information.
#' 
#' The \strong{win.info} slot contains a data.frame containing information about windows.
#' For example, window, start, end, length, A, C, G, T, N, other, variants and genic fields are stored here.
#' 
#' The \strong{seq.info} slot is a list containing two matrices.
#' The first describes rectangles for called nucleotides and the second describes rectangles for 'N' calls.
#' Within each matrix, the first column indicates tha start position and the second column indicates the end position of each rectangle.
#' 
#' 
#' 
#' @seealso \code{\link{vcfR-class}}, \code{\link[ape]{DNAbin}},
#' \href{http://www.1000genomes.org/wiki/analysis/variant\%20call\%20format/vcf-variant-call-format-version-41}{vcf format}, 
#' \href{http://www.sequenceontology.org/gff3.shtml}{gff3 format}
#' 
#' 
#' 
#' @export
#' @import methods
setClass(
  Class="chromR",
  representation=representation(
    name = "character",
    len = "integer",
    window_size = "integer",
    vcf = "vcfR",
    seq = "DNAbin",
    ann = "data.frame",
    #
    var.info = "data.frame",
    win.info = "data.frame",
    seq.info = "list",
    #
    gt.m = "matrix"
  ),
  prototype=prototype(
    name = "Chromosome",
    len = as.integer(0),
    window_size = as.integer(1e3),
    vcf = new(Class="vcfR"),
    seq = as.DNAbin('n'),
    seq = NULL,
    ann = data.frame(matrix(ncol=9, nrow=0,
                            dimnames=list(c(),
                                          c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
    ),
    stringsAsFactors=FALSE)
  )
)





