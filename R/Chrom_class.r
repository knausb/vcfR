# chromR.

#### Class definition. ####

setOldClass("DNAbin")

#' @title Chrom class
#'
#' @name Chrom-class
#' @rdname Chrom-class
#'
#' @description
#' A class for representing chromosomes (or contigs).
#'
#' @details Defines a class for chromosomal or contig data.
#' 
#' This object has a numbr of slots.
#' 
#' \itemize{
#'   \item name name of the object  
#'   \item seq object of class DNAbin (ape)
#'   \item len length of sequence (integer)
#'   \item vcf.meta vcf meta data
#'   \item vcf.fix vcf fixed data
#'   \item vcf.gt vcf genotype data
#'   \item ann annotation data in a gff-like data.frame
#'
#'   \item var.info a data.frame containing information on variants
#'   \item win.info a data.frame containing information on windows
#'   \item seq.info a list containing information on the sequence
#'      
#   \item gt.m matrix of genotypes
#   \item sfs matrix for the site frequency spectrum
#   \item link matrix for linkages
#   
#   \item mask a logical vector to indicate masked variants
#' }
#' 
#' More descriptions can be put here.
#' 
#' The \strong{vcf.meta} slot is a list containing the meta information from the top of the vcf file.
#' 
#' The \strong{vcf.fix} slot is a data.frame containing the fixed data (the first eight columns) from the vcf file.
#' 
#' The \strong{vcf.gt} slot is a data.frame containing information about the samples.  The number of rows is the number of samples plus one, where the first row describes the format of the data in each cell.  The number of rows is equal to the numnber of variants.
#' 
#' The \strong{ann} slot is a data.frame containing gff format data.
#' 
#' The \strong{var.info} slot contains a data.frame containing information about variants
#' 
#' The \strong{win.info} slot contains a data.frame containing information about windows.  For example, window, start, end, length, A, C, G, T, N, other, variants and genic fields are stored here.
#' 
#' The \strong{seq.info} slot is a list containing two data.frames.  The first describes rectangles for called nucleotides and the second describes rectangles for 'N' calls.
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
  Class="Chrom",
  representation=representation(
    name = "character",
    seq = "DNAbin",
    len = "integer",
    vcf.meta = "character",
    vcf.fix = "data.frame",
    vcf.gt = "data.frame",
    ann = "data.frame",
    #
    var.info = "data.frame",
    win.info = "data.frame",
    seq.info = "list",
    #
    gt.m = "matrix",
    vcf.stat = "data.frame",
    sfs = "matrix",
    link = "matrix",
    #
    mask = "logical"
  ),
  prototype=prototype(
    vcf.fix = data.frame(matrix(ncol=8, nrow=0,
                                dimnames=list(c(),
                                              c('CHROM','POS','ID','REF','ALT','QUAL','FILTER','INFO'))
#                                              c('chrom','pos','id','ref','alt','qual','filter','info'))
    ),
    stringsAsFactors=FALSE),
    vcf.stat = data.frame(matrix(ncol=11, nrow=0, 
                                 dimnames=list(c(),
                                               c('Allele_num','R_num','A_num','Ho','He','Ne','theta_pi','theta_w','theta_h','tajimas_d','fw_h'))
    ),
    stringsAsFactors=FALSE),
    ann = data.frame(matrix(ncol=9, nrow=0,
                            dimnames=list(c(),
                                          c("seqid", "source", "type", "start", "end", "score", "strand", "phase", "attributes"))
    ),
    stringsAsFactors=FALSE)
  )
)



#### EOF ####
