

#' @title Convert vcfR to DNAbin
#' @name vcfR2DNAbin
#' 
#' @rdname vcfR2DNAbin
#' @aliases vcfR2DNAbin
#' 
#' @description 
#' Convert objects of class vcfR to objects of class ape::DNAbin
#' 
#' @param x an object of class chromR or vcfR
#' @param extract.indels logical, at present, the only option is TRUE
#' @param consensus logical, at present, the only option is TRUE
#' @param extract.haps logical specifying whether to separate each genotype into alleles based on a delimiting character
#' @param gt.split character to delimit alleles within genotypes
#' @param ref.seq reference sequence (DNAbin) for the region being converted
#' @param start.pos chromosomal position for the start of the ref.seq
#' @param verbose logical specifying whether to produce verbose output
#' 
#' @details
#' Objects of class \strong{DNAbin}, from the package ape, store nucleotide sequence information.
#' Typically, nucleotide sequence information contains all the nucleotides within a region, for example, a gene.
#' Because most sites are typically invariant, this results in a large amount of redundant data.
#' This is why files in the vcf format only contain information on variant sites, it results in a smaller file.
#' Nucleotide sequences can be generated which only contain variant sites.
#' However, some applications require the invariant sites.
#' For example, inference of phylogeny based on maximum likelihood or Bayesian methods requires invariant sites.
#' The function vcfR2DNAbin therefore includes a number of options in attempt to accomodate various scenarios.
#' 
#' 
#' The presence of indels (insertions or deletions)in a sequence typically presents a data analysis problem.
#' Mutation models typically do not accomodate this data well.
#' For now, the only option is for indels to be omitted from the conversion of vcfR to DNAbin objects.
#' The option \strong{extract.indels} was included to remind us of this, and to provide a placeholder in case we wish to address this in the future.
#' 
#' 
#' The \strong{ploidy} of the samples is inferred from the first non-missing genotype.
#' The option \code{gt.split} is used to split this genotype into alleles and these are counted.
#' Values for \code{gt.split} are typically '|' for phased data or '/' for unphased data.
#' Note that this option is an exact match and not used in a regular expression, as the 'sep' parameter in \code{\link{vcfR2genind}} is used.
#' All samples and all variants within each sample are assumed to be of the same ploid.
#' 
#' 
#' Conversion of \strong{haploid data} is fairly straight forward.
#' The options \code{consensus}, \code{extract.haps} and \code{gt.split} are not relevant here.
#' When vcfR2DNAbin encounters missing data in the vcf data (NA) it is coded as an ambiguous nucleotide (n) in the DNAbin object.
#' When no reference sequence is provided (option \code{ref.seq}), a DNAbin object consisting only of variant sites is created.
#' When a reference sequence and a starting position are provided the entire sequence, including invariant sites, is returned.
#' The reference sequence is used as a starting point and variable sitees are added to this.
#' Because the data in the vcfR object will be using a chromosomal coordinate system, we need to tell the function where on this chromosome the reference sequence begins.
#' 
#' 
#' Conversion of \strong{diploid data} presents a number of scenarios.
#' When the option \code{consensus} is TRUE, each genotype is split into two alleles using gt.split and the two alleles are converted into their IUPAC ambiguity code.
#' This results in one sequence for each diploid sample.
#' This may be an appropriate path when you have unphased data.
#' Note that functions called downstream of this choice may handle IUPAC ambiguity codes in unexpected manners.
#' When extract.haps is set to TRUE, each genotype is split into two alleles using gt.split.
#' These alleles are inserted into two sequences.
#' Thsi results in two sequences per diploid sample.
#' Note that this really only makes sense if you have phased data.
#' The options ref.seq and start.pos are used as in halpoid data.
#' 
#' 
#' 
#' Conversion of \strong{polyploid data} is currently not supported.
#' However, I have made some attempts at accomodating polyploid data.
#' If you have polyploid data and are interested in giving this a try, feel free.
#' But be prepared to scrutinize the output to make sure it appears reasonable.
#' 
#' 
#' Creation of DNAbin objects from large chromosomal regions may result in objects which occupy large amounts of memory.
#' If in doubt, begin by subsetting your data and the scale up to ensure you do not run out of memory.
#' 
#' 
#' 
#' 
#' @seealso 
#' \href{http://cran.r-project.org/package=ape}{ape}
#' 
#' @export
vcfR2DNAbin <- function( x, extract.indels = TRUE , consensus = FALSE,
                         extract.haps = TRUE, gt.split="|",
                         ref.seq = NULL, start.pos = NULL,
                         verbose = TRUE )
{
  # Sanitize input.
  if( class(x) == 'chromR' ){ x <- x@vcf }
  if( class(x) != 'vcfR' ){ stop( "Expecting an object of class chromR or vcfR" ) }
  if( consensus == TRUE & extract.haps == TRUE){
    stop("consensus and extract_haps both set to TRUE. These options are incompatible. A haplotype should not be ambiguous.")
  }
  if( !is.null(start.pos) & class(start.pos) == "character" ){
    start.pos <- as.integer(start.pos)
  }
  
  # Check and sanitize ref.seq.
  if( class(ref.seq) != 'DNAbin' & !is.null(ref.seq) ){
    stop( paste("expecting ref.seq to be of class DNAbin but it is of class", class(ref.seq)) )
  }
  if( is.list(ref.seq) ){
    ref.seq <- as.matrix(ref.seq)
  }
  
  # Extract indels.
  # Currently the only option is TRUE.
  if( extract.indels == TRUE ){
    x <- extract.indels(x)
  } else {
    stop("extract.indels == FALSE is not currently implemented.")
  }
  
  # Save POS in case we need it.
  # i.e., for inserting variants into a matrix.
  pos <- as.numeric(x@fix[,'POS'])
  
  # If we think we have variants we should extract them.
  # Our GT matrix may contain zero rows, all NA or data.
  
  # Check for zero rows.
  if( nrow(x@fix) == 0 ){
    # Create an empty matrix.
    x <- x@gt[ 0, -1 ]
  } else if( sum(!is.na(x@gt[,-1])) == 0 ){
    # Check for all NA.
    # Case of zero rows will sum to zero here.
    # Create an empty matrix.
    x <- x@gt[ 0, -1 ]    
  } else {
    # If x is still of class vcfR, we should process it.
#  if( class(x) == "vcfR" ){
#    first.gt <- x@gt[ ,-1 ][ !is.na(x@gt[,-1]) ][1]
    if( consensus == TRUE & extract.haps == FALSE ){
      x <- extract.gt( x, return.alleles = TRUE, allele.sep = gt.split )       
      x <- alleles2consensus( x, sep = gt.split )
    } else {
      x <- extract.haps( x, gt.split = gt.split, verbose = verbose )
    }
  }

  
  

  
  # The only reason to NOT extract haplotypes is if
  # consensus = TRUE
#  if( extract.haps == TRUE ){
#    x <- extract.haps( x, gt.split = gt.split, verbose = verbose )
#  }
  
  # Determine if extract.haps returned something successfully.
  # nrow of x could be because x had no variants,
  # or because gt.split was not properly specified.
#  if( nrow(x) == 0 ){
#    warning( paste('extract.haps returned nrow = 0, gt.split = ', gt.split, '. Does this seem appropriate?', sep="") )
#  }
  

  # Data could be haploid, diploid or higher ploid.
  
  
  
  
  # x should be a matrix of variants by here.
  
  # Return full sequence when ref.seq is not NULL
  if( is.null(ref.seq) == FALSE ){
    # Create a matrix of nucleotides.
    # The number of columns should match the number
    # of columns in x (i.e., number of haplotypes).
    # The number of rows should match the reference
    # sequence length.
    # The matrix will be initialized with the 
    # reference and will have no variants.
    variants <- x
    x <- matrix( as.character(ref.seq),
                     nrow = length(ref.seq),
                     ncol = length(colnames(x)),
                     byrow = FALSE
    )
    colnames(x) <- colnames(variants)

    # Populate with variants.
    # We need to subset the variant data to the region of interest.
    # First we remove variants above the region.
    # Then we remove variants below this region.
    # Then we rescale the region to be one-based.
    variants <- variants[ pos <= start.pos + dim(ref.seq)[2], ]
    pos <- pos[ pos <= start.pos + dim(ref.seq)[2] ]    
    variants <- variants[ pos >= start.pos, ]
    pos <- pos[ pos >= start.pos ]
    pos <- pos - start.pos + 1
    x[pos,] <- variants
  }
  
  # DNAbin characters must be lower case.
  # tolower requires dim(X) to be positive.
  if( nrow(x) > 0 ){
#    x <- apply( x, MARGIN=2, tolower )
    x <- tolower( x )
  }
  
  # Convert matrix to DNAbin
  x <- ape::as.DNAbin(t(x))
  
  return(x)
}

  
  
  
  
  
  
  
  
vcfR2DNAbin1 <- function( x, extract.indels = TRUE , consensus = TRUE,
                           extract.haps = FALSE, gt.split="|",
                           ref.seq = NULL, start.pos = NULL,
                           verbose = TRUE )
    
  {
  # Sanitize input.
  if( class(x) == 'chromR' )
  {
    x <- x@vcf
  }
  if( class(x) != 'vcfR' )
  {
    stop( "Expecting an object of class chromR or vcfR" )
  }
  if( consensus == TRUE & extract.haps == TRUE){
    stop("consensus and extract_haps both set to true. These options are incompatible. A haplotype should not be ambiguous.")
  }
  if( !is.null(start.pos) & class(start.pos) == "character" ){
    start.pos <- as.integer(start.pos)
  }
  
  # Check and sanitize ref.seq.
  if( class(ref.seq) != 'DNAbin' & !is.null(ref.seq) ){
    stop( paste("expecting ref.seq to be of class DNAbin but it is of class", class(ref.seq)) )
  }
  if( is.list(ref.seq) ){
    ref.seq <- as.matrix(ref.seq)
  }
  
  # Extract indels.
  # Currently the only option is TRUE.
  if( extract.indels == TRUE ){
    x <- extract.indels(x)
  } else {
    stop("extract.indels == FALSE is not currently implemented.")
  }
  
  # If we removed all variants, return.
#  if( nrow(x@fix) < 1 & is.null(ref.seq) ){
#    return( x )
#    return( NA )
#  }
  
  # Save POS in case we need it.
  # i.e., for inserting variants into a matrix.
  pos <- as.numeric(x@fix[,'POS'])
  
  if( extract.haps == FALSE ){
    x <- extract.gt( x, return.alleles = TRUE, allele.sep = gt.split )
  } else if( extract.haps == TRUE & nrow(x@gt) > 0 ){
#    x <- extract_haps(x, gt_split=gt.split, verbose = verbose)
    x <- extract.haps( x, gt.split = gt.split, verbose = verbose )
  } else if( extract.haps == TRUE & nrow(x@gt) == 0 ){
#    ref.seq <- as.matrix(ref.seq)
    out.seq <- matrix( as.character(ref.seq),
                       nrow = dim(ref.seq)[2],
                       ncol = ncol(x) * 2,
                       byrow = FALSE
    )
    colnames(out.seq) <- paste(rep(colnames(x@gt)[-1], each=2), rep(c(0,1), times=ncol(x@gt[,-1]), sep="_"))
    out.seq <- ape::as.DNAbin(t(out.seq))
    return( out.seq )
  } else {
    stop( "Invalid specification of extract.haps.\nShould be a logical." )
  }
  
  # If we have no variants (rows) return NA.
#  if( nrow(x) < 1 )
#  {
#    return( NA )
#  }

  # Strategies to convert genotypes (with a delimiter) to nucleotides.
  # If extract.haps was set to TRUE, then our data is effectively haploid now.
  ploid <- length(unlist( strsplit( x[!is.na(x)][1], split=gt.split, fixed=TRUE ) ))
  
#  if( verbose == TRUE ){
#    message( paste("Ploidy detected to be:", ploid) )
#  }
  
  if( ploid == 1 ){
    # Haploid case
    x[is.na(x)] <- 'n'
    if( nrow(x) > 1 ){
      x <- apply(x, MARGIN=2, tolower)
    } else if( nrow(x) == 1 ){
      x <- apply(x, MARGIN=2, tolower)
      x <- matrix(x, nrow=1, dimnames = list( NULL, names(x)))
    }
  } else if ( ploid == 2 & consensus == TRUE ){
    # Diploid case
    x <- alleles2consensus( x, sep = gt.split, NA_to_n = TRUE )
    x <- apply(x, MARGIN=2, tolower)
  } else if ( ploid == 2 & consensus == FALSE ){
    stop( "Ploidy is 2 but consensus is set to FALSE.\nGenotypes must be converted to a nucleotide.\nConsider setting consensus to TRUE." )
#    stop( "function only valid for diploid and haploid genotypes." )
  } else if ( ploid > 2 & extract_haps == FALSE ){
    stop( "Ploidy is greater than 2 but extract_haps is set to FALSE.\nGenotypes must be converted to a nucleotide.\nConsider setting extract_haps to TRUE." )
  }
  
  # Fill with reference sequence.
  if( !is.null(ref.seq) ){
    if( class(ref.seq) == "list" & inherits(ref.seq, "DNAbin") ){
      ref.seq <- as.matrix(ref.seq)
    }
    
    # Subset the vcf data to our region of interest.
    end.pos <- start.pos + dim(ref.seq)[2] - 1
    x <- x[ pos >= start.pos & pos <= end.pos, , drop = FALSE ]
    
    # Subset pos to our region of interest
    # and set to new coordinate system.
    pos <- pos[ pos >= start.pos & pos <= end.pos ]
    pos <- pos - start.pos + 1
    
    # Create a matrix with the same number of columns as vcf data
    # and the number of rows to match the sequence length.
    out.seq <- matrix( as.character(ref.seq),
                       nrow = dim(ref.seq)[2],
                       ncol = ncol(x),
                       byrow = FALSE
    )
    colnames(out.seq) <- colnames(x)
    
#    pos <- pos - start.pos + 1
    out.seq[pos,] <- x
    x <- out.seq
  }

  # DNAbin characters must be lower case.
  x <- ape::as.DNAbin(t(x))
  x
}




