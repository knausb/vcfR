\name{extract.gt.allR}
\alias{extract.gt}
\alias{extract.gt.allR}
\alias{extract_haps}
\alias{extract_indels}
\alias{extract_info}
\title{Extract elements from the GT section of a vcf format object}
\usage{
extract.gt.allR(x, element = "GT", mask = logical(0), as.matrix = FALSE)

extract.gt(x, element = "GT", mask = FALSE, as.numeric = FALSE)

extract_indels(x, return_indels = FALSE)

extract_info(x, element, as.numeric = FALSE, mask = FALSE)

extract_haps(x, mask = FALSE, verbose = TRUE)
}
\arguments{
  \item{x}{An object of class Chrom, vcfR or data.frame}

  \item{element}{element to extract from vcf genotype data.
  Common options include "DP", "GT" and "GQ"}

  \item{mask}{a logical indicating whether to apply the
  mask (TRUE) or return all variants (FALSE).
  Alternatively, a vector of logicals may be provided.}

  \item{as.matrix}{attempt to recast as a numeric matrix}

  \item{verbose}{should verbose output be generated}

  \item{as.numeric}{Logical, should the matrix be converted
  to numerics}

  \item{return_indels}{logical indicating whether to return
  indels or not}
}
\description{
Extract elements from the GT section of a vcf format object
}
\details{
Note that when 'as.numeric' is set to 'TRUE' but the data
are not actually numeric, unexpected results will likely
occur.
}

