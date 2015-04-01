\name{extract.gt.allR}
\alias{extract.gt}
\alias{extract.gt.allR}
\title{Extract elements from the GT section of a vcf format object}
\usage{
extract.gt.allR(x, element = "GT", mask = logical(0), as.matrix = FALSE)

extract.gt(x, element = "GT", mask = logical(0), as.numeric = FALSE)
}
\arguments{
  \item{x}{An object of class Chrom, vcfR or data.frame}

  \item{element}{element to extract from vcf genotype data.
  Common options include "DP", "GT" and "GQ"}

  \item{mask}{a logical vector indicating which variants
  (rows) to include}

  \item{as.matrix}{attempt to recast as a numeric matrix}

  \item{as.numeric}{Logical, should the matrix be converted
  to numerics}
}
\description{
Extract elements from the GT section of a vcf format object
}
\details{
Note that when 'as.numeric' is set to 'TRUE' but the data
are not actually numeric, unexpected results will likely
occur.
}

