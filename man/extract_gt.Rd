\name{extract.gt}
\alias{extract.gt}
\alias{extract.gt2}
\title{Extract elements from the GT section of a vcf format object}
\usage{
extract.gt(x, element = "GT", mask = logical(0), as.matrix = FALSE)

extract.gt2(x, element = "GT", mask = logical(0), as.matrix = FALSE)
}
\arguments{
  \item{x}{An object of class chrom, vcfR or data.frame}

  \item{element}{element to extract from vcf genotype data.
  Common options include "DP", "GT" and "GQ"}

  \item{mask}{a logical vector indicating which variants
  (rows) to include}

  \item{as.matrix}{attempt to recast as a numeric matrix}
}
\description{
Extract elements from the GT section of a vcf format object
}

