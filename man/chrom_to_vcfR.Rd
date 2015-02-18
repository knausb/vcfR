\name{chrom_to_vcfR}
\alias{chrom_to_vcfR}
\title{Convert chrom objects to vcfR objects}
\usage{
chrom_to_vcfR(x, use.mask = TRUE)
}
\arguments{
  \item{x}{Object of class chrom}

  \item{use.mask}{Logical, determine if mask from chrom
  object should be used to subset vcf data}
}
\value{
Returns an object of class vcfR.
}
\description{
Convert chrom objects to vcfR objects.
}
\details{
The chrom object is subset and recast as a vcfR object.
When use.mask is set to TRUE (the default), the object is
subset to only the variants (rows) indicated to include by
the mask.  When use.mask is set to FALSE, all variants
(rows) from the chrom object are included in the new vcfR
object.
}

