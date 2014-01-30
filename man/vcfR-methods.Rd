\docType{methods}
\name{plot,vcfR-method}
\alias{plot,vcfR-method}
\alias{plot.vcfR}
\alias{read.vcf}
\alias{write.vcf}
\title{vcfR methods}
\usage{
\S4method{plot}{vcfR}(x, y, ...)

read.vcf(x)

write.vcf(xvcf, vfile, mask = logical(0))
}
\arguments{
  \item{...}{Arguments to be passed to methods}

  \item{x}{variant call format (vcf) file}

  \item{xvcf}{a vcfR object}

  \item{vfile}{an output filename}

  \item{mask}{logical vector indicating rows to use}

  \item{y}{the y coordinates of points in the plot,
  \emph{optional} if \code{x} is an appropriate structure.}
}
\description{
Reads in a vcf file and stores it in a vcf class.
}
\details{
Reads in a vcf file and stores it in a vcf class.  Once the
number of lines the meta information contains the data is
divided into three tables: meta data, fixed data and
genotype data.
}
\examples{
library(vcfR)
data(vcfR_example)
plot(pinf_vcf)
}

