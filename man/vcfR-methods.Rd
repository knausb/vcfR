\docType{methods}
\name{head,vcfR-method}
\alias{head,vcfR-method}
\alias{head.vcfR}
\alias{plot,vcfR-method}
\alias{plot.vcfR}
\alias{read.vcf}
\alias{subset,vcfR-method}
\alias{subset.vcf}
\alias{write.vcf}
\title{vcfR methods}
\usage{
\S4method{head}{vcfR}(x, y, n = 6, ...)

\S4method{plot}{vcfR}(x, y, ...)

\S4method{subset}{vcfR}(x, regex, ...)

read.vcf(x)

write.vcf(xvcf, vfile, mask = logical(0))
}
\arguments{
  \item{...}{Arguments to be passed to methods}

  \item{regex}{a regular expression to search for}

  \item{x}{variant call format (vcf) file}

  \item{xvcf}{a vcfR object}

  \item{vfile}{an output filename}

  \item{mask}{logical vector indicating rows to use}

  \item{n}{a single integer. If positive, size for the
  resulting object: number of elements for a vector
  (including lists), rows for a matrix or data frame or
  lines for a function. If negative, all but the \code{n}
  last/first number of elements of \code{x}.}

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
head(pinf_vcf)
plot(pinf_vcf)
}

