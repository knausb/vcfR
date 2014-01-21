\name{read.vcf}
\alias{read.vcf}
\alias{write.vcf}
\title{vcf-methods}
\usage{
read.vcf(x)

write.vcf(x.vcf, vfile)
}
\arguments{
  \item{x}{variant call format (vcf) file}

  \item{x.vcf}{a vcf object}

  \item{vfile}{a character string naming a file

  #}
}
\description{
Reads in a vcf file and stores it in a vcf class.

#Takes a vcf object and writes it to file.
}
\details{
Reads in a vcf file and stores it in a vcf class.  Once the
number of lines the meta information contains the data is
divided into three tables: meta data, fixed data and
genotype data.

#

#Creates a file from the meta data and then uses cbind to
append the fixed and genotype data.

#
}

