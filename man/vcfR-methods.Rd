\name{read.vcf}
\alias{read.vcf}
\alias{write.vcf}
\title{vcfR methods}
\usage{
read.vcf(x)

write.vcf(xvcf, vfile)
}
\arguments{
  \item{x}{variant call format (vcf) file}

  \item{xvcf}{a vcfR object}

  \item{vfile}{an output filename}
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

