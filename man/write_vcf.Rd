\name{write.vcf}
\alias{write.vcf}
\title{write.vcf

#}
\usage{
write.vcf(xvcf, vfile)
}
\arguments{
  \item{xvcf}{a vcf object}

  \item{vfile}{a character string naming a file

  #}
}
\description{
#Takes a vcf object and writes it to file.
}
\details{
#Creates a file from the meta data and then uses cbind to
append the fixed and genotype data.
}

