\name{read.vcf}
\alias{read.vcf}
\alias{read.vcf.devel}
\alias{write.vcf}
\title{Read and write vcf format files}
\usage{
read.vcf(file)

write.vcf(x, file = "", mask = logical(0), APPEND = FALSE)

read.vcf.devel(file)
}
\arguments{
  \item{file}{A filename for a variant call format (vcf)
  file}

<<<<<<< HEAD
  \item{xvcf}{A vcfR object}
=======
\item{x}{A vcfR object}
>>>>>>> 4c072b5db59e6b12ee92f3c0f3b6ee2fafaff731

  \item{mask}{logical vector indicating rows to use}

  \item{APPEND}{logical indicating whether to append to
  existing vcf file or write a new file}
}
\description{
Read and write files in the vcf format.
}
\details{
Reads in a vcf file and stores it in a vcf class.  Once the
number of lines the meta information contains the data is
divided into three tables: meta data, fixed data and
genotype data.

\code{read.vcf.devel} is a custom C++ implementation which
may be faster than read.vcf.
}
\examples{
library(vcfR)
data(vcfR_example)
head(pinf_vcf)
plot(pinf_vcf)
pinf_vcf[1:6,]
}
\seealso{
\code{\link[PopGenome]{readVCF}} \link[pegas]{read.vcf}

Bioconductor:
\href{www.bioconductor.org/packages/3.0/bioc/html/VariantAnnotation.html}{VariantAnnotation::readVcf}
}

