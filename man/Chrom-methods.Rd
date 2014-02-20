\docType{methods}
\name{plot,Chrom-method}
\alias{ann2chrom}
\alias{chrom-methods}
\alias{chromo}
\alias{chromopop}
\alias{chromoqc}
\alias{create.chrom}
\alias{extract.gt}
\alias{head,Chrom-method}
\alias{head.chrom}
\alias{masker}
\alias{names,Chrom-method}
\alias{names.chrom}
\alias{plot,Chrom-method}
\alias{plot.chrom}
\alias{proc.chrom}
\alias{set.pop1}
\alias{set.pop2}
\alias{variant.table}
\alias{vcf2chrom}
\alias{windowize}
\title{Chrom methods}
\usage{
\S4method{plot}{Chrom}(x, y, ...)

\S4method{names}{Chrom}(x)

\S4method{head}{Chrom}(x)

vcf2chrom(x, y, ...)

ann2chrom(x, y, ...)

create.chrom(name, seq, vcf = NULL, ann = NULL)

set.pop1(x, pop1)

masker(x, QUAL = 999, mindp = 0.25, maxdp = 0.75, minmq = 0.25,
  maxmq = 0.75, ...)

windowize(x, win.size = 1000, max.win = 10000)

proc.chrom(x, pop1 = NA, pop2 = NA, verbose = TRUE)

chromoqc(x, nsum = FALSE, ...)

chromopop(x, ...)

chromo(x, verbose = TRUE, nsum = TRUE, DP = FALSE, QUAL = FALSE,
  MQ = FALSE, NE = FALSE, TPI = FALSE, TAJD = FALSE, FWH = FALSE,
  SNPDEN = FALSE, NUC = FALSE, ANN = FALSE, ...)

extract.gt(x, element = "GT", mask = logical(0), as.matrix = FALSE)

variant.table(x)
}
\arguments{
  \item{x}{an object of class chrom}

  \item{y}{an object}

  \item{...}{arguments}

  \item{name}{a name for the object}

  \item{seq}{a sequence as a DNAbin object}

  \item{vcf}{a vcfR object}

  \item{ann}{an annotation file (gff-like)}

  \item{pop1}{a numeric vector indicating the samples in
  population 1}

  \item{pop2}{a numeric vector indicating the samples in
  population 2}

  \item{mindp}{minimum cumulative depth}

  \item{maxdp}{maximum cumulative depth}

  \item{minmq}{minimum mapping quality}

  \item{maxmq}{maximum mapping quality}

  \item{win.size}{window size, in base pairs}

  \item{max.win}{maximum window size}

  \item{verbose}{logical stating whether to produce verbose
  output}

  \item{nsum}{logical for whether nsum will be displayed}

  \item{DP}{logical for whether cumulative depth will be
  displayed}

  \item{QUAL}{logical for whether variant quality will be
  displayed}

  \item{MQ}{logical for whether mapping quality will be
  displayed}

  \item{NE}{logical for whether effective size will be
  displayed}

  \item{TPI}{logical for whether Theta sub pi will be
  displayed}

  \item{TAJD}{logical for whether Tajima's D will be
  displayed}

  \item{FWH}{logical for whether Fay and Wu's H will be
  displayed}

  \item{SNPDEN}{logical for whether variant density will be
  displayed}

  \item{NUC}{logical for whether nucleotide content will be
  displayed}

  \item{ANN}{logical for whether annotation will be
  displayed}

  \item{element}{element to extract from vcf genotype data.
  Common options include "DP", "GT" and "GQ"}

  \item{mask}{a logical vector indicating which variants
  (rows) to include}

  \item{as.matrix}{attempt to recast as a numeric matrix}
}
\description{
Methods to work with objects of the chrom class
}
\details{
Reads in a vcf file and stores it in a Chrom class.

Reads in a vcf file and stores it in a Chrom class.

Creates and names a chrom object from a name, a sequence (a
DNA bin object), vcf data (a vcfR object) and annotation
data (gff-like).

Reads in a vcf file and stores it in a Chrom class.
}
\examples{
library(vcfR)
data(vcfR_example)
pinf_mt <- create.chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)
head(pinf_mt)
pinf_mt
names(pinf_mt)
plot(pinf_mt)
pinf_mt <- masker(pinf_mt)
pinf_mt <- proc.chrom(pinf_mt)
plot(pinf_mt)
chromoqc(pinf_mt)
chromopop(pinf_mt)
gt <- extract.gt(pinf_mt)
head(gt)
tab <- variant.table(pinf_mt)
head(tab)
hist(tab$Ho - tab$He, col=5)
# Note that this example is a mitochondrion, so this is a bit silly.
}

