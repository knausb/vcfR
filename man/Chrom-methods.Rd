\docType{methods}
\name{names,Chrom-method}
\alias{ann2chrom}
\alias{chrom-methods}
\alias{create.chrom}
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
\alias{variant_table}
\alias{vcf2chrom}
\alias{window_table}
\alias{windowize}
\title{Chrom methods}
\usage{
\S4method{names}{Chrom}(x)

\S4method{head}{Chrom}(x)

\S4method{plot}{Chrom}(x, y, ...)

vcf2chrom(x, y, ...)

ann2chrom(x, y, ...)

create.chrom(name, seq, vcf = NULL, ann = NULL, verbose = TRUE)

set.pop1(x, pop1)

set.pop2(x, pop2)

masker(x, QUAL = 999, mindp = 0.25, maxdp = 0.75, minmq = 20,
  maxmq = 50, ...)

windowize(x, win.size = 1000, max.win = 10000)

proc.chrom(x, verbose = TRUE, ...)

variant_table(x)

window_table(x)
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

  \item{QUAL}{variant quality}

  \item{mindp}{minimum cumulative depth}

  \item{maxdp}{maximum cumulative depth}

  \item{minmq}{minimum mapping quality}

  \item{maxmq}{maximum mapping quality}

  \item{win.size}{window size, in base pairs}

  \item{max.win}{maximum window size}

  \item{verbose}{logical stating whether to produce verbose
  output}
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
pinf_mt <- proc.chrom(pinf_mt, win.size=1000)
plot(pinf_mt)

chromoqc(pinf_mt)
chromoqc(pinf_mt, xlim=c(25e+03, 3e+04), dot.alpha=99)

set.seed(10)
x1 <- as.integer(runif(n=20, min=1, max=39000))
y1 <- runif(n=length(x1), min=1, max=100)
chromodot(pinf_mt, x1=x1, y1=y1)
chromodot(pinf_mt, x1=x1, y1=y1, label1='My data', x2=x1, y2=y1, label2='More of my data', dot.alpha='ff')

chromohwe(pinf_mt, dot.alpha='ff')

chromopop(pinf_mt)
gt <- extract.gt(pinf_mt)
head(gt)
tab <- variant_table(pinf_mt)
head(tab)
win <- window_table(pinf_mt)
head(win)
}

