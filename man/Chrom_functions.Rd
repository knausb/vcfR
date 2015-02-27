\name{Chrom functions}
\alias{Chrom functions}
\alias{masker}
\alias{proc_chrom}
\alias{variant_table}
\alias{window_table}
\title{Chrom_functions}
\usage{
masker(x, QUAL = 999, mindp = 0.25, maxdp = 0.75, minmq = 20,
  maxmq = 50, ...)

proc_chrom(x, verbose = TRUE, ...)

variant_table(x)

window_table(x)
}
\arguments{
  \item{QUAL}{minimum variant quality}

  \item{mindp}{minimum cumulative depth}

  \item{maxdp}{maximum cumulative depth}

  \item{minmq}{minimum mapping quality}

  \item{maxmq}{maximum mapping quality}

  \item{...}{arguments to be passed to methods}

  \item{verbose}{logical stating whether to produce verbose
  output}

  \item{x}{object of class Chrom}
}
\description{
Functions which act on Chrom objects
}
\details{
The function \strong{masker} creates a logical vector that
determines which variants are masked. This vector is stored
in the var.info slot of a Chrom object.

The function \strong{variant_table} creates a data.frame
containing information about variants.

The funciton \strong{window_table}
}

