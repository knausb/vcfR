\name{Proc_chrom}
\alias{Proc_chrom}
\alias{gt2popsum}
\alias{proc_chrom}
\alias{regex.win}
\alias{var.win}
\title{Process chrom object}
\usage{
proc_chrom(x, win.size = 1000, verbose = TRUE, ...)

regex.win(x, max.win = 1000, regex = "[acgtwsmkrybdhv]")

var.win(x, win.size = 1000)

gt2popsum(x)
}
\arguments{
  \item{x}{oject of class chrom}

  \item{win.size}{integer indicating size for windowing
  processes}

  \item{verbose}{logical indicating whether verbose output
  should be reported}

  \item{...}{arguments to be passed to methods}
}
\description{
Functions which process Chrom objects
}

