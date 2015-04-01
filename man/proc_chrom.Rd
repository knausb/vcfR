\name{Process chrom objects}
\alias{Process chrom objects}
\alias{gt2popsum}
\alias{gt_to_popsum}
\alias{proc_chrom}
\alias{proc_chrom2}
\alias{regex.win}
\alias{seq_to_rects}
\alias{var.win}
\title{Process chrom object}
\usage{
proc_chrom(x, win.size = 1000, verbose = TRUE)

proc_chrom2(x, win.size = 1000, verbose = TRUE)

regex.win(x, max.win = 1000, regex = "[acgtwsmkrybdhv]")

seq_to_rects(x, chars = "acgtwsmkrybdhv", lower = TRUE)

var.win(x, win.size = 1000)

gt2popsum(x)

gt_to_popsum(x)
}
\arguments{
  \item{x}{object of class Chrom}

  \item{win.size}{integer indicating size for windowing
  processes}

  \item{verbose}{logical indicating whether verbose output
  should be reported}

  \item{...}{arguments to be passed to methods}

  \item{max.win}{maximum window size}

  \item{regex}{a regular expression to indicate nucleotides
  to be searched for}

  \item{chars}{a vector of characters to be used as a
  database for inclusion in rectangles}

  \item{lower}{converts the sequence and database to lower
  case, making the search case insensitive}
}
\description{
Functions which process Chrom objects

Create representation of a sequence. Begining and end
points are determined for stretches of nucleotides.
Stretches are determined by querying each nucleotides in a
sequence to determine if it is represented in the database
of characters (chars).
}
\details{
The function \strong{proc_chrom()} calls helper functions
to process the data present in a Chrom object into
summaries statistics.

The function \strong{regex.win()} is used to generate
coordinates to define rectangles to represent regions of
the chromosome containing called nucleotides
(acgtwsmkrybdhv). It is then called a second time to
generate coordinates to define rectangles to represent
regions called as uncalled nucleotides (n, but not gaps).

The function \strong{gt2popsum} is called to create
summaries of the variant data.

The function \strong{var.win} is called to create
windowized summaries of the Chrom object.
}

