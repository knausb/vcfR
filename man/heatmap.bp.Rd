\name{heatmap-barplots}
\alias{heatmap-barplots}
\alias{heatmap.bp}
\title{Heatmap with barplots}
\usage{
heatmap.bp(x, cbarplot = TRUE, rbarplot = TRUE, legend = TRUE, ...)
}
\arguments{
  \item{x}{a numeric matrix}

  \item{rbarplot}{logical indicating whether the rows
  should be summarized with a barplot}

  \item{cbarplot}{logical indicating whether the columns
  should be summarized with a barplot}

  \item{legend}{logical indicating whether a legend should
  be plotted}

  \item{...}{additional arguments to be passed on}
}
\description{
Heatmap of a numeric matrix with barplots summarizing
columns and rows.
}
\details{
The function heatmap.bp creates a heatmap from a numeric
matrix with optional barplots to summarize the rows and
columns.
}
\examples{
library(vcfR)
data(vcfR_example)
pinf_mt <- create.chrom('pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)
pinf_mt <- masker(pinf_mt)
pinf_gq <- extract.gt(pinf_mt, element="GQ", as.matrix=TRUE)
heatmap.bp(pinf_gq)
}
\seealso{
heatmap, image, heatmap2 in library(gplots)
}

