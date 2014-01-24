\docType{class}
\name{Chrom-class}
\alias{Chrom-class}
\title{Chrom class}
\description{
A class for representing chromosomes (or contigs).
}
\details{
Defines a class for chromosomal or contig data.

This object has quite a few slots.

\itemize{ \item name name of the object \item seq object of
class DNAbin (ape) \item len length of sequence (integer)
\item vcf.meta vcf meta data \item vcf.fix vcf fixed data
\item vcf.gt vcf genotype data \item vcf.info data
extracted from the vcf.gt \item ann annotation data in a
gff-like data.frame

\item pop1 vector indicating members of pop1 \item pop2
vector indicating members of pop2

\item acgt matrix of nucleotide compositions \item n.w
matrix of \item windows matrix of windows \item nuccomp.w
data.frame of nucleotide composition windows \item snpden.w
data.frame of snp density windows

\item gt.m matrix of genotypes \item vcf.stat data.frame of
variant stats \item sfs matrix for the site frequency
spectrum \item link matrix for linkages

\item mask a logical vector to indicate masked variants }

More descriptions can be put here.
}

