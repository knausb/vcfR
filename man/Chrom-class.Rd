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
\item vcf.gt vcf genotype data \item ann annotation data in
a gff-like data.frame

\item var.info a data.frame containing information on
variants \item win.info a data.frame containing information
on windows \item seq.info a list containing information on
the sequence

\item pop1 vector indicating members of pop1 \item pop2
vector indicating members of pop2

\item windows matrix of windows \item nuccomp.w data.frame
of nucleotide composition windows \item snpden.w data.frame
of snp density windows

\item gt.m matrix of genotypes \item sfs matrix for the
site frequency spectrum \item link matrix for linkages

}

More descriptions can be put here.
}

