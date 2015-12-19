

context("chromoqc functions")

library(vcfR)
data("vcfR_example")

#chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
#chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, verbose=FALSE)
#chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
#chrom <- masker(chrom, min_QUAL=0, min_DP=350, max_DP=650, min_MQ=59.5, max_MQ=60.5)
#chrom <- proc.chromR(chrom, verbose = FALSE)


##### ##### ##### ##### #####

# Base: only dot plots.
chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
chromoqc( chrom )
chromoqc( chrom, boxp = FALSE )

# Dot plots with variants.
chrom <- create.chromR(name="Supercontig", vcf=vcf, verbose=FALSE)
chrom <- proc.chromR(chrom, verbose = FALSE)
chromoqc( chrom )
chromoqc( chrom, boxp = FALSE )

# Dot plots with nucleotides
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, verbose=FALSE)
chromoqc( chrom )
chromoqc( chrom, boxp = FALSE )

# Dot plots with nucleotides, processed
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, verbose=FALSE)
chrom <- proc.chromR(chrom, verbose = FALSE)
chromoqc( chrom )
chromoqc( chrom, boxp = FALSE )

# Dot plots with annotations
chrom <- create.chromR(name="Supercontig", vcf=vcf, ann=gff, verbose=FALSE)
chromoqc( chrom )
chromoqc( chrom, boxp = FALSE )

# Dot plots with annotations, processed
chrom <- create.chromR(name="Supercontig", vcf=vcf, ann=gff, verbose=FALSE)
chrom <- proc.chromR(chrom, verbose = FALSE)
chromoqc( chrom )
chromoqc( chrom, boxp = FALSE )

# Dot plots with nucleotides and annotations
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chromoqc( chrom )
chromoqc( chrom, boxp = FALSE )

# Dot plots with nucleotides and annotations, processed
chrom <- create.chromR(name="Supercontig", vcf=vcf, seq=dna, ann=gff, verbose=FALSE)
chrom <- proc.chromR(chrom, verbose = FALSE)
chromoqc( chrom )
chromoqc( chrom, boxp = FALSE )


##### ##### ##### ##### #####
# EOF.