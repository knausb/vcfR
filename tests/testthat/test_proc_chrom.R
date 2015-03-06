# detach(package:vcfR, unload=T)
library(vcfR)
context("proc_chrom functions")

data(vcfR_example)


pinf_mt <- create_chrom(name='pinf_mt', seq=pinf_dna, vcf=pinf_vcf, ann=pinf_gff)
#head(pinf_mt)
#pinf_mt
#names(pinf_mt)
#plot(pinf_mt)
#pinf_mt <- masker(pinf_mt)
#pinf_mt <- proc_chrom(pinf_mt, win.size=1000)



# proc_chrom2
#wins <- .Call('vcfR_window_init', PACKAGE = 'vcfR', window_size=1e3, max_bp=length(pinf_dna))
#win2 <- .Call('vcfR_windowize_fasta', PACKAGE = 'vcfR', wins=wins, seq=as.character(pinf_dna)[1,])
#pinf_mt2 <- proc_chrom2(pinf_mt, win.size=1000)



