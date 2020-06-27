# DNAbin and Annotation objects with more than one chromosome

# Example VCF data
system("cd data-raw; tabix -h ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20101123/interim_phase1_release/ALL.chr1.phase1.projectConsensus.genotypes.vcf.gz 1:10000-50000 > vcfex1.vcf")
system("cd data-raw; tabix -h ftp://ftp-trace.ncbi.nih.gov/1000genomes/ftp/release/20101123/interim_phase1_release/ALL.chr2.phase1.projectConsensus.genotypes.vcf.gz 2:10000-50000 > vcfex2.vcf")
#system("gunzip data-raw/vcfex1.vcf data-raw/vcfex2.vcf")
prepbcf <- function(filenm) {
  system(paste0("bcftools view ",filenm,".vcf -Oz -o ",filenm,".vcf.gz"))
  system(paste0("bcftools index ",filenm,".vcf.gz"))
}
prepbcf("data-raw/vcfex1")
prepbcf("data-raw/vcfex2")
system("bcftools merge --force-samples data-raw/vcfex1.vcf.gz data-raw/vcfex2.vcf.gz -Ov -o data-raw/merged.vcf")
vcf_multi <- vcfR::read.vcfR("data-raw/merged.vcf")
system("rm data-raw/vcfex* data-raw/*.tbi data-raw/*.vcf")
###

download.file("ftp://ftp.ensembl.org/pub/release-55/gtf/homo_sapiens/Homo_sapiens.GRCh37.55.gtf.gz", "data-raw/Homo_sapiens.GRCh37.55.gtf.gz")
gtfnames <- c("sequence", "source", "feature", "start", "end", "score", "strand","phase","attributes")
#tab <- read.table("data-raw/Homo_sapiens.GRCh37.55.gtf.gz",sep = "\t", col.names = gtfnames)
tab <- as.data.frame(vroom::vroom("data-raw/Homo_sapiens.GRCh37.55.gtf.gz",col_names = gtfnames))
tab_multi <- tab[tab$sequence %in% c("1","2") & (tab$start > 10000) & (tab$end < 50000),]
system("rm data-raw/Homo_sapiens.GRCh37.55.gtf.gz")
##

library(BSgenome.Hsapiens.UCSC.hg19)
seqs <- getSeq(
  BSgenome.Hsapiens.UCSC.hg19,
  names = c("chr1","chr2"),start = c(10001, 10001), end = c(50000, 50000)
)
seq_multi <- ape::as.DNAbin(seqs)

##

save(vcf_multi, tab_multi, seq_multi, file = "data/vcfR_example_extra.RData")
