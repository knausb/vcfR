
library(vcfR)

data(vcfR_example)

head(pinf_vcf)


working_dir <- getwd()

setwd(tempdir())

write.vcf(pinf_vcf, "pinf_mt.vcf")
x <- read.vcf("pinf_mt.vcf")

head(x)

unlink("pinf_mt.vcf")

setwd(working_dir)


