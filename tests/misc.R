
library(vcfR)

data(vcfR_example)

#head(pinf_vcf)


original_wd <- getwd()

setwd(tempdir())

write.vcf(pinf_vcf, "pinf_mt.vcf")
x <- read.vcf("pinf_mt.vcf")

x <- .Call('vcfR_readVcfBody', PACKAGE = 'vcfR', "pinf_mt.vcf")

x <- .Call('vcfR_readVcfBody2', PACKAGE = 'vcfR', "pinf_mt.vcf")

#head(x)

unlink("pinf_mt.vcf")

setwd(original_wd)





