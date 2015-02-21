# extractgt devel

library(vcfR)

data(vcfR_example)

#head(pinf_vcf)


ncol(pinf_vcf@gt)

outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="GQ")
outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="GT")
outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="PL")
outm <- .Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="DP")


head(pinf_vcf@gt)
head(outm)



