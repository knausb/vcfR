# extractgt devel

library(vcfR)

data(vcfR_example)

#head(pinf_vcf)


.Call('vcfR_extract_GT_to_DF', PACKAGE = 'vcfR', pinf_vcf@gt, element="GQ")


