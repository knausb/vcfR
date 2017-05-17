


#
library(testthat)
#detach(package:vcfR, unload=TRUE)
#
library(vcfR)


#
context("genetic_diff")

##### ##### ##### ##### #####
# Jost's D

test_that("Jost's example works",{
  data("vcfR_test")
  
  # Create VCF data.
  jost <- vcfR_test[1,]
  jost@gt <- matrix(nrow=1, ncol=221)
  jost@gt[1,1] <- "GT"
  jost@gt[,2:11] <- "0/1"
  jost@gt[,12:21] <- "2/3"
  jost@gt[,22:121] <- "2/3"
  jost@gt[,122:221] <- "4/5"
  colnames(jost@gt) <- c("FORMAT", paste("sample", 1:220, sep="_"))
  
  # Pop factor
  myPops <- rep("b", times=220)
  myPops[1:20] <- "a"
  myPops <- as.factor(myPops)
  
  # genetic_diff
  tmp <- genetic_diff(jost, myPops, method = "jost")
  
  expect_equal(trunc(1e2*tmp$a), 25)
  expect_equal(trunc(1e7*tmp$b), 4788895)
  expect_equal(trunc(1e7*tmp$Dest_Chao), 4779589)
  expect_equal(trunc(1e7*tmp$Db), 13333333)
})


##### ##### ##### ##### #####
# Nei's Gst


test_that("Nei's method works",{
  #  devtools::load_all(".")
  #debug("calc_nei")
  data("vcfR_test")
  vcfR_test@gt <- cbind(vcfR_test@gt, vcfR_test@gt[,2:4])
  myPops <- as.factor(rep(c('a','b'), each=3))
  
  tmp <- genetic_diff(vcfR_test, myPops, method = "nei")
  
  expect_equal(10*tmp$Ht[1], 5)
  expect_equal(tmp$n_a[1], 6)
  expect_equal(tmp$Gprimest[1], 0)
})


##### ##### ##### ##### #####
# P. rubi example


vcf <- new("vcfR"
    , meta = c("##fileformat=VCFv4.1", "##FILTER=<ID=LowQual,Description=\"Low quality\">", 
"##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">", 
"##FORMAT=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth (reads with MQ=255 or with bad mates are filtered)\">", 
"##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">", 
"##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">", 
"##FORMAT=<ID=PL,Number=G,Type=Integer,Description=\"Normalized, Phred-scaled likelihoods for genotypes as defined in the VCF specification\">", 
"##INFO=<ID=AC,Number=A,Type=Integer,Description=\"Allele count in genotypes, for each ALT allele, in the same order as listed\">", 
"##INFO=<ID=AF,Number=A,Type=Float,Description=\"Allele Frequency, for each ALT allele, in the same order as listed\">", 
"##INFO=<ID=AN,Number=1,Type=Integer,Description=\"Total number of alleles in called genotypes\">", 
"##INFO=<ID=BaseQRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt Vs. Ref base qualities\">", 
"##INFO=<ID=ClippingRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref number of hard clipped bases\">", 
"##INFO=<ID=DP,Number=1,Type=Integer,Description=\"Approximate read depth; some reads may have been filtered\">", 
"##INFO=<ID=DS,Number=0,Type=Flag,Description=\"Were any of the samples downsampled?\">", 
"##INFO=<ID=FS,Number=1,Type=Float,Description=\"Phred-scaled p-value using Fisher's exact test to detect strand bias\">", 
"##INFO=<ID=HaplotypeScore,Number=1,Type=Float,Description=\"Consistency of the site with at most two segregating haplotypes\">", 
"##INFO=<ID=InbreedingCoeff,Number=1,Type=Float,Description=\"Inbreeding coefficient as estimated from the genotype likelihoods per-sample when compared against the Hardy-Weinberg expectation\">", 
"##INFO=<ID=MLEAC,Number=A,Type=Integer,Description=\"Maximum likelihood expectation (MLE) for the allele counts (not necessarily the same as the AC), for each ALT allele, in the same order as listed\">", 
"##INFO=<ID=MLEAF,Number=A,Type=Float,Description=\"Maximum likelihood expectation (MLE) for the allele frequency (not necessarily the same as the AF), for each ALT allele, in the same order as listed\">", 
"##INFO=<ID=MQ,Number=1,Type=Float,Description=\"RMS Mapping Quality\">", 
"##INFO=<ID=MQ0,Number=1,Type=Integer,Description=\"Total Mapping Quality Zero Reads\">", 
"##INFO=<ID=MQRankSum,Number=1,Type=Float,Description=\"Z-score From Wilcoxon rank sum test of Alt vs. Ref read mapping qualities\">", 
"##INFO=<ID=QD,Number=1,Type=Float,Description=\"Variant Confidence/Quality by Depth\">", 
"##INFO=<ID=ReadPosRankSum,Number=1,Type=Float,Description=\"Z-score from Wilcoxon rank sum test of Alt vs. Ref read position bias\">", 
"##reference=file:///nfs1/BPP/Grunwald_Lab/home/tabimaj/GBS/barcoded/rubi/Pr4671.fa"
)
    , fix = structure(c("scaffold_1", "scaffold_10", "scaffold_10", "scaffold_10", 
"98978", "51595", "51602", "51611", NA, NA, NA, NA, "G", "A", 
"G", "G", "A", "AGAC", "A", "A", "36.56", "13818.27", "13827.33", 
"13827.36", "PASS", "PASS", "PASS", "PASS", "AC=1;AF=2.924e-03;AN=342;BaseQRankSum=-3.672;ClippingRankSum=0.574;DP=1911;FS=0.000;InbreedingCoeff=-0.0286;MLEAC=1;MLEAF=2.924e-03;MQ=43.91;MQ0=0;MQRankSum=0.622;QD=2.03;ReadPosRankSum=0.033", 
"AC=80;AF=0.229;AN=350;BaseQRankSum=6.043;ClippingRankSum=0.027;DP=1997;FS=0.000;InbreedingCoeff=0.8732;MLEAC=75;MLEAF=0.214;MQ=44.10;MQ0=0;MQRankSum=0.740;QD=31.78;ReadPosRankSum=-2.538", 
"AC=80;AF=0.229;AN=350;BaseQRankSum=9.395;ClippingRankSum=0.221;DP=1992;FS=0.000;InbreedingCoeff=0.8732;MLEAC=75;MLEAF=0.214;MQ=44.10;MQ0=0;MQRankSum=0.690;QD=27.36;ReadPosRankSum=0.349", 
"AC=80;AF=0.229;AN=350;BaseQRankSum=4.568;ClippingRankSum=-0.194;DP=1990;FS=0.000;InbreedingCoeff=0.8730;MLEAC=75;MLEAF=0.214;MQ=44.11;MQ0=0;MQRankSum=0.755;QD=28.46;ReadPosRankSum=0.031"
), .Dim = c(4L, 8L), .Dimnames = list(NULL, c("CHROM", "POS", 
"ID", "REF", "ALT", "QUAL", "FILTER", "INFO")))
    , gt = structure(c("GT:AD:DP:GQ:PL", "GT:AD:DP:GQ:PL", "GT:AD:DP:GQ:PL", 
"GT:AD:DP:GQ:PL", "0/1:19,0:19:57:0,57,531", "0/0:17,0:17:51:0,51,765", 
"0/0:17,0:17:51:0,51,765", "0/0:17,0:17:51:0,51,765", NA, "0/0:18,0:18:54:0,54,810", 
"0/0:18,0:18:54:0,54,810", "0/0:18,0:18:54:0,54,810", "0/0:6,0:6:18:0,18,148", 
"0/0:11,0:11:33:0,33,495", "0/0:11,0:11:33:0,33,495", "0/0:10,0:10:30:0,30,450", 
"0/0:12,0:12:36:0,36,424", "0/0:14,0:14:42:0,42,630", "0/0:14,0:14:42:0,42,630", 
"0/0:14,0:14:42:0,42,630", "0/0:5,0:5:15:0,15,126", "0/0:4,0:4:12:0,12,180", 
"0/0:4,0:4:12:0,12,180", "0/0:4,0:4:12:0,12,180", "0/0:25,0:25:75:0,75,798", 
"0/0:11,0:11:33:0,33,495", "0/0:11,0:11:33:0,33,495", "0/0:11,0:11:33:0,33,495", 
"0/0:25,0:25:75:0,75,843", "0/0:14,0:14:42:0,42,630", "0/0:14,0:14:42:0,42,630", 
"0/0:14,0:14:42:0,42,630", "0/0:17,0:17:51:0,51,535", "0/0:12,0:12:36:0,36,540", 
"0/0:12,0:12:36:0,36,540", "0/0:12,0:12:36:0,36,540", "0/0:4,0:4:12:0,12,145", 
"0/0:10,0:10:30:0,30,450", "0/0:10,0:10:30:0,30,450", "0/0:10,0:10:30:0,30,450", 
"0/0:22,1:23:57:0,57,496", "0/0:15,0:15:45:0,45,675", "0/0:15,0:15:45:0,45,675", 
"0/0:15,0:15:45:0,45,675", "0/0:13,0:13:39:0,39,321", "0/0:10,0:10:30:0,30,450", 
"0/0:10,0:10:30:0,30,450", "0/0:10,0:10:30:0,30,450", NA, "0/0:9,0:9:27:0,27,405", 
"0/0:9,0:9:27:0,27,405", "0/0:9,0:9:27:0,27,405", "0/0:6,1:7:11:0,11,144", 
"0/0:13,0:13:39:0,39,585", "0/0:13,0:13:39:0,39,585", "0/0:13,0:13:39:0,39,585", 
"0/0:14,0:14:42:0,42,368", "0/0:18,0:18:54:0,54,810", "0/0:18,0:18:54:0,54,810", 
"0/0:18,0:18:54:0,54,810", "0/0:10,0:10:30:0,30,262", "0/0:19,0:19:57:0,57,855", 
"0/0:19,0:19:57:0,57,855", "0/0:19,0:19:57:0,57,855", "0/0:11,0:11:33:0,33,377", 
"0/0:6,0:6:18:0,18,270", "0/0:6,0:6:18:0,18,270", "0/0:6,0:6:18:0,18,270", 
"0/0:20,0:20:60:0,60,661", "0/0:18,0:18:54:0,54,810", "0/0:18,0:18:54:0,54,810", 
"0/0:18,0:18:54:0,54,810", NA, NA, NA, NA, "0/0:7,0:7:21:0,21,207", 
"0/0:4,0:4:12:0,12,180", "0/0:4,0:4:12:0,12,180", "0/0:4,0:4:12:0,12,180", 
"0/0:10,0:10:30:0,30,275", "0/0:13,0:13:39:0,39,585", "0/0:13,0:13:39:0,39,585", 
"0/0:13,0:13:39:0,39,585", "0/0:9,0:9:27:0,27,310", "0/0:13,0:13:39:0,39,585", 
"0/0:13,0:13:39:0,39,585", "0/0:13,0:13:39:0,39,585", "0/0:10,0:10:30:0,30,264", 
"0/0:18,0:18:54:0,54,810", "0/0:18,0:18:54:0,54,810", "0/0:18,0:18:54:0,54,810", 
"0/0:8,0:8:24:0,24,198", "0/0:15,0:15:45:0,45,675", "0/0:15,0:15:45:0,45,675", 
"0/0:15,0:15:45:0,45,675", "0/0:8,0:8:24:0,24,198", "0/0:30,0:30:90:0,90,1350", 
"0/0:30,0:30:90:0,90,1350", "0/0:30,0:30:90:0,90,1350", "0/0:15,0:15:45:0,45,451", 
"0/0:39,0:39:99:0,117,1755", "0/0:39,0:39:99:0,117,1755", "0/0:39,0:39:99:0,117,1755", 
"0/0:4,0:4:12:0,12,84", "0/0:11,0:11:33:0,33,495", "0/0:11,0:11:33:0,33,495", 
"0/0:11,0:11:33:0,33,495", "0/0:9,0:9:27:0,27,277", "0/0:4,0:4:12:0,12,180", 
"0/0:4,0:4:12:0,12,180", "0/0:4,0:4:12:0,12,180", "0/0:10,0:10:30:0,30,232", 
"0/0:7,0:7:21:0,21,315", "0/0:7,0:7:21:0,21,315", "0/0:7,0:7:21:0,21,315", 
"0/0:16,0:16:48:0,48,425", "0/0:21,0:21:63:0,63,945", "0/0:21,0:21:63:0,63,945", 
"0/0:21,0:21:63:0,63,945", "0/0:13,0:13:39:0,39,412", "0/0:9,0:9:27:0,27,405", 
"0/0:9,0:9:27:0,27,405", "0/0:9,0:9:27:0,27,405", "0/0:5,0:5:15:0,15,184", 
"0/0:11,0:11:33:0,33,495", "0/0:11,0:11:33:0,33,495", "0/0:11,0:11:33:0,33,495", 
"0/0:27,0:27:81:0,81,751", "0/0:37,0:37:99:0,111,1665", "0/0:37,0:37:99:0,111,1665", 
"0/0:37,0:37:99:0,111,1665", "0/0:19,0:19:57:0,57,531", "0/0:27,0:27:81:0,81,1215", 
"0/0:27,0:27:81:0,81,1215", "0/0:27,0:27:81:0,81,1215", "0/0:5,0:5:15:0,15,184", 
"0/0:10,0:10:30:0,30,450", "0/0:10,0:10:30:0,30,450", "0/0:10,0:10:30:0,30,450", 
NA, NA, NA, NA, "0/0:12,0:12:36:0,36,402", "0/0:7,0:7:21:0,21,315", 
"0/0:7,0:7:21:0,21,315", "0/0:7,0:7:21:0,21,315", "0/0:5,0:5:15:0,15,192", 
NA, NA, NA, "0/0:23,0:23:69:0,69,715", "0/0:23,0:23:69:0,69,1035", 
"0/0:23,0:23:69:0,69,1035", "0/0:23,0:23:69:0,69,1035", "0/0:15,0:15:45:0,45,580", 
"0/0:4,0:4:12:0,12,180", "0/0:4,0:4:12:0,12,180", "0/0:4,0:4:12:0,12,180", 
"0/0:10,0:10:30:0,30,273", "0/0:4,0:4:12:0,12,180", "0/0:4,0:4:12:0,12,180", 
"0/0:4,0:4:12:0,12,180", "0/0:23,0:23:69:0,69,619", "0/0:10,0:10:30:0,30,450", 
"0/0:10,0:10:30:0,30,450", "0/0:10,0:10:30:0,30,450", "0/0:15,0:15:45:0,45,500", 
"0/0:8,0:8:24:0,24,360", "0/0:8,0:8:24:0,24,360", "0/0:8,0:8:24:0,24,360", 
"0/0:8,0:8:24:0,24,198", "0/0:11,0:11:33:0,33,495", "0/0:11,0:11:33:0,33,495", 
"0/0:11,0:11:33:0,33,495", "0/0:4,0:4:12:0,12,99", "0/0:17,0:17:51:0,51,765", 
"0/0:17,0:17:51:0,51,765", "0/0:17,0:17:51:0,51,765", "0/0:5,0:5:15:0,15,158", 
"0/0:6,0:6:18:0,18,270", "0/0:6,0:6:18:0,18,270", "0/0:6,0:6:18:0,18,270", 
"0/0:18,0:18:54:0,54,659", "0/0:24,0:24:72:0,72,1080", "0/0:24,0:24:72:0,72,1080", 
"0/0:24,0:24:72:0,72,1080", "0/0:13,0:13:39:0,39,410", "0/0:8,0:8:24:0,24,360", 
"0/0:8,0:8:24:0,24,360", "0/0:8,0:8:24:0,24,360", NA, "0/0:4,0:4:12:0,12,180", 
"0/0:4,0:4:12:0,12,180", "0/0:4,0:4:12:0,12,180", "0/0:5,0:5:15:0,15,124", 
"0/0:8,0:8:24:0,24,360", "0/0:8,0:8:24:0,24,360", "0/0:8,0:8:24:0,24,360", 
NA, "0/0:13,0:13:39:0,39,585", "0/0:13,0:13:39:0,39,585", "0/0:13,0:13:39:0,39,585", 
"0/0:47,0:47:99:0,140,1366", "0/0:23,0:23:69:0,69,1035", "0/0:23,0:23:69:0,69,1035", 
"0/0:23,0:23:69:0,69,1035", "0/0:7,0:7:21:0,21,266", "0/0:12,0:12:36:0,36,540", 
"0/0:12,0:12:36:0,36,540", "0/0:12,0:12:36:0,36,540", "0/0:5,0:5:15:0,15,193", 
"0/0:6,0:6:18:0,18,259", "0/0:5,0:5:18:0,18,259", "0/0:5,0:5:15:0,15,225", 
"0/0:25,0:25:74:0,74,701", "0/0:35,0:35:99:0,105,1575", "0/0:35,0:35:99:0,105,1575", 
"0/0:35,0:35:99:0,105,1575", "0/0:7,0:7:21:0,21,173", "0/0:7,0:7:21:0,21,315", 
"0/0:7,0:7:21:0,21,315", "0/0:7,0:7:21:0,21,315", "0/0:6,0:6:18:0,18,229", 
"0/0:9,0:9:27:0,27,405", "0/0:9,0:9:27:0,27,405", "0/0:9,0:9:27:0,27,405", 
"0/0:23,0:23:69:0,69,649", "0/0:22,0:22:66:0,66,990", "0/0:22,0:22:66:0,66,990", 
"0/0:22,0:22:66:0,66,990", "0/0:5,0:5:15:0,15,191", "0/0:4,0:4:12:0,12,180", 
"0/0:4,0:4:12:0,12,180", "0/0:4,0:4:12:0,12,180", "0/0:12,0:12:36:0,36,297", 
"0/0:17,0:17:51:0,51,765", "0/0:17,0:17:51:0,51,765", "0/0:17,0:17:51:0,51,765", 
"0/0:14,0:14:42:0,42,381", "0/0:12,0:12:36:0,36,540", "0/0:12,0:12:36:0,36,540", 
"0/0:12,0:12:36:0,36,540", "0/0:7,0:7:20:0,20,159", "0/0:5,0:5:15:0,15,225", 
"0/0:5,0:5:15:0,15,225", "0/0:5,0:5:15:0,15,225", "0/0:12,0:12:36:0,36,327", 
"0/0:15,0:15:45:0,45,675", "0/0:15,0:15:45:0,45,675", "0/0:15,0:15:45:0,45,675"
), .Dim = c(4L, 63L), .Dimnames = list(NULL, c("FORMAT", "4291", 
"4813", "4815", "4821", "4822", "4823", "4824", "4826", "4827", 
"4828", "4829", "4830", "4831", "4832", "4833", "4834", "4835", 
"4839", "4840", "4845", "4976", "4977", "4978", "4979", "4993", 
"5079", "5080", "5081", "5082", "5083", "5084", "5085", "5086", 
"5089", "5090", "5091", "5093", "5094", "5097", "5098", "5099", 
"5114", "5115", "5116", "5117", "5118", "5119", "5120", "5121", 
"5122", "5123", "5126", "5140", "5141", "5142", "5293", "5294", 
"5295", "5296", "5353", "5354", "5357")))
)
  
  


myPops <- structure(c(2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 
1L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 1L, 
1L, 1L, 1L, 1L, 1L, 2L, 2L, 2L, 1L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 
2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 2L, 1L, 
1L), .Label = c("OR", "WA"), class = "factor")


#myDiff <- genetic_diff(vcf, myPops, method = "nei")


#gt <- extract.gt(vcf)
#table(gt[1,])


test_that("GBS Jost's method works",{
#  devtools::load_all(".")
#
  myDiff <- genetic_diff(vcf, myPops, method = "jost")

  expect_equal(myDiff$Hs_OR[1], 0)
  expect_equal(trunc(1e10*myDiff$Hs_WA[1]), 289792388)
  expect_equal(trunc(1e10*myDiff$a[1]), 19705882352)
  expect_equal(trunc(1e10*myDiff$b[1]), 19705882352)
  expect_equal(myDiff$Dest_Chao[1], 0)
  
  expect_equal(trunc(1e10*myDiff$Da[1]), 10147026552)
  expect_equal(trunc(1e10*myDiff$Dg[1]), 10148140019)
  expect_equal(trunc(1e10*myDiff$Db[1]), 10001097333)
})


