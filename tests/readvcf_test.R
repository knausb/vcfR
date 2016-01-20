

library(Rcpp)

cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}')


add(2, 3, 4)

##### ##### ##### ##### #####


library(Rcpp)

sourceCpp("tests/readvcf2.cpp")


library(vcfR)
data("vcfR_example")
test_dir <- tempdir()
ex_file <- paste(test_dir, "/test.vcf.gz", sep="")
write.vcf(vcf, file=ex_file)

stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)

read_body_gz2(ex_file, stats, nrows = -1, skip = 0, cols = 1:stats['columns'], verbose = 1)

read_body_gz2(ex_file, stats, nrows = 10, skip = 0, cols = 1:stats['columns'], verbose = 1)
read_body_gz2(ex_file, stats, nrows =  1, skip = 0, cols = 1:stats['columns'], verbose = 1)
read_body_gz2(ex_file, stats, nrows =  0, skip = 0, cols = 1:stats['columns'], verbose = 1)




