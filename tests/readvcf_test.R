

library(Rcpp)

cppFunction('int add(int x, int y, int z) {
  int sum = x + y + z;
  return sum;
}')


add(2, 3, 4)

##### ##### ##### ##### #####


library(Rcpp)

cppFunction('Rcpp::StringMatrix read_body_gz2(std::string x,
                                   Rcpp::NumericVector stats,
                                   int nrows = -1,
                                   int skip = 0,
                                   Rcpp::IntegerVector cols = 0,
                                   int verbose = 1) {

//  Rcpp::StringMatrix mymatrix(1,1);
//  mymatrix(0,0) = NA_STRING;
//  return mymatrix;

  cols = cols - 1; // R is 1-based
  int row_num = 0;

  if( nrows == -1 & skip == 0 ){
    nrows = stats[2];
  } else if ( nrows != -1 & skip == 0 ){
    // nrows = nrows;
  } else if ( nrows == -1 & skip > 0){
    nrows = stats[2] - skip;
  } else if ( nrows != -1 & skip > 0){
    // nrows = nrows;
  } else {
    Rcpp::Rcerr << "failed to calculate return matrix geometry.";
  }
  Rcpp::StringMatrix gt( nrows, cols.size() );

  return gt;
}')


library(vcfR)
data("vcfR_example")
test_dir <- tempdir()
ex_file <- paste(test_dir, "/test.vcf.gz", sep="")
write.vcf(vcf, file=ex_file)

stats <- .Call('vcfR_vcf_stats_gz', PACKAGE = 'vcfR', ex_file)

read_body_gz2(ex_file, stats, nrows = 10, skip = 0, cols = 1:stats['columns'], verbose = 1)
  
