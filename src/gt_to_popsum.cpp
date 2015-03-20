#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::DataFrame gt_to_popsum(Rcpp::DataFrame var_info, Rcpp::CharacterMatrix gt) {
  // Calculate popgen summaries for the sample.
  // var_info should contain columns named 'CHROM', 'POS', 'mask' and possibly others.
  
  Rcpp::LogicalVector mask = var_info["mask"];
  Rcpp::IntegerVector nsample(mask.size());
  Rcpp::CharacterVector alleles(mask.size());
  Rcpp::CharacterVector allele_counts(mask.size());
  Rcpp::NumericVector Ho(mask.size());
  Rcpp::NumericVector Ne(mask.size());
  
  int i = 0;
  int j = 0;
  int cols = gt.ncol();
  
  for(i=0; i<mask.size(); i++){
    if(mask[i] == TRUE){
      for(j=0; j < cols; j++){
        if(gt(i, j) != NA_STRING){
          nsample[i]++;
        }
      }
    } else {
      nsample[i] = NA_INTEGER;
    }
  }
  
  
//  return var_info;
  return Rcpp::DataFrame::create(var_info, _["n"]=nsample);
}
