#include <Rcpp.h>
using namespace Rcpp;




Rcpp::CharacterVector gtsplit(std::string line, std::string delimiter = "\t"){
  // Based on:
  // http://stackoverflow.com/a/14266139
  std::string token;
  std::vector<std::string> tempv;
  
  size_t pos = 0;
  while ((pos = line.find(delimiter)) != std::string::npos) {
    token = line.substr(0, pos);
    tempv.push_back(token);
    line.erase(0, pos + delimiter.length());
  }
  tempv.push_back(line);
  
  Rcpp::CharacterVector charvec(tempv.size());
  for(int i=0; i<tempv.size(); i++){charvec[i] = tempv[i];}
  return charvec;
}



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
          std::vector<int> allele_cnt(2,0);
          nsample[i]++;
//          strsplit(gt(i, j), delimiter = "/")
        }
      }
    } else {
      nsample[i] = NA_INTEGER;
    }
  }
  
  
//  return var_info;
  return Rcpp::DataFrame::create(var_info, _["n"]=nsample);
}
