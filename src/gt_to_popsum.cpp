#include <Rcpp.h>
using namespace Rcpp;




std::vector < int > gtsplit(std::string line){
//Rcpp::IntegerVector gtsplit(std::string line){
//  Rcpp::IntegerVector intv;
  std::vector < int > intv;
//  Rcout << "Genotype: " << line << "\n";

  // Case of a single digit.
  if(line.size() == 1){
    intv.push_back(atoi(line.c_str()));
  }
  
  int start=0;
  for(int i=1; i<line.size(); i++){
    if( line[i] == '/' || line[i] == '|' ){
      std::string temp = line.substr(start, i);
//      Rcout << "  i: " << i << ", temp: " << temp << "\n";
      intv.push_back(atoi(temp.c_str()));
      start = i+1;
      i = i+2;
    }
  }

  // Handle last element.
  std::string temp = line.substr(start, line.size());
//  Rcout << "  temp: " << temp << "\n";
  intv.push_back(atoi(temp.c_str()));
  
  return intv;
}



// [[Rcpp::export]]
Rcpp::DataFrame gt_to_popsum(Rcpp::DataFrame var_info, Rcpp::CharacterMatrix gt) {
  // Calculate popgen summaries for the sample.
  // var_info should contain columns named 'CHROM', 'POS', 'mask' and possibly others.
  
  Rcpp::LogicalVector mask = var_info["mask"];
  Rcpp::IntegerVector nsample(mask.size());
  Rcpp::CharacterVector alleles(mask.size());
  
  Rcpp::StringVector allele_counts(mask.size());
  Rcpp::NumericVector Ho(mask.size());
  Rcpp::NumericVector Ne(mask.size());
  
  int i = 0;
  int j = 0;
  int k = 0;
  int cols = gt.ncol();
  
  for(i=0; i<mask.size(); i++){ // Iterate over variants (rows)
    if(mask[i] == TRUE){
      for(j=0; j < cols; j++){ // Iterate over samples (columns)
        if(gt(i, j) != NA_STRING){
//          std::vector<int> allele_cnt(2,0);
          nsample[i]++;  // Increment sample count.
//          Rcpp::IntegerVector intv = gtsplit(as<std::string>(gt(i, j)));
          std::vector < int > intv = gtsplit(as<std::string>(gt(i, j)));
//          std::string allele_string = 
          char buffer [33];
//          itoa(intv[0], buffer, 10);
//          allele_counts(i) = intv[0];
          for(k=1; k<intv.size(); k++){
//            allele_counts(i) = strcat(allele_counts(i), ',', intv[k]);
//            allele_counts(i) = strcat(allele_counts(i), ',', intv[k]);
          }
          Rcout << "Allele_counts: " << allele_counts(i) << "\n";
//          string tmp_alleles = intv[0];
        }
      }
    } else { // Missing data (NA)
      nsample[i] = NA_INTEGER;
    }
  }
  
//  return var_info;
  return Rcpp::DataFrame::create(var_info, _["n"]=nsample);
}
