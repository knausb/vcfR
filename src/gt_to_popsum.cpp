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
  
  Rcpp::LogicalVector   mask = var_info["mask"];
  Rcpp::IntegerVector   nsample(mask.size());
  Rcpp::CharacterVector alleles(mask.size());
  
  Rcpp::StringVector    allele_counts(mask.size());
  Rcpp::NumericVector   Ho(mask.size());
  Rcpp::NumericVector   Ne(mask.size());
  
  int i = 0;
  int j = 0;
  int k = 0;
  int cols = gt.ncol();
  
  for(i=0; i<mask.size(); i++){ // Iterate over variants (rows)
    if(mask[i] == TRUE){
//      int myints[] = {0,0};
      std::vector<int> myints (1,0);
//      std::array<int,5> myints;
      for(j=0; j < cols; j++){ // Iterate over samples (columns)
        if(gt(i, j) != NA_STRING){
          nsample[i]++;  // Increment sample count.

          // Count alleles.
          std::vector < int > intv = gtsplit(as<std::string>(gt(i, j)));
          while(myints.size() - 1 < intv[0]){myints.push_back(0);}
          myints[intv[0]]++;
          for(k=1; k<intv.size(); k++){
            while(myints.size() - 1 < intv[0]){myints.push_back(0);}
            myints[intv[k]]++;
          }
        }
      }
      
      Rcout << "Sample count: " << nsample[i] << "\n";

      int n;
      char buffer [50];
      n=sprintf (buffer, "%d", myints[0]);
//      for(j=1; j < sizeof(myints); j++){
      for(j=1; j < myints.size(); j++){
        n=sprintf (buffer, "%s,%d", buffer, myints[j]);
      }

      Rcout << "Allele_counts: " << myints[0] << "," << myints[1] << "\n";
      allele_counts[i] = buffer;      
    } else { // Missing variant (row=NA)
      nsample[i] = NA_INTEGER;        
    }
 }
  
//  return var_info;
//  return Rcpp::DataFrame::create(var_info, _["n"]=nsample);
  return Rcpp::DataFrame::create(var_info, _["n"]=nsample, _["Allele_counts"]=allele_counts);
}
