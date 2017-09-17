#include <Rcpp.h>
#include "vcfRCommon.h"
#include <string>
//#include <vector>



using namespace Rcpp;


/* 2017-09-15 This function does not appear to be used by anything. */
/* Should be deleted. */

//' @export
// [[Rcpp::export(name=".extract_GT_to_CM_B")]]
Rcpp::CharacterMatrix extract_GT_to_CM_B(Rcpp::DataFrame x, std::string element="DP", int depr = 1) {
  int i = 0;
  int j = 0;
  Rcpp::StringVector column = x(0);   // Vector to check out DataFrame columns to
  std::vector<int> positions(column.size());  // Vector to hold position data
  Rcpp::CharacterMatrix cm(column.size(), x.size() - 1);  // CharacterMatrix for output
  
  // Swap column names, minus the first, from x to cm
  Rcpp::StringVector colnames = x.attr("names");
  colnames.erase(0);
  cm.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(), colnames);
  
  if( depr == 1 ){
    Rcpp::Rcerr << "The function extract_GT_to_CM was deprecated in vcfR 1.6.0" << std::endl;
    Rcpp::Rcerr << "If you use this function and you would like to advocate its persistence, please contact the maintainer." << std::endl;
    Rcpp::Rcerr << "The maintainer of this package can be found with" << std::endl;
    Rcpp::Rcerr << "maintainer('vcfR')" << std::endl;
    Rcpp::stop("");
  }
  
  // Determine the position where the query element is 
  // located in each row (variant)
  for(i=0; i<column.size(); i++){
    positions[i] = elementNumber(column(i), element);
  }
  
  // Process the input DataFrame
  for(i = 1; i < x.size(); i++){ // Sample (column) counter
    column = x(i);
    for(j=0; j<column.size(); j++){ // Variant (row) counter
      Rcpp::checkUserInterrupt();
//      Rcout << column(j) << "\tposition: " << positions[j] << "\n";
      cm(j, i-1) = extractElementS(column(j), positions[j]);
//      Rcout << "Returned value: " << cm(j, i-1) << "\n";
    }
  }

  return cm;
}
















