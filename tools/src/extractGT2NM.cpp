#include <Rcpp.h>
#include "vcfRCommon.h"
#include <string>
//#include <vector>



using namespace Rcpp;





double extractElementD(String x, int number=1){
  //
  // x is a string similar to:
  // GT:GQ:DP:RO:QR:AO:QA:GL
  //
  // number is the position in the colon delimited 
  // string which needs to be extracted.
  //
//  int count = 0;
  int start = 0;
  int pos = 1;
  std::string istring = x;
  std::string teststring;
  unsigned int i = 0;
  
  for(i=1; i <= istring.size(); i++){
    if(istring[i] == ':'){
      if(pos == number){
        teststring = istring.substr(start, i-start);
        double teststring2 = atof(teststring.c_str());
        return teststring2;
//        return std::stod(teststring);
      } else {
        start = i+1;
        pos++;
        i++;
      }
    }
  }
  // If we get here we did not find the element.
  return(0);
}


// [[Rcpp::export]]
Rcpp::CharacterMatrix extract_GT_to_CM(Rcpp::DataFrame x, std::string element="DP") {
  int i = 0;
  int j = 0;
  Rcpp::StringVector column = x(0);   // Vector to check out DataFrame columns to
  std::vector<int> positions(column.size());  // Vector to hold position data
  Rcpp::CharacterMatrix cm(column.size(), x.size() - 1);  // CharacterMatrix for output
  
  // Swap column names, minus the first, from x to cm
  Rcpp::StringVector colnames = x.attr("names");
  colnames.erase(0);
  cm.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(), colnames);
  
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













