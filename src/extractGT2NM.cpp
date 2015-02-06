// [[Rcpp::plugins(cpp11)]]
#include <Rcpp.h>
#include <string>

using namespace Rcpp;

// For more on using Rcpp click the Help button on the editor toolbar


int elementNumber(String x, std::string element = "DP"){
  //
  //  x is a string similar to:
  //  GT:GQ:DP:RO:QR:AO:QA:GL
  //
  //  element is which element 
  //  we're looking for.

  int eNum = 0;
  int start = 0;
  int len = 0;
  int pos = 1;
  int i;
  std::string istring = x;
  istring = istring + ":";
  std::string teststring;

  // Scroll through string.
  // Every time we find a delimiter (:)
  // check to see if we have our field.
  for(i=1; i <= istring.size(); i++){
    if(istring[i] == ':'){
      teststring = istring.substr(start, i-start);
      if(teststring == element){
        return pos;
      } else {
        pos = pos + 1;
        start = i+1;
        i++;
      }
    }
  }
  // If we get here then we did not observe our element.
  return 0;
}


// Don't really need this here.
// Perhaps for CharacterMatrix function?
std::string extractElementS(String x, int number=1){
  //
  // x is a string similar to:
  // GT:GQ:DP:RO:QR:AO:QA:GL
  //
  // number is the position in the colon delimited 
  // string which needs to be extracted.
  //
  int count = 0;
  int start = 0;
  int pos = 1;
  std::string istring = x;
  std::string teststring;
  
  for(int i=1; i <= istring.size(); i++){
    if(istring[i] == ':'){
//      Rcout << "Pos: " << pos << ",\tNumber: " << number << "\n";
      if(pos == number){
        teststring = istring.substr(start, i-start);
        return teststring;
      } else {
        start = i+1;
        pos++;
        i++;
      }
    }
  }
//  std::string ostring = istring.substr(3,2);
//  return "yup\n";
//  return ostring;
}


double extractElementD(String x, int number=1){
  //
  // x is a string similar to:
  // GT:GQ:DP:RO:QR:AO:QA:GL
  //
  // number is the position in the colon delimited 
  // string which needs to be extracted.
  //
  int count = 0;
  int start = 0;
  int pos = 1;
  std::string istring = x;
  std::string teststring;
  
  for(int i=1; i <= istring.size(); i++){
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
}


//' Extract numeric data from genotype field of VCF
//' 
//' @param x A dataframe
//' @param element A string matching the element to be retrieved
//' @export
// [[Rcpp::export]]
NumericMatrix extractGT2NM(DataFrame x, std::string element="DP") {
  int i = 0;
  int j = 0;
  
/*  
  Rcout << "Inside of extractGT2NV";
  Rcout << "\n";
  int dfsize = x.size();
  Rcout << "DF size: ";
  Rcout << dfsize;
  Rcout << "\n";
//  CharacterVector format = x[0];
  StringVector format = x[0];
  Rcout << "format 0: ";
  Rcout << format(0);
  Rcout << "\t";
  Rcout << format[1];
  Rcout << "\t";
  Rcout << format[1][0];
  Rcout << "\n";
  */

  NumericMatrix outM(x.nrows(), x.size()-1);
  StringVector format = x(0);

  std::vector<int> positions(format.size());
  
  for(i=0; i<format.size(); i++){
    positions[i] = elementNumber(format(i), element);
  }

  Rcpp::NumericVector outV(format.size());
  
  for(i=1; i<x.size(); i++){ // Sample counter
    StringVector sample = x(i);
    for(j=0; j<format.size(); j++){ // Variant counter
      outM(j, i-1) = extractElementD(sample(j), positions[j]);
    }
  }

  return outM;
}




