#include <Rcpp.h>
// #include <string>
#include "vcfRCommon.h"


using namespace Rcpp;


int elementNumber(String x, std::string element = "GT"){
  //
  //  Determine the position of a query element
  //  in a colon delimited string.
  //
  //  x is a string similar to:
  //  GT:GQ:DP:RO:QR:AO:QA:GL
  //
  //  element is which element 
  //  we're looking for.

//  int eNum = 0;
  int start = 0;
//  int len = 0;
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


Rcpp::String extractElementS(Rcpp::String x, int position=0){
  //
  // x is a colon delimited string similar to:
  // GT:GQ:DP:RO:QR:AO:QA:GL
  //
  // position is the position in the colon delimited 
  // string which needs to be extracted.
  //

  int start = 0;
//  int pos = 1;
  int current_position = 1; // One-based so that 0 means did not observe.
  std::string istring = x; // Convert Rcpp::String to std::string
  istring.push_back(':');
  std::string teststring;
  int i;

  for(i=1; i <= istring.size(); i++){
    if(istring[i] == ':'){
//      if(pos == number){
//      Rcout << "Current position: " << current_position << ", Desired position: " << position << "\n";
//      Rcout << "Test string: " << istring.substr(start, i-start) << "\n";
      if(position == current_position){
        teststring = istring.substr(start, i-start);
        return teststring;
      } else {
        start = i+1;
        current_position++;
//        pos++;
        i++;
      }
    }
  }
  // If we get here we did not find the element.
  return NA_STRING;
}


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
  // If we get here we did not find the element.
  return(0);
}


// [[Rcpp::export]]
CharacterMatrix extract_GT_to_CM(DataFrame x, std::string element="DP") {
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
  // located in each row (varioant)
  for(i=0; i<column.size(); i++){
    positions[i] = elementNumber(column(i), element);
  }
  
  // Process the input DataFrame
  for(i = 1; i < x.size(); i++){ // Sample (column) counter
    column = x(i);
    for(j=0; j<column.size(); j++){ // Variant (row) counter
//      Rcout << column(j) << "\tposition: " << positions[j] << "\n";
      cm(j, i-1) = extractElementS(column(j), positions[j]);
//      Rcout << "Returned value: " << cm(j, i-1) << "\n";
    }
  }

  return cm;
}


// [[Rcpp::export]]
NumericMatrix CM_to_NM(CharacterMatrix x) {
  int i = 0;
  int j = 0;
  Rcpp::NumericMatrix nm(x.nrow(), x.ncol());  // NumericMatrix for output
  nm.attr("dimnames") = x.attr("dimnames");

  for(i=0; i<x.ncol(); i++){
    for(j=0; j<x.nrow(); j++){
//      Rcpp::String element = x(j, i);
//      Rcout << x(j,i) << " is NA? "  << CharacterMatrix::is_na(x(j, i)) << "\n";
//      if(CharacterMatrix::is_na(x(j, i)) == TRUE){
      if(x(j,i) == NA_STRING){
//      if(ISNA(x(j, i))){
        nm(j, i) = NA_REAL;
      } else {
        nm(j, i) = atof(x(j, i));
      }
    }
  }

  return nm;
}



// [[Rcpp::export]]
Rcpp::StringMatrix extract_haps(Rcpp::StringVector ref,
                                Rcpp::StringVector alt,
                                Rcpp::StringMatrix gt,
                                int vebosity) {
  // Vcf files are typically of one ploidy.
  
  std::string temp = Rcpp::as< std::string >(gt(1,1));
  int ploidy = 1;
  int i = 0;
  int j = 0;
  
  for(i=0; i<temp.length(); i++){
//    char allele_split = '|';
//    if(temp[i] == allele_split){ploidy++;}
    if(temp[i] == '|'){ploidy++;}
  }

  if(ploidy == 1){
    // Is either haploid or is not phased.
    return gt;
  }

  Rcpp::StringMatrix haps(gt.nrow(), gt.ncol() * ploidy);

  // Iterate over variants (rows of gt)
  for(i=0; i<gt.nrow(); i++){
    // Create a vector where the reference allele is at position 0
    // and alternate alleles are pushed on.
    std::vector<std::string> alleles;
    alleles[0] = ref(i);

    // The alternate alleles.
    std::vector < std::string > alt_vec;
    char alt_split = ','; // Must be single quotes!
    std::string line = Rcpp::as< std::string >(alt(i));
    vcfRCommon::strsplit(line, alt_vec, alt_split);
    for(j=0; j<alt_vec.size(); j++){
      alleles.push_back(alt_vec[i]);
    }

    // Process the genotypes (columns) into haplotypes.
    int hap_col = 0;
    for(j=0; j<gt.ncol(); j++){
      std::vector < std::string > al_vec;
      char al_split = '|'; // Must be single quotes!
      std::string line = Rcpp::as< std::string >(gt(i, j));
      vcfRCommon::strsplit(line, al_vec, al_split);
      int hap_num = 0;
      while(hap_num < ploidy){
        int al_num = stoi(al_vec[hap_num]);
        haps(1, hap_col) = alleles[al_num];
        hap_num++;
        hap_col++;
      }
    }
  }
  
  return(haps);
}
