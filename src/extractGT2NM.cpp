#include <Rcpp.h>
#include <string>
#include "vcfRCommon.h"
//#include <vector>

// Number of records to report progress at.
const int nreport = 1000;

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
      Rcpp::checkUserInterrupt();
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
      Rcpp::checkUserInterrupt();
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
                                char gt_split,
                                int verbose) {
  // Vcf files are typically of one ploidy.

  std::string temp = Rcpp::as< std::string >(gt(1,1));
  int ploidy = 1;
  int i = 0;
  int j = 0;
  int hap_col = 0;
  int hap_num = 0;

  for(i=0; i<temp.length(); i++){
//    if( temp[i] == '/' ){ploidy++;}
    if( temp[i] == gt_split ){ploidy++;}
  }
//  Rcout << "Ploidy = " << ploidy << "\n";

  if(ploidy == 1){
    // Is either haploid or is not phased.
//    return gt;
//    Rcout << "Ploidy = 1.\n";
    Rcpp::StringMatrix haps(1, 1);
    haps(0, 0) = NA_STRING;
    return haps;
  }

  // Initialize return structure
  Rcpp::StringMatrix haps(gt.nrow(), gt.ncol() * ploidy);

  Rcpp::List gt_names = gt.attr("dimnames");
  Rcpp::StringVector sample_names = gt_names(1);
  Rcpp::StringVector haplo_names(gt.ncol() * ploidy);
  
  j = 0;
  while(j < sample_names.size()){
    hap_num = 0;
    std::string sname = Rcpp::as< std::string >(sample_names[j]);
    while(hap_num < ploidy){
//      sname = sname + "_" + std::to_string(hap_num);
//      haplo_names(hap_col) = sname;
      haplo_names(hap_col) = sname + "_" + std::to_string(hap_num);
//      Rcout << haplo_names(hap_col) << "\n";
      hap_num++;
      hap_col++;
    }
    j++;
  }
  
  haps.attr("dimnames") = Rcpp::List::create(gt_names(0), haplo_names);
  
//  Rcpp::StringVector sample_names(gt.colnames());
//  Rcpp::StringVector haplo_names(gt.ncol() * ploidy);
//  int name_num = 0;



  // Iterate over variants (rows of gt)
  for(i=0; i<gt.nrow(); i++){

    // Split the alternate alleles string into alleles.
    std::vector < std::string > alleles_vec;
    char alleles_split = ','; // Must be single quotes!
    std::string line = Rcpp::as< std::string >(alt(i));
//    Rcout << i << " line: " << line << "\n";
    vcfRCommon::strsplit(line, alleles_vec, alleles_split);

    // Insert reference allele at begining of vector.
    std::string ref_allele = Rcpp::as< std::string >(ref(i));
    alleles_vec.insert(alleles_vec.begin(), ref_allele);

/*
    for(j=0; j<alleles_vec.size(); j++){
      Rcout << "alles_vec " << j << " value " << alleles_vec[j] << "\n";
    }
    Rcout << "\n";
    Rcout << "\n";
*/

    // Process the genotypes (columns) into haplotypes.
    hap_col = 0;
    for(j=0; j<gt.ncol(); j++){
      Rcpp::checkUserInterrupt();
      std::vector < std::string > al_vec;
      char al_split = gt_split; // Must be single quotes!
//      char al_split = '/'; // Must be single quotes!
      
      std::string line = Rcpp::as< std::string >(gt(i, j));
      vcfRCommon::strsplit(line, al_vec, al_split);
      hap_num = 0;
      while(hap_num < ploidy){
        int al_num = stoi(al_vec[hap_num]);
        haps(i, hap_col) = alleles_vec[al_num];
        hap_num++;
        hap_col++;
      }
    }
    if(i % nreport == 0 && verbose == 1){
      Rcout << "\rVariant " << i << " processed";
    }
  }
  if(verbose == 1){
    Rcout << "\rVariant " << i << " processed\n";
  }
  
  return(haps);
}
