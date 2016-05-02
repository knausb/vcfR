#include <Rcpp.h>
#include "vcfRCommon.h"
#include <string>
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


Rcpp::String extractElementS(Rcpp::String x, int position=0, int extract=1){
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
        if(extract == 1){
          teststring = istring.substr(start, i-start);
          return teststring;
        } else {
//          Rcpp::Rcout << "istring: " << istring << ", position: " << position << "\n";
          if(position == 1){
//            Rcpp::Rcout << "i: " << i << ", istring: " << istring << "\n";
            teststring = istring.substr( i + 1, istring.size() - i );
//            Rcpp::Rcout << "  teststring: " << teststring << "\n";
            teststring = teststring.substr(0, teststring.size() - 1); // Remove terminating : added above.
            return teststring;
          } else {
            teststring = istring.substr(0, start) + istring.substr( i + 1, istring.size());
            teststring = teststring.substr(0, teststring.size() - 1); // Remove terminating : added above.
            return teststring;
          }
        }
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


std::vector < std::string > get_allele_vector( Rcpp::String ref,
                                               Rcpp::String alt )
{
  std::string allele_string = alt;
  std::vector < std::string > allele_vector;
  char alleles_split = ','; // Must be single quotes!
  vcfRCommon::strsplit(allele_string, allele_vector, alleles_split);
  std::string ref2 = ref;
  allele_vector.insert( allele_vector.begin(), ref2 );

  return allele_vector;
}


std::string gt2alleles( Rcpp::String gt, std::vector< std::string > allele_vector, char allele_sep )
{
  // gt is a genotpye (e.g., 0/0, 1/2, 1/.).
  // allele_vector is a concatenation of the REF and ALT data.
  // allele_sep specifies the allele delimiter.
  
  // Recast allele_sep (char) to sep (std::string).
  std::string sep = std::string(1, allele_sep);
  
  // Create a std::string NA character.
  std::string na_allele = ".";

  // Split the genotype into alleles delimited by allele_sep.
  // This results in the vector gt_vector.
  std::string gt2 = gt;
  std::vector < std::string > gt_vector;
  vcfRCommon::strsplit( gt2, gt_vector, allele_sep );
  
//  Rcpp::Rcout << "Made it.!\n";
//  Rcpp::Rcout << "  gt_vector[0]" << gt_vector[0] << "\n";

  // stoi is not supported in MinGW so should not be used.
  // Instead we use std:istringstream to set (>>) the int value.
  
  //  int allele_number = std::stoi( gt_vector[0] );
  int allele_number; // Initialize allele counter.


  // Initialize our return value, gt3.
  std::string gt3;
  if( gt_vector[0].compare( na_allele ) == 0 ){
    gt3.append( na_allele );
  } else if ( ! (std::istringstream( gt_vector[0] ) >> allele_number) ){
    gt3.append( na_allele );
  } else {
    std::istringstream(gt_vector[0]) >> allele_number;
    gt3.append( allele_vector[ allele_number ] );
  }
  
  // If we're greater than haploid.
  if( gt_vector.size() > 1 )
  {
    for( int i=1; i<gt_vector.size(); i++ )
    {
//      if ( ! (std::istringstream(gt_vector[i]) >> allele_number) ) allele_number = 0;
//      gt3.append( sep );
//      gt3.append( allele_vector[ allele_number ] );
      if( gt_vector[i].compare( na_allele ) == 0 ){
        gt3.append( sep );
        gt3.append( na_allele );
      } else if ( ! (std::istringstream( gt_vector[i] ) >> allele_number) ){
        //  
        Rcpp::Rcout << "Couldn't convert string to int!\n";
        gt3.append( sep );
        gt3.append( na_allele );
      } else {
        gt3.append( sep );
        std::istringstream(gt_vector[i]) >> allele_number;
        gt3.append( allele_vector[ allele_number ] );
      }
    }
  }
    
  return gt3;
}



// [[Rcpp::export]]
Rcpp::StringMatrix extract_GT_to_CM2( Rcpp::StringMatrix fix,
                                         Rcpp::StringMatrix gt,
                                         std::string element="DP",
                                         char allele_sep = '/',
                                         int alleles = 0,
                                         int extract = 1 ) {
  int i = 0;
  int j = 0;

  // Initialize a return matrix.
  // The first column of gt is FORMAT, 
  // so the return_matrix will have one less column than in gt.
  // We'll preserve the column names of gt in return_matrix.
  Rcpp::StringMatrix return_matrix( gt.nrow(), gt.ncol() - 1 );
  Rcpp::List matrix_names = gt.attr("dimnames");
  Rcpp::StringVector colnames = matrix_names(1);
  colnames.erase(0);
  return_matrix.attr("dimnames") = Rcpp::List::create(matrix_names(0), colnames);
  
  for( i=0; i<gt.nrow(); i++) // Row (variant) counter,
  {
    Rcpp::checkUserInterrupt();
    // Determine the position where the query element is located.
    int position = elementNumber(gt(i,0), element);
    std::vector< std::string > allele_vector;
    
    // If we are to return alleles, 
    // we'll need to create an array for them.
    if( alleles == 1 )
    {
      allele_vector = get_allele_vector( fix(i,3), fix(i,4) );
    }
    
    // Process columns (samples).
    for( j=1; j<gt.ncol(); j++)
    {
      if( gt(i, j) == NA_STRING ){
        return_matrix(i, j-1) = NA_STRING;
      } else {
        return_matrix(i, j-1) = extractElementS( gt(i, j), position, extract );
        // Manage NAs.
        if( return_matrix(i, j-1) == "." ){ return_matrix(i, j-1) = NA_STRING; }
        // Convert to alleles
        if( alleles == 1 )
        {
          std::string gt_string = Rcpp::as< std::string >( return_matrix(i, j-1) );
        
//        Rcpp::Rcout << "gt_string: " << gt_string << "\n";
//        Rcpp::Rcout << "test ic_naC: " << gt_string == na_string << "\n";
        
          gt_string = gt2alleles( gt_string, allele_vector, allele_sep );
          return_matrix(i, j-1) = gt_string;
        }
      }
    }

  }
  
  return return_matrix;
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

  int ploidy = 1;
  int i = 0;
  int j = 0;
  int hap_col = 0;
  int hap_num = 0;

  // Determine ploidy.
  // Vcf files are typically of one ploidy.
  while( gt(i,1) == NA_STRING ){
    i++;
  }
  std::string temp = Rcpp::as< std::string >(gt(i,1));
//  Rcpp:Rcout << "GT to test ploidy: " << temp << "\n";

  // Count elements to determine ploidy.
  for(i=0; i<temp.length(); i++){
    if( temp[i] == gt_split ){ploidy++;}
  }

  if(ploidy == 1){
    // Is either haploid or is not phased.
    Rcpp::StringMatrix haps(1, 1);
    haps(0, 0) = NA_STRING;
    return haps;
  }

  // Initialize return structure
  Rcpp::StringMatrix haps(gt.nrow(), gt.ncol() * ploidy);

  Rcpp::List gt_names = gt.attr("dimnames");
  Rcpp::StringVector sample_names = gt_names(1);

  // Manage haplotype names with postfixed number.
  Rcpp::StringVector haplo_names(gt.ncol() * ploidy);    
  j = 0;
  while(j < sample_names.size()){
    hap_num = 0;
    std::string sname = Rcpp::as< std::string >(sample_names[j]);
    while(hap_num < ploidy){
      std::ostringstream stm;
      stm << hap_num ;
      haplo_names(hap_col) = sname + "_" + stm.str();
      hap_num++;
      hap_col++;
    }
    j++;
  }
  haps.attr("dimnames") = Rcpp::List::create(gt_names(0), haplo_names);


  // Iterate over variants (rows of gt)
  // Each variant has a REF and ALT, so this can't be by sample.
  // Create a vector where position zero is REF and subsequent positions are ALT.
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

    // Process the genotypes (columns) into haplotypes.
    // hap_num counts haplotypes per sample.
    // hap_col counts columns in the return matrix
    hap_col = 0;
    for(j=0; j<gt.ncol(); j++){
      Rcpp::checkUserInterrupt();
      std::vector < std::string > al_vec;
      char al_split = gt_split; // Must be single quotes!

      if( gt(i, j) == NA_STRING ){
        hap_num = 0;
        while(hap_num < ploidy){
          haps(i, hap_col) = NA_STRING;
          hap_num++;
          hap_col++;
        }
      } else {
        std::string line = Rcpp::as< std::string >(gt(i, j));
        vcfRCommon::strsplit(line, al_vec, al_split);
        hap_num = 0;
        while(hap_num < ploidy){
          // Manage missing alleles.
          if( al_vec[hap_num] == "." ){
            haps(i, hap_col) = NA_STRING;
          } else {
            int al_num = atoi(al_vec[hap_num].c_str());
            haps(i, hap_col) = alleles_vec[al_num];
          }
          hap_num++;
          hap_col++;
        }
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
