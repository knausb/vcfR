#include <Rcpp.h>
#include "vcfRCommon.h"


//' @rdname is_het
//' @name is_het
//' 
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix is_het(Rcpp::StringMatrix gt,
                           Rcpp::LogicalVector na_is_false = true
){

  // NA matrix to return in case of unexpected results.
  Rcpp::LogicalMatrix nam( 1, 1 );
  nam(0,0) = NA_LOGICAL;
  
  // Create return data matrix.
  Rcpp::LogicalMatrix hets( gt.nrow(), gt.ncol() );
  hets.attr("dimnames") = gt.attr("dimnames");
  
  int i;
  int j;
  
  for( i=0; i<gt.nrow(); i++){
    for( j=0; j<gt.ncol(); j++){

//      char my_split = '|'; // Must be single quotes!
      std::string my_string;
      my_string = gt(i,j);
      std::vector < std::string > allele_vec;
//      vcfRCommon::strsplit(my_string, allele_vec, my_split);
      int unphased_as_na = 0; // 0 == FALSE
      vcfRCommon::gtsplit( my_string, allele_vec, unphased_as_na );

      // Remove missing alleles from start of vector.
      while( allele_vec[0] == "." ){
        allele_vec.erase( allele_vec.begin() );
      }
      
      // Case of all missing data.
      if( allele_vec.size() == 0 ){
        if( na_is_false ){
          hets(i,j) = false;
        } else {
          hets(i,j) = NA_LOGICAL;
        }
      }
      
      // Case of one allele.
      if( allele_vec.size() == 1 ){
        hets(i,j) = false;
      }
      
      // Case of more than one allele.
      if( allele_vec.size() > 1 ){
        std::vector < std::string > allele_vec2;
        allele_vec2.push_back( allele_vec[0] );
        allele_vec.erase( allele_vec.begin() );
        int k;
        for(k=0; k<allele_vec.size(); k++){
          if( allele_vec[k] == "." ){
          }
          if( allele_vec[k] != allele_vec2[0] ){
            allele_vec2.push_back( allele_vec[k] );
          }
        }
        if( allele_vec2.size() == 1){
          hets(i,j) = false;
        }        
        if( allele_vec2.size() > 1){
          hets(i,j) = true;
        }
      }
    }
  }
  
  return( hets );
}

