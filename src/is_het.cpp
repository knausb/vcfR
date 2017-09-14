#include <Rcpp.h>
#include "vcfRCommon.h"



//' @rdname is_het
//' @name is_het
//' 
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::LogicalMatrix is_het(Rcpp::StringMatrix x,
                           Rcpp::LogicalVector na_is_false = true
){

  // NA matrix to return in case of unexpected results.
//  Rcpp::LogicalMatrix nam( 1, 1 );
//  nam(0,0) = NA_LOGICAL;
  
  // Initialize return data matrix.
  Rcpp::LogicalMatrix hets( x.nrow(), x.ncol() );
  hets.attr("dimnames") = x.attr("dimnames");
  
  int i;
  int j;
  unsigned int k;
  for( i=0; i<x.nrow(); i++ ){
    for( j=0; j<x.ncol(); j++){

      // Parse genotype string into alleles.
      std::string my_string;
      if( x(i,j) == NA_STRING ){
        my_string = ".";
      } else {
        my_string = x(i,j);
      }


      std::vector < std::string > allele_vec;
//      vcfRCommon::strsplit(my_string, allele_vec, my_split);
      int unphased_as_na = 0; // 0 == FALSE
      vcfRCommon::gtsplit( my_string, allele_vec, unphased_as_na );

//      Rcpp::Rcout << "gtsplit returned: " << allele_vec[0];
//      for( k=1; k<allele_vec.size(); k++){
//        Rcpp::Rcout << "," << allele_vec[k];
//      }
//      Rcpp::Rcout << "\n";


      // Initialize new vector of alleles with first element of allele_vec.
      std::vector < std::string > allele_vec2;

      // Scroll through vector looking for alleles.
      for(k=0; k<allele_vec.size(); k++){
        if( allele_vec[k] == "." ){
          // Found missing value.
          // Delete and bail out.
          while( allele_vec2.size() > 0 ){
            allele_vec2.erase( allele_vec2.begin() );
          }
          k = allele_vec.size();
        } else if( allele_vec2.size() == 0 ){
          // Initialize.
          allele_vec2.push_back( allele_vec[k] );
        } else if( allele_vec2[0] != allele_vec[k] ){
          allele_vec2.push_back( allele_vec[k] );
        }
      }
      
//      Rcpp::Rcout << "allele_vec2.size(): " << allele_vec2.size();
//      Rcpp::Rcout << "\n";
//      Rcpp::Rcout << "\n";
      
      // Score return value.
      if( allele_vec2.size() == 0){
        if( na_is_false[0] == true ){
          hets(i,j) = false;
        } else if( na_is_false[0] == false ){
          hets(i,j) = NA_LOGICAL;
        }
      } else if( allele_vec2.size() == 1){
        hets(i,j) = false;
      } else if( allele_vec2.size() > 1){
        hets(i,j) = true;
      }

    }
  }
  
  return( hets );
}

