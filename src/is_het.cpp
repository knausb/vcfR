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
//  Rcpp::LogicalMatrix nam( 1, 1 );
//  nam(0,0) = NA_LOGICAL;
  
  // Initialize return data matrix.
  Rcpp::LogicalMatrix hets( gt.nrow(), gt.ncol() );
  hets.attr("dimnames") = gt.attr("dimnames");
  
  int i;
  int j;
  
  for( i=0; i<gt.nrow(); i++){
    for( j=0; j<gt.ncol(); j++){

      // Parse genotype string into alleles.
      std::string my_string;
      my_string = gt(i,j);
      Rcpp::Rcout << my_string << "\n";
      
      std::vector < std::string > allele_vec;
//      vcfRCommon::strsplit(my_string, allele_vec, my_split);
      int unphased_as_na = 0; // 0 == FALSE
      vcfRCommon::gtsplit( my_string, allele_vec, unphased_as_na );

      Rcpp::Rcout << "allele_vec.size(): " << allele_vec.size() << "\n";
      Rcpp::Rcout << "allele_vec[0]: " << allele_vec[0] << "\n";
      
      // Remove missing alleles from start of vector.
      while( allele_vec[0] == "." ){
        allele_vec.erase( allele_vec.begin() );
      }
      Rcpp::Rcout << "allele_vec.size(): " << allele_vec.size() << "\n";
      
      // Case of all missing alleles.
      if( allele_vec.size() == 0 ){
        if( na_is_false ){
          hets(i,j) = false;
        } else {
          hets(i,j) = NA_LOGICAL;
        }
      }
      
      // Case of one allele.
//      if( allele_vec.size() == 1 ){
//        hets(i,j) = false;
//      }
      
      // Case of alleles present.
      // Initialize new vector of alleles with first element of allele_vec.
      std::vector < std::string > allele_vec2;
      allele_vec2.push_back( allele_vec[0] );
//      if( allele_vec.size() > 0 ){
//        allele_vec.erase( allele_vec.begin() );
        
      // Scroll through vector looking for a second allele.
      int k;
      for(k=1; k<allele_vec.size(); k++){
        if( allele_vec[k] == "." ){
          // Ignore missing alleles.
        }
        if( allele_vec2[0] != allele_vec[k] ){
          Rcpp::Rcout << "New allele: " << allele_vec[k] << "\n";
          allele_vec2.push_back( allele_vec[k] );
        }
//        }
      }
      
      Rcpp::Rcout << "allele_vec2 size: " << allele_vec2.size() << "\n";
      Rcpp::Rcout << allele_vec2[0];
//      int k;
      for(k=1; k<allele_vec2.size(); k++){
        Rcpp::Rcout << "," << allele_vec2[k];
      }
      Rcpp::Rcout << "\n";
      Rcpp::Rcout << "\n";
        
      if( allele_vec2.size() == 1){
        hets(i,j) = false;
      }        
      if( allele_vec2.size() > 1){
        hets(i,j) = true;
      }
        
//      }
    }
  }
  
  return( hets );
}

