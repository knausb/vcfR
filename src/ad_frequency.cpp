#include <Rcpp.h>
#include "vcfRCommon.h"

// Helper for std::sort.
struct greater
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};


// Convert vectors of strings to floats.
std::vector<float> str_vec_to_float_vec( std::vector<std::string> str_vec ){
  // Initialize return vector.
  std::vector<float> float_vec( str_vec.size(), 0 );
  
  int i;
  for( i=0 ; i < str_vec.size() ; i++ ){
  
//    Rcpp::Rcout << "  " << str_vec[i] << "\n";
    std::istringstream ss0(str_vec[i]);
    if ( !( ss0 >> float_vec[i] ) ){
      // error: didn't convert to a float
      Rcpp::Rcerr << "Failed to convert to a float.\n";
    }
  }
  return float_vec;
}


//' @title AD_frequency
//' @name AD_frequency
//' @rdname AD_frequency
//'
//' @description
//' Create allele frequencies from matrices of allelic depths (AD)
//'
//' @param ad a matrix of allele depths (e.g., "7,2")
//' @param allele which (1-based) allele to report frequency for
//' @param sum_type type of sum to calculate, see details
//' @param delim character that delimits values
//'
//' @details
//' To do.
//'
//' There are two methods to calculate a sum for the denominator of the frequency.
//' When sum_type = 0 the alleles are sorted decendingly and the first two allele counts are used for the sum.
//' This may be useful when a state of diploidy may be known to be appropriate and other alleles may be interpreted as erroneous.
//' When sum_type = 1 a sum is taken over all the observed alleles for a variant.
//'
//' @return a numeric matrix of frequencies
//' [[Rcpp::export]]
Rcpp::NumericMatrix AD_frequency(Rcpp::StringMatrix ad,
                                 int allele = 1,
                                 int sum_type = 0,
                                 char delim = ','
                                 ) {

  // Initialize return data structure.
  Rcpp::NumericMatrix adf( ad.nrow(), ad.ncol() );
  adf.attr("dimnames") = ad.attr("dimnames");
  
  // NA matrix to return in case of unexpected results.
  Rcpp::NumericMatrix nam( 1, 1 );
  nam(0,0) = NA_REAL;
    
  if( allele < 1){
    Rcpp::Rcerr << "Specified allele number is less than one.\n";
    return nam;
  }
  
  // R is one based, C++ zero based.
  allele = allele - 1;

  int i;
  int j;

  for(i=0; i<adf.nrow(); i++){   // Count rows (variants).
    for(j=0; j<adf.ncol(); j++){ // Count columns (samples).

      if( ad(i,j) != NA_STRING ){
        std::vector < std::string > col_vec;
        char my_split = ','; // Must be single quotes!
        std::string my_string;
        my_string = ad(i,j);
        vcfRCommon::strsplit(my_string, col_vec, my_split);

        // Recast vector of string to vector of floats.
        std::vector < float > float_vec;( col_vec.size(), 0);
        float_vec = str_vec_to_float_vec(col_vec);

        // Sort the vector.
        std::sort ( float_vec.begin(), float_vec.end(), greater() );
        
        // Calculate a sum.
        float my_sum = 0;
        if( sum_type == 0 ){
          my_sum = float_vec[0] + float_vec[1];
        } else if ( sum_type == 1 ){
          for(std::vector<float>::iterator it = float_vec.begin(); it != float_vec.end(); ++it)
            my_sum += *it;
        } else {
          Rcpp::Rcerr << "Undefined sum type.\n";
          return nam;
        }
        
        // Calculate the frequency.
        if( allele < float_vec.size() ){
          adf(i,j) = float_vec[ allele ]/my_sum;
        } else {
          adf(i,j) = NA_REAL;
        }
        
      } else {
        adf(i,j) = NA_REAL;
      }
    }
  }

  return adf;
}

