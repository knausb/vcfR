#include <Rcpp.h>
#include "vcfRCommon.h"

//using namespace Rcpp;

// Helper for std::sort.
struct greater
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

struct lesser
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a < b; }
};


// Convert vectors of strings to floats.
std::vector<float> str_vec_to_float_vec2( std::vector<std::string> str_vec ){
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
//' @param decreasing should the values be sorted decreasing (1) or increasing (0)?
//'
//' @details
//' Files containing VCF data frequently include data on allelic depth (e.g., AD).
//' This is the number of times each allele has been sequenced.
//' Our naive assumption for diploids is that these alleles should be observed at a frequency of 1 or zero for homozygous positions and near 0.5 for heterozygous positions.
//' Deviations from this expectation may indicate allelic imbalance or ploidy differences.
//' This function is intended to facilitate the exploration of allele frequencies for all positions in a sample.
//'
//' The alleles are sorted by their frequency within the function.
//' The user can then specify is the would like to calculate the frequency of the most frequent allele (allele = 1), the second most frequent allele (allele = 2) and so one.
//' If an allele is requested that does not exist it should result in NA for that position and sample.
//'
//' There are two methods to calculate a sum for the denominator of the frequency.
//' When sum_type = 0 the alleles are sorted decendingly and the first two allele counts are used for the sum.
//' This may be useful when a state of diploidy may be known to be appropriate and other alleles may be interpreted as erroneous.
//' When sum_type = 1 a sum is taken over all the observed alleles for a variant.
//'
//' @return A numeric matrix of frequencies
//' 
//' @examples
//' set.seed(999)
//' x1 <- round(rnorm(n=9, mean=10, sd=2))
//' x2 <- round(rnorm(n=9, mean=20, sd=2))
//' ad <- matrix(paste(x1, x2, sep=","), nrow=3, ncol=3)
//' colnames(ad) <- paste('Sample', 1:3, sep="_")
//' rownames(ad) <- paste('Variant', 1:3, sep="_")
//' ad[1,1] <- "9,23,12"
//' AD_frequency(ad=ad)
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix AD_frequency(Rcpp::StringMatrix ad,
                                 std::string delim = ",",
                                 int allele = 1,
                                 int sum_type = 0,
                                 int decreasing = 1
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
//        char my_split = ','; // Must be single quotes!
        char my_split = delim[0];
        std::string my_string;
        my_string = ad(i,j);
        vcfRCommon::strsplit(my_string, col_vec, my_split);

        // Recast vector of string to vector of floats.
        std::vector < float > float_vec;( col_vec.size(), 0);
        float_vec = str_vec_to_float_vec2(col_vec);

        // Sort the vector.
        if( decreasing == 1 ){
          std::sort ( float_vec.begin(), float_vec.end(), greater() );
        } else if ( decreasing == 0 ){
          std::sort ( float_vec.begin(), float_vec.end(), lesser() );
        } else {
          Rcpp::Rcerr << "Specification of 'decreasing' should be either 0 or 1.\n";
          return nam;
        }

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
          if( my_sum == 0 ){
            adf(i,j) = 0;
          } else {
            adf(i,j) = float_vec[ allele ]/my_sum;
          }
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




