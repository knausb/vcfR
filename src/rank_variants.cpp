#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::DataFrame rank_variants(Rcpp::DataFrame variants,
                              Rcpp::NumericVector ends,
                              Rcpp::NumericVector score){
  // Rank variants by window and by a vector of scores.
  // Input is a DataFrame of variant information which includes
  // the chromosomal position of each variant.
  // This may be originally derived from a vcf file.  
  // It should also contain information on whether each variant has been masked.
  // Also input is a vector of window ending positions and a vector of scores for each
  // variant.  The scores may be extracted from a vcf file, such
  // as the genotype quality filed (GQ) or any sort of composite.
  // The variants with the greatest score are given the lowest
  // rank order.
  
  Rcpp::NumericVector pos = variants["POS"];
  Rcpp::NumericVector mask = variants["mask"];
  
  Rcpp::NumericVector win_num(score.size());  // Vector of window numbers.
  Rcpp::NumericVector win_rank(score.size()); // Vector to be returned as per window rank.
  Rcpp::NumericVector window;   // Vector to hold information for a single 
                                // window which will be sorted.



  return DataFrame::create(variants, _["window_number"]=win_num, _["window_rank"]=win_rank);
//  return variants;
}
