#include <Rcpp.h>
using namespace Rcpp;
//#include <algorithm>    // std::sort

const int na_int = -9999;


// Create type and function to help sort scores.
// This funciton is used by rank_variants.
// http://stackoverflow.com/a/1577627
//
typedef std::pair<int,int> mypair;
bool comparator ( const mypair& l, const mypair& r)
   { return l.first > r.first; }
   
//typedef std::pair<int,int> mypair;
bool minimize ( const mypair& l, const mypair& r)
   { return l.first < r.first; }


//' @export
// [[Rcpp::export(name=".rank_variants")]]
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
  // as the genotype quality filed (GQ), or any sort of composite score.
  // The variants with the greatest score are given the lowest
  // rank order.
  
  // Sorting on pairs from:
  // http://stackoverflow.com/a/1577627
  
  Rcpp::NumericVector pos = variants["POS"];
  Rcpp::LogicalVector mask = variants["mask"];
  
  Rcpp::NumericVector win_num(score.size());  // Vector of window numbers.
  Rcpp::NumericVector win_rank(score.size()); // Vector to be returned as per window rank.
//  Rcpp::NumericVector rank;   // Vector to hold information for a single 
                              // window which will be sorted.

  // Create a vector of vectors of pairs.
  // There will be a vector of windows,
  // where each window will be a vector of
  // a pair containing a variant's score and position in the window.
//  std::vector< std::vector< std::pair< int, int > > > vec_vec_pair;
  // or
  std::vector< std::pair< float, int > > vec_pair;
  std::vector< std::pair< int, int > > rank_pair;
//  std::vector < int > ranks;
  Rcpp::IntegerVector ranks;
  
  int win = 0;
  int i = 0;  // Variant counter.
  int j = -1;  // Variant within window counter.

  // Iterate to first window.
  while(pos(0) > ends(win)){
    win++;
    }

  for(i=0; i<score.size(); i++){
//    Rcout << "Counter: " << i << " Position: " << pos(i) << "\n";
    Rcpp::checkUserInterrupt();
    if( pos(i) < ends(win) ){
      // Position is within current window.
      win_num(i) = win;
      j++;
      // Handle missing data.
      if(score(i) == NA_REAL || mask(i) == FALSE){
        vec_pair.push_back( std::make_pair( na_int, j) );
      } else {
        vec_pair.push_back( std::make_pair( score(i), j) );
      }

    } else {
      // Begin a new window.

      // Sort vector of pairs for present window.
      // vec_pair contains score and position.
      std::sort(vec_pair.begin(), vec_pair.end(), comparator);
      
      // Rank within the window.
      unsigned int k = 0;
      for(k=0; k<vec_pair.size(); k++){
//        Rcout << "Window number: " << win << ", pair1: " << vec_pair[k].first << ", pair2: " << vec_pair[k].second  << "\n";
          rank_pair.push_back( std::make_pair( vec_pair[k].second, k + 1 ) );
//        Rcout << "Window number: " << win << ", Rpair1: " << rank_pair[k].first << ", Rpair2: " << rank_pair[k].second  << "\n";
      }
      
      // Restore original order.
      // rank_pair contains position and rank.
      std::sort(rank_pair.begin(), rank_pair.end(), minimize);
      
      for(k=0; k<vec_pair.size(); k++){
//        Rcout << "Window number: " << win << ", Rpair1: " << rank_pair[k].first << ", Rpair2: " << rank_pair[k].second  << "\n";
//        Rcout << "Window number: " << win << " Pair number: " << k;
//        Rcout << " Rpair1: " << rank_pair[k].first << " Rpair2: " << rank_pair[k].second << "\n";
//        Rcout << "Pair: " << rank_pair[k].second << "\n";
        if( rank_pair[k].second == na_int ){
          ranks.push_back( NA_INTEGER );
        } else {
          ranks.push_back( rank_pair[k].second );
        }
      }
//      Rcout << "\n";
      // Increment to next window.
      while(pos(i) > ends(win)){win++;}
      win_num(i) = win;
      rank_pair.erase(rank_pair.begin(), rank_pair.end());
      vec_pair.erase(vec_pair.begin(), vec_pair.end());
      j = 0;
      vec_pair.push_back( std::make_pair( score(i), j) );
    }
//    Rcout << "Counter: " << i << " Window: " << win_num(i) << " Position:" << pos(i) << " Score: " << score(i);
//    Rcout << " Rank: " << ranks[i];
//    Rcout << "\n";
  }
  
  // Handle the last window.
  std::sort(vec_pair.begin(), vec_pair.end(), comparator);
  unsigned int k = 0;
  for(k=0; k<vec_pair.size(); k++){
    rank_pair.push_back( std::make_pair( vec_pair[k].second, k + 1) );
  }
  std::sort(rank_pair.begin(), rank_pair.end(), minimize);
  
  for(k=0; k<vec_pair.size(); k++){
    if( rank_pair[k].second == na_int){
      ranks.push_back( NA_INTEGER );
    } else {
      ranks.push_back( rank_pair[k].second );
    }
  }
  
//  for(i=0; i<ranks.size(); i++){
//    Rcout << i << "\n";
//    Rcout << "Ranks: " << ranks[i] << "\n";
//  }


//  Rcpp::IntegerVector rank = wrap(ranks);
  return Rcpp::DataFrame::create(variants, _["window_number"]=win_num, _["rank"]=ranks);
//  return Rcpp::DataFrame::create(variants, _["window_number"]=win_num);
}




