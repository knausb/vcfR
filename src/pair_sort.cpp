#include <Rcpp.h>
using namespace Rcpp;
// Created for testing.

typedef std::pair<int,int> mypair;
bool comparator2 ( const mypair& l, const mypair& r)
   { return l.first > r.first; }

// [[Rcpp::export]]
Rcpp::DataFrame pair_sort(){
  std::vector < std::pair < int,int > > v_prime;
  v_prime.push_back(std::make_pair(5, 0));
  v_prime.push_back(std::make_pair(8, 1));
  v_prime.push_back(std::make_pair(7, 2));
  v_prime.push_back(std::make_pair(3, 3));
  
  Rcpp::NumericVector pair1;
  Rcpp::NumericVector pair2;
  
  std::sort(v_prime.begin(), v_prime.end(), comparator2);
  
  unsigned int i = 0;
  for(i=0; i<v_prime.size(); i++){
    pair1.push_back(v_prime[i].first);
    pair2.push_back(v_prime[i].second);
  }
  
  return Rcpp::DataFrame::create(_["pair1"]=pair1, _["pair2"]=pair2);
}


