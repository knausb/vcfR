#include <Rcpp.h>
using namespace Rcpp;

//' Multiply a number by two
//'
//' @param x A single integer.
//' @export
// [[Rcpp::export(name=".timesTwo")]]
int timesTwo(int x) {
   return x * 2;
}
