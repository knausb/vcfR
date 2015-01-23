#include <Rcpp.h>
#include <iostream>

using namespace Rcpp;


// To compile the C++ code, use sourceCpp("path/to/file.cpp")


// [[Rcpp::export]]
int one() {
//  std::cout << "In a function.\n"; 
  return 1;
}


