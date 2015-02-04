#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar



//' Extract numeric data from genotype field of VCF
//' 
//' @param x A NumericMatrix
//' @param maxbp max
//' @param winsize
//' @export
// [[Rcpp::export]]
NumericMatrix NM2winNM(NumericMatrix x, std::vector<int> pos, int maxbp, int winsize=100) {
//  int nwins = 100;
//  int nwins = x.nrow();
  int nwins = ceil(maxbp/winsize);
//  std::vector<std::string> nm = x.attributeNames();
//  std::vector<std::string> nm = x.attr("dimnames");
//  for(int i = 0; i < nm.size(); i++){
//    Rcout << nm[i] << "\n";
//  }
//  Rcout << x.dimnames[0] << "\n";
  
//  std::vector<std::string> nm = x.rowNames();
///winsize;
//  NumericMatrix outM(nwins, x.size());
//  NumericMatrix outM(x.nrows(), x.size()-1);
  NumericMatrix outM(nwins, x.ncol());
  
  return outM;
}


