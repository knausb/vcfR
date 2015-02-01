#include <Rcpp.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar


//' Extract numeric data from genotype field of VCF
//' 
//' @param x A dataframe.
//' @export
// [[Rcpp::export]]
NumericMatrix extractGT2NM(DataFrame x) {
  Rcout << "Inside of extractGT2NV";
  Rcout << "\n";
  int dfsize = x.size();
  Rcout << "DF size: ";
  Rcout << dfsize;
  Rcout << "\n";
  CharacterVector format = x[0];
//  StringVector format2 = x[0];
  Rcout << "format 0: ";
  Rcout << format[0];
  Rcout << "\n";
  int vsize = format.size();
  Rcout << "V size: ";
  Rcout << vsize;
  Rcout << "\n";
//  NumericMatrix outM(vsize, dfsize-1);
  // Allocate return matrix
//  NumericMatrix outM(x.nrows(), x.ncol()-1);
  NumericMatrix outM(x.nrows(), x.size()-1);

  /*  In the genotype portion of a vcf file
      we have samples in columns and variants
      (i.e., loci) in rows.  The first column
      is a format specifier as opposed to a
      sample.  Each variant may have a
      different format, so we want to process
      by row.  Our input DataFrame is actually
      a vector of columnns (samples).  If we
      transform our DataFrame we should be able
      to process it by variant instead.
  */
//  x = std::transform(x);

//  return 1;
  return outM;
}
