#include <Rcpp.h>
#include <fstream>
#include <zlib.h>
#include <cerrno>

// Number of records to report progress at.
const int nreport = 1000;

/* Size of the block of memory to use for reading. */
#define LENGTH 0x1000 // hexadecimel for 4096.




// [[Rcpp::export]]
Rcpp::CharacterMatrix read_body_gz2(std::string x,
                                   Rcpp::NumericVector stats,
                                   int nrows = -1,
                                   int skip = 0,
                                   Rcpp::IntegerVector cols = 0,
                                   int verbose = 1) {

  // Sort the column numbers.
  cols.sort();
  
  // Validate that column numbers have been specified.
  if(cols[0] == 0){
    Rcpp::Rcerr << "User must specify which (positive integer) columns to extract from the file.\n";
    Rcpp::StringMatrix mymatrix(1,1);
    mymatrix(0,0) = NA_STRING;
    return mymatrix;
  }
  cols = cols - 1; // R is 1-based

    
  // Initialize matrix for body data.
  if( nrows == -1 & skip == 0 ){
    nrows = stats[2];
  } else if ( nrows != -1 & skip == 0 ){
    // nrows = nrows;
  } else if ( nrows == -1 & skip > 0){
    nrows = stats[2] - skip;
  } else if ( nrows != -1 & skip > 0){
    // nrows = nrows;
  } else {
    Rcpp::Rcerr << "failed to calculate return matrix geometry.";
  }
  Rcpp::CharacterMatrix gt( nrows, cols.size() );
  

  int row_num = 0;
  
  
  // Create filehandle and open.
  gzFile file;
  file = gzopen (x.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << x << " failed: " << strerror (errno) << ".\n";
    return Rcpp::CharacterMatrix(1);
  }

  // Because the last line may be incomplete,
  // We'll typically omit it from processing and
  // concatenate it to the first line.
  // But first we'll have to initialize it.
  std::string lastline = "";
  
  // String vector to store the header (^#CHROM...).
  std::vector<std::string> header_vec;
  
  // variant counter.  
  int var_num = 0;

  
  
  
  
  
  
  // Close filehandle.
  gzclose (file);

  if( stats[1] == 0 ){
    if( verbose == 1 ){
      Rcpp::Rcout << "Warning: no header information was found! Data contains no sample names!\n";
    }
  } else {
    header_vec[0] = "CHROM";
    if( header_vec.size() == gt.ncol() ){
      gt.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(), header_vec);
    } else {
      if( verbose == 1 ){
        Rcpp::Rcout << "Warning: no header information found!\n";
      }
    }
  }

  if(verbose == 1){
    Rcpp::Rcout << "\rProcessed variant: " << var_num;
    Rcpp::Rcout << "\nAll variants processed\n";
  }
  
  return gt;
}