#include <Rcpp.h>
#include <fstream>
#include <cerrno>
#include <zlib.h>

#include "../src/vcfRCommon.h"


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

  
  // Scroll through buffers.
  while (1) {
    Rcpp::checkUserInterrupt();
    int err;
    
    // Slurp in a buffer.
    int bytes_read;
    char buffer[LENGTH];
    bytes_read = gzread (file, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0'; // Terminate the buffer.
    
    std::string mystring(reinterpret_cast<char*>(buffer));  // Recast buffer as a string.
    mystring = lastline + mystring; // Concatenate last line to the buffer
    
    // Delimit into lines.
    std::vector < std::string > svec;  // Initialize vector of strings for parsed buffer.
    char split = '\n'; // Must be single quotes!
    vcfRCommon::strsplit(mystring, svec, split);

    
    // Scroll through lines of buffer.
    for(int i = 0; i < svec.size() - 1; i++){
      
      if(svec[i][0] == '#' && svec[i][1] == '#'){
        // Meta line, ignore.
      } else if(svec[i][0] == '#' && svec[i][1] == 'C'){
        // Process header.
        char header_split = '\t';
        vcfRCommon::strsplit(svec[i], header_vec, header_split);
        
        // Subset the header to select columns.
        std::vector<std::string> header_vec2( cols.size() );
        for(int j=0; j<cols.size(); j++){
          header_vec2[j] = header_vec[ cols[j] ];
        }
        header_vec = header_vec2;
      } else {
        // Variant line.
        
        var_num++; // Input row number.
      }
    }
    
    
    
    
        
    // Check for EOF or errors.
    if (bytes_read < LENGTH - 1) {
      if (gzeof (file)) {
        break;
      }
      else {
        const char * error_string;
        error_string = gzerror (file, & err);
        if (err) {
          Rcpp::Rcerr << "Error: " << error_string << ".\n";
//          return Rcpp::StringVector(1);
          return Rcpp::CharacterMatrix(1);
        }
      }
    }

  }
  
  
  
  // Close filehandle.
  gzclose (file);


  if(verbose == 1){
    Rcpp::Rcout << "\rProcessed variant: " << var_num;
    Rcpp::Rcout << "\nAll variants processed\n";
  }
 
//  Rcpp::CharacterMatrix gt( 2, 4 ); 
  return gt;
}