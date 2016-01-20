#include <Rcpp.h>
#include <fstream>
#include <cerrno>
#include <zlib.h>

#include "../src/vcfRCommon.h"


// Number of records to report progress at.
const int nreport = 1000;

/* Size of the block of memory to use for reading. */
#define LENGTH 0x1000 // hexadecimel for 4096.


void proc_body_line(Rcpp::CharacterMatrix gt,
                    int var_num,
                    std::string myline,
                    Rcpp::IntegerVector cols){
  
  char split = '\t'; // Must be single quotes!
  std::vector < std::string > data_vec;
  
  vcfRCommon::strsplit(myline, data_vec, split);
  
  Rcpp::Rcout << "\n";
  Rcpp::Rcout << "    Inside proc_body_line" << "\n";
  Rcpp::Rcout << "    proc_body_line::var_num: " << var_num << "\n";
  Rcpp::Rcout << "    cols.size(): " << cols.size() << "\n";
  Rcpp::Rcout << "\n";
  
  for(int i = 0; i < cols.size(); i++){
//    if( data_vec[ cols[i] ] == "." ){
//      gt(var_num, i) = NA_STRING;
//    }
//    } else if( data_vec[ cols[i] ][0] == '.' & data_vec[ cols[i] ][2] == '.' & data_vec[ cols[i] ].size() == 3 ){
//      gt(var_num, i) = NA_STRING;
//    } else {
//      gt(var_num, i) = data_vec[ cols[i] ];
//    }
  }
  
}



// [[Rcpp::export]]
Rcpp::CharacterMatrix read_body_gz2(std::string x,
                                   Rcpp::NumericVector stats,
                                   int nrows = -1,
                                   int skip = 0,
                                   Rcpp::IntegerVector cols = 0,
                                   int verbose = 1) {

  // NA matrix for unexpected results.
  Rcpp::StringMatrix na_matrix(1,1);
  na_matrix(0,0) = NA_STRING;
  
  // Sort the column numbers.
  cols.sort();

  // Validate that column numbers have been specified.
  if(cols[0] == 0){
    Rcpp::Rcerr << "User must specify which (positive integer) columns to extract from the file.\n";
    return na_matrix;
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

  // Create filehandle and open.
  gzFile file;
  file = gzopen (x.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << x << " failed: " << strerror (errno) << ".\n";
    return na_matrix;
  }

  // Because the last line may be incomplete,
  // We'll typically omit it from processing and
  // concatenate it to the first line.
  // But first we'll have to initialize it.
  std::string lastline = "";
  
  // String vector to store the header (^#CHROM...).
  std::vector<std::string> header_vec;

  int var_num = 0; // input variant counter.
  int row_num = 0; // Retained variant counter.
  int err;
  
  
  if( verbose == 1 ){
    Rcpp::Rcout << "Character matrix gt created.\n";
    Rcpp::Rcout << "Character matrix gt rows: ";  Rcpp::Rcout << gt.rows();
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "Character matrix gt cols: ";  Rcpp::Rcout << gt.cols();
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "skip: ";  Rcpp::Rcout << skip;
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "nrows: ";  Rcpp::Rcout << nrows;
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "row_num: ";  Rcpp::Rcout << row_num;
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "\n";
  }
  
  
  // Scroll through buffers.
  while (1) {
    Rcpp::checkUserInterrupt();

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
        Rcpp::Rcout << row_num << "|" << var_num << " ";
        
        if ( var_num >= skip & row_num < nrows ){
          proc_body_line(gt, row_num, svec[i], cols);
          row_num++; // Return matrix row number.
        }
        var_num++; // Input row number.
      }
    }
    
    if( nrows != -1 & row_num >= nrows ){
      break;
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
//          return Rcpp::CharacterMatrix(1);
          return na_matrix;
        }
      }
    }

  }
  
  
  
  // Close filehandle.
  gzclose (file);


  if(verbose == 1){
//    Rcpp::Rcout << "\rProcessed variant: " << var_num;
//    Rcpp::Rcout << "\nAll variants processed\n";
  }
 
//  Rcpp::CharacterMatrix gt( 2, 4 ); 
  return gt;
}