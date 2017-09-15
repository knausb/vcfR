#include <Rcpp.h>
#include <zlib.h>
#include <fstream>

// using namespace Rcpp;

/*  Write vcf body  */

//' @export
// [[Rcpp::export(name=".write_vcf_body")]]
void write_vcf_body( Rcpp::CharacterMatrix fix,
                     Rcpp::CharacterMatrix gt,
                     std::string filename="myFile.vcf.gz",
                     int mask=0 ) {
  // http://stackoverflow.com/a/5649224
  
//  
int verbose = 0;
//  int verbose = 1;
  
  if( verbose == 1 ){
    Rcpp::Rcout << "Made it into the function!\n";
  }
  
  int i = 0; // Rows
  int j = 0; // Columns
  std::string tmpstring;  // Assemble each line before writing

  // Initialize filehandle.
  gzFile fi;
  
  // Initialize file.
  // Note that gzfile does not tolerate initializing an empty file.
  // Use ofstream instead.
  if ( ! std::ifstream( filename ) ){
    if( verbose == 1 ){
      Rcpp::Rcout << "File does not exist." << std::endl;
    }
    
    std::ofstream myfile;
    myfile.open (filename, std::ios::out | std::ios::binary);
    myfile.close();
    
    // This should make valgrind hang.
    // Or not???
//    fi = gzopen( filename.c_str(), "ab" );
//    gzclose(fi);
  }

  // In order for APPEND=TRUE to work the header
  // should not be printed here.

  if( verbose == 1 ){
    Rcpp::Rcout << "Matrix fix has " << fix.nrow() << " rows (variants).\n";
  }
  
  // Manage body
  if( fix.nrow() >= 1 ){
    if( verbose == 1 ){
      Rcpp::Rcout << "Processing the body (variants).\n";
    }
    // There is at least one variant.
    fi = gzopen( filename.c_str(), "ab" );
    if (! fi) {
      Rcpp::Rcerr << "gzopen of " << filename << " failed: " << strerror (errno) << ".\n";
    }

    for(i = 0; i < fix.nrow(); i++){
      Rcpp::checkUserInterrupt();

      if(mask == 1 && fix(i,6) != "PASS" ){
        // Don't print variant.
      } else {
        // Print variant.
        j = 0;
        tmpstring = fix(i,j);
        for(j = 1; j < fix.ncol(); j++){
          if(fix(i,j) == NA_STRING){
            tmpstring = tmpstring + "\t" + ".";
          } else {
            tmpstring = tmpstring + "\t" + fix(i,j);
          }
        }

        // gt portion
        for(j = 0; j < gt.ncol(); j++){
          if(gt(i, j) == NA_STRING){
            tmpstring = tmpstring + "\t" + "./.";
          } else {
            tmpstring = tmpstring + "\t" + gt(i, j);
          }
        }

        gzwrite(fi, tmpstring.c_str(), tmpstring.size());
        gzwrite(fi,"\n",strlen("\n"));
      }
    }
    if( verbose == 1 ){
      Rcpp::Rcout << "Finished processing the body (variants).\n";
    }
    gzclose(fi);
  } else {
    if( verbose == 1 ){
      Rcpp::Rcout << "No rows (variants).\n";
    }
  }
  
//  return void;
}



