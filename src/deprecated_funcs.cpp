#include <Rcpp.h>
#include <fstream>
#include <zlib.h>
#include "vcfRCommon.h"

// These are deprecated functions I hope to remove.


Rcpp::StringMatrix DataFrame_to_StringMatrix( Rcpp::DataFrame df, int depr = 1 ){
  
  Rcpp::StringVector sv = df(0);
  Rcpp::StringMatrix sm(sv.size(), df.size());

  sm.attr("col.names") = df.attr("col.names");
  sm.attr("row.names") = df.attr("row.names");

  if( depr == 1 ){
    Rcpp::Rcerr << "The function DataFrame_to_StringMatrix was deprecated in vcfR 1.6.0" << std::endl;
    Rcpp::Rcerr << "If you use this function and you would like to advocate its persistence, please contact the maintainer." << std::endl;
    Rcpp::Rcerr << "The maintainer of this package can be found with" << std::endl;
    Rcpp::Rcerr << "maintainer('vcfR')" << std::endl;
    return sm;
  }
  
  for(int i=0; i < df.size(); i++){
    sv = df(i);
    for(int j=0; j < sv.size(); j++){
      sm(j, i) = sv(j);
    }
  }

  return sm;
}


/* Write data to fasta file */

//' @export
// [[Rcpp::export(name=".write_fasta")]]
void write_fasta( Rcpp::CharacterVector seq,
                  std::string seqname, 
                  std::string filename, 
                  int rowlength=80,
                  int verbose=1, int depr = 1) {
//  rowlength=rowlength-1;
  FILE * pFile;
//  pFile=fopen(filename.c_str(),"wt");
  pFile=fopen(filename.c_str(),"at");
  int i = 0;
//  unsigned int i = 0;
  
  if( depr == 1 ){
    Rcpp::Rcerr << "The function write_fasta was deprecated in vcfR 1.6.0" << std::endl;
    Rcpp::Rcerr << "If you use this function and you would like to advocate its persistence, please contact the maintainer." << std::endl;
    Rcpp::Rcerr << "The maintainer of this package can be found with" << std::endl;
    Rcpp::Rcerr << "maintainer('vcfR')" << std::endl;
    Rcpp::stop("");
  }
  
  if(verbose == 1){
    Rcpp::Rcout << "Processing sample: " << seqname << "\n";
  }

  putc ('>' , pFile);
  for(i=0; (unsigned)i<seqname.size(); i++){
    putc (seqname[i] , pFile);
  }
  putc ('\n' , pFile);

  // Number of records to report progress at.
  const int nreport = 1000;

  putc (Rcpp::as< char >(seq[0]) , pFile);
  for(i=1; i<seq.size(); i++){
    Rcpp::checkUserInterrupt();
//    putc (seq[i][0] , pFile);
    if( i % rowlength == 0){
      putc('\n', pFile);
    }
    putc (Rcpp::as< char >(seq[i]) , pFile);
    if(i % nreport == 0 && verbose == 1){
      Rcpp::Rcout << "\rNucleotide " << i << " processed";
    }
  }
  putc('\n', pFile);
  fclose (pFile);
  if(verbose == 1){
    Rcpp::Rcout << "\rNucleotide " << i << " processed\n";
  }
//  return 0;
}


double extractElementD(Rcpp::String x, int number=1, int depr = 1){
  //
  // x is a string similar to:
  // GT:GQ:DP:RO:QR:AO:QA:GL
  //
  // number is the position in the colon delimited 
  // string which needs to be extracted.
  //
//  int count = 0;
  int start = 0;
  int pos = 1;
  std::string istring = x;
  std::string teststring;
  unsigned int i = 0;

  if( depr == 1 ){
    Rcpp::Rcerr << "The function extractElementD was deprecated in vcfR 1.6.0" << std::endl;
    Rcpp::Rcerr << "If you use this function and you would like to advocate its persistence, please contact the maintainer." << std::endl;
    Rcpp::Rcerr << "The maintainer of this package can be found with" << std::endl;
    Rcpp::Rcerr << "maintainer('vcfR')" << std::endl;
    Rcpp::stop("");
  }
    
  for(i=1; i <= istring.size(); i++){
    if(istring[i] == ':'){
      if(pos == number){
        teststring = istring.substr(start, i-start);
        double teststring2 = atof(teststring.c_str());
        return teststring2;
//        return std::stod(teststring);
      } else {
        start = i+1;
        pos++;
        i++;
      }
    }
  }
  // If we get here we did not find the element.
  return(0);
}



