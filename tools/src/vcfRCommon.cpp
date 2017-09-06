#include <Rcpp.h>
#include "vcfRCommon.h"
#include <string>


// From:
// https://github.com/RcppCore/Rcpp/issues/636#issuecomment-280985661
/*
void R_init_vcfR(DllInfo* info) {
	R_registerRoutines(info, NULL, NULL, NULL, NULL);
	R_useDynamicSymbols(info, TRUE);
}
*/


//static 
//static 
void vcfRCommon::strsplit(std::string& mystring, std::vector<std::string>& vec_o_strings, char& split){
//void strsplit(std::string& mystring, std::vector<std::string>& vec_o_strings, char& split){
  // mystring is a string to be split on the character 'split'.
  // vec_o_strings is empty and will be pushed on to.

  int start = 0;
  unsigned int i=0;

  for(i = 1; i < mystring.size(); i++){
    if( mystring[i] == split){
      std::string temp = mystring.substr(start, i - start);
      vec_o_strings.push_back(temp);
      start = i+1;
      i = i+1;
    }
  }

  // Handle the last element.
  std::string temp = mystring.substr(start, i - start);
  vec_o_strings.push_back(temp);
  
}


void vcfRCommon::gtsplit(std::string& mystring,
                         std::vector<std::string>& vec_o_strings,
                         int& unphased_as_na){

  // mystring is a string of genotypes to be split to a character.
  // Genotypes may be delimited as | or /.
  // vec_o_strings is empty and will be pushed on to.
  // Sometimes, genotypes delimited with / may be undesireable.
  // In this case missing data should be returned.
  
//  Rcpp::Rcout << "In gtsplit\n";
  
  int start = 0;
  unsigned int i = 0;
  int is_phased = 0;
  
  char split1 = '|'; // Must be single quotes!
  char split2 = '/'; // Must be single quotes!

  // Iterate through genotype string looking for delimiters.  
  for(i = 0; i < mystring.size(); i++){
    if( mystring[i] == split1 ){
      // Found a delimiter.
      is_phased = 1;
      std::string temp = mystring.substr(start, i - start);
      vec_o_strings.push_back(temp);
      start = i+1;
      i = i+1;
    } else if ( mystring[i] == split2 ){
      // Found a delimiter.
      is_phased = 0;
      if( unphased_as_na == 1 ){
        vec_o_strings.push_back( "." );
      } else {
        std::string temp = mystring.substr(start, i - start);
        vec_o_strings.push_back(temp);
        start = i+1;
        i = i+1;
      }
    }
  }
  
  // Handle the last element.
  if( is_phased == 0 && unphased_as_na == 1 ){
    vec_o_strings.push_back( "." );
  } else {
    std::string temp = mystring.substr(start, i - start);
    vec_o_strings.push_back(temp);
  }

  // Debug
  /*
  Rcpp::Rcout << "mystring: " << mystring << " ";
  Rcpp::Rcout << "vec_o_strings: ";
  for(i=0; i<vec_o_strings.size(); i++){
    Rcpp::Rcout << " " << vec_o_strings[i];
  }
  Rcpp::Rcout << "\n";
  */
}




void vcfRCommon::gtdelim(std::string& mystring,
                         std::vector<std::string>& vec_o_strings){
  // Collect a genotype's allelic delimitors.
  
  char split1 = '|'; // Must be single quotes!
  char split2 = '/'; // Must be single quotes!
  
  unsigned int i = 0;
  
  for(i = 0; i < mystring.size(); i++){
    if( ( mystring[i] == split1 ) | ( mystring[i] == split2 ) ){
      std::stringstream ss;
      ss << mystring[i];
//      vec_o_strings.push_back( ss );
      std::string myDelim;
      ss >> myDelim;
      vec_o_strings.push_back( myDelim );
//      vec_o_strings.push_back( split1 );
    }
  }
}



