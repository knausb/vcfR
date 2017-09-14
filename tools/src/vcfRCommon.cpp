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



