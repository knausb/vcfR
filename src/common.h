// #include <Rcpp.h>
// #include <string>

// using namespace Rcpp;

#ifndef __COMMON_H_INCLUDED__   // if x.h hasn't been included yet...
#define __COMMON_H_INCLUDED__   //   #define this so the compiler knows it has been included

class common
{
public:
  static void strsplit(std::string&, std::vector<std::string>&, char&);
//  static void strsplit(unsigned char[]&, std::vector<std::string>&, char&);
};


void common::strsplit(std::string& mystring, std::vector<std::string>& vec_o_strings, char& split){
  // mystring is a string to be split on the character 'split'.
  // vec_o_strings is empty and will be pushed on to.

  int start = 0;
  int i=0;

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


#endif 