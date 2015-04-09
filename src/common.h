#include <Rcpp.h>
using namespace Rcpp;


#ifndef __COMMON_H_INCLUDED__   // if x.h hasn't been included yet...
#define __COMMON_H_INCLUDED__   //   #define this so the compiler knows it has been included

class common
{
public:
  std::vector < std::string > strsplit();
//  void foo();
//  int bar;
};



std::vector < std::string > strsplit(std::string line, char split) {
  std::vector < std::string > stringv;

  int start=0;
  int j = 0;

  for(int i=1; i<line.size(); i++){
    if( line[i] == '\t'){
      std::string temp = line.substr(start, i - start);
//      Rcout << "  i: " << i << ", temp: " << temp << "\n";
//      stringv[j] = temp;
      stringv.push_back(temp);
      j++;
      start = i+1;
      i = i+1;
    }
  }

  // Handle last element.
  std::string temp = line.substr(start, line.size());
//  stringv[j] = temp;
  stringv.push_back(temp);
      
  return stringv;
  
//   return 1;
}




#endif 