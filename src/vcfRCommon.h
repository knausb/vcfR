// #include <Rcpp.h>
#include <string>
#include <vector>

// using namespace Rcpp;

#ifndef __VCFRCOMMON_H_INCLUDED__   // if x.h hasn't been included yet...
#define __VCFRCOMMON_H_INCLUDED__   //   #define this so the compiler knows it has been included

class vcfRCommon
{
public:
  static void strsplit(std::string&, std::vector<std::string>&, char&);
//  static void strsplit(unsigned char[]&, std::vector<std::string>&, char&);
};



#endif 