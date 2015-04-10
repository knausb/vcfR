// #include <Rcpp.h>
// #include <string>

// using namespace Rcpp;


#ifndef __COMMON_H_INCLUDED__   // if x.h hasn't been included yet...
#define __COMMON_H_INCLUDED__   //   #define this so the compiler knows it has been included

class common
{
public:
//  static std::vector < std::string > strsplit();
  static void strsplit_old(std::string);
  static void strsplit(std::string&, std::vector<std::string>&, char&);
  static void foo();
  static void foo2(char[], std::vector<std::string>&, std::string&);
//  static void foo2(char);
  int bar();
};


void common::foo(){
//  Rcpp::Rcout << "In foo!\n";
}


void common::foo2(char[], std::vector<std::string>&, std::string&){
//  Rcpp::Rcout << "In foo!\n";
}


int common::bar(){
  return 1;
}

void fun(std::vector<int>& v){
//fill the vector v;
  v.push_back(1);
  v.push_back(2);
}


void common::strsplit(std::string& mystring, std::vector<std::string>& vec_o_strings, char& split){
  Rcpp::Rcout << "Inside strsplit!\n";
  vec_o_strings.push_back("Element 1");
  
}

void strsplit_old(std::string s) {
  std::string question2 = "Where do you live? ";
//  Rcpp::Rcout << "In foo!\n";
//  stringv.push_back("a");
//  stringv.push_back("b");
}

/*
void strsplit2(std::vector<std::string>& stringv) {
//std::vector < std::string > strsplit(std::string line) {
//std::vector < std::string > strsplit(std::string line, char split) {
//  std::vector < std::string > stringv;

  std::string line = "phrase with a tab\tto split on";

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
      
//  return stringv;
  
//   return 1;
}
*/



#endif 