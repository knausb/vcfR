#include <Rcpp.h>
#include <fstream>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::NumericVector vcf_stats(std::string x) {
  Rcpp::NumericVector stats(4);
  stats.names() = Rcpp::StringVector::create("meta", "header", "variants", "columns");

  std::string line;  // String for reading file into
  long int i;

  std::ifstream myfile;
  myfile.open (x.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }
  
  // Loop over the file.
  while ( getline (myfile,line) ){
    if(line[0] == '#' && line[1] == '#'){
      stats[0]++;
    } else if (line[0] == '#'){
      stats[1] = stats[0] + 1;
    } else {
      stats[2]++;
    }
    i++;
  }
  myfile.close();

  // Reopen the file to count columns for first variant.
  myfile.open (x.c_str(), std::ios::in);
  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }

  i = 0;
  while(i <= stats[1]){
    getline (myfile,line);
    i++;
  }

  myfile.close();

//  Rcout << "Line: " << line << "\n";
  for(int j = 0; j < line.size(); j++ ){
//    Rcout << line[j] << "\n";
    if(line[j] == '\t'){
      stats[3]++;
    }
  }
  stats[3]++;
  
  return stats;
}



// [[Rcpp::export]]
Rcpp::List vcf_meta(std::string x, Rcpp::NumericVector stats) {
  // stats consists of elements ("meta", "header", "variants", "columns");
  Rcpp::List meta(stats[0]);
  std::string line;  // String for reading file into
  
  std::ifstream myfile;
  myfile.open (x.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }
  
  // Loop over the file.
  int i = 0;
  while ( i < stats[0] ){
//    Rcout << i << "\n";
    getline (myfile,line);
    meta(i) = line;
    i++;
  }
  myfile.close();

  return meta;
}