#include <Rcpp.h>
#include <fstream>
#include <string>

using namespace Rcpp;


// Read vcf header
// 
// [[Rcpp::export]]
std::vector<std::string> readVcfHeader(String x) {
  std::vector<std::string> header;
  std::string x2 = x;
  std::string line;
  std::string temp;


  std::ifstream myfile;
  myfile.open (x2.c_str(), std::ios::in);

  if (myfile.is_open()){
    getline (myfile, line);
    temp = line.substr(0,2);
    if(temp == "##"){
      header.push_back(line);
      while(temp == "##"){
        getline (myfile, line);
        temp = line.substr(0,2);
        if(temp == "##"){
          header.push_back(line);
        }
      }
    }
    myfile.close();
  }

  else Rcout << "Unable to open file";

  return header;
}


std::vector<std::string> splitTab(std::string line){
  // Based on:
  // http://stackoverflow.com/a/14266139
  std::string delimiter = "\t";  // Delimiting string
  std::string token;
  std::vector<std::string> tempv;
  
  size_t pos = 0;
  while ((pos = line.find(delimiter)) != std::string::npos) {
    token = line.substr(0, pos);
    tempv.push_back(token);
    line.erase(0, pos + delimiter.length());
  }
  tempv.push_back(line);
  
  return(tempv);
}


Rcpp::CharacterVector strsplit(std::string line, std::string delimiter = "\t"){
  // Based on:
  // http://stackoverflow.com/a/14266139
  std::string token;
  std::vector<std::string> tempv;
  
  size_t pos = 0;
  while ((pos = line.find(delimiter)) != std::string::npos) {
    token = line.substr(0, pos);
    tempv.push_back(token);
    line.erase(0, pos + delimiter.length());
  }
  tempv.push_back(line);
  
  Rcpp::CharacterVector charvec(tempv.size());
  for(int i=0; i<tempv.size(); i++){charvec[i] = tempv[i];}
  return charvec;
}



CharacterMatrix addRow(CharacterMatrix x){
  CharacterMatrix cm2(x.nrow()+1, x.ncol());
  for(int i=0; i<x.nrow(); i++){
    cm2(i, _) = x(i, _);
  }
  return(cm2);  
}



// Read vcf body
// 
// [[Rcpp::export]]
CharacterMatrix readVcfBody(String x) {
//  CharacterMatrix body;
  std::string x2 = x;
  std::string line;
  std::string temp;
  std::vector<std::string> tempv;
  std::vector<std::string> header;
//  Rcpp::CharacterVector header;
  std::vector<std::string> linesv;
  int i;
  int j;


  std::ifstream myfile;
  myfile.open (x2.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }

  // Iterate past meta lines.
  getline (myfile, line);
  temp = line.substr(0,2);
  while(temp == "##"){
    getline (myfile, line);
    temp = line.substr(0,2);
  }

  // Process header line.
  line = line.substr(1,line.size());
    
  // Split on tab.
  header = splitTab(line);

  // Loop over the rest of the file.
  while ( getline (myfile,line) ){
    linesv.push_back(line);
  }

  myfile.close();

//  Rcout << "linesv size: " << linesv.size() << "\n";

  CharacterMatrix body( linesv.size(), header.size() );

  for(i=0; i<linesv.size(); i++){
    tempv = splitTab(linesv[i]);
    for(j = 0; j < header.size(); j++){
      body(i, j) = tempv[j];
    }
  }

  // Add header to first row.
//  for(j = 0; j < header.size(); j++){
//    body(0, j) = header[j];
//  }
  return body;
}





// [[Rcpp::export]]
Rcpp::DataFrame readVcfBody2(std::string x) {
  std::string line;  // String for reading file into
  std::string temp;  // Temp for subsetting line
  Rcpp::CharacterVector header;  // Vector for the header
  Rcpp::CharacterVector tempv;  // Vector lines after split, before insertion into matrix
  std::vector<std::string> vector_of_lines;  // Vector of lines for the body of the file

  std::ifstream myfile;
  myfile.open (x.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }

  // Iterate past meta lines.
  getline (myfile, line);
  temp = line.substr(0,2);
  while(temp == "##"){
    getline (myfile, line);
    temp = line.substr(0,2);
  }

  // Process header line.
  line = line.substr(1,line.size());

  // Split on tab.
  header = strsplit(line);
  
  // Loop over the rest of the file.
  while ( getline (myfile,line) ){
    vector_of_lines.push_back(line);
  }

  myfile.close();

  // Create matrix for output
  CharacterMatrix body( vector_of_lines.size(), header.size() );
  body.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(), header);

  // Process body of file
  for(int i=0; i<vector_of_lines.size(); i++){
    tempv = strsplit(vector_of_lines[i]);
    for(int j = 0; j < header.size(); j++){
      body(i, j) = tempv[j];
    }
  }

  return body;
}




