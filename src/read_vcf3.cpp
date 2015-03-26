#include <Rcpp.h>
#include <fstream>
#include <iostream>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>

using namespace Rcpp;

const int nreport = 1000;


// [[Rcpp::export]]
Rcpp::NumericVector vcf_stats(std::string x) {
  Rcpp::NumericVector stats(4);
  stats.names() = Rcpp::StringVector::create("meta", "header", "variants", "columns");

  std::string line;  // String for reading file into
  long int i = 0;

  std::ifstream myfile;
  myfile.open (x.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }
  
  // Loop over the file.
  while ( getline (myfile,line) ){
    Rcpp::checkUserInterrupt();
    if(line[0] == '#' && line[1] == '#'){
      stats[0]++;
    } else if (line[0] == '#'){
      stats[1] = stats[0] + 1;
    } else {
      stats[2]++;
    }
    i++;
    
    if( i % nreport == 0){
      Rcout << "\rProcessed line: " << i;
    }
  }
  myfile.close();

  Rcout << "\rProcessed line: " << i;
  Rcout << "\nAll lines processed.\n";

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
Rcpp::StringVector vcf_meta(std::string x, Rcpp::NumericVector stats) {
//Rcpp::List vcf_meta(std::string x, Rcpp::NumericVector stats) {
  // stats consists of elements ("meta", "header", "variants", "columns");
//  Rcpp::List meta(stats[0]);
  Rcpp::StringVector meta(stats[0]);
  std::string line;  // String for reading file into
  
  std::ifstream myfile;
  myfile.open (x.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }
  
  // Loop over the file.
  int i = 0;
  while ( i < stats[0] ){
    Rcpp::checkUserInterrupt();
//    Rcout << i << "\n";
    getline (myfile,line);
    meta(i) = line;
    i++;
    
    if( i % nreport == 0){
      Rcout << "\rProcessed meta line: " << i;
    }

  }
  myfile.close();
  Rcout << "\rProcessed meta line: " << i;
  Rcout << "\nMeta lines processed.\n";

  return meta;
}



std::vector < std::string > tabsplit(std::string line, int elements){
  std::vector < std::string > stringv(elements);
//  Rcout << "Genotype: " << line << "\n";

  int start=0;
  int j = 0;
//  char c;
  for(int i=1; i<line.size(); i++){
//    c = line[i];
    if( line[i] == '\t'){
      std::string temp = line.substr(start, i - start);
//      Rcout << "  i: " << i << ", temp: " << temp << "\n";
      stringv[j] = temp;
      j++;
      start = i+1;
      i = i+1;
    }
  }

  // Handle last element.
  std::string temp = line.substr(start, line.size());
//  Rcout << "  temp: " << temp << "\n";
//  stringv.push_back(temp);
  stringv[j] = temp;
  
//  for(j=0; j<stringv.size(); j++){ Rcout << "j: " << j << " " << stringv[j] << "\n"; }
  
  return stringv;
}



// [[Rcpp::export]]
Rcpp::DataFrame vcf_body_old1(std::string x, Rcpp::NumericVector stats) {
//  Rcpp::StringMatrix body(stats[2], stats[3]);
  Rcpp::CharacterMatrix body(stats[2], stats[3]);
  
  std::string line;  // String for reading file into
  
  std::ifstream myfile;
  myfile.open (x.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }

  // Loop over the file.
  int i = 0;
  int j = 0;
  while ( i < stats[0] ){
    getline (myfile,line);
    i++;
  }
  // Get header.
  getline (myfile,line);
  std::string header = line;
  
  // Get body.
  i = 0;
  while ( getline (myfile,line) ){
    Rcpp::checkUserInterrupt();
//    body(i, _) = tabsplit(header);
    std::vector < std::string > temps = tabsplit(line, stats[3]);
    for(j=0; j<stats[3]; j++){
      body(i, j) = temps[j];
    }
    i++;
  }

  myfile.close();

  Rcpp::DataFrame df1(body);
  std::vector < std::string > temps = tabsplit(header, stats[3]);
  temps[0].erase(0,1);
//  df1.names() = Rcpp::StringVector::create(temps);
  df1.names() = temps;

  return df1;
}



// [[Rcpp::export]]
Rcpp::DataFrame vcf_body(std::string x, Rcpp::NumericVector stats) {
//  Rcpp::StringVector   chrom(stats[2]);
  Rcpp::CharacterVector   chrom(stats[2]);
  Rcpp::IntegerVector  pos(stats[2]);
  Rcpp::StringVector   id(stats[2]);
  Rcpp::StringVector   ref(stats[2]);
  Rcpp::StringVector   alt(stats[2]);
  Rcpp::NumericVector  qual(stats[2]);
  Rcpp::StringVector   filter(stats[2]);
  Rcpp::StringVector   info(stats[2]);
  
  Rcpp::CharacterMatrix gt(stats[2], stats[3] - 8);
  
  std::string line;  // String for reading file into

  // Open file.
  std::ifstream myfile;
  myfile.open (x.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }

  // Iterate past meta.
  int i = 0;
  int j = 0;
  while ( i < stats[0] ){
    getline (myfile,line);
    i++;
  }
  
  // Get header.
  getline (myfile,line);
  std::string header = line;
  
  // Get body.
  i = 0;
  char buffer [50];
  while ( getline (myfile,line) ){
    Rcpp::checkUserInterrupt();
    std::vector < std::string > temps = tabsplit(line, stats[3]);

    if(temps[0] == "."){
      chrom[i] = NA_STRING;
    } else {
      chrom[i] = temps[0];
    }
    if(temps[1] == "."){
      pos[i] = NA_INTEGER;
    } else {
        pos[i] = atoi(temps[1].c_str());
    }
    if(temps[2] == "."){
      id[i] = NA_STRING;
    } else {
      id[i] = temps[2];
    }
    if(temps[3] == "."){
      ref[i] = NA_STRING;
    } else {
      ref[i] = temps[3];
    }
    if(temps[4] == "."){
      alt[i] = NA_STRING;
    } else {
      alt[i] = temps[4];
    }
    if(temps[5] == "."){
      qual[i] = NA_REAL;
    } else {
      qual[i] = atof(temps[5].c_str());
    }
    if(temps[6] == "."){
      filter[i] = NA_STRING;
    } else {
      filter[i] = temps[6];
    }
    if(temps[7] == "."){
      info[i] = NA_STRING;
    } else {
      info[i] = temps[7];
    }

    for(j=8; j<stats[3]; j++){
      gt(i, j-8) = temps[j];
//      body(i, j-8) = temps[j];
    }
    i++;
    
    if( i % nreport == 0){
      Rcout << "\rProcessed variant: " << i;
    }

  }
  Rcout << "\rProcessed variant: " << i;
  Rcout << "\nAll variants processed\n";
      
  myfile.close();

  Rcpp::DataFrame df1 = Rcpp::DataFrame::create(
    _["CHROM"]= chrom,
    _["POS"]= pos,
    _["ID"] = id,
    _["REF"] = ref,
    _["ALT"] = alt,
    _["QUAL"] = qual,
    _["FILTER"] = filter,
    _["INFO"] = info,
    gt);

  std::vector < std::string > temps = tabsplit(header, stats[3]);
  temps[0].erase(0,1);
  df1.names() = temps;
  
  return df1;
}




// [[Rcpp::export]]
int read_to_line(std::string x) {
  std::string line;  // String for reading file into
  long int i = 0;

  std::ifstream myfile;
  myfile.open (x.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }
  
  // Loop over the file.
  while ( getline (myfile,line) ){
    Rcpp::checkUserInterrupt();
//    Rcout << line << "\n";
    i++;
  }

  myfile.close();
  return i;
}


// Function requires boost for gz compression.
// [[Rcpp::depends(BH)]]
// [[Rcpp::export]]
int read_gz_to_line(std::string x) {
  std::string line;  // String for reading file into
  long int i = 0;

//  std::ifstream myfile(x.c_str(), std::ios_base::in | std::ios_base::binary);
//  boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
//  inbuf.push(boost::iostreams::gzip_decompressor());
//  inbuf.push(file);

//    std::ifstream myfile(x.c_str(), std::ios_base::in | std::ios_base::binary);
    std::ifstream myfile("hello.gz", std::ios_base::in | std::ios_base::binary);
//    boost::iostreams::filtering_streambuf<boost::iostreams::input> in;
//    in.push(boost::iostreams::gzip_decompressor());
//    in.push(myfile);
//    boost::iostreams::copy(in, Rcout);



   //Read from the first command line argument, assume it's gzipped

//  std::ifstream myfile;
//  myfile.open (x.c_str(), std::ios::in);

  if (!myfile.is_open()){
    Rcout << "Unable to open file";
  }
  
  // Loop over the file.
  while ( getline (myfile,line) ){
    Rcout << line << "\n";
    i++;
  }

  myfile.close();
  return i;
}




// [[Rcpp::export]]
Rcpp::CharacterMatrix ram_test(int nrow=1, int ncol=1) {
  Rcpp::CharacterMatrix gt(nrow, ncol);
//Rcpp::StringMatrix ram_test() {
//  Rcpp::StringMatrix gt(1, 1);
  
  return(gt);
}

