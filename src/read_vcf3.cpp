#include <Rcpp.h>
#include <fstream>
#include <zlib.h>
#include "vcfRCommon.h"

// #include <iostream>
//using namespace Rcpp;



// Number of records to report progress at.
const int nreport = 1000;

/* Size of the block of memory to use for reading. */
#define LENGTH 0x1000


/*  Helper functions */

void stat_line(Rcpp::NumericVector stats, std::string line){
  if(line[0] == '#' && line[1] == '#'){
    // Meta
    stats(0)++;
  } else if (line[0] == '#' && line[1] != '#'){
    // Header
    stats(1) = stats(0) + 1;
    std::vector < std::string > col_vec;
    char col_split = '\t'; // Must be single quotes!
    vcfRCommon::strsplit(line, col_vec, col_split);
    stats(3) = col_vec.size();
  } else {
    // Variant
    stats(2)++;
  }
}



/*  Single pass of vcf file to get statistics */

// [[Rcpp::export]]
Rcpp::NumericVector vcf_stats_gz(std::string x) {
  Rcpp::NumericVector stats(4);  // 4 elements, all zero.  Zero is default.
  stats.names() = Rcpp::StringVector::create("meta", "header", "variants", "columns");
  
  gzFile file;
  file = gzopen (x.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << x << " failed: " << strerror (errno) << ".\n";
    return stats;
  }

  // Scroll through buffers.
  std::string lastline = "";
  while (1) {
    Rcpp::checkUserInterrupt();
    int err;
    int bytes_read;
    char buffer[LENGTH];
    bytes_read = gzread (file, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0';

    std::string mystring(reinterpret_cast<char*>(buffer));  // Recast buffer as a string.
    mystring = lastline + mystring;
    std::vector < std::string > svec;  // Initialize vector of strings for parsed buffer.
    
    char split = '\n'; // Must be single quotes!
    vcfRCommon::strsplit(mystring, svec, split);
        
    // Scroll through lines derived from the buffer.
    for(int i=0; i < svec.size() - 1; i++){
      stat_line(stats, svec[i]);
    }
    // Manage the last line.
    lastline = svec[svec.size() - 1];

    // Check for EOF or errors.
    if (bytes_read < LENGTH - 1) {
      if (gzeof (file)) {
        break;
      }
      else {
        const char * error_string;
        error_string = gzerror (file, & err);
        if (err) {
          Rcpp::Rcerr << "Error: " << error_string << ".\n";
          return stats;
        }
      }
    }
  }
  gzclose (file);

  return stats;
}




/*  Read vcf meta region  */

// [[Rcpp::export]]
Rcpp::StringVector read_meta_gz(std::string x, Rcpp::NumericVector stats, int verbose) {
  // Read in the meta lines.
  // stats consists of elements ("meta", "header", "variants", "columns");

  Rcpp::StringVector meta(stats[0]);
  std::string line;  // String for reading file into
  int meta_row = 0;

  gzFile file;
  file = gzopen (x.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << x << " failed: " << strerror (errno) << ".\n";
    return Rcpp::StringVector(1);
  }

  // Scroll through buffers.
  std::string lastline = "";
  while (1) {
    Rcpp::checkUserInterrupt();
    int err;
    int bytes_read;
    char buffer[LENGTH];
    bytes_read = gzread (file, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0';

    std::string mystring(reinterpret_cast<char*>(buffer));  // Recast buffer as a string.
    mystring = lastline + mystring;
    std::vector < std::string > svec;  // Initialize vector of strings for parsed buffer.
    char split = '\n'; // Must be single quotes!
    vcfRCommon::strsplit(mystring, svec, split);


    int i = 0;
    while(meta_row < stats(0) && i < svec.size() - 1){
      meta(meta_row) = svec[i];
      meta_row++;
      i++;
      if(meta_row % nreport == 0 && verbose == 1){
        Rcpp::Rcout << "\rMeta line " << meta_row << " read in.";
      }
    }
    lastline = svec[svec.size() - 1];


    // Check for EOF or errors.
    if (bytes_read < LENGTH - 1) {
      if (gzeof (file)) {
        break;
      }
      else {
        const char * error_string;
        error_string = gzerror (file, & err);
        if (err) {
          Rcpp::Rcerr << "Error: " << error_string << ".\n";
          return Rcpp::StringVector(1);
        }
      }
    }
  }
  gzclose (file);

  if(verbose == 1){
    Rcpp::Rcout << "\rMeta line " << meta_row << " read in.";
    Rcpp::Rcout << "\nAll meta lines processed.\n";
  }

  return meta;
}




/*  Read in the body of the vcf file  */

void proc_body_line(Rcpp::CharacterMatrix gt, int var_num, std::string myline){
  char split = '\t'; // Must be single quotes!
  std::vector < std::string > data_vec;
  
  vcfRCommon::strsplit(myline, data_vec, split);

  for(int i = 0; i < data_vec.size(); i++){
    if(data_vec[i] == "."){
      gt(var_num, i) = NA_STRING;
    } else if(data_vec[i] == "./."){
      gt(var_num, i) = NA_STRING;
    } else {
      gt(var_num, i) = data_vec[i];
    }
  }
}


// [[Rcpp::export]]
Rcpp::CharacterMatrix read_body_gz(std::string x, Rcpp::NumericVector stats, int verbose = 1) {
//Rcpp::DataFrame read_body_gz(std::string x, Rcpp::NumericVector stats, int verbose = 1) {
  // Read in the fixed and genotype portion of the file.

  // Stats contains:
  // "meta", "header", "variants", "columns"

  // Matrix for body data.
  Rcpp::CharacterMatrix gt(stats[2], stats[3]);

  gzFile file;
  file = gzopen (x.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << x << " failed: " << strerror (errno) << ".\n";
//    return Rcpp::StringVector(1);
    return Rcpp::CharacterMatrix(1);
  }

  // Scroll through buffers.
  std::string lastline = "";
  std::vector<std::string> header_vec;
  int var_num = 0;
  while (1) {
    Rcpp::checkUserInterrupt();
    int err;
    int bytes_read;
    char buffer[LENGTH];
    bytes_read = gzread (file, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0';

    std::string mystring(reinterpret_cast<char*>(buffer));  // Recast buffer as a string.
    mystring = lastline + mystring;
    std::vector < std::string > svec;  // Initialize vector of strings for parsed buffer.
    char split = '\n'; // Must be single quotes!
    vcfRCommon::strsplit(mystring, svec, split);

    for(int i = 0; i < svec.size() - 1; i++){
      if(svec[i][0] == '#' && svec[i][1] == '#'){
        // Meta line, ignore.
      } else if(svec[i][0] == '#' && svec[i][1] != '#'){
        // Process header
        char header_split = '\t';
        vcfRCommon::strsplit(svec[i], header_vec, header_split);
      } else {
        // Variant line.
        proc_body_line(gt, var_num, svec[i]);
        var_num++;
        
        if(var_num % nreport == 0 && verbose == 1){
          Rcpp::Rcout << "\rProcessed variant " << var_num;
        }
      }   
    }
    lastline = svec[svec.size() - 1];


    // Check for EOF or errors.
    if (bytes_read < LENGTH - 1) {
      if (gzeof (file)) {
        break;
      }
      else {
        const char * error_string;
        error_string = gzerror (file, & err);
        if (err) {
          Rcpp::Rcerr << "Error: " << error_string << ".\n";
//          return Rcpp::StringVector(1);
          return Rcpp::CharacterMatrix(1);
        }
      }
    }
  }
  gzclose (file);
  
  header_vec[0] = "CHROM";

  gt.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(), header_vec);

  if(verbose == 1){
    Rcpp::Rcout << "\rProcessed variant: " << var_num;
    Rcpp::Rcout << "\nAll variants processed\n";
  }

//  Rcpp::DataFrame df1 = Rcpp::DataFrame::create(gt);
//  Rcpp::DataFrame df1(gt);
//  df1.names() = header_vec;
//  if(verbose == 1){
//    Rcpp::Rcout << "Rcpp::DataFrame created.\n";
//  }
  
  return gt;
//  return df1;
}



/*  Memory test  */



Rcpp::StringMatrix DataFrame_to_StringMatrix( Rcpp::DataFrame df ){
  Rcpp::StringVector sv = df(0);
  Rcpp::StringMatrix sm(sv.size(), df.size());
  
  sm.attr("col.names") = df.attr("col.names");
  sm.attr("row.names") = df.attr("row.names");

  for(int i=0; i < df.size(); i++){
    sv = df(i);
    for(int j=0; j < sv.size(); j++){
      sm(j, i) = sv(j);
    }
  }

  return sm;
}


/*  Write vcf body  */

// [[Rcpp::export]]
void write_vcf_body( Rcpp::DataFrame fix, Rcpp::DataFrame gt, std::string filename , int mask=0 ) {
//int write_vcf_body( Rcpp::DataFrame fix, Rcpp::DataFrame gt, std::string filename , int mask=0 ) {

  // fix DataFrame
  Rcpp::StringVector chrom  = fix["CHROM"];
  Rcpp::StringVector pos    = fix["POS"];
  Rcpp::StringVector id     = fix["ID"];
  Rcpp::StringVector ref    = fix["REF"];
  Rcpp::StringVector alt    = fix["ALT"];
  Rcpp::StringVector qual   = fix["QUAL"];
  Rcpp::StringVector filter = fix["FILTER"];
  Rcpp::StringVector info   = fix["INFO"];

  // gt DataFrame
  Rcpp::StringMatrix gt_cm = DataFrame_to_StringMatrix(gt);
  Rcpp::StringVector column_names(gt.size());
  column_names = gt.attr("names");
//  column_names = gt_cm.attr("col.names");
//  delete gt;
  
  int i = 0;
  int j = 0;

  // Uncompressed.
  std::ofstream myfile;
  myfile.open (filename.c_str(), std::ios::out | std::ios::app | std::ios::binary);
  
//  gzFile *fi = (gzFile *)gzopen("file.gz","wb");
  

  for(i=0; i<chrom.size(); i++){
    Rcpp::checkUserInterrupt();
    if(mask == 1 && filter(i) == "PASS" ){
      // Don't print variant.
    } else {
      myfile << chrom(i);
      myfile << "\t";
      myfile << pos(i);
      myfile << "\t";
      if(id(i) == NA_STRING){
        myfile << ".";
        myfile << "\t";
      } else {
        myfile << id(i);
        myfile << "\t";
      }
      myfile << ref(i);
      myfile << "\t";
      myfile << alt(i);
      myfile << "\t";
      if(qual(i) == NA_STRING){
        myfile << ".";
        myfile << "\t";
      } else {
        myfile << qual(i);
        myfile << "\t";
      }
      if(filter(i) == NA_STRING){
        myfile << ".";
        myfile << "\t";
      } else {
        myfile << filter(i);
        myfile << "\t";
      }
      if(info(i) == NA_STRING){
        myfile << ".";
        myfile << "\t";
      } else {
        myfile << info(i);
      }
      
      // gt region.
      myfile << "\t";
      myfile << gt_cm(i, 0);
      for(j=1; j<column_names.size(); j++){
        myfile << "\t";
        myfile << gt_cm(i, j);
      }

      myfile << "\n";
    }
  }

  myfile.close();
  
  return;
}



// [[Rcpp::export]]
void write_vcf_body_gz( Rcpp::DataFrame fix, Rcpp::DataFrame gt, std::string filename , int mask=0 ) {
  // http://stackoverflow.com/a/5649224
  
  // fix DataFrame
  Rcpp::StringVector chrom  = fix["CHROM"];
  Rcpp::StringVector pos    = fix["POS"];
  Rcpp::StringVector id     = fix["ID"];
  Rcpp::StringVector ref    = fix["REF"];
  Rcpp::StringVector alt    = fix["ALT"];
  Rcpp::StringVector qual   = fix["QUAL"];
  Rcpp::StringVector filter = fix["FILTER"];
  Rcpp::StringVector info   = fix["INFO"];
  
  // gt DataFrame
  Rcpp::StringMatrix gt_cm = DataFrame_to_StringMatrix(gt);
  Rcpp::StringVector column_names(gt.size());
  column_names = gt.attr("names");
  
  int i = 0;
  int j = 0;
  
  
  gzFile fi = gzopen( filename.c_str(), "ab" );
//  gzFile *fi = (gzFile *)gzopen( filename.c_str(), "ab" );
//  gzFile *fi = (gzFile *)gzopen(filename.c_str(),"abw");
  for(i=0; i<chrom.size(); i++){
    Rcpp::checkUserInterrupt();
    if(mask == 1 && filter(i) != "PASS" ){
      // Don't print variant.
    } else {
      std::string tmpstring;
      tmpstring = chrom(i);
      tmpstring = tmpstring + "\t" + pos(i) + "\t";
      if(id(i) == NA_STRING){
        tmpstring = tmpstring + ".";
      } else {
        tmpstring = tmpstring + id(i);
      }
      tmpstring = tmpstring + "\t" + ref(i) + "\t" + alt(i) + "\t";
      if(qual(i) == NA_STRING){
        tmpstring = tmpstring + "." + "\t";
      } else {
        tmpstring = tmpstring + qual(i) + "\t";
      }
      if(filter(i) == NA_STRING){
        tmpstring = tmpstring + "." + "\t";
      } else {
        tmpstring = tmpstring + filter(i) + "\t";
      }
      tmpstring = tmpstring + info(i);

      // gt portion
      for(j=0; j<column_names.size(); j++){
        if(gt_cm(i, j) == NA_STRING){
          tmpstring = tmpstring + "\t" + "./.";
        } else {
          tmpstring = tmpstring + "\t" + gt_cm(i, j);
        }
      }


//      gzwrite(fi,"my decompressed data",strlen("my decompressed data"));
//      gzwrite(fi,"\n",strlen("\n"));
//      std::string tmpstring = "test string\n";
      gzwrite(fi, (char *)tmpstring.c_str(), tmpstring.size());
      
      gzwrite(fi,"\n",strlen("\n"));
    }
  }
  gzclose(fi);
  
  return;
}




// [[Rcpp::export]]
void write_fasta( Rcpp::CharacterVector seq,
                  std::string seqname, 
                  std::string filename, 
                  int rowlength=80,
                  int verbose=1) {
//  rowlength=rowlength-1;
  FILE * pFile;
//  pFile=fopen(filename.c_str(),"wt");
  pFile=fopen(filename.c_str(),"at");
  int i=0;

  if(verbose == 1){
    Rcpp::Rcout << "Processing sample: " << seqname << "\n";
  }

  putc ('>' , pFile);
  for(i=0; i<seqname.size(); i++){
    putc (seqname[i] , pFile);
  }
  putc ('\n' , pFile);

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



