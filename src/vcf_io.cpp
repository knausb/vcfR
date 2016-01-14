#include <Rcpp.h>
#include <fstream>
#include <zlib.h>
#include "vcfRCommon.h"


// Number of records to report progress at.
const int nreport = 1000;

/* Size of the block of memory to use for reading. */
#define LENGTH 0x1000


/*  Helper functions */

/* 
Called by vcf_stats_gz.
Processes lines from vcf files.
Counts meta (^##), header (^#C), columns in the header and remaining lines.
*/
void stat_line(Rcpp::NumericVector stats, std::string line){
  if(line[0] == '#' && line[1] == '#'){
    // Meta
    stats(0)++;
  } else if (line[0] == '#' && line[1] == 'C'){
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
      if ( gzeof (file) ) {
        if( stats(3) == 0 ){
          // Count columns from last line.
          std::vector < std::string > col_vec;
          char col_split = '\t'; // Must be single quotes!
          vcfRCommon::strsplit(svec[0], col_vec, col_split);
          stats(3) = col_vec.size();
        }
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

/*  Helper function to process one line  */
void proc_body_line(Rcpp::CharacterMatrix gt,
                    int var_num,
                    std::string myline,
                    Rcpp::IntegerVector cols){
  
  char split = '\t'; // Must be single quotes!
  std::vector < std::string > data_vec;
  
  vcfRCommon::strsplit(myline, data_vec, split);

//  for(int i = 0; i < data_vec.size(); i++){
  for(int i = 0; i < cols.size(); i++){
//    if(data_vec[i] == "."){
    if(data_vec[ cols[i] ] == "."){
      gt(var_num, i) = NA_STRING;
//    } else if(data_vec[i] == "./."){
//    } else if( data_vec[ i ][0] == '.' & data_vec[ i ][2] == '.' ){
    } else if( data_vec[ cols[i] ][0] == '.' & data_vec[ cols[i] ][2] == '.' & data_vec[ cols[i] ].size() == 3 ){
      gt(var_num, i) = NA_STRING;
    } else {
      gt(var_num, i) = data_vec[ cols[i] ];
    }
  }
}

/* 
 Read in the fixed and genotype portion of the file.
 x is the file name.
 Stats is a vector containing information about the file contents (see below).
 varbose indicates whether verbose output should be generated.
 
 Stats contains:
 "meta", "header", "variants", "columns"
 
*/
// [[Rcpp::export]]
Rcpp::CharacterMatrix read_body_gz(std::string x,
                                   Rcpp::NumericVector stats,
                                   int nrows,
                                   int skip,
                                   Rcpp::IntegerVector cols,
                                   int verbose = 1) {

  // Sort the column numbers.
  cols.sort();
  
  // Validate that column numbers have been specified.
  if(cols[0] == 0){
    Rcpp::Rcerr << "User must specify which (positive integer) columns to extract from the file.\n";
    Rcpp::StringMatrix mymatrix(1,1);
    mymatrix(0,0) = NA_STRING;
    return mymatrix;
  }
  cols = cols - 1; // R is 1-based
  
  // Initialize matrix for body data.
//  Rcpp::CharacterMatrix gt(stats[2], stats[3]);
  Rcpp::CharacterMatrix gt( stats[2], cols.size() );


  // Create filehandle and open.
  gzFile file;
  file = gzopen (x.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << x << " failed: " << strerror (errno) << ".\n";
    return Rcpp::CharacterMatrix(1);
  }


  // Because the last line may be incomplete,
  // We'll typically omit it from processing and
  // concatenate it to the first line.
  // But first we'll have to initialize it.
  std::string lastline = "";
  
  // String vector to store the header (^#CHROM...).
  std::vector<std::string> header_vec;
  
  // variant counter.  
  int var_num = 0;

  // Scroll through buffers.
  while (1) {
    Rcpp::checkUserInterrupt();
    int err;

    // Slurp in a buffer.
    int bytes_read;
    char buffer[LENGTH];
    bytes_read = gzread (file, buffer, LENGTH - 1);
    buffer[bytes_read] = '\0'; // Terminate the buffer.

    std::string mystring(reinterpret_cast<char*>(buffer));  // Recast buffer as a string.
    mystring = lastline + mystring; // Concatenate last line to the buffer
    
    // Delimit into lines.
    std::vector < std::string > svec;  // Initialize vector of strings for parsed buffer.
    char split = '\n'; // Must be single quotes!
    vcfRCommon::strsplit(mystring, svec, split);
    
    /* 
    svec should now contain a vector of strings,
    one string for each line
    where the last line may be incomplete.
    We can now process each line except the last.
    */

    // Scroll through lines.
    for(int i = 0; i < svec.size() - 1; i++){
      if(svec[i][0] == '#' && svec[i][1] == '#'){
        // Meta line, ignore.
      } else if(svec[i][0] == '#' && svec[i][1] == 'C'){
        // Process header.
        char header_split = '\t';
        vcfRCommon::strsplit(svec[i], header_vec, header_split);
        
        // Subset the header to select columns.
        std::vector<std::string> header_vec2( cols.size() );
        for(int j=0; j<cols.size(); j++){
          header_vec2[j] = header_vec[ cols[j] ];
        }
        header_vec = header_vec2;
        
      } else {
        // Variant line.

//        Rcpp::Rcout << "var_num: " << var_num << "\n";
//        Rcpp::Rcout << "skip: " << skip << "\n";
        if( var_num >= skip ){
          proc_body_line(gt, var_num, svec[i], cols);
        }
        var_num++;        

        if(var_num % nreport == 0 && verbose == 1){
          Rcpp::Rcout << "\rProcessed variant " << var_num;
        }
      }   
    }
    // Keep the last line so we can append it to 
    //the beginning of the next buffer
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
    
  // Return to top of loop and process another buffer.
  }
  
  // Close filehandle.
  gzclose (file);


  if( stats[1] == 0 ){
    if( verbose == 1 ){
      Rcpp::Rcout << "Warning: no header information was found! Data contains no sample names!\n";
    }
  } else {
    header_vec[0] = "CHROM";
    if( header_vec.size() == gt.ncol() ){
      gt.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(), header_vec);
    } else {
      if( verbose == 1 ){
        Rcpp::Rcout << "Warning: no header information found!\n";
      }
    }
  }

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
}





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
void write_vcf_body( Rcpp::CharacterMatrix fix, Rcpp::CharacterMatrix gt, std::string filename , int mask=0 ) {
//void write_vcf_body( Rcpp::DataFrame fix, Rcpp::DataFrame gt, std::string filename , int mask=0 ) {
//int write_vcf_body( Rcpp::DataFrame fix, Rcpp::DataFrame gt, std::string filename , int mask=0 ) {
  // http://stackoverflow.com/a/5649224
  
  int i = 0; // Rows
  int j = 0; // Columns
  std::string tmpstring;  // Assemble each line before writing
  
  
  gzFile fi = gzopen( filename.c_str(), "ab" );

  // In order for APPEND=TRUE to work the header
  // should not be printed here.

  // Manage header.
/*  Rcpp::List matrix_names = fix.attr("dimnames");
  Rcpp::StringVector head_names = matrix_names(1);
  tmpstring = "#" + head_names(0);
  for(i = 1; i < head_names.size(); i++){
    tmpstring = tmpstring + "\t" + head_names(i);
  }

  matrix_names = gt.attr("dimnames");
  head_names = matrix_names(1);
  for(i = 0; i < head_names.size(); i++){
    tmpstring = tmpstring + "\t" + head_names(i);
  }
*/
  // Write header.
//  gzwrite(fi, (char *)tmpstring.c_str(), tmpstring.size());
//  gzwrite(fi,"\n",strlen("\n"));

  
  
  // Manage body
  for(i = 0; i < fix.nrow(); i++){
    Rcpp::checkUserInterrupt();
    if(mask == 1 && fix(i,6) != "PASS" ){
      // Don't print variant.
    } else {
      j = 0;
      tmpstring = fix(i,j);
      for(j = 1; j < fix.ncol(); j++){
        if(fix(i,j) == NA_STRING){
          tmpstring = tmpstring + "\t" + ".";
        } else {
          tmpstring = tmpstring + "\t" + fix(i,j);
        }
      }

      // gt portion
      for(j = 0; j < gt.ncol(); j++){
        if(gt(i, j) == NA_STRING){
          tmpstring = tmpstring + "\t" + "./.";
        } else {
          tmpstring = tmpstring + "\t" + gt(i, j);
        }
      }

      gzwrite(fi, (char *)tmpstring.c_str(), tmpstring.size());
      gzwrite(fi,"\n",strlen("\n"));
    }
  }
  gzclose(fi);
  
  return;
}



/* Write data to fasta file */

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



