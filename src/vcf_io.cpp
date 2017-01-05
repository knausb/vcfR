#include <Rcpp.h>
#include <fstream>
#include <zlib.h>
#include "vcfRCommon.h"


// Number of records to report progress at.
const int nreport = 1000;

/* Size of the block of memory to use for reading. */
#define LENGTH 0x1000 // hexadecimel for 4096.


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
    unsigned int i = 0;
    for(i=0; i < svec.size() - 1; i++){
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


    unsigned int i = 0;
    while(meta_row < stats(0) && i < svec.size() - 1){
      
      // Check and remove carriage returns (Windows).
      if( svec[i][ svec[i].size()-1] == '\r' ){
        svec[i].erase( svec[i].size() - 1 );
      }

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

/*  
Helper function to process one line.
gt is a pointer to the data matrix output by read_vcf.
var_num
myline is a tab delimited string row from the orriginal vcf file.
cols is a vector indicating which columns to select from the file data.
*/
void proc_body_line(Rcpp::CharacterMatrix gt,
                    int var_num,
                    std::string myline,
                    Rcpp::IntegerVector cols,
                    int convertNA){
  
  char split = '\t'; // Must be single quotes!
  std::vector < std::string > data_vec;
  
  vcfRCommon::strsplit(myline, data_vec, split);

  for(int i = 0; i < cols.size(); i++){
    if( convertNA == 1 ){
      // The VCF specification encodes missing data as ".".
      // Missing genotypes may be encoded as "./.", ".|.", etc.
      // Here we convert is to NA.
      if( data_vec[ cols[i] ] == "." ){
        gt(var_num, i) = NA_STRING;
      } else {
        // Possible genotype missing.
        std::vector < std::string > allele_vec;
        int unphased_as_na = 0; // 0 == FALSE
        std::string my_string;
//        if( data_vec[ cols[i] ] == NA_STRING ){
//          my_string = ".";
//        } else {
          my_string = data_vec[ cols[i] ];
//        }
        
        vcfRCommon::gtsplit( my_string, allele_vec, unphased_as_na );
        int gtNA = 1;
        unsigned int k = 0;
        for( k = 0; k < allele_vec.size(); k++ ){
//            Rcpp::Rcout << "allele_vec[k]: " << allele_vec[k] << "\n";
          if( allele_vec[k] != "." ){ gtNA = 0; }
        }
        if( gtNA == 1 ){
          gt(var_num, i) = NA_STRING;
        } else {
          gt(var_num, i) = data_vec[ cols[i] ];
        }
      }
//        int gtNA = 1;
//        for(int j = 0; j < data_vec[ cols[i] ].size(); j++){
//          if( data_vec[ cols[i] ][j] != '.' ){
//            gtNA = 0;
//          }
//          j++; // Every other character should be a delimiter.
//        }
//        if( gtNA == 1 ){
//          gt(var_num, i) = NA_STRING;
//        } else {
//          gt(var_num, i) = data_vec[ cols[i] ];
//        }
//      }
//      } else if( data_vec[ cols[i] ][0] == '.' & data_vec[ cols[i] ][2] == '.' & 
//               data_vec[ cols[i] ].size() == 3 & convertNA == 1 ){
      // We can also convert diploid genotypes where both alleles are "." to NA.
//      gt(var_num, i) = NA_STRING;
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
                                   int nrows = -1,
                                   int skip = 0,
                                   Rcpp::IntegerVector cols = 0,
                                   int convertNA = 1,
                                   int verbose = 1) {

  // NA matrix for unexpected results.
  Rcpp::StringMatrix na_matrix(1,1);
  na_matrix(0,0) = NA_STRING;
  
  
  /*
   * Manage cols vector.
   * The first eight (1-based) columns are mandatory.
   * We can ensure they are there by adding them,
   * sorting and removing adjacent non-identical values.
   */
//  for( int i=9; i >= 1; i-- ){
  for( int i=8; i >= 1; i-- ){
    cols.push_front(i);
  }
  cols.sort();

  // Remove duplicate values using a set.
  std::set<int> s( cols.begin(), cols.end() );
  cols.assign( s.begin(), s.end() );

  cols = cols - 1; // R is 1-based, C is 0-based.

  
  // Initialize matrix for body data.
  // old: Rcpp::CharacterMatrix gt(stats[2], stats[3]);
  int row_num = 0;
  

  if( ( nrows == -1 ) & ( skip == 0 ) ){
    nrows = stats[2];
  } else if ( ( nrows != -1 ) & ( skip == 0 ) ){
    // nrows = nrows;
  } else if ( ( nrows == -1 ) & ( skip > 0) ){
    nrows = stats[2] - skip;
  } else if ( ( nrows != -1 ) & ( skip > 0) ){
    // nrows = nrows;
  } else {
    Rcpp::Rcerr << "failed to calculate return matrix geometry.";
    return na_matrix;
  }
  Rcpp::CharacterMatrix gt( nrows, cols.size() );
  
//  if ( nrows > -1 & skip == 0 ){
//    row_num = nrows;
//  } else if ( nrows == -1 & skip > 0 ){
//    row_num = stats[2] - skip;
//  } else {
//    row_num = stats[2];    
//  }
//  Rcpp::CharacterMatrix gt( row_num, cols.size() );

  row_num = 0;
  
  if( verbose == 1 ){
    Rcpp::Rcout << "Character matrix gt created.\n";
    Rcpp::Rcout << "Character matrix gt rows: ";  Rcpp::Rcout << gt.rows();
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "Character matrix gt cols: ";  Rcpp::Rcout << gt.cols();
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "skip: ";  Rcpp::Rcout << skip;
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "nrows: ";  Rcpp::Rcout << nrows;
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "row_num: ";  Rcpp::Rcout << row_num;
    Rcpp::Rcout << "\n";
    Rcpp::Rcout << "\n";
  }

  
  // Create filehandle and open.
  gzFile file;
  file = gzopen (x.c_str(), "r");
  if (! file) {
    Rcpp::Rcerr << "gzopen of " << x << " failed: " << strerror (errno) << ".\n";
    return na_matrix;
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
    unsigned int i = 0;
    for(i = 0; i < svec.size() - 1; i++){
      
      // Check and remove carriage returns (Windows).
      if( svec[i][ svec[i].size()-1] == '\r' ){
        svec[i].erase( svec[i].size() - 1 );
      }

      if(svec[i][0] == '#' && svec[i][1] == '#'){
        // Meta line, ignore.
      } else if(svec[i][0] == '#' && svec[i][1] == 'C'){
        // Process header.
//        Rcpp::Rcout << svec[i].substr(0,40) << "\n\n";
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

        if ( ( var_num >= skip ) & ( row_num < nrows ) ){
          proc_body_line(gt, row_num, svec[i], cols, convertNA);
          row_num++; // Return matrix row number.
        }
        var_num++; // Input row number.


        if(var_num % nreport == 0 && verbose == 1){
          Rcpp::Rcout << "\rProcessed variant " << var_num;
        }
      }   
    }

    
    // Processed all lines of current buffer.
    // Keep the last line so we can append it to 
    //the beginning of the next buffer.
    lastline = svec[svec.size() - 1];

//    Rcpp::Rcout << "line-2:" << svec[svec.size() - 2].substr(0,40) << "|<-\n";
//    Rcpp::Rcout << "line-1:" << svec[svec.size() - 1].substr(0,40) << "|<-\n";
//    Rcpp::Rcout << "\n";

      
    /*
     * We can bail out early if we have read nrows.
     * Before we do we need to check that:
     * 1) we have read in nrows
     * 2) we have processed the header
     * (important when nrows is small)
     * 3) we actually have a line (when buffer ends at the end of a line).
     */
    if( ( row_num >= nrows ) & ( lastline[0] != '#' ) & ( lastline.size() > 0 ) ){
//        Rcpp::Rcout << "\nBreaking!\n";
//        Rcpp::Rcout << "lastline: " << lastline.substr(0,40) << "\n";
        break;
    }


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
          return na_matrix;
        }
      }
    }
    
  // Return to top of loop and process another buffer.
  } // Close while.
  
  // Close filehandle.
  gzclose (file);

//  Rcpp::Rcout << "\n\n>>---<< Made it: file close! >>---<<\n\n";
//  Rcpp::Rcout << "header_vec.size(): " << header_vec.size() << "\n";
  
  
  if( stats[1] == 0 ){
    if( verbose == 1 ){
      Rcpp::Rcout << "Warning: no header information was found! Data contains no sample names!\n";
    }
  } else {

    if( header_vec.size() == (unsigned)gt.ncol() ){
      header_vec[0] = "CHROM";
      gt.attr("dimnames") = Rcpp::List::create(Rcpp::CharacterVector::create(), header_vec);
    } else {
      if( verbose == 1 ){
        Rcpp::Rcout << "Warning: no header information found!\n";
      }
    }
  }

//  Rcpp::Rcout << "\n\n>>---<< Made it! >>---<<\n\n";

    
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
void write_vcf_body( Rcpp::CharacterMatrix fix,
                     Rcpp::CharacterMatrix gt,
                     std::string filename,
                     int mask=0 ) {
  // http://stackoverflow.com/a/5649224
  
//  
int verbose = 0;
//  int verbose = 1;
  
  if( verbose == 1 ){
    Rcpp::Rcout << "Made it into the function!\n";
  }
  
  int i = 0; // Rows
  int j = 0; // Columns
  std::string tmpstring;  // Assemble each line before writing

  // Initialize filehandle.
  gzFile fi;
  
  // Initialize file.
  // Note that gzfile does not tolerate initializing an empty file.
  // Use ofstream instead.
  if ( ! std::ifstream( filename ) ){
    if( verbose == 1 ){
      Rcpp::Rcout << "File does not exist." << std::endl;
    }
    
    std::ofstream myfile;
    myfile.open (filename, std::ios::out | std::ios::binary);
    myfile.close();
    
    // This should make valgrind hang.
    // Or not???
//    fi = gzopen( filename.c_str(), "ab" );
//    gzclose(fi);
  }

  // In order for APPEND=TRUE to work the header
  // should not be printed here.

  if( verbose == 1 ){
    Rcpp::Rcout << "Matrix fix has " << fix.nrow() << " rows (variants).\n";
  }
  
  // Manage body
  if( fix.nrow() >= 1 ){
    if( verbose == 1 ){
      Rcpp::Rcout << "Processing the body (variants).\n";
    }
    // There is at least one variant.
    fi = gzopen( filename.c_str(), "ab" );
    if (! fi) {
      Rcpp::Rcerr << "gzopen of " << filename << " failed: " << strerror (errno) << ".\n";
    }

    for(i = 0; i < fix.nrow(); i++){
      Rcpp::checkUserInterrupt();

      if(mask == 1 && fix(i,6) != "PASS" ){
        // Don't print variant.
      } else {
        // Print variant.
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

        gzwrite(fi, tmpstring.c_str(), tmpstring.size());
        gzwrite(fi,"\n",strlen("\n"));
      }
    }
    if( verbose == 1 ){
      Rcpp::Rcout << "Finished processing the body (variants).\n";
    }
    gzclose(fi);
  } else {
    if( verbose == 1 ){
      Rcpp::Rcout << "No rows (variants).\n";
    }
  }
  
//  return void;
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
  int i = 0;
//  unsigned int i = 0;
  
  if(verbose == 1){
    Rcpp::Rcout << "Processing sample: " << seqname << "\n";
  }

  putc ('>' , pFile);
  for(i=0; (unsigned)i<seqname.size(); i++){
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



