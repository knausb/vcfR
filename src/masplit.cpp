#include <Rcpp.h>
#include "vcfRCommon.h"

// Helper for std::sort.
struct greater
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a > b; }
};

struct lesser
{
    template<class T>
    bool operator()(T const &a, T const &b) const { return a < b; }
};


// Convert vectors of strings to floats.
std::vector<float> str_vec_to_float_vec( std::vector<std::string> str_vec ){
  // Initialize return vector.
  std::vector<float> float_vec( str_vec.size(), 0 );
  
  unsigned int i = 0;
  for( i=0 ; i < str_vec.size() ; i++ ){
  
//    Rcpp::Rcout << "  " << str_vec[i] << "\n";
    std::istringstream ss0(str_vec[i]);
    if( str_vec[i] == "." ){
      float_vec[i] = -99999;
    } else if ( !( ss0 >> float_vec[i] ) ){
      // error: didn't convert to a float
      Rcpp::Rcout << "ss0: " << ss0.str() << "\n";
      Rcpp::Rcerr << "Failed to convert to a float.\n";
    }
  }
  return float_vec;
}


// Convert vectors of strings to NumericVector.
Rcpp::NumericVector str_vec_to_NumericVector( std::vector<std::string> str_vec ){
  // Initialize return vector.
//  std::vector<float> float_vec( str_vec.size(), 0 );
  Rcpp::NumericVector num_vec( str_vec.size(), 0 );
  
  unsigned int i = 0;
  for( i=0 ; i < str_vec.size() ; i++ ){
    
    //    Rcpp::Rcout << "  " << str_vec[i] << "\n";
    std::istringstream ss0(str_vec[i]);
    
    if( str_vec[i] == "." ){
      num_vec[i] = NA_REAL;
    } else if ( !( ss0 >> num_vec[i] ) ){
      // error: didn't convert to a float
      Rcpp::Rcout << "ss0: " << ss0.str() << "\n";
      Rcpp::Rcerr << "Failed to convert to a float.\n";
    }
  }
  return num_vec;
}



//' 
//' @rdname masplit
//' 
//' @title masplit
//' @description Split a matrix of delimited strings.
//' 
//' @param myMat a matrix of delimited strings (e.g., "7,2").
//' @param delim character that delimits values.
//' @param count return the count of delimited records.
//' @param record which (1-based) record to return.
//' @param sort should the records be sorted prior to selecting the element (0,1)?
//' @param decreasing should the values be sorted decreasing (1) or increasing (0)?
//' 
//' 
//' @details 
//' Split a matrix of delimited strings that represent numerics into numerics.
//' The parameter \strong{count} returns a matrix of integers indicating how many delimited records exist in each element.
//' This is intended to help if you do not know how many records are in each element particularly if there is a mixture of numbers of records.
//' The parameter \strong{record} indicates which record to return (first, second, third, ...).
//' The parameter \strong{sort} indicates whether the records in each element should be sorted (1) or not (0) prior to selection.
//' When sorting has been selected \strong{decreasing} indicates if the sorting should be performed in a decreasing (1) or increasing (0) manner prior to selection.
//' 
//' 
//' 
//' 
//' @return A numeric matrix
//' 
//' 
//' @examples
//' set.seed(999)
//' x1 <- round(rnorm(n=9, mean=10, sd=2))
//' x2 <- round(rnorm(n=9, mean=20, sd=2))
//' ad <- matrix(paste(x1, x2, sep=","), nrow=3, ncol=3)
//' colnames(ad) <- paste('Sample', 1:3, sep="_")
//' rownames(ad) <- paste('Variant', 1:3, sep="_")
//' ad[1,1] <- "9,23,12"
//' is.na(ad[3,1]) <- TRUE
//' 
//' ad
//' masplit(ad, count = 1)
//' masplit(ad, sort = 0)
//' masplit(ad, sort = 0, record = 2)
//' masplit(ad, sort = 0, record = 3)
//' masplit(ad, sort = 1, decreasing = 0)
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::NumericMatrix masplit(Rcpp::StringMatrix myMat,
                                 std::string delim = ",",
                                 int count = 0,
                                 int record = 1,
                                 int sort = 1,
                                 int decreasing = 1
                                 ) {

  // Initialize return data structure.
  Rcpp::NumericMatrix retMat( myMat.nrow(), myMat.ncol() );
  retMat.attr("dimnames") = myMat.attr("dimnames");
  
  // NA matrix to return in case of unexpected results.
  Rcpp::NumericMatrix naMat( 1, 1 );
  naMat(0,0) = NA_REAL;
    
  if( record < 1){
    Rcpp::Rcerr << "Specified record number is less than one.\n";
    return naMat;
  }
  
  // R is one based, C++ zero based.
  record = record - 1;

  int i;
  int j;

  for(i=0; i<retMat.nrow(); i++){   // Count rows (variants).
    for(j=0; j<retMat.ncol(); j++){ // Count columns (samples).

      if( myMat(i,j) != NA_STRING ){
        std::vector < std::string > col_vec;
//        char my_split = ','; // Must be single quotes!
        char my_split = delim[0];
        std::string my_string;
        my_string = myMat(i,j);
        vcfRCommon::strsplit(my_string, col_vec, my_split);

        // Recast vector of string to vector of floats.
        std::vector < float > float_vec;( col_vec.size(), 0);
        float_vec = str_vec_to_float_vec(col_vec);
//        Rcpp::NumericVector col_vec2( col_vec.size(), 0);
//        col_vec2 = str_vec_to_NumericVector(col_vec);
        
        // Process the vector.
        if( count == 1 ){
          // Return the length of the vector.
          retMat(i,j) = float_vec.size();
//          retMat(i,j) = col_vec2.size();
        } else {

          // Sort the vector.
          if( sort == 1 ){
            if( decreasing == 1 ){
              std::sort ( float_vec.begin(), float_vec.end(), greater() );
            } else if ( decreasing == 0 ){
              std::sort ( float_vec.begin(), float_vec.end(), lesser() );
            } else {
              Rcpp::Rcerr << "Specification of 'decreasing' should be either 0 or 1.\n";
              return naMat;
            }
          }  
          
          // Select the record.
          if( (unsigned)record + 1 > float_vec.size() ){
            retMat(i,j) = NA_REAL;
          } else if( float_vec[ record ] == -99999 ){
            retMat(i,j) = NA_REAL;
          } else {
            retMat(i,j) = float_vec[ record ];
          }
        }
      } else if( myMat(i,j) == NA_STRING ){
        retMat(i,j) = NA_REAL;
//        Rcpp::Rcout << "NA input to NA out.\n";
      }

//      } else {
//        retMat(i,j) = NA_REAL;
//      }
    }
  }

  return retMat;
}



