#include <Rcpp.h>

using namespace Rcpp;


NumericVector win_mean(std::vector< std::vector<double> > win){
  unsigned int i = 0;
  unsigned int j = 0;
  int k = 0;
  NumericVector means(win.size());

  // Sum over samples.
  for(i=0; i<win.size(); i++){
    if(win[i].size() > 0){
      means[i] = win[i][0];
      if(win[i].size() > 1){
        for(j=1; j<win[i].size(); j++){
          means[i] = means[i] + win[i][j];
        }
      }
    } else {
      means[i] = 0;
    }
  }

  for(k=0; k < means.size(); k++){
    if(win[k].size() > 1){
      means[k] = means[k] / win[k].size();
    }
  }

  return means;
}






// [[Rcpp::export]]
NumericMatrix NM2winNM(NumericMatrix x, std::vector<int> pos, int maxbp, int winsize=100) {
  int nwins;
  nwins = maxbp % winsize;
  if(maxbp % winsize == 0){
    nwins = maxbp/winsize;
  } else {
    nwins = 1 + maxbp/winsize;
  }
  
//  Rcout << "nwins: " << nwins << "\n";

  // Output matrix.
  NumericMatrix outM(nwins, x.ncol());

  // Window boundaries and counter.
  int winmin = 1;
  winsize = winsize - 1;
  int winmax = + winsize;
  int winnum = 0;
  int i = 0;
  unsigned int j = 0;

  // Hold data for each window and sample.
//  std::vector<NumericVector> temp(x.ncol());
  std::vector< std::vector<double> > temp(x.ncol());


  // Loop over data.
//  for(i=0; i <= outM.nrow(); i++){
  for(i=0; i <= x.nrow(); i++){
//    if(pos[0]  > winmin){
//    } else 

//    Rcout << "Record: " << pos[i] << ",\t" << x(i, 0) << ",\t" << x(i, 1) << "\n";

    if (pos[i] > winmax){
      // Proc window.
      outM(winnum,_) = win_mean(temp);
      
      // Clear window.
      for(j=0; j<temp.size(); j++){
        temp[j].clear();
      }

      // Increment window
      winmin = winmax + 1;
      winmax = winmin + winsize;
      winnum++;
      
    } else {
      // Add to window.
      for(j=0; j<temp.size(); j++){
        if( x(i, j) > 0){
          temp[j].push_back(x(i,j));
        }
      }
    }
  }
  
  // Proc last window.

  return outM;
}


double vector_count(std::vector<double> x){
//  Rcout << "In vector_count\n";
//  for(int i=0; i<x.size(); i++){
//    Rcout << x[i] << "; ";
//    sum = sum + x[i];
//  }
//  Rcout << "Count: " << x.size() << "\n";
//  Rcout << "\nEnd vector_count\n\n";
  return x.size();
}


double vector_sum(std::vector<double> x){
//  Rcout << "In vector_sum\n";
  double sum = 0;
  unsigned int i=0;
  for(i=0; i<x.size(); i++){
//    Rcout << x[i] << "; ";
    sum = sum + x[i];
  }
//  Rcout << "Sum: " << sum << "\n";
//  Rcout << "\nEnd vector_sum\n\n";
  return sum;
}


double vector_mean(std::vector<double> x){
  double mean = 0;
  double count = 0;
  unsigned int i = 0;
  for( i=0; i<x.size(); i++ ){
    if( ( x[i] != NA_REAL ) | ( x[i] != NA_INTEGER ) ){
      mean = mean + x[i];
      count++;
    }
  }
//  mean = mean / x.size();
  mean = mean / count;
  return mean;
}


double vector_median(std::vector<double> x){
  // http://stackoverflow.com/a/2114817
  double median = 0;
  size_t x_size = x.size();
  
  sort(x.begin(), x.end());

  if (x_size  % 2 == 0){
      median = (x[x_size / 2 - 1] + x[x_size / 2]) / 2;
  } else {
      median = x[x_size / 2];
  }

  return median;
}


// [[Rcpp::export]]
NumericMatrix windowize_NM(Rcpp::NumericMatrix x, Rcpp::NumericVector pos,
                           Rcpp::NumericVector starts, Rcpp::NumericVector ends,
                           Rcpp::String summary="mean") {

  // Declare a matrix for output.
  Rcpp::NumericMatrix outM(starts.size(), x.ncol());
  Rcpp::List dimnames = x.attr("dimnames");
  dimnames(0) = Rcpp::CharacterVector::create();
  outM.attr("dimnames") = dimnames;

  // Declare a vector of a vector of doubles to hold each
  // window prior to summarization.
  std::vector < std::vector < double > > window_tmp(x.ncol());

  // Counters
  int window_num = 0;
  int i; // Variant counter
  int j; // Sample counter

  // Manage the first variant.
  while(pos(0) > ends(window_num)){
    window_num++;
  }
  
  // Scroll over variants.
  for(i = 0; i < pos.size(); i++){
    if(pos(i) > ends(window_num)){
      // Summarize window.
//      Rcout << "  Sumarizing window " << window_num << ", size is: " << window_tmp[0].size() << "\n";
      if(summary == "mean"){
        for(j=0; j<x.ncol(); j++){
          outM(window_num, j) = vector_mean(window_tmp[j]);
          window_tmp[j].clear();
          window_tmp[j].push_back(x(i,j));
        }
      }
      if(summary == "median"){
        for(j=0; j<x.ncol(); j++){
          outM(window_num, j) = vector_median(window_tmp[j]);
          window_tmp[j].clear();
          window_tmp[j].push_back(x(i,j));
        }
      }
      if(summary == "sum"){
        for(j=0; j<x.ncol(); j++){
          outM(window_num, j) = vector_sum(window_tmp[j]);
          window_tmp[j].clear();
          window_tmp[j].push_back(x(i,j));
        }
      }
      if(summary == "count"){
        for(j=0; j<x.ncol(); j++){
          outM(window_num, j) = vector_count(window_tmp[j]);
          window_tmp[j].clear();
          window_tmp[j].push_back(x(i,j));
        }
      }
      window_num++;
    } else {
      // Add values to current window.
      for(j=0; j<x.ncol(); j++){
        if(x(i,j) != NA_REAL){
          window_tmp[j].push_back(x(i,j));
        }
      }
    }
    
//    Rcout << "i: " << i << ", pos: " << pos(i) << ", window_num:" << window_num << "\n";
  }
  
  // Summarize the last window.
  if(summary == "mean"){
    for(j=0; j<x.ncol(); j++){
      outM(window_num, j) = vector_mean(window_tmp[j]);
      window_tmp[j].clear();
    }
  }
  if(summary == "median"){
    for(j=0; j<x.ncol(); j++){
      outM(window_num, j) = vector_median(window_tmp[j]);
      window_tmp[j].clear();
    }
  }
  if(summary == "sum"){
    for(j=0; j<x.ncol(); j++){
      outM(window_num, j) = vector_sum(window_tmp[j]);
      window_tmp[j].clear();
    }
  }      
  if(summary == "count"){
    for(j=0; j<x.ncol(); j++){
      outM(window_num, j) = vector_count(window_tmp[j]);
      window_tmp[j].clear();
    }
  }

  return outM;
}


