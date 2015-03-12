#include <Rcpp.h>
using namespace Rcpp;



NumericVector win_mean(std::vector< std::vector<double> > win){
  int i;
  NumericVector means(win.size());

  // Sum over samples.
  for(i=0; i<win.size(); i++){
    if(win[i].size() > 0){
      means[i] = win[i][0];
      if(win[i].size() > 1){
        for(int j=1; j<win[i].size(); j++){
          means[i] = means[i] + win[i][j];
        }
      }
    } else {
      means[i] = 0;
    }
  }

  for(i=0; i < means.size(); i++){
    if(win[i].size() > 1){
      means[i] = means[i] / win[i].size();
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
  int j = 0;

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


double vector_sum(std::vector<double> x){
  Rcout << "In vector_sum\n";
  double sum = 0;
  for(int i=0; i<x.size(); i++){
    Rcout << x[i] << "; ";
    sum = sum + x[i];
  }
  Rcout << "End vector_sum\n\n";
  return sum;
}


double vector_mean(std::vector<double> x){
  double mean = 0;
  int i = 0;
  for(i=0; i<x.size(); i++){
    mean = mean + x[i];
  }
  mean = mean / x.size();
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
NumericMatrix windowize_NM(Rcpp::NumericMatrix x,
                           Rcpp::NumericVector pos,
                           Rcpp::NumericVector starts,
                           Rcpp::NumericVector ends,
                           Rcpp::String centrality="mean") {
                             

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



//  std::string fname = Rcpp::as<std::string>(centrality);

  Rcout << "centrality set to: ";
//  Rcout << Rcpp::as<std::string>(centrality);
  Rcout << "\n";


//Rcout << "Made it.\n";
//  for(i=0; i<cnames.size(); i++){Rcout << cnames(i) << "\n";}

  // Manage the first variant.
  while(pos(0) > ends(window_num)){
    window_num++;
  }
  
  // Scroll over variants.
  for(i = 0; i < pos.size(); i++){
    if(pos(i) > ends(window_num)){
      // Summarize window.
      if(centrality == "mean"){
        for(j=0; j<x.ncol(); j++){
          outM(window_num, j) = vector_mean(window_tmp[j]);
          window_tmp[j].clear();
        }
      }
      if(centrality == "median"){
        for(j=0; j<x.ncol(); j++){
          outM(window_num, j) = vector_median(window_tmp[j]);
          window_tmp[j].clear();
        }
      }
      if(centrality == "sum"){
        for(j=0; j<x.ncol(); j++){
          outM(window_num, j) = vector_sum(window_tmp[j]);
          window_tmp[j].clear();
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
  }
  // Summarize the last window.
  if(centrality == "mean"){
    for(j=0; j<x.ncol(); j++){
      outM(window_num, j) = vector_mean(window_tmp[j]);
      window_tmp[j].clear();
    }
  }
  if(centrality == "median"){
    for(j=0; j<x.ncol(); j++){
      outM(window_num, j) = vector_median(window_tmp[j]);
      window_tmp[j].clear();
    }
  }
  if(centrality == "sum"){
    for(j=0; j<x.ncol(); j++){
      outM(window_num, j) = vector_sum(window_tmp[j]);
      window_tmp[j].clear();
    }
  }      

  return outM;
}


