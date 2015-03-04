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



// Extract windows of numeric data from genotype field of VCF
// 
// @param x A NumericMatrix
// @param pos A vector of chromosomal positions
// @param maxbp Length of chromosome
// @param winsize Size (in bp) for windows
// @export
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


