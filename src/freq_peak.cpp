#include <Rcpp.h>




// Slice a window of rows out of a matrix.
Rcpp::NumericMatrix mat_to_win( Rcpp::NumericMatrix myMat, 
                                int start_row, 
                                int end_row
){
  Rcpp::NumericMatrix retMatrix( end_row - start_row + 1, myMat.ncol() );
  
  int i = 0;
  int j = 0;
  
  for(i=0; i<myMat.ncol(); i++){
    for(j=start_row; j<end_row; j++){
      retMatrix(j - start_row + 1, i) = myMat(j,i);
    }
  }
  
  return(retMatrix);
}

// Count non missing values in a matrix.
Rcpp::NumericVector count_nonNA( Rcpp::NumericMatrix myMat ){
  Rcpp::NumericVector myCounts( myMat.ncol() );

  int i = 0;
  int j = 0;
  
  // Initialize to zero.
  for(i=0; i<myCounts.size(); i++){
    myCounts(i) = 0;
  }

  for(i=0; i<myMat.ncol(); i++){
    for(j=0; j<myMat.nrow(); j++){
      if( myMat(i,j) != NA_REAL){
        myCounts(i) = myCounts(i) + 1;
      }
    }
  }
  
  return(myCounts);  
}



// Find peaks from frequency values [0-1].
Rcpp::NumericVector find_peaks( Rcpp::NumericMatrix myMat, float bin_width ){
  Rcpp::NumericVector myPeaks( myMat.ncol() );
  int i = 0;
  int j = 0;
  int k = 0;
  
  // Initialize to zero.
  for(i=0; i<myPeaks.size(); i++){
    myPeaks(i) = 0;
  }

  int nbins = 1 / bin_width;
  Rcpp::NumericVector breaks( nbins + 1 );
  Rcpp::NumericVector mids( nbins );
  Rcpp::NumericVector counts( nbins );
  
  // Test thet 1/bin_width does not have a remainder.
  if( 1/bin_width - nbins > 0 ){
    Rcpp::Rcerr << "1/bin_width has a remainder, please try another bin_width.\n";
    return(myPeaks);
  }
  
  // Initialize vectors.
  breaks(0) = 0;
  for(i=1; i<mids.size(); i++ ){
    breaks(i) = breaks(i-1) + bin_width;
    mids(i-1) = breaks(i) - bin_width/2;
    counts(i-1) = 0;
  }

  for(i=0; i<myMat.ncol(); i++){ // Column (sample) counter.
    for(j=0; j<myMat.nrow(); j++){ // Row (variant) counter.
      for(k=0; k<counts.size(); k++){ // Bin counter.
        if( myMat(j,i) >= breaks(i) & myMat(j,i) < breaks(i + 1) ){
          counts(i)++;
        }
      }
      
    }
  }

  return(myPeaks);
}


//' 
//' @rdname freq_peak
//' 
//' @title freq_peak
//' @description Find peaks in frequency data.
//' 
//' @param myMat a matrix of frequencies [0-1].
//' @param pos a numeric vector describing the position of variants in myMat.
//' @param winsize sliding window size.
//' @param bin_width Width of bins to summarize ferequencies in (0-1].
//' @param count logical specifying to count the number of non-NA values intead of reporting peak.
//' 
//' @details
//' More to come.
//' 
//' @return 
//' A list
//' 
//' @examples
//' freqs <- matrix(runif(n=9), ncol=3, nrow=3)
//' pos <- 1:3
//' myPeaks <- freq_peak(freqs, pos)
//' 
//' data(vcfR_example)
//' ad <- extract.gt(vcf, element = "AD")
//' ad1 <- masplit(ad, record = 1)
//' ad2 <- masplit(ad, record = 2)
//' freqs <- ad1/(ad1+ad2)
//' myPeaks <- freq_peak(freqs, getPOS(vcf))
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List freq_peak(Rcpp::NumericMatrix myMat,
                     Rcpp::NumericVector pos,
                     int winsize = 10000,
                     float bin_width = 0.02,
                     Rcpp::LogicalVector count = false
                     ){
  int i = 0;
  int j = 0;
  
  // NA matrix to return in case of unexpected results.
  Rcpp::NumericMatrix naMat( 1, 1 );
  naMat(0,0) = NA_REAL;
  
  // Create a matrix of windows.
// Rcpp::Rcout << "pos.size() is: " << pos.size() << ".\n"; 
  int max_pos = pos[ pos.size() - 1 ] / winsize + 1;
// Rcpp::Rcout << "max_pos is: " << max_pos << ".\n"; 
  Rcpp::NumericMatrix wins( max_pos, 6);
  Rcpp::StringVector rownames( max_pos );
  for(i=0; i<max_pos; i++){
    wins(i,0) = i * winsize + 1;
    wins(i,1) = i * winsize + winsize;
    rownames(i) = "win" + std::to_string(i+1);
  }
//  Rcpp::Rcout << "wins initialized!\n";
  Rcpp::StringVector colnames(6);
  colnames(0) = "START";
  colnames(1) = "END";
  colnames(2) = "START_row";
  colnames(3) = "END_row";
  colnames(4) = "START_pos";
  colnames(5) = "END_pos";
  wins.attr("dimnames") = Rcpp::List::create(rownames, colnames);

  // Initialize a freq matrix.
  Rcpp::NumericMatrix freqs( max_pos, myMat.ncol() );
//  Rcpp::Rcout << "Trying dimnames.\n";
  Rcpp::StringVector myColNames = Rcpp::colnames(myMat);
//  Rcpp::Rcout << "myColNames.size(): " << myColNames.size() << "\n";

  Rcpp::rownames(freqs) = rownames;
  if( myColNames.size() > 0 ){
    Rcpp::colnames(freqs) = myColNames;
  }
  
//  Rcpp::Rcout << "Finished dimnames.\n";
  
  // Find windows in pos.
  int win_num = 0;
  i = 0;

  // First row.
//  Rcpp::Rcout << "First row.\n";
  while( pos(i) < wins(win_num,0) ){
    win_num++;
  }
  wins(win_num,2) = i + 1;
  wins(win_num,4) = pos(0);
  
  // Remaining rows.
//  Rcpp::Rcout << "Windowing.\n";
  for(i=1; i<myMat.nrow(); i++){
    
    if( pos(i) > wins(win_num,1) ){
      // Increment window.
//      Rcpp::Rcout << "  New window, pos(i): " << pos(i) << " wins(win_num,0): " << wins(win_num,0) << " wins(win_num,1): " << wins(win_num,1) << "\n";
      wins(win_num,3) = i;
      wins(win_num,5) = pos(i-1);
      
      while( pos(i) > wins(win_num,1) ){
//        win_num++;
        win_num = win_num + 1;
//        Rcpp::Rcout << "    Incrementing win_num: " << win_num << "\n";
      }
//      Rcpp::Rcout << "    win_num: " << win_num << "\n";
      wins(win_num,2) = i + 1;
      wins(win_num,4) = pos(i);
    }
  }
  
  // Last row.
  wins(win_num,3) = i;
  wins(win_num,5) = pos(i-1);

  // Windowize and process.

  // Check bin_width validity.
  if( !count(0) ){
    // Positive bin width.
    if( bin_width <= 0 ){
      Rcpp::Rcerr << "bin_width must be greater than zero, please try another bin_width.\n";
      Rcpp::List myList = Rcpp::List::create(
        Rcpp::Named("wins") = wins,
        Rcpp::Named("peaks") = naMat
      );
      return( myList );
    }
    
    // bin width less than one.
    if( bin_width > 1 ){
      Rcpp::Rcerr << "bin_width must be no greater than one, please try another bin_width.\n";
      Rcpp::List myList = Rcpp::List::create(
        Rcpp::Named("wins") = wins,
        Rcpp::Named("peaks") = naMat
      );
      return( myList );
    }
    
    // No remainder to bin width.
    int nbins = 1 / bin_width;
    if( 1/bin_width - nbins > 0 ){
      Rcpp::Rcerr << "1/bin_width has a remainder, please try another bin_width.\n";
      Rcpp::List myList = Rcpp::List::create(
        Rcpp::Named("wins") = wins,
        Rcpp::Named("peaks") = naMat
      );
    return( myList );
    }
  }
  
  // Window counter.
  for(i=0; i<freqs.nrow(); i++){
    Rcpp::NumericMatrix myWin(wins(i,3) - wins(i,2) + 1, freqs.ncol());
    // Remember, R=1-based, C++=0-based!
    myWin = mat_to_win(myMat, wins(i,2) - 1, wins(i,3) - 1 );

    
//    Rcpp::Rcout << "count(0):" << count(0) << "\n";
    if( count(0) ){
//          Rcpp::Rcout << "count(0):" << count(0) << " must be true!\n";
      freqs(i,Rcpp::_) = count_nonNA( myWin );
    } else {
      freqs(i,Rcpp::_) = find_peaks( myWin, bin_width );
    }
  }
  
  // Create the return List.
  Rcpp::List myList = Rcpp::List::create(
    Rcpp::Named("wins") = wins,
    Rcpp::Named("peaks") = freqs
  );
  
  return(myList);
}



