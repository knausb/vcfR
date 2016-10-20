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
    for(j=start_row; j<=end_row; j++){
      retMatrix(j - start_row + 0, i) = myMat(j,i);
    }
  }
  
  return(retMatrix);
}


// Count non missing values in a NumericMatrix.
Rcpp::NumericVector count_nonNA( Rcpp::NumericMatrix myMat ){
  Rcpp::NumericVector myCounts( myMat.ncol() );

  int i = 0;
  int j = 0;
  
  // Initialize to zero.
  for(i=0; i<myCounts.size(); i++){
    myCounts(i) = 0;
  }

  for(i=0; i<myMat.ncol(); i++){
    for(j=0; j<myMat.nrow() - 1; j++){
      if( !Rcpp::NumericVector::is_na( myMat(j,i) ) ){
        myCounts(i) = myCounts(i) + 1;
      }
    }
  }

  return(myCounts);  
}

// Sort data into bins based on breaks (boundaries).
Rcpp::NumericMatrix bin_data( Rcpp::NumericVector myFreqs,
                              float bin_width
){
  int i = 0;
  int j = 0;
  int nbins = 1/bin_width;
  Rcpp::NumericMatrix breaks( nbins, 4 );
//  int multiplier = 1000; // Convert floating points to ints.
  int multiplier = 10000000; // Convert floating points to ints, graphics::hist uses 1e7.
  Rcpp::IntegerMatrix intBreaks( nbins, 4 );
  
  Rcpp::StringVector colnames(4);
  colnames(0) = "START";
  colnames(1) = "MID";
  colnames(2) = "END";
  colnames(3) = "COUNT";
  Rcpp::colnames(breaks) = colnames;
  Rcpp::colnames(intBreaks) = colnames;
  
//  Rcpp::IntegerVector counts(nbins);
  
  breaks(0,0) = 0;
  breaks(0,1) = bin_width/2;
  breaks(0,2) = bin_width;
  intBreaks(0,0) = 0;
  intBreaks(0,1) = (bin_width/2) * multiplier;
  intBreaks(0,2) = bin_width * multiplier;
  
  for(i=1; i<breaks.nrow(); i++){
    breaks(i,0) = breaks(i-1,0) + bin_width;
    breaks(i,1) = breaks(i-1,1) + bin_width;
    breaks(i,2) = breaks(i-1,2) + bin_width;
    
    intBreaks(i,0) = intBreaks(i-1,0) + bin_width * multiplier;
    intBreaks(i,1) = intBreaks(i-1,1) + bin_width * multiplier;
    intBreaks(i,2) = intBreaks(i-1,2) + bin_width * multiplier;
  }
  
//  Rcpp::Rcout << "\nBinning!\n\n";
  
  for(i=0; i<myFreqs.size(); i++){
    int intQuery = myFreqs(i) * multiplier;
    j = 0;
    if( intQuery >= intBreaks(j,0) & intQuery <= intBreaks(j,2) ){
//      Rcpp::Rcout << "Binned: " <<  myFreqs(i) << " is >= " << breaks(j,0) << " & <= " << breaks(j,2) << "\n";
//      breaks(j,3) = breaks(j,3) + 1;
    }
    for(j=1; j<breaks.nrow(); j++){
      if( intQuery > intBreaks(j,0) & intQuery <= intBreaks(j,2) ){
//        Rcpp::Rcout << "Binned: " <<  myFreqs(i) << " is > " << breaks(j,0) << " & <= " << breaks(j,2) << "\n";
        breaks(j,3) = breaks(j,3) + 1;
      }
    }
    
  }
  
  return(breaks);  
}


double find_one_peak( Rcpp::NumericMatrix binned_data, 
                      Rcpp::LogicalVector lhs
                    ){
  double myPeak = 0;
  int max_peak = 0;
  int i = 0;

  for(i=1; i<binned_data.nrow(); i++){
    if( lhs(0) == 1 ){
      if( binned_data(i,3) > binned_data(max_peak,3) ){
        max_peak = i;
      }
    } else {
      if( binned_data(i,3) >= binned_data(max_peak,3) ){
        max_peak = i;
      }
    }
  }
    
  myPeak = binned_data(max_peak,1);
  return( myPeak );  
}



void dput_NumericMatrix( Rcpp::NumericMatrix myMat){
  int i = 0;
  int j = 0;
  
  Rcpp::StringVector myRowNames( myMat.nrow() );
  Rcpp::StringVector myColNames( myMat.ncol() );
  
  if( !Rf_isNull(rownames(myMat)) && Rf_length(rownames(myMat)) > 1 ){
    myRowNames = Rcpp::rownames(myMat);    
  }
//  Rcpp::Rcout << "Checking for colnames!\n";
  if( !Rf_isNull(colnames(myMat)) && Rf_length(colnames(myMat)) > 1 ){
//    Rcpp::Rcout << "Found colnames!\n";
    myColNames = Rcpp::colnames(myMat);
  }
  
  
  Rcpp::Rcout << "\n";
  Rcpp::Rcout << "structure(c(";
  
  // Print the first column.
  if( Rcpp::NumericVector::is_na( myMat(0,0) ) ){
    Rcpp::Rcout << "NA";
  } else {
    Rcpp::Rcout << myMat(0,0);
  }
  for(i=1; i<myMat.nrow(); i++){
//    Rcpp::Rcout << ", " << myMat(i,j);
    if( Rcpp::NumericVector::is_na( myMat(i,j) ) ){
      Rcpp::Rcout << ",NA";
    } else {
      Rcpp::Rcout << "," << myMat(i,j);
    }
  }
  
  // Print the remaining columns.
  for(j=1; j<myMat.ncol(); j++){
    for(i=0; i<myMat.nrow(); i++){
      if( Rcpp::NumericVector::is_na( myMat(i,j) ) ){
        Rcpp::Rcout << ",NA";
      } else {
        Rcpp::Rcout << "," << myMat(i,j);
      }
    }
  }
  
  Rcpp::Rcout << "),";
  
  // Dimensions
  Rcpp::Rcout << " .Dim = c(" << myMat.nrow() << "L, " << myMat.ncol() << "L)";
  
  Rcpp::Rcout << ", .Dimnames = list(";
  // Row names.
  if( !Rf_isNull(rownames(myMat)) && Rf_length(rownames(myMat)) > 1 ){
    //  myRowNames = Rcpp::rownames(myMat);
    Rcpp::Rcout << "c(\"" << myRowNames(0);
    for(i=1; i<myRowNames.size(); i++){
      Rcpp::Rcout << "\", \"" << myRowNames(i);
    }
    Rcpp::Rcout << "\")";
  } else {
    Rcpp::Rcout << "NULL";
  }
  
  Rcpp::Rcout << ",";
  // Column names.
  if( !Rf_isNull(colnames(myMat)) && Rf_length(colnames(myMat)) > 1 ){
    Rcpp::Rcout << "c(\"" << myColNames(0);
    for(i=1; i<myColNames.size(); i++){
      Rcpp::Rcout << "\", \"" << myColNames(i);
    }
    Rcpp::Rcout << "\")";
  } else {
    Rcpp::Rcout << "NULL";
  }
  
  //
  Rcpp::Rcout << ")"; // Close .Dimnames list.
  Rcpp::Rcout << ")\n"; // Close structure.
  Rcpp::Rcout << "\n\n"; // Delimit dput statements.
}



// Find peaks from frequency values [0-1]
// from a single window (matrix of columns) of data.
//
Rcpp::NumericVector find_peaks( Rcpp::NumericMatrix myMat, 
                                float bin_width,
                                Rcpp::LogicalVector lhs
                                ){
  // myMat is a matrix that consists of samples in rows
  // and one window's length of frequencies in rows.
  
  int i = 0;
//  int j = 0;
//  int k = 0;
  
  // Create return vector and initialize to zero.
  // Return vector contains peaks and is as long
  // as the number of samples which is the same as
  // the number of columns in myMat.
  Rcpp::NumericVector myPeaks( myMat.ncol() );
  for(i=0; i<myPeaks.size(); i++){
    myPeaks(i) = 0;
  }

  int nbins = 1 / bin_width;
  std::vector<double> breaks ( nbins + 1, 0 );
  Rcpp::NumericVector mids( nbins );

  // Test thet 1/bin_width does not have a remainder.
  if( 1/bin_width - nbins > 0 ){
    Rcpp::Rcerr << "1/bin_width has a remainder.\nThis will result in uneven bins.\nPlease try another bin_width.\n";
    return(myPeaks);
  }
  
  // Initialize vectors.
  breaks[0] = 0;
  for(i=0; i<mids.size(); i++ ){
    breaks[i+1] = breaks[i] + bin_width;
    mids(i) = breaks[i] + bin_width/2;
  }

  // Process by sample (column) each set of frequencies in the matrix myMat.
  // First bin the data using 'bin_data'.
  // Then find the peak in the data using 'find_one_peak'.
  //
  // Debugging info:
  // dput_NumericMatrix will take an Rcpp::NumericMatrix and print it to stdout in
  // a format similar to R's dput.
  //
  for(i=0; i<myMat.ncol(); i++){ // Column (sample) counter.
    Rcpp::NumericMatrix binned_data;
    binned_data = bin_data( myMat( Rcpp::_, i), bin_width );

// dput_NumericMatrix( myMat );
// dput_NumericMatrix(binned_data);
//    myPeaks(i) = find_one_peak( myMat( Rcpp::_, i), breaks, mids, lhs );
    myPeaks(i) = find_one_peak( binned_data, lhs );
  }

  return(myPeaks);
}


//' 
//' @rdname freq_peak
//' 
//' @title freq_peak
//' @description Find density peaks in frequency data.
//' 
//' @param myMat a matrix of frequencies [0-1].
//' @param pos a numeric vector describing the position of variants in myMat.
//' @param winsize sliding window size.
//' @param bin_width Width of bins to summarize ferequencies in (0-1].
//' @param lhs logical specifying whether the search for the bin of greatest density should favor values from the left hand side.
//' 
//' @details
//' Noisy data, such as genomic data, lack a clear consensus.
//' Summaries may be made in an attempt to 'clean it up.'
//' Common summaries, such as the mean, rely on an assumption of normalicy.
//' An assumption that frequently can be violated.
//' This leaves a conundrum as to how to effectively summarize these data.
//' 
//' 
//' Here we implement an attempt to summarize noisy data through binning the data and selecting the bin containing the greatest density of data.
//' The data are first divided into parameter sized windows.
//' Next the data are categorized by parameterizable bin widths.
//' Finally, the bin with the greatest density, the greatest count of data, is used as a summary.
//' Because this method is based on binning the data it does not rely on a distributional assumption.
//' 
//' 
//' The parameter `lhs` specifyies whether the search for the bin of greatest density should be performed from the left hand side.
//' The default value of TRUE starts at the left hand side, or zero, and selects a new bin as having the greatest density only if a new bin has a greater density.
//' If the new bin has an equal density then no update is made.
//' This causees the analysis to select lower frequencies.
//' When this parameter is set to FALSE ties result in an update of the bin of greatest density.
//' This causes the analysis to select higher frequencies.
//' It is recommended that when testing the most abundant allele (typically [0.5-1]) to use the default of TURE so that a low value is preferred.
//' Similarly, when testing the less abundant alleles it is recommended to set this value at FALSE to preferentially select high values.
//' 
//' 
//' @return 
//' A list containing:
//' \itemize{
//'   \item a matrix containing window coordinates
//'   \item a matrix containing peak locations
//'   \item a matrix containing the counts of variants for each sample in each window
//' }
//' 
//' The window matrix contains start and end coordinates for each window, the rows of the original matrix that demarcate each window and the position of the variants that begin and end each window.
//' 
//' The matrix of peak locations contains the midpoint for the bin of greatest density for each sample and each window.
//' Alternatively, if `count = TRUE` the number of non-missing values in each window is reported.
//' The number of non-mising values in each window may be used to censor windows containing low quantities of data.
//' 
//' 
//' @examples
//' data(vcfR_example)
//' gt <- extract.gt(vcf)
//' hets <- is_het(gt)
//' # Censor non-heterozygous positions.
//' is.na(vcf@gt[,-1][!hets]) <- TRUE
//' # Extract allele depths.
//' ad <- extract.gt(vcf, element = "AD")
//' ad1 <- masplit(ad, record = 1)
//' ad2 <- masplit(ad, record = 2)
//' freq1 <- ad1/(ad1+ad2)
//' freq2 <- ad2/(ad1+ad2)
//' myPeaks1 <- freq_peak(freq1, getPOS(vcf))
//' is.na(myPeaks1$peaks[myPeaks1$counts < 20]) <- TRUE
//' myPeaks2 <- freq_peak(freq2, getPOS(vcf), lhs = FALSE)
//' is.na(myPeaks2$peaks[myPeaks2$counts < 20]) <- TRUE
//' #myPeaks <- freq_peak(freqs[1:115,], getPOS(vcf)[1:115])
//' 
//' # Visualize
//' mySample <- "P17777us22"
//' myWin <- 2
//' hist(freq1[myPeaks1$wins[myWin,'START_row']:myPeaks1$wins[myWin,'END_row'], mySample], 
//'      breaks=seq(0,1,by=0.02), col="#A6CEE3", main="", xlab="", xaxt="n")
//' hist(freq2[myPeaks2$wins[myWin,'START_row']:myPeaks2$wins[myWin,'END_row'], mySample], 
//'      breaks=seq(0,1,by=0.02), col="#1F78B4", main="", xlab="", xaxt="n", add = TRUE)
//' axis(side=1, at=c(0,0.25,0.333,0.5,0.666,0.75,1), 
//'      labels=c(0,'1/4','1/3','1/2','2/3','3/4',1), las=3)
//' abline(v=myPeaks1$peaks[myWin,mySample], col=2, lwd=2)
//' abline(v=myPeaks2$peaks[myWin,mySample], col=2, lwd=2)
//' 
//' # Visualize #2
//' mySample <- "P17777us22"
//' plot(getPOS(vcf), freq1[,mySample], ylim=c(0,1), type="n", yaxt='n', 
//'      main = mySample, xlab = "POS", ylab = "Allele balance")
//' axis(side=2, at=c(0,0.25,0.333,0.5,0.666,0.75,1), 
//'      labels=c(0,'1/4','1/3','1/2','2/3','3/4',1), las=1)
//' abline(h=c(0.25,0.333,0.5,0.666,0.75), col=8)
//' points(getPOS(vcf), freq1[,mySample], pch = 20, col= "#A6CEE3")
//' points(getPOS(vcf), freq2[,mySample], pch = 20, col= "#1F78B4")
//' segments(x0=myPeaks1$wins[,'START_pos'], y0=myPeaks1$peaks[,mySample],
//'          x1=myPeaks1$wins[,'END_pos'], lwd=3)
//' segments(x0=myPeaks1$wins[,'START_pos'], y0=myPeaks2$peaks[,mySample],
//'          x1=myPeaks1$wins[,'END_pos'], lwd=3)
//' 
//' 
//' 
//' @export
// [[Rcpp::export]]
Rcpp::List freq_peak(Rcpp::NumericMatrix myMat,
                     Rcpp::NumericVector pos,
                     int winsize = 10000,
                     float bin_width = 0.02,
                     Rcpp::LogicalVector lhs = true
                     ){
//  Rcpp::Rcout << "In freq_peak" << "\n"; 
  int i = 0;
//  int j = 0;
  
  // NA matrix to return in case of unexpected results.
  Rcpp::NumericMatrix naMat( 1, 1 );
  naMat(0,0) = NA_REAL;
  
  //                                 //
  // Initialize a matrix of windows. //
  //                                 //

//  Rcpp::Rcout << "Initializing win matrix: " << "\n";
  //Rcpp::Rcout << "pos.size() is: " << pos.size() << ".\n";
  // May have zero variants.
  int max_pos = 0;
  if( pos.size() > 0 ){
     max_pos = pos[ pos.size() - 1 ] / winsize + 1;
  }
  //
//Rcpp::Rcout << "max_pos is: " << max_pos << ".\n"; 
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
  
  
  //                             //
  // Initialize a freq matrix.   //
  //                             //

  Rcpp::NumericMatrix freqs( max_pos, myMat.ncol() );
  Rcpp::NumericMatrix cnts( max_pos, myMat.ncol() );
//  Rcpp::Rcout << "Trying dimnames.\n";
  Rcpp::StringVector myColNames = Rcpp::colnames(myMat);
//  Rcpp::Rcout << "myColNames.size(): " << myColNames.size() << "\n";

  Rcpp::rownames(freqs) = rownames;
  Rcpp::rownames(cnts) = rownames;
  if( myColNames.size() > 0 ){
    Rcpp::colnames(freqs) = myColNames;
    Rcpp::colnames(cnts) = myColNames;
  }
//  Rcpp::Rcout << "Finished dimnames.\n";

//  Rcpp::NumericMatrix cnts = freqs;

  //                                 //
  // Find windows in pos.            //
  // Assign rows and POS to windows. //
  //                                 //
    
  int win_num = 0;
  i = 0;

  if( pos.size() > 0 ){
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
      R_CheckUserInterrupt();
    
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

  
    //                //
    // Sanity checks. //
    //                //
  
    // Check bin_width validity.
//    if( !count(0) ){
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
//    }

    //                    //
    // Process by window. //
    //                    //
  
    // Window counter.
    for(i=0; i<freqs.nrow(); i++){
      R_CheckUserInterrupt();
      Rcpp::NumericMatrix myWin(wins(i,3) - wins(i,2) + 1, freqs.ncol());
      // Remember, R=1-based, C++=0-based!
      myWin = mat_to_win(myMat, wins(i,2) - 1, wins(i,3) - 1 );

//      if( count(0) ){
//          Rcpp::Rcout << "count(0):" << count(0) << " must be true!\n";
      cnts(i,Rcpp::_) = count_nonNA( myWin );
//      } else {
      freqs(i,Rcpp::_) = find_peaks( myWin, bin_width, lhs );
//      }
    }

    
  }
  
  
  // Create the return List.
  Rcpp::List myList = Rcpp::List::create(
    Rcpp::Named("wins") = wins,
    Rcpp::Named("peaks") = freqs,
    Rcpp::Named("counts") = cnts
  );
  
  return(myList);
}



