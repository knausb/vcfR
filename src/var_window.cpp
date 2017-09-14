#include <Rcpp.h>

using namespace Rcpp;

//' @export
// [[Rcpp::export(name=".window_init")]]
Rcpp::DataFrame window_init(int window_size, int max_bp) {
  // Initialize windows.
  // Input is the window size and the maximum coordinate
  // to include in the windows.
  // A DataFrame is returned which contains a vector of 
  // window numbers, start, end and length of each window.
  
  int max_window = max_bp / window_size;
  max_window++;
  Rcpp::NumericVector window(max_window);
  Rcpp::NumericVector start(max_window);
  Rcpp::NumericVector end(max_window);
  Rcpp::NumericVector length(max_window);

  for(int i=0; i<window.size(); i++){
    window(i) = i+1;
    start(i) = i * window_size + 1;
    end(i) = (i + 1) * window_size;
    length(i) = window_size;
  }

  return DataFrame::create(_["window"]=window, _["start"]=start, _["end"]=end, _["length"]=length);
}


//' @export
// [[Rcpp::export(name=".windowize_fasta")]]
Rcpp::DataFrame windowize_fasta(Rcpp::DataFrame wins, Rcpp::CharacterVector seq) {
  // Windowizes the nucleotide information from a DNA sequence.
  // Input is a DataFrame containing a vector of end points for
  // each window.  This could come from vcfR::window_init.
  // The other input is a vector of chatacters from a sequence.
  // Output is a DataFrame that includes the input DataFrame and
  // appends to it the windowed nucleotide counts.
  
  // The sequence must begin at position one.

  Rcpp::NumericVector ends = wins["end"];
  
  Rcpp::NumericVector A(ends.size());
  Rcpp::NumericVector C(ends.size());
  Rcpp::NumericVector G(ends.size());
  Rcpp::NumericVector T(ends.size());
  Rcpp::NumericVector N(ends.size());
  Rcpp::NumericVector other(ends.size());

  int window_num = 0;

//Rcout << "seq.size: " << seq.size() << "\n";
//Rcout << "ends.size: " << ends.size() << "\n";

  // Vectors are zero based.
  // Sequences are one based.
  for(int i=0; i < seq.size(); i++){
    if(i + 1 > ends(window_num)){
      window_num++;
    }
//Rcout << "Made it to: " << i << "\n";

    if(seq(i) == "A")
      A(window_num)++;
    else if(seq(i) == "a")
      A(window_num)++;
    else if(seq(i) == "C")
      C(window_num)++;
    else if(seq(i) == "c")
      C(window_num)++;
    else if(seq(i) == "G")
      G(window_num)++;
    else if(seq(i) == "g")
      G(window_num)++;
    else if(seq(i) == "T")
      T(window_num)++;
    else if(seq(i) == "t")
      T(window_num)++;
    else if(seq(i) == "N")
      N(window_num)++;
    else if(seq(i) == "n")
      N(window_num)++;
    else
      other(window_num)++;
  }

  return DataFrame::create(wins, _["A"]=A, _["C"]=C, _["G"]=G, _["T"]=T, _["N"]=N, _["other"]=other);
}



// Windowize variants
//
//' @export
// [[Rcpp::export(name=".windowize_variants")]]
Rcpp::DataFrame windowize_variants(Rcpp::DataFrame windows, Rcpp::DataFrame variants) {
  Rcpp::NumericVector ends = windows["end"];
  Rcpp::NumericVector pos = variants["POS"];
  Rcpp::LogicalVector mask = variants["mask"];
  Rcpp::NumericVector var_counts(ends.size());
  int window_num = 0;
  
  // Vectors are zero based.
  // Positions are one based.
  for(int i=0; i < pos.size(); i++){
//    Rcout << "i: " << i << ", pos(i): " << pos(i) << ", ends(window_num): " << ends(window_num) <<"\n";
    while(pos(i) > ends(window_num)){
      window_num++;
    }
//    Rcout << "i: " << i << ", pos(i): " << pos(i) << ", ends(window_num): " << ends(window_num) <<"\n";
    if (mask(i) == TRUE){
      var_counts(window_num)++;
    }
  }

  return DataFrame::create(windows, _["variants"]=var_counts);
}




// Windowize annotated nucleotides
//
//' @export
// [[Rcpp::export(name=".windowize_annotations")]]
Rcpp::DataFrame windowize_annotations(Rcpp::DataFrame wins,
                                      Rcpp::NumericVector ann_starts,
                                      Rcpp::NumericVector ann_ends,
                                      int chrom_length) {
                                        
  // Tally the number of annotated nucleotides in windows.

  Rcpp::NumericVector win_ends = wins["end"];
  Rcpp::NumericVector chrom(chrom_length);
  Rcpp::NumericVector window_tally(win_ends.size());
  int i;
  int j;
  int tmp;

  // We create a chromosome of zeros and replace
  // each genic position with a one.
  // Some positions may be annotated more than
  // once, so we need to handle this too.
  
//  for(i=0; i<chrom_length; i++){ Rcout << chrom(i); }
//  for(i=0; i < window_tally.size(); i++){ Rcout << window_tally(i) << "\n"; }
  
  // Reorient reverse strand features.
  for(i = 0; i < ann_starts.size(); i++){
    if(ann_starts(i) > ann_ends(i)){
      tmp = ann_starts(i);
      ann_starts(i) = ann_ends(i);
      ann_ends(i) = tmp;
    }

    // Mark genic bases.
//    Rcpp::Rcout << ann_starts(i) << "\t" << ann_ends(i) << "\n";
    for(j = ann_starts(i); j < ann_ends(i); j++){
      chrom(j) = 1;
    }

  }

  // Now tally the number of genic positions in each window.
  tmp = 0;
  for(i = 0; i < chrom_length; i++){
//    Rcout << chrom(i); // << "\n";
    if(i + 1 > win_ends(tmp)){
      tmp++;
      }
    if(chrom(i) == 1){window_tally(tmp)++;}
  }

  return DataFrame::create(wins, _["genic"]=window_tally);
//  return wins;
}



