#include <Rcpp.h>
using namespace Rcpp;


// [[Rcpp::export]]
Rcpp::DataFrame window_init(int window_size, int max_bp) {
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


// [[Rcpp::export]]
Rcpp::DataFrame windowize_fasta(Rcpp::DataFrame wins, Rcpp::CharacterVector seq) {
  // The sequence must begin at position one.

  Rcpp::NumericVector ends = wins["end"];
  
  Rcpp::NumericVector A(ends.size());
  Rcpp::NumericVector C(ends.size());
  Rcpp::NumericVector G(ends.size());
  Rcpp::NumericVector T(ends.size());
  Rcpp::NumericVector N(ends.size());
  Rcpp::NumericVector other(ends.size());

  int window_num = 0;
  
  // Vectors are zero based.
  // Sequences are one based.
  for(int i; i < seq.size(); i++){
    if(i + 1 > ends(window_num)){
      window_num++;
    }

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



// Windowize variant number
// [[Rcpp::export]]
Rcpp::DataFrame windowize_variants(Rcpp::DataFrame wins, Rcpp::NumericVector pos) {
  Rcpp::NumericVector ends = wins["end"];
  Rcpp::NumericVector var_counts(ends.size());
  int window_num = 0;
  
  // Vectors are zero based.
  // Positions are one based.
  for(int i; i < pos.size(); i++){
    if(pos(i) > ends(window_num)){
      window_num++;
    } else {
      var_counts(window_num)++;
    }
  }

  return DataFrame::create(wins, _["variants"]=var_counts);
}




// Windowize genic nucleotides
// [[Rcpp::export]]
Rcpp::DataFrame windowize_annotations(Rcpp::DataFrame wins,
                                      Rcpp::NumericVector ann_starts,
                                      Rcpp::NumericVector ann_ends,
                                      int chrom_length) {
  Rcpp::NumericVector win_ends = wins["end"];
//  Rcpp::NumericVector ann_starts = ann(3);
//  Rcpp::NumericVector ann_ends = ann(4);
    
  Rcpp::NumericVector chrom(chrom_length);
  Rcpp::NumericVector window_tally(win_ends);
  int i;
  int j;
  int tmp;

  // We create a chromosome of zeros and replace
  // each genic position with a one.
  // Some positions may be annotated more than
  // once, so we need to handle this too.
  
  // Reorient reverse strand features.
  for(i = 0; i < ann_starts.size(); i++){
    if(ann_starts(i) > ann_ends(i)){
      tmp = ann_starts(i);
      ann_starts(i) = ann_ends(i);
      ann_ends(i) = tmp;
    }
    
    // Mark genic bases.
    for(j = ann_starts(i); j < ann_ends(i); j++){
      Rcout << ann_starts(i) << "\t" << ann_ends(i) << "\n";
      chrom(i) = 1;
    }
  }

  // Now tally the number of genic positions in each window.
  tmp = 0;
  window_tally(tmp) = 0;
  for(i = 0; i < chrom_length; i++){
//    Rcout << chrom(i); // << "\n";
//    if(i > win_ends(tmp)){
//      tmp++;
//      window_tally(tmp) = 0;
//      }
//    if(chrom(i) == 1){window_tally(tmp)++;}
  }

  return DataFrame::create(wins, _["genic"]=window_tally);
//  return wins;
}



