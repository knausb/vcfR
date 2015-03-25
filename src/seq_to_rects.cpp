#include <Rcpp.h>
using namespace Rcpp;



int minimal_regex(char query, std::string db){
  int test = 0;
  int i = 0;
  
  while(test == 0){
//    Rcout << "Query: " << query << " db: " << db << "\n";
    if(query == db[i]){
      test = 1;
      } else {
        i++;
      }
  }
//  Rcout << "Done" << "\n\n";
  
  return test;
}



// [[Rcpp::export]]
Rcpp::IntegerMatrix seq_to_rects(Rcpp::CharacterVector seq, std::string targets) {
  std::vector < int > starts;
  std::vector < int > ends;
//  Rcpp::IntegerVector starts;
//  Rcpp::IntegerVector ends;

  int i=0;
  
/*
  for(i=0; i<targets.size(); i++){
    Rcout << targets[i] << " ";
  }
  Rcout << "\n";
*/

  // Test first nucleotide.
//  int test = minimal_regex(Rcpp::as<char>(seq[0]), targets);
  
  int in_rect = 0;
  int test = 0;
  
  for(i=0; i<seq.size(); i++){
//    Rcout << seq[i] << "\n";
    test = minimal_regex(Rcpp::as<char>(seq[i]), targets);
    if(in_rect == 0 && test == 0){
      // Not in rectangle and nucleotide does not belong in one.
    } else if(in_rect == 0 && test == 1){
      // Not in rectangle but need to initiate one.
      starts.push_back(i + 1);
      in_rect = 1;
    } else if(in_rect == 1 && test == 0){
      // In rectangle but need to end it.
      ends.push_back(i);
      in_rect = 0;      
    } else if(in_rect == 1 && test == 1){
      // In rectangle and nucleotide belongs in one.
    }
  }

  // Handle the last nucleotide.
  if(in_rect == 0 && test == 0){
    // Not in rectangle and nucleotide does not belong in one.
  } else if(in_rect == 0 && test == 1){
    // Not in rectangle but need to initiate one.
    // A one bp feature.
//    starts.push_back(i +1);
    ends.push_back(i + 1);
  } else if(in_rect == 1 && test == 0){
    // In rectangle but need to end it.
    // Should have been finalized above.
  } else if(in_rect == 1 && test == 1){
    // In rectangle and nucleotide belongs in one.
    // Need to end rectangle.
    ends.push_back(i + 1);
  }

//  Rcpp::IntegerMatrix rects = Rcpp::IntegerMatrix::create(_["starts"] = starts, _["ends"] = ends);
  Rcpp::IntegerMatrix rects(starts.size(), 2);
  for(i=0; i<starts.size(); i++){
    rects(i,0) = starts[0];
    rects(i,1) = ends[0];
  }
//  rects.names() = Rcpp::StringVector::create("starts", "ends");
  rects.attr("dimnames") = Rcpp::List::create(Rcpp::StringVector::create(), 
                                   Rcpp::StringVector::create("starts", "ends"));

  return rects;
}


