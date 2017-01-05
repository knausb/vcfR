#include <Rcpp.h>
using namespace Rcpp;



int minimal_regex(char query, std::string db){
  int test = 0;
  unsigned int i = 0;
  
//  Rcout << "Query: " << query << " db: " << db << "\n";
  for(i=0; i<db.size(); i++){
    if(query == db[i]){
      test = 1;
      break;
    }
  }
//  Rcout << "Done, test is: " << test << "\n\n";
  
  return test;
}



// [[Rcpp::export]]
Rcpp::IntegerMatrix seq_to_rects(Rcpp::CharacterVector seq, std::string targets) {
  std::vector < int > starts;
  std::vector < int > ends;
//  Rcpp::IntegerVector starts;
//  Rcpp::IntegerVector ends;

  int i=0;
//  unsigned int i=0;
  int in_rect = 0;
  int test = 0;
  
  for(i = 0; i < seq.size(); i++){
//    Rcout << seq[i] << "\n";
    test = minimal_regex(Rcpp::as<char>(seq[i]), targets);
//    Rcout << i << ": in_rect: " << in_rect << ", test: " << test << "\n";
    
    if(in_rect == 0 && test == 0){
      // Not in rectangle and nucleotide does not belong in one.
    } else if(in_rect == 0 && test == 1){
      // Not in rectangle but need to initiate one.
//      Rcout << "New rectangle.\n";
      starts.push_back(i + 1);
      in_rect = 1;
    } else if(in_rect == 1 && test == 0){
      // In rectangle but need to end it.
//      Rcout << "End of rectangle.\n";
      ends.push_back(i);
      in_rect = 0;      
    } else if(in_rect == 1 && test == 1){
      // In rectangle and nucleotide belongs in one.
    }
  }
//  Rcout << "Ending i: " << i << "\n";

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
    ends.push_back(i);
  }

  Rcpp::IntegerMatrix rects(starts.size(), 2);
  for(i=0; (unsigned)i<starts.size(); i++){
//    Rcout << "starts: " << starts[i] << ", ends: " << ends[i] << "\n";
    rects(i,0) = starts[i];
    rects(i,1) = ends[i];
  }

  rects.attr("dimnames") = Rcpp::List::create(Rcpp::StringVector::create(), 
                                   Rcpp::StringVector::create("starts", "ends"));

  return rects;
}


