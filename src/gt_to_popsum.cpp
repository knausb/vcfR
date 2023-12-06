#include <Rcpp.h>
#include "vcfRCommon.h"

using namespace Rcpp;


std::vector < int > gtsplit(std::string line){
//Rcpp::IntegerVector gtsplit(std::string line){
//  Rcpp::IntegerVector intv;
  std::vector < int > intv;
//  Rcout << "Genotype: " << line << "\n";

  // Case of a single digit.
  if(line.size() == 1){
    intv.push_back(atoi(line.c_str()));
  }
  
  unsigned int start=0;
  unsigned int i = 0;
  for(i=1; i<line.size(); i++){
    if( line[i] == '/' || line[i] == '|' ){
      std::string temp = line.substr(start, i);
//      Rcout << "  i: " << i << ", temp: " << temp << "\n";
      intv.push_back(atoi(temp.c_str()));
      start = i+1;
      i = i+2;
    }
  }

  // Handle last element.
  std::string temp = line.substr(start, line.size());
//  Rcout << "  temp: " << temp << "\n";
  intv.push_back(atoi(temp.c_str()));
  
  return intv;
}


// ' @export
// [[Rcpp::export(name=".gt_to_popsum")]]
Rcpp::DataFrame gt_to_popsum(Rcpp::DataFrame var_info, Rcpp::CharacterMatrix gt) {
  // Calculate popgen summaries for the sample.
  // var_info should contain columns named 'CHROM', 'POS', 'mask' and possibly others.
  Rcpp::LogicalVector   mask = var_info["mask"];
  Rcpp::IntegerVector   nsample(mask.size());
  Rcpp::StringVector    allele_counts(mask.size());
  Rcpp::NumericVector   Hes(mask.size());
  Rcpp::NumericVector   Nes(mask.size());
  
  int i = 0;
  int j = 0;
//  unsigned int j = 0;
  unsigned int k = 0;
  
  for(i=0; i < gt.nrow(); i++){
    // Iterate over variants (rows)
    if(mask[i] == TRUE){
      // Initialize a vector to count alleles
      // for this variant (row).
      std::vector<int> myalleles (1,0);
      for(j=0; j < gt.ncol(); j++){
        // Iterate over samples (columns)
        if(gt(i, j) != NA_STRING){
          // Process a single genotype.
          nsample[i]++;  // Increment sample count.

//          Rcout << "gt: " << gt(i, j) << "\n";

          // Count alleles per sample.
          int unphased_as_na = 0; // 0 == FALSE
          std::vector < std::string > gt_vector;
          std::string gt2 = as<std::string>(gt(i,j));
          vcfRCommon::gtsplit( gt2, gt_vector, unphased_as_na );
          
//          Rcout << "  gt_vector.size: " << gt_vector.size() << "\n";

          for(k=0; k<gt_vector.size(); k++){
//            Rcout << "  gt_vector[k]: " << gt_vector[k] << ", ";
            if( gt_vector[k] != "." ){
              unsigned int myAllele = std::stoi(gt_vector[k]);
//              Rcout << "myAllele: " << myAllele;
//            // If this genotype had an allele we did not previously observe
            // we'll have to grow the vector.
              while(myalleles.size() - 1 < myAllele){
                myalleles.push_back(0);
              }
            myalleles[myAllele]++;
            }
          }
//          Rcout << "\n";
        }
      }

      // Report allele vector.
//      Rcpp::Rcout << "Allele vector size: " << myalleles.size();
//      Rcpp::Rcout << "; Allele vector: " << myalleles[0];
//      for(k=1; k<myalleles.size(); k++){
//        Rcpp::Rcout << ", " << myalleles[k];
//      }
//      Rcpp::Rcout << "\n";

      // Concatenate allele counts into a comma delimited string.
      std::string allele_string = std::to_string(myalleles[0]);
      for(k=1; k<myalleles.size(); k++){
        allele_string = allele_string + ",";
        allele_string = allele_string + std::to_string(myalleles[k]);
      }
//      Rcpp::Rcout << "allele_string: " << allele_string << "\n";
      allele_counts[i] = allele_string;


      // Sum all alleles.
      int nalleles = myalleles[0];
      for(j=1; (unsigned)j < myalleles.size(); j++){
        nalleles = nalleles + myalleles[j];
      }
//      Rcpp::Rcout << "nalleles: " << nalleles << "\n";
//      Rcpp::Rcout << "\n";

      // Stats.
      double He = 1;
      He = He - pow(double(myalleles[0])/double(nalleles), myalleles.size());
      for(j=1; (unsigned)j < myalleles.size(); j++){
        He = He - pow(double(myalleles[j])/double(nalleles), myalleles.size());
      }
      Hes[i] = He;
      Nes[i] = 1/(1-He);
    } else { // Missing variant (row=NA)
      nsample[i] = NA_INTEGER;     
    }
  }

  return Rcpp::DataFrame::create(var_info, 
      _["n"]=nsample, 
      _["Allele_counts"]=allele_counts,
      _["He"]=Hes,
      _["Ne"]=Nes
  );
}
