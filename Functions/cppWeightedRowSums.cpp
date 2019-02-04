#include <Rcpp.h>
using namespace Rcpp;

// [[Rcpp::export]]
NumericMatrix cppWeightedRowSums(NumericMatrix X, IntegerVector itmax, double alpha){
  int NT = X.ncol();
  NumericMatrix result(X.nrow(), NT);
  result(_,0) = X(_,0);
  for(int j=1; j<NT; j++){
    result(_,j) =  X(_,j) + alpha * result(_,j-1);
  } 
  return result;
}
