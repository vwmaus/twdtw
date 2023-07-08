#ifndef TWDTW_H
#define TWDTW_H

#include <Rcpp.h>

// Define TWDTW distance
void twdtw_cpp(const Rcpp::NumericMatrix& XM, const Rcpp::NumericMatrix& YM, Rcpp::NumericMatrix& CM, Rcpp::IntegerMatrix& DM, Rcpp::IntegerMatrix& VM,
               int N, int M, int D, const Rcpp::NumericVector& TW, double LB, Rcpp::IntegerVector& JB, double CL);

#endif
