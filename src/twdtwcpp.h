#ifndef TWDTW_H
#define TWDTW_H

#include <Rcpp.h>

// Compute TWDTW distance using logistic weight
double distancecpp(const Rcpp::NumericVector& YM, const Rcpp::NumericVector& XM, int N, int M, int D, int I, int J, const Rcpp::NumericVector& TW);

// Compute elapsed time in days
double ellapsedcpp(double X);

// Compute TWDTW distance using logistic weight
void twdtw_cpp(const Rcpp::NumericMatrix& XM, const Rcpp::NumericMatrix& YM, Rcpp::NumericMatrix& CM, Rcpp::IntegerMatrix& DM, Rcpp::IntegerMatrix& VM,
               const Rcpp::IntegerMatrix& SM, int N, int M, int D, int NS, const Rcpp::NumericVector& TW, bool LB, Rcpp::IntegerVector& JB);

#endif
