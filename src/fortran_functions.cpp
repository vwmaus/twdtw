#include <Rcpp.h>
using namespace Rcpp;

// Forward declaration of the Fortran functions
extern "C" {

  void twdtwf90gt_(double* XM, double* YM, double* CM, int* DM, int* VM, int* SM,
                   int* N, int* M, int* D, int* NS, double* TW, bool* LB, int* JB);

  void twdtwf90_(double* XM, double* YM, double* CM, int* DM, int* VM, int* SM,
                 int* N, int* M, int* D, int* NS, double* TW, bool* LB, int* JB);

}

// Wrapper functions
// [[Rcpp::export]]
void twdtw_f90gt(NumericMatrix XM, NumericMatrix YM, NumericMatrix CM, IntegerMatrix DM,
                 IntegerMatrix VM, IntegerMatrix SM, int N, int M, int D, int NS,
                 NumericVector TW, bool LB, IntegerVector JB) {
  twdtwf90gt_(XM.begin(), YM.begin(), CM.begin(), DM.begin(), VM.begin(), SM.begin(),
              &N, &M, &D, &NS, TW.begin(), &LB, JB.begin());
}

// [[Rcpp::export]]
void twdtw_f90(NumericMatrix XM, NumericMatrix YM, NumericMatrix CM, IntegerMatrix DM,
               IntegerMatrix VM, IntegerMatrix SM, int N, int M, int D, int NS,
               NumericVector TW, bool LB, IntegerVector JB) {
  twdtwf90_(XM.begin(), YM.begin(), CM.begin(), DM.begin(), VM.begin(), SM.begin(),
            &N, &M, &D, &NS, TW.begin(), &LB, JB.begin());
}
