#include <Rcpp.h>
#include <functional>
using namespace Rcpp;

// Forward declaration of the Fortran functions
extern "C" {

  void twdtwf90gt_(double* XM, double* YM, double* CM, int* DM, int* VM, int* SM,
                   int* N, int* M, int* D, int* NS, double* TW, bool* LB, int* JB);

  void twdtwf90(double* XM, double* YM, double* CM, int* DM, int* VM, int* SM,
                 int* N, int* M, int* D, int* NS, double* TW, bool* LB, int* JB, double (*callback_func)(double, double));
}

// Define the callback function pointer type for time weight function
typedef double (*CallbackFunc)(double, double);

// Define a function object for the callback function
struct CallbackFuncObject {
  Function tw_r;

  CallbackFuncObject(Function input_tw_r) : tw_r(input_tw_r) {}

  double call(double x, double y) const {
    NumericVector result = tw_r(x, y);  // Call the R function
    return result[0];  // Extract the first element as the result
  }
};

typedef std::function<double(double, double)> FuncType;

// Global pointer to CallbackFuncObject
CallbackFuncObject* gCallbackFuncObject;

// Bridge function to convert the callback
extern "C" double callback_bridge(double x, double y) {
  return gCallbackFuncObject->call(x, y);
}

// Wrapper Fortran functions
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
               NumericVector TW, bool LB, IntegerVector JB, Function tw_r) {

  // Create an instance of the callback function object
  CallbackFuncObject callback_obj(tw_r);
  gCallbackFuncObject = &callback_obj;

  // Call the Fortran wrapper function
  twdtwf90(XM.begin(), YM.begin(), CM.begin(), DM.begin(), VM.begin(), SM.begin(),
          &N, &M, &D, &NS, TW.begin(), &LB, JB.begin(), callback_bridge);
}
