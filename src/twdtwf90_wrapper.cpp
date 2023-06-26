#include <Rcpp.h>
#include <functional>
using namespace Rcpp;

// Forward declaration of the Fortran functions
extern "C" {

  double logistic_tw(double* DIST, double* TD, double* TW1, double* TW2);

  double tw_distance(double* Y, double* X, double* PAR, int* D, int* NPAR);

  void twdtwf90gt_(double* XM, double* YM, double* CM, int* DM, int* VM, int* SM,
                   int* N, int* M, int* D, int* NS, double* TW, double* LB, int* JB);

  void twdtwf90(double* XM, double* YM, double* CM, int* DM, int* VM, int* SM,
                int* N, int* M, int* D, int* NS, double* TW, double* LB, int* JB,
                int* NPAR, double (*callback_func)(double*, double*, double*, int*, int*));
}

// Define the callback function pointer type for time weight function
typedef double (*CallbackFunc)(double*, double*, double*, double*);

typedef double (*CallbackFunc2)(double*, double*, double*, int*, int*);

// Define a function object for the callback function
struct CallbackFuncObject {
  Function tw_r;

  CallbackFuncObject(Function input_tw_r) : tw_r(input_tw_r) {}

  double call(double* x, double* y, double* z, double* w) const {
    NumericVector result = tw_r(*x, *y, *z, *w);  // Call the R function
    return result[0];  // Extract the first element as the result
  }
};

struct CallbackFuncObject2 {
  Function tw_r;

  CallbackFuncObject2(Function input_tw_r) : tw_r(input_tw_r) {}

  double call(double* x, double* y, double* z, int* w, int* k) const {
    NumericVector result = tw_r(*x, *y, *z, *w, *k);  // Call the R function
    return result[0];  // Extract the first element as the result
  }
};

typedef std::function<double(double*, double*, double*, double*)> FuncType;

// Global pointer to CallbackFuncObject
CallbackFuncObject* gCallbackFuncObject = nullptr; // initialize it to nullptr
CallbackFuncObject2* gCallbackFuncObject2 = nullptr; // initialize it to nullptr

// Bridge function to convert the callback
extern "C" double callback_bridge(double* x, double* y, double* z, double* w) {
  return gCallbackFuncObject->call(x, y, z, w);
}

extern "C" double callback_bridge2(double* x, double* y, double* z, int* w, int* k) {
  return gCallbackFuncObject2->call(x, y, z, w, k);
}

// Wrapper Fortran functions
// [[Rcpp::export]]
void twdtw_f90gt(NumericMatrix XM, NumericMatrix YM, NumericMatrix CM, IntegerMatrix DM,
                 IntegerMatrix VM, IntegerMatrix SM, int N, int M, int D, int NS,
                 NumericVector TW, double LB, IntegerVector JB) {
  twdtwf90gt_(XM.begin(), YM.begin(), CM.begin(), DM.begin(), VM.begin(), SM.begin(),
              &N, &M, &D, &NS, TW.begin(), &LB, JB.begin());
}

// [[Rcpp::export]]
void twdtw_f90(NumericMatrix XM, NumericMatrix YM, NumericMatrix CM, IntegerMatrix DM,
               IntegerMatrix VM, IntegerMatrix SM, int N, int M, int D, int NS,
               NumericVector TW, double LB, int NPAR, IntegerVector JB, Rcpp::Nullable<Rcpp::Function> tw_r = R_NilValue) {

  // Check if the R function is null or not
  bool is_tw_r_null = tw_r.isNull();

  // If the global callback object exists, delete it before assigning a new one
  if (gCallbackFuncObject2 != nullptr) {
    delete gCallbackFuncObject2;
    gCallbackFuncObject2 = nullptr;
  }

  if (is_tw_r_null) {
    // TW is NULL, handle this situation
    twdtwf90(XM.begin(), YM.begin(), CM.begin(), DM.begin(), VM.begin(), SM.begin(),
             &N, &M, &D, &NS, TW.begin(), &LB, JB.begin(), &NPAR, tw_distance);
  } else {
    Function tw_r_func(tw_r);
    // Allocate the CallbackFuncObject on heap
    gCallbackFuncObject2 = new CallbackFuncObject2(tw_r_func);
    twdtwf90(XM.begin(), YM.begin(), CM.begin(), DM.begin(), VM.begin(), SM.begin(),
             &N, &M, &D, &NS, TW.begin(), &LB, JB.begin(), &NPAR, callback_bridge2);

    // Delete the CallbackFuncObject immediately after using it
    delete gCallbackFuncObject2;
    gCallbackFuncObject2 = nullptr;
  }
}
