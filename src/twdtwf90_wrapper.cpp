#include <Rcpp.h>
#include <functional>
using namespace Rcpp;

// Forward declaration of the Fortran functions
extern "C" {

  double logistic_tw(double* DIST, double* TD, double* TW1, double* TW2);

  void twdtwf90gt_(double* XM, double* YM, double* CM, int* DM, int* VM, int* SM,
                   int* N, int* M, int* D, int* NS, double* TW, double* LB, int* JB, double* CL);

  void twdtwf90(double* XM, double* YM, double* CM, int* DM, int* VM, int* SM,
                int* N, int* M, int* D, int* NS, double* TW, double* LB, int* JB,
                double* CL, double (*callback_func)(double*, double*, double*, double*));
}

// Define the callback function pointer type for time weight function
typedef double (*CallbackFunc)(double*, double*, double*, double*);

// Define a function object for the callback function
struct CallbackFuncObject {
  Function tw_r_fun;

  CallbackFuncObject(Function input_tw_r_fun) : tw_r_fun(input_tw_r_fun) {}

  double call(double* x, double* y, double* z, double* w) const {
    NumericVector result = tw_r_fun(*x, *y, *z, *w);  // Call the R function
    return result[0];  // Extract the first element as the result
  }
};

typedef std::function<double(double*, double*, double*, double*)> FuncType;

// Global pointer to CallbackFuncObject
CallbackFuncObject* gCallbackFuncObject = nullptr; // initialize it to nullptr

// Bridge function to convert the callback
extern "C" double callback_bridge(double* x, double* y, double* z, double* w) {
  return gCallbackFuncObject->call(x, y, z, w);
}

// Wrapper Fortran functions
// [[Rcpp::export]]
void twdtw_f90gt(NumericMatrix XM, NumericMatrix YM, NumericMatrix CM, IntegerMatrix DM,
                 IntegerMatrix VM, IntegerMatrix SM, int N, int M, int D, int NS,
                 NumericVector TW, double LB, IntegerVector JB, double CL) {
  twdtwf90gt_(XM.begin(), YM.begin(), CM.begin(), DM.begin(), VM.begin(), SM.begin(),
              &N, &M, &D, &NS, TW.begin(), &LB, JB.begin(), &CL);
}

// [[Rcpp::export]]
void twdtw_f90(NumericMatrix XM, NumericMatrix YM, NumericMatrix CM, IntegerMatrix DM,
               IntegerMatrix VM, IntegerMatrix SM, int N, int M, int D, int NS,
               NumericVector TW, double LB, IntegerVector JB, double CL,
               Rcpp::Nullable<Rcpp::Function> tw_r_fun = R_NilValue) {

  // Check if the R function is null or not
  bool is_tw_r_fun_null = tw_r_fun.isNull();

  // If the global callback object exists, delete it before assigning a new one
  if (gCallbackFuncObject != nullptr) {
    delete gCallbackFuncObject;
    gCallbackFuncObject = nullptr;
  }

  if (is_tw_r_fun_null) {
    // TW is NULL, handle this situation
    twdtwf90(XM.begin(), YM.begin(), CM.begin(), DM.begin(), VM.begin(), SM.begin(),
             &N, &M, &D, &NS, TW.begin(), &LB, JB.begin(), &CL, logistic_tw);
  } else {
    Function tw_r_fun_func(tw_r_fun);
    // Allocate the CallbackFuncObject on heap
    gCallbackFuncObject = new CallbackFuncObject(tw_r_fun_func);
    twdtwf90(XM.begin(), YM.begin(), CM.begin(), DM.begin(), VM.begin(), SM.begin(),
             &N, &M, &D, &NS, TW.begin(), &LB, JB.begin(), &CL, callback_bridge);

    // Delete the CallbackFuncObject immediately after using it
    delete gCallbackFuncObject;
    gCallbackFuncObject = nullptr;
  }
}
