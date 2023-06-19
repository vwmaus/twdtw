#include <R.h>
#include <Rinternals.h>
#include <R_ext/Rdynload.h>
#include <stdbool.h>
#include "twdtwcpp.h"

/* Declare Fortran functions */
extern void F77_NAME(twdtw)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(twdtw90gt)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);
extern void F77_NAME(twdtw90)(void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *, void *);

/* Declare static array with information about the Fortran and C++ functions */
static const R_FortranMethodDef FortranEntries[] = {
  {"twdtw_f77", (DL_FUNC) &F77_NAME(twdtw), 13},
  {"twdtw_f90gt", (DL_FUNC) &F77_NAME(twdtw90gt), 13},
  {"twdtw_f90", (DL_FUNC) &F77_NAME(twdtw90), 13},
  {NULL, NULL, 0}
};

/* Declare C++ functions */
void twdtwcpp(const double* XM, const double* YM, double* CM, int* DM, int* VM, const int* SM,
              const int N, const int M, const int D, const int NS, const double* TW, const bool LB,
              int* JB);

/* Register C++ function as callable from R */
void register_twdtwcpp() {
  R_RegisterCCallable("twdtw_cpp", "twdtwcpp", (DL_FUNC) &twdtwcpp);
}

/* Initialize the package DLL */
void R_init_twdtw(DllInfo *dll) {
  R_registerRoutines(dll, NULL, NULL, FortranEntries, NULL);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, FALSE);

  register_twdtwcpp();  // Register the C++ function
}
