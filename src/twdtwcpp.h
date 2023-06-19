#ifndef TWDTW_FUNCTIONS_H
#define TWDTW_FUNCTIONS_H

double ellapsedcpp(double X);
double distancecpp(const double* YM, const double* XM, int N, int M, int D, int I, int J,
                   const double* TW, double TD);
void twdtwcpp(const double* XM, const double* YM, double* CM, int* DM, int* VM, const int* SM,
              const int N, const int M, const int D, const int NS, const double* TW, const bool LB,
              int* JB);

#endif
