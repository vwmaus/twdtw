#include <Rcpp.h>
#include "twdtwcpp.h"
using namespace Rcpp;


/**
 * Compute elapsed time in days
 *
 * @param X Time difference in days
 * @return Elapsed time in days
 */
double ellapsedcpp(double X) {
  const double PC = 366.0;
  const double HPC = PC / 2.0;
  double elapsed = std::abs(X);

  if (elapsed > HPC) {
    elapsed = PC - elapsed;
  }

  return elapsed;
}

/**
 * Compute TWDTW distance using logistic weight
 *
 * @param YM Matrix with the temporal profile (N,D)
 * @param XM Matrix with the time series (M,D)
 * @param N Number of rows in CM, DM, and VM - time series
 * @param M Number of columns CM, DM, and VM - temporal profile
 * @param D Number of spectral dimensions including time in XM and YM
 * @param I Single point in the time series to calculate the local distance
 * @param J Single point in the temporal profile to calculate the local distance
 * @param TW Time-Weight parameter alpha and beta
 * @param TD Time difference
 * @return TWDTW distance
 */
double distancecpp(const NumericMatrix& YM, const NumericMatrix& XM, int N, int M, int D, int I, int J, const NumericVector& TW) {
  double BD, dist = 0.0;

  for (int K = 1; K < D; ++K) {
    BD = YM((I-1), K) - XM((J-1), K);
    dist += BD * BD;
  }

  dist = std::sqrt(dist);
  dist += 1.0 / (1.0 + std::exp(TW(0) * (ellapsedcpp(YM((I-1), 0) - XM((J-1), 0)) - TW(1))));

  return dist;
}

/**
 * Compute TWDTW distance using logistic weight
 *
 * @param XM   Matrix with the time series (N,D)
 * @param YM   Matrix with the temporal profile (M,D)
 * @param CM   Output cumulative cost matrix
 * @param DM   Direction matrix
 * @param VM   Starting points matrix
 * @param SM   Matrix of step patterns
 * @param N    Number of rows in CM, DM, and VM - time series
 * @param M    Number of columns CM, DM, and VM - temporal profile
 * @param D    Number of spectral dimensions including time in XM and YM
 * @param NS   Number of rows in SM
 * @param TW   Time-Weight parameters alpha and beta
 * @param LB   Constrain TWDTW calculation to band given by TW(2)
 * @param JB   Output array of starting points
 */
// [[Rcpp::export]]
void twdtw_cpp(const NumericMatrix& XM, const NumericMatrix& YM, NumericMatrix& CM, IntegerMatrix& DM, IntegerMatrix& VM,
               const IntegerMatrix& SM, int N, int M, int D, int NS, const NumericVector& TW, double LB, IntegerVector& JB) {
  const double INF = std::numeric_limits<double>::infinity();
  const int ZERO = 0;
  const int ONE = 1;

  NumericVector CP(NS);
  int JM = 0, ILMIN, JLMIN;
  IntegerVector IL(NS), JL(NS);

  // Initialize VM matrix
  VM(0, 0) = 1;

  // Initialize the first row and column of the matrices
  for (int I = 2; I <= N + 1; ++I) {
    CM(I - 1, 0) = CM(I - 2, 0) + distancecpp(YM, XM, N, M, D, I - 1, 1, TW);
    DM(I - 1, 0) = 3;
    VM(I - 1, 0) = 1;
  }

  for (int J = 2; J <= M; ++J) {
    CM(1, J - 1) = CM(0, J - 2) + distancecpp(YM, XM, N, M, D, 1, J, TW);
    DM(0, J - 1) = 2;
    VM(0, J - 1) = J;
  }

  // Compute cumulative cost matrix
  int J = 2;
  while (J <= M) {
    int I = 2;
    while (I <= N + 1) {
      double TD = ellapsedcpp(YM(I - 1, 0) - XM(J - 1, 0));
      if (TD > LB) {
        CM(I - 1, J - 1) = INF;
        DM(I - 1, J - 1) = -ONE;
        VM(I - 1, J - 1) = ZERO;
      } else {
        CM(I - 1, J - 1) = distancecpp(YM, XM, N, M, D, I - 1, J, TW);
      }
      CP.fill(NA_REAL);
      for (int K = 0; K < NS; ++K) {
        int PK = SM(K, 0);
        IL(K) = I - SM(K, 1);
        JL(K) = J - SM(K, 2);
        if (IL(K) > ZERO && JL(K) > ZERO) {
          double W = SM(K, 3);
          if (W == -ONE) {
            CP(PK - 1) = CM(IL(K) - 1, JL(K) - 1);
          } else {
            CP(PK - 1) += CM(IL(K) - 1, JL(K) - 1) * W;
          }
        }
      }
      int KMIN = -ONE;
      double VMIN = INF;
      for (int K = 0; K < NS; ++K) {
        int PK = SM(K, 0);
        if (std::isfinite(CP(PK - 1)) && CP(PK - 1) < VMIN) {
          KMIN = PK;
          VMIN = CP(PK - 1);
          ILMIN = IL(K);
          JLMIN = JL(K);
        }
      }
      if (KMIN > -ONE) {
        CM(I - 1, J - 1) = VMIN;
        DM(I - 1, J - 1) = KMIN;
        VM(I - 1, J - 1) = VM(ILMIN - 1, JLMIN - 1);
      }
      ++I;
    }
    ++J;
  }

  J = 1;
  int K = ZERO;
  while (J <= M) {
    if (VM(N, J - 1) != ZERO) {
      if (K == ZERO) {
        K = 1;
        JB(K - 1) = J;
        JM = VM(N, J - 1);
      } else if (VM(N, J - 1) != JM) {
        ++K;
        JB(K - 1) = J;
        JM = VM(N, J - 1);
      } else if (CM(N, J - 1) < CM(N, JB(K - 1) - 1)) {
        JB(K - 1) = J;
      }
    }
    ++J;
  }
}

