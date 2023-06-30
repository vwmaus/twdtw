#include <Rcpp.h>
#include <algorithm> // for std::min
#include <cmath>     // for std::exp and std::sqrt
#include "twdtwcpp.h"
using namespace Rcpp;

/**
 * This function applies a logistic transformation on an input time difference (TD) and distance (DIST)
 * using a pair of weights (TW1, TW2).
 *
 * @param DIST : A double representing the distance parameter in the logistic transformation.
 * @param TD : A double representing the time difference parameter in the logistic transformation.
 * @param TW1 : A double representing the first weight parameter in the logistic transformation.
 * @param TW2 : A double representing the second weight parameter in the logistic transformation.
 *
 * @return : A double that is the result of applying the logistic transformation on DIST and TD
 * using the weights TW1 and TW2.
 */
double logistic_tw_cpp(double DIST, double TD, double TW1, double TW2) {
  return DIST + 1.0 / (1.0 + std::exp(-TW1 * (TD - TW2)));
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
 * @param CL   The length of the time cycle
 */
// [[Rcpp::export]]
void twdtw_cpp(const NumericMatrix& XM, const NumericMatrix& YM, NumericMatrix& CM, IntegerMatrix& DM, IntegerMatrix& VM,
               const IntegerMatrix& SM, int N, int M, int D, int NS, const NumericVector& TW, double LB, IntegerVector& JB, double CL) {
  const double INF = std::numeric_limits<double>::infinity();
  const int ZERO = 0;
  const int ONE = 1;

  NumericVector CP(NS);
  int JM = 0, ILMIN, JLMIN;
  IntegerVector IL(NS), JL(NS);

  // Initialize VM matrix
  VM(0, 0) = 1;

  // Initialize the first row and column of the matrices
  for (int I = 1; I <= N; ++I) {
    double TD = std::fabs(YM(I-1, 0) - XM(0, 0));
    TD = std::min(TD, CL - TD);
    double DIST = 0.0;
    for (int K = 1; K < D; ++K) {
      DIST += std::pow(YM(I-1, K) - XM(0, K), 2);
    }
    CM(I, 0) = CM(I-1, 0) + sqrt(DIST) + 1.0 / (1.0 + std::exp(-TW[0] * (TD - TW[1])));
    DM(I, 0) = 3;
    VM(I, 0) = 1;
  }

  for (int J = 1; J < M; ++J) {
    double TD = std::fabs(YM(1, 0) - XM(J, 0));
    TD = std::min(TD, CL - TD);
    double DIST = 0.0;
    for (int K = 1; K < D; ++K) {
      DIST += std::pow(YM(0, K) - XM(J, K), 2);
    }
    CM(1, J) = CM(1, J-1) + sqrt(DIST) + 1.0 / (1.0 + std::exp(-TW[0] * (TD - TW[1])));
    DM(0, J) = 2;
    VM(0, J) = J;
  }

  // Compute cumulative cost matrix
  int J = 2;
  while (J <= M) {
    int I = 2;
    while (I <= N + 1) {
      double TD = std::fabs(YM(I - 1, 0) - XM(J - 1, 0));
      TD = std::min(TD, CL - TD);
      if (TD > LB) {
        CM(I - 1, J - 1) = INF;
        DM(I - 1, J - 1) = -ONE;
        VM(I - 1, J - 1) = ZERO;
      } else {
        double DIST = 0.0;
        for (int K = 1; K < D; ++K) {
          DIST += std::pow(YM(I - 1, K) - XM(J - 1, K), 2);
        }
        CM(I - 1, J - 1) = sqrt(DIST) + 1.0 / (1.0 + std::exp(-TW[0] * (TD - TW[1])));
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

