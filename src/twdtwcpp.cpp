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
 * @param N    Number of rows in CM, DM, and VM - time series
 * @param M    Number of columns CM, DM, and VM - temporal profile
 * @param D    Number of spectral dimensions including time in XM and YM
 * @param TW   Time-Weight parameters alpha and beta
 * @param LB   Constrain TWDTW calculation to band given by TW(2)
 * @param JB   Output array of starting points
 * @param CL   The length of the time cycle
 */
// [[Rcpp::export]]
void twdtw_cpp(const NumericMatrix& XM, const NumericMatrix& YM, NumericMatrix& CM, IntegerMatrix& DM, IntegerMatrix& VM,
               int N, int M, int D, const NumericVector& TW, double LB, IntegerVector& JB, double CL) {
  const double INF = std::numeric_limits<double>::infinity();
  const int ZERO = 0;
  const int ONE = 1;

  int JM = 0;

  // Initialize VM matrix
  VM(0, 0) = 1; // R index

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


  // Compute cumulative cost matrix
  int J = 2;
  while (J <= M) {
    VM(0,J-1) = J + 1;
    int I = 2;
    while (I <= N + 1) {
      double TD = std::fabs(YM(I - 2, 0) - XM(J - 1, 0));
      TD = std::min(TD, CL - TD);
      if (TD > LB) {
        CM(I - 1, J - 1) = INF;
        DM(I - 1, J - 1) = -ONE;
        VM(I - 1, J - 1) = ZERO;
      } else {
        double DIST = 0.0;
        for (int K = 1; K < D; ++K) {
          DIST += std::pow(YM(I - 2, K) - XM(J - 1, K), 2);
        }
        double CP = sqrt(DIST) + 1.0 / (1.0 + std::exp(-TW[0] * (TD - TW[1])));
        CM(I - 1, J - 1) = CP + CM(I - 2, J - 2);
        DM(I - 1, J - 1) = ONE;
        VM(I - 1, J - 1) = VM(I - 2, J - 2);

        double ST = CP + CM(I - 1, J - 2);
        if (ST < CM(I - 1, J - 1)) {
          DM(I - 1, J - 1) = 2;
          CM(I - 1, J - 1) = ST;
          VM(I - 1, J - 1) = VM(I - 1, J - 2);
        }

        ST = CP + CM(I - 2, J - 1);
        if (ST < CM(I - 1, J - 1)){
            DM(I - 1, J - 1) = 3;
            CM(I - 1, J - 1) = ST;
            VM(I - 1, J - 1) = VM(I - 2, J - 1);
        }
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

