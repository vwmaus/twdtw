#include <cmath>
#include <limits>
#include <algorithm>


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
double distancecpp(const double* YM, const double* XM, int N, int M, int D, int I, int J, const double* TW, double TD) {
  double BD, CD, dist = 0.0;

  CD = 0.0;
  for (int K = 0; K < D; K++) {
    BD = YM[(I-1)*D + K] - XM[(J-1)*D + K];
    CD += BD * BD;
  }

  dist = std::sqrt(CD);
  dist += 1.0 / (1.0 + std::exp(TW[0] * (TD - TW[1])));

  return dist;
}



/**
 * Compute elapsed time in days
 *
 * @param X Time difference in days
 * @return Elapsed time in days
 */
double ellapsedcpp(double X) {
  const double PC = 366.0;
  const double HPC = PC / 2;
  double elapsed = 0.0;

  // Compute elapsed time difference
  elapsed = std::sqrt(X * X);

  // Correct elapsed time with year cycle
  if (elapsed > HPC) {
    elapsed = PC - elapsed;
  }

  elapsed = std::abs(elapsed);

  return elapsed;
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
void twdtwcpp(const double* XM, const double* YM, double* CM, int* DM, int* VM, const int* SM,
              const int N, const int M, const int D, const int NS, const double* TW, const bool LB,
              int* JB) {
  const double INF = std::numeric_limits<double>::infinity();
  const int ZERO = 0;
  const int ONE = 1;
  // const double NAN = std::numeric_limits<double>::quiet_NaN();

  double* CP = new double[NS];
  int JM = 0;

  VM[0] = 1;

  // Initialize the first row and col of the matrices
  for (int I = 2; I <= N + 1; ++I) {
    double TD = ellapsedcpp(YM[(I - 1) * D] - XM[0]);
    CM[I * M] = CM[(I - 1) * M] + distancecpp(YM, XM, N, M, D, I - 1, 1, TW, TD);
    DM[I * M] = 3;
    VM[I * M] = 1;
  }

  for (int J = 2; J <= M; ++J) {
    double TD = ellapsedcpp(YM[D] - XM[J * D]);
    CM[2 * M + J - 1] = CM[2 * M + J - 2] + distancecpp(YM, XM, N, M, D, 1, J, TW, TD);
    DM[J - 1] = 2;
    VM[J - 1] = J;
  }

  // Compute cumulative cost matrix
  int J = 2;
  while (J <= M) {
    int I = 2;
    while (I <= N + 1) {
      // Calculate local distance
      // the call takes I-1 because local matrix has an additional row at the beginning
      double TD = ellapsedcpp(YM[(I - 1) * D] - XM[J * D]);
      if (LB && (TD > TW[1])) {
        CM[I * M + J - 1] = INF;
        DM[I * M + J - 1] = -ONE;
        VM[I * M + J - 1] = ZERO;
      } else {
        CM[I * M + J - 1] = distancecpp(YM, XM, N, M, D, I - 1, J, TW, TD);
      }
      // Initialize list of step cost
      std::fill(CP, CP + NS, NAN);
      for (int K = 0; K < NS; ++K) {
        int PK = SM[K];
        int IL = I - SM[K + 2 * NS];
        int JL = J - SM[K + 3 * NS];
        if (IL > ZERO && JL > ZERO) {
          double W = SM[K + NS];
          if (W == -ONE) {
            CP[PK - 1] = CM[(IL - 1) * M + JL - 1];
          } else {
            CP[PK - 1] += CM[(IL - 1) * M + JL - 1] * W;
          }
        }
      }
      int KMIN = -ONE;
      double VMIN = INF;
      for (int K = 0; K < NS; ++K) {
        int PK = SM[K];
        if (std::isfinite(CP[PK - 1]) && CP[PK - 1] < VMIN) {
          KMIN = PK;
          VMIN = CP[PK - 1];
        }
      }
      if (KMIN > -ONE) {
        CM[I * M + J - 1] = VMIN;
        DM[I * M + J - 1] = KMIN;
        VM[I * M + J - 1] = VM[(I - SM[KMIN - 1 + NS]) * M + J - SM[KMIN - 1 + 2 * NS] - 1];
      }
      ++I;
    }
    ++J;
  }

  J = 1;
  int K = ZERO;
  while (J <= M) {
    if (VM[(N + 1) * M + J - 1] != ZERO) {
      if (K == ZERO) {
        K = 1;
        JB[K - 1] = J;
        JM = VM[(N + 1) * M + J - 1];
      } else if (VM[(N + 1) * M + J - 1] != JM) {
        ++K;
        JB[K - 1] = J;
        JM = VM[(N + 1) * M + J - 1];
      } else if (CM[(N + 1) * M + J - 1] < CM[(N + 1) * M + JB[K - 1] - 1]) {
        JB[K - 1] = J;
      }
    }
    ++J;
  }

  delete[] CP;
}

