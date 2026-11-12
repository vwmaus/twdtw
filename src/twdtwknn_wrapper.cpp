#include <Rcpp.h>
#include <algorithm>
#include <vector>
using namespace Rcpp;

// Forward declaration of the Fortran functions
extern "C" {

  void twdtwknn_weights(double* TX, int* NT, double* TY, int* NY, double* TW,
                        double* LB, double* CL, double* W);

  void twdtwknn(double* XM, int* NS, int* NT, int* ND, int* P1, int* P2,
                double* YM, int* YO, int* NY, int* NP, double* WU, double* WL,
                int* K, int* KI, double* KD, int* IERR);
}

/**
 * Compute the TWDTW k-nearest neighbours of a set of time series
 *
 * @param XM   Array with the time series (NS,NT,ND), missing values are NA
 * @param XD   Dimensions of XM: number of time series, dates, and spectral dimensions
 * @param TX   Dates of the time series in units of the time cycle (NT)
 * @param YM   Matrix with the temporal patterns stacked by rows (NY,ND), without missing values
 * @param TY   Dates of the observations of the temporal patterns in units of the time cycle (NY)
 * @param YO   Row offset of each temporal pattern in YM (NP+1)
 * @param TW   Time-Weight parameters alpha and beta
 * @param LB   Maximum elapsed time to constrain TWDTW calculation
 * @param CL   The length of the time cycle
 * @param K    Number of nearest neighbours
 *
 * @return A list with the index of the K nearest patterns (NS,K) and their
 *         distance (NS,K). Time series without observations have NA.
 */
// [[Rcpp::export]]
List twdtw_knn_f90(NumericVector XM, IntegerVector XD, NumericVector TX,
                   NumericMatrix YM, NumericVector TY, IntegerVector YO,
                   NumericVector TW, double LB, double CL, int K) {

  // The Fortran code does not check array bounds
  if (XD.size() != 3) {
    stop("XD must have the number of time series, dates, and dimensions");
  }
  int NS = XD[0], NT = XD[1], ND = XD[2];
  int NY = YM.nrow(), NP = YO.size() - 1;

  if (NS < 0 || NT < 1 || ND < 1 ||
      XM.size() != static_cast<R_xlen_t>(NS) * NT * ND || TX.size() != NT) {
    stop("The dimensions of the time series do not match");
  }

  if (NP < 1 || YM.ncol() != ND || TY.size() != NY || YO[0] != 0 || YO[NP] != NY) {
    stop("The dimensions of the temporal patterns do not match");
  }

  for (int i = 0; i < NP; ++i) {
    if (YO[i + 1] <= YO[i]) {
      stop("Temporal patterns must have at least one observation");
    }
  }

  if (TW.size() != 2) {
    stop("TW must have the two parameters of the logistic time weight");
  }

  if (K < 1 || K > NP) {
    stop("K must be between one and the number of temporal patterns");
  }

  // Time weights are the same for all time series. The weights constrained by
  // the maximum elapsed time are only needed if the constraint is active
  double INF = R_PosInf;
  std::vector<double> WU(static_cast<std::size_t>(NT) * NY);
  twdtwknn_weights(TX.begin(), &NT, TY.begin(), &NY, TW.begin(), &INF, &CL, WU.data());

  std::vector<double> WC;
  double* WL = WU.data();
  if (LB != R_PosInf) {
    WC.resize(WU.size());
    twdtwknn_weights(TX.begin(), &NT, TY.begin(), &NY, TW.begin(), &LB, &CL, WC.data());
    WL = WC.data();
  }

  // Time series without observations are not modified by the Fortran code
  IntegerMatrix KI(NS, K);
  NumericMatrix KD(NS, K);
  std::fill(KI.begin(), KI.end(), NA_INTEGER);
  std::fill(KD.begin(), KD.end(), NA_REAL);

  // Process blocks of time series of about the same computational cost
  // to check for user interrupts between the blocks
  double cost = static_cast<double>(NY) * NT * ND;
  R_xlen_t block = static_cast<R_xlen_t>(std::max(1.0, std::min(static_cast<double>(NS), 2e7 / cost)));

  int IERR = 0;
  for (R_xlen_t first = 1; first <= NS; first += block) {
    int P1 = static_cast<int>(first);
    int P2 = static_cast<int>(std::min<R_xlen_t>(NS, first + block - 1));
    twdtwknn(XM.begin(), &NS, &NT, &ND, &P1, &P2, YM.begin(), YO.begin(), &NY, &NP,
             WU.data(), WL, &K, KI.begin(), KD.begin(), &IERR);
    if (IERR != 0) {
      stop("Could not allocate memory to compute TWDTW");
    }
    checkUserInterrupt();
  }

  return List::create(Named("index") = KI, Named("distance") = KD);
}
