#include <Rcpp.h>
#include <algorithm>
#include <cmath>
#include <vector>
using namespace Rcpp;

// Forward declaration of the Fortran functions
extern "C" {

  void twdtwknn_weights(double* TX, int* NT, double* TY, int* NY, double* TW,
                        double* LB, double* CL, double* W);

  void twdtwknn(double* XM, int* NS, int* NT, int* ND, int* P1, int* P2,
                double* YM, int* YO, int* NY, int* NP, double* WU, double* WL,
                int* K, int* KI, double* KD, int* IERR);

  void twdtwknn_dist(double* XM, int* NS, int* NT, int* ND, int* P1, int* P2,
                     double* YM, int* NY, int* Y1, int* Y2, double* DM, int* OK);

  void twdtwknn_cost(double* CM, int* NB, int* NT, int* NI, int* OK, double* EL,
                     double* LB, int* IP, int* K, int* NS, int* P1,
                     int* KN, int* KI, double* KD, int* IERR);
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

/**
 * Compute the TWDTW k-nearest neighbours of a set of time series using a
 * time-weight function defined in R
 *
 * The time-weight function is called once for each temporal pattern and block
 * of time series. It receives two vectors with the local distance and the
 * elapsed time between each observation of the time series and each
 * observation of the pattern, and must return a vector with their cost. The
 * function must be vectorised, i.e. it must return the same values as one
 * call for each pair of distance and elapsed time.
 *
 * @param XM   Array with the time series (NS,NT,ND), missing values are NA
 * @param XD   Dimensions of XM: number of time series, dates, and spectral dimensions
 * @param TX   Dates of the time series in units of the time cycle (NT)
 * @param YM   Matrix with the temporal patterns stacked by rows (NY,ND), without missing values
 * @param TY   Dates of the observations of the temporal patterns in units of the time cycle (NY)
 * @param YO   Row offset of each temporal pattern in YM (NP+1)
 * @param TW   Time-weight function of the local distance and the elapsed time
 * @param LB   Maximum elapsed time to constrain TWDTW calculation
 * @param CL   The length of the time cycle
 * @param K    Number of nearest neighbours
 *
 * @return A list with the index of the K nearest patterns (NS,K) and their
 *         distance (NS,K). Time series without observations have NA.
 */
// [[Rcpp::export]]
List twdtw_knn_fun_f90(NumericVector XM, IntegerVector XD, NumericVector TX,
                       NumericMatrix YM, NumericVector TY, IntegerVector YO,
                       Function TW, double LB, double CL, int K) {

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

  int NI_MAX = 0;
  for (int i = 0; i < NP; ++i) {
    if (YO[i + 1] <= YO[i]) {
      stop("Temporal patterns must have at least one observation");
    }
    NI_MAX = std::max(NI_MAX, YO[i + 1] - YO[i]);
  }

  if (K < 1 || K > NP) {
    stop("K must be between one and the number of temporal patterns");
  }

  // Elapsed time between each date of the time series and each observation of
  // the temporal patterns, calculated as in twdtwf90
  std::vector<double> EL(static_cast<std::size_t>(NT) * NY);
  for (int i = 0; i < NY; ++i) {
    for (int j = 0; j < NT; ++j) {
      double TD = std::fabs(TY[i] - TX[j]);
      EL[static_cast<std::size_t>(i) * NT + j] = std::min(TD, CL - TD);
    }
  }

  // Time series without observations are not modified by the Fortran code
  IntegerMatrix KI(NS, K);
  NumericMatrix KD(NS, K);
  std::fill(KI.begin(), KI.end(), NA_INTEGER);
  std::fill(KD.begin(), KD.end(), NA_REAL);
  std::vector<int> KN(NS, 0);

  // Process blocks of time series, such that the vectors passed to the
  // time-weight function have at most about two million values
  double cells = static_cast<double>(NT) * NI_MAX;
  R_xlen_t block = static_cast<R_xlen_t>(std::max(1.0, std::min(static_cast<double>(NS), 2e6 / cells)));

  int IERR = 0;
  for (R_xlen_t first = 1; first <= NS; first += block) {
    int P1 = static_cast<int>(first);
    int P2 = static_cast<int>(std::min<R_xlen_t>(NS, first + block - 1));
    int NB = P2 - P1 + 1;
    std::vector<int> OK(static_cast<std::size_t>(NB) * NT);

    for (int ip = 0; ip < NP; ++ip) {
      int Y1 = YO[ip] + 1, Y2 = YO[ip + 1];
      int NI = Y2 - Y1 + 1;
      int IP = ip + 1;
      R_xlen_t n = static_cast<R_xlen_t>(NB) * NT * NI;

      // Local distance and elapsed time of each pair of observations
      NumericVector DM(n), EM(n);
      twdtwknn_dist(XM.begin(), &NS, &NT, &ND, &P1, &P2, YM.begin(), &NY, &Y1, &Y2,
                    DM.begin(), OK.data());

      const double* EP = EL.data() + static_cast<std::size_t>(Y1 - 1) * NT;
      R_xlen_t pos = 0;
      for (int i = 0; i < NI; ++i) {
        for (int j = 0; j < NT; ++j) {
          double TD = EP[static_cast<std::size_t>(i) * NT + j];
          for (int p = 0; p < NB; ++p) {
            EM[pos++] = TD;
          }
        }
      }

      // Call the time-weight function once for all pairs of observations
      Shield<SEXP> res(TW(DM, EM));
      if (!Rf_isNumeric(res) || Rf_xlength(res) != n) {
        stop("The time-weight function must return a numeric vector with one value for each pair of distance and elapsed time");
      }
      NumericVector CM(res);

      // The cost of the observations that are used cannot be missing
      pos = 0;
      for (int i = 0; i < NI; ++i) {
        for (int j = 0; j < NT; ++j) {
          for (int p = 0; p < NB; ++p) {
            if (OK[static_cast<std::size_t>(j) * NB + p] == 1 && ISNAN(CM[pos])) {
              stop("The time-weight function returned missing values");
            }
            ++pos;
          }
        }
      }

      twdtwknn_cost(CM.begin(), &NB, &NT, &NI, OK.data(), const_cast<double*>(EP), &LB,
                    &IP, &K, &NS, &P1, KN.data(), KI.begin(), KD.begin(), &IERR);
      if (IERR != 0) {
        stop("Could not allocate memory to compute TWDTW");
      }
    }
    checkUserInterrupt();
  }

  return List::create(Named("index") = KI, Named("distance") = KD);
}
