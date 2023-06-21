#'
#' Calculate TWDTW Distance between Time Series
#'
#' This function calculates the Time-Weighted Dynamic Time Warping (TWDTW) distance
#' between two time series.
#'
#' @param x A tibble, data.frame, or data.table representing the first time series.
#' @param y A tibble, data.frame, or data.table representing the second time series.
#' @param tw Numeric vector of length 2 representing the time-weight parameters alpha and beta.
#'   Defaults to c(100, 1) (set weight to zero).
#' @param step_matrix A matrix specifying the step pattern for the TWDTW algorithm.
#'   Defaults to symmetric1.
#' @param index_column The column name of the time index. Defaults to "date".
#' @param lower_band Logical indicating whether to constrain the TWDTW calculation
#'   to the lower band given by the time-weight parameter beta. Defaults to TRUE.
#' @param all_matches Logical indicating whether to find all matches within the
#'   TWDTW matrix. Defaults to FALSE.
#' @param twdtw_version A string identifying the version of twdtw implementation.
#' Options are 'f77' for Fortran 77, 'f90' for Fortran 90, 'cpp' for C++. Defaults to 'f77'.
#'
#' @return A numeric value representing the TWDTW distance between the two time series.
#'
#' @examples
#'
#' n <- 20
#' t <- seq(0, pi, length.out = n)
#'
#' x <- data.frame(date = seq(as.Date("2020-01-01"), by = "day", length.out = n),
#'              value = sin(t)*2 + runif(n))
#'
#' y <- data.frame(date = seq(as.Date("2020-01-05"), by = "day", length.out = n),
#'              value = sin(t)*2 + runif(n))
#'
#' rbenchmark::benchmark(
#' f90=twdtw(x, y, tw = c(-.1,50), all_matches = TRUE, twdtw_version = 'f90'),
#' fgt=twdtw(x, y, tw = c(-.1,50), all_matches = TRUE, twdtw_version = 'f90gt'),
#' cpp=twdtw(x, y, tw = c(-.1,50), all_matches = TRUE, twdtw_version = 'cpp'),
#' replications = 10000
#' )
#'
#' @export
twdtw <- function(x, y, tw = c(100, 1), step_matrix = symmetric1, twdtw_version = 'f77',
                  index_column = 'date', lower_band = TRUE, all_matches = FALSE) {

  # The dimensions of the time series must match
  if (!setequal(names(y), names(x))) {
    stop("The dimensions (columns) of the time series 'x' and 'y' do not match.")
  }

  # Convert dates to numeric
  if(inherits(x[,index_column], "Date")){
    x[,index_column] <- as.numeric(format(x[,index_column], "%j"))
  }
  if(inherits(y[,index_column], "Date")){
    y[,index_column] <- as.numeric(format(y[,index_column], "%j"))
  }

  # Position time index at the first column
  x <- x[, c(index_column, setdiff(names(x), index_column)), drop = FALSE]

  # Sort dimensions of y according to x
  y <- y[, names(x), drop = FALSE]

  # Initialize data with correct dimensions and types
  N <- as.integer(nrow(y))
  M <- as.integer(nrow(x))
  D <- as.integer(ncol(y))
  XM <- matrix(as.double(as.matrix(x)), N, D)
  YM <- matrix(as.double(as.matrix(y)), N, D)
  CM <- matrix(0, nrow = N+1, ncol = M)
  DM <- matrix(0L, nrow = N+1, ncol = M)
  VM <- matrix(0L, nrow = N+1, ncol = M)
  JB <- as.integer(rep(0, N))
  SM <- matrix(as.integer(step_matrix), nrow(step_matrix), ncol(step_matrix))
  NS <- as.integer(nrow(SM))
  TW = as.double(tw)
  LB = as.logical(lower_band)

  # Get the function version using its name
  fn <- get(c('twdtw_f90',
              'twdtw_f90gt',
              'twdtw_cpp')[twdtw_version == c('f90', 'f90gt', 'cpp')])

  # Call twdtw
  fn(XM, YM, CM, DM, VM, SM, N, M, D, NS, TW, LB, JB)

  b <- JB[JB!=0]
  a <- VM[-1,][N,b]
  d <- CM[-1,][N,b]
  candidates <- matrix(c(a, b, d), ncol = 3, byrow = F)

  if (all_matches) {
    return(candidates)
  }

  candidates

  return(candidates)
}
