#' @title Calculate Time-Weighted Dynamic Time Warping (TWDTW) distance
#'
#' @description
#' This function calculates the Time-Weighted Dynamic Time Warping (TWDTW) distance between two time series.
#'
#' @param x A data.frame or matrix representing time series.
#' @param y A data.frame or matrix representing a labelled time series (reference).
#' @param time_weight A numeric vector with lenght two (steepness and midpoint of logistic weight) or a function. See details.
#' @param cycle_length A character string or a numeric indicating the larger unit of time.
#' It must be one of "year", "month", "day", "hour", or "minute". It can also receive a numeric value when \code{time_scale} is numeric.
#' @param time_scale A character string or numeric indicating the smaller unit of time,
#' which is a division of the \code{cycle_length}. If \code{cycle_length} is "year", \code{time_scale} can be one of
#' "month", "day", "hour", "minute", "second". If \code{cycle_length} is "month", \code{time_scale} can be "day",
#' "hour", "minute", "second", and so on. It can also receive a numeric value when \code{cycle_length} is numeric.
#' @param index_column (optional) The column name of the time index for data.frame inputs. Defaults to "time".
#' For matrix input, an integer indicating the column with the time index. Defaults to 1.
#' @param dtw_step_matrix A matrix specifying the step pattern for the DTW algorithm. Defaults to \code{\link[dtw]{symmetric1}}.
#' @param max_elapsed Numeric value constraining the TWDTW calculation to the lower band given by a maximum elapsed time. Defaults to Inf.
#' @param output A character string defining the output. It must be one of 'distance', 'matches', 'internals'. Defaults to 'distance'.
#' 'distance' will return the lowest TWDTW distance between \code{x} and \code{y}.
#' 'matches' will return all matches within the TWDTW matrix. 'internals' will return all TWDTW internal data.
#' @param version A string identifying the version of TWDTW implementation.
#' Options are 'f90' for Fortran 90, 'f90goto' for Fortran 90 with goto statements,
#' or 'cpp' for C++ version. Defaults to 'f90'. See details.
#' @param ... ignore
#'
#'
#' @details TWDTW calculates a time-weighted version of DTW by modifying each element of the
#' DTW's local cost matrix (see details in Maus et al. (2016) and Maus et al. (2019)).
#' The default time weight is calculated using a logistic function
#' that adds a weight to each pair of observations in the time series \code{x} and \code{y}
#' based on the time difference between observations, such that
#'
#' \deqn{tw(dist_{i,j}) = dist_{i,j} + \frac{1}{1 + e^{-{\alpha} (el_{i,j} - {\beta})}}}
#'
#'Where:
#' \itemize{
#'  \item{\eqn{tw} is the time-weight function}
#'  \item{\eqn{dist_{i,j}} is the Euclidean distance between the i-th element of \code{x} and the j-th element of \code{y} in a multi-dimensional space}
#'  \item{\eqn{el_{i,j}} is the time elapsed between the i-th element of \code{x} and the j-th element of \code{y}}
#'  \item{\eqn{\alpha} and \eqn{\beta} are the steepness and midpoint of the logistic function, respectively}
#' }
#'
#' The logistic function is implemented as the default option in the C++ and Fortran versions of the code.
#' To use the native implementation, \eqn{\alpha} and \eqn{\beta} must be provided as a numeric vector of
#' length two using the argument \code{time_weight}. This implementation provides high processing performance.
#'
#' The \code{time_weight} argument also accepts a function defined in R, allowing the user to define a different
#' weighting scheme. However, passing a function to \code{time_weight} can degrade the processing performance,
#' i.e., it can be up to 3x slower than using the default logistic time-weight.
#'
#' A time-weight function passed to \code{time_weight} must receive two numeric arguments and return a
#' single numeric value. The first argument received is the Euclidean \eqn{dist_{i,j}} and the second
#' is the elapsed time \eqn{el_{i,j}}. For example,
#' \code{time_weight = function(dist, el) dist + 0.1*el} defines a linear weighting scheme with a slope of 0.1.
#'
#' The Fortran 90 versions of \code{twdtw} typically outperform the C++ version.
#' Additionally, the '`f90goto`' version, which utilizes `goto` statements,
#' tends to be slightly faster than the '`f90`' version, which uses only while and for loops.
#'
#' @return
#' If output = 'distance', a numeric value representing the TWDTW distance between the two time series.
#' If output = 'matches', a numeric matrix of all TWDTW matches.
#' For each match the starting index, ending index, and distance are returned.
#' If output = 'internals' a list of all TWDTW internal data is returned.
#'
#' @references
#' Maus, V., Camara, G., Cartaxo, R., Sanchez, A., Ramos, F. M., & de Moura, Y. M. (2016).
#' A Time-Weighted Dynamic Time Warping Method for Land-Use and Land-Cover Mapping.
#' IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 9(8), 3729-3739.
#' \url{https://doi.org/10.1109/JSTARS.2016.2517118}
#'
#' Maus, V., Camara, G., Appel, M., & Pebesma, E. (2019).
#' dtwSat: Time-Weighted Dynamic Time Warping for Satellite Image Time Series Analysis in R.
#' Journal of Statistical Software, 88(5), 1-31.
#' \url{https://doi.org/10.18637/jss.v088.i05}
#'
#' @examples
#'
#' # Create a time series
#' n <- 23
#' t <- seq(0, pi, length.out = n)
#' d <- seq(as.Date('2020-09-01'), length.out = n, by = "15 day")
#'
#' x <- data.frame(time = d,      v1 = sin(t)*2 + runif(n))
#'
#' # shift time by 30 days
#' y <- data.frame(time = d + 30, v1 = sin(t)*2 + runif(n))
#'
#' plot(x, type = "l", xlim = range(c(d, d + 5)))
#' lines(y, col = "red")
#'
#' # Calculate TWDTW distance between x and y using logistic weight
#' twdtw(x, y,
#'       cycle_length = 'year',
#'       time_scale = 'day',
#'       time_weight = c(steepness = 0.1, midpoint = 50))
#'
#' # Pass a generic time-weight function
#' twdtw(x, y,
#'       cycle_length = 'year',
#'       time_scale = 'day',
#'       time_weight = function(x,y) x + 1.0 / (1.0 + exp(-0.1 * (y - 50))))
#'
#' # Test other version
#' twdtw(x, y,
#'       cycle_length = 'year',
#'       time_scale = 'day',
#'       time_weight = c(steepness = 0.1, midpoint = 50),
#'       version = 'f90goto')
#'
#' twdtw(x, y,
#'       cycle_length = 'year',
#'       time_scale = 'day',
#'       time_weight = c(steepness = 0.1, midpoint = 50),
#'       version = 'cpp')
#'
#' @include convert_date_to_numeric.R
#'
#' @export
twdtw <- function(x, y, time_weight, cycle_length, time_scale, ...) {
  UseMethod("twdtw")
}

#' @export
#' @rdname twdtw
twdtw.data.frame <- function(x, y,
                             time_weight,
                             cycle_length,
                             time_scale,
                             index_column = 'time',
                             dtw_step_matrix = symmetric1,
                             max_elapsed = Inf,
                             output = 'distance',
                             version = 'f90', ...) {

  # Check that 'x' and 'y' are data.frames
  if (!inherits(y, "data.frame")) {
    stop("Both x and y need to be data.frames")
  }

  if (is.null(cycle_length)) stop("The 'cycle_length' argument is missing.")
  if (is.null(time_scale)) stop("The 'time_scale' argument is missing.")

  # The dimensions of the time series must match
  if (!setequal(names(y), names(x))) {
    stop("The dimensions (columns) of the time series 'x' and 'y' do not match.")
  }

  # The dimensions of the time series must match and 'index_column' must be present in both 'x' and 'y'
  if (!setequal(names(y), names(x)) & !index_column %in% names(x)) {
    stop("The dimensions (columns) of the time series 'x' and 'y' do not match.")
  }

  # Position 'index_column' at the first column
  x <- x[, c(index_column, setdiff(names(x), index_column)), drop = FALSE]

  # Sort dimensions of y according to x
  y <- y[, names(x), drop = FALSE]

  # Convert 'index_column' to numeric based on 'cycle_length' and 'time_scale'
  x[, index_column] <- convert_date_to_numeric(x[, index_column], cycle_length, time_scale)
  y[, index_column] <- convert_date_to_numeric(y[, index_column], cycle_length, time_scale)

  # call .twdtw function
  twdtw(x = as.matrix(x),
        y = as.matrix(y),
        time_weight = time_weight,
        cycle_length = cycle_length,
        time_scale = time_scale,
        index_column = 1,
        dtw_step_matrix = dtw_step_matrix,
        max_elapsed = max_elapsed,
        output = output,
        version = version)
}

#' @export
#' @rdname twdtw
twdtw.matrix <- function(x, y,
                         time_weight,
                         cycle_length,
                         time_scale = NULL,
                         index_column = 1,
                         dtw_step_matrix = symmetric1,
                         max_elapsed = Inf,
                         output = 'distance',
                         version = 'f90', ...) {

  # Check that 'x' and 'y' are matrix
  if (!inherits(y, "matrix")) {
    stop("Both x and y need to be data.frames")
  }

  # The dimensions of the time series must match
  if (ncol(x) != ncol(y)) {
    stop("The dimensions (columns) of the time series 'x' and 'y' do not match.")
  }

  if (!(is.function(time_weight) || (is.numeric(time_weight) && length(time_weight) == 2))) {
    stop("'time_weight' should be either a function or a numeric vector with length two")
  }

  # Check if 'cycle_length' is declared
  if (is.null(cycle_length)) stop("The 'cycle_length' argument is missing.")

  # Get maximum possible value that a specific time cycle and scale
  if (is.character(cycle_length)){
    if (is.null(time_scale)) stop("The 'time_scale' argument is missing for 'cycle_length' type character.")
    cycle_length <- calculate_max_cycle_length(cycle_length, time_scale)
  }

  # Position 'index_column' at the first column
  new_order <- c(index_column, setdiff(1:ncol(x), index_column))
  x <- x[, new_order, drop = FALSE]
  y <- y[, new_order, drop = FALSE]

  # call .twdtw function
  .twdtw(x = x,
         y = y,
         time_weight = time_weight,
         dtw_step_matrix = dtw_step_matrix,
         cycle_length = cycle_length,
         max_elapsed = max_elapsed,
         output = output,
         version = version)
}


# Internal function
.twdtw <- function(x, y,
                   time_weight,
                   cycle_length,
                   dtw_step_matrix = symmetric1,
                   max_elapsed = Inf,
                   output = 'distance',
                   version = 'f90', ...) {

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
  SM <- matrix(as.integer(dtw_step_matrix), nrow(dtw_step_matrix), ncol(dtw_step_matrix))
  NS <- as.integer(nrow(SM))
  CL <- as.double(cycle_length)
  LB <- as.double(max_elapsed)
  if (is.function(time_weight)){
    TW <- as.double(c(0.0, 0.0))
  } else {
    TW <- as.double(time_weight)
  }

  # Get the function version using its name
  fn <- get(c('twdtw_f90',
              'twdtw_f90gt',
              'twdtw_cpp')[version == c('f90', 'f90goto', 'cpp')])

  # Prepare arguments
  args <- list(XM, YM, CM, DM, VM, SM, N, M, D, NS, TW, LB, JB, CL)
  if (version == 'f90' & is.function(time_weight)){
    args$tw_r <- time_weight
  }

  # Call twdtw
  do.call(fn, args)
  b <- JB[JB!=0]

  switch (output,
    'distance' = min(CM[-1,][N,b]),
    'matches' = matrix(c(VM[-1,][N,b], b, CM[-1,][N,b]), ncol = 3, byrow = F),
    'internals' = list(XM, YM, CM, DM, VM, SM, N, M, D, NS, TW, LB, JB, CL)
  )

}
