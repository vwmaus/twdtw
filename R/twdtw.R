#' @title Calculate Time-Weighted Dynamic Time Warping (TWDTW) distance
#'
#' @description
#' This function calculates the Time-Weighted Dynamic Time Warping (TWDTW) distance between two time series.
#'
#' @param x A data.frame or matrix representing time series.
#' @param y A data.frame or matrix representing a labelled time series (reference).
#' @param time_weight A numeric vector with lenght two (steepness and midpoint of logistic weight) or a function. See details.
#' @param time_cycle_length A character string or a numeric indicating the larger unit of time.
#' It must be one of "year", "month", "day", "hour", or "minute". It can also receive a numeric value when \code{time_cycle_scale} is numeric.
#' @param time_cycle_scale A character string or numeric indicating the smaller unit of time,
#' which is a division of the \code{time_cycle_length}. If \code{time_cycle_length} is "year", \code{time_cycle_scale} can be one of
#' "month", "day", "hour", "minute", "second". If \code{time_cycle_length} is "month", \code{time_cycle_scale} can be "day",
#' "hour", "minute", "second", and so on. It can also receive a numeric value when \code{time_cycle_length} is numeric.
#' @param index_column_name (optional) The column name of the time index for data.frame inputs. Defaults to "time".
#' @param index_column_position (optional) The column position of the time index for matrix inputs. Defaults to 1.
#' @param dtw_step_matrix A matrix specifying the step pattern for the DTW algorithm. Defaults to \code{\link[dtw]{symmetric1}}.
#' @param max_elapsed Numeric value constraining the TWDTW calculation to the lower band given by a maximum elapsed time. Defaults to Inf.
#' @param return_distance Logical indicating whether to return the TWDTW distance. Defaults to TRUE.
#' @param return_all_matches Logical indicating whether to find all matches within the TWDTW matrix. Defaults to FALSE.
#' @param version A string identifying the version of TWDTW implementation. Options are 'f90' for Fortran 90, 'f90goto' for Fortran 90 with goto statements, or 'cpp' for C++ version. Defaults to 'f90'.
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
#' @return
#' If return_distance = TRUE, a numeric value representing the TWDTW distance between the two time series.
#' If return_all_matches = TRUE, a matrix of all matches within the TWDTW matrix.
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
#'       time_cycle_length = 'year',
#'       time_cycle_scale = 'day',
#'       time_weight = c(steepness = 0.1, midpoint = 50))
#'
#' # Pass a generic time-weight fucntion
#' twdtw(x, y,
#'       time_cycle_length = 'year',
#'       time_cycle_scale = 'day',
#'       time_weight = function(x,y) x + 1.0 / (1.0 + exp(-0.1 * (y - 50))))
#'
#' # Test other version
#' twdtw(x, y,
#'       time_cycle_length = 'year',
#'       time_cycle_scale = 'day',
#'       time_weight = c(steepness = 0.1, midpoint = 50),
#'       version = 'f90goto')
#'
#' twdtw(x, y,
#'       time_cycle_length = 'year',
#'       time_cycle_scale = 'day',
#'       time_weight = c(steepness = 0.1, midpoint = 50),
#'       version = 'cpp')
#'
#' @include convert_date_to_numeric.R
#'
#' @export
twdtw <- function(x, y, time_weight, time_cycle_length, time_cycle_scale, ...) {
  UseMethod("twdtw")
}

#' @export
#' @rdname twdtw
twdtw.data.frame <- function(x, y,
                             time_weight,
                             time_cycle_length,
                             time_cycle_scale,
                             index_column_name = 'time',
                             dtw_step_matrix = symmetric1,
                             max_elapsed = Inf,
                             return_distance = TRUE,
                             return_all_matches = FALSE,
                             version = 'f90', ...) {

  # Check that 'x' and 'y' are data.frames
  if (!inherits(y, "data.frame")) {
    stop("Both x and y need to be data.frames")
  }

  if (is.null(time_cycle_length)) stop("The 'time_cycle_length' argument is missing.")
  if (is.null(time_cycle_scale)) stop("The 'time_cycle_scale' argument is missing.")

  # The dimensions of the time series must match
  if (!setequal(names(y), names(x))) {
    stop("The dimensions (columns) of the time series 'x' and 'y' do not match.")
  }

  # The dimensions of the time series must match and 'index_column_name' must be present in both 'x' and 'y'
  if (!setequal(names(y), names(x)) & !index_column_name %in% names(x)) {
    stop("The dimensions (columns) of the time series 'x' and 'y' do not match.")
  }

  # Position 'index_column_name' at the first column
  x <- x[, c(index_column_name, setdiff(names(x), index_column_name)), drop = FALSE]

  # Sort dimensions of y according to x
  y <- y[, names(x), drop = FALSE]

  # Convert 'index_column_name' to numeric based on 'time_cycle_length' and 'time_cycle_scale'
  x[, index_column_name] <- convert_date_to_numeric(x[, index_column_name], time_cycle_length, time_cycle_scale)
  y[, index_column_name] <- convert_date_to_numeric(y[, index_column_name], time_cycle_length, time_cycle_scale)

  # call .twdtw function
  twdtw(x = as.matrix(x),
        y = as.matrix(y),
        time_weight = time_weight,
        time_cycle_length = time_cycle_length,
        time_cycle_scale = time_cycle_scale,
        index_column_position = 1,
        dtw_step_matrix = dtw_step_matrix,
        max_elapsed = max_elapsed,
        return_distance = return_distance,
        return_all_matches = return_all_matches,
        version = version)
}

#' @export
#' @rdname twdtw
twdtw.matrix <- function(x, y,
                         time_weight,
                         time_cycle_length,
                         time_cycle_scale = NULL,
                         index_column_position = 1,
                         dtw_step_matrix = symmetric1,
                         max_elapsed = Inf,
                         return_distance = TRUE,
                         return_all_matches = FALSE,
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

  # Check if 'time_cycle_length' is declared
  if (is.null(time_cycle_length)) stop("The 'time_cycle_length' argument is missing.")

  # Get maximum possible value that a specific time cycle and scale
  if (is.character(time_cycle_length)){
    if (is.null(time_cycle_scale)) stop("The 'time_cycle_scale' argument is missing for 'time_cycle_length' type character.")
    time_cycle_length <- calculate_max_cycle_length(time_cycle_length, time_cycle_scale)
  }

  # Position 'index_column_name' at the first column
  new_order <- c(index_column_position, setdiff(1:ncol(x), index_column_position))
  x <- x[, new_order, drop = FALSE]
  y <- y[, new_order, drop = FALSE]

  # call .twdtw function
  .twdtw(x = x,
         y = y,
         time_weight = time_weight,
         dtw_step_matrix = dtw_step_matrix,
         time_cycle_length = time_cycle_length,
         max_elapsed = max_elapsed,
         return_distance = return_distance,
         return_all_matches = return_all_matches,
         version = version)
}


# Internal function
.twdtw <- function(x, y,
                   time_weight,
                   time_cycle_length,
                   dtw_step_matrix = symmetric1,
                   max_elapsed = Inf,
                   return_distance = TRUE,
                   return_all_matches = FALSE,
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
  CL <- as.double(time_cycle_length)
  LB <- as.double(max_elapsed)
  if (is.function(time_weight)){
    TW <- as.double(c(0.0, 0.0))
    time_weight <- time_weight_wrapper(time_weight)
  } else {
    TW <- as.double(time_weight)
  }

  # Get the function version using its name
  fn <- get(c('twdtw_f90',
              'twdtw_f90gt',
              'twdtw_cpp')[version == c('f90', 'f90goto', 'cpp')])

  # Prepare arguments
  args <- list(XM, YM, CM, DM, VM, SM, N, M, D, NS, TW, LB, JB, CL)
  if (version == 'f90'){
    args$tw_r <- time_weight
  }

  # Call twdtw
  do.call(fn, args)

  b <- JB[JB!=0]
  a <- VM[-1,][N,b]
  d <- CM[-1,][N,b]
  candidates <- matrix(c(a, b, d), ncol = 3, byrow = F)

  if (return_all_matches) {
    return(candidates)
  }

  return(min(d))

}

# Define the time-weight wrapper function
time_weight_wrapper <- function(fn) {
  function(x, y, z = 0.0, w = 0.0) {
    fn(x, y)
  }
}

