#' @title Calculate Time-Weighted Dynamic Time Warping (TWDTW) distance
#'
#' @description
#' This function calculates the Time-Weighted Dynamic Time Warping (TWDTW) distance between two time series.
#'
#' @param x A data.frame or matrix representing the first time series.
#' @param y A data.frame or matrix representing the second time series.
#' @param time_weight_par (optional) A numeric vector of parameters to be used with the time weighting function.
#' @param time_weight_fun (optional) A user-defined function for time weighting.
#' @param time_cycle_length Required for data.frame inputs. A character string indicating the larger unit of time.
#' It must be one of "year", "month", "day", "hour", "minute". It can also receive a numeric value when time_cycle_scale is numeric.
#' @param time_cycle_scale Required for data.frame inputs. A character string indicating the smaller unit of time,
#' which is a division of the time_cycle_length. If time_cycle_length is "year", time_cycle_scale can be one of
#' "month", "day", "hour", "minute", "second". If time_cycle_length is "month", time_cycle_scale can be "day",
#' "hour", "minute", "second", and so on. It can also receive a numeric value when time_cycle_length is numeric.
#' @param time_cycle_scale Required for data.frame inputs. A character string specifying the time cycle scale or an integer value.
#' @param index_column_name (optional) The column name of the time index for data.frame inputs. Defaults to "time".
#' @param index_column_position (optional) The column position of the time index for matrix inputs. Defaults to 1.
#' @param dtw_step_matrix A matrix specifying the step pattern for the TWDTW algorithm. Defaults to symmetric1.
#' @param max_elapsed Numeric value constraining the TWDTW calculation to the lower band given by a maximum elapsed time. Defaults to Inf.
#' @param return_distance Logical indicating whether to return the TWDTW distance. Defaults to TRUE.
#' @param return_all_matches Logical indicating whether to find all matches within the TWDTW matrix. Defaults to FALSE.
#' @param version A string identifying the version of TWDTW implementation. Options are 'f90' for Fortran 90, 'f90goto' for Fortran 90 with goto statements, or 'cpp' for C++ version. Defaults to 'f90'.
#' @param ... Additional parameters to pass to the time weighting function.
#'
#' @return
#' If return_distance = TRUE, a numeric value representing the TWDTW distance between the two time series.
#' If return_all_matches = TRUE, a matrix of all matches within the TWDTW matrix.
#'
#' @examples
#'
#' n <- 23
#' t <- seq(0, pi, length.out = n)
#' d <- seq(as.Date('2020-09-01'), length.out = n, by = "15 day")
#'
#' x <- data.frame(time = d, value = sin(t)*2 + runif(n))
#'
#' y <- data.frame(time = d + 30, value = sin(t)*2 + runif(n))
#'
#' plot(x, type = "l", xlim = range(c(d, d + 5)))
#' lines(y, col = "red")
#'
#' # Calculate TWDTW distance between x and y
#' twdtw(x, y,
#'       time_cycle_length = 'year',
#'       time_cycle_scale = 'day',
#'       time_weight_par = c(steepness = -0.1, midpoint = 50))
#'
#' #twdtw(x, y, max_elapsed = 20,
#' #   tw_r = function(dist,td,tw1,tw2) dist + 1.0 / (1.0 + exp(-0.1 * (td - 50))))
#'
#' # twdtw(x, y, tw = c(-.1, 50), max_elapsed = 50, version = 'f90goto')
#'
#' # twdtw(x, y, tw = c(-.1, 50), max_elapsed = 50, version = 'cpp')
#'
#' @include convert_date_to_numeric.R
#'
#' @export
twdtw <- function(x, y, ...) {
  UseMethod("twdtw")
}

#' @export
#' @rdname twdtw
twdtw.data.frame <- function(x, y,
                             time_weight_par = NULL,
                             time_weight_fun = NULL,
                             time_cycle_length = NULL,
                             time_cycle_scale = NULL,
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
        time_weight_fun = time_weight_fun,
        time_weight_par = time_weight_par,
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
                         time_weight_par = NULL,
                         time_weight_fun = NULL,
                         time_cycle_length = NULL,
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

  # Check if 'time_weight_fun' is a function or NULL
  if(!is.null(time_weight_fun) & !is.function(time_weight_fun)) stop("'time_weight_fun' should be either a function or NULL")

  # Check if 'time_weight_par' is a numeric vector or NULL
  if(!is.null(time_weight_par) & !(is.numeric(time_weight_par))) stop("'time_weight_par' should be either a numeric vector of length 2 or NULL")

  # If both 'time_weight_fun' and 'time_weight_par' are NULL, stop the function and return an error message
  if(is.null(time_weight_fun) & is.null(time_weight_par)) stop("Either 'time_weight_fun' or 'time_weight_par' or both should be provided")

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
         time_weight_fun = time_weight_fun,
         time_weight_par = time_weight_par,
         dtw_step_matrix = dtw_step_matrix,
         time_cycle_length = time_cycle_length,
         max_elapsed = max_elapsed,
         return_distance = return_distance,
         return_all_matches = return_all_matches,
         version = version)
}


# Internal function
.twdtw <- function(x, y,
                   time_weight_par = NULL,
                   time_cycle_length = NULL,
                   time_weight_fun = NULL,
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
  TW <- as.double(time_weight_par)
  CL <- as.double(time_cycle_length)
  LB <- as.double(max_elapsed)

  # Get the function version using its name
  fn <- get(c('twdtw_f90',
              'twdtw_f90gt',
              'twdtw_cpp')[version == c('f90', 'f90goto', 'cpp')])

  # Prepare arguments
  args <- list(XM, YM, CM, DM, VM, SM, N, M, D, NS, TW, LB, JB, CL)
  if (version == 'f90') args$tw_r <- time_weight_fun

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

