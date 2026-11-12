#' @title Find the k-nearest neighbours of time series using TWDTW
#'
#' @description
#' This function calculates the Time-Weighted Dynamic Time Warping (TWDTW) distance between
#' each time series in \code{x} and each temporal pattern in \code{y}, and returns for each time series
#' the \code{k} patterns with the lowest distance.
#'
#' @param x A numeric array with three dimensions: series, time, and band, i.e. \code{x[i, , ]} has the
#' observations of the i-th time series with one row for each date and one column for each band.
#' All time series share the same dates, missing observations are \code{NA}.
#' @param y A list of temporal patterns (references). Each pattern is either a data.frame or a matrix with
#' one row for each observation, a column with the time index, and one column for each band in \code{x}.
#' The patterns can have different number of observations.
#' @param time A vector with the dates of the time dimension of \code{x}. Its length must be equal to \code{dim(x)[2]}.
#' @param time_weight A numeric vector with length two (steepness and midpoint of logistic weight) or a function.
#' See details in \code{\link{twdtw}}. A function is much slower than the logistic weight, see details.
#' @param k The number of nearest neighbours. Defaults to 1, i.e. only the nearest pattern.
#' @param origin For numeric cycle_length, the origin must be specified. This is the point
#' from which the elapsed time is computed. Must be of the same class as \code{time}.
#' @param index_column The column name or the column position of the time index in the patterns in \code{y}.
#' Defaults to "time".
#' @inheritParams twdtw
#'
#' @details The distances are the same as the distances calculated with \code{\link{twdtw}} for each pair
#' of time series and pattern. However, all time series are processed in a single call to the 'Fortran' routines
#' and only the distances are calculated. This is much faster than calling \code{\link{twdtw}} for each pair
#' and suitable for a large number of time series, e.g. the pixels of satellite images.
#' The matches and the internal matrices are not available, use \code{\link{twdtw}} for that.
#'
#' The band names of \code{x}, i.e. \code{dimnames(x)[[3]]}, are matched to the column names of the patterns in \code{y}.
#' If either of them has no names, the bands are matched by their position.
#'
#' Time indices of class \code{Date} or \code{POSIXct}, or character strings, are converted to
#' the units of the cycle using \code{\link{date_to_numeric_cycle}}. Numeric time indices are
#' assumed to be already in the units of the cycle and are not converted.
#' This applies to both \code{time} and the time index of each pattern in \code{y}.
#'
#' An observation with a missing value in any band is removed from the respective time series or pattern.
#' A time series without any complete observation has \code{NA} for all neighbours. Patterns must
#' have at least one complete observation.
#'
#' Patterns with the same distance to a time series are returned in the order of their position in \code{y}.
#'
#' If \code{time_weight} is a function, \code{\link{twdtw}} is called for each pair of time series
#' and pattern, which is much slower than using the logistic time weight.
#'
#' @return A list with two matrices with one row for each time series in \code{x} and \code{k} columns
#' ordered from the nearest to the farthest neighbour:
#' \item{index}{The position in \code{y} of the nearest patterns.}
#' \item{distance}{The TWDTW distance to the nearest patterns.}
#'
#' @references
#' Maus, V., Camara, G., Cartaxo, R., Sanchez, A., Ramos, F. M., & de Moura, Y. M. (2016).
#' A Time-Weighted Dynamic Time Warping Method for Land-Use and Land-Cover Mapping.
#' IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 9(8), 3729-3739.
#' \doi{10.1109/JSTARS.2016.2517118}
#'
#' @seealso \code{\link{twdtw}}
#'
#' @examples
#'
#' # Satellite image time series of 999 pixels and temporal patterns of 5 land use classes
#' ts_x <- readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))
#' ts_y <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))
#'
#' # Stack the time series into an array with dimensions: series, time, and band
#' bands <- c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")
#' x <- sapply(ts_x$observations, function(ts) as.matrix(ts[bands]), simplify = "array")
#' x <- aperm(x, c(3, 1, 2))
#' dim(x)
#'
#' # Find the nearest pattern for each time series
#' knn <- twdtw_knn(x = x,
#'                  y = ts_y$observations,
#'                  time = ts_x$observations[[1]]$time,
#'                  cycle_length = 'year',
#'                  time_scale = 'day',
#'                  time_weight = c(steepness = 0.1, midpoint = 50))
#'
#' table(ts_y$label[knn$index[, 1]])
#'
#' # Find the two nearest patterns and their distances
#' knn <- twdtw_knn(x = x,
#'                  y = ts_y$observations,
#'                  time = ts_x$observations[[1]]$time,
#'                  cycle_length = 'year',
#'                  time_scale = 'day',
#'                  time_weight = c(steepness = 0.1, midpoint = 50),
#'                  k = 2)
#'
#' head(knn$index)
#' head(knn$distance)
#'
#' @export
twdtw_knn <- function(x, y,
                      time,
                      time_weight,
                      cycle_length,
                      time_scale = NULL,
                      k = 1,
                      origin = NULL,
                      index_column = 'time',
                      max_elapsed = Inf) {

  # Check that 'x' is an array with the time series
  if (!is.numeric(x) || length(dim(x)) != 3) {
    stop("'x' must be a numeric array with three dimensions: series, time, and band")
  }

  if (length(time) != dim(x)[2]) {
    stop("The length of 'time' must be equal to the second dimension of 'x'")
  }

  if (!(is.function(time_weight) || (is.numeric(time_weight) && length(time_weight) == 2))) {
    stop("'time_weight' should be either a function or a numeric vector with length two")
  }

  # Check if 'cycle_length' is declared
  if (is.null(cycle_length)) stop("The 'cycle_length' argument is missing.")

  if (!is.numeric(max_elapsed) || length(max_elapsed) != 1 || is.na(max_elapsed)) {
    stop("'max_elapsed' must be a single numeric value")
  }

  # A single pattern does not need to be in a list
  if (is.data.frame(y) || is.matrix(y)) {
    y <- list(y)
  }

  if (!is.list(y) || length(y) < 1) {
    stop("'y' must be a list of data.frames or matrices")
  }

  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != round(k) || k < 1 || k > length(y)) {
    stop("'k' must be an integer between one and the number of patterns in 'y'")
  }

  # Convert time index to numeric based on 'cycle_length' and 'time_scale'
  # Numeric time indices are already in the units of the cycle
  to_cycle <- function(index) {
    if (is.numeric(index)) {
      return(as.double(index))
    }
    if (is.null(time_scale)) stop("The 'time_scale' argument is missing.")
    as.double(date_to_numeric_cycle(index, cycle_length, time_scale, origin))
  }

  tx <- to_cycle(time)

  if (anyNA(tx)) {
    stop("'time' cannot have missing values")
  }

  # Get time index, bands, and complete observations of each pattern
  y <- lapply(seq_along(y), function(i) {
    .get_pattern(y[[i]], i, index_column, dimnames(x)[[3]], dim(x)[3], to_cycle)
  })

  # Get maximum possible value that a specific time cycle and scale
  if (is.character(cycle_length)) {
    if (is.null(time_scale)) stop("The 'time_scale' argument is missing for 'cycle_length' type character.")
    cycle_length <- max_cycle_length(cycle_length, time_scale)
  }

  if (is.function(time_weight)) {

    # The batch routine does not call R functions
    out <- .twdtw_knn_loop(x = x,
                           y = y,
                           tx = tx,
                           time_weight = time_weight,
                           cycle_length = cycle_length,
                           max_elapsed = max_elapsed,
                           k = k)

  } else {

    # Stack patterns by rows and get the row offset of each pattern
    YM <- do.call(rbind, lapply(y, function(p) p$values))
    YM <- matrix(as.double(YM), nrow(YM), ncol(YM))
    TY <- as.double(unlist(lapply(y, function(p) p$time)))
    YO <- as.integer(c(0, cumsum(sapply(y, function(p) nrow(p$values)))))

    # Avoids a copy if the time series are already double
    if (!is.double(x)) {
      storage.mode(x) <- "double"
    }

    # Call twdtw for all time series and patterns
    out <- twdtw_knn_f90(XM = x,
                         XD = dim(x),
                         TX = tx,
                         YM = YM,
                         TY = TY,
                         YO = YO,
                         TW = as.double(time_weight),
                         LB = as.double(max_elapsed),
                         CL = as.double(cycle_length),
                         K = as.integer(k))

  }

  rownames(out$index) <- dimnames(x)[[1]]
  rownames(out$distance) <- dimnames(x)[[1]]

  out

}

# Internal function
# Get time index in units of the cycle and the observations of a pattern as a
# matrix with the bands in the same order as the time series
.get_pattern <- function(pattern, i, index_column, bands, n_bands, to_cycle) {

  if (!(is.data.frame(pattern) || is.matrix(pattern))) {
    stop(sprintf("Pattern %d in 'y' is neither a data.frame nor a matrix", i))
  }

  # Find 'index_column' by its name or position
  if (is.character(index_column)) {
    index <- match(index_column, colnames(pattern))
  } else {
    index <- index_column
  }

  if (length(index) != 1 || is.na(index) || index < 1 || index > ncol(pattern)) {
    stop(sprintf("The 'index_column' was not found in pattern %d in 'y'", i))
  }

  time <- to_cycle(if (is.data.frame(pattern)) pattern[[index]] else pattern[, index])
  values <- as.matrix(pattern[, -index, drop = FALSE])

  if (!is.numeric(values)) {
    stop(sprintf("Pattern %d in 'y' has non-numeric columns", i))
  }

  # The dimensions of the patterns must match the time series
  if (!is.null(bands) && !is.null(colnames(values))) {
    if (!setequal(colnames(values), bands)) {
      stop(sprintf("The dimensions (columns) of pattern %d in 'y' do not match the bands of 'x'", i))
    }
    # Sort dimensions of the pattern according to x
    values <- values[, bands, drop = FALSE]
  } else if (ncol(values) != n_bands) {
    stop(sprintf("The dimensions (columns) of pattern %d in 'y' do not match the bands of 'x'", i))
  }

  # Patterns should not contain NAs
  keep <- complete.cases(time, values)

  if (!any(keep)) {
    stop(sprintf("Pattern %d in 'y' has no complete observations", i))
  }

  list(time = time[keep], values = values[keep, , drop = FALSE])

}

# Internal function
# Find the k-nearest neighbours calling twdtw for each time series and pattern
.twdtw_knn_loop <- function(x, y, tx, time_weight, cycle_length, max_elapsed, k) {

  index <- matrix(NA_integer_, nrow = dim(x)[1], ncol = k)
  distance <- matrix(NA_real_, nrow = dim(x)[1], ncol = k)

  # Position time index at the first column
  y <- lapply(y, function(p) cbind(p$time, p$values))

  for (i in seq_len(dim(x)[1])) {

    ts <- cbind(tx, matrix(x[i, , ], nrow = length(tx)))

    # Time series should not contain NAs
    ts <- ts[complete.cases(ts), , drop = FALSE]

    # Time series without observations have NA
    if (nrow(ts) < 1) next

    d <- sapply(y, function(p) {
      .twdtw(x = ts,
             y = p,
             time_weight = time_weight,
             cycle_length = cycle_length,
             max_elapsed = max_elapsed)[1]
    })

    # Sorting is stable, i.e. ties are in the order of the patterns
    nearest <- order(d)[seq_len(k)]
    index[i, ] <- nearest
    distance[i, ] <- d[nearest]

  }

  list(index = index, distance = distance)

}
