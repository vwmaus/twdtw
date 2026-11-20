#' @title Create a k-nearest neighbour model using TWDTW
#'
#' @description
#' This function creates a k-nearest neighbour (KNN) model to classify time series using the
#' Time-Weighted Dynamic Time Warping (TWDTW) distance. The model keeps labelled time series as
#' temporal patterns. \code{\link[=predict.twdtw_knn]{predict}} classifies new time series with
#' the labels of their nearest temporal patterns.
#'
#' @param x The labelled time series (samples) in one of three forms: a list of data.frames or matrices,
#' a data.frame with one row for each time series, or a numeric array with the dimensions
#' series, time, and band. See details.
#' @param y A vector with the label of each time series in \code{x}, in the same order.
#' @param cycle_length The length of the cycle, e.g. \code{'year'}. See details in \code{\link{twdtw}}.
#' @param time_scale The time scale of the observations, e.g. \code{'day'}. See details in \code{\link{twdtw}}.
#' It is not needed if both the time indices and \code{cycle_length} are numeric.
#' @param resampling_fun A function that creates one temporal pattern for each label from all time series
#' with this label. If \code{NULL} (default), the model keeps all time series as temporal patterns. See details.
#' @param resampling_freq The time between the observations of the temporal patterns created with
#' \code{resampling_fun}, e.g. \code{16} for 16 days if the time indices are dates. If \code{NULL} (default),
#' the temporal patterns have the time indices of the observations of the time series.
#' @param keep_models Whether the model keeps the fitted models that \code{resampling_fun} returns for
#' each label and band, e.g. the \code{gam} objects. If \code{FALSE} (default), it keeps only their root
#' mean squared error. The temporal patterns and the predictions do not depend on it. See details.
#' @param start,end The time interval in which the label of each time series is valid. Either one value
#' for all time series or a vector with one value for each time series, of the same class as the time indices.
#' Only the observations from \code{start} to \code{end}, including both, are used.
#' If \code{NULL} (default), all observations are used.
#' @param index_column The column name or the column position of the time index of the time series in a list.
#' Defaults to "time".
#' @param time A vector with the dates of the time dimension of an array \code{x}.
#' Its length must be equal to \code{dim(x)[2]}. It is only used if \code{x} is an array.
#' @param ... Arguments of the methods. An argument that a method does not know is an error.
#' @inheritParams twdtw
#'
#' @details The time series in \code{x} can have three forms:
#'
#' \itemize{
#'   \item A list of data.frames or matrices. Each of them is a time series with one row for each observation,
#'   a column with the time index, and one column for each band, as the time series of \code{\link{twdtw}}.
#'   The time series can have different dates and different number of observations.
#'   \item A data.frame with one row for each time series and one column for each band and date,
#'   named with the band and the date separated by a period, e.g. \code{EVI.2011-09-14}. Names with the date
#'   first, e.g. \code{2011-09-14.EVI}, and dates written with periods, e.g. \code{EVI.2011.09.14}, are also
#'   recognised. Columns with only the date in their name are the observations of a single band.
#'   Columns with other names are not used.
#'   \item A numeric array with three dimensions: series, time, and band, i.e. \code{x[i, , ]} has the
#'   observations of the i-th time series with one row for each date and one column for each band.
#'   The dates of the time dimension are given in \code{time}.
#' }
#'
#' All time series must have the same bands. Observations with dates are used in time order.
#' Numeric time indices are in the units of the cycle, as in \code{\link{twdtw}}, and their
#' observations are used in the order in which they are given.
#' An observation with a missing value in any band is not used to calculate the TWDTW distance,
#' and time series without any complete observation are removed with a warning.
#'
#' If \code{resampling_fun} is \code{NULL}, the model keeps all time series as temporal patterns.
#' Such a model can have several patterns for each label and can be classified with more than one
#' nearest neighbour, see the argument \code{k} of \code{\link{predict.twdtw_knn}}.
#'
#' If \code{resampling_fun} is a function, the model has one temporal pattern for each label. This
#' reduces the number of TWDTW distances that \code{predict} calculates from the number of labelled
#' time series to the number of labels. The function is called for each label and band with a
#' data.frame that has two columns: \code{x}, the time index as a number, and \code{y}, the values
#' of all time series with this label. It must return a model for which
#' \code{predict(model, newdata = data.frame(x = ...))} returns one numeric value for each row of
#' \code{newdata}, e.g. \code{function(data) lm(y ~ factor(x), data = data)} for the mean of each date, or
#' \code{function(data) mgcv::gam(y ~ s(x), data = data)} for a smoothed pattern. Before the function is
#' called, time series with dates are shifted to start on the same day of the year in the year 1970,
#' so that time series of different years overlap. This aligns the time series by the year. For another
#' cycle, e.g. \code{cycle_length = 'day'}, use numeric time indices in the units of the cycle, which are
#' not shifted.
#'
#' Most of the time to create the temporal patterns is used to fit these models. With thousands of
#' time series, \code{function(data) mgcv::bam(y ~ s(x), data = data, discrete = TRUE)} gives nearly
#' the same patterns as \code{mgcv::gam} and is about two times faster.
#'
#' The models that \code{resampling_fun} returns, one for each label and band, are only used to create
#' the temporal patterns: \code{predict} uses the patterns and not these models. By default the model
#' keeps the patterns and the root mean squared error of each fitted model, which \code{print} shows,
#' and not the fitted models. A fitted model usually contains the data it was fitted to, i.e. the
#' observations of all time series of a label. Without the fitted models the model is much smaller,
#' e.g. to save it to a file or to send it to other processes. With \code{keep_models = TRUE} the fitted
#' models are kept in \code{smooth_models}, e.g. to plot them or to see their summary:
#' \code{summary(m$smooth_models[[label]][[band]])}.
#'
#' The temporal patterns of different labels can have a different number of observations, e.g. if the
#' time series of the labels are from different years and their dates are one day apart after the shift.
#' The TWDTW distance is not divided by the number of observations of the pattern, so a pattern with fewer
#' observations tends to be nearer. With \code{resampling_freq} all patterns have observations at the
#' same interval.
#'
#' Creating the model does not calculate TWDTW distances. The distances between new time series and
#' the temporal patterns are calculated by \code{\link{predict.twdtw_knn}} in a single call to
#' \code{\link{twdtw_nearest}}.
#'
#' \code{twdtw_knn} is a generic function. Other packages can add methods for their classes of
#' time series, e.g. for the pixels of a satellite image data cube.
#'
#' @return An object of class \code{twdtw_knn}, which is a list with:
#' \item{call}{The call of the function.}
#' \item{smooth_models}{The models returned by \code{resampling_fun}, if \code{keep_models} is
#' \code{TRUE}: a list with one element for each label, which is a list with the model of each band.
#' Otherwise \code{NULL}.}
#' \item{rmse}{The root mean squared error of the models returned by \code{resampling_fun}: a matrix with
#' one row for each band and one column for each label, or \code{NULL} if the model keeps all time series.}
#' \item{data}{A data.frame with one row for each temporal pattern and the columns \code{label} and
#' \code{observations}. The latter is a list of data.frames with the column \code{time} and one column for each band.}
#' \item{twdtw_args}{A list with the arguments of the TWDTW distance: \code{time_weight}, \code{cycle_length},
#' \code{time_scale}, \code{origin}, and \code{max_elapsed}.}
#'
#' @references
#' Maus, V., Camara, G., Cartaxo, R., Sanchez, A., Ramos, F. M., & de Queiroz, G. R. (2016).
#' A Time-Weighted Dynamic Time Warping Method for Land-Use and Land-Cover Mapping.
#' IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 9(8), 3729-3739.
#' \doi{10.1109/JSTARS.2016.2517118}
#'
#' @seealso \code{\link{predict.twdtw_knn}}, \code{\link{twdtw_nearest}}, \code{\link{twdtw}}
#'
#' @examples
#'
#' # Temporal patterns of 5 land use classes and satellite image time series of 999 pixels
#' ts_y <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))
#' ts_x <- readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))
#'
#' # Create a model from a list of labelled time series
#' m <- twdtw_knn(x = ts_y$observations,
#'                y = ts_y$label,
#'                cycle_length = 'year',
#'                time_scale = 'day',
#'                time_weight = c(steepness = 0.1, midpoint = 50))
#'
#' m
#'
#' # Classify the time series of the pixels
#' lu <- predict(m, newdata = ts_x$observations)
#'
#' head(lu)
#'
#' table(lu$prediction)
#'
#' # Use the first 300 classified time series as samples and create a model
#' # with one temporal pattern for each label: the mean of the samples on each date
#' samples <- ts_x$observations[1:300]
#' labels <- lu$prediction[1:300]
#'
#' m <- twdtw_knn(x = samples,
#'                y = labels,
#'                resampling_fun = function(data) lm(y ~ factor(x), data = data),
#'                cycle_length = 'year',
#'                time_scale = 'day',
#'                time_weight = c(steepness = 0.1, midpoint = 50))
#'
#' m
#'
#' # A model without 'resampling_fun' keeps all samples as temporal patterns
#' # and can be classified with more than one nearest neighbour
#' m <- twdtw_knn(x = samples,
#'                y = labels,
#'                cycle_length = 'year',
#'                time_scale = 'day',
#'                time_weight = c(steepness = 0.1, midpoint = 50))
#'
#' lu <- predict(m, newdata = ts_x$observations[301:310], k = 5)
#'
#' lu
#'
#' @export
twdtw_knn <- function(x, y, time_weight, cycle_length, time_scale, ..., time = NULL) {
  UseMethod("twdtw_knn")
}

#' @export
#' @rdname twdtw_knn
twdtw_knn.list <- function(x, y,
                           time_weight,
                           cycle_length,
                           time_scale,
                           resampling_fun = NULL,
                           resampling_freq = NULL,
                           start = NULL,
                           end = NULL,
                           origin = NULL,
                           max_elapsed = Inf,
                           index_column = 'time', ...,
                           keep_models = FALSE) {

  .check_unknown_args(list(...), "twdtw_knn", "x")

  if (missing(time_scale)) {
    time_scale <- NULL
  }

  ts <- .time_series_list(x, index_column, "x")

  .twdtw_knn_model(ts = ts,
                   y = y,
                   time_weight = time_weight,
                   cycle_length = cycle_length,
                   time_scale = time_scale,
                   resampling_fun = resampling_fun,
                   resampling_freq = resampling_freq,
                   keep_models = keep_models,
                   start = start,
                   end = end,
                   origin = origin,
                   max_elapsed = max_elapsed,
                   call = .generic_call(match.call(), "twdtw_knn"))

}

#' @export
#' @rdname twdtw_knn
twdtw_knn.data.frame <- function(x, y,
                                 time_weight,
                                 cycle_length,
                                 time_scale,
                                 resampling_fun = NULL,
                                 resampling_freq = NULL,
                                 start = NULL,
                                 end = NULL,
                                 origin = NULL,
                                 max_elapsed = Inf, ...,
                                 keep_models = FALSE) {

  .check_unknown_args(list(...), "twdtw_knn", "x")

  if (missing(time_scale)) {
    time_scale <- NULL
  }

  # Columns with only the date in their name are the observations of a single band
  x <- .wide_table_array(x, single_band = "value", arg = "x")
  ts <- .array_time_series(x$values, x$time, "x")

  .twdtw_knn_model(ts = ts,
                   y = y,
                   time_weight = time_weight,
                   cycle_length = cycle_length,
                   time_scale = time_scale,
                   resampling_fun = resampling_fun,
                   resampling_freq = resampling_freq,
                   keep_models = keep_models,
                   start = start,
                   end = end,
                   origin = origin,
                   max_elapsed = max_elapsed,
                   call = .generic_call(match.call(), "twdtw_knn"))

}

#' @export
#' @rdname twdtw_knn
twdtw_knn.array <- function(x, y,
                            time_weight,
                            cycle_length,
                            time_scale,
                            resampling_fun = NULL,
                            resampling_freq = NULL,
                            start = NULL,
                            end = NULL,
                            origin = NULL,
                            max_elapsed = Inf,
                            time = NULL, ...,
                            keep_models = FALSE) {

  .check_unknown_args(list(...), "twdtw_knn", "x")

  if (missing(time_scale)) {
    time_scale <- NULL
  }

  x <- .check_time_series_array(x, time, "x")
  ts <- .array_time_series(x$values, x$time, "x")

  .twdtw_knn_model(ts = ts,
                   y = y,
                   time_weight = time_weight,
                   cycle_length = cycle_length,
                   time_scale = time_scale,
                   resampling_fun = resampling_fun,
                   resampling_freq = resampling_freq,
                   keep_models = keep_models,
                   start = start,
                   end = end,
                   origin = origin,
                   max_elapsed = max_elapsed,
                   call = .generic_call(match.call(), "twdtw_knn"))

}

#' @export
twdtw_knn.default <- function(x, y, time_weight, cycle_length, time_scale, ...) {

  stop("No method of twdtw_knn for 'x' of class '", class(x)[1], "'. ",
       "'x' must be a list of data.frames or matrices, a data.frame with one row for each time series, ",
       "or a numeric array with the dimensions series, time, and band")

}

# Internal function
# Create a model of class 'twdtw_knn' from a list of time series, which are
# data.frames with the column 'time' followed by one column for each band
.twdtw_knn_model <- function(ts, y,
                             time_weight,
                             cycle_length,
                             time_scale,
                             resampling_fun,
                             resampling_freq,
                             keep_models,
                             start,
                             end,
                             origin,
                             max_elapsed,
                             call) {

  # One label for each time series
  if (is.factor(y)) {
    y <- as.character(y)
  }

  if (!is.atomic(y) || is.null(y) || length(y) != length(ts)) {
    stop(sprintf("'y' must be a vector with one label for each time series in 'x' (%d)", length(ts)))
  }

  if (anyNA(y)) {
    stop("'y' cannot have missing values")
  }

  labels <- as.character(y)

  if (!is.null(resampling_fun) && !is.function(resampling_fun)) {
    stop("'resampling_fun' must be a function")
  }

  if (!is.logical(keep_models) || length(keep_models) != 1 || is.na(keep_models)) {
    stop("'keep_models' must be TRUE or FALSE")
  }

  # All time series have the same kind of time index
  time_class <- vapply(ts, function(s) .time_class(s$time), character(1))
  if (length(unique(time_class)) > 1) {
    stop("The time indices of the time series in 'x' must be all numeric, all dates, or all date-times")
  }

  .check_twdtw_args(time_weight, cycle_length, time_scale, origin, max_elapsed, ts[[1]]$time)

  # Keep only the observations in the time interval of the label
  ts <- .cut_time_series(ts, start, end)

  # A time series without complete observations cannot be compared
  has_observations <- vapply(ts, .has_complete_observation, logical(1))
  interval <- if (is.null(start) && is.null(end)) "" else " between 'start' and 'end'"

  if (!any(has_observations)) {
    stop("No time series in 'x' has a complete observation", interval)
  }

  if (!all(has_observations)) {
    warning("Time series without complete observations", interval, " were removed: ",
            sum(!has_observations), " of ", length(has_observations), call. = FALSE)
    ts <- ts[has_observations]
    labels <- labels[has_observations]
  }

  # Create one temporal pattern for each label
  smooth_models <- NULL
  rmse <- NULL
  if (!is.null(resampling_fun)) {
    patterns <- .resample_patterns(ts, labels, resampling_fun, resampling_freq, keep_models)
    ts <- patterns$patterns
    labels <- patterns$labels
    smooth_models <- patterns$smooth_models
    rmse <- patterns$rmse

    # A pattern without complete observations cannot be compared
    has_observations <- vapply(ts, .has_complete_observation, logical(1))
    if (!all(has_observations)) {
      stop("The temporal pattern created with 'resampling_fun' for the label '", labels[!has_observations][1],
           "' has no complete observation")
    }
  }

  data <- data.frame(label = labels, stringsAsFactors = FALSE)
  data$observations <- ts

  model <- list(call = call,
                smooth_models = smooth_models,
                rmse = rmse,
                data = data,
                twdtw_args = list(time_weight = time_weight,
                                  cycle_length = cycle_length,
                                  time_scale = time_scale,
                                  origin = origin,
                                  max_elapsed = max_elapsed))

  # A model with a single class can be used by predict functions of other
  # packages that look for a predict method by the class of the model
  class(model) <- "twdtw_knn"

  model

}

# Internal function
# The call of a method as a call of its generic function
.generic_call <- function(call, generic) {

  call[[1L]] <- as.name(generic)

  call

}

# Internal function
# Arguments that a function does not know are an error, so that they are not
# ignored without a message, e.g. an argument with a spelling mistake
.check_unknown_args <- function(args, fun, data_arg) {

  if (length(args) == 0) {
    return(invisible(NULL))
  }

  arg_names <- names(args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(args))
  }

  msg <- paste0("Unknown arguments in ", fun, ": ",
                paste(ifelse(arg_names == "", "an argument without name", paste0("'", arg_names, "'")),
                      collapse = ", "), ".")

  if ("time" %in% arg_names) {
    msg <- paste0(msg, " 'time' is only used if '", data_arg, "' is an array.")
  }

  stop(msg, call. = FALSE)

}

# Internal function
# Check the arguments of the TWDTW distance of a model with the messages of
# twdtw and twdtw_nearest, so that a model that was created can be used
.check_twdtw_args <- function(time_weight, cycle_length, time_scale, origin, max_elapsed, time) {

  if (!(is.function(time_weight) || (is.numeric(time_weight) && length(time_weight) == 2))) {
    stop("'time_weight' should be either a function or a numeric vector with length two")
  }

  if (is.null(cycle_length)) stop("The 'cycle_length' argument is missing.")

  if (!is.numeric(max_elapsed) || length(max_elapsed) != 1 || is.na(max_elapsed)) {
    stop("'max_elapsed' must be a single numeric value")
  }

  if (is.character(cycle_length)) {
    if (is.null(time_scale)) stop("The 'time_scale' argument is missing for 'cycle_length' type character.")
    max_cycle_length(cycle_length, time_scale)
  }

  # Dates are converted to the units of the cycle
  if (!is.numeric(time)) {
    if (is.null(time_scale)) stop("The 'time_scale' argument is missing.")
    date_to_numeric_cycle(time, cycle_length, time_scale, origin)
  }

  invisible(NULL)

}

#' Print method for twdtw_knn class
#'
#' @param x An object of class `twdtw_knn`, created by \code{\link{twdtw_knn}}.
#' @param ... Not used.
#'
#' @return The object `x`, invisibly. The function prints the call that created the model,
#' the number of temporal patterns and observations for each label, the root mean squared error
#' of the models returned by `resampling_fun`, and the arguments of the TWDTW distance.
#'
#' @export
print.twdtw_knn <- function(x, ...) {

  cat("\nModel of class 'twdtw_knn'\n")
  cat("----------------------------\n")

  # A call that has the data instead of the names of the objects is shortened
  cat("Call:\n")
  call_lines <- deparse(x$call, width.cutoff = 80L, nlines = 13L)
  if (length(call_lines) > 12) {
    call_lines <- c(call_lines[1:12], "...")
  }
  cat(call_lines, sep = "\n")

  patterns <- x$data$observations
  labels <- as.character(x$data$label)
  bands <- setdiff(names(patterns[[1]]), "time")

  # The root mean squared error of the models returned by 'resampling_fun'. A
  # model that was created before it was kept has the models
  rmse <- x$rmse
  if (is.null(rmse) && !is.null(x$smooth_models)) {
    rmse <- sapply(x$smooth_models, function(label_models) sapply(label_models, .model_rmse))
  }

  cat("\nData:\n")
  if (is.null(rmse)) {
    cat(" ", length(patterns), "temporal patterns, the labelled time series\n")
  } else {
    cat(" ", length(patterns), "temporal patterns, one for each label created with 'resampling_fun'\n")
  }
  cat("  Bands:", paste(bands, collapse = ", "), "\n\n")

  # Number of patterns and observations and the time interval of each label
  label_summary <- do.call(rbind, lapply(sort(unique(labels)), function(label) {
    label_patterns <- patterns[labels == label]
    n_observations <- range(vapply(label_patterns, nrow, integer(1)))
    time <- do.call(c, lapply(label_patterns, function(p) p$time))
    data.frame(label = label,
               patterns = length(label_patterns),
               observations = paste(unique(n_observations), collapse = "-"),
               from = format(min(time)),
               to = format(max(time)),
               stringsAsFactors = FALSE)
  }))
  print(label_summary, row.names = FALSE)

  if (!is.null(rmse)) {
    cat("\nRoot Mean Squared Error (RMSE) of smooth models:\n")
    print(rmse)
  }

  cat("\nTWDTW Arguments:\n")
  for (name in names(x$twdtw_args)) {
    cat(paste0(" - ", name, ": ", .format_argument(x$twdtw_args[[name]]), "\n"))
  }

  invisible(x)

}

# Internal function
# The value of an argument as a character string to print
.format_argument <- function(value) {

  if (is.null(value)) {
    return("NULL")
  }

  # The arguments of a function, e.g. 'function (dist, el)'
  if (is.function(value)) {
    fun_args <- args(value)
    return(if (is.null(fun_args)) "function" else trimws(deparse(fun_args)[1]))
  }

  if (is.atomic(value) && !is.null(names(value))) {
    return(paste0("c(", paste(names(value), value, sep = "=", collapse = ", "), ")"))
  }

  if (is.atomic(value) && length(value) > 1) {
    return(paste0("c(", paste(value, collapse = ", "), ")"))
  }

  format(value)

}
