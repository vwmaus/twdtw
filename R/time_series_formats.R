# Internal functions that convert time series between the forms accepted by
# twdtw_knn and predict.twdtw_knn: a list of data.frames or matrices, a data.frame
# with one column for each band and date, and an array with the dimensions series,
# time, and band, which is the input of twdtw_nearest


# Internal function
# Check a time index and convert character strings to dates. Numeric time
# indices are kept, they are already in the units of the cycle
.standard_time <- function(x, what) {

  if (is.factor(x)) {
    x <- as.character(x)
  }

  if (is.character(x)) {
    x <- to_date_time(x)
  }

  if (inherits(x, "POSIXlt")) {
    x <- as.POSIXct(x)
  }

  if (!(is.numeric(x) || inherits(x, c("Date", "POSIXct")))) {
    stop(sprintf("%s must be numeric, a date, a date-time, or a character string with dates", what))
  }

  x

}

# Internal function
# Values of a band are numeric. A band without any value can be logical, which
# is the type of NA
.is_band_values <- function(x) {

  is.numeric(x) || (is.logical(x) && all(is.na(x)))

}

# Internal function
# The kind of a time index, to check that two time indices can be compared
.time_class <- function(x) {

  if (inherits(x, "Date")) {
    return("Date")
  }

  if (inherits(x, "POSIXct")) {
    return("POSIXct")
  }

  "numeric"

}

# Internal function
# Dates and date-times that are not in time order. Numeric time indices are in
# the units of the cycle and can start again within a time series, e.g. the day
# of the year of a time series from September to August, so they are used in
# the order in which they are given
.is_unsorted_time <- function(time) {

  !is.numeric(time) && is.unsorted(time)

}

# Internal function
# Get the time index and the observations of a time series given as a data.frame
# or a matrix. The observations are returned as a numeric matrix with one column
# for each band, without observations that have no time index. Observations with
# dates are in time order
.split_time_series <- function(ts, i, index_column, arg) {

  if (!(is.data.frame(ts) || is.matrix(ts))) {
    stop(sprintf("Element %d of '%s' is neither a data.frame nor a matrix", i, arg))
  }

  # Find 'index_column' by its name or position
  if (is.character(index_column)) {
    index <- match(index_column, colnames(ts))
  } else {
    index <- index_column
  }

  if (length(index) != 1 || is.na(index) || index < 1 || index > ncol(ts)) {
    stop(sprintf("The 'index_column' was not found in element %d of '%s'", i, arg))
  }

  time <- .standard_time(if (is.data.frame(ts)) ts[[index]] else ts[, index],
                         sprintf("The time index of element %d of '%s'", i, arg))

  if (is.data.frame(ts)) {
    values <- ts[-index]
    is_numeric <- all(vapply(values, .is_band_values, logical(1)))
    values <- as.matrix(values)
  } else {
    values <- ts[, -index, drop = FALSE]
    is_numeric <- .is_band_values(values)
  }

  if (!is_numeric) {
    stop(sprintf("Element %d of '%s' has non-numeric columns", i, arg))
  }

  # A time series without observations is a matrix without rows
  storage.mode(values) <- "double"

  if (is.null(colnames(values))) {
    colnames(values) <- paste0("V", seq_len(ncol(values)))
  }
  rownames(values) <- NULL

  if (anyNA(time)) {
    keep <- !is.na(time)
    time <- time[keep]
    values <- values[keep, , drop = FALSE]
  }

  # Observations are used in time order
  if (.is_unsorted_time(time)) {
    time_order <- order(time)
    time <- time[time_order]
    values <- values[time_order, , drop = FALSE]
  }

  list(time = time, values = values)

}

# Internal function
# Sort the bands of a time series according to 'bands'
.sort_bands <- function(values, bands, i, arg, of) {

  if (identical(colnames(values), bands)) {
    return(values)
  }

  if (!setequal(colnames(values), bands)) {
    stop(sprintf("The bands of element %d of '%s' (%s) are not the bands of %s (%s)", i, arg,
                 paste(colnames(values), collapse = ", "), of, paste(bands, collapse = ", ")))
  }

  values[, bands, drop = FALSE]

}

# Internal function
# Convert a list of time series to a list of data.frames with the column 'time'
# followed by one column for each band, with the bands in the same order
.time_series_list <- function(x, index_column, arg) {

  if (length(x) < 1) {
    stop(sprintf("'%s' has no time series", arg))
  }

  bands <- NULL

  lapply(seq_along(x), function(i) {

    ts <- .split_time_series(x[[i]], i, index_column, arg)

    if (is.null(bands)) {
      bands <<- colnames(ts$values)
      .check_band_names(bands, arg)
    }

    values <- .sort_bands(ts$values, bands, i, arg, "the first time series")

    data.frame(time = ts$time, values, check.names = FALSE)

  })

}

# Internal function
# Convert a list of time series to arrays with the dimensions series, time, and
# band. Time series with the same time index are in the same array, so that a
# list of time series with the same dates gives a single array. 'bands' has the
# band names in the order of the third dimension. Without band names in the time
# series the bands are matched by their position
.time_series_arrays <- function(x, index_column, bands, arg) {

  if (length(x) < 1) {
    stop(sprintf("'%s' has no time series", arg))
  }

  ts <- lapply(seq_along(x), function(i) {

    ts <- .split_time_series(x[[i]], i, index_column, arg)

    if (is.null(colnames(x[[i]]))) {
      if (ncol(ts$values) != length(bands)) {
        stop(sprintf("Element %d of '%s' has %d bands and the model has %d",
                     i, arg, ncol(ts$values), length(bands)))
      }
      colnames(ts$values) <- bands
    }

    ts$values <- .sort_bands(ts$values, bands, i, arg, "the model")

    ts

  })

  # Group time series with the same time index
  first_time <- ts[[1]]$time
  if (all(vapply(ts, function(s) identical(s$time, first_time), logical(1)))) {
    groups <- list(seq_along(ts))
  } else {
    key <- vapply(ts, function(s) {
      paste(c(.time_class(s$time), as.numeric(s$time)), collapse = " ")
    }, character(1))
    groups <- unname(split(seq_along(ts), factor(key, levels = unique(key))))
  }

  lapply(groups, function(rows) {

    time <- ts[[rows[1]]]$time
    values <- vapply(ts[rows], function(s) s$values, matrix(0, length(time), length(bands)))
    values <- aperm(array(values, dim = c(length(time), length(bands), length(rows))), c(3, 1, 2))
    dimnames(values) <- list(NULL, NULL, bands)

    list(rows = rows, time = time, values = values)

  })

}

# Internal function
# Get the date and the band from the names of the columns of a data.frame with
# one column for each band and date. The names have the band and the date
# separated by a period, 'band.YYYY-MM-DD', or the date first, 'YYYY-MM-DD.band'
# and 'XYYYY-MM-DD.band'. The date can also be written with periods, 'YYYY.MM.DD',
# which is how data.frame() writes the names. Names with only the date,
# 'YYYY-MM-DD' or 'XYYYY-MM-DD', have a missing band. Other names have a missing
# date and a missing band
.parse_band_date <- function(var_names) {

  date_format <- "[0-9]{4}[.-][0-9]{2}[.-][0-9]{2}"

  name_formats <- list(
    list(pattern = paste0("^(.+)\\.(", date_format, ")$"), band = "\\1", date = "\\2"),
    list(pattern = paste0("^X?(", date_format, ")\\.(.+)$"), band = "\\2", date = "\\1"),
    list(pattern = paste0("^X?(", date_format, ")$"), band = NA, date = "\\1")
  )

  band <- rep(NA_character_, length(var_names))
  date <- rep(NA_character_, length(var_names))
  not_parsed <- rep(TRUE, length(var_names))

  for (name_format in name_formats) {
    found <- not_parsed & grepl(name_format$pattern, var_names)
    date[found] <- sub(name_format$pattern, name_format$date, var_names[found])
    if (!is.na(name_format$band)) {
      band[found] <- sub(name_format$pattern, name_format$band, var_names[found])
    }
    not_parsed <- not_parsed & !found
  }

  date <- as.Date(gsub(".", "-", date, fixed = TRUE), format = "%Y-%m-%d")
  band[is.na(date)] <- NA_character_

  list(date = date, band = band)

}

# Internal function
# Convert a data.frame with one row for each time series and one column for
# each band and date to an array with the dimensions series, time, and band.
# Columns with other names, such as the coordinates of a pixel, are not used.
# Columns with only the date in their name have the values of 'single_band', or
# are not used if 'single_band' is NULL
.wide_table_array <- function(x, single_band = NULL, arg = "x") {

  x <- as.data.frame(x)

  band_date <- .parse_band_date(names(x))

  if (!is.null(single_band)) {
    band_date$band[!is.na(band_date$date) & is.na(band_date$band)] <- single_band
  }

  is_observation <- !is.na(band_date$date) & !is.na(band_date$band)

  if (!any(is_observation)) {
    stop(sprintf(paste("'%s' has no columns named with a band and a date, e.g. 'EVI.2011-09-14'.",
                       "Time series with a column for the time index must be in a list"), arg))
  }

  not_numeric <- names(x)[is_observation][!vapply(x[is_observation], .is_band_values, logical(1))]
  if (length(not_numeric) > 0) {
    stop(sprintf("The column '%s' of '%s' is not numeric", not_numeric[1], arg))
  }

  # Observations are used in time order
  ts_time <- sort(unique(band_date$date[is_observation]))
  ts_bands <- unique(band_date$band[is_observation])
  .check_band_names(ts_bands, arg)

  values <- array(NA_real_,
                  dim = c(nrow(x), length(ts_time), length(ts_bands)),
                  dimnames = list(NULL, NULL, ts_bands))

  time_position <- match(band_date$date, ts_time)
  band_position <- match(band_date$band, ts_bands)

  for (j in which(is_observation)) {
    values[, time_position[j], band_position[j]] <- as.numeric(x[[j]])
  }

  list(values = values, time = ts_time)

}

# Internal function
# Check an array of time series and its time index, and sort the time dimension
# if its dates are not in time order
.check_time_series_array <- function(x, time, arg) {

  if (!is.numeric(x) || length(dim(x)) != 3) {
    stop(sprintf("'%s' must be a numeric array with three dimensions: series, time, and band", arg))
  }

  if (is.null(time)) {
    stop(sprintf("The 'time' argument is missing: the dates of the second dimension of '%s'", arg))
  }

  time <- .standard_time(time, "'time'")

  if (length(time) != dim(x)[2]) {
    stop(sprintf("The length of 'time' must be equal to the second dimension of '%s'", arg))
  }

  if (anyNA(time)) {
    stop("'time' cannot have missing values")
  }

  if (.is_unsorted_time(time)) {
    time_order <- order(time)
    time <- time[time_order]
    x <- x[, time_order, , drop = FALSE]
  }

  list(values = x, time = time)

}

# Internal function
# Convert an array with the dimensions series, time, and band to a list of
# data.frames with the column 'time' followed by one column for each band
.array_time_series <- function(values, time, arg) {

  bands <- dimnames(values)[[3]]
  if (is.null(bands)) {
    bands <- paste0("V", seq_len(dim(values)[3]))
  }
  .check_band_names(bands, arg)

  lapply(seq_len(dim(values)[1]), function(i) {
    ts <- matrix(values[i, , ], nrow = length(time), dimnames = list(NULL, bands))
    data.frame(time = time, ts, check.names = FALSE)
  })

}

# Internal function
# The time index of the temporal patterns of a model is in the column 'time',
# so a band cannot have this name, and each band needs its own name
.check_band_names <- function(bands, arg) {

  if ("time" %in% bands) {
    stop(sprintf("A band of '%s' is named 'time', which is the name of the time index", arg))
  }

  if (anyDuplicated(bands) > 0) {
    stop(sprintf("The bands of '%s' do not have unique names", arg))
  }

  invisible(NULL)

}

# Internal function
# Keep the observations of each time series from 'start' to 'end', including
# both. 'start' and 'end' have one value for all time series or one value for
# each time series
.cut_time_series <- function(ts, start, end) {

  if (is.null(start) && is.null(end)) {
    return(ts)
  }

  time_class <- .time_class(ts[[1]]$time)

  time_limit <- function(limit, name) {

    if (is.null(limit)) {
      return(NULL)
    }

    limit <- .standard_time(limit, sprintf("'%s'", name))

    if (!length(limit) %in% c(1, length(ts))) {
      stop(sprintf("'%s' must have one value or one value for each time series in 'x' (%d)", name, length(ts)))
    }

    if (anyNA(limit)) {
      stop(sprintf("'%s' cannot have missing values", name))
    }

    if (.time_class(limit) != time_class) {
      stop(sprintf("'%s' must have the class of the time index of the time series (%s)", name, time_class))
    }

    rep(limit, length.out = length(ts))

  }

  start <- time_limit(start, "start")
  end <- time_limit(end, "end")

  lapply(seq_along(ts), function(i) {

    keep <- rep(TRUE, nrow(ts[[i]]))

    if (!is.null(start)) {
      keep <- keep & ts[[i]]$time >= start[i]
    }

    if (!is.null(end)) {
      keep <- keep & ts[[i]]$time <= end[i]
    }

    ts[[i]][keep, , drop = FALSE]

  })

}

# Internal function
# Shift dates to start on the same day of the year in the year of 'origin',
# keeping the number of days between the observations. Time series of different
# years then overlap, and their observations can be combined in one temporal
# pattern. A time series that starts on day 366 of a leap year starts on the
# last day of the origin year if this year has 365 days. Date-times keep their
# time of the day. Numeric time indices are not shifted
.shift_dates <- function(x, origin = "1970-01-01") {

  if (is.numeric(x) || length(x) < 1) {
    return(x)
  }

  # Calendar date of the first observation, for a date-time in its own time zone
  first_date <- as.Date(format(x[1], "%Y-%m-%d"))
  day_of_year <- as.numeric(format(first_date, "%j"))

  origin_year <- format(as.Date(origin), "%Y")
  origin_days <- as.numeric(format(as.Date(paste0(origin_year, "-12-31")), "%j"))
  day_of_year <- min(day_of_year, origin_days)

  new_start <- as.Date(paste0(origin_year, sprintf("-%03d", day_of_year)), format = "%Y-%j")
  day_diff <- as.numeric(first_date - new_start)

  if (inherits(x, "Date")) {
    return(x - day_diff)
  }

  x - day_diff * 86400

}

# Internal function
# Create one temporal pattern for each label. The function 'resampling_fun' is
# called for each label and band with a data.frame that has the time index 'x'
# and the values 'y' of all time series with this label. The model that it
# returns gives the values of the temporal pattern
.resample_patterns <- function(ts, labels, resampling_fun, resampling_freq) {

  bands <- setdiff(names(ts[[1]]), "time")

  # Align time series of different years
  ts <- lapply(ts, function(s) {
    s$time <- .shift_dates(s$time)
    s
  })

  # One pattern for each label, in the order in which the labels first appear
  pattern_labels <- unique(labels)
  smooth_models <- vector("list", length(pattern_labels))
  names(smooth_models) <- pattern_labels
  patterns <- vector("list", length(pattern_labels))

  for (l in seq_along(pattern_labels)) {

    label_ts <- ts[labels == pattern_labels[l]]
    label_time <- do.call(c, lapply(label_ts, function(s) s$time))

    # Fit a model to the observations of each band
    smooth_models[[l]] <- lapply(setNames(bands, bands), function(band) {
      values <- unlist(lapply(label_ts, function(s) s[[band]]), use.names = FALSE)
      resampling_fun(data.frame(x = as.numeric(label_time), y = values))
    })

    # Time index of the temporal pattern: the dates of the observations in time
    # order, or the numeric time indices in the order of the observations
    if (is.null(resampling_freq)) {
      pattern_time <- unique(label_time)
      if (!is.numeric(pattern_time)) {
        pattern_time <- sort(pattern_time)
      }
    } else {
      pattern_time <- seq(min(label_time), max(label_time), by = resampling_freq)
    }

    pattern_values <- lapply(setNames(bands, bands), function(band) {
      .predict_pattern(smooth_models[[l]][[band]], pattern_time, pattern_labels[l], band)
    })

    patterns[[l]] <- data.frame(time = pattern_time, pattern_values, check.names = FALSE)

  }

  list(patterns = patterns, labels = pattern_labels, smooth_models = smooth_models)

}

# Internal function
# Values of a temporal pattern at the time index 'pattern_time', predicted with
# a model returned by 'resampling_fun'
.predict_pattern <- function(model, pattern_time, label, band) {

  pattern_x <- as.numeric(pattern_time)

  # A model that uses the time index as a factor, e.g. y ~ factor(x), needs the
  # levels of the factor to predict
  model_x <- tryCatch(model.frame(model)[[2]], error = function(e) NULL)
  if (is.factor(model_x)) {
    pattern_x <- factor(pattern_x, levels = levels(model_x))
  }

  values <- tryCatch(predict(model, newdata = data.frame(x = pattern_x)), error = function(e) e)

  if (inherits(values, "error") || !is.numeric(values) || length(values) != length(pattern_time)) {
    stop(sprintf(paste("The model returned by 'resampling_fun' for the label '%s' and the band '%s' does not",
                       "predict one numeric value for each row of 'newdata', see the details of twdtw_knn%s"),
                 label, band, if (inherits(values, "error")) paste0(": ", conditionMessage(values)) else ""),
         call. = FALSE)
  }

  as.numeric(values)

}
