# twdtw_knn() creates a model from labelled time series in three forms: a list of
# data.frames or matrices, a data.frame with one column for each band and date,
# and an array with the dimensions series, time, and band.

test_that("the model keeps the labelled time series as temporal patterns", {
  m <- knn_model(knn_samples, knn_sample_labels)

  # A single class, so that predict functions of other packages find the method
  expect_identical(class(m), "twdtw_knn")
  expect_named(m, c("call", "smooth_models", "rmse", "data", "twdtw_args"))
  expect_null(m$smooth_models)
  expect_null(m$rmse)
  expect_identical(m$call[[1]], as.name("twdtw_knn"))

  expect_s3_class(m$data, "data.frame")
  expect_named(m$data, c("label", "observations"))
  expect_equal(nrow(m$data), length(knn_samples))
  expect_identical(m$data$label, knn_sample_labels)
  for (i in c(1, length(knn_samples))) {
    expect_equal(m$data$observations[[i]], knn_samples[[i]])
  }

  expect_equal(m$twdtw_args, list(time_weight = c(steepness = 0.1, midpoint = 50),
                                  cycle_length = 'year',
                                  time_scale = 'day',
                                  origin = NULL,
                                  max_elapsed = Inf))

  # Labels given as a factor are kept as character strings
  m_factor <- knn_model(knn_samples, factor(knn_sample_labels))
  expect_identical(m_factor$data, m$data)

  # Additional arguments of the TWDTW distance
  m_args <- knn_model(knn_samples, knn_sample_labels, max_elapsed = 60)
  expect_equal(m_args$twdtw_args$max_elapsed, 60)
})

test_that("the three forms of x give the same model", {
  m_list <- knn_model(knn_samples, knn_sample_labels)
  m_array <- knn_model(knn_array(knn_samples), knn_sample_labels, time = knn_time)
  m_wide <- knn_model(knn_wide(knn_samples), knn_sample_labels)

  expect_equal(m_array$data, m_list$data)
  expect_equal(m_wide$data, m_list$data)
  expect_equal(m_array$twdtw_args, m_list$twdtw_args)
  expect_equal(m_wide$twdtw_args, m_list$twdtw_args)

  # Matrices with the time index in the first column
  ts_matrix <- lapply(knn_samples, function(ts) cbind(as.numeric(format(ts$time, "%j")), as.matrix(ts[knn_bands])))
  m_matrix <- knn_model(ts_matrix, knn_sample_labels, index_column = 1)
  expect_equal(lapply(m_matrix$data$observations, function(ts) ts[knn_bands]),
               lapply(m_list$data$observations, function(ts) ts[knn_bands]))
  expect_type(m_matrix$data$observations[[1]]$time, "double")

  # An array without band names
  x <- unname(knn_array(knn_samples))
  m_unnamed <- knn_model(x, knn_sample_labels, time = knn_time)
  expect_named(m_unnamed$data$observations[[1]], c("time", paste0("V", 1:6)))
})

test_that("time series in a list can have different dates, bands in another order, and rows in any order", {
  samples <- knn_samples[1:6]
  samples[[2]] <- samples[[2]][23:1, ]
  samples[[3]] <- samples[[3]][c("NDVI", "time", "MIR", "EVI", "RED", "NIR", "BLUE")]
  samples[[4]] <- samples[[4]][seq(1, 23, by = 2), ]
  samples[[5]]$time <- samples[[5]]$time + 365
  samples[[6]]$time <- as.character(samples[[6]]$time)

  m <- knn_model(samples, knn_sample_labels[1:6])

  # Observations in time order and bands in the order of the first time series
  expect_equal(m$data$observations[[2]], knn_samples[[2]], ignore_attr = "row.names")
  expect_equal(m$data$observations[[3]], knn_samples[[3]])
  expect_equal(nrow(m$data$observations[[4]]), 12)
  expect_equal(m$data$observations[[5]]$time, knn_time + 365)
  expect_equal(m$data$observations[[6]], knn_samples[[6]])
})

test_that("a data.frame has one column for each band and date in several name formats", {
  samples <- knn_samples[1:8]
  labels <- knn_sample_labels[1:8]
  m <- knn_model(samples, labels)

  # Dates written with periods, as data.frame() writes the names
  wide <- knn_wide(samples)
  names(wide) <- gsub("-", ".", names(wide), fixed = TRUE)
  expect_equal(knn_model(wide, labels)$data, m$data)

  # The date first, with and without the X that data.frame() adds
  wide <- knn_wide(samples)
  date_first <- sub("^(.+)\\.([0-9-]+)$", "\\2.\\1", names(wide))
  names(wide) <- date_first
  expect_equal(knn_model(wide, labels)$data, m$data)
  names(wide) <- paste0("X", gsub("-", ".", date_first, fixed = TRUE))
  expect_equal(knn_model(wide, labels)$data, m$data)

  # Dates that are not in time order, and columns that are not observations
  wide <- knn_wide(samples)
  later_first <- unlist(lapply(seq_along(knn_bands) - 1, function(band) band * 23 + 23:1))
  wide <- cbind(x = 1:8, wide[later_first], y = 8:1)
  expect_equal(knn_model(wide, labels)$data, m$data)

  # The bands are in the order of the columns
  wide <- knn_wide(samples)
  expect_named(knn_model(wide[rev(seq_along(wide))], labels)$data$observations[[1]], c("time", rev(knn_bands)))

  # Columns with only the date are the observations of a single band
  wide <- knn_wide(samples, bands = "EVI")
  names(wide) <- sub("^EVI\\.", "", names(wide))
  m_single <- knn_model(wide, labels)
  expect_named(m_single$data$observations[[1]], c("time", "value"))
  expect_equal(m_single$data$observations[[3]]$value, samples[[3]]$EVI)

  # A band name with a hyphen is kept
  wide <- knn_wide(samples, bands = "EVI")
  names(wide) <- sub("^EVI\\.", "EVI-2.", names(wide))
  expect_named(knn_model(wide, labels)$data$observations[[1]], c("time", "EVI-2"))

  # A data.frame with a column for the time index is a single time series and must be in a list
  expect_error(knn_model(samples[[1]], labels[1]),
               "'x' has no columns named with a band and a date", fixed = TRUE)
  wide <- knn_wide(samples)
  wide[[3]] <- as.character(wide[[3]])
  expect_error(knn_model(wide, labels), "is not numeric", fixed = TRUE)
})

test_that("start and end keep the observations in the time interval of each label", {
  samples <- knn_samples[1:5]
  labels <- knn_sample_labels[1:5]

  # One interval for all time series, including the start and end dates
  m <- knn_model(samples, labels, start = knn_time[3], end = knn_time[20])
  for (ts in m$data$observations) {
    expect_equal(ts$time, knn_time[3:20])
  }
  expect_equal(m$data$observations[[2]], samples[[2]][3:20, ], ignore_attr = "row.names")

  # Only start or only end, given as character strings
  m <- knn_model(samples, labels, start = as.character(knn_time[10]))
  expect_equal(m$data$observations[[1]]$time, knn_time[10:23])
  m <- knn_model(samples, labels, end = as.character(knn_time[10]))
  expect_equal(m$data$observations[[1]]$time, knn_time[1:10])

  # One interval for each time series, in all forms of x
  start <- knn_time[c(1, 2, 3, 4, 5)]
  end <- knn_time[c(23, 22, 21, 20, 19)] + 3
  m <- knn_model(samples, labels, start = start, end = end)
  expect_equal(vapply(m$data$observations, nrow, integer(1)), c(23L, 21L, 19L, 17L, 15L))
  expect_equal(knn_model(knn_array(samples), labels, time = knn_time, start = start, end = end)$data, m$data)
  expect_equal(knn_model(knn_wide(samples), labels, start = start, end = end)$data, m$data)

  expect_error(knn_model(samples, labels, start = knn_time[1:2]),
               "'start' must have one value or one value for each time series in 'x' (5)", fixed = TRUE)
  expect_error(knn_model(samples, labels, end = as.Date(NA)),
               "'end' cannot have missing values", fixed = TRUE)
  expect_error(knn_model(samples, labels, start = 15000),
               "'start' must have the class of the time index of the time series (Date)", fixed = TRUE)
})

test_that("time series without complete observations are removed with a warning", {
  samples <- knn_samples[1:5]
  labels <- paste("label", 1:5)

  # No observation in the time interval
  start <- knn_time[c(1, 1, 1, 1, 1)]
  end <- c(knn_time[23], knn_time[1] - 10, knn_time[23], knn_time[1] - 1, knn_time[1])
  expect_warning(m <- knn_model(samples, labels, start = start, end = end),
                 "Time series without complete observations between 'start' and 'end' were removed: 2 of 5",
                 fixed = TRUE)
  expect_identical(m$data$label, paste("label", c(1, 3, 5)))
  expect_equal(nrow(m$data$observations[[3]]), 1)

  # A missing value in every observation
  samples[[4]]$NDVI <- NA
  samples[[2]]$EVI[1:22] <- NA
  expect_warning(m <- knn_model(samples, labels),
                 "Time series without complete observations were removed: 1 of 5", fixed = TRUE)
  expect_identical(m$data$label, paste("label", c(1, 2, 3, 5)))

  # Observations with missing values are kept in the model and not used by predict
  expect_equal(sum(complete.cases(m$data$observations[[2]])), 1)

  expect_error(knn_model(samples, labels, start = knn_time[23] + 1),
               "No time series in 'x' has a complete observation between 'start' and 'end'", fixed = TRUE)
})

test_that("numeric time indices are in the units of the cycle", {
  samples <- lapply(knn_samples[1:6], function(ts) {
    ts$time <- as.numeric(format(ts$time, "%j"))
    ts
  })

  # The time scale is not needed for a numeric cycle
  m <- twdtw_knn(samples, knn_sample_labels[1:6], time_weight = c(0.1, 50), cycle_length = 366)
  expect_null(m$twdtw_args$time_scale)

  # The time series are from September to August: the day of the year starts again
  # within the time series, and the observations are kept in their order
  expect_equal(m$data$observations[[1]]$time, as.numeric(format(knn_time, "%j")))
  expect_true(is.unsorted(m$data$observations[[1]]$time))

  # An array with numeric time and the TWDTW arguments given by position
  x <- knn_array(knn_samples[1:6])
  m_array <- twdtw_knn(x, knn_sample_labels[1:6], c(0.1, 50), 366, time = seq_along(knn_time))
  expect_equal(m_array$data$observations[[1]]$time, seq_along(knn_time))
  expect_null(m_array$twdtw_args$time_scale)
})

test_that("arguments that cannot give a model are an error", {
  samples <- knn_samples[1:4]
  labels <- knn_sample_labels[1:4]

  # Labels
  expect_error(knn_model(samples, labels[1:3]),
               "'y' must be a vector with one label for each time series in 'x' (4)", fixed = TRUE)
  expect_error(knn_model(samples, c(labels[1:3], NA)), "'y' cannot have missing values", fixed = TRUE)

  # Arguments that are not known are not ignored
  expect_error(knn_model(samples, labels, resample = 16),
               "Unknown arguments in twdtw_knn: 'resample'.", fixed = TRUE)
  expect_error(knn_model(samples, labels, time = knn_time),
               "Unknown arguments in twdtw_knn: 'time'. 'time' is only used if 'x' is an array.", fixed = TRUE)
  expect_error(knn_model(knn_wide(samples), labels, index_column = "time"),
               "Unknown arguments in twdtw_knn: 'index_column'.", fixed = TRUE)

  # The number of nearest neighbours is an argument of predict. It is not taken
  # as a short form of 'keep_models', which must be given with its full name
  expect_error(knn_model(samples, labels, k = 3), "Unknown arguments in twdtw_knn: 'k'.", fixed = TRUE)
  expect_error(knn_model(samples, labels, keep = FALSE), "Unknown arguments in twdtw_knn: 'keep'.", fixed = TRUE)
  for (value in list(NA, "no", c(TRUE, FALSE), NULL)) {
    expect_error(knn_model(samples, labels, keep_models = value), "'keep_models' must be TRUE or FALSE", fixed = TRUE)
  }

  # Time series
  expect_error(knn_model(list(), character(0)), "'x' has no time series", fixed = TRUE)
  expect_error(knn_model(c(samples[1:3], list(1:10)), labels),
               "Element 4 of 'x' is neither a data.frame nor a matrix", fixed = TRUE)
  expect_error(knn_model(samples, labels, index_column = "date"),
               "The 'index_column' was not found in element 1 of 'x'", fixed = TRUE)
  expect_error(knn_model(c(samples[1:3], list(samples[[4]][c("time", "EVI")])), labels),
               "The bands of element 4 of 'x' (EVI) are not the bands of the first time series", fixed = TRUE)
  samples_character <- samples
  samples_character[[2]]$EVI <- as.character(samples_character[[2]]$EVI)
  expect_error(knn_model(samples_character, labels), "Element 2 of 'x' has non-numeric columns", fixed = TRUE)
  samples_numeric <- samples
  samples_numeric[[3]]$time <- as.numeric(samples_numeric[[3]]$time)
  expect_error(knn_model(samples_numeric, labels),
               "The time indices of the time series in 'x' must be all numeric, all dates, or all date-times",
               fixed = TRUE)

  # Arrays
  x <- knn_array(samples)
  expect_error(knn_model(x, labels), "The 'time' argument is missing", fixed = TRUE)
  expect_error(knn_model(x, labels, time = knn_time[-1]),
               "The length of 'time' must be equal to the second dimension of 'x'", fixed = TRUE)
  expect_error(knn_model(x[, , 1], labels, time = knn_time),
               "'x' must be a numeric array with three dimensions: series, time, and band", fixed = TRUE)
  dimnames(x)[[3]][2] <- "time"
  expect_error(knn_model(x, labels, time = knn_time),
               "A band of 'x' is named 'time', which is the name of the time index", fixed = TRUE)
  expect_error(knn_model(1:10, labels), "No method of twdtw_knn for 'x' of class 'integer'", fixed = TRUE)

  # Arguments of the TWDTW distance
  expect_error(twdtw_knn(samples, labels, time_weight = 1, cycle_length = 'year', time_scale = 'day'),
               "'time_weight' should be either a function or a numeric vector with length two", fixed = TRUE)
  expect_error(twdtw_knn(samples, labels, time_weight = c(0.1, 50), cycle_length = 'year'),
               "The 'time_scale' argument is missing for 'cycle_length' type character.", fixed = TRUE)
  expect_error(twdtw_knn(samples, labels, time_weight = c(0.1, 50), cycle_length = 'year', time_scale = 'week'),
               "Invalid time_scale for the provided cycle_length", fixed = TRUE)
  expect_error(twdtw_knn(samples, labels, time_weight = c(0.1, 50), cycle_length = 365, time_scale = 'day'),
               "For numeric cycle_length, origin must be given.", fixed = TRUE)
  expect_error(knn_model(samples, labels, max_elapsed = NA), "'max_elapsed' must be a single numeric value", fixed = TRUE)
  expect_error(knn_model(samples, labels, resampling_fun = "mean"), "'resampling_fun' must be a function", fixed = TRUE)
})

test_that("the internal functions called for each time series give the results of the functions they replace", {
  time <- as.Date("2020-01-01") + c(0, 16, 32, 48)
  values <- matrix(c(1.5, 2, NA, 4, 5, 6, 7, 8), ncol = 2, dimnames = list(NULL, c("RED", "NIR")))

  # A data.frame with the time index and one column for each band
  inputs <- list(
    dates = list(time, values),
    date_times = list(as.POSIXct("2020-01-01 10:00:00", tz = "America/Sao_Paulo") + (0:3) * 86400, values),
    numeric_time = list(c(1, 17, 33, 49), values),
    integer_values = list(time, matrix(1:8, ncol = 2, dimnames = list(NULL, c("RED", "NIR")))),
    one_band = list(time, values[, "NIR", drop = FALSE]),
    one_observation = list(time[1], values[1, , drop = FALSE]),
    no_observation = list(time[0], values[0, , drop = FALSE]),
    named_time = list(setNames(time, c("a", "b", "c", "d")), values))

  for (case in names(inputs)) {
    expect_true(identical(twdtw:::.time_series_frame(inputs[[case]][[1]], inputs[[case]][[2]]),
                          data.frame(time = inputs[[case]][[1]], inputs[[case]][[2]], check.names = FALSE)),
                info = case)
  }

  # The observations of a time series that are kept
  ts <- data.frame(time = time, values, check.names = FALSE)
  keep <- list(all = rep(TRUE, 4), none = rep(FALSE, 4), first = c(TRUE, TRUE, FALSE, FALSE),
               last = c(FALSE, FALSE, FALSE, TRUE), some = c(FALSE, TRUE, FALSE, TRUE))

  for (case in names(keep)) {
    expect_true(identical(twdtw:::.keep_observations(ts, keep[[case]]), ts[keep[[case]], , drop = FALSE]),
                info = case)
  }

  # A time series with at least one complete observation
  series <- list(complete = ts[c(1, 2, 4), ], some_missing = ts, none_complete = ts[3, ],
                 no_observation = ts[0, ], band_missing = transform(ts, RED = NA_real_))

  for (case in names(series)) {
    expect_identical(twdtw:::.has_complete_observation(series[[case]]), any(complete.cases(series[[case]])),
                     info = case)
  }
})

test_that("the model is printed with its patterns and the arguments of the distance", {
  m <- knn_model(knn_samples, knn_sample_labels)
  out <- capture.output(printed <- print(m))
  expect_identical(printed, m)

  expect_match(out, "Model of class 'twdtw_knn'", fixed = TRUE, all = FALSE)
  expect_match(out, "twdtw_knn(x = x, y = y, time_weight = c(steepness = 0.1, midpoint = 50)", fixed = TRUE, all = FALSE)
  expect_match(out, paste(length(knn_samples), "temporal patterns, the labelled time series"), fixed = TRUE, all = FALSE)
  expect_match(out, "Bands: EVI, NDVI, RED, BLUE, NIR, MIR", fixed = TRUE, all = FALSE)
  for (label in unique(knn_sample_labels)) {
    expect_match(out, paste0(label, " +", sum(knn_sample_labels == label), " +23 +2011-09-14 +2012-08-28"), all = FALSE)
  }
  expect_match(out, " - time_weight: c(steepness=0.1, midpoint=50)", fixed = TRUE, all = FALSE)
  expect_match(out, " - cycle_length: year", fixed = TRUE, all = FALSE)
  expect_match(out, " - origin: NULL", fixed = TRUE, all = FALSE)
  expect_match(out, " - max_elapsed: Inf", fixed = TRUE, all = FALSE)
  expect_false(any(grepl("RMSE", out, fixed = TRUE)))

  # A time-weight function and patterns created with a resampling function
  m <- twdtw_knn(knn_samples, knn_sample_labels,
                 time_weight = function(dist, el) dist + 0.1 * el,
                 cycle_length = 'year', time_scale = 'day',
                 resampling_fun = function(data) lm(y ~ factor(x), data = data))
  out <- capture.output(print(m))
  expect_match(out, "temporal patterns, one for each label created with 'resampling_fun'", fixed = TRUE, all = FALSE)
  expect_match(out, "Root Mean Squared Error (RMSE) of smooth models:", fixed = TRUE, all = FALSE)
  expect_match(out, " - time_weight: function (dist, el)", fixed = TRUE, all = FALSE)

  # A time weight without names
  m <- twdtw_knn(knn_samples, knn_sample_labels, time_weight = c(0.1, 50), cycle_length = 'year', time_scale = 'day')
  expect_match(capture.output(print(m)), " - time_weight: c(0.1, 50)", fixed = TRUE, all = FALSE)

  # A call that has the data instead of the names of the objects is shortened
  m <- do.call(twdtw_knn, list(x = knn_samples[1:20], y = knn_sample_labels[1:20],
                               time_weight = c(0.1, 50), cycle_length = 'year', time_scale = 'day'))
  out <- capture.output(print(m))
  expect_lt(length(out), 40)
  expect_match(out, "^\\.\\.\\.$", all = FALSE)
})
