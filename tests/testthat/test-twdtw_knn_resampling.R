# twdtw_knn() with a resampling function creates one temporal pattern for each label
# from all time series with this label.

# Dates of the example time series shifted to start on the same day of the year in 1970
shifted_time <- knn_time - as.numeric(knn_time[1] - as.Date("1970-09-14"))

test_that("the resampling function creates one temporal pattern for each label", {
  m <- knn_model(knn_samples, knn_sample_labels,
                 resampling_fun = function(data) lm(y ~ factor(x), data = data), keep_models = TRUE)

  # One pattern for each label, in the order in which the labels first appear
  labels <- unique(knn_sample_labels)
  expect_identical(m$data$label, labels)
  expect_equal(nrow(m$data), length(labels))

  # The model of each label and band
  expect_named(m$smooth_models, labels)
  for (label_models in m$smooth_models) {
    expect_named(label_models, knn_bands)
    expect_s3_class(label_models$EVI, "lm")
  }

  # The root mean squared error of each model, with one row for each band and one column for each label
  expect_identical(dimnames(m$rmse), list(knn_bands, labels))
  for (l in seq_along(labels)) {
    expect_equal(unname(m$rmse[, l]),
                 unname(sapply(m$smooth_models[[l]], function(model) sqrt(mean(residuals(model)^2)))))
  }

  # With the dates as a factor, the pattern is the mean of the samples on each date
  x <- knn_array(knn_samples)
  for (l in seq_along(labels)) {
    pattern <- m$data$observations[[l]]
    expect_named(pattern, c("time", knn_bands))
    expect_equal(pattern$time, shifted_time)
    expected <- apply(x[knn_sample_labels == labels[l], , , drop = FALSE], c(2, 3), mean)
    expect_equal(as.matrix(pattern[knn_bands]), expected, ignore_attr = TRUE)
  }

  # The three forms of x give the same patterns
  mean_by_date <- function(data) lm(y ~ factor(x), data = data)
  m_array <- knn_model(knn_array(knn_samples), knn_sample_labels, time = knn_time, resampling_fun = mean_by_date)
  m_wide <- knn_model(knn_wide(knn_samples), knn_sample_labels, resampling_fun = mean_by_date)
  expect_equal(m_array$data, m$data)
  expect_equal(m_wide$data, m$data)
})

test_that("the model keeps the root mean squared error of the fitted models, and the models only with keep_models = TRUE", {
  mean_by_date <- function(data) lm(y ~ factor(x), data = data)
  m <- knn_model(knn_samples, knn_sample_labels, resampling_fun = mean_by_date, keep_models = TRUE)
  m_small <- knn_model(knn_samples, knn_sample_labels, resampling_fun = mean_by_date)

  # The default is keep_models = FALSE, for the three forms of x
  for (class in c("list", "data.frame", "array")) {
    expect_identical(formals(getS3method("twdtw_knn", class))$keep_models, FALSE)
  }
  m_false <- knn_model(knn_samples, knn_sample_labels, resampling_fun = mean_by_date, keep_models = FALSE)
  expect_identical(m_small[names(m_small) != "call"], m_false[names(m_false) != "call"])
  expect_named(m_small, names(m))
  expect_null(m_small$smooth_models)
  expect_named(m$smooth_models, m$data$label)
  expect_identical(m_small$rmse, m$rmse)
  expect_identical(m_small$data, m$data)
  expect_lt(as.numeric(object.size(m_small)), as.numeric(object.size(m)) / 10)

  # The same predictions, and the same summary after the call
  new_series <- knn_series[knn_new_ids]
  expect_identical(predict(m_small, newdata = new_series), predict(m, newdata = new_series))
  printed_data <- function(model) {
    out <- capture.output(print(model))
    out[seq(which(out == "Data:"), length(out))]
  }
  expect_identical(printed_data(m_small), printed_data(m))
  expect_match(printed_data(m_small), "one for each label created with 'resampling_fun'", fixed = TRUE, all = FALSE)
  expect_match(printed_data(m_small), "Root Mean Squared Error (RMSE) of smooth models:", fixed = TRUE, all = FALSE)

  # A model that was created before the error was kept has only the models
  m_before <- m
  m_before$rmse <- NULL
  expect_identical(printed_data(m_before), printed_data(m))

  # The other two forms of x
  m_array <- knn_model(knn_array(knn_samples), knn_sample_labels, time = knn_time,
                       resampling_fun = mean_by_date, keep_models = FALSE)
  m_wide <- knn_model(knn_wide(knn_samples), knn_sample_labels, resampling_fun = mean_by_date)
  for (model in list(m_array, m_wide)) {
    expect_null(model$smooth_models)
    expect_equal(model$rmse, m$rmse)
    expect_equal(model$data, m$data)
  }
  m_array <- knn_model(knn_array(knn_samples), knn_sample_labels, time = knn_time,
                       resampling_fun = mean_by_date, keep_models = TRUE)
  m_wide <- knn_model(knn_wide(knn_samples), knn_sample_labels, resampling_fun = mean_by_date, keep_models = TRUE)
  for (model in list(m_array, m_wide)) {
    expect_named(model$smooth_models, m$data$label)
    expect_s3_class(model$smooth_models[[1]]$EVI, "lm")
  }

  # A model that keeps all time series has no models to keep
  m_all <- knn_model(knn_samples, knn_sample_labels, keep_models = TRUE)
  expect_null(m_all$smooth_models)
  expect_null(m_all$rmse)
  expect_identical(m_all$data, knn_model(knn_samples, knn_sample_labels)$data)

  # A model without residuals has no root mean squared error
  registerS3method("predict", "knn_test_mean", function(object, newdata, ...) rep(object$mean, nrow(newdata)))
  m_mean <- knn_model(knn_samples, knn_sample_labels,
                      resampling_fun = function(data) structure(list(mean = mean(data$y)), class = "knn_test_mean"))
  expect_identical(dimnames(m_mean$rmse), dimnames(m$rmse))
  expect_true(all(is.na(m_mean$rmse)))
  expect_equal(m_mean$data$observations[[1]]$EVI,
               rep(mean(knn_array(knn_samples)[knn_sample_labels == m_mean$data$label[1], , "EVI"]), length(knn_time)))
})

test_that("the resampling function receives the time index and the values of all samples of a label", {
  samples <- knn_samples[1:6]
  labels <- c("a", "b", "a", "b", "a", "b")
  received <- list()
  resampling_fun <- function(data) {
    received[[length(received) + 1]] <<- data
    lm(y ~ poly(x, 3), data = data)
  }

  m <- knn_model(samples, labels, resampling_fun = resampling_fun, keep_models = TRUE)

  # One call for each label and band, with the samples of the label one after the other
  expect_length(received, 2 * length(knn_bands))
  expect_named(received[[1]], c("x", "y"))
  expect_equal(received[[1]]$x, rep(as.numeric(shifted_time), 3))
  expect_equal(received[[1]]$y, c(samples[[1]]$EVI, samples[[3]]$EVI, samples[[5]]$EVI))
  expect_equal(received[[length(knn_bands) + 2]]$y, c(samples[[2]]$NDVI, samples[[4]]$NDVI, samples[[6]]$NDVI))

  # The pattern has the values that the model predicts for its dates
  pattern <- m$data$observations[[2]]
  expected <- predict(m$smooth_models$b$RED, newdata = data.frame(x = as.numeric(shifted_time)))
  expect_equal(pattern$RED, as.numeric(expected))
})

test_that("resampling_freq sets the time between the observations of the patterns", {
  smooth <- function(data) loess(y ~ x, data = data, span = 0.5)

  m <- knn_model(knn_samples, knn_sample_labels, resampling_fun = smooth, resampling_freq = 8, keep_models = TRUE)
  pattern_time <- seq(min(shifted_time), max(shifted_time), by = 8)
  for (pattern in m$data$observations) {
    expect_equal(pattern$time, pattern_time)
    expect_false(anyNA(pattern))
  }
  expected <- predict(m$smooth_models[[1]]$NIR, newdata = data.frame(x = as.numeric(pattern_time)))
  expect_equal(m$data$observations[[1]]$NIR, as.numeric(expected))

  # Without resampling_freq the patterns have the dates of the observations
  m <- knn_model(knn_samples, knn_sample_labels, resampling_fun = smooth)
  expect_equal(m$data$observations[[1]]$time, shifted_time)

  # A model that uses the dates as a factor has no value between the dates of the observations
  m <- knn_model(knn_samples, knn_sample_labels, resampling_freq = 8,
                 resampling_fun = function(data) lm(y ~ factor(x), data = data))
  expect_equal(is.na(m$data$observations[[1]]$EVI), !pattern_time %in% shifted_time)
})

test_that("time series of different years overlap in the patterns", {
  # The dates of the example time series, which start in 2011, on the same days of the
  # year starting in another year
  time_from <- function(first_year) {
    year <- first_year + as.numeric(format(knn_time, "%Y")) - 2011
    as.Date(paste(year, format(knn_time, "%j")), format = "%Y %j")
  }
  expect_equal(time_from(2011), knn_time)

  # The same time series in three years
  ts <- knn_samples[[1]]
  samples <- lapply(c(2011, 2013, 2014), function(first_year) {
    ts$time <- time_from(first_year)
    ts
  })

  m <- knn_model(samples, rep("a", 3), resampling_fun = function(data) lm(y ~ factor(x), data = data))
  expect_equal(nrow(m$data), 1)
  expect_equal(m$data$observations[[1]]$time, shifted_time)
  expect_equal(m$data$observations[[1]][knn_bands], ts[knn_bands], ignore_attr = TRUE)

  # The number of days between the observations is kept: a time series that starts in a
  # leap year has its observations of the next year one day later
  leap <- ts
  leap$time <- time_from(2016)
  next_year <- format(knn_time, "%Y") == "2012"
  m <- knn_model(list(ts, leap), rep("a", 2), resampling_fun = function(data) lm(y ~ factor(x), data = data))
  expect_equal(m$data$observations[[1]]$time, sort(c(shifted_time, shifted_time[next_year] + 1)))

  # The labels of an array are valid in a different time interval for each time series
  first_years <- c(2011, 2013, 2014)
  cube_time <- do.call(c, lapply(first_years, time_from))
  x <- knn_array(lapply(knn_samples[1:6], function(ts) {
    do.call(rbind, lapply(first_years, function(first_year) ts))
  }))
  labels <- rep(c("a", "b"), 3)
  start <- do.call(c, lapply(first_years, function(first_year) time_from(first_year)[c(1, 1)]))
  end <- do.call(c, lapply(first_years, function(first_year) time_from(first_year)[c(23, 23)]))

  m <- knn_model(x, labels, time = cube_time, start = start, end = end,
                 resampling_fun = function(data) lm(y ~ factor(x), data = data))
  m_one_year <- knn_model(knn_samples[1:6], labels,
                          resampling_fun = function(data) lm(y ~ factor(x), data = data))
  expect_equal(m$data, m_one_year$data)

  # A time series that starts on day 366 of a leap year starts on the last day of 1970
  expect_equal(twdtw:::.shift_dates(as.Date(c("2012-12-31", "2013-01-16"))),
               as.Date(c("1970-12-31", "1971-01-16")))
  expect_equal(twdtw:::.shift_dates(as.Date(c("2011-09-14", "2011-09-30"))),
               as.Date(c("1970-09-14", "1970-09-30")))

  # Date-times keep their time of the day
  date_time <- as.POSIXct(c("2011-09-14 10:30:00", "2011-09-15 22:00:00"), tz = "UTC")
  expect_equal(twdtw:::.shift_dates(date_time),
               as.POSIXct(c("1970-09-14 10:30:00", "1970-09-15 22:00:00"), tz = "UTC"))
})

test_that("each time series is shifted according to its own first observation", {
  # Four time series with one observation every 16 days. The first two start at the same
  # time, which is 31 December in one time zone and 1 January in the other. The third
  # starts one year later, and the fourth is equal to the first
  first <- as.POSIXct("2020-12-31 20:00:00", tz = "UTC")
  series <- function(first, tz) {
    time <- first + (0:5) * 16 * 86400
    attr(time, "tzone") <- tz
    data.frame(time = time, v = as.numeric(1:6))
  }
  samples <- list(series(first, "UTC"), series(first, "Asia/Tokyo"),
                  series(first + 365 * 86400, "UTC"), series(first, "UTC"))

  received <- NULL
  m <- knn_model(samples, rep("a", 4), resampling_fun = function(data) {
    received <<- data
    lm(y ~ x, data = data)
  })

  # The model of the label receives the observations of all time series, each
  # shifted as by .shift_dates
  expected <- unlist(lapply(samples, function(ts) as.numeric(twdtw:::.shift_dates(ts$time))))
  expect_identical(received$x, expected)
  expect_identical(received$y, as.numeric(rep(1:6, 4)))

  # The first time series starts on the last day of 1970 and the second on the first
  time_1970 <- as.POSIXct(received$x[c(1, 7, 13, 19)], origin = "1970-01-01", tz = "UTC")
  expect_equal(format(time_1970, "%Y-%m-%d %H:%M"),
               c("1970-12-31 20:00", "1969-12-31 20:00", "1970-12-31 20:00", "1970-12-31 20:00"))
})

test_that("numeric time indices are not shifted and keep their order", {
  samples <- lapply(knn_samples[1:6], function(ts) {
    ts$time <- as.numeric(format(ts$time, "%j"))
    ts
  })
  day_of_year <- samples[[1]]$time

  m <- twdtw_knn(samples, rep(c("a", "b"), 3), time_weight = c(0.1, 50), cycle_length = 366,
                 resampling_fun = function(data) lm(y ~ factor(x), data = data))

  # The patterns start in September, as the time series
  expect_equal(m$data$observations[[1]]$time, day_of_year)
  expected <- (samples[[1]]$EVI + samples[[3]]$EVI + samples[[5]]$EVI) / 3
  expect_equal(m$data$observations[[1]]$EVI, expected)
})

test_that("a resampling function must return a model that predicts the values of the pattern", {
  expect_error(knn_model(knn_samples, knn_sample_labels,
                         resampling_fun = function(data) smooth.spline(data$x, data$y)),
               "The model returned by 'resampling_fun' for the label", fixed = TRUE)
  expect_error(knn_model(knn_samples, knn_sample_labels, resampling_fun = function(data) mean(data$y)),
               "The model returned by 'resampling_fun' for the label", fixed = TRUE)

  # A model that predicts only missing values gives a pattern that cannot be compared
  registerS3method("predict", "knn_test_missing", function(object, newdata, ...) rep(NA_real_, nrow(newdata)),
                   envir = asNamespace("stats"))
  expect_error(knn_model(knn_samples, knn_sample_labels,
                         resampling_fun = function(data) structure(list(), class = "knn_test_missing")),
               "The temporal pattern created with 'resampling_fun' for the label", fixed = TRUE)
})
