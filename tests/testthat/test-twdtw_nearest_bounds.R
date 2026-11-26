# twdtw_nearest() skips patterns using a lower bound of their distance if k is lower than
# the number of patterns. With k equal to the number of patterns no pattern can be skipped,
# which gives the distance of all pairs to compare with
ts_x <- readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))
ts_y <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))

bands <- c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")
time <- ts_x$observations[[1]]$time
tw <- c(steepness = 0.1, midpoint = 50)

# Array of time series with dimensions: series, time, and band
x_all <- sapply(ts_x$observations, function(ts) as.matrix(ts[bands]), simplify = "array")
x_all <- aperm(x_all, c(3, 1, 2))

# 60 complete time series are the patterns, 300 other time series are classified
complete <- which(!sapply(ts_x$observations, anyNA))
pattern_ids <- complete[seq(1, length(complete), length.out = 60)]
patterns <- lapply(pattern_ids, function(i) as.data.frame(ts_x$observations[[i]]))
x <- x_all[setdiff(seq_len(dim(x_all)[1]), pattern_ids)[seq(1, 900, by = 3)], , ]

nearest <- function(x, y = patterns, k = 1, time_weight = tw, ...) {
  twdtw_nearest(x = x, y = y, time = time, cycle_length = 'year', time_scale = 'day',
                time_weight = time_weight, k = k, ...)
}

# The k nearest patterns must be the first k of all patterns
expect_first_of_all <- function(x, y = patterns, k = c(1, 3, 10), ...) {
  all_patterns <- nearest(x, y, k = length(y), ...)
  for (ki in k) {
    knn <- nearest(x, y, k = ki, ...)
    expect_identical(knn$index, all_patterns$index[, seq_len(ki), drop = FALSE])
    expect_identical(knn$distance, all_patterns$distance[, seq_len(ki), drop = FALSE])
  }
  invisible(all_patterns)
}

test_that("skipping patterns does not change the nearest neighbours", {
  all_patterns <- expect_first_of_all(x)

  # The distance is the same as calling twdtw for each pair
  ids <- c(1, 150, 300)
  distance <- t(sapply(ids, function(i) {
    sapply(patterns, function(pattern) {
      twdtw(x = data.frame(time = time, x[i, , ]), y = pattern, cycle_length = 'year',
            time_scale = 'day', time_weight = tw)[1]
    })
  }))
  expect_equal(all_patterns$distance[ids, ], t(apply(distance, 1, sort)))
  expect_equal(nearest(x[ids, , ])$index[, 1], apply(distance, 1, which.min))
})

test_that("skipping patterns does not change the nearest neighbours for other time weights", {
  # A steep weight, a weak weight, and a weight that decreases with the elapsed time
  expect_first_of_all(x, time_weight = c(steepness = 5, midpoint = 30))
  expect_first_of_all(x, time_weight = c(steepness = 0.01, midpoint = 100))
  expect_first_of_all(x, time_weight = c(steepness = -0.1, midpoint = 50))

  # Without time weight
  expect_first_of_all(x, time_weight = c(steepness = 1, midpoint = 1e12))
})

test_that("skipping patterns does not change the nearest neighbours using max_elapsed", {
  expect_first_of_all(x, max_elapsed = 60)

  # Patterns with dates from another season cannot be matched with a short max_elapsed
  shifted <- lapply(patterns[1:20], function(pattern) {
    pattern$time <- pattern$time + 120
    pattern
  })
  all_patterns <- expect_first_of_all(x, y = c(shifted, patterns[21:60]), max_elapsed = 20)
  expect_true(all(is.infinite(all_patterns$distance[, 41:60])))

  # Fewer patterns with a finite distance than neighbours
  expect_first_of_all(x, y = c(shifted, patterns[21:25]), k = c(3, 10), max_elapsed = 20)
})

test_that("skipping patterns does not change the nearest neighbours of time series with missing values", {
  x_na <- x
  set.seed(3)
  x_na[sample(length(x_na), 1000)] <- NA
  x_na[5, , ] <- NA
  x_na[6, -1, ] <- NA
  all_patterns <- expect_first_of_all(x_na)
  expect_true(all(is.na(all_patterns$index[5, ])))
  expect_false(anyNA(all_patterns$index[-5, ]))
})

test_that("patterns with the same distance are in the order of their position", {
  # Each pattern is three times in the list
  repeated <- c(patterns[1:20], patterns[20:1], patterns[1:20])
  all_patterns <- expect_first_of_all(x, y = repeated, k = c(1, 2, 3, 4, 7))
  knn <- nearest(x, y = repeated, k = 3)
  expect_true(all(knn$distance[, 1] == knn$distance[, 3]))
  expect_true(all(knn$index[, 1] < knn$index[, 2] & knn$index[, 2] < knn$index[, 3]))
})

test_that("patterns can have other dates and numbers of observations than the time series", {
  # Patterns with observations every eight days, with every second, and with a part of the observations
  resampled <- lapply(patterns[1:20], function(pattern) {
    dates <- seq(min(pattern$time), max(pattern$time), by = 8)
    out <- data.frame(time = dates)
    for (band in bands) out[[band]] <- approx(as.numeric(pattern$time), pattern[[band]], as.numeric(dates))$y
    out
  })
  thinned <- lapply(patterns[21:40], function(pattern) pattern[seq(1, nrow(pattern), by = 2), ])
  parts <- lapply(patterns[41:60], function(pattern) pattern[8:15, ])
  expect_first_of_all(x, y = c(resampled, thinned, parts, ts_y$observations))
})

test_that("the nearest neighbours do not depend on the order and the number of the time series", {
  knn <- nearest(x, k = 3)

  set.seed(7)
  shuffled <- sample(dim(x)[1])
  knn_shuffled <- nearest(x[shuffled, , ], k = 3)
  expect_identical(unname(knn_shuffled$index), unname(knn$index[shuffled, ]))
  expect_identical(unname(knn_shuffled$distance), unname(knn$distance[shuffled, ]))

  knn_single <- nearest(x[17, , , drop = FALSE], k = 3)
  expect_identical(unname(knn_single$index), unname(knn$index[17, , drop = FALSE]))
  expect_identical(unname(knn_single$distance), unname(knn$distance[17, , drop = FALSE]))
})
