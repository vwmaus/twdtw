# k-nearest neighbours with a time-weight function defined in R

ts_x <- readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))
ts_y <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))

bands <- c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")
time <- ts_x$observations[[1]]$time
patterns <- ts_y$observations

# Array of time series with dimensions: series, time, and band
x <- sapply(ts_x$observations, function(ts) as.matrix(ts[bands]), simplify = "array")
x <- aperm(x, c(3, 1, 2))

# Time series with missing values and a sample of the remaining time series
has_na <- which(sapply(ts_x$observations, anyNA))
ids <- sort(unique(c(has_na, seq(1, dim(x)[1], by = 25))))

# Linear time weight
tw_linear <- function(dist, el) dist + 0.02 * el

# The cost can be negative
tw_negative <- function(dist, el) dist - 0.5 + 0.02 * el

# Linear time weight up to 60 days: not vectorised and vectorised
tw_if <- function(dist, el) if (el > 60) dist + 5 else dist + 0.02 * el
tw_ifelse <- function(dist, el) ifelse(el > 60, dist + 5, dist + 0.02 * el)

# No match of observations that are more than 45 days apart
tw_infinite <- function(dist, el) ifelse(el > 45, Inf, dist + 0.02 * el)

# The distance without a time weight
tw_none <- function(dist, el) dist

# Define TWDTW k-nearest neighbours call
twdtw_knn_call <- function(x, time_weight, y = patterns, ...) {
  twdtw_knn(x = x,
            y = y,
            time = time,
            cycle_length = 'year',
            time_scale = 'day',
            time_weight = time_weight, ...)
}

# Distance between each pair of time series and pattern calling twdtw for each pair
twdtw_pairs <- function(ids, time_weight, ...) {
  t(sapply(ids, function(i) {
    sapply(patterns, function(pattern) {
      twdtw(x = as.data.frame(ts_x$observations[[i]]),
            y = as.data.frame(pattern),
            cycle_length = 'year',
            time_scale = 'day',
            time_weight = time_weight, ...)[1]
    })
  }))
}

# The k lowest distances in each row and the respective columns
k_lowest <- function(distances, k) {
  list(index = t(apply(distances, 1, order))[, seq_len(k), drop = FALSE],
       distance = t(apply(distances, 1, sort))[, seq_len(k), drop = FALSE])
}

# Record the calls of a time-weight function: the number of calls, the length of
# the arguments in the first 100 calls, and the largest length after the first call
record_calls <- function(time_weight) {
  calls <- new.env()
  calls$n <- 0L
  calls$lengths <- integer(100)
  calls$largest_later <- 0L
  list(time_weight = function(dist, el) {
         calls$n <- calls$n + 1L
         if (calls$n <= 100L) calls$lengths[calls$n] <- length(dist)
         if (calls$n > 1L) calls$largest_later <- max(calls$largest_later, length(dist))
         time_weight(dist, el)
       },
       calls = calls)
}

test_that("a vectorised time-weight function is called once for each pattern", {
  recorded <- record_calls(tw_linear)
  knn <- twdtw_knn_call(x[ids, , , drop = FALSE], recorded$time_weight)

  # One call for each pattern with all pairs of observations, after the calls that test the function
  n_pairs <- length(ids) * length(time) * sapply(patterns, nrow)
  expect_lt(recorded$calls$n, 100)
  expect_equal(tail(recorded$calls$lengths[seq_len(recorded$calls$n)], length(patterns)), n_pairs)
})

test_that("a time-weight function that is not vectorised is called for each pair of observations", {
  x_small <- x[ids[1:12], , , drop = FALSE]

  recorded <- record_calls(tw_if)
  knn <- twdtw_knn_call(x_small, recorded$time_weight)

  # The first call tests the function with vectors, all other calls have one pair of observations
  expect_gt(recorded$calls$lengths[1], 1)
  expect_equal(recorded$calls$largest_later, 1)
  expect_gt(recorded$calls$n, dim(x_small)[1] * length(time))

  # The result is the same as for the vectorised function
  expect_equal(knn, twdtw_knn_call(x_small, tw_ifelse))
})

test_that("k-nearest neighbours are the same as calling twdtw for each pair", {
  for (time_weight in list(tw_linear, tw_negative, tw_ifelse)) {
    distances <- twdtw_pairs(ids, time_weight)

    for (k in c(1, 3, 5)) {
      knn <- twdtw_knn_call(x[ids, , , drop = FALSE], time_weight, k = k)
      expect_equal(knn, k_lowest(distances, k))
    }
  }
})

test_that("k-nearest neighbours are the same as calling twdtw for each pair using max_elapsed", {
  # Without a time weight, the distances depend on the maximum elapsed time
  unconstrained <- twdtw_knn_call(x[ids, , , drop = FALSE], tw_none, k = 5)

  for (max_elapsed in c(20, 30, 90)) {
    distances <- suppressWarnings(twdtw_pairs(ids, tw_none, max_elapsed = max_elapsed))

    for (k in c(1, 5)) {
      knn <- twdtw_knn_call(x[ids, , , drop = FALSE], tw_none, k = k, max_elapsed = max_elapsed)
      expect_equal(knn, k_lowest(distances, k))
    }

    expect_true(all(knn$distance >= unconstrained$distance, na.rm = TRUE))
    expect_false(isTRUE(all.equal(knn$distance, unconstrained$distance)))
  }
})

test_that("a time-weight function can return infinite values", {
  distances <- suppressWarnings(twdtw_pairs(ids, tw_infinite))
  knn <- twdtw_knn_call(x[ids, , , drop = FALSE], tw_infinite, k = 5)

  expect_equal(knn, k_lowest(distances, 5))
})

test_that("patterns can be matrices with different number of observations", {
  # Numeric time index in units of the cycle
  doy <- date_to_numeric_cycle(time, 'year', 'day')

  # Bands of the patterns in a different order than in the time series
  y <- lapply(patterns, function(pattern) {
    cbind(time = date_to_numeric_cycle(pattern$time, 'year', 'day'), as.matrix(pattern[rev(bands)]))
  })
  y[[2]] <- y[[2]][1:15, ]
  y[[3]][4, "NDVI"] <- NA

  distances <- t(sapply(ids, function(i) {
    sapply(y, function(pattern) {
      ts <- cbind(time = doy, x[i, , ])
      twdtw(x = ts, y = pattern[, colnames(ts)], cycle_length = 366, time_weight = tw_linear)[1]
    })
  }))

  knn <- twdtw_knn(x = x[ids, , , drop = FALSE], y = y, time = doy, cycle_length = 366, time_weight = tw_linear, k = 3)
  expect_equal(knn, k_lowest(distances, 3))
})

test_that("results do not depend on the number and the position of the time series", {
  # The time series are processed in more than one block
  n <- dim(x)[1]
  id <- rep(seq_len(n), 4)
  empty <- c(1, 5, 2000, length(id))

  x_rep <- x[id, , , drop = FALSE]
  x_rep[empty, , ] <- NA
  x_rep[5, , ] <- x[5, , ]
  x_rep[5, , "EVI"] <- NA

  knn <- twdtw_knn_call(x_rep, tw_linear, k = 2)
  ref <- twdtw_knn_call(x, tw_linear, k = 2)

  expect_equal(knn$index[-empty, ], ref$index[id[-empty], ])
  expect_equal(knn$distance[-empty, ], ref$distance[id[-empty], ])

  # Time series without complete observations
  expect_true(all(is.na(knn$index[empty, ])))
  expect_true(all(is.na(knn$distance[empty, ])))
})

test_that("patterns with the same distance are in the order of their position", {
  n <- dim(x)[1]

  knn <- twdtw_knn_call(x, tw_linear, y = patterns[c(1, 1)], k = 2)
  expect_equal(knn$index, cbind(rep(1L, n), rep(2L, n)))
  expect_equal(knn$distance[, 1], knn$distance[, 2])
})

test_that("a time-weight function that returns missing values is an error", {
  expect_error(twdtw_knn_call(x[ids, , , drop = FALSE], function(dist, el) dist + NA_real_),
               "returned missing values")
})
