ts_x <- readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))
ts_y <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))
expected_vector <- readRDS(system.file("sits/modis_mod13q1_results.rds", package = "twdtw"))

bands <- c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")
time <- ts_x$observations[[1]]$time
patterns <- ts_y$observations
tw <- c(steepness = 0.1, midpoint = 50)

# Array of time series with dimensions: series, time, and band
x <- sapply(ts_x$observations, function(ts) as.matrix(ts[bands]), simplify = "array")
x <- aperm(x, c(3, 1, 2))

# Time series with missing values and a sample of the remaining time series
has_na <- which(sapply(ts_x$observations, anyNA))
ids <- sort(unique(c(has_na, seq(1, dim(x)[1], by = 10))))

# Define TWDTW k-nearest neighbours call
twdtw_knn_call <- function(x, y = patterns, ...) {
  twdtw_knn(x = x,
            y = y,
            time = time,
            cycle_length = 'year',
            time_scale = 'day',
            time_weight = tw, ...)
}

# Distance between each pair of time series and pattern calling twdtw for each pair
twdtw_pairs <- function(ids, ...) {
  t(sapply(ids, function(i) {
    sapply(patterns, function(pattern) {
      twdtw(x = as.data.frame(ts_x$observations[[i]]),
            y = as.data.frame(pattern),
            cycle_length = 'year',
            time_scale = 'day',
            time_weight = tw, ...)[1]
    })
  }))
}

# The k lowest distances in each row and the respective columns
k_lowest <- function(distances, k) {
  list(index = t(apply(distances, 1, order))[, seq_len(k), drop = FALSE],
       distance = t(apply(distances, 1, sort))[, seq_len(k), drop = FALSE])
}

test_that("Nearest neighbor calculation is correct", {
  knn <- twdtw_knn_call(x)

  expect_equal(dim(knn$index), c(dim(x)[1], 1))
  expect_equal(dim(knn$distance), c(dim(x)[1], 1))
  expect_equal(knn$index[, 1], expected_vector)
})

test_that("k-nearest neighbours are the same as calling twdtw for each pair", {
  distances <- twdtw_pairs(ids)

  for (k in c(1, 2, 3, 5)) {
    knn <- twdtw_knn_call(x[ids, , , drop = FALSE], k = k)
    expect_equal(knn, k_lowest(distances, k))
  }
})

test_that("k-nearest neighbours are the same as calling twdtw for each pair using max_elapsed", {
  for (max_elapsed in c(30, 90)) {
    distances <- suppressWarnings(twdtw_pairs(ids, max_elapsed = max_elapsed))

    for (k in c(1, 5)) {
      knn <- twdtw_knn_call(x[ids, , , drop = FALSE], k = k, max_elapsed = max_elapsed)
      expect_equal(knn, k_lowest(distances, k))
    }
  }
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
      twdtw(x = ts, y = pattern[, colnames(ts)], cycle_length = 366, time_weight = tw)[1]
    })
  }))

  knn <- twdtw_knn(x = x[ids, , , drop = FALSE], y = y, time = doy, cycle_length = 366, time_weight = tw, k = 3)
  expect_equal(knn, k_lowest(distances, 3))

  # Time index given by its position
  knn <- twdtw_knn(x = x[ids, , , drop = FALSE], y = y, time = doy, cycle_length = 366, time_weight = tw, k = 3, index_column = 1)
  expect_equal(knn, k_lowest(distances, 3))
})

test_that("results do not depend on the number and the position of the time series", {
  n <- dim(x)[1]
  id <- rep(seq_len(n), 4)
  empty <- c(1, 5, 2000, length(id))

  x_rep <- x[id, , , drop = FALSE]
  x_rep[empty, , ] <- NA
  x_rep[5, , ] <- x[5, , ]
  x_rep[5, , "EVI"] <- NA

  knn <- twdtw_knn_call(x_rep, k = 2)
  ref <- twdtw_knn_call(x, k = 2)

  expect_equal(knn$index[-empty, ], ref$index[id[-empty], ])
  expect_equal(knn$distance[-empty, ], ref$distance[id[-empty], ])

  # Time series without complete observations
  expect_true(all(is.na(knn$index[empty, ])))
  expect_true(all(is.na(knn$distance[empty, ])))
})

test_that("patterns with the same distance are in the order of their position", {
  n <- dim(x)[1]

  knn <- twdtw_knn_call(x, y = patterns[c(1, 1)], k = 2)
  expect_equal(knn$index, cbind(rep(1L, n), rep(2L, n)))
  expect_equal(knn$distance[, 1], knn$distance[, 2])

  knn <- twdtw_knn_call(x, y = patterns[c(1, 1)])
  expect_equal(knn$index[, 1], rep(1L, n))
})

test_that("time-weight function returns the same value as the logistic time-weight", {
  tw_r_fun <- function(x, y) x + 1.0 / (1.0 + exp(-0.1 * (y - 50)))

  x_small <- x[c(has_na[1:2], 1:4), , , drop = FALSE]
  x_small[6, , ] <- NA

  knn <- twdtw_knn_call(x_small, k = 2)
  knn_fun <- twdtw_knn(x = x_small,
                       y = patterns,
                       time = time,
                       cycle_length = 'year',
                       time_scale = 'day',
                       time_weight = tw_r_fun,
                       k = 2)

  expect_equal(knn_fun, knn)
  expect_true(all(is.na(knn$index[6, ])))
})

test_that("row names are the names of the time series", {
  x_named <- x[1:3, , , drop = FALSE]
  dimnames(x_named)[[1]] <- c("a", "b", "c")

  knn <- twdtw_knn_call(x_named)
  expect_equal(rownames(knn$index), c("a", "b", "c"))
  expect_equal(rownames(knn$distance), c("a", "b", "c"))
})

test_that("incorrect inputs return an error", {
  expect_error(twdtw_knn_call(x[1, , ]), "three dimensions")
  expect_error(twdtw_knn_call(x, k = 0), "'k' must be")
  expect_error(twdtw_knn_call(x, k = length(patterns) + 1), "'k' must be")
  expect_error(twdtw_knn_call(x, index_column = "date"), "'index_column' was not found")
  expect_error(twdtw_knn_call(x, y = lapply(patterns, function(p) p[c("time", "EVI")])), "do not match")

  expect_error(twdtw_knn(x = x, y = patterns, time = time[-1], cycle_length = 'year',
                         time_scale = 'day', time_weight = tw), "length of 'time'")

  empty_pattern <- patterns
  empty_pattern[[2]]$EVI <- NA
  expect_error(twdtw_knn_call(x, y = empty_pattern), "no complete observations")
})
