# twdtw removes the observations with missing values and compares the remaining
# observations. If every observation of a time series has a missing value, none
# is left to compare: the distance is missing and there are no matches.

n <- 23
t <- seq(0, pi, length.out = n)
d <- seq(as.Date('2020-09-01'), length.out = n, by = "15 day")

ts_y <- data.frame(time = d, v1 = sin(t) * 2, v2 = cos(t) * 3)

# Time series shifted by 30 days
ts_x <- data.frame(time = d + 30, v1 = sin(t) * 2 + 0.5, v2 = cos(t) * 3 - 0.5)

# All values of one dimension are missing, so no observation is complete
ts_na <- ts_x
ts_na$v1 <- NA

# Time series without rows
ts_empty <- ts_x[0, ]

twdtw_call <- function(x, y = ts_y, time_weight = c(steepness = 0.1, midpoint = 50), ...) {
  twdtw(x = x,
        y = y,
        cycle_length = 'year',
        time_scale = 'day',
        time_weight = time_weight, ...)
}

test_that("the distance is missing if a time series has no complete observation", {
  for (version in c('f90', 'f90goto', 'cpp')) {
    # No warning
    expect_warning(distance <- twdtw_call(ts_na, version = version), NA)
    expect_s3_class(distance, "twdtw")
    expect_true(is.na(distance))

    expect_true(is.na(twdtw_call(ts_empty, version = version)))

    # The same if the second time series has no complete observation
    expect_true(is.na(twdtw_call(ts_x, y = ts_na, version = version)))
    expect_true(is.na(twdtw_call(ts_na, y = ts_na, version = version)))
  }

  # Time-weight function
  expect_true(is.na(twdtw_call(ts_na, time_weight = function(dist, el) dist + 0.1 * el)))

  # Matrix input
  x <- cbind(date_to_numeric_cycle(ts_na$time, 'year', 'day'), as.matrix(ts_na[c("v1", "v2")]))
  y <- cbind(date_to_numeric_cycle(ts_y$time, 'year', 'day'), as.matrix(ts_y[c("v1", "v2")]))
  expect_true(is.na(twdtw(x = x, y = y, cycle_length = 366, time_weight = c(0.1, 50))))
})

test_that("there are no matches if a time series has no complete observation", {
  matches <- twdtw_call(ts_na, output = 'matches')
  expect_s3_class(matches, "twdtw")
  expect_equal(dim(matches), c(0, 3))

  internals <- twdtw_call(ts_na, output = 'internals')
  expect_s3_class(internals, "twdtw")
  expect_equal(internals$M, 0L)
  expect_equal(internals$N, n)
  expect_equal(dim(internals$CM), c(n + 1, 0))
  expect_equal(dim(internals$matches), c(0, 3))
})

test_that("a time series with some missing values is compared with its complete observations", {
  missing <- c(2, 7, 8, 15, 20)
  ts_some <- ts_x
  ts_some$v1[missing] <- NA

  for (version in c('f90', 'f90goto', 'cpp')) {
    distance <- twdtw_call(ts_some, version = version)
    expect_false(is.na(distance))
    expect_equal(distance[1], twdtw_call(ts_x[-missing, ], version = version)[1])
  }
})

test_that("a time series with a single complete observation is compared", {
  ts_one <- ts_x
  ts_one$v1[-5] <- NA

  distance <- twdtw_call(ts_one)
  expect_false(is.na(distance))
  expect_equal(distance[1], twdtw_call(ts_x[5, ])[1])
})

test_that("the distance registered in proxy is missing if a time series has no complete observation", {
  distance <- proxy::dist(x = ts_na, y = ts_y, method = "twdtw",
                          cycle_length = 'year', time_scale = 'day',
                          time_weight = c(steepness = 0.1, midpoint = 50))
  expect_true(is.na(distance[1]))
})
