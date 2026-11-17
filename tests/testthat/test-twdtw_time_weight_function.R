# A time-weight function defined in R is only called by the 'f90' version of twdtw.

n <- 23
t <- seq(0, pi, length.out = n)
d <- seq(as.Date('2020-09-01'), length.out = n, by = "16 day")

ts_y <- data.frame(time = d, v1 = sin(t) * 2, v2 = cos(t) * 3)
ts_x <- data.frame(time = d + 30, v1 = sin(t) * 2 + 0.3, v2 = cos(t) * 3 - 0.2)

# Linear time weight
tw_linear <- function(dist, el) dist + 0.1 * el

twdtw_call <- function(time_weight, ...) {
  twdtw(x = ts_x,
        y = ts_y,
        cycle_length = 'year',
        time_scale = 'day',
        time_weight = time_weight, ...)
}

test_that("a time-weight function is used by the default version", {
  distance <- twdtw_call(tw_linear)

  # The result of the function is not the result of a logistic time weight
  expect_false(isTRUE(all.equal(distance[1], twdtw_call(c(0, 0))[1])))
  expect_false(isTRUE(all.equal(distance[1], twdtw_call(c(0.1, 50))[1])))

  expect_equal(twdtw_call(tw_linear, version = 'f90')[1], distance[1])
})

test_that("a time-weight function with a version that cannot call it is an error", {
  expected <- "A time-weight function can only be used with version = 'f90'"

  expect_error(twdtw_call(tw_linear, version = 'f90goto'), expected, fixed = TRUE)
  expect_error(twdtw_call(tw_linear, version = 'cpp'), expected, fixed = TRUE)

  # Matrix input
  x <- cbind(date_to_numeric_cycle(ts_x$time, 'year', 'day'), as.matrix(ts_x[c("v1", "v2")]))
  y <- cbind(date_to_numeric_cycle(ts_y$time, 'year', 'day'), as.matrix(ts_y[c("v1", "v2")]))
  expect_error(twdtw(x = x, y = y, cycle_length = 366, time_weight = tw_linear, version = 'cpp'),
               expected, fixed = TRUE)

  # The versions accept the parameters of the logistic time weight
  expect_equal(twdtw_call(c(0.1, 50), version = 'f90goto')[1], twdtw_call(c(0.1, 50))[1], tolerance = 1e-5)
  expect_equal(twdtw_call(c(0.1, 50), version = 'cpp')[1], twdtw_call(c(0.1, 50))[1], tolerance = 1e-5)
})

test_that("a logistic weight with a midpoint above any elapsed time gives the distance without time weight", {
  no_weight <- c(steepness = 1, midpoint = 1e6)

  # The same distance as a function that returns the distance of the observations
  distance <- twdtw_call(no_weight)
  expect_identical(distance[1], twdtw_call(function(dist, el) dist)[1])
  expect_lt(distance[1], twdtw_call(c(0.1, 50))[1])

  # A steepness of zero is not a weight of zero: it adds 0.5 for each pair of observations
  expect_gt(twdtw_call(c(steepness = 0, midpoint = 50))[1], distance[1])

  # Dynamic time warping in which all observations of y are matched, and the match can
  # start and end at any observation of x
  open_dtw <- function(x, y) {
    local <- as.matrix(proxy::dist(as.matrix(y[-1]), as.matrix(x[-1])))
    cost <- matrix(0, nrow(local) + 1, ncol(local))
    for (i in seq_len(nrow(local))) cost[i + 1, 1] <- cost[i, 1] + local[i, 1]
    for (j in 2:ncol(local)) {
      for (i in seq_len(nrow(local))) {
        cost[i + 1, j] <- local[i, j] + min(cost[i, j - 1], cost[i + 1, j - 1], cost[i, j])
      }
    }
    min(cost[nrow(cost), ])
  }
  expect_equal(distance[1], open_dtw(ts_x, ts_y))

  # A pattern with a part of the observations
  pattern <- ts_y[8:15, ]
  expect_equal(twdtw(ts_x, pattern, cycle_length = 'year', time_scale = 'day', time_weight = no_weight)[1],
               open_dtw(ts_x, pattern))

  # The same in a model
  m <- twdtw_knn(list(ts_y, pattern), c("a", "b"), cycle_length = 'year', time_scale = 'day', time_weight = no_weight)
  nearest <- predict(m, newdata = list(ts_x), type = "neighbors", k = 2)
  expect_equal(c(nearest$distance_1, nearest$distance_2), sort(c(open_dtw(ts_x, ts_y), open_dtw(ts_x, pattern))))
})
