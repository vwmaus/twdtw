# Reference calculation of the TWDTW distance with a maximum elapsed time, without a time weight.
# All observations of 'y' are matched and the match can start and end at any observation of 'x'.
# Pairs of observations with an elapsed time greater than 'max_elapsed' cannot be matched
reference_distance <- function(x, y, cycle_length, max_elapsed) {
  n <- nrow(y)
  m <- nrow(x)
  elapsed <- abs(outer(y[, 1], x[, 1], `-`))
  elapsed <- pmin(elapsed, cycle_length - elapsed)
  cost <- abs(outer(y[, 2], x[, 2], `-`))
  cost[elapsed > max_elapsed] <- Inf
  cm <- matrix(0, n + 1, m)
  for (j in seq_len(m)) {
    for (i in seq_len(n)) {
      previous <- if (j == 1) cm[i, j] else min(cm[i, j - 1], cm[i + 1, j - 1], cm[i, j])
      cm[i + 1, j] <- cost[i, j] + previous
    }
  }
  min(cm[n + 1, ])
}

no_weight <- c(steepness = 1, midpoint = 1e12)
no_weight_fun <- function(dist, el) dist

x <- cbind(time = 1:5, value = c(3, 0, 10, 0, 10))
y <- cbind(time = 1:5, value = rep(3, 5))

twdtw_versions <- function(x, y, ...) {
  suppressWarnings(c(
    f90 = twdtw(x, y, time_weight = no_weight, version = 'f90', ...)[1],
    f90_fun = twdtw(x, y, time_weight = no_weight_fun, version = 'f90', ...)[1],
    f90goto = twdtw(x, y, time_weight = no_weight, version = 'f90goto', ...)[1],
    cpp = twdtw(x, y, time_weight = no_weight, version = 'cpp', ...)[1]))
}

test_that("max_elapsed constrains the first observation of x", {
  # Without a limit all observations of y match the first observation of x
  expect_equal(unname(twdtw_versions(x, y, cycle_length = 20)), rep(0, 4))

  # With max_elapsed = 0 each observation of y can only match the observation of x with the same time
  expect_equal(unname(twdtw_versions(x, y, cycle_length = 20, max_elapsed = 0)), rep(sum(abs(y[, 2] - x[, 2])), 4))

  expect_equal(unname(twdtw_versions(x, y, cycle_length = 20, max_elapsed = 1)),
               rep(reference_distance(x, y, 20, 1), 4))
})

test_that("twdtw with max_elapsed is the same as the reference calculation", {
  set.seed(11)
  for (i in 1:50) {
    m <- sample(2:15, 1)
    n <- sample(1:10, 1)
    ts_x <- cbind(time = sort(runif(m, 0, 100)), value = rnorm(m))
    ts_y <- cbind(time = sort(runif(n, 0, 100)), value = rnorm(n))
    max_elapsed <- runif(1, 5, 60)
    expected <- reference_distance(ts_x, ts_y, 100, max_elapsed)
    expect_equal(unname(twdtw_versions(ts_x, ts_y, cycle_length = 100, max_elapsed = max_elapsed)), rep(expected, 4))

    # The k-nearest neighbours have the same distance
    for (tw in list(no_weight, no_weight_fun)) {
      knn <- twdtw_nearest(x = array(ts_x[, 2], c(1, m, 1)), y = list(ts_y), time = ts_x[, 1], index_column = 1,
                           cycle_length = 100, time_weight = tw, max_elapsed = max_elapsed)
      expect_equal(knn$distance[1, 1], expected)
    }
  }
})

test_that("time series that cannot be matched within max_elapsed have an infinite distance and no matches", {
  # All observations of y are more than two time units away from the observations of x
  y_far <- cbind(time = 11:15, value = rep(3, 5))

  expect_no_warning(distance <- twdtw(x, y_far, time_weight = no_weight, cycle_length = 40, max_elapsed = 2))
  expect_equal(distance[1], Inf)
  expect_equal(unname(twdtw_versions(x, y_far, cycle_length = 40, max_elapsed = 2)), rep(Inf, 4))

  matches <- twdtw(x, y_far, time_weight = no_weight, cycle_length = 40, max_elapsed = 2, output = 'matches')
  expect_equal(dim(matches), c(0, 3))

  internals <- twdtw(x, y_far, time_weight = no_weight, cycle_length = 40, max_elapsed = 2, output = 'internals')
  expect_equal(dim(internals$matches), c(0, 3))
  expect_true(all(is.infinite(internals$CM[-1, ])))

  # The pattern that cannot be matched is the farthest neighbour
  for (tw in list(no_weight, no_weight_fun)) {
    knn <- twdtw_nearest(x = array(x[, 2], c(1, 5, 1)), y = list(y_far, y), time = x[, 1], index_column = 1,
                         cycle_length = 40, time_weight = tw, max_elapsed = 2, k = 2)
    expect_equal(knn$index[1, ], c(2L, 1L))
    expect_equal(knn$distance[1, ], c(reference_distance(x, y, 40, 2), Inf))
  }
})
