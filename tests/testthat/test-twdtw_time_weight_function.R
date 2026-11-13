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
