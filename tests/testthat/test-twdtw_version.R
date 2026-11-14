# The versions 'f90goto' and 'cpp' of the twdtw routine are deprecated:
# twdtw gives a warning once per session if one of them is used.

n <- 23
t <- seq(0, pi, length.out = n)
d <- seq(as.Date('2020-09-01'), length.out = n, by = "16 day")

ts_y <- data.frame(time = d, v1 = sin(t) * 2, v2 = cos(t) * 3)
ts_x <- data.frame(time = d + 30, v1 = sin(t) * 2 + 0.3, v2 = cos(t) * 3 - 0.2)

twdtw_call <- function(...) {
  twdtw(x = ts_x,
        y = ts_y,
        cycle_length = 'year',
        time_scale = 'day',
        time_weight = c(steepness = 0.1, midpoint = 50), ...)
}

# Set whether the warning was already given in this session
set_warned <- function(warned) {
  assign("deprecated_version_warned", warned, envir = twdtw:::.twdtw_session)
}

test_that("a deprecated version gives a warning once per session", {
  expected <- "The versions 'f90goto' and 'cpp' of twdtw are deprecated and will be removed in version 1.2-0"

  for (version in c('cpp', 'f90goto')) {
    set_warned(FALSE)

    expect_warning(distance <- twdtw_call(version = version), expected, fixed = TRUE)

    # The same result as the default version
    expect_equal(distance[1], twdtw_call()[1])

    # No warning for the next calls, with any version
    expect_warning(twdtw_call(version = 'cpp'), NA)
    expect_warning(twdtw_call(version = 'f90goto'), NA)
  }

  # The other tests do not show the warning
  set_warned(TRUE)
})

test_that("the default version gives no warning", {
  set_warned(FALSE)

  expect_warning(twdtw_call(), NA)
  expect_warning(twdtw_call(version = 'f90'), NA)

  # The warning is still to be given
  expect_warning(twdtw_call(version = 'cpp'), "deprecated")

  set_warned(TRUE)
})

test_that("an unknown version is an error", {
  expected <- "Incorrect version. Must be one of 'f90', 'f90goto', 'cpp'"

  expect_error(twdtw_call(version = 'f77'), expected, fixed = TRUE)
  expect_error(twdtw_call(version = c('f90', 'cpp')), expected, fixed = TRUE)
  expect_error(twdtw_call(version = 1), expected, fixed = TRUE)
  expect_error(twdtw_call(version = NULL), expected, fixed = TRUE)
})
