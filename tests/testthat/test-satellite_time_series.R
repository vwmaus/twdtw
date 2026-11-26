ts_x <- read.csv(system.file("sits/timeseries.csv", package = "twdtw"))
ts_y <- read.csv(system.file("sits/reference.csv", package = "twdtw"))

# Compute TWDTW distance - result should be
twdtw(x = ts_x,
      y = ts_y,
      cycle_length = 'year',
      time_scale = 'day',
      time_weight = c(steepness = 0.1, midpoint = 50), output = 'distance')

twdtw(x = ts_x,
      y = ts_y,
      cycle_length = 'year',
      time_scale = 'day',
      time_weight = c(steepness = 0.1, midpoint = 50), output = 'matches')

# Check matrices
m <- twdtw(x = ts_x,
           y = ts_y,
           cycle_length = 'year',
           time_scale = 'day',
           time_weight = c(steepness = 0.1, midpoint = 50), output = 'internals')

# Check plot cost matrix
plot_cost_matrix(m)

test_that("twdtw returns the correct value", {
  expected_value <- 5.93395
  actual_value <- twdtw(x = ts_x,
                        y = ts_y,
                        cycle_length = 'year',
                        time_scale = 'day',
                        time_weight = c(steepness = 0.1, midpoint = 50), output = 'distance')[1]
  expect_equal(actual_value, expected_value, tolerance = 0.0001)
})


test_that("multiple functions return the same value", {
  results <- list(
    twdtw_f90 = twdtw(x = ts_x,
                      y = ts_y,
                      cycle_length = 'year',
                      time_scale = 'day',
                      time_weight = c(steepness = 0.1, midpoint = 50), output = 'matches', version = 'f90'),

    twdtw_f90goto = twdtw(x = ts_x,
                          y = ts_y,
                          cycle_length = 'year',
                          time_scale = 'day',
                          time_weight = c(steepness = 0.1, midpoint = 50), output = 'matches', version = 'f90goto'),

    twdtw_cpp = twdtw(x = ts_x,
                      y = ts_y,
                      cycle_length = 'year',
                      time_scale = 'day',
                      time_weight = c(steepness = 0.1, midpoint = 50), output = 'matches', version = 'cpp'),

    twdtw_f90_me = twdtw(x = ts_x,
                         y = ts_y,
                         cycle_length = 'year',
                         time_scale = 'day',
                         time_weight = c(steepness = 0.1, midpoint = 50), output = 'matches', version = 'f90', max_elapsed = 60),

    twdtw_f90goto_me = twdtw(x = ts_x,
                             y = ts_y,
                             cycle_length = 'year',
                             time_scale = 'day',
                             time_weight = c(steepness = 0.1, midpoint = 50), output = 'matches', version = 'f90goto', max_elapsed = 60),

    twdtw_cpp_me = twdtw(x = ts_x,
                         y = ts_y,
                         cycle_length = 'year',
                         time_scale = 'day',
                         time_weight = c(steepness = 0.1, midpoint = 50), output = 'matches', version = 'cpp', max_elapsed = 60)

  )

  # The versions return the same matches, without and with a maximum elapsed time
  for (i in 2:3) {
    expect_equal(results[[1]], results[[i]], tolerance = 1e-5)
  }

  for (i in 5:6) {
    expect_equal(results[[4]], results[[i]], tolerance = 1e-5)
  }

  # Without a maximum elapsed time there is a match of all observations of the pattern to
  # the first observation of the time series. The maximum elapsed time excludes this match
  # and does not change the other matches
  all_to_first <- results[[1]][, 1] == 1 & results[[1]][, 2] == 1
  expect_equal(sum(all_to_first), 1)
  expect_equal(unclass(results[[4]]), unclass(results[[1]])[!all_to_first, , drop = FALSE], tolerance = 1e-5)
})
