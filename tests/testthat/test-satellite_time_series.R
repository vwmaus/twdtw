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

  for (i in 2:length(results)) {
    expect_equal(results[[1]], results[[i]], tolerance = 1e-5)
  }
})
