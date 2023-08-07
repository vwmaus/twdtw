n <- 23
t <- seq(0, pi, length.out = n)
d <- seq(as.Date('2020-09-01'), length.out = n, by = "15 day")

# Create time series with 4 dimensions
ts_x <- data.frame(time = d,
                   v1 = sin(t)*2 + runif(n),
                   v2 = sin(t)*3 + runif(n),
                   v3 = sin(t)*4 + runif(n),
                   v4 = sin(t)*5 + runif(n))

# Create a second time series with 4 dimensions an shift time by 30 days
ts_y <- data.frame(time = d + 30,
                   v1 = sin(t)*2 + runif(n),
                   v2 = sin(t)*3 + runif(n),
                   v3 = sin(t)*4 + runif(n),
                   v4 = sin(t)*5 + runif(n))

# Define TWDTW call
twdtw_call <- function(x = ts_x,
                       y = ts_y,
                       cycle_length = 'year',
                       time_scale = 'day',
                       time_weight = c(steepness = 0.1, midpoint = 50),
                       output = 'matches', ...){
  twdtw(x = x,
        y = y,
        cycle_length = cycle_length,
        time_scale = time_scale,
        time_weight = time_weight,
        output = output, ...)
}

tw_r_fun <- function(x,y) x + 1.0 / (1.0 + exp(-0.1 * (y - 50)))

# All calls must return the same result
test_that("multiple functions return the same value", {
  results <- list(
    twdtw_f90     = twdtw_call(version = 'f90'),
    twdtw_f90_fun = twdtw_call(version = 'f90', time_weight = tw_r_fun),
    twdtw_f90goto = twdtw_call(version = 'f90goto'),
    twdtw_cpp     = twdtw_call(version = 'cpp')
  )

  for (i in 2:length(results)) {
    expect_equal(results[[1]], results[[i]], tolerance = 1e-5)
  }
})

# Check NA inputs
insert_random_NAs <- function(df, num_NAs) {
  n <- nrow(df)
  m <- ncol(df)
  row_indices <- c(1, sample(1:n, num_NAs-1, replace = TRUE))
  col_indices <- sample(1:m, num_NAs, replace = TRUE)
  for(i in 1:num_NAs){
    df[row_indices[i], col_indices[i]] <- NA
  }
  return(df)
}

ts_x_na <- insert_random_NAs(ts_x, 20)
ts_y_na <- insert_random_NAs(ts_y, 20)
test_that("multiple functions return the same value", {
  results <- list(
    twdtw_f90     = twdtw_call(ts_x_na, ts_y_na, version = 'f90'),
    twdtw_f90_fun = twdtw_call(ts_x_na, ts_y_na, version = 'f90', time_weight = tw_r_fun),
    twdtw_f90goto = twdtw_call(ts_x_na, ts_y_na, version = 'f90goto'),
    twdtw_cpp     = twdtw_call(ts_x_na, ts_y_na, version = 'cpp')
  )

  for (i in 2:length(results)) {
    expect_equal(results[[1]], results[[i]], tolerance = 1e-5)
  }
})

# Test output types
twdtw_call(output = 'distance')

twdtw_call(output = 'matches')

twdtw_call(output = 'internals')

# Test proxy call
test_that("twdtw_call returns the same result as proxy::dist with given parameters", {
  expected <- proxy::dist(x = ts_x, y = ts_y, method = "twdtw", cycle_length = 'year', time_scale = 'day', time_weight = c(steepness = 0.1, midpoint = 50))
  actual <- twdtw_call(output = 'distance')
  expect_equal(actual[1], expected[1])
})