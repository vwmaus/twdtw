library(twdtw)
library(dtw)
library(ggplot2)
library(rbenchmark)

n <- 23
t <- seq(0, pi, length.out = n)
d <- seq(as.Date('2020-09-01'), length.out = n, by = "15 day")

# Define funcrion to normalize time series
fun_range <- function(x) {
  (x - min(x)) / (max(x) - min(x))
}

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
                       time_cycle_length = 'year',
                       time_cycle_scale = 'day',
                       time_weight_par = c(steepness = 0.1, midpoint = 50), ...){
  twdtw(x = x,
        y = y,
        time_cycle_length = time_cycle_length,
        time_cycle_scale = time_cycle_scale,
        time_weight = time_weight_par, ...)
}

tw_r_fun <- function(x,y) x + 1.0 / (1.0 + exp(-0.1 * (y - 50)))

# Benchmark default TWDTW call
benchmark(
  twdtw_f90        = twdtw_call(version = 'f90'),
  twdtw_f90_fun    = twdtw_call(version = 'f90', time_weight = tw_r_fun),
  twdtw_f90goto    = twdtw_call(version = 'f90goto'),
  twdtw_cpp        = twdtw_call(version = 'cpp'),
  twdtw_f90_lb     = twdtw_call(version = 'f90', max_elapsed = 30),
  twdtw_f90_fun_lb = twdtw_call(version = 'f90', max_elapsed = 30, time_weight = tw_r_fun),
  twdtw_f90goto_lb = twdtw_call(version = 'f90goto', max_elapsed = 30),
  twdtw_cpp_lb     = twdtw_call(version = 'cpp', max_elapsed = 30),
  dtw_cpp          = dtw(ts_x[,c(2,2:4)], ts_y[,c(2,2:4)], distance.only = TRUE), # does not support time dimension
  replications = 1000
)

# All calls must return the same result
twdtw_f90     = twdtw_call(version = 'f90')
twdtw_f90_fun = twdtw_call(version = 'f90', time_weight = tw_r_fun)
twdtw_f90goto = twdtw_call(version = 'f90goto')
twdtw_cpp     = twdtw_call(version = 'cpp')
twdtw_f90_lb     = twdtw_call(version = 'f90', max_elapsed = 30)
twdtw_f90_fun_lb = twdtw_call(version = 'f90', max_elapsed = 30, time_weight = tw_r_fun)
twdtw_f90goto_lb = twdtw_call(version = 'f90goto', max_elapsed = 30)
twdtw_cpp_lb     = twdtw_call(version = 'cpp', max_elapsed = 30)

# Check results
identical(twdtw_f90, twdtw_f90_fun)
identical(twdtw_f90, twdtw_f90goto)
identical(twdtw_f90, twdtw_cpp)
identical(twdtw_f90, twdtw_f90_lb)
identical(twdtw_f90, twdtw_f90_fun_lb)
identical(twdtw_f90, twdtw_f90goto_lb)
identical(twdtw_f90, twdtw_cpp_lb)
