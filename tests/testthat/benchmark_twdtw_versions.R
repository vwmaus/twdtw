library(twdtw)
library(rbenchmark)

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
                       time_weight = c(steepness = 0.1, midpoint = 50), ...){
  twdtw(x = x,
        y = y,
        cycle_length = cycle_length,
        time_scale = time_scale,
        time_weight = time_weight, ...)
}

tw_r_fun <- function(x,y) x + 1.0 / (1.0 + exp(-0.1 * (y - 50)))

# Benchmark twdtw versions
# For the same task the dtw package has an elapsed time between 0.8 - 1.0 in comparison
benchmark(
  twdtw_f90        = twdtw_call(version = 'f90'),
  twdtw_f90_fun    = twdtw_call(version = 'f90', time_weight = tw_r_fun),
  twdtw_f90goto    = twdtw_call(version = 'f90goto'),
  twdtw_cpp        = twdtw_call(version = 'cpp'),
  twdtw_f90_lb     = twdtw_call(version = 'f90', max_elapsed = 30),
  twdtw_f90_fun_lb = twdtw_call(version = 'f90', max_elapsed = 30, time_weight = tw_r_fun),
  twdtw_f90goto_lb = twdtw_call(version = 'f90goto', max_elapsed = 30),
  twdtw_cpp_lb     = twdtw_call(version = 'cpp', max_elapsed = 30)
)

