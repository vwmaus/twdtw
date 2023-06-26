library(twdtw)
library(dtw)

n <- 500
t <- seq(0, 4*pi, length.out = n)

x <- data.frame(date = seq(as.Date("2020-01-01"), by = "day", length.out = n),
             value = sin(t)*2 + runif(n))

y <- data.frame(date = seq(as.Date("2020-01-15"), by = "day", length.out = n),
             value = sin(t)*2 + runif(n))

plot(x, type = "l")
lines(y, col = "red")

rbenchmark::benchmark(
  twdtw_f90     = twdtw(x, y, tw = c(-.1, 50), version = 'f90'),
  twdtw_f90goto = twdtw(x, y, tw = c(-.1, 50), version = 'f90goto'),
  twdtw_cpp     = twdtw(x, y, tw = c(-.1, 50), version = 'cpp'),
  dtw_cpp       = dtw(x$value, y$value, distance.only = TRUE), # does not support time dimension
  replications = 100
)

# Check results
f90     = twdtw(x, y, tw = c(-.1, 50), version = 'f90')
f90goto = twdtw(x, y, tw = c(-.1, 50), version = 'f90goto')
cpp     = twdtw(x, y, tw = c(-.1, 50), version = 'cpp')

identical(f90, f90goto)
identical(f90, cpp)


rbenchmark::benchmark(
  #twdtw_f90     = twdtw(x, y, tw = c(-.1, 50, 50),
  #                      tw_r = function(dist,td,tw1,tw2) dist + 1.0 / (1.0 + exp(-0.1 * (td - 50)))),
  twdtw_f90_tw  = twdtw(x, y, tw = c(-0.1, 50, 50), NPAR = as.integer(3)),
  twdtw_f90goto = twdtw(x, y, tw = c(-.1, 50), version = 'f90goto'),
  twdtw_cpp     = twdtw(x, y, tw = c(-.1, 50), version = 'cpp'),
  dtw_cpp       = dtw(x$value, y$value, distance.only = TRUE), # does not support time dimension
  replications = 100
)
