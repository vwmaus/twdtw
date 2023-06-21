library(twdtw)

n <- 200
t <- seq(0, 4*pi, length.out = n)

x <- data.frame(date = seq(as.Date("2020-01-01"), by = "day", length.out = n),
             value = sin(t)*2 + runif(n))

y <- data.frame(date = seq(as.Date("2020-01-15"), by = "day", length.out = n),
             value = sin(t)*2 + runif(n))

plot(x, type = "l")
lines(y, col = "red")

rbenchmark::benchmark(
  f90     = twdtw(x, y, tw = c(-.1, 50), version = 'f90'),
  f90goto = twdtw(x, y, tw = c(-.1, 50), version = 'f90goto'),
  cpp     = twdtw(x, y, tw = c(-.1, 50), version = 'cpp'),
  replications = 100
)

# Check results
f90     = twdtw(x, y, tw = c(-.1, 50), version = 'f90')
f90goto = twdtw(x, y, tw = c(-.1, 50), version = 'f90goto')
cpp     = twdtw(x, y, tw = c(-.1, 50), version = 'cpp')

identical(f90, f90goto)
identical(f90, cpp)
