library(twdtw)

n <- 23
t <- seq(0, pi, length.out = n)
d <- seq(as.Date('2020-09-01'), length.out = n, by = "15 day")

x <- data.frame(time = d,      v1 = sin(t)*2 + runif(n))

# shift time by 30 days
y <- data.frame(time = d + 30, v1 = sin(t)*2 + runif(n))

plot(x, type = "l", xlim = range(c(d, d + 5)))
lines(y, col = "red")

# Test call output
twdtw(x, y,
      cycle_length = 'year',
      time_scale = 'day',
      time_weight = c(steepness = 0.1, midpoint = 50),
      output = 'distance')

twdtw(x, y,
      cycle_length = 'year',
      time_scale = 'day',
      time_weight = c(steepness = 0.1, midpoint = 50),
      output = 'matches')

twdtw(x, y,
      cycle_length = 'year',
      time_scale = 'day',
      time_weight = c(steepness = 0.1, midpoint = 50),
      output = 'internals')
