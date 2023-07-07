library(twdtw)

ts_x <- read.csv(system.file("sits/timeseries.csv", package = "twdtw"))
ts_y <- read.csv(system.file("sits/reference.csv", package = "twdtw"))

# Compute TWDTW distance
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


