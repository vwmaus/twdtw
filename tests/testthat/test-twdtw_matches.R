# A pattern that is much shorter than the time series can match more often than
# it has observations: the number of matches is bounded by the length of the
# time series (one per starting point), not by the length of the pattern.
per_year <- 23
years <- 20
doy <- rep(round(seq(1, 353, length.out = per_year)), years)
season <- sin(seq(0, 2 * pi, length.out = per_year))

# 20 seasonal cycles with a small deterministic disturbance
ts_long <- cbind(time = doy, v1 = rep(season, years) + 0.05 * cos(0.7 * seq_along(doy)))

# 5 observations from one cycle
ts_short <- cbind(time = doy[6:10], v1 = season[6:10])

twdtw_call <- function(...) {
  twdtw(x = ts_long,
        y = ts_short,
        cycle_length = 366,
        time_weight = c(steepness = 0.1, midpoint = 50), ...)
}

test_that("all matches are returned when they outnumber the pattern observations", {
  n <- nrow(ts_short)

  for (version in c('f90', 'f90goto', 'cpp')) {
    internals <- twdtw_call(output = 'internals', version = version)

    # Last row of the cost and starting point matrices, one column per observation in x
    cost <- internals$CM[n + 1, ]
    start <- internals$VM[n + 1, ]

    # One match for each run of columns sharing the same starting point
    n_matches <- length(rle(start[start != 0])$lengths)
    expect_gt(n_matches, n)
    expect_equal(nrow(internals$matches), n_matches)

    # The lowest distance is the minimum of the last row of the cost matrix
    expect_equal(min(internals$matches[, 3]), min(cost))
    expect_equal(unclass(twdtw_call(output = 'distance', version = version)), min(cost))
    expect_equal(nrow(twdtw_call(output = 'matches', version = version)), n_matches)
  }
})
