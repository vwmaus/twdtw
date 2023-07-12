# Your setup
date_sys <- "2023-07-12 09:59:25 EEST"

cycle_options <- c('year', 'month', 'day', 'hour', 'minute', 'second')
time_options <-  c('year', 'month', 'day', 'hour', 'minute', 'second')

# Expected results matrix
expected_results <- matrix(
  list(
    NULL, NULL, NULL, NULL, NULL, NULL,
    c(7), NULL, NULL, NULL, NULL, NULL,
    c(193), c(12), NULL, NULL, NULL, NULL,
    c(4641), c(297), c(9), NULL, NULL, NULL,
    c(278519), c(17879), c(599), c(59), NULL, NULL,
    c(16711165), c(1072765), c(35965), c(3565), c(25), NULL
  ),
  nrow = length(time_options),
  ncol = length(cycle_options),
  dimnames = list(time_options, cycle_options)
)

# Testing
test_that("sapply with date_to_numeric_cycle returns expected results", {
  actual_results <- sapply(time_options, function(y) {
    sapply(cycle_options, function(z) {
      date_to_numeric_cycle(date_sys, z, y)
    })
  })
  expect_equal(actual_results, expected_results)
})
