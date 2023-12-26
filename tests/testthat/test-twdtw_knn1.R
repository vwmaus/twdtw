ts_x <- readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))
ts_y <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))
expected_vector <- readRDS(system.file("sits/modis_mod13q1_results.rds", package = "twdtw"))

test_that("Nearest neighbor calculation is correct", {
  # Compute TWDTW distances
  distances <- sapply(ts_y[[2]], function(pattern) {
    sapply(ts_x[[3]], function(ts) {
      twdtw(x = as.data.frame(ts), y = as.data.frame(pattern), cycle_length = 'year', time_scale = 'day', time_weight = c(steepness = 0.1, midpoint = 50))
    })
  })

  # Find the nearest neighbor for each observation in newdata
  nearest_neighbor <- apply(distances, 1, which.min)

  # Test if nearest_neighbor is equal to the expected vector
  expect_equal(nearest_neighbor, expected_vector)
})


