# Example data and functions for the tests of twdtw_knn() and its predict method

# Time series of 999 pixels with 23 dates and six bands, the temporal patterns of
# five land use classes, and the position of the nearest pattern of each time series
knn_series <- lapply(readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))$observations,
                     as.data.frame)
knn_patterns <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))
knn_nearest <- readRDS(system.file("sits/modis_mod13q1_results.rds", package = "twdtw"))

knn_bands <- c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")
knn_time <- knn_series[[1]]$time

# Labelled time series: a sample of the time series with the label of their nearest pattern
knn_sample_ids <- seq(5, 999, by = 12)
knn_samples <- knn_series[knn_sample_ids]
knn_sample_labels <- knn_patterns$label[knn_nearest[knn_sample_ids]]

# Time series to classify: the time series with missing values and a sample of the others
knn_new_ids <- sort(unique(c(which(sapply(knn_series, anyNA)), seq(2, 999, by = 45))))

# Create a model with the TWDTW arguments of the examples
knn_model <- function(x, y, ...) {
  twdtw_knn(x = x,
            y = y,
            cycle_length = 'year',
            time_scale = 'day',
            time_weight = c(steepness = 0.1, midpoint = 50), ...)
}

# Array with the dimensions series, time, and band from a list of time series with the same dates
knn_array <- function(ts_list, bands = knn_bands) {
  x <- sapply(ts_list, function(ts) as.matrix(ts[bands]), simplify = "array")
  aperm(x, c(3, 1, 2))
}

# Data frame with one row for each time series and one column for each band and date
knn_wide <- function(ts_list, bands = knn_bands, time = knn_time, sep = ".") {
  x <- knn_array(ts_list, bands)
  wide <- as.data.frame(matrix(x, nrow = dim(x)[1]))
  names(wide) <- paste(rep(bands, each = length(time)), rep(format(time), length(bands)), sep = sep)
  wide
}

# TWDTW distance between each time series and each pattern, calculated with twdtw()
# one pair at a time
knn_pair_distance <- function(ts_list, patterns, ...) {
  t(sapply(ts_list, function(ts) {
    sapply(patterns, function(pattern) {
      twdtw(x = as.data.frame(ts),
            y = as.data.frame(pattern),
            cycle_length = 'year',
            time_scale = 'day',
            time_weight = c(steepness = 0.1, midpoint = 50), ...)[1]
    })
  }))
}

# Vote of the k nearest patterns of one time series, written without matrices:
# the label with the most votes and the share of the votes of each label
knn_expected_vote <- function(distance, labels, k, weighted = FALSE) {
  nearest <- order(distance)[seq_len(k)]
  weight <- if (weighted) 1 / distance[nearest] else rep(1, k)
  votes <- tapply(weight, labels[nearest], sum)
  most_voted <- names(votes)[votes == max(votes)]
  first <- nearest[labels[nearest] %in% most_voted][1]
  share <- setNames(rep(0, length(unique(labels))), sort(unique(labels)))
  share[names(votes)] <- votes / sum(votes)
  list(label = labels[first], distance = distance[first], share = share)
}
