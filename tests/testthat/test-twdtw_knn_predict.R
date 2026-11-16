# predict() classifies time series with a model created by twdtw_knn(): the k nearest
# temporal patterns, the vote of their labels, and the forms of the time series.

# Model with one pattern for each label and model that keeps the labelled time series
pattern_model <- knn_model(knn_patterns$observations, knn_patterns$label)
sample_model <- knn_model(knn_samples, knn_sample_labels)

# Time series to classify and their distance to each sample, calculated one pair at a time
new_series <- knn_series[knn_new_ids]
sample_distance <- knn_pair_distance(new_series, knn_samples)

test_that("the nearest pattern is the pattern with the lowest distance calculated with twdtw", {
  lu <- predict(pattern_model, newdata = knn_series)

  expect_s3_class(lu, "data.frame")
  expect_named(lu, c("prediction", "distance"))
  expect_equal(nrow(lu), length(knn_series))

  # The predictions are factors with the labels of the model in alphabetical order
  expect_s3_class(lu$prediction, "factor")
  expect_identical(levels(lu$prediction), sort(knn_patterns$label))
  expect_identical(as.character(lu$prediction), knn_patterns$label[knn_nearest])

  # Time series with missing values are classified with their remaining observations
  pattern_distance <- knn_pair_distance(new_series, knn_patterns$observations)
  expect_equal(lu$distance[knn_new_ids], apply(pattern_distance, 1, min))
  expect_identical(as.character(lu$prediction[knn_new_ids]),
                   knn_patterns$label[apply(pattern_distance, 1, which.min)])
  expect_false(anyNA(lu))
})

test_that("the forms of newdata give the same result", {
  lu <- predict(sample_model, newdata = new_series, k = 3)

  # An array, with and without band names, and with the bands in another order
  x <- knn_array(new_series)
  expect_equal(predict(sample_model, newdata = x, time = knn_time, k = 3), lu)
  expect_equal(predict(sample_model, newdata = unname(x), time = knn_time, k = 3), lu)
  expect_equal(predict(sample_model, newdata = x[, , 6:1], time = knn_time, k = 3), lu)
  expect_equal(predict(sample_model, newdata = x, time = as.character(knn_time), k = 3), lu)

  # An array with dates that are not in time order
  expect_equal(predict(sample_model, newdata = x[, 23:1, ], time = rev(knn_time), k = 3), lu)

  # A data.frame with one column for each band and date, and columns that are not observations
  wide <- knn_wide(new_series)
  expect_equal(predict(sample_model, newdata = wide, k = 3), lu)
  expect_equal(predict(sample_model, newdata = cbind(x = 1, y = 2, wide[rev(seq_along(wide))]), k = 3), lu)

  # A list of matrices with the time index in the last column
  ts_matrix <- lapply(new_series, function(ts) cbind(as.matrix(ts[knn_bands]), time = as.numeric(format(ts$time, "%j"))))
  expect_equal(predict(sample_model, newdata = ts_matrix, k = 3, index_column = "time"), lu)
  expect_equal(predict(sample_model, newdata = lapply(ts_matrix, unname), k = 3, index_column = 7), lu)

  # A single time series in a list
  expect_equal(predict(sample_model, newdata = new_series[4], k = 3), lu[4, ], ignore_attr = "row.names")
})

test_that("a list of time series with the same dates is converted to one array", {
  x <- knn_array(new_series)
  lu <- predict(sample_model, newdata = x, time = knn_time, k = 3, type = "neighbors")

  # Data frames with the same columns and the same dates are converted together
  converted <- twdtw:::.same_time_series_array(new_series, "time", knn_bands)
  expect_equal(converted$values, x)
  expect_identical(converted$time, knn_time)
  expect_identical(converted$rows, seq_along(new_series))
  expect_equal(predict(sample_model, newdata = new_series, k = 3, type = "neighbors"), lu)

  # Rows and bands in another order, and dates as character strings
  reordered <- lapply(new_series, function(ts) ts[23:1, c("MIR", "NIR", "time", "BLUE", "RED", "NDVI", "EVI")])
  expect_equal(twdtw:::.same_time_series_array(reordered, "time", knn_bands), converted)
  expect_equal(predict(sample_model, newdata = reordered, k = 3, type = "neighbors"), lu)
  as_text <- lapply(new_series, function(ts) {
    ts$time <- as.character(ts$time)
    ts
  })
  expect_equal(twdtw:::.same_time_series_array(as_text, "time", knn_bands), converted)

  # Matrices with the time index in the first column, with and without column names
  ts_matrix <- lapply(new_series, function(ts) cbind(time = as.numeric(format(ts$time, "%j")), as.matrix(ts[knn_bands])))
  converted_matrix <- twdtw:::.same_time_series_array(ts_matrix, "time", knn_bands)
  expect_equal(converted_matrix$values, x)
  expect_equal(converted_matrix$time, as.numeric(format(knn_time, "%j")))
  expect_equal(twdtw:::.same_time_series_array(lapply(ts_matrix, unname), 1, knn_bands), converted_matrix)
  expect_equal(predict(sample_model, newdata = ts_matrix, k = 3, type = "neighbors"), lu)

  # Time series that are not as the first one are converted one by one, with the same result
  one_differs <- new_series
  one_differs[[5]] <- one_differs[[5]][c("time", rev(knn_bands))]
  expect_null(twdtw:::.same_time_series_array(one_differs, "time", knn_bands))
  expect_equal(predict(sample_model, newdata = one_differs, k = 3, type = "neighbors"), lu)
  other_dates <- new_series
  other_dates[[2]] <- other_dates[[2]][-1, ]
  expect_null(twdtw:::.same_time_series_array(other_dates, "time", knn_bands))
  missing_date <- new_series
  missing_date[[1]]$time[3] <- NA
  expect_null(twdtw:::.same_time_series_array(missing_date, "time", knn_bands))
  expect_null(twdtw:::.same_time_series_array(c(new_series, list(ts_matrix[[1]])), "time", knn_bands))
  expect_null(twdtw:::.same_time_series_array(lapply(new_series, function(ts) ts[c("time", "EVI")]), "time", knn_bands))

  # The conversion one by one gives the error messages
  not_numeric <- new_series
  not_numeric[[3]]$EVI <- as.character(not_numeric[[3]]$EVI)
  expect_null(twdtw:::.same_time_series_array(not_numeric, "time", knn_bands))
  expect_error(predict(sample_model, newdata = not_numeric),
               "Element 3 of 'newdata' has non-numeric columns", fixed = TRUE)
  expect_error(predict(sample_model, newdata = new_series, index_column = "date"),
               "The 'index_column' was not found in element 1 of 'newdata'", fixed = TRUE)
  no_dates <- lapply(new_series[1:3], function(ts) {
    ts$time <- "first"
    ts
  })
  expect_null(twdtw:::.same_time_series_array(no_dates, "time", knn_bands))
  expect_error(predict(sample_model, newdata = no_dates))
})

test_that("time series with different dates are classified in the order of newdata", {
  # Time series without some of their observations, with rows and bands in another order
  different <- new_series[1:12]
  different[[2]] <- different[[2]][-(1:4), ]
  different[[5]] <- different[[5]][-(1:4), ]
  different[[3]] <- different[[3]][seq(1, 23, by = 2), ]
  different[[7]] <- different[[7]][23:1, c("MIR", "time", "EVI", "NDVI", "RED", "BLUE", "NIR")]
  different[[9]]$time <- different[[9]]$time + 365
  different[[11]] <- different[[11]][0, ]

  for (type in c("class", "prob", "neighbors")) {
    lu <- predict(sample_model, newdata = different, k = 3, vote = "weighted", type = type)
    one_by_one <- do.call(rbind, lapply(different, function(ts) {
      predict(sample_model, newdata = list(ts), k = 3, vote = "weighted", type = type)
    }))
    expect_equal(lu, one_by_one, ignore_attr = "row.names")
  }

  # The distances are the distances calculated with twdtw for the observations in time order
  lu <- predict(sample_model, newdata = different)
  in_time_order <- lapply(different[-11], function(ts) ts[order(ts$time), ])
  distance <- knn_pair_distance(in_time_order, knn_samples)
  expect_equal(lu$distance[-11], apply(distance, 1, min))
  expect_identical(as.character(lu$prediction[-11]), knn_sample_labels[apply(distance, 1, which.min)])

  # A time series without observations has no prediction
  expect_true(is.na(lu$prediction[11]) && is.na(lu$distance[11]))
})

test_that("the labels of the k nearest patterns are combined by vote", {
  labels <- c("a", "a", "b", "b", "c")

  # Four time series and their three nearest patterns. The third has no neighbours
  index <- rbind(c(3, 1, 2),
                 c(1, 3, 5),
                 c(NA, NA, NA),
                 c(5, 3, 4))
  distance <- rbind(c(1, 4, 4),
                    c(2, 3, 5),
                    c(NA, NA, NA),
                    c(0, 2, 2))

  # Majority: two of three, a tie that goes to the nearest pattern, and two of three
  majority <- twdtw:::.knn_result(index, distance, labels, vote = "majority", type = "class")
  expect_named(majority, c("prediction", "distance"))
  expect_s3_class(majority$prediction, "factor")
  expect_identical(levels(majority$prediction), c("a", "b", "c"))
  expect_equal(as.character(majority$prediction), c("a", "a", NA, "b"))
  expect_equal(majority$distance, c(4, 2, NA, 2))

  # Weights 1/distance: the nearest pattern outweighs two patterns that are four times
  # as far, and a pattern at distance zero gets all the votes
  weighted <- twdtw:::.knn_result(index, distance, labels, vote = "weighted", type = "class")
  expect_equal(as.character(weighted$prediction), c("b", "a", NA, "c"))
  expect_equal(weighted$distance, c(1, 2, NA, 0))

  # Share of the votes of each label
  share <- twdtw:::.knn_result(index, distance, labels, vote = "majority", type = "prob")
  expect_named(share, c("a", "b", "c"))
  expect_equal(share$a, c(2/3, 1/3, NA, 0))
  expect_equal(share$b, c(1/3, 1/3, NA, 2/3))
  expect_equal(share$c, c(0, 1/3, NA, 1/3))

  share <- twdtw:::.knn_result(index, distance, labels, vote = "weighted", type = "prob")
  expect_equal(share$a, c(1/3, (1/2) / (1/2 + 1/3 + 1/5), NA, 0))
  expect_equal(share$b, c(2/3, (1/3) / (1/2 + 1/3 + 1/5), NA, 0))
  expect_equal(share$c, c(0, (1/5) / (1/2 + 1/3 + 1/5), NA, 1))

  # No vote: the neighbours from the nearest to the farthest
  neighbors <- twdtw:::.knn_result(index, distance, labels, vote = "majority", type = "neighbors")
  expect_named(neighbors, c(paste0("prediction_", 1:3), paste0("distance_", 1:3)))
  expect_equal(as.character(neighbors$prediction_1), c("b", "a", NA, "c"))
  expect_equal(as.character(neighbors$prediction_2), c("a", "b", NA, "b"))
  expect_equal(as.character(neighbors$prediction_3), c("a", "c", NA, "b"))
  expect_identical(levels(neighbors$prediction_3), c("a", "b", "c"))
  expect_equal(neighbors$distance_2, c(4, 3, NA, 2))

  # A single neighbour gives its label with any vote
  for (vote in c("majority", "weighted")) {
    nearest <- twdtw:::.knn_result(index[, 1, drop = FALSE], distance[, 1, drop = FALSE], labels, vote = vote, type = "class")
    expect_equal(as.character(nearest$prediction), c("b", "a", NA, "c"))
    expect_equal(nearest$distance, c(1, 2, NA, 0))
    expect_identical(levels(nearest$prediction), c("a", "b", "c"))
    share <- twdtw:::.knn_result(index[, 1, drop = FALSE], distance[, 1, drop = FALSE], labels, vote = vote, type = "prob")
    expect_equal(share$b, c(1, 0, NA, 0))
  }

  # Several patterns at distance zero share the votes, and the label with the most of them is predicted.
  # Patterns at infinite distance have one vote each
  index <- rbind(c(1, 3, 4), c(2, 5, 3))
  distance <- rbind(c(0, 0, 0), c(Inf, Inf, Inf))
  weighted <- twdtw:::.knn_result(index, distance, labels, vote = "weighted", type = "class")
  expect_equal(as.character(weighted$prediction), c("b", "a"))
  expect_equal(weighted$distance, c(0, Inf))
  share <- twdtw:::.knn_result(index, distance, labels, vote = "weighted", type = "prob")
  expect_equal(share$a, c(1/3, 1/3))
  expect_equal(share$b, c(2/3, 1/3))

  # Negative distances of a time-weight function cannot be weights
  expect_error(twdtw:::.knn_result(index, -1 - distance[c(1, 1), ], labels, vote = "weighted", type = "class"),
               "vote = 'weighted' cannot be used with negative TWDTW distances", fixed = TRUE)
  expect_no_error(twdtw:::.knn_result(index, -1 - distance[c(1, 1), ], labels, vote = "majority", type = "class"))
})

test_that("k nearest samples classify time series as a vote calculated one pair at a time", {
  labels <- knn_sample_labels
  n <- length(new_series)
  k <- 7

  for (vote in c("majority", "weighted")) {
    expected <- lapply(seq_len(n), function(i) {
      knn_expected_vote(sample_distance[i, ], labels, k, weighted = vote == "weighted")
    })

    # The label with the most votes and the distance to its nearest sample
    lu <- predict(sample_model, newdata = new_series, k = k, vote = vote)
    expect_named(lu, c("prediction", "distance"))
    expect_identical(as.character(lu$prediction), vapply(expected, function(e) e$label, character(1)))
    expect_equal(lu$distance, vapply(expected, function(e) e$distance, numeric(1)))

    # The share of the votes of each label
    share <- predict(sample_model, newdata = new_series, k = k, vote = vote, type = "prob")
    expect_named(share, sort(unique(labels)))
    expect_equal(as.matrix(share), t(vapply(expected, function(e) e$share, numeric(length(unique(labels))))),
                 ignore_attr = TRUE)
    expect_equal(rowSums(share), rep(1, n), ignore_attr = TRUE)

    # The predicted label has the largest share of the votes
    predicted_share <- as.matrix(share)[cbind(seq_len(n), as.integer(lu$prediction))]
    expect_equal(predicted_share, apply(share, 1, max), ignore_attr = TRUE)
  }

  # The votes of some time series are divided between labels
  share <- predict(sample_model, newdata = new_series, k = k, type = "prob")
  expect_true(any(apply(share, 1, max) < 1))

  # No vote: the labels and the distances of the k nearest samples
  neighbors <- predict(sample_model, newdata = new_series, k = k, type = "neighbors")
  expect_named(neighbors, c(paste0("prediction_", seq_len(k)), paste0("distance_", seq_len(k))))
  for (j in seq_len(k)) {
    nearest <- apply(sample_distance, 1, function(d) order(d)[j])
    expect_identical(as.character(neighbors[[paste0("prediction_", j)]]), labels[nearest])
    expect_equal(neighbors[[paste0("distance_", j)]], sample_distance[cbind(seq_len(n), nearest)])
    expect_identical(levels(neighbors[[paste0("prediction_", j)]]), sort(unique(labels)))
  }

  # The default is the nearest sample
  lu <- predict(sample_model, newdata = new_series)
  expect_equal(predict(sample_model, newdata = new_series, k = 1), lu)
  expect_identical(as.character(lu$prediction), labels[apply(sample_distance, 1, which.min)])
  expect_equal(lu$distance, apply(sample_distance, 1, min))
  share <- predict(sample_model, newdata = new_series, type = "prob")
  expect_equal(as.matrix(share)[cbind(seq_len(n), as.integer(lu$prediction))], rep(1, n))
  expect_equal(rowSums(share), rep(1, n), ignore_attr = TRUE)
})

test_that("with one pattern for each label a vote gives the label of the nearest pattern", {
  lu <- predict(pattern_model, newdata = new_series)
  expect_equal(predict(pattern_model, newdata = new_series, k = 3), lu)
  expect_equal(predict(pattern_model, newdata = new_series, k = 3, vote = "weighted"), lu)
  expect_equal(predict(pattern_model, newdata = new_series, k = 5), lu)

  # Weighted votes of all patterns give shares that are proportional to 1/distance
  distance <- knn_pair_distance(new_series, knn_patterns$observations)
  share <- predict(pattern_model, newdata = new_series, k = 5, vote = "weighted", type = "prob")
  expected <- (1 / distance) / rowSums(1 / distance)
  expect_equal(as.matrix(share), expected[, order(knn_patterns$label)], ignore_attr = TRUE)

  # Without a vote the nearest pattern is the prediction with k = 1
  neighbors <- predict(pattern_model, newdata = new_series, k = 2, type = "neighbors")
  expect_identical(neighbors$prediction_1, lu$prediction)
  expect_equal(neighbors$distance_1, lu$distance)
  expect_true(all(neighbors$distance_2 >= neighbors$distance_1))
  expect_false(any(neighbors$prediction_2 == neighbors$prediction_1))
})

test_that("a time series without complete observations has no prediction for any k and type", {
  x <- knn_array(new_series[1:3])
  x[2, , "BLUE"] <- NA

  for (k in c(1, 5)) {
    for (vote in c("majority", "weighted")) {
      lu <- predict(sample_model, newdata = x, time = knn_time, k = k, vote = vote)
      expect_equal(is.na(lu$prediction), c(FALSE, TRUE, FALSE))
      expect_equal(is.na(lu$distance), c(FALSE, TRUE, FALSE))
      share <- predict(sample_model, newdata = x, time = knn_time, k = k, vote = vote, type = "prob")
      expect_equal(rowSums(is.na(share)), c(0, ncol(share), 0), ignore_attr = TRUE)
    }
    neighbors <- predict(sample_model, newdata = x, time = knn_time, k = k, type = "neighbors")
    expect_equal(rowSums(is.na(neighbors)), c(0, 2 * k, 0), ignore_attr = TRUE)
  }

  # The same time series in a data.frame and in a list
  wide <- knn_wide(new_series[1:3])
  wide[2, paste0("BLUE.", knn_time)] <- NA
  expect_equal(is.na(predict(sample_model, newdata = wide, k = 5)$prediction), c(FALSE, TRUE, FALSE))
  ts_list <- new_series[1:3]
  ts_list[[2]]$BLUE <- NA
  expect_equal(is.na(predict(sample_model, newdata = ts_list, k = 5)$prediction), c(FALSE, TRUE, FALSE))
})

test_that("arguments of the TWDTW distance given to predict replace the arguments of the model", {
  other_model <- twdtw_knn(knn_samples, knn_sample_labels,
                           time_weight = c(steepness = 0.2, midpoint = 30),
                           cycle_length = 'year', time_scale = 'day', max_elapsed = 80)

  lu <- predict(sample_model, newdata = new_series, k = 3,
                time_weight = c(steepness = 0.2, midpoint = 30), max_elapsed = 80)
  expect_equal(lu, predict(other_model, newdata = new_series, k = 3))
  expect_false(isTRUE(all.equal(lu$distance, predict(sample_model, newdata = new_series, k = 3)$distance)))

  # The model is not changed
  expect_equal(sample_model$twdtw_args$time_weight, c(steepness = 0.1, midpoint = 50))

  # A time-weight function
  lu_fun <- predict(sample_model, newdata = new_series, k = 3,
                    time_weight = function(dist, el) dist + 1 / (1 + exp(-0.1 * (el - 50))))
  expect_equal(lu_fun, predict(sample_model, newdata = new_series, k = 3))
})

test_that("arguments that cannot be used to classify are an error", {
  x <- knn_array(new_series[1:4])
  n_patterns <- length(knn_samples)
  msg <- paste0("'k' must be an integer between one and the number of temporal patterns of the model (", n_patterns, ")")

  expect_error(predict(sample_model, newdata = new_series, k = 0), msg, fixed = TRUE)
  expect_error(predict(sample_model, newdata = new_series, k = n_patterns + 1), msg, fixed = TRUE)
  expect_error(predict(sample_model, newdata = new_series, k = 1.5), msg, fixed = TRUE)
  expect_error(predict(sample_model, newdata = new_series, k = NA), msg, fixed = TRUE)
  expect_error(predict(sample_model, newdata = new_series, k = "2"), msg, fixed = TRUE)
  expect_error(predict(sample_model, newdata = new_series, vote = "mean"), "'arg' should be one of")
  expect_error(predict(sample_model, newdata = new_series, vote = "none"), "'arg' should be one of")
  expect_error(predict(sample_model, newdata = new_series, type = "label"), "'arg' should be one of")
  expect_error(predict(sample_model), "The 'newdata' argument is missing", fixed = TRUE)

  # Arguments that are not known are not ignored
  expect_error(predict(sample_model, newdata = new_series, K = 3),
               "Unknown arguments in predict: 'K'.", fixed = TRUE)
  expect_error(predict(sample_model, newdata = new_series, time = knn_time),
               "Unknown arguments in predict: 'time'. 'time' is only used if 'newdata' is an array.", fixed = TRUE)
  expect_error(predict(sample_model, newdata = knn_wide(new_series), time = knn_time),
               "Unknown arguments in predict: 'time'. 'time' is only used if 'newdata' is an array.", fixed = TRUE)

  # Time series
  expect_error(predict(sample_model, newdata = x), "The 'time' argument is missing", fixed = TRUE)
  expect_error(predict(sample_model, newdata = x, time = knn_time[-1]),
               "The length of 'time' must be equal to the second dimension of 'newdata'", fixed = TRUE)
  expect_error(predict(sample_model, newdata = x[, , 1], time = knn_time),
               "'newdata' must be a numeric array with three dimensions: series, time, and band", fixed = TRUE)
  expect_error(predict(sample_model, newdata = 1:10), "No method of twdtw_knn_predict for 'newdata' of class 'integer'", fixed = TRUE)
  expect_error(predict(sample_model, newdata = list()), "'newdata' has no time series", fixed = TRUE)

  # A data.frame with a column for the time index is a single time series and must be in a list
  expect_error(predict(sample_model, newdata = new_series[[1]]),
               "'newdata' has no columns named with a band and a date", fixed = TRUE)

  # The time series must have the bands of the model
  expect_error(predict(sample_model, newdata = x[, , 1:5], time = knn_time),
               "The bands of newdata (EVI, NDVI, RED, BLUE, NIR) are not the bands of the model (EVI, NDVI, RED, BLUE, NIR, MIR)",
               fixed = TRUE)
  expect_error(predict(sample_model, newdata = unname(x[, , 1:5]), time = knn_time),
               "newdata has 5 bands and the model has 6", fixed = TRUE)
  expect_error(predict(sample_model, newdata = lapply(new_series[1:3], function(ts) ts[c("time", "EVI")])),
               "The bands of element 1 of 'newdata' (EVI) are not the bands of the model", fixed = TRUE)
  expect_error(predict(sample_model, newdata = knn_wide(new_series, bands = knn_bands[1:2])),
               "The bands of newdata (EVI, NDVI) are not the bands of the model", fixed = TRUE)

  # An object that is not a model
  expect_error(twdtw_knn_predict(x, list(data = 1), time = knn_time), "'object' is not a model created by twdtw_knn", fixed = TRUE)
})

test_that("a model with a single band classifies a data.frame with only the date in its column names", {
  evi <- function(ts_list) lapply(ts_list, function(ts) ts[c("time", "EVI")])
  m <- knn_model(evi(knn_samples), knn_sample_labels)
  lu <- predict(m, newdata = evi(new_series), k = 3)

  wide <- knn_wide(new_series, bands = "EVI")
  expect_equal(predict(m, newdata = wide, k = 3), lu)
  names(wide) <- sub("^EVI\\.", "", names(wide))
  expect_equal(predict(m, newdata = wide, k = 3), lu)
  names(wide) <- paste0("X", gsub("-", ".", names(wide), fixed = TRUE))
  expect_equal(predict(m, newdata = wide, k = 3), lu)

  # With more than one band the name of the band is needed
  expect_error(predict(sample_model, newdata = wide),
               "'newdata' has no columns named with a band and a date", fixed = TRUE)
})

test_that("other packages can add methods to classify their classes of time series", {
  # A class that keeps the values of all time series in an array and their dates
  pixels <- structure(list(values = knn_array(new_series), dates = knn_time), class = "knn_test_pixels")

  expect_error(predict(sample_model, newdata = pixels),
               "No method of twdtw_knn_predict for 'newdata' of class 'knn_test_pixels'", fixed = TRUE)

  # The method converts the time series to an array and arranges the result for its class
  registerS3method("twdtw_knn_predict", "knn_test_pixels", function(newdata, object, ...) {
    result <- twdtw_knn_predict(newdata$values, object, time = newdata$dates, ...)
    structure(list(result = result), class = "knn_test_map")
  }, envir = asNamespace("twdtw"))

  # predict passes k, vote, type, and the arguments of the TWDTW distance to the method
  for (type in c("class", "prob", "neighbors")) {
    map <- predict(sample_model, newdata = pixels, k = 4, vote = "weighted", type = type, max_elapsed = 90)
    expect_s3_class(map, "knn_test_map")
    expect_equal(map$result, predict(sample_model, newdata = new_series, k = 4, vote = "weighted", type = type, max_elapsed = 90))
  }
  expect_equal(predict(sample_model, newdata = pixels)$result, predict(sample_model, newdata = new_series))
})
