#' @title Classify time series using a TWDTW k-nearest neighbour model
#'
#' @description
#' This function classifies time series with a model created by \code{\link{twdtw_knn}}.
#' Each time series gets a label from the labels of its \code{k} nearest temporal patterns,
#' i.e. the patterns of the model with the lowest Time-Weighted Dynamic Time Warping (TWDTW) distance.
#' By default \code{k = 1} and each time series gets the label of the nearest pattern.
#'
#' @param object A model of class \code{twdtw_knn}, created by \code{\link{twdtw_knn}}.
#' @param newdata The time series to classify, in one of the three forms of the time series of
#' \code{\link{twdtw_knn}}: a list of data.frames or matrices, a data.frame with one row for each time series
#' and one column for each band and date, or a numeric array with the dimensions series, time, and band.
#' @param k The number of nearest temporal patterns used for each time series. An integer
#' between one and the number of temporal patterns of the model. Defaults to 1, i.e. only the
#' nearest pattern is used.
#' @param vote The weight of the vote of each of the \code{k} nearest patterns. Either
#' \code{"majority"} (default), where each pattern has one vote, or \code{"weighted"}, where each
#' pattern votes with the weight 1/distance. It can change the result only if \code{k} is larger
#' than one. See details.
#' @param type What the function returns for each time series. One of:
#' \itemize{
#'   \item \code{"class"} (default): the predicted label, i.e. the label with the most votes,
#'   and the distance to the nearest pattern with this label.
#'   \item \code{"prob"}: a score from 0 to 1 for each label of the model, i.e. the share of the
#'   votes that the label gets.
#'   \item \code{"neighbors"}: the labels and the distances of the \code{k} nearest patterns,
#'   without a vote.
#' }
#' See details.
#' @param time A vector with the dates of the time dimension of an array \code{newdata}.
#' Its length must be equal to \code{dim(newdata)[2]}. It is only used if \code{newdata} is an array.
#' @param ... Arguments of the TWDTW distance that replace the arguments of the model:
#' \code{time_weight}, \code{cycle_length}, \code{time_scale}, \code{origin}, and \code{max_elapsed}.
#' For a list \code{newdata} also \code{index_column}, the column name or the column position of the
#' time index, which defaults to "time". An argument that is not known is an error.
#'
#' @details A time series is classified in three steps, which are set by \code{k}, \code{vote},
#' and \code{type}:
#'
#' \enumerate{
#'   \item The \code{k} temporal patterns of the model with the lowest TWDTW distance to the time series
#'   are found.
#'   \item Each of these patterns votes for its label. With \code{vote = "majority"} each pattern has
#'   one vote. With \code{vote = "weighted"} each pattern votes with the weight 1/distance, so that a
#'   pattern twice as far has half the weight.
#'   \item \code{type} sets what is returned: the label with the most votes (\code{"class"}), the share
#'   of the votes of each label (\code{"prob"}), or the \code{k} nearest patterns without a vote
#'   (\code{"neighbors"}).
#' }
#'
#' With the defaults \code{k = 1}, \code{vote = "majority"}, and \code{type = "class"}, each time series
#' gets the label of its nearest pattern.
#'
#' As an example, let the three nearest patterns of a time series (\code{k = 3}) have the labels
#' soybean, maize, and maize, and the distances 1, 4, and 4:
#'
#' \itemize{
#'   \item With \code{vote = "majority"}, soybean has one vote and maize has two. \code{type = "class"}
#'   returns maize and the distance 4, and \code{type = "prob"} returns 1/3 for soybean and 2/3 for maize.
#'   \item With \code{vote = "weighted"}, soybean has the vote 1 and maize has the votes 1/4 + 1/4.
#'   \code{type = "class"} returns soybean and the distance 1, and \code{type = "prob"} returns
#'   2/3 for soybean and 1/3 for maize.
#'   \item \code{type = "neighbors"} returns the three labels and the three distances, with any \code{vote}.
#' }
#'
#' The three types of result in more detail:
#'
#' \itemize{
#'   \item \code{"class"}: if several labels have the same number of votes, the prediction is the one of
#'   them whose pattern is nearest. The distance is the TWDTW distance to the nearest pattern with the
#'   predicted label.
#'   \item \code{"prob"}: the share of a label is its votes divided by the votes of all \code{k} patterns.
#'   The shares of a time series sum to one, a label without a pattern among the \code{k} nearest has
#'   the share zero, and the label of \code{type = "class"} always has the largest share. With
#'   \code{k = 1} the share is one for the label of the nearest pattern and zero for the other labels.
#'   The name follows other classifiers in R, where \code{type = "prob"} returns a score for each class.
#'   \item \code{"neighbors"}: the patterns are ordered from the nearest to the farthest. This is the
#'   result to use for another way of combining the labels and the distances.
#' }
#'
#' With \code{vote = "weighted"}, patterns at distance zero get all the votes. Negative distances,
#' which are possible only with a time-weight function, are an error with this vote.
#'
#' A vote can only change the prediction if the model has more than one temporal pattern for a label,
#' as a model created without \code{resampling_fun}. With one pattern for each label, every label has
#' one vote: the prediction is the label of the nearest pattern for any \code{k}, and
#' \code{vote = "weighted"} with \code{k} equal to the number of patterns gives shares that are
#' proportional to 1/distance.
#'
#' The votes use the nearest patterns that are found for any \code{type}, so the shares of the
#' votes need no additional TWDTW distances. A larger \code{k} takes more time, because the calculation
#' of the distance to a pattern stops as soon as the pattern cannot be among the \code{k} nearest.
#'
#' The TWDTW distances between all time series and all temporal patterns are calculated in
#' a single call to \code{\link{twdtw_nearest}}. An array \code{newdata} is passed to this function
#' as it is and is the fastest form. The other forms are converted to arrays.
#'
#' The time series must have the bands of the temporal patterns of the model. An observation with a missing
#' value in any band is not used, and the time series is classified with its remaining observations.
#' A time series without any complete observation has \code{NA} in all columns of the result.
#'
#' The time series are converted by \code{\link{twdtw_knn_predict}}. Other packages can add methods to
#' this function to classify their classes of time series.
#'
#' @return A data.frame with one row for each time series in \code{newdata}, in the same order.
#' The columns depend on \code{type}:
#' \item{\code{"class"}}{\code{prediction}, a factor with the predicted label, and \code{distance},
#' the TWDTW distance to the nearest temporal pattern with this label.}
#' \item{\code{"prob"}}{One column for each label of the model, named with the label, with its share of the votes.}
#' \item{\code{"neighbors"}}{\code{prediction_1}, ..., \code{prediction_k} and \code{distance_1}, ..., \code{distance_k},
#' where 1 is the nearest temporal pattern.}
#' The predictions are factors with the labels of the model as levels, in alphabetical order.
#'
#' @seealso \code{\link{twdtw_knn}}, \code{\link{twdtw_nearest}}
#'
#' @examples
#'
#' # Temporal patterns of 5 land use classes and satellite image time series of 999 pixels
#' ts_y <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))
#' ts_x <- readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))
#'
#' m <- twdtw_knn(x = ts_y$observations,
#'                y = ts_y$label,
#'                cycle_length = 'year',
#'                time_scale = 'day',
#'                time_weight = c(steepness = 0.1, midpoint = 50))
#'
#' # With the defaults each time series gets the label of its nearest pattern
#' lu <- predict(m, newdata = ts_x$observations)
#'
#' head(lu)
#'
#' # The model has one pattern for each label: weighted votes of all patterns
#' # give shares that are proportional to 1/distance
#' head(predict(m, newdata = ts_x$observations, k = 5, vote = "weighted", type = "prob"))
#'
#' # A model that keeps 300 labelled time series as temporal patterns:
#' # the labels have several patterns, and the votes of the nearest patterns can be divided
#' m_samples <- twdtw_knn(x = ts_x$observations[1:300],
#'                        y = lu$prediction[1:300],
#'                        cycle_length = 'year',
#'                        time_scale = 'day',
#'                        time_weight = c(steepness = 0.1, midpoint = 50))
#'
#' new_ts <- ts_x$observations[c(301, 307, 791)]
#'
#' # type = "class": the label with the most votes among the 5 nearest patterns
#' predict(m_samples, newdata = new_ts, k = 5)
#'
#' # type = "prob": the share of the votes of each label
#' predict(m_samples, newdata = new_ts, k = 5, type = "prob")
#'
#' # Votes with the weight 1/distance change the shares,
#' # and the label of the third time series
#' predict(m_samples, newdata = new_ts, k = 5, vote = "weighted", type = "prob")
#'
#' predict(m_samples, newdata = new_ts, k = 5, vote = "weighted")
#'
#' # type = "neighbors": the labels and the distances of the nearest patterns
#' predict(m_samples, newdata = new_ts, k = 3, type = "neighbors")
#'
#' # The fastest form is an array with the dimensions series, time, and band
#' bands <- c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")
#' x <- sapply(ts_x$observations, function(ts) as.matrix(ts[bands]), simplify = "array")
#' x <- aperm(x, c(3, 1, 2))
#'
#' lu_array <- predict(m, newdata = x, time = ts_x$observations[[1]]$time)
#'
#' all.equal(lu, lu_array)
#'
#' @export
predict.twdtw_knn <- function(object, newdata,
                              k = 1,
                              vote = c("majority", "weighted"),
                              type = c("class", "prob", "neighbors"),
                              time = NULL, ...) {

  if (missing(newdata)) {
    stop("The 'newdata' argument is missing: the time series to classify")
  }

  vote <- match.arg(vote)
  type <- match.arg(type)

  # 'time' is passed only if it is given, so that methods of other packages
  # for time series that have their own dates do not receive it
  if (is.null(time)) {
    return(twdtw_knn_predict(newdata, object, k = k, vote = vote, type = type, ...))
  }

  twdtw_knn_predict(newdata, object, k = k, vote = vote, type = type, time = time, ...)

}

#' @title Classify time series of a given class with a TWDTW k-nearest neighbour model
#'
#' @description
#' This generic function is called by \code{\link{predict.twdtw_knn}} to classify the time series in
#' \code{newdata}, and the method is chosen by the class of \code{newdata}. It is exported for packages
#' that add methods for their classes of time series, e.g. the pixels of a satellite image data cube.
#' To classify time series call \code{\link[=predict.twdtw_knn]{predict}}.
#'
#' @param newdata The time series to classify. The methods of this package are for a list of data.frames
#' or matrices, a data.frame with one row for each time series and one column for each band and date,
#' and a numeric array with the dimensions series, time, and band, see \code{\link{twdtw_knn}}.
#' @param object A model of class \code{twdtw_knn}, created by \code{\link{twdtw_knn}}.
#' @param index_column The column name or the column position of the time index of the time series in a list.
#' Defaults to "time".
#' @param ... Arguments passed by \code{\link{predict.twdtw_knn}} to the methods: \code{k}, \code{vote},
#' \code{type}, and the arguments of the TWDTW distance that replace the arguments of the model.
#' \code{time} is passed only if it is given to \code{predict}.
#' @inheritParams predict.twdtw_knn
#'
#' @details \code{predict} chooses its method by the class of the model, which is its first argument.
#' This function chooses its method by the class of the time series, so that a package can add a method
#' for its class without changing the model.
#'
#' A method for another class usually converts \code{newdata} to a numeric array with the dimensions
#' series, time, and band, calls \code{twdtw_knn_predict} with this array, the dates of its time dimension
#' in \code{time}, and the other arguments in \code{...}, and then arranges the rows of the result in the
#' form of its class, see the examples. The method for an array passes the array to
#' \code{\link{twdtw_nearest}} as it is, so that missing values remain missing.
#'
#' @return The methods of this package return a data.frame with one row for each time series in
#' \code{newdata}, in the same order, see \code{\link{predict.twdtw_knn}}.
#'
#' @seealso \code{\link{predict.twdtw_knn}}, \code{\link{twdtw_knn}}, \code{\link{twdtw_nearest}}
#'
#' @examples
#'
#' # Temporal patterns of 5 land use classes and satellite image time series of 999 pixels
#' ts_y <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))
#' ts_x <- readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))
#'
#' m <- twdtw_knn(x = ts_y$observations,
#'                y = ts_y$label,
#'                cycle_length = 'year',
#'                time_scale = 'day',
#'                time_weight = c(steepness = 0.1, midpoint = 50))
#'
#' # A class of time series that keeps the values of all time series in an array
#' # and the dates of the observations
#' bands <- c("EVI", "NDVI", "RED", "BLUE", "NIR", "MIR")
#' values <- sapply(ts_x$observations, function(ts) as.matrix(ts[bands]), simplify = "array")
#' pixels <- structure(list(values = aperm(values, c(3, 1, 2)),
#'                          dates = ts_x$observations[[1]]$time),
#'                     class = "pixel_series")
#'
#' # A method for this class: pass the array and its dates to the method for arrays
#' twdtw_knn_predict.pixel_series <- function(newdata, object, ...) {
#'   twdtw_knn_predict(newdata$values, object, time = newdata$dates, ...)
#' }
#'
#' # predict() now classifies objects of this class
#' head(predict(m, newdata = pixels, k = 2, type = "neighbors"))
#'
#' @export
twdtw_knn_predict <- function(newdata, object, ...) {
  UseMethod("twdtw_knn_predict")
}

#' @export
#' @rdname twdtw_knn_predict
twdtw_knn_predict.array <- function(newdata, object,
                                    k = 1,
                                    vote = c("majority", "weighted"),
                                    type = c("class", "prob", "neighbors"),
                                    time = NULL, ...) {

  vote <- match.arg(vote)
  type <- match.arg(type)

  newdata <- .check_time_series_array(newdata, time, "newdata")

  .twdtw_knn_classify(object = object,
                      x = newdata$values,
                      time = newdata$time,
                      k = k,
                      vote = vote,
                      type = type,
                      new_twdtw_args = list(...))

}

#' @export
#' @rdname twdtw_knn_predict
twdtw_knn_predict.data.frame <- function(newdata, object,
                                         k = 1,
                                         vote = c("majority", "weighted"),
                                         type = c("class", "prob", "neighbors"), ...) {

  vote <- match.arg(vote)
  type <- match.arg(type)

  # Columns with only the date in their name have the values of the band of
  # a model with a single band
  bands <- .model_bands(object)
  newdata <- .wide_table_array(newdata,
                               single_band = if (length(bands) == 1) bands else NULL,
                               arg = "newdata")

  .twdtw_knn_classify(object = object,
                      x = newdata$values,
                      time = newdata$time,
                      k = k,
                      vote = vote,
                      type = type,
                      new_twdtw_args = list(...))

}

#' @export
#' @rdname twdtw_knn_predict
twdtw_knn_predict.list <- function(newdata, object,
                                   k = 1,
                                   vote = c("majority", "weighted"),
                                   type = c("class", "prob", "neighbors"),
                                   index_column = 'time', ...) {

  vote <- match.arg(vote)
  type <- match.arg(type)

  # Time series with the same dates are classified in the same call
  groups <- .time_series_arrays(newdata, index_column, .model_bands(object), "newdata")

  results <- lapply(groups, function(group) {
    .twdtw_knn_classify(object = object,
                        x = group$values,
                        time = group$time,
                        k = k,
                        vote = vote,
                        type = type,
                        new_twdtw_args = list(...))
  })

  if (length(groups) == 1) {
    return(results[[1]])
  }

  # Arrange the results in the order of the time series in newdata
  result <- results[[1]][rep(NA_integer_, length(newdata)), , drop = FALSE]
  rownames(result) <- NULL
  for (i in seq_along(groups)) {
    result[groups[[i]]$rows, ] <- results[[i]]
  }

  result

}

#' @export
twdtw_knn_predict.default <- function(newdata, object, ...) {

  stop("No method of twdtw_knn_predict for 'newdata' of class '", class(newdata)[1], "'. ",
       "'newdata' must be a list of data.frames or matrices, a data.frame with one row for each time series, ",
       "or a numeric array with the dimensions series, time, and band")

}

# Internal function
# Band names of the temporal patterns of a model
.model_bands <- function(object) {

  .check_model(object)

  setdiff(names(object$data$observations[[1]]), "time")

}

# Internal function
# Check that an object has the components of a model created by twdtw_knn
.check_model <- function(object) {

  if (!is.list(object) || !is.list(object$twdtw_args) || is.null(object$data$label) ||
      length(object$data$observations) < 1) {
    stop("'object' is not a model created by twdtw_knn")
  }

  invisible(NULL)

}

# Internal function
# Classify the time series in an array with the dimensions series, time, and
# band: find the k nearest temporal patterns of the model and combine their labels
.twdtw_knn_classify <- function(object, x, time, k, vote, type, new_twdtw_args) {

  .check_model(object)

  # Arguments of the TWDTW distance that replace the arguments of the model
  twdtw_args <- object$twdtw_args
  arg_names <- names(new_twdtw_args)
  if (is.null(arg_names)) {
    arg_names <- rep("", length(new_twdtw_args))
  }

  known_args <- c('time_weight', 'cycle_length', 'time_scale', 'origin', 'max_elapsed')
  unknown <- !arg_names %in% known_args
  if (any(unknown)) {
    .check_unknown_args(new_twdtw_args[unknown], "predict", "newdata")
  }

  for (name in arg_names) {
    twdtw_args[name] <- list(new_twdtw_args[[name]])
  }

  if (is.null(twdtw_args$max_elapsed)) {
    twdtw_args$max_elapsed <- Inf
  }

  # Temporal patterns of the model and their bands
  patterns <- lapply(object$data$observations, as.data.frame)
  labels <- as.character(object$data$label)
  pattern_bands <- setdiff(names(patterns[[1]]), 'time')

  # The neighbours are temporal patterns of the model
  if (!is.numeric(k) || length(k) != 1 || is.na(k) || k != round(k) || k < 1 || k > length(patterns)) {
    stop(paste0("'k' must be an integer between one and the number of temporal patterns of the model (",
                length(patterns), ")"))
  }

  # The time series must have the bands of the temporal patterns
  # Without band names the bands are matched by their position
  newdata_bands <- dimnames(x)[[3]]
  if (is.null(newdata_bands)) {
    if (dim(x)[3] != length(pattern_bands)) {
      stop(paste0("newdata has ", dim(x)[3], " bands and the model has ", length(pattern_bands),
                  " (", paste(pattern_bands, collapse = ", "), ")"))
    }
  } else if (!setequal(newdata_bands, pattern_bands)) {
    stop(paste0("The bands of newdata (", paste(newdata_bands, collapse = ", "),
                ") are not the bands of the model (", paste(pattern_bands, collapse = ", "), ")"))
  }

  if (dim(x)[2] > 0) {

    # Find the k nearest temporal patterns of each time series
    nearest <- twdtw_nearest(x = x,
                             y = patterns,
                             time = time,
                             time_weight = twdtw_args$time_weight,
                             cycle_length = twdtw_args$cycle_length,
                             time_scale = twdtw_args$time_scale,
                             k = k,
                             origin = twdtw_args$origin,
                             index_column = 'time',
                             max_elapsed = twdtw_args$max_elapsed)

  } else {

    # Time series without observations have no neighbours
    nearest <- list(index = matrix(NA_integer_, nrow = dim(x)[1], ncol = k),
                    distance = matrix(NA_real_, nrow = dim(x)[1], ncol = k))

  }

  result <- .knn_result(index = nearest$index,
                        distance = nearest$distance,
                        labels = labels,
                        vote = vote,
                        type = type)

  as.data.frame(result, check.names = FALSE)

}

# Internal function
# Combine the labels of the k nearest temporal patterns of each time series.
# 'index' and 'distance' are matrices with one row for each time series and k
# columns, ordered from the nearest to the farthest pattern, and 'labels' has
# the labels of the temporal patterns of the model. Returns a list of vectors
# with one element for each time series. A time series without neighbours,
# i.e. without any complete observation, has NA
.knn_result <- function(index, distance, labels, vote, type) {

  n <- nrow(index)
  k <- ncol(index)

  # The labels of the model in alphabetical order are the levels of the predictions
  classes <- sort(unique(labels))
  as_label <- function(class_id) {
    structure(as.integer(class_id), levels = classes, class = "factor")
  }

  # Label of the k nearest patterns as its position in 'classes'
  neighbor_class <- matrix(match(labels, classes)[index], nrow = n, ncol = k)

  if (type == "neighbors") {
    result <- c(lapply(seq_len(k), function(j) as_label(neighbor_class[, j])),
                lapply(seq_len(k), function(j) distance[, j]))
    names(result) <- c(paste0("prediction_", seq_len(k)), paste0("distance_", seq_len(k)))
    return(result)
  }

  # A single neighbour gives its label
  if (type == "class" && k == 1) {
    return(list(prediction = as_label(neighbor_class[, 1]), distance = distance[, 1]))
  }

  # Weight of the vote of each of the k nearest patterns
  weight <- matrix(1, nrow = n, ncol = k)
  if (vote == "weighted") {
    # A time-weight function can give negative distances, which cannot be used as weights
    if (any(distance < 0, na.rm = TRUE)) {
      stop("vote = 'weighted' cannot be used with negative TWDTW distances, use vote = 'majority'")
    }
    weight <- 1 / distance
    # Patterns at distance zero get all the votes
    zero_distance <- which(rowSums(distance == 0, na.rm = TRUE) > 0)
    weight[zero_distance, ] <- as.numeric(distance[zero_distance, ] == 0)
    # If no pattern has a weight, e.g. if all distances are infinite, each pattern has one vote
    no_weight <- which(rowSums(weight, na.rm = TRUE) == 0)
    weight[no_weight, ] <- 1
  }

  # Votes for each label
  votes <- matrix(0, nrow = n, ncol = length(classes))
  for (j in seq_len(k)) {
    has_neighbor <- which(!is.na(neighbor_class[, j]))
    id <- cbind(has_neighbor, neighbor_class[has_neighbor, j])
    votes[id] <- votes[id] + weight[has_neighbor, j]
  }

  if (type == "prob") {
    all_votes <- rowSums(votes)
    all_votes[is.na(neighbor_class[, 1])] <- NA
    result <- lapply(seq_along(classes), function(i) votes[, i] / all_votes)
    names(result) <- classes
    return(result)
  }

  most_votes <- votes[, 1]
  for (i in seq_along(classes)[-1]) {
    most_votes <- pmax(most_votes, votes[, i])
  }

  # The neighbours are ordered by distance, so the first of them that has a label with
  # the most votes gives the prediction, and its distance is the lowest for this label
  prediction <- rep(NA_integer_, n)
  prediction_distance <- rep(NA_real_, n)
  for (j in seq_len(k)) {
    not_predicted <- which(is.na(prediction) & !is.na(neighbor_class[, j]))
    not_predicted <- not_predicted[votes[cbind(not_predicted, neighbor_class[not_predicted, j])] == most_votes[not_predicted]]
    prediction[not_predicted] <- neighbor_class[not_predicted, j]
    prediction_distance[not_predicted] <- distance[not_predicted, j]
  }

  list(prediction = as_label(prediction), distance = prediction_distance)

}
