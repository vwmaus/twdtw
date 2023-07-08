#' Plot TWDTW cost matrix
#'
#' This function visualizes the Time-Weighted Dynamic Time Warping cost matrix.
#'
#' @param twdtw_obj An object of class 'twdtw' including internal data.
#' @param ... Additional arguments passed to \code{\link[graphics]{image}}.
#'
#' @examples
#' \dontrun{
#' # Create a time series
#' n <- 23
#' t <- seq(0, pi, length.out = n)
#' d <- seq(as.Date('2020-09-01'), length.out = n, by = "15 day")
#'
#' x <- data.frame(time = d,      v1 = sin(t)*2 + runif(n))
#'
#' # shift time by 30 days
#' y <- data.frame(time = d + 30, v1 = sin(t)*2 + runif(n))
#'
#' plot(x, type = "l", xlim = range(c(d, d + 5)))
#' lines(y, col = "red")
#'
#' Call twdtw using "output = 'internals'"
#' twdtw_obj <- twdtw(x, y,
#'        cycle_length = 'year',
#'        time_scale = 'day',
#'        time_weight = c(steepness = 0.1, midpoint = 50), output = 'internals')
#'
#' plot_cost_matrix(twdtw_obj)
#'
#' }
#' @export
plot_cost_matrix <- function(twdtw_obj, ...) {

  if (!inherits(twdtw_obj, "twdtw")) {
    stop("twdtw_obj must be of class 'twdtw'")
  }

  if (!"CM" %in% names(twdtw_obj)) {
    stop("Cost matrix (CM) not found in twdtw_obj")
  }

  # Get the cost matrix
  cost_matrix <- t(twdtw_obj$CM[-1,])

  # Create the plot
  image(cost_matrix, ..., xlab = "Time Series X", ylab = "Time Series Y", main = "TWDTW Cost Matrix")

}
