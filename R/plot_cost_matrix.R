#' Plot TWDTW cost matrix
#'
#' This function visualizes the Time-Weighted Dynamic Time Warping cost matrix.
#'
#' @param x An object of class 'twdtw' including internal data.
#' @param ... Additional arguments passed to \code{\link[graphics]{image}}.
#'
#' @return An image plot of the TWDTW cost matrix. The x-axis represents the time series x,
#'         and the y-axis represents the time series y. The cost matrix is color-coded,
#'         with darker shades indicating higher costs and lighter shades indicating lower costs.
#'         No object is returned by this function; the plot is directly outputted to the active device.
#'
#' @examples
#'
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
#' # Call twdtw using "output = 'internals'
#' twdtw_obj <- twdtw(x, y,
#'        cycle_length = 'year',
#'        time_scale = 'day',
#'        time_weight = c(steepness = 0.1, midpoint = 50), output = 'internals')
#'
#' plot_cost_matrix(twdtw_obj)
#'
#' @export
plot_cost_matrix <- function(x, ...) {

  if (!inherits(x, "twdtw")) {
    stop("x must be of class 'twdtw'")
  }

  if (!"CM" %in% names(x)) {
    stop("Cost matrix (CM) not found in x")
  }

  # Get the cost matrix
  cost_matrix <- t(x$CM[-1,])

  # Create the plot
  image(cost_matrix, ..., xlab = "Time Series x", ylab = "Time Series y", main = "TWDTW Cost Matrix")
}
