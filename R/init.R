#' @useDynLib twdtw, .registration = TRUE
#' @import Rcpp
NULL

#' @importFrom dtw symmetric1
#' @export
dtw::symmetric1

#' @importFrom dtw symmetric2
#' @export
dtw::symmetric2

#' @importFrom dtw asymmetric
#' @export
dtw::asymmetric

#' @importFrom dtw rabinerJuangStepPattern
#' @export
dtw::rabinerJuangStepPattern
