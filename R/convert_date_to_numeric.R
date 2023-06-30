#' Convert date to numeric value based on time cycle
#'
#' @description This function converts a date or datetime column to a numeric format
#' based on a specified time cycle length and scale.
#'
#' @param date_column (Date or POSIXct) The column of dates to be converted.
#' @param cycle_length (character or integer) The length of the time cycle.
#' If character, it should be one of: 'year', 'month', 'day', 'hour', 'minute', 'second'.
#' If numeric, it represents the length of the cycle in arbitrary units.
#' @param time_scale (character or integer) The scale of the time cycle.
#' If character, it should be one of: 'year', 'month', 'day', 'hour', 'minute', 'second'.
#' If numeric, it represents the scale of the cycle in arbitrary units.
#'
#' @return A numeric vector with the same length as `date_column`,
#' where each date has been converted to a numeric value based on the specified time cycle and scale.
#'
#'
#' @examples
#'
#' convert_date_to_numeric(Sys.Date(), "year", "day") # Day of the year
#' convert_date_to_numeric(Sys.time(), "day", "hour") # Hour of the day
#'
#' @export
convert_date_to_numeric <- function(date_column, cycle_length, time_scale) {

  # Check if the date_column is of type Date or POSIXct, if not try to convert it
  if (!inherits(date_column, c("Date", "POSIXt"))) {
    date_column <- try(as.Date(date_column), silent = TRUE)
    if (inherits(date_column, "try-error")) {
      date_column <- try(as.POSIXct(date_column), silent = TRUE)
    }
    if (inherits(date_column, "try-error")) {
      stop("index_column_name could not be converted to a date or datetime format")
    }
  }

  # Check if cycle_length and time_scale are one of the valid options
  valid_options <- c("year", "month", "day", "hour", "minute", "second")
  if (!(cycle_length %in% valid_options) | !(time_scale %in% valid_options)) {
    if (!is.numeric(cycle_length) | !is.numeric(time_scale)) {
      stop("cycle_length and time_scale must be one of: 'year', 'month', 'day', 'hour', 'minute', 'second' or numeric")
    }
  }

  # Convert the dates
  if (is.numeric(cycle_length) && is.numeric(time_scale)) {
    # Ensure cycle_length > time_scale
    if (cycle_length <= time_scale) {
      stop("For numeric cycle_length and time_scale, cycle_length must be greater than time_scale")
    }
    return(as.numeric(date_column) %% cycle_length / time_scale)
  } else {
    switch(cycle_length,
           "year" = switch(time_scale,
                           "month" = month(date_column),
                           "day" = yday(date_column),
                           "hour" = hour(date_column),
                           "minute" = minute(date_column),
                           "second" = second(date_column)),
           "month" = switch(time_scale,
                            "day" = mday(date_column),
                            "hour" = hour(date_column),
                            "minute" = minute(date_column),
                            "second" = second(date_column)),
           "day" = switch(time_scale,
                          "hour" = hour(date_column),
                          "minute" = minute(date_column),
                          "second" = second(date_column)),
           "hour" = switch(time_scale,
                           "minute" = minute(date_column),
                           "second" = second(date_column)),
           "minute" = switch(time_scale,
                             "second" = second(date_column))
    )
  }
}

#' Calculate the Maximum Possible Value of a Time Cycle
#'
#' This function returns the maximum possible value that a specific time
#' component can take, given a cycle length and scale.
#'
#' @param cycle_length A character string indicating the larger unit of time.
#' It must be one of "year", "month", "day", "hour", "minute".
#'
#' @param time_scale A character string indicating the smaller unit of time,
#' which is a division of the \code{cycle_length}. If \code{cycle_length} is "year",
#' \code{time_scale} can be one of "month", "day", "hour", "minute", "second".
#' If \code{cycle_length} is "month", \code{time_scale} can be "day", "hour", "minute",
#' "second", and so on.
#'
#' @return The function returns the maximum possible value that the \code{time_scale}
#' can take within one \code{cycle_length}.
#'
#' @examples
#'
#' calculate_max_cycle_length("year", "month")  # Maximum months is a year 12
#' calculate_max_cycle_length("day", "minute")  # Maximum minutes in a day 1440
#' calculate_max_cycle_length("year", "day")    # Maximum days in a year 366
#'
#' @export
calculate_max_cycle_length <- function(cycle_length, time_scale) {

  units_in_larger <- list(
    "year" = list("month" = 12, "day" = 366, "hour" = 24 * 366, "minute" = 60 * 24 * 366, "second" = 60 * 60 * 24 * 366),
    "month" = list("day" = 31, "hour" = 24 * 31, "minute" = 60 * 24 * 31, "second" = 60 * 60 * 24 * 31),
    "day" = list("hour" = 24, "minute" = 24 * 60, "second" = 24 * 60 * 60),
    "hour" = list("minute" = 60, "second" = 60 * 60),
    "minute" = list("second" = 60)
  )

  if (!(cycle_length %in% names(units_in_larger))) {
    stop("Invalid cycle_length")
  }

  if (!(time_scale %in% names(units_in_larger[[cycle_length]]))) {
    stop("Invalid time_scale for the provided cycle_length")
  }

  max_value <- units_in_larger[[cycle_length]][[time_scale]]
  return(max_value)
}


