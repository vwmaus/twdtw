#' Convert date to numeric value based on time cycle
#'
#' @description This function converts a date or datetime column to a numeric format
#' based on a specified time cycle length and scale.
#'
#' @param date_column (Date or POSIXct) The column of dates to be converted.
#' @param time_cycle_length (character or integer) The length of the time cycle. If character, it should be one of: 'year', 'month', 'day', 'hour', 'minute', 'second'. If numeric, it represents the length of the cycle in arbitrary units.
#' @param time_cycle_scale (character or integer) The scale of the time cycle. If character, it should be one of: 'year', 'month', 'day', 'hour', 'minute', 'second'. If numeric, it represents the scale of the cycle in arbitrary units.
#'
#' @return A numeric vector with the same length as `date_column`, where each date has been converted to a numeric value based on the specified time cycle.
#'
#'
#' @examples
#'
#' convert_date_to_numeric(Sys.Date(), "year", "day") # Day of the year
#' convert_date_to_numeric(Sys.time(), "day", "hour") # Hour of the day
#'
#' @export
convert_date_to_numeric <- function(date_column, time_cycle_length, time_cycle_scale) {

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

  # Check if time_cycle_length and time_cycle_scale are one of the valid options
  valid_options <- c("year", "month", "day", "hour", "minute", "second")
  if (!(time_cycle_length %in% valid_options) | !(time_cycle_scale %in% valid_options)) {
    if (!is.numeric(time_cycle_length) | !is.numeric(time_cycle_scale)) {
      stop("time_cycle_length and time_cycle_scale must be one of: 'year', 'month', 'day', 'hour', 'minute', 'second' or numeric")
    }
  }

  # Convert the dates
  if (is.numeric(time_cycle_length) && is.numeric(time_cycle_scale)) {
    # Ensure time_cycle_length > time_cycle_scale
    if (time_cycle_length <= time_cycle_scale) {
      stop("For numeric time_cycle_length and time_cycle_scale, time_cycle_length must be greater than time_cycle_scale")
    }
    return(as.numeric(date_column) %% time_cycle_length / time_cycle_scale)
  } else {
    switch(time_cycle_length,
           "year" = switch(time_cycle_scale,
                           "month" = month(date_column),
                           "day" = yday(date_column),
                           "hour" = hour(date_column),
                           "minute" = minute(date_column),
                           "second" = second(date_column)),
           "month" = switch(time_cycle_scale,
                            "day" = mday(date_column),
                            "hour" = hour(date_column),
                            "minute" = minute(date_column),
                            "second" = second(date_column)),
           "day" = switch(time_cycle_scale,
                          "hour" = hour(date_column),
                          "minute" = minute(date_column),
                          "second" = second(date_column)),
           "hour" = switch(time_cycle_scale,
                           "minute" = minute(date_column),
                           "second" = second(date_column)),
           "minute" = switch(time_cycle_scale,
                             "second" = second(date_column))
    )
  }
}

#' Calculate the Maximum Possible Value of a Time Cycle
#'
#' This function returns the maximum possible value that a specific time
#' component can take, given a cycle length and scale.
#'
#' @param time_cycle_length A character string indicating the larger unit of time.
#' It must be one of "year", "month", "day", "hour", "minute".
#'
#' @param time_cycle_scale A character string indicating the smaller unit of time,
#' which is a division of the time_cycle_length. If time_cycle_length is "year",
#' time_cycle_scale can be one of "month", "day", "hour", "minute", "second".
#' If time_cycle_length is "month", time_cycle_scale can be "day", "hour", "minute",
#' "second", and so on.
#'
#' @return The function returns the maximum possible value that the time_cycle_scale
#' can take within one time_cycle_length.
#'
#' @examples
#'
#' calculate_max_cycle_length("year", "month")  # Maximum months is a year 12
#' calculate_max_cycle_length("day", "minute")  # Maximum minutes in a day 1440
#' calculate_max_cycle_length("year", "day")    # Maximum days in a year 366
#'
#' @export
calculate_max_cycle_length <- function(time_cycle_length, time_cycle_scale) {

  units_in_larger <- list(
    "year" = list("month" = 12, "day" = 366, "hour" = 24 * 366, "minute" = 60 * 24 * 366, "second" = 60 * 60 * 24 * 366),
    "month" = list("day" = 31, "hour" = 24 * 31, "minute" = 60 * 24 * 31, "second" = 60 * 60 * 24 * 31),
    "day" = list("hour" = 24, "minute" = 24 * 60, "second" = 24 * 60 * 60),
    "hour" = list("minute" = 60, "second" = 60 * 60),
    "minute" = list("second" = 60)
  )

  if (!(time_cycle_length %in% names(units_in_larger))) {
    stop("Invalid time_cycle_length")
  }

  if (!(time_cycle_scale %in% names(units_in_larger[[time_cycle_length]]))) {
    stop("Invalid time_cycle_scale for the provided time_cycle_length")
  }

  max_value <- units_in_larger[[time_cycle_length]][[time_cycle_scale]]
  return(max_value)
}


