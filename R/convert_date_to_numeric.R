#' Convert Date/POSIXct to a Numeric Cycle
#'
#' This function takes a date or datetime and converts it to a numeric cycle.
#' The cycle can be specified in units of years, months, days, hours, minutes, or seconds.
#' When cycle_length is a string, time_scale only changes the unit in which the result is expressed.
#' When cycle_length is numeric, time_scale and origin are used to compute the elapsed time.
#'
#' @param x A vector of dates or datetimes to convert. If not of type Date or POSIXct,
#' the function attempts to convert it.
#' @param cycle_length The length of the cycle. Can be a numeric value or a string
#' specifying the units ('year', 'month', 'day', 'hour', 'minute', 'second').
#' When numeric, the cycle length is in the same units as time_scale. When a string,
#' it specifies the time unit of the cycle.
#' @param time_scale Specifies the time scale for the conversion. Must be one of
#' 'year', 'month', 'day', 'hour', 'minute', 'second'. When cycle_length is a string,
#' time_scale changes the unit in which the result is expressed.
#' When cycle_length is numeric, time_scale is used to compute the elapsed time in seconds.
#' @param origin For numeric cycle_length, the origin must be specified. This is the point
#' from which the elapsed time is computed. Must be of the same class as x.
#' @return The numeric cycle value(s) corresponding to x.
#' @examples
#' date_to_numeric_cycle(Sys.time(), "year", "day") # Returns the day of the year
#' date_to_numeric_cycle(Sys.time(), "day", "hour") # Returns the hour of the day
#' @export
date_to_numeric_cycle <- function(x, cycle_length, time_scale, origin = NULL) {

  # Check if the x is of type Date or POSIXct, if not try to convert it
  x <- to_date_time(x)

  # Define the time scale in seconds
  time_scale_seconds <- switch(time_scale,
                               "year" = 60 * 60 * 24 * 365.25, # Note: Using 365.25 to account for leap years
                               "month" = 60 * 60 * 24 * 30, # Note: This approximation does not account for months with 28, 29, or 31 days
                               "day" = 60 * 60 * 24,
                               "hour" = 60 * 60,
                               "minute" = 60,
                               "second" = 1,
                               stop("Invalid time_scale, must be one of 'year', 'month', 'day', 'hour', 'minute', 'second'")
  )

  if (is.numeric(cycle_length)) {
    if (is.null(origin)) {
      stop("For numeric cycle_length, origin must be given.")
    }
    origin <- to_date_time(origin)
    if (!identical(class(x), class(origin))) {
      stop("x and origin must be of the same class Date or POSIXct")
    }

    # Compute elapsed time from origin in the appropriate scale
    elapsed_time <- as.numeric(difftime(x, origin, units = "secs"))

    # Compute cycle
    cycle_value <- (elapsed_time / time_scale_seconds) %% cycle_length + 1

    return(cycle_value)

  } else {

    # Check if cycle_length is one of the valid options
    valid_options <- c("year", "month", "day", "hour", "minute", "second")
    if (!(cycle_length %in% valid_options)) {
      stop("cycle_length must be one of: 'year', 'month', 'day', 'hour', 'minute', or 'second'")
    }

    # Compute cycle
    switch(cycle_length,
           "year" = switch(time_scale,
                           "month" = as.numeric(format(x, "%m")),
                           "day" = as.numeric(format(x, "%j")),
                           "hour" = as.numeric(format(x, "%j"))*24 + as.numeric(format(x, "%H")),
                           "minute" = (as.numeric(format(x, "%j"))*24 + as.numeric(format(x, "%H")))*60 + as.numeric(format(x, "%M")),
                           "second" = ((as.numeric(format(x, "%j"))*24 + as.numeric(format(x, "%H")))*60 + as.numeric(format(x, "%M")))*60 + as.numeric(format(x, "%S"))),
           "month" = switch(time_scale,
                            "day" = as.numeric(format(x, "%d")),
                            "hour" = as.numeric(format(x, "%d"))*24 + as.numeric(format(x, "%H")),
                            "minute" = (as.numeric(format(x, "%d"))*24 + as.numeric(format(x, "%H")))*60 + as.numeric(format(x, "%M")),
                            "second" = ((as.numeric(format(x, "%d"))*24 + as.numeric(format(x, "%H")))*60 + as.numeric(format(x, "%M")))*60 + as.numeric(format(x, "%S"))),
           "day" = switch(time_scale,
                          "hour" = as.numeric(format(x, "%H")),
                          "minute" = as.numeric(format(x, "%H"))*60 + as.numeric(format(x, "%M")),
                          "second" = (as.numeric(format(x, "%H"))*60 + as.numeric(format(x, "%M")))*60 + as.numeric(format(x, "%S"))),
           "hour" = switch(time_scale,
                           "minute" = as.numeric(format(x, "%M")),
                           "second" = as.numeric(format(x, "%M"))*60 + as.numeric(format(x, "%S"))),
           "minute" = switch(time_scale,
                             "second" = as.numeric(format(x, "%S")))
    )
  }
}

to_date_time <- function(x){
  if (!inherits(x, c("Date", "POSIXt"))) {
    # check if all strings in the vector include hours, minutes, and seconds
    if (all(grepl("\\d{4}-\\d{2}-\\d{2} \\d{2}:\\d{2}:\\d{2}", x))) {
      x <- try(as.POSIXct(x), silent = TRUE)
    } else {
      x <- try(as.Date(x), silent = TRUE)
    }
    if (inherits(x, "try-error")) {
      stop("Some elements of x could not be converted to a date or datetime format")
    }
  }
  return(x)
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
#' max_cycle_length("year", "month")  # Maximum months is a year 12
#' max_cycle_length("day", "minute")  # Maximum minutes in a day 1440
#' max_cycle_length("year", "day")    # Maximum days in a year 366
#'
#' @export
max_cycle_length <- function(cycle_length, time_scale) {

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


