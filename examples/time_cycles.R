z <- "2023-01-24 23:59:00 CEST"
lubridate::second(z) # 0-60  # second
lubridate::minute(z) # 0-59  # minute
lubridate::hour(z)   # 0-23  # hour
lubridate::wday(z)   # 1-7   # week
lubridate::mday(z)   # 1-31  # month
lubridate::month(z)  # 1-12  # year # e.g. monthly composite images
lubridate::yday(z)   # 1-366 # year
