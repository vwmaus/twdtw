# version 1.1-0, November 13, 2026

* New function `twdtw_knn()`. It calculates the TWDTW distance between a set of time series and a set of temporal patterns in a single call to the Fortran routines, and returns for each time series the k patterns with the lowest distance. The distances are the same as the distances calculated with `twdtw()`, and the calculation is much faster than calling `twdtw()` for each pair.

* `twdtw_knn()` calls a time-weight function defined in R once for each temporal pattern, and not once for each pair of observations, if the function is vectorised.

* The function `to_date_time()` is now exported.

* Fixes an out-of-bounds write in `twdtw()` if a pattern matched a longer time series more often than the pattern has observations. This could return an incorrect distance, omit matches, or crash R.

* `twdtw()` returns `NA` if every observation of a time series has a missing value, so that no observation is left to compare. Before, it wrote out of bounds in this case, which could crash R. A time series with some missing values is compared as before, using its observations without missing values.

# version 1.0-1, August 7, 2023

* Minor fixes. Package usage and results are not affected.

# version 1.0-0, Jul 8, 2023

* Efficient generalized implementation of Time-Weighted Dynamic Time Warping (TWDTW)
