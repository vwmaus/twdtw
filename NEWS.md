# version 1.1-0, November 14, 2026

* New function `twdtw_knn()`. It calculates the TWDTW distance between a set of time series and a set of temporal patterns in a single call to the Fortran routines, and returns for each time series the k patterns with the lowest distance. The distances are the same as the distances calculated with `twdtw()`. In a benchmark with 10,000 time series (23 observations, 6 bands) and 5 temporal patterns, `twdtw_knn()` was about 140 times faster than calling `twdtw()` for each pair of data frames, and about 30 times faster than calling `twdtw()` for each pair of numeric matrices.

* `twdtw_knn()` calls a time-weight function defined in R once for each temporal pattern, and not once for each pair of observations, if the function is vectorised. In the same benchmark, this was about 25 times faster than a function that is not vectorised. The logistic time weight is the fastest option, about 7 times faster than a vectorised function.

* The function `to_date_time()` is now exported.

* Fixes an out-of-bounds write in `twdtw()` if a pattern matched a longer time series more often than the pattern has observations. This could return an incorrect distance, omit matches, or crash R.

* `twdtw()` returns `NA` if every observation of a time series has a missing value, so that no observation is left to compare. Before, it wrote out of bounds in this case and returned `Inf`, stopped with an error, or crashed R. A time series with some missing values is compared as before, using its observations without missing values.

* `twdtw()` stops with an error if a time-weight function is used with `version = 'f90goto'` or `version = 'cpp'`. These two versions do not call the function: they ignored it without a message and added a constant time weight of 0.5.

* The versions `'f90goto'` and `'cpp'` of `twdtw()` are deprecated and will be removed in version 1.2-0. They give the same result as the default version `'f90'`. `twdtw()` gives a warning once per session if one of them is used.

# version 1.0-1, August 7, 2023

* Minor fixes. Package usage and results are not affected.

# version 1.0-0, Jul 8, 2023

* Efficient generalized implementation of Time-Weighted Dynamic Time Warping (TWDTW)
