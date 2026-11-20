# version 1.1-0, November 20, 2026

* Faster classification of time series. Version 1.1-0 is 20 to 140 times faster than version 1.0-1 on classification tasks with the new model `twdtw_knn()`. Assigning 10,000 time series (23 observations, 6 bands) to the nearest of 5 temporal patterns took 13 seconds in version 1.0-1. In version 1.1-0, it took 0.17 seconds starting from the same data frames, and 0.09 seconds for time series that are already in an array.

* New function `twdtw_knn()`. It creates a k-nearest neighbour model from labelled time series, and `predict()` classifies time series with the labels of their nearest temporal patterns. The time series can be a list of data frames or matrices, a data frame with one column for each band and date, or an array. The model keeps all labelled time series as temporal patterns, or creates one pattern for each label with a resampling function.

* A model with one pattern for each label keeps the patterns and the root mean squared error of the models fitted by the resampling function, which `print()` shows. The fitted models, e.g. the `gam` objects, contain the time series they were fitted to and are not kept by default; `keep_models = TRUE` keeps them in the model.

* `predict()` uses the `k` nearest temporal patterns of each time series, by default only the nearest. `vote` sets how their votes are counted: one vote for each pattern (`"majority"`, the default) or votes with the weight 1/distance (`"weighted"`). `type` sets the result: the label with the most votes (`"class"`, the default), the share of the votes of each label (`"prob"`), or the labels and the distances of the nearest patterns (`"neighbors"`).

* New function `twdtw_nearest()`. It calculates the TWDTW distance between a set of time series and a set of temporal patterns in a single call to the Fortran routines, and returns for each time series the k patterns with the lowest distance. The distances are the same as the distances calculated with `twdtw()`. `predict()` calls this function to find the nearest patterns.

* New generic function `twdtw_knn_predict()`, which `predict()` calls with the time series. Other packages can add methods to it to classify their own classes of time series.

* Details of the benchmark: version 1.0-1 was timed with a call to `twdtw()` for each time series and each pattern, which took 13 seconds with data frames and 2.8 seconds with numeric matrices. In version 1.1-0, `twdtw_knn()` and `predict()` gave the same labels and distances in 0.17 seconds for the list of data frames (about 75 times faster), in 0.14 seconds for the list of matrices (about 20 times faster), and in 0.09 seconds for an array (about 140 and 30 times faster). `twdtw()` has the same speed as in version 1.0-1.

* `predict()` and `twdtw_nearest()` call a time-weight function defined in R once for each temporal pattern, and not once for each pair of observations, if the function is vectorised. In the same benchmark, this was about 25 times faster than a function that is not vectorised. The logistic time weight is the fastest option, about 7 times faster than a vectorised function.

* The function `to_date_time()` is now exported.

* Fixes an out-of-bounds write in `twdtw()` if a pattern matched a longer time series more often than the pattern has observations. This could return an incorrect distance, omit matches, or crash R.

* `twdtw()` returns `NA` if every observation of a time series has a missing value, so that no observation is left to compare. Before, it wrote out of bounds in this case and returned `Inf`, stopped with an error, or crashed R. A time series with some missing values is compared as before, using its observations without missing values.

* `twdtw()` stops with an error if a time-weight function is used with `version = 'f90goto'` or `version = 'cpp'`. These two versions do not call the function: they ignored it without a message and added a constant time weight of 0.5.

* The versions `'f90goto'` and `'cpp'` of `twdtw()` are deprecated and will be removed in version 1.2-0. They give the same result as the default version `'f90'`. `twdtw()` gives a warning once per session if one of them is used.

# version 1.0-1, August 7, 2023

* Minor fixes. Package usage and results are not affected.

# version 1.0-0, Jul 8, 2023

* Efficient generalized implementation of Time-Weighted Dynamic Time Warping (TWDTW)
