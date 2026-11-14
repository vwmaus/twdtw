# twdtw-1.1-0

This version adds new features and fixes two memory errors. No function or argument
of version 1.0-1 was removed or changed.

## Changes

* New function `twdtw_knn()`. It calculates the TWDTW distance between a set of
  time series and a set of temporal patterns in a single call to the 'Fortran'
  routines, and returns for each time series the k patterns with the lowest
  distance. The distances are the same as the distances calculated with
  `twdtw()`. In our benchmark with 10,000 time series and 5 patterns it was
  about 140 times faster than calling `twdtw()` for each pair of data frames,
  and about 30 times faster than calling `twdtw()` for each pair of numeric
  matrices. A time-weight function defined in 'R' is called once for each
  temporal pattern, and not once for each pair of observations, if the
  function is vectorised. In the benchmark, this was about 25 times faster
  than calling the function for each pair of observations.

* The function `to_date_time()` is now exported. It was an internal function in
  version 1.0-1.

* Fixes an out-of-bounds write in `twdtw()`. The array that stores the matches
  had the length of the pattern (`y`), but there can be as many matches as
  observations in the time series (`x`). If a pattern matched a longer time
  series more often than the pattern has observations, the 'Fortran' and 'C++'
  routines wrote beyond the end of the array. This could return an incorrect
  distance, omit matches, or crash R. The array now has the length of the time
  series.

* Fixes a second out-of-bounds write in `twdtw()`. The function removes the
  observations with missing values and compares the remaining observations. If
  every observation of a time series had a missing value, no observation was
  left, but the 'Fortran' and 'C++' routines still read and wrote the first
  element of arrays without elements. `twdtw()` then returned `Inf`, stopped
  with an error, or crashed R. It now returns `NA` for the distance and no
  matches in this case, and does not call the compiled routines. Time series
  with some missing values are compared as before.

* `twdtw()` stops with an error if a time-weight function defined in 'R' is
  used with `version = 'f90goto'` or `version = 'cpp'`. Only the default
  version, 'f90', calls the function. The other two versions ignored it
  without a message and added a constant time weight of 0.5 to the distance of
  each pair of observations.

* The versions 'f90goto' and 'cpp' of `twdtw()` are deprecated. They give the
  same result as the default version 'f90' and will be removed in version
  1.2-0. `twdtw()` gives a warning once per session if one of them is used,
  and stops with an error for an unknown version.

* Adds tests for the new function and for the memory errors, and a larger test
  data set with 999 time series.

## Segmentation fault in the tests of 'dtwSat'

The first memory error fixed in this version caused the segmentation fault in the
tests of 'dtwSat' 1.0.0, which was archived on 2023-10-31
(<https://cran-archive.r-project.org/web/checks/2023/2023-10-31_check_results_dtwSat.html>).
The test that failed compares time series with 23 observations to patterns with
6 observations. We reproduced the segmentation fault with twdtw 1.0-1, and it
does not occur with this version. We plan to resubmit 'dtwSat' after this
update.

## Test environments

* Local Ubuntu 24.04.5 LTS, R 4.6.2, GCC and GNU Fortran 13.3.0
* GitHub Actions, Ubuntu 26.04.1 LTS, R-devel (2026-11-12 r90881)
* GitHub Actions, Ubuntu 26.04.1 LTS, R 4.6.2
* GitHub Actions, Ubuntu 26.04.1 LTS, R 4.5.3
* GitHub Actions, macOS Tahoe 26.7.2 (arm64), R 4.6.2
* GitHub Actions, Windows Server 2022, R 4.6.2

## R CMD check results

There were no ERRORs or WARNINGs.

There were no NOTEs on GitHub Actions. The local check with `--as-cran` returned
two NOTEs that are caused by the local installation of R and not by the
package: the compiler flag `-mno-omit-leaf-frame-pointer`, which is set in the
R configuration of the system, and 'HTML Tidy' not being installed.

## Reverse dependencies

There is one reverse dependency on CRAN, 'rasterdiv' 0.3.8, which imports
'twdtw' and calls `twdtw()` for pairs of time series with the same dates. The
interface of `twdtw()` did not change. In `R CMD check` with 'twdtw' 1.0-1 and
with this version, 'rasterdiv' installs and its examples run without error.

The check of 'rasterdiv' does not run its call to `twdtw()`, because its
examples, tests and vignettes do not evaluate code that uses it. We therefore
compared the results of this call between version 1.0-1 and this version for
4,300 random pairs of time series:

* Time series without missing values: the results are identical.

* Time series with some missing values: the results are identical, unless
  version 1.0-1 wrote out of bounds because of the first memory error. This
  occurred if there were more matches than observations without missing values
  in the second time series (317 pairs), and version 1.0-1 then returned an
  incorrect distance for 185 of them.

* Time series in which every observation has a missing value: version 1.0-1
  returned `Inf`, stopped with an error, or crashed R. This version returns
  `NA`.


# twdtw-1.0-1

* Fixes error: 
    error #7723: RESULT name must not be same as FUNCTION name.   [LOGISTIC_TW]
    function logistic_tw(DIST, TD, TW1, TW2) bind(C, name="logistic_tw_") 
    result(logistic_tw)

The f90 files compile locally using intel-basekit on Ubuntu 22.04.2 LTS.

## Test environments

* Local Ubuntu 22.04.2 LTS, R-4.3.1, GCC GNU Fortran
* macOS Monterey 12.6.7, R-4.3.1, GCC GNU Fortran
* Ubuntu 22.04.2 LTS, R-devel, GCC GNU Fortran
* Ubuntu 22.04.2 LTS, R-4.2.3, GCC GNU Fortran
* Windows Server 2022, R-4.3.1, GCC GNU Fortran
* Windows Server 2022, R-devel, GCC GNU Fortran
* Ubuntu Linux 20.04.1 LTS, R-release, GCC
* Fedora Linux, R-devel, clang, gfortran
* Debian Linux, R-devel, GCC ASAN/UBSAN

## R CMD check results

* Debian Linux, R-devel, GCC ASAN/UBSAN returned PREPERROR

The compilation error is caused by the -std=gnu++98 compiler flag. This flag tells g++ to use the GNU dialect of the C++98 standard, which predates C++11. As a result, the C++11 features (e.g., std::function, nullptr) used in my code aren't recognized in the system. All other test environments work.
