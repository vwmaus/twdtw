# twdtw-1.1-0

This version adds new features and fixes a memory error. No function or argument
of version 1.0-1 was removed or changed.

## Changes

* New function `twdtw_knn()`. It calculates the TWDTW distance between a set of
  time series and a set of temporal patterns in a single call to the 'Fortran'
  routines, and returns for each time series the k patterns with the lowest
  distance. The distances are the same as the distances calculated with
  `twdtw()`. In our benchmark with 10,000 time series and 5 patterns it was
  about 120 times faster than calling `twdtw()` for each pair.

* The function `to_date_time()` is now exported. It was an internal function in
  version 1.0-1.

* Fixes an out-of-bounds write in `twdtw()`. The array that stores the matches
  had the length of the pattern (`y`), but there can be as many matches as
  observations in the time series (`x`). If a pattern matched a longer time
  series more often than the pattern has observations, the 'Fortran' and 'C++'
  routines wrote beyond the end of the array. This could return an incorrect
  distance, omit matches, or crash R. The array now has the length of the time
  series.

* Adds tests for the new function and for the memory error, and a larger test
  data set with 999 time series.

## Segmentation fault in the tests of 'dtwSat'

The memory error fixed in this version caused the segmentation fault in the
tests of 'dtwSat' 1.0.0, which was archived on 2023-10-31
(<https://cran-archive.r-project.org/web/checks/2023/2023-10-31_check_results_dtwSat.html>).
The test that failed compares time series with 23 observations to patterns with
6 observations. We reproduced the segmentation fault with twdtw 1.0-1, and it
does not occur with this version. We plan to resubmit 'dtwSat' after this
update.

## Test environments

* Local Ubuntu 24.04.5 LTS, R 4.6.2, GCC and GNU Fortran 13.3.0
* GitHub Actions, Ubuntu 26.04.1 LTS, R-devel (2026-11-11 r90876)
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

There is one reverse dependency on CRAN, 'rasterdiv', which imports 'twdtw' and
calls `twdtw()` for pairs of time series with the same dates. The interface of
`twdtw()` did not change, and for this kind of input the results of this
version are identical to the results of version 1.0-1.


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
