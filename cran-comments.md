This is a new release.

The twdtw package provides a generalized implementation of Time-Weighted Dynamic Time Warping (TWDTW) for any type of time series data.
Previously, TWDTW was only available for remote sensing time series from the package dtwSat.

## Test environnement

* Local Ubuntu 22.04, R-4.3.1, GCC GNU Fortran
* macOS Monterey 12.6.7, R-4.3.1, GCC GNU Fortran
* Ubuntu 22.04.2 LTS, R-devel, GCC GNU Fortran
* Ubuntu 22.04.2 LTS, R-4.2.3, GCC GNU Fortran
* Windows Server 2022, R-4.3.1, GCC GNU Fortran
* Windows Server 2022, R-devel, GCC GNU Fortran

## R CMD check results

* checking CRAN incoming feasibility ... NOTE
  Maintainer: 'Victor Maus <vwmaus1@gmail.com>'
  
  New submission

* Possibly misspelled words in DESCRIPTION: TWDTW (14:61, 14:109, 14:418)
  
The word is correctly spelled

* Debian Linux, R-devel, GCC ASAN/UBSAN returned PREPERROR

The compilation error is caused by the -std=gnu++98 compiler flag. This flag tells g++ to use the GNU dialect of the C++98 standard, which predates C++11. As a result, the C++11 features (e.g., std::function, nullptr) used in my code aren't recognized in the system.
