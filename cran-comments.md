# twdtw-1.0-1

* Release requested to fix error in the Fortran code https://www.stats.ox.ac.uk/pub/bdr/Intel/twdtw.log

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
