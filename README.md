<!-- badges: start -->
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/vwmaus/twdtw/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vwmaus/twdtw/actions/workflows/R-CMD-check.yaml)
[![Coverage Status](https://img.shields.io/codecov/c/github/vwmaus/twdtw/main.svg)](https://app.codecov.io/gh/vwmaus/twdtw)
[![CRAN](https://www.r-pkg.org/badges/version/twdtw)](https://cran.r-project.org/package=twdtw)
[![Downloads](https://cranlogs.r-pkg.org/badges/twdtw?color=brightgreen)](https://www.r-pkg.org/pkg/twdtw)
<!-- badges: end -->
  
# twdtw

Implements Time-Weighted Dynamic Time Warping (TWDTW), 
a measure for quantifying time series similarity. The TWDTW algorithm, 
described in Maus et al. (2016) and 
Maus et al. (2019), is applicable to multi-dimensional 
time series of various resolutions. It is particularly suitable for comparing 
time series with seasonality for environmental and ecological data analysis, 
covering domains such as remote sensing imagery, climate data, hydrology, 
and animal movement. The 'twdtw' package offers a user-friendly 'R' interface, 
efficient 'Fortran' routines for TWDTW calculations, flexible time weighting 
definitions, as well as utilities for time series preprocessing and visualization.

# References

Maus, V., Camara, G., Cartaxo, R., Sanchez, A., Ramos, F. M., & de Moura, Y. M. (2016).
A Time-Weighted Dynamic Time Warping Method for Land-Use and Land-Cover Mapping.
IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 9(8), 3729-3739.
[10.1109/JSTARS.2016.2517118](https://doi.org/10.1109/JSTARS.2016.2517118)

Maus, V., Camara, G., Appel, M., & Pebesma, E. (2019).
dtwSat: Time-Weighted Dynamic Time Warping for Satellite Image Time Series Analysis in R.
Journal of Statistical Software, 88(5), 1-31. [10.18637/jss.v088.i05](https://doi.org/10.18637/jss.v088.i05)
