<!-- badges: start -->
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/vwmaus/twdtw/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vwmaus/twdtw/actions/workflows/R-CMD-check.yaml)
[![Coverage Status](https://img.shields.io/codecov/c/github/vwmaus/twdtw/main.svg)](https://app.codecov.io/gh/vwmaus/twdtw)
<!-- badges: end -->
  
# twdtw

Implements Time-Weighted Dynamic Time Warping (TWDTW), a measure for time series similarities. 
TWDTW is applicable to multi-dimensional time series with various resolutions, making it highly suited 
for environmental and ecological data analysis, such as remote sensing imagery, climate data, hydrology, 
or animal movement. The package provides a user-friendly R interface, efficient Fortran routines for TWDTW 
calculations, flexible time weighting definition, utilities for time series preprocessing, and visualization.

# References

Maus, V., Camara, G., Cartaxo, R., Sanchez, A., Ramos, F. M., & de Moura, Y. M. (2016).
A Time-Weighted Dynamic Time Warping Method for Land-Use and Land-Cover Mapping.
IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 9(8), 3729-3739.
\doi{10.1109/JSTARS.2016.2517118}

Maus, V., Camara, G., Appel, M., & Pebesma, E. (2019).
dtwSat: Time-Weighted Dynamic Time Warping for Satellite Image Time Series Analysis in R.
Journal of Statistical Software, 88(5), 1-31. \doi{10.18637/jss.v088.i05}
