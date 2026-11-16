<!-- badges: start -->
[![License](https://img.shields.io/badge/license-GPL%20%28%3E=%202%29-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html)
[![R-CMD-check](https://github.com/vwmaus/twdtw/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/vwmaus/twdtw/actions/workflows/R-CMD-check.yaml)
[![Coverage Status](https://img.shields.io/codecov/c/github/vwmaus/twdtw/main.svg)](https://app.codecov.io/gh/vwmaus/twdtw)
[![CRAN](https://www.r-pkg.org/badges/version/twdtw)](https://cran.r-project.org/package=twdtw)
[![Downloads](https://cranlogs.r-pkg.org/badges/twdtw?color=brightgreen)](https://www.r-pkg.org/pkg/twdtw)
[![total](https://cranlogs.r-pkg.org/badges/grand-total/twdtw)](https://www.r-pkg.org/pkg/twdtw)
[![cran checks](https://badges.cranchecks.info/worst/twdtw.svg)](https://cran.r-project.org/web/checks/check_results_twdtw.html)
[![status](https://tinyverse.netlify.com/badge/twdtw)](https://CRAN.R-project.org/package=twdtw)
<!-- badges: end -->
  
# twdtw

Provides a fast implementation of Time-Weighted Dynamic Time Warping (TWDTW), 
a measure for quantifying time series similarity. The TWDTW algorithm, 
described in Maus et al. (2016) and 
Maus et al. (2019), is applicable to multi-dimensional 
time series of various resolutions. It is particularly suitable for comparing 
time series with seasonality for environmental and ecological data analysis, 
covering domains such as remote sensing imagery, climate data, hydrology, 
and animal movement. The 'twdtw' package offers a user-friendly 'R' interface, 
efficient 'Fortran' routines for TWDTW calculations, a k-nearest neighbour model 
to classify time series, flexible time weighting definitions, as well as utilities 
for time series preprocessing and visualization. With one temporal pattern for 
each class, the classification of time series has a processing time comparable 
to that of other machine learning models available in 'R', such as random forest 
and support vector machines.

# Installation

```r
# Released version from CRAN
install.packages("twdtw")

# Development version from GitHub
# install.packages("remotes")
remotes::install_github("vwmaus/twdtw")
```

# Usage

The TWDTW distance between two time series:

```r
library(twdtw)

n <- 23
t <- seq(0, pi, length.out = n)
d <- seq(as.Date("2020-09-01"), length.out = n, by = "15 day")

x <- data.frame(time = d, v1 = sin(t) * 2 + runif(n))
y <- data.frame(time = d + 30, v1 = sin(t) * 2 + runif(n))

twdtw(x, y,
      cycle_length = "year",
      time_scale = "day",
      time_weight = c(steepness = 0.1, midpoint = 50))
```

Classification of time series with a k-nearest neighbour model, available from version 1.1-0:

```r
# Temporal patterns of 5 land use classes and satellite image time series of 999 pixels
ts_y <- readRDS(system.file("sits/modis_mod13q1_patterns.rds", package = "twdtw"))
ts_x <- readRDS(system.file("sits/modis_mod13q1_999ts.rds", package = "twdtw"))

# Create the model from the labelled time series
m <- twdtw_knn(x = ts_y$observations,
               y = ts_y$label,
               cycle_length = "year",
               time_scale = "day",
               time_weight = c(steepness = 0.1, midpoint = 50))

# Classify each time series with the label of its nearest temporal pattern
lu <- predict(m, newdata = ts_x$observations)

table(lu$prediction)
```

The help pages of `twdtw_knn()` and `predict.twdtw_knn()` describe models with several
temporal patterns for each label, the number of nearest patterns `k`, the votes, and the
types of prediction.

# References

Maus, V., Camara, G., Cartaxo, R., Sanchez, A., Ramos, F. M., & de Queiroz, G. R. (2016).
A Time-Weighted Dynamic Time Warping Method for Land-Use and Land-Cover Mapping.
IEEE Journal of Selected Topics in Applied Earth Observations and Remote Sensing, 9(8), 3729-3739.
[10.1109/JSTARS.2016.2517118](https://doi.org/10.1109/JSTARS.2016.2517118)

Maus, V., Camara, G., Appel, M., & Pebesma, E. (2019).
dtwSat: Time-Weighted Dynamic Time Warping for Satellite Image Time Series Analysis in R.
Journal of Statistical Software, 88(5), 1-31. [10.18637/jss.v088.i05](https://doi.org/10.18637/jss.v088.i05)
