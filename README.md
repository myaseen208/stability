
## `stability`: Stability Analysis of Genotype by Environment Interaction (GEI)

###### Version : [0.5.0](https://myaseen208.com/stability/); Copyright (C) 2018-2024: License: [GPL-2\|GPL-3](https://www.r-project.org/Licenses/)

##### *Muhammad Yaseen<sup>1,2</sup>, and Kent M. Eskridge<sup>3</sup>*

1.  School of Mathematical & Statistical Sciences, Clemson University,
    Clemson, South Carolina, USA

2.  Department of Mathematics and Statistics, University of Agriculture
    Faisalabad, Pakistan

3.  Department of Statistics, University of Nebraska Lincoln, NE, USA

------------------------------------------------------------------------

[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.5.0-6666ff.svg)](https://cran.r-project.org/)
[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version-last-release/stability)](https://cran.r-project.org/package=stability)
[![rstudio mirror
downloads](https://cranlogs.r-pkg.org/badges/grand-total/stability?color=green)](https://CRAN.R-project.org/package=stability)
<!-- [![packageversion](https://img.shields.io/badge/Package%20version-0.2.3.3-orange.svg)](https://github.com/myaseen208/stability) -->

[![develVersion](https://img.shields.io/badge/devel%20version-0.5.0-orange.svg)](https://github.com/myaseen208/stability)

<!-- [![GitHub Download Count](https://github-basic-badges.herokuapp.com/downloads/myaseen208/stability/total.svg)] -->

[![Project Status:
WIP](https://www.repostatus.org/badges/latest/inactive.svg)](https://www.repostatus.org/#inactive)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--09--28-yellowgreen.svg)](https://github.com/myaseen208/stability)
\*\*\*

## Description

Functionalities to perform Stability Analysis of Genotype by Environment
Interaction (GEI) to identify superior and stable genotypes under
diverse environments. It performs Eberhart & Russel’s ANOVA (1966)
([doi:10.2135/cropsci1966.0011183X000600010011x\>), Finlay and Wilkinson
(1963) Joint Linear Regression (\<doi:10.1071/AR9630742\>), Wricke
(1962, 1964) Ecovalence, Shukla’s stability variance parameter (1972)
(\<doi:10.1038/hdy.1972.87\>) and Kang’s (1991)
(\<doi:10.2134/agronj1991.00021962008300010037x](https://doi.org/10.2135/cropsci1966.0011183X000600010011x%3E),
Finlay and Wilkinson (1963) Joint Linear Regression
(<doi:10.1071/AR9630742>), Wricke (1962, 1964) Ecovalence, Shukla’s
stability variance parameter (1972) (<doi:10.1038/hdy.1972.87>) and
Kang’s (1991) (\<<doi:10.2134/agronj1991.00021962008300010037x>))
simultaneous selection for high yielding and stable parameter.

## Installation

The package can be installed from CRAN as follows:

``` r
install.packages("stability", dependencies = TRUE)
```

The development version can be installed from github as follows:

``` r
if (!require("remotes")) install.packages("remotes")
remotes::install_github("myaseen208/stability")
```

## What’s new

To know whats new in this version type:

``` r
news(package = "stability")
```

## Links

[CRAN page](https://cran.r-project.org/package=stability)

[Github page](https://github.com/myaseen208/stability)

[Documentation website](https://myaseen208.com/stability/)

## Citing `stability`

To cite the methods in the package use:

``` r
citation("stability")
To cite package 'stability' in publications use:

  Muhammad Yaseen, Kent M. Eskridge, Ghulam Murtaza (2018). _stability:
  Stability Analysis of Genotype by Environment Interaction (GEI)_. R
  package version 0.5.0,
  <https://CRAN.R-project.org/package=stability>.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {stability: Stability Analysis of Genotype by Environment Interaction (GEI)},
    author = {{Muhammad Yaseen} and {Kent M. Eskridge} and {Ghulam Murtaza}},
    year = {2018},
    note = {R package version 0.5.0},
    url = {https://CRAN.R-project.org/package=stability},
  }
```
