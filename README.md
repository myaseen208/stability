
## `stability`: Stability Analysis of Genotype by Environment Interaction (GEI)

###### Version : [0.7.0](https://myaseen208.com/stability/); Copyright (C) 2018-2024: License: [GPL-2\|GPL-3](https://www.r-project.org/Licenses/)

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
[![Last-changedate](https://img.shields.io/badge/last%20change-2024--09--29-yellowgreen.svg)](https://github.com/myaseen208/stability)
\*\*\*

## Description

Provides functionalities for performing stability analysis of genotype
by environment interaction (GEI) to identify superior and stable
genotypes across diverse environments. It implements Eberhart and
Russell’s ANOVA method
(1966)([doi:10.2135/cropsci1966.0011183X000600010011x\>), Finlay and
Wilkinson’s Joint Linear Regression method (1963)
(\<doi:10.1071/AR9630742\>), Wricke’s Ecovalence (1962, 1964), Shukla’s
stability variance parameter (1972) (\<doi:10.1038/hdy.1972.87\>),
Kang’s simultaneous selection for high yield and stability (1991)
(\<doi:10.2134/agronj1991.00021962008300010037x](https://doi.org/10.2135/cropsci1966.0011183X000600010011x%3E),
Finlay and Wilkinson’s Joint Linear Regression method (1963)
(<doi:10.1071/AR9630742>), Wricke’s Ecovalence (1962, 1964), Shukla’s
stability variance parameter (1972) (<doi:10.1038/hdy.1972.87>), Kang’s
simultaneous selection for high yield and stability (1991)
(\<<doi:10.2134/agronj1991.00021962008300010037x>)), Additive Main
Effects and Multiplicative Interaction (AMMI) method and Genotype plus
Genotypes by Environment (GGE) Interaction methods.

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
Please, support this project by citing it in your publications!

  Yaseen M, Eskridge KM (2018). _stability: Stability Analysis of
  Genotype by Environment Interaction (GEI)_.

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {stability: Stability Analysis of Genotype by Environment Interaction (GEI)},
    author = {Muhammad Yaseen and Kent M. Eskridge},
    year = {2018},
    journal = {The Comprehensive R Archive Network (CRAN)},
  }
```
