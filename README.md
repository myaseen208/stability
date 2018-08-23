# stability: Stability Analysis of Genotype by Environment Interaction (GEI)

## Introduction
 The **stability** package provides functionalities to perform
 Stability Analysis of Genotype by Environment Interaction (GEI)
 to identify superior and stable genotypes under diverse environments.
 It performs  Eberhart & Russel's ANOVA (1966) (https://dl.sciencesocieties.org/publications/cs/abstracts/6/1/CS0060010036),
 Finlay and Wilkinson (1963) Joint Linear Regression (http://www.publish.csiro.au/cp/AR9630742),
 Wricke (1962, 1964) Ecovalence, Shukla's stability variance parameter (1972) (https://www.nature.com/articles/hdy197287)
 and  Kang's (1991) simultaneous selection for high yielding and stable parameter (https://dl.sciencesocieties.org/publications/aj/abstracts/83/1/AJ0830010161).

## Authors
1. Muhammad Yaseen (\email{myaseen208@gmail.com}
2. Kent M. Edkridge (\email{keskridge1@unl.edu}

## Installation
Use **devtools** to install the development version from Github:

```{r}
if(!require("devtools")) install.packages("devtools")
devtools::install_github('MYaseen208/stability', build_vignettes = TRUE)
```

## License
This package is free and open source software, licensed under GPL.
