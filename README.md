
<!-- README.md is generated from README.Rmd. Please edit that file -->

# bayesRCM

<!-- badges: start -->
<!-- badges: end -->

The goal of bayesRCM is to …

## Installation

You can install the development version of bayesRCM like so:

``` r
devtools::install_github("nevilleq/bayesRCM")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(bayesRCM)
## basic example code

data_list <- bayesRCM::sim_data()
model     <- bayesRCM::rcm(y = data_list)
```
