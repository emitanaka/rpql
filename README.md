
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rpql: Regularized PQL for Joint Selection in GLMMs

<!-- badges: start -->

[![CRAN
status](https://www.r-pkg.org/badges/version/rpql)](https://cran.r-project.org/package=rpql)
[![Lifecycle:
maturing](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
<!-- badges: end -->

`rpql` offers fast joint selection of fixed and random effects in
Generalized Linear Mixed Model (GLMMs) via ularization. The penalized
quasi-likelihood (PQL) is used as a loss ction, and penalties are added
on to perform fixed and random effects ection. This method of joint
selection in GLMMs, referred to regularized, is fast compared to
information criterion and hypothesis testing i et al., 2016).

Please note `rpql` is the core workshops function that performed
regularized PQL on a single set of tuning parameters. `rpqlseq` is a
wrapper to permit a sequence of tuning parameter values. The latter is
often what users want to use.

## Installation

You can install the released version of rpql from
[CRAN](https://CRAN.R-project.org)
with:

``` r
install.packages("rpql")
```

<!-- And the development version from [GitHub](https://github.com/) with: -->

<!-- ``` r -->

<!-- # install.packages("devtools") -->

<!-- devtools::install_github("emitanaka/rpql") -->

<!-- ``` -->
