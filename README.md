
<!-- README.md is generated from README.Rmd. Please edit that file -->

# KSLE

<!-- badges: start -->
<!-- badges: end -->

The goal of KSLE is to …

## Installation

You can install the development version of KSLE from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("kexiin/BIOSTAT615_KSLE")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(KSLE)

set.seed(123)
n <- 1000  # matrix size
# Some random data
M <- matrix(rnorm(n^2), n)
# Make it symmetric
A <- crossprod(M)
```

``` r
k <- 5    # number of eigenvalues to calculate
system.time({eigen(A)})
#>    user  system elapsed 
#>   1.172   0.014   1.192
system.time({eigen_IRAM(A, k)})
#>    user  system elapsed 
#>   1.224   0.036   1.267
system.time({eigen_KSLE(A, k)})
#>    user  system elapsed 
#>   1.282   0.006   1.313
```

``` r
k <- 10    # number of eigenvalues to calculate
system.time({eigen_IRAM(A, k)})
#>    user  system elapsed 
#>   7.081   0.067   7.189
system.time({eigen_KSLE(A, k)})
#>    user  system elapsed 
#>   3.833   0.022   3.877
```
