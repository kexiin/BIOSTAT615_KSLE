
<!-- README.md is generated from README.Rmd. Please edit that file -->

# KSLE

<!-- badges: start -->

[![Author](https://img.shields.io/badge/Author-KexinLi-red.svg "Author")](https://kexiin.gitee.io "Author")
<!-- badges: end -->

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
m <- 500  # matrix size
n <- 50
k <- 2    # number of eigenvalues to calculate
# Some random data
M <- matrix(rnorm(n * m), n, m)
# Make it symmetric
A <- crossprod(M)
```

``` r
result_eigen <- bench::mark(
  round(eigen(A, symmetric = TRUE)$values[1:k], 0),
  round(eigen_IRAM(A, k)$eigenVal, 0),
  round(eigen_KSLE(A, k)$eigenVal, 0)
)
result_eigen
#> # A tibble: 3 × 6
#>   expression                            min  median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                        <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl>
#> 1 round(eigen(A, symmetric = TRUE)… 156.3ms 159.3ms      6.32    6.87MB     2.11
#> 2 round(eigen_IRAM(A, k)$eigenVal,…  84.8ms    87ms     11.1    13.27MB    11.1 
#> 3 round(eigen_KSLE(A, k)$eigenVal,…  61.2ms  72.2ms     13.4     8.96MB     2.24
```

``` r
result_SVD <- bench::mark(
  round(svd(M)$d[1:k], 0),
  round(truncatedSVD(M, k, method = "IRAM")$D, 0),
  round(truncatedSVD(M, k, method = "KSLE")$D, 0),
)
result_SVD
#> # A tibble: 3 × 6
#>   expression                            min  median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                        <bch:t> <bch:t>     <dbl> <bch:byt>    <dbl>
#> 1 "round(svd(M)$d[1:k], 0)"          2.85ms  3.04ms     319.     1.11MB     8.81
#> 2 "round(truncatedSVD(M, k, method… 88.87ms 88.87ms      11.3   19.27MB    45.0 
#> 3 "round(truncatedSVD(M, k, method… 90.62ms 97.15ms      10.3   12.91MB    15.4
```
