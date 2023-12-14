
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

Example on EVD:

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
#> 1 round(eigen(A, symmetric = TRUE)…   163ms 164.1ms      6.10    6.87MB     2.03
#> 2 round(eigen_IRAM(A, k)$eigenVal,…  96.4ms  97.3ms     10.1    13.27MB     6.76
#> 3 round(eigen_KSLE(A, k)$eigenVal,…  69.5ms  84.9ms     11.6     8.96MB     4.64
```

Example on SVD:

``` r
result_SVD <- bench::mark(
  round(svd(M)$d[1:k], 0),
  round(truncatedSVD(M, k, method = "IRAM")$D, 0),
  round(truncatedSVD(M, k, method = "KSLE")$D, 0),
)
result_SVD
#> # A tibble: 3 × 6
#>   expression                            min  median `itr/sec` mem_alloc `gc/sec`
#>   <bch:expr>                       <bch:tm> <bch:t>     <dbl> <bch:byt>    <dbl>
#> 1 "round(svd(M)$d[1:k], 0)"          2.99ms   3.1ms    312.      1.11MB     8.65
#> 2 "round(truncatedSVD(M, k, metho… 119.64ms 119.6ms      8.36   20.36MB    25.1 
#> 3 "round(truncatedSVD(M, k, metho… 118.51ms 118.5ms      8.44   14.05MB    25.3
```
