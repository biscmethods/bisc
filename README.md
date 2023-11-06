# biclust
Cell clustering wrapper for [scregclust](https://github.com/sven-nelander/scregclust)

## Installation

You can install the development version of biclust from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("sebfranz/biclust")
```

## Example

Load package

``` r
library(devtools)
load_all(".")
```

Then run demo/lm_example.R.

## Development

Load package

``` r
library(devtools)
load_all(".")
```

To generate documentation load all and

``` r
library(roxygen2) #  Read in the roxygen2 R package
roxygenise() #  Builds the help files
```

To check compliance of package load all and

``` r
check()
```