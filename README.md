# bisc
Bisc is an abbrevation for "BI(clust) SC(regclust)". It's an iterative CEM-based (Classification Expectation Maximization) cell clusterer that wraps a method clustering genes called [scregclust](https://github.com/scmethods/scregclust).

# WIP
This is current not working as library you can install.

## Installation

You can install the development version of biclust from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("biscmethods/bisc")
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
