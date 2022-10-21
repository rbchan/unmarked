# R package unmarked

<!-- badges: start -->

[![R build
status](https://github.com/rbchan/unmarked/workflows/R-CMD-check/badge.svg)](https://github.com/rbchan/unmarked/actions)
[![CRAN
status](https://www.r-pkg.org/badges/version/unmarked)](https://cran.r-project.org/package=unmarked)
<!-- badges: end -->

`unmarked` is an [R](https://www.r-project.org/) package for analyzing
ecological data arising from several popular sampling techniques. The
sampling methods include point counts, occurrence sampling, distance
sampling, removal, double observer, and many others. `unmarked` uses
hierarchical models to incorporate covariates of the latent abundance
(or occupancy) and imperfect detection processes.

## Installation

The latest stable version of unmarked can be downloaded from
[CRAN](https://cran.r-project.org/package=unmarked):

``` r
install.packages("unmarked")
```

The latest development version can be installed from Github:

``` r
install.packages("remotes")
remotes::install_github("rbchan/unmarked")
```

## Support

Support is provided through the [unmarked Google
group](http://groups.google.com/group/unmarked). The package
[website](https://rbchan.github.io/unmarked) has more information. You
can report bugs [here](https://github.com/rbchan/unmarked/issues), by
posting to the Google group, or by emailing [the current
maintainer](https://kenkellner.com).

## Example analysis

Below we demonstrate a simple single-season occupancy analysis using
`unmarked`. First, load in a dataset from a CSV file and format:

``` r
library(unmarked)
wt <- read.csv(system.file("csv","widewt.csv", package="unmarked"))

# Presence/absence matrix
y <- wt[,2:4]

# Site and observation covariates
siteCovs <-  wt[,c("elev", "forest", "length")]
obsCovs <- list(date=wt[,c("date.1", "date.2", "date.3")]) 
```

Create an `unmarkedFrame`, a special type of `data.frame` for `unmarked`
analyses:

``` r
umf <- unmarkedFrameOccu(y = y, siteCovs = siteCovs, obsCovs = obsCovs)
summary(umf)
```

    ## unmarkedFrame Object
    ## 
    ## 237 sites
    ## Maximum number of observations per site: 3 
    ## Mean number of observations per site: 2.81 
    ## Sites with at least one detection: 79 
    ## 
    ## Tabulation of y observations:
    ##    0    1 <NA> 
    ##  483  182   46 
    ## 
    ## Site-level covariates:
    ##       elev               forest              length      
    ##  Min.   :-1.436125   Min.   :-1.265352   Min.   :0.1823  
    ##  1st Qu.:-0.940726   1st Qu.:-0.974355   1st Qu.:1.4351  
    ##  Median :-0.166666   Median :-0.064987   Median :1.6094  
    ##  Mean   : 0.007612   Mean   : 0.000088   Mean   :1.5924  
    ##  3rd Qu.: 0.994425   3rd Qu.: 0.808005   3rd Qu.:1.7750  
    ##  Max.   : 2.434177   Max.   : 2.299367   Max.   :2.2407  
    ## 
    ## Observation-level covariates:
    ##       date         
    ##  Min.   :-2.90434  
    ##  1st Qu.:-1.11862  
    ##  Median :-0.11862  
    ##  Mean   :-0.00022  
    ##  3rd Qu.: 1.30995  
    ##  Max.   : 3.80995  
    ##  NA's   :42

Fit a null occupancy model and a model with covariates, using the `occu`
function:

``` r
(mod_null <- occu(~1~1, data=umf))
```

    ## 
    ## Call:
    ## occu(formula = ~1 ~ 1, data = umf)
    ## 
    ## Occupancy:
    ##  Estimate    SE     z  P(>|z|)
    ##    -0.665 0.139 -4.77 1.82e-06
    ## 
    ## Detection:
    ##  Estimate    SE    z  P(>|z|)
    ##      1.32 0.174 7.61 2.82e-14
    ## 
    ## AIC: 528.987

``` r
(mod_covs <- occu(~date~elev, data=umf))
```

    ## 
    ## Call:
    ## occu(formula = ~date ~ elev, data = umf)
    ## 
    ## Occupancy:
    ##             Estimate    SE     z  P(>|z|)
    ## (Intercept)   -0.738 0.157 -4.71 2.45e-06
    ## elev           0.885 0.174  5.10 3.49e-07
    ## 
    ## Detection:
    ##             Estimate    SE     z  P(>|z|)
    ## (Intercept)   1.2380 0.180 6.869 6.47e-12
    ## date          0.0603 0.121 0.497 6.19e-01
    ## 
    ## AIC: 498.158

Rank them using AIC:

``` r
fl <- fitList(null=mod_null, covs=mod_covs)
modSel(fl)
```

    ##      nPars    AIC delta AICwt cumltvWt
    ## covs     4 498.16  0.00 1e+00     1.00
    ## null     2 528.99 30.83 2e-07     1.00

Estimate occupancy probability using the top-ranked model at the first
six sites:

``` r
head(predict(mod_covs, type='state'))
```

    ##   Predicted         SE      lower     upper
    ## 1 0.1448314 0.03337079 0.09080802 0.2231076
    ## 2 0.1499962 0.03351815 0.09535878 0.2280473
    ## 3 0.2864494 0.03346270 0.22555773 0.3562182
    ## 4 0.3035399 0.03371489 0.24175619 0.3733387
    ## 5 0.1607798 0.03374307 0.10502635 0.2382512
    ## 6 0.1842147 0.03392277 0.12669813 0.2600662

Predict occupancy probability at a new site with given covariate values:

``` r
nd <- data.frame(elev = 1.2)
predict(mod_covs, type="state", newdata=nd)
```

    ##   Predicted         SE     lower     upper
    ## 1 0.5803085 0.06026002 0.4598615 0.6918922
