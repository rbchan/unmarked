\name{unmarkedFrameMPois}

\title{Organize data for the multinomial-Poisson mixture model of Royle (2004) 
    fit by multinomPois}

\alias{unmarkedFrameMPois}

\usage{unmarkedFrameMPois(y, siteCovs=NULL, obsCovs=NULL, type, obsToY, 
    mapInfo, piFun)}

\description{Organizes count data along with the covariates. 
    This S4 class is required by the data argument of 
    \code{\link{multinomPois}}}

\arguments{
    \item{y}{An RxJ matrix of count data, where R is the 
        number of sites (transects) and J is the maximum number of observations
        per site.}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the 
        site level. This should have R rows and one column per covariate}
    \item{obsCovs}{Either a named list of RxJ \code{\link{data.frame}}s or 
        a \code{data.frame} with RxJ rows and one column per covariate. 
        For the latter format, the covariates should be in site-major order.}
    \item{type}{Either "removal" or "double" for removal sampling or double 
        observer sampling. If this argument not specified, the user must
        provide an \code{obsToY} matrix. See details.}
    \item{obsToY}{A matrix describing the relationship between \code{obsCovs} 
        and \code{y}. This is necessary because under some sampling designs 
        the dimensions of \code{y} do not equal the dimensions of each 
        observation level covariate. For example, in double observer sampling
        there are 3 observations (seen only by observer A, detected only by 
        observer B, and detected by both), but each observation-level covariate
        can only have 2 columns, one for each observer. This matrix is created 
        automatically if \code{type} is either "removal" or "double". } 
    \item{mapInfo}{Currently ignored}
    \item{piFun}{Function used to compute the multinomial cell probabilities 
        from a matrix of detection probabilities. This is created automatically
        if \code{type} is either "removal" or "double".}}

\details{
    unmarkedFrameMPois is the S4 class that holds data to be passed 
    to the \code{\link{multinomPois}} model-fitting function.}

\value{an object of class unmarkedFrameMPois}

\references{

Royle, J. A. (2004). Generalized estimators of avian abundance from
count survey data. Animal Biodiversity and Conservation, 27(1), 375-386.

}


\seealso{\code{\link{unmarkedFrame-class}}, \code{\link{unmarkedFrame}}, 
    \code{\link{multinomPois}}, \code{\link{piFuns}}}

\examples{

# Fake doulbe observer data
R <- 4 # number of sites
J <- 2 # number of observers

y <- matrix(c(
   1,0,3,
   0,0,0,
   2,0,1,
   0,0,2), nrow=R, ncol=J+1, byrow=TRUE)
y

site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
site.covs

obs.covs <- list(
   x3 = matrix(c(
      -1,0,
      -2,0,
      -3,1,
       0,0), 
      nrow=R, ncol=J, byrow=TRUE),
   x4 = matrix(c(
      'a','b',
      'a','b',
      'a','b',
      'a','b'), 
      nrow=R, ncol=J, byrow=TRUE))
obs.covs


# Create unmarkedFrame
umf <- unmarkedFrameMPois(y=y, siteCovs=site.covs, obsCovs=obs.covs,
    type="double")
    
# The above is the same as:
o2y <- matrix(1, 2, 3)
pifun <- function(p)
{
    M <- nrow(p)
    pi <- matrix(NA, M, 3)
    pi[, 1] <- p[, 1] * (1 - p[, 2])
    pi[, 2] <- p[, 2] * (1 - p[, 1])
    pi[, 3] <- p[, 1] * p[, 2]
    return(pi)
}

umf <- unmarkedFrameMPois(y=y, siteCovs=site.covs, obsCovs=obs.covs,
    obsToY=o2y, piFun="pifun")


# Fit a model
fm <- multinomPois(~1 ~1, umf)



}
