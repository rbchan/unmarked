\name{unmarkedFrameGDR}

\title{Organize data for the combined distance and removal point-count model of 
  Amundson et al. (2014) fit by \code{gdistremoval}}

\alias{unmarkedFrameGDR}
\alias{unmarkedFrameGDR-class}

\usage{
  unmarkedFrameGDR(yDistance, yRemoval, numPrimary=1, siteCovs=NULL, obsCovs=NULL, 
                   yearlySiteCovs=NULL, dist.breaks, unitsIn, period.lengths=NULL)
}

\description{Organize data for the combined distance and removal point-count model of 
  Amundson et al. (2014) fit by \code{gdistremoval}}

\arguments{
  \item{yDistance}{An MxTJ matrix of count data, where M is the number of sites
    (points), T is the number of primary periods (can be 1) and J is the number of 
    distance classes}
  \item{yRemoval}{An MxTJ matrix of count data, where M is the number of sites
    (points), T is the number of primary periods (can be 1) and J is the number of 
    time removal periods}
  \item{numPrimary}{Number of primary periods in the dataset}
  \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the
    site level. This should have M rows and one column per covariate}
  \item{obsCovs}{A \code{\link{data.frame}} of covariates that vary at the
    site level. This should have MxTJ rows and one column per covariate.
    These covariates are used only by the removal part of the model}
  \item{yearlySiteCovs}{A \code{\link{data.frame}} of covariates that vary
   by site and primary period. This should have MxT rows and one column per covariate}
  \item{dist.breaks}{vector of distance cut-points delimiting the
    distance classes. It must be of length J+1}
  \item{unitsIn}{Either "m" or "km" defining the measurement units for
    \code{dist.breaks}}
  \item{period.lengths}{Optional vector of time lengths of each removal period. 
    Each value in the vector must be a positive integer, and the total length
    of the vector must be equal to the number of removal periods J. If this is 
    not provided (the default), then all periods are assumed to have an equal 
    length of 1 time unit}
}

\details{
  unmarkedFrameGDR is the S4 class that holds data to be passed
  to the \code{\link{gdistremoval}} model-fitting function.}

\value{an object of class \code{unmarkedFrameGDR}}

\note{If you have continuous distance data, they must be "binned" into
  discrete distance classes, which are delimited by dist.breaks.}

\references{
  Amundson, C.L., Royle, J.A. and Handel, C.M., 2014. A hierarchical model 
  combining distance sampling and time removal to estimate detection probability 
  during avian point counts. The Auk 131: 476-494.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarkedFrame-class}}, \code{\link{unmarkedFrame}},
  \code{\link{gdistremoval}}}
