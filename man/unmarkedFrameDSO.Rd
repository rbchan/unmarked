\name{unmarkedFrameDSO}
\title{Create an object of class unmarkedFrameDSO that contains data used by distsampOpen.}
\alias{unmarkedFrameDSO}
\usage{unmarkedFrameDSO(y, siteCovs=NULL, yearlySiteCovs=NULL, numPrimary, 
       primaryPeriod, dist.breaks, tlength, survey, unitsIn)}

\description{Organizes distance sampling data and experimental design information 
  from multiple primary periods along with associated covariates. This S4 class 
  is required by the data argument of \code{\link{distsampOpen}}}

\arguments{
    \item{y}{An MxJT matrix of the repeated count data, where M is the number
        of sites (i.e., points or transects), J is the number of distance classes
        and T is the maximum number of primary sampling periods per site}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the
        site level. This should have M rows and one column per covariate}
    \item{yearlySiteCovs}{Either a named list of MxT \code{\link{data.frame}}s,
        or a site-major \code{\link{data.frame}} with MT rows and 1 column per
        covariate}
    \item{numPrimary}{Maximum number of observed primary periods for each site}
    \item{primaryPeriod}{An MxJT matrix of integers indicating the primary period 
      of each observation}
    \item{dist.breaks}{vector of distance cut-points delimiting the
        distance classes. It must be of length J+1}
    \item{tlength}{A vector of length R containing the transect lengths. This is
        ignored when survey="point"}
    \item{survey}{Either "point" or "line" for point- and line-transects}
    \item{unitsIn}{Either "m" or "km" defining the measurement units for
        \emph{both} \code{dist.breaks} and \code{tlength}}
}
\details{
  \code{unmarkedFrameDSO} is the S4 class that holds data to be passed
  to the \code{\link{distsampOpen}} model-fitting function. Unlike
  most unmarked functions, \code{obsCovs} cannot be supplied.

  If you have continuous distance data, they must be "binned" into
  discrete distance classes, which are delimited by dist.breaks.

  When gamma or omega are modeled using year-specific covariates, the
  covariate data for the final year will be ignored; however,
  they must be supplied.

  If the time gap between primary periods is not constant, an M by T
  matrix of integers should be supplied using the \code{primaryPeriod} argument.
  
}

\value{an object of class unmarkedFrameDSO}
\seealso{\code{\link{unmarkedFrame-methods}}, \code{\link{distsampOpen}}}

\examples{

# Fake data
M <- 4 # number of sites
J <- 3 # number of distance classes
T <- 2 # number of primary periods

db <- c(0, 10, 20, 30) # distance break points

y <- matrix(c(
   5,4,3, 6,2,1, # In bin 1: 5 detections in primary period 1, 6 in period 2
   0,0,0, 0,1,0,
   2,1,1, 0,0,0,
   1,1,0, 1,1,1), nrow=M, ncol=J*T, byrow=TRUE)
y

# Primary periods of observations
# In this case there are no gaps
primPer <- matrix(as.integer(c(
    1,2,
    1,2,
    1,2,
    1,2)), nrow=M, ncol=T, byrow=TRUE)

#Site covs: M rows and 1 column per covariate
site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
site.covs

#Yearly site covs on gamma/omega
ysc <- list(
    x3 = matrix(c(
        1,2,
        1,2,
        1,2,
        1,2), nrow=M, ncol=T, byrow=TRUE))

umf <- unmarkedFrameDSO(y=y, siteCovs=site.covs, yearlySiteCovs=ysc,
                        numPrimary=T, primaryPeriod=primPer,
                        dist.breaks=db, survey="point", unitsIn="m")            

umf                         # look at data
summary(umf)                # summarize


}
