\name{unmarkedFrameOccu}

\title{Organize data for the single season occupancy models fit by occu and 
    occuRN}

\alias{unmarkedFrameOccu}

\usage{unmarkedFrameOccu(y, siteCovs=NULL, obsCovs=NULL, mapInfo)}

\description{Organizes detection, non-detection data along with the covariates. 
    This S4 class is required by the data argument of \code{\link{occu}} and
    \code{\link{occuRN}}}

\arguments{
    \item{y}{An RxJ matrix of the detection, non-detection data, where R is the 
        number of sites, J is the maximum number of sampling periods per site.}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the 
        site level. This should have M rows and one column per covariate}
    \item{obsCovs}{Either a named list of \code{\link{data.frame}}s of 
        covariates that vary within sites, or a \code{\link{data.frame}} with 
        RxJ rows in site-major order.}
    \item{mapInfo}{Currently ignored}
  
}

\details{
    unmarkedFrameOccu is the S4 class that holds data to be passed 
    to the \code{\link{occu}} and \code{\link{occuRN}} model-fitting function.
}

\value{an object of class unmarkedFrameOccu}

\seealso{\code{\link{unmarkedFrame-class}}, \code{\link{unmarkedFrame}}, 
    \code{\link{occu}}, \code{\link{occuRN}}}

\examples{

# Fake data
R <- 4 # number of sites
J <- 3 # number of visits
y <- matrix(c(
   1,1,0,
   0,0,0,
   1,1,1,
   1,0,1), nrow=R, ncol=J, byrow=TRUE)
y

site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
site.covs

obs.covs <- list(
   x3 = matrix(c(
      -1,0,1,
      -2,0,0,
      -3,1,0,
      0,0,0), nrow=R, ncol=J, byrow=TRUE),
   x4 = matrix(c(
      'a','b','c',
      'd','b','a',
      'a','a','c',
      'a','b','a'), nrow=R, ncol=J, byrow=TRUE))
obs.covs

umf <- unmarkedFrameOccu(y=y, siteCovs=site.covs, 
    obsCovs=obs.covs)   # organize data
umf                     # look at data
summary(umf)            # summarize      
fm <- occu(~1 ~1, umf)  # fit a model


}
