\name{unmarkedFramePCount}

\title{Organize data for the N-mixture model fit by pcount}
\alias{unmarkedFramePCount}

\usage{unmarkedFramePCount(y, siteCovs=NULL, obsCovs=NULL, mapInfo)}

\description{Organizes repeated count data along with the covariates. 
    This S4 class is required by the data argument of \code{\link{pcount}}}

\arguments{
    \item{y}{An RxJ matrix of the repeated count data, where R is the number 
        of sites, J is the maximum number of sampling periods per site.}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the 
        site level. This should have M rows and one column per covariate}
    \item{obsCovs}{Either a named list of \code{\link{data.frame}}s of 
        covariates that vary within sites, or a \code{\link{data.frame}} with 
        RxJ rows in site-major order.}
    \item{mapInfo}{Currently ignored}}

\details{
    unmarkedFramePCount is the S4 class that holds data to be passed 
    to the \code{\link{pcount}} model-fitting function.}

\value{an object of class unmarkedFramePCount}

\seealso{\code{\link{unmarkedFrame-class}}, \code{\link{unmarkedFrame}}, 
    \code{\link{pcount}}}

\examples{

# Fake data
R <- 4 # number of sites
J <- 3 # number of visits
y <- matrix(c(
   1,2,0,
   0,0,0,
   1,1,1,
   2,2,1), nrow=R, ncol=J, byrow=TRUE)
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

umf <- unmarkedFramePCount(y=y, siteCovs=site.covs, 
    obsCovs=obs.covs)          # organize data
umf                            # take a l
summary(umf)                   # summarize data
fm <- pcount(~1 ~1, umf, K=10) # fit a model
      

}
