\name{unmarkedFrameDS}

\title{Organize data for the distance sampling model of Royle et al. (2004)
    fit by distsamp}

\alias{unmarkedFrameDS}

\usage{unmarkedFrameDS(y, siteCovs=NULL, dist.breaks, tlength, survey,
    unitsIn, mapInfo)}

\description{Organizes count data along with the covariates and metadata.
    This S4 class is required by the data argument of \code{\link{distsamp}}}

\arguments{
    \item{y}{An RxJ matrix of count data, where R is the
        number of sites (transects) and J is the number of distance
    classes.}
    \item{siteCovs}{A \code{\link{data.frame}} of covariates that vary at the
        site level. This should have R rows and one column per covariate}
    \item{dist.breaks}{vector of distance cut-points delimiting the
        distance classes. It must be of length J+1.}
    \item{tlength}{A vector of length R containing the trasect lengths. This is
        ignored when survey="point".}
    \item{survey}{Either "point" or "line" for point- and line-transects.}
    \item{unitsIn}{Either "m" or "km" defining the measurement units for
        \emph{both} \code{dist.breaks} and \code{tlength}}.
    \item{mapInfo}{Currently ignored}}

\details{
    unmarkedFrameDS is the S4 class that holds data to be passed
    to the \code{\link{distsamp}} model-fitting function.}

  \value{an object of class unmarkedFrameDS}

  \note{If you have continuous distance data, they must be "binned" into
  discrete distance classes, which are delimited by dist.breaks.}

\references{

Royle, J. A., D. K. Dawson, and S. Bates (2004) Modeling
abundance effects in distance sampling. \emph{Ecology} 85, pp. 1591-1597.

}


\seealso{\code{\link{unmarkedFrame-class}}, \code{\link{unmarkedFrame}},
    \code{\link{distsamp}}}

\examples{

# Fake data
R <- 4 # number of sites
J <- 3 # number of distance classes

db <- c(0, 10, 20, 30) # distance break points

y <- matrix(c(
   5,4,3, # 5 detections in 0-10 distance class at this transect
   0,0,0,
   2,1,1,
   1,1,0), nrow=R, ncol=J, byrow=TRUE)
y

site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
site.covs

umf <- unmarkedFrameDS(y=y, siteCovs=site.covs, dist.breaks=db, survey="point",
    unitsIn="m")            # organize data
umf                         # look at data
summary(umf)                # summarize
fm <- distsamp(~1 ~1, umf)  # fit a model


}
