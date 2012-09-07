\name{issj}
\alias{issj}
\docType{data}
\title{
  Distance-sampling data for the Island Scrub Jay (\emph{Aphelocoma
    insularis})
}
\description{
  Data were collected at 307 survey locations ("point transects") on
  Santa Cruz Island, California during the Fall of 2008. The distance
  data are binned into 3 distance intervals
  [0-100], (100-200], and (200-300]. The coordinates of the survey
  locations as well as 3 habitat covariates are also included.
}
\usage{data(issj)}
\format{
  A data frame with 307 observations on the following 8 variables.
  \describe{
    \item{\code{issj[0-100]}}{Number of individuals detected within 100m}
    \item{\code{issj(100-200]}}{Detections in the interval (100-200m]}
    \item{\code{issj(200-300]}}{Detections in the interval (200-300m]}
    \item{\code{x}}{Easting (meters)}
    \item{\code{y}}{Northing (meters)}
    \item{\code{elevation}}{Elevation in meters}
    \item{\code{forest}}{Forest cover}
    \item{\code{chaparral}}{Chaparral cover}
  }
}
\seealso{Island-wide covariates are also available \code{\link{cruz}}}
\references{
Sillett, S. and Chandler, R.B. and Royle, J.A. and Kery, M. and
  Morrison, S.A. In Press. Hierarchical distance sampling models to
  estimate population size and habitat-specific abundance of an island
  endemic. \emph{Ecological Applications}
}
\examples{
data(issj)
str(issj)
head(issj)

umf <- unmarkedFrameDS(y=as.matrix(issj[,1:3]), siteCovs=issj[,6:8],
    dist.breaks=c(0,100,200,300), unitsIn="m", survey="point")
summary(umf)

}
\keyword{datasets}
