\name{cruz}
\alias{cruz}
\docType{data}
\title{
  Landscape data for Santa Cruz Island
}
\description{
  Spatially-referenced elevation, forest cover, and vegetation data for
  Santa Cruz Island.
}
\usage{data(cruz)}
\format{
  A data frame with 2787 observations on the following 5 variables.
  \describe{
    \item{\code{x}}{Easting (meters)}
    \item{\code{y}}{Northing (meters)}
    \item{\code{elevation}}{a numeric vector, meters}
    \item{\code{forest}}{a numeric vector, proportion cover}
    \item{\code{chaparral}}{a numeric vector, proportion cover}
  }
}
\details{
  The resolution is 300x300 meters.

  The Coordinate system is EPSG number 26911

  NAD_1983_UTM_Zone_11N
  Projection: Transverse_Mercator
  False_Easting: 500000.000000
  False_Northing: 0.000000
  Central_Meridian: -117.000000
  Scale_Factor: 0.999600
  Latitude_Of_Origin: 0.000000
  Linear Unit: Meter
  GCS_North_American_1983
  Datum: D_North_American_1983

}
\source{
  Brian Cohen of the Nature Conservancy helped prepare the data
}
\references{
Sillett, S. and Chandler, R.B. and Royle, J.A. and Kery, M. and
  Morrison, S.A. In Press. Hierarchical distance sampling models to
  estimate population size and habitat-specific abundance of an island
  endemic. \emph{Ecological Applications}
}
\examples{
data(cruz)
str(cruz)

levelplot(elevation ~ x + y, cruz, aspect="iso",
    col.regions=terrain.colors(100))

if(require(raster)) {
elev <- rasterFromXYZ(cruz[,1:3],
     crs="+proj=utm +zone=11 +ellps=GRS80 +datum=NAD83 +units=m +no_defs")
elev
plot(elev)
}



}
\keyword{datasets}
