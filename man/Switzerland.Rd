\name{Switzerland}
\alias{Switzerland}
\docType{data}
\title{
Swiss landscape data
}
\description{
Spatially-referenced data on elevation, forest cover, and water at a
1km-sq resolution.
}
\usage{data(Switzerland)}
\format{
  A data frame with 42275 observations on the following 5 variables.
  \describe{
    \item{\code{x}}{Easting (m)}
    \item{\code{y}}{Northing (m)}
    \item{\code{elevation}}{a numeric vector (m)}
    \item{\code{forest}}{a numeric vector (percent cover)}
    \item{\code{water}}{a numeric vector (percent cover)}
  }
}
\details{
Forest and water coverage (in percent area) was computed using the
1992-97 landcover dataset of the Swiss Federal Statistical Office
(http://www.bfs.admin.ch). Median elevation (in metres) was computed
using a median aggregation of the digital elevation model of the Swiss
Federal Statistical Office.

x and y are the coordinates of the center of each 1km2 pixel.

The coordinate reference system is CH1903/LV03 (EPSG number: 21781).

These data can only be used for non-profit projects. Otherwise, written
permission must be obtained from the Swiss Federal Statistical Office
}
\source{
Swiss Federal Statistical Office (http://www.bfs.admin.ch)
}
\examples{
data(Switzerland)
str(Switzerland)

levelplot(elevation ~ x + y, Switzerland, aspect="iso",
    col.regions=terrain.colors(100))

\dontrun{
library(raster)
el.r <- rasterFromXYZ(Switzerland[,c("x","y","elevation")], crs =
"+proj=somerc +lat_0=46.95240555555556 +lon_0=7.439583333333333
+k_0=1 +x_0=600000 +y_0=200000 +ellps=bessel
+towgs84=674.374,15.056,405.346,0,0,0,0 +units=m +no_defs")
plot(el.r)
spplot(el.r)
}

}
\keyword{datasets}
