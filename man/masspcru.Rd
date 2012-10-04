\name{masspcru}
\alias{masspcru}
\docType{data}
\title{ Massachusetts North American Amphibian Monitoring Program Data }
\description{
\code{masspcru} contains NAAMP data for Pseudacris
crucifer (pcru) in Massachusetts from 2001 to 2007 in the raw long format.
}
\usage{data(masspcru)}
\format{
  Data frame with
  	\describe{
	  \item{SurveyYear}{Year of data collection.}
	  \item{RouteNumStopNum}{Stop number.}
	  \item{JulianDate}{Day of year.}
	  \item{Pcru}{Observed calling index.}
	  \item{MinAfterSunset}{Minutes after sunset of the observation.}
	  \item{Temperature}{Temperature measured during observation.}
	}

}
\details{
These data come from the North American Amphibian Monitoring Program.
Please see the reference below for more details.
}
\source{
https://www.pwrc.usgs.gov/naamp/
}
\references{
Mossman MJ, Weir LA. North American Amphibian Monitoring Program (NAAMP). Amphibian Declines: the conservation status of United States species. University of California Press, Berkeley, California, USA. 2005:307-313.
}
\examples{
data(masspcru)
str(masspcru)
}
\keyword{datasets}
