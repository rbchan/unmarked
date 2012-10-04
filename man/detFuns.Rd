\name{detFuns}
\title{Distance-sampling detection functions and associated density functions}
\usage{
gxhn(x, sigma)
gxexp(x, rate)
gxhaz(x, shape, scale)

dxhn(x, sigma)
dxexp(x, rate)
dxhaz(x, shape, scale)
drhn(r, sigma)
drexp(r, rate)
drhaz(r, shape, scale)
}
\alias{detFuns}
\alias{gxhn}
\alias{gxexp}
\alias{gxhaz}
\alias{grhn}
\alias{grexp}
\alias{grhaz}
\alias{dxhn}
\alias{dxexp}
\alias{dxhaz}
\alias{drhn}
\alias{drexp}
\alias{drhaz}
\description{These functions represent the currently available detection
    functions used for modeling line and point transect data with
	\code{\link{distsamp}. Detection functions begin with "g", and density
    functions begin with a "d".}
	}
\arguments{
	\item{x}{Perpendicular distance}
	\item{r}{Radial distance}
	\item{sigma}{Shape parameter of half-normal detection function}
	\item{rate}{Shape parameter of negative-exponential detection function}
	\item{shape}{Shape parameter of hazard-rate detection function}
	\item{scale}{Scale parameter of hazard-rate detection function}
	}
\seealso{\code{\link{distsamp}} for example of using these for plotting
detection function}
\examples{
# Detection probabilities at 25m for range of half-normal sigma values.
round(gxhn(25, 10:15), 2)

# Plot negative exponential distributions
plot(function(x) gxexp(x, rate=10), 0, 50, xlab="distance",
    ylab="Detection probability")
plot(function(x) gxexp(x, rate=20), 0, 50, add=TRUE, lty=2)
plot(function(x) gxexp(x, rate=30), 0, 50, add=TRUE, lty=3)

# Plot half-normal probability density functions for line- and point-transects
par(mfrow=c(2, 1))
plot(function(x) dxhn(x, 20), 0, 50, xlab="distance",
    ylab="Probability density", main="Line-transect")
plot(function(x) drhn(x, 20), 0, 50, xlab="distance",
    ylab="Probability density", main="Point-transect")
}