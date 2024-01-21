\name{unmarkedFrame-methods}
\alias{unmarkedFrame-methods}
\alias{numSites}
\alias{numSites,unmarkedFrame-method}
\alias{numY}
\alias{numY,unmarkedFrame-method}
\alias{obsNum}
\alias{obsNum,unmarkedFrame-method}
\alias{plot,unmarkedFrame,missing-method}
\alias{summary,unmarkedFrame-method}

\title{Methods for unmarkedFrames}
\description{Methods for summarizing and viewing unmarkedFrames}

\usage{
\S4method{numSites}{unmarkedFrame}(object)
\S4method{numY}{unmarkedFrame}(object)
\S4method{obsNum}{unmarkedFrame}(object)
\S4method{plot}{unmarkedFrame,missing}(x, y, panels = 1, colorkey, strip = FALSE,
  ylab = "Site", xlab = "Observation", ...)
\S4method{summary}{unmarkedFrame}(object, ...)
}

\arguments{
	\item{object}{An object of class \code{unmarkedFrame}}
	\item{x}{An object of class \code{unmarkedFrame}}
	\item{y}{Unused}
  \item{panels}{How many panels plotted data should be split over}
	\item{colorkey}{Colors to use in plot, see \code{\link[lattice]{levelplot}}}
  \item{strip}{See \code{\link[lattice]{levelplot}}}
  \item{ylab}{Y-axis label for plot}
  \item{xlab}{X-axis label for plot}
  \item{...}{Other arguments, currently unused}
}

\section{Methods}{
  \describe{
    \item{numSites:}{Number of sites in the \code{unmarkedFrame}}
    \item{numY:}{Number of columns in the y-matrix of the \code{unmarkedFrame}}
    \item{obsNum:}{Total number of observations per site (not necessarily equal to
      \code{numY} depending on the \code{unmarkedFrame} type)}
    \item{plot:}{Summary plot of \code{unmarkedFrame}}
    \item{summary:}{Summary statistics for \code{unmarkedFrame}}
  }
}

\seealso{\code{\link{unmarkedFrame-[-methods}}, \code{\link{unmarkedFrame-extractors}}}

\examples{

# Organize data for pcount()
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
	obsCovs = mallard.obs)

mallardUMF

# Vizualize it
plot(mallardUMF)

# Summary stats
summary(mallardUMF)

numSites(mallardUMF)
numY(mallardUMF)
obsNum(mallardUMF)
}
