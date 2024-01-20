\name{[-methods for unmarkedFrames}
\alias{[,unmarkedFrame,numeric,missing,missing-method}
\alias{[,unmarkedFrame,missing,numeric,missing-method}
\alias{[,unmarkedFrame,numeric,numeric,missing-method}
\alias{[,unmarkedFrame,list,missing,missing-method}
\alias{[,unmarkedFrame,logical,missing,missing-method}

\title{Methods for bracket extraction [ of unmarkedFrames}

\description{
  These methods extract subsets of sites, primary periods, or observations 
  from \code{unmarkedFrame} objects using numeric vectors of indices,
  a list of indices, or a logical vector.
}

\usage{
\S4method{[}{unmarkedFrame,numeric,missing,missing}(x, i, j)
\S4method{[}{unmarkedFrame,missing,numeric,missing}(x, i, j)
\S4method{[}{unmarkedFrame,numeric,numeric,missing}(x, i, j)
\S4method{[}{unmarkedFrame,list,missing,missing}(x, i, j)
\S4method{[}{unmarkedFrame,logical,missing,missing}(x, i, j)
}

\arguments{
	\item{x}{An object of class \code{unmarkedFrame}}
	\item{i}{Sites to select}
	\item{j}{Observation numbers (eg occasions, distance classes, etc.) to so select
    for single-season models, and primary periods to select for multi-season and 
    temporary emigration models}
}

\value{
  A new \code{unmarkedFrame} object of the class as the input.
}

\examples{

data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
	obsCovs = mallard.obs)
summary(mallardUMF)

mallardUMF[1:5,]
mallardUMF[,1:2]
mallardUMF[1:5, 1:2]

}
