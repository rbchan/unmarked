\name{simulate-methods}
\docType{methods}
\alias{simulate-methods}
\alias{simulate,unmarkedFitColExt-method}
\alias{simulate,unmarkedFitDS-method}
\alias{simulate,unmarkedFitMPois-method}
\alias{simulate,unmarkedFitOccu-method}
\alias{simulate,unmarkedFitOccuRN-method}
\alias{simulate,unmarkedFitPCount-method}
\alias{simulate,unmarkedFitPCO-method}
\alias{simulate,unmarkedFitGMM-method}
\alias{simulate,unmarkedFitGDS-method}
\title{Methods for Function simulate in Package `unmarked'}
\description{
Simulate data from a fitted model.
}
\usage{
\S4method{simulate}{unmarkedFitColExt}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitDS}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitMPois}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitOccu}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitOccuRN}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitPCount}(object, nsim, seed, na.rm)
}
\arguments{
\item{object}{Fitted model of appropriate S4 class}
\item{nsim}{Number of simulations}
\item{seed}{Seed for random number generator. Not currently implemented}
\item{na.rm}{Logical, should missing values be removed?}
}
\section{Methods}{
\describe{
\item{object = "unmarkedFitColExt"}{A model fit by \code{\link{colext}}}
\item{object = "unmarkedFitDS"}{A model fit by \code{\link{distsamp}}}
\item{object = "unmarkedFitMPois"}{A model fit by \code{\link{multinomPois}}}
\item{object = "unmarkedFitOccu"}{A model fit by \code{\link{occu}}}
\item{object = "unmarkedFitOccuRN"}{A model fit by \code{\link{occuRN}}}
\item{object = "unmarkedFitPCount"}{A model fit by \code{\link{pcount}}}
}}
\keyword{methods}

