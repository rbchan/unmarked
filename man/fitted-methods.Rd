\name{fitted-methods}
\docType{methods}
\alias{fitted-methods}
%%\alias{fitted,unmarkedFit-method}
%%\alias{fitted,unmarkedFitColExt-method}
%%\alias{fitted,unmarkedFitOccu-method}
%%\alias{fitted,unmarkedFitOccuRN-method}
%%\alias{fitted,unmarkedFitPCount-method}
\title{Methods for Function fitted in Package `unmarked'}
\description{Extracted fitted values from a fitted model.
}
\usage{
\S4method{fitted}{unmarkedFit}(object, na.rm = FALSE)
\S4method{fitted}{unmarkedFitColExt}(object, na.rm = FALSE)
\S4method{fitted}{unmarkedFitOccu}(object, na.rm = FALSE)
\S4method{fitted}{unmarkedFitOccuRN}(object, K, na.rm = FALSE)
\S4method{fitted}{unmarkedFitPCount}(object, K, na.rm = FALSE)
}
\arguments{
\item{object}{A fitted model of appropriate S4 class}
\item{K}{Integer specifying upper bound of integration.}
\item{na.rm}{Logical. Should missing values be removed from data?}
}
\section{Methods}{
\describe{
\item{object = "unmarkedFit"}{A fitted model}
\item{object = "unmarkedFitColExt"}{A model fit by \code{\link{colext}}}
\item{object = "unmarkedFitOccu"}{A model fit by \code{\link{occu}}}
\item{object = "unmarkedFitOccuRN"}{A model fit by \code{\link{occuRN}}}
\item{object = "unmarkedFitPCount"}{A model fit by \code{\link{pcount}}}
}}
\value{Returns a matrix of expected values}
\keyword{methods}
