\name{SSE}
\alias{SSE}
\alias{SSE-methods}
\alias{SSE,unmarkedFit-method}
\alias{SSE,unmarkedFitOccuMulti-method}
\alias{SSE,unmarkedFitGDR-method}
\alias{SSE,unmarkedFitIDS-method}
\alias{SSE,unmarkedFitOccuRNMulti-method}
\title{Compute Sum of Squared Residuals for a Model Fit.}
\description{
Compute the sum of squared residuals for an unmarked fit object.  This
is useful for a \code{\link{parboot}}.
}
\usage{
SSE(fit, ...)
}
\arguments{
  \item{fit}{An unmarked fit object.}
  \item{...}{Additional arguments to be passed to statistic}
}
\value{
A numeric value for the models SSE.}
\seealso{
\code{\link{parboot}}
}


