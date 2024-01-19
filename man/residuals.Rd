\name{residuals}
\alias{residuals-methods}
\alias{residuals,unmarkedFit-method}

\title{Calculate residuals from a model}

\description{Calculates residuals from a model. The approach varies based
on the model type, but is usually the observed data (\code{y}) minus the 
fitted values. See \code{\link{fitted}}.}

\usage{
\S4method{residuals}{unmarkedFit}(object, ...)
}

\arguments{
  \item{object}{A \code{unmarkedFit} object}
  \item{...}{Other arguments, currently unused.}
}

\value{In most cases, \code{residuals} returns a matrix of residual values with the 
  same dimensions as the \code{y} matrix (M x J or M x JT). For some models with 
  multiple observation processes or multiple species, such as \code{gdistremoval}
  or \code{occuMulti}, the output may be a named list of matrices.
}
