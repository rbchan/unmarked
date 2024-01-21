\name{fitted}
\alias{fitted}
\alias{fitted-methods}
\alias{fitted,unmarkedFit-method}

\title{Calculate fitted (expected) values from a model}

\description{Calculates expected values from a model. The approach varies based
on the model type but is usually the product of the state and detection parameters.}

\usage{
\S4method{fitted}{unmarkedFit}(object, na.rm = FALSE, ...)
}

\arguments{
  \item{object}{A \code{unmarkedFit} object}
  \item{na.rm}{Logical. Should missing values be removed from data? Not all model
               types support this.}
  \item{...}{Other arguments, currently ignored.}
}

\value{In most cases, \code{fitted} returns a matrix of expected values with the 
  same dimensions as the \code{y} matrix (M x J or M x JT). For some models with 
  multiple observation processes or multiple species, such as \code{gdistremoval}
  or \code{occuMulti}, the output may be a named list of matrices.
}
