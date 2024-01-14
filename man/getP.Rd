\name{getP}
\alias{getP}
\alias{getP-methods}
\alias{getP,unmarkedFit-method}

\title{Get detection parameter estimates from a fitted 'unmarked' model}

\description{
The \code{getP} function returns a matrix of the back-transformed detection parameter (\eqn{p} the detection probability or \eqn{\lambda} the detection rate, depending on the model).
}

\usage{
\S4method{getP}{unmarkedFit}(object, na.rm = TRUE)
}

\arguments{
  \item{object}{An \code{unmarkedFit} object.}
  \item{na.rm}{Logical. Should missing values be removed?}
}

\value{A matrix of dimension MxJ, with M the number of sites and J the number of sampling periods; or of dimension MxJT for models with multiple primary periods T. Some fit types with multiple detection processes or multiple species return a named list of matrices, one per detection process or species.
} 
