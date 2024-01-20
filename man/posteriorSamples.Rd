\name{posteriorSamples}
\alias{posteriorSamples}
\alias{posteriorSamples-methods}
\alias{posteriorSamples,unmarkedRanef-method}
\alias{posteriorSamples,unmarkedFit-method}
\alias{unmarkedPostSamples-class}
\alias{[,unmarkedPostSamples,ANY,ANY,ANY-method}

\title{Draw samples from the posterior predictive distribution}

\description{
  Draw samples from the empirical Bayes posterior predictive distribution 
  derived from unmarked models or ranef objects
}

\usage{
\S4method{posteriorSamples}{unmarkedRanef}(object, nsims=100, ...)
\S4method{posteriorSamples}{unmarkedFit}(object, nsims=100, ...)
}

\arguments{
  \item{object}{An object inheriting class \code{unmarkedRanef} or 
    \code{unmarkedFit}}
  \item{nsims}{Number of draws to make from the posterior predictive distribution}
  \item{...}{Other arguments}
} 

\value{\code{unmarkedPostSamples} object containing the draws from the
 posterior predictive distribution. The draws are in the \code{@samples} slot. 
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{
  \code{\link{ranef}},
  \code{\link{predict}}
}

\examples{

# Simulate data under N-mixture model
set.seed(4564)
R <- 20
J <- 5
N <- rpois(R, 10)
y <- matrix(NA, R, J)
y[] <- rbinom(R*J, N, 0.5)

# Fit model
umf <- unmarkedFramePCount(y=y)
fm <- pcount(~1 ~1, umf, K=50)

# Estimates of conditional abundance distribution at each site
(re <- ranef(fm))

#Draw from the posterior predictive distribution
(ppd <- posteriorSamples(re, nsims=100))

}
