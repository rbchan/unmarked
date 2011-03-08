\name{multinomPois}

\alias{multinomPois}

\title{Multinomial-Poisson Mixtures Model}

\usage{multinomPois(formula, data, starts, method = "BFGS", control = list(), se = TRUE)}

\description{Fit the multinomial-Poisson mixture model to data collected using
    survey methods such as removal sampling or double observer sampling.}

\details{
This function takes advantage of the closed form of the integrated
likelihood when a latent Poisson distribution is assumed for abundance
at each site and a multinomial distribution is taken for the observation
state. Many common sampling methods can be framed in this context.  For
example, double-observer point counts and removal sampling
can be analyzed with this function by specifying the proper
multinomial cell probablilities.  This is done with by supplying the
appropriate function (piFun) argument.  \link{removalPiFun} and \link{doublePiFun} are supplied as example cell probability functions.}

\value{unmarkedFit object describing the model fit.}

\author{Ian Fiske}

\seealso{\link{piFuns}, \code{\link{unmarkedFrameMPois}}}

\references{

Royle, J. A. (2004). Generalized estimators of avian abundance from
count survey data. Animal Biodiversity and Conservation, 27(1), 375-386.

Royle, J. A., & Dorazio, R. M. (2006). Hierarchical Models of Animal Abundance and Occurrence. Journal Of Agricultural Biological And Environmental Statistics, 11(3), 249.}

\arguments{
    \item{formula}{double right-hand side formula for detection and abundance 
        covariates, in that order.}
    \item{data}{unmarkedFrame supplying data.}
    \item{starts}{vector of starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{control}{Other arguments passed to \code{\link{optim}}.}
    \item{se}{logical specifying whether or not to compute standard errors.}}

\examples{

data(ovendata)
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
    siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])), 
    type = "removal")
(fm1 <- multinomPois(~ 1 ~ ufp + trba, ovenFrame))

}

\keyword{models}
