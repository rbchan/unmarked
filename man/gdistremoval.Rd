\name{gdistremoval}
\alias{gdistremoval}
\title{
Fit the combined distance and removal model of Amundson et al. (2014).
}
\description{

Fit the model of Amundson et al. (2014) to point count datasets containing
both distance and time of observation data. The Amundson et al. (2014) model
is extended to account for temporary emigration by estimating an additional
availability probability if multiple counts at a site are available.
Abundance can be modeled as a Poisson, negative binomial, or Zero-inflated Poisson.
Multiple distance sampling key functions are also available.
}
\usage{
gdistremoval(lambdaformula=~1, phiformula=~1, removalformula=~1,
  distanceformula=~1, data, keyfun=c("halfnorm", "exp", "hazard", "uniform"),
  output=c("abund", "density"), unitsOut=c("ha", "kmsq"), mixture=c('P', 'NB', 'ZIP'), 
  K, starts, method = "BFGS", se = TRUE, engine=c("C","TMB"), threads=1, ...)
}

\arguments{
  \item{lambdaformula}{A right-hand side formula describing the abundance covariates}
  \item{phiformula}{A right-hand side formula describing the availability covariates}
  \item{removalformula}{A right-hand side formula describing removal probability covariates}
  \item{distanceformula}{A right-hand side formula describing the detection function covariates}
  \item{data}{An object of class \code{unmarkedFrameGDR}}
  \item{keyfun}{One of the following detection functions: "halfnorm", "hazard", "exp",
    or "uniform"}
  \item{output}{Model either "abund" or "density"}
  \item{unitsOut}{Units of density. Either "ha" or "kmsq" for hectares and square
    kilometers, respectively}
  \item{mixture}{Either "P", "NB", or "ZIP" for the Poisson, negative binomial,
    and Zero-inflated Poisson models of abundance}
  \item{K}{An integer value specifying the upper bound used in the integration}
  \item{starts}{A numeric vector of starting values for the model parameters}
  \item{method}{Optimization method used by \code{\link{optim}}}
  \item{se}{logical specifying whether or not to compute standard errors}
  \item{engine}{Either "C" to use C++ code or "TMB" to use TMB for optimization}
  \item{threads}{Set the number of threads to use for optimization in C++, if
      OpenMP is available on your system. Increasing the number of threads
      may speed up optimization in some cases by running the likelihood 
      calculation in parallel. If \code{threads=1} (the default), OpenMP is disabled} 
  \item{\dots}{Additional arguments to optim, such as lower and upper bounds}

}

\value{An object of class unmarkedFitGDR}

\references{
  Amundson, C.L., Royle, J.A. and Handel, C.M., 2014. A hierarchical model 
  combining distance sampling and time removal to estimate detection probability 
  during avian point counts. The Auk 131: 476-494.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarkedFrameGDR}}, \code{\link{gdistsamp}}, \code{\link{gmultmix}}}

\keyword{models}
