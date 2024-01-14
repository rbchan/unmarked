\name{nonparboot}
\alias{nonparboot}
\alias{nonparboot-methods}
\alias{nonparboot,unmarkedFit-method}

\title{Nonparametric bootstrapping in unmarked}

\description{
Call \code{nonparboot} on an \code{unmarkedFit} to obtain non-parametric
bootstrap samples.  These can then be used by \code{vcov} in order to
get bootstrap estimates of standard errors. The function can be run
repeatedly on an \code{unmarkedFit} object to add additional bootstrap samples.
}

\usage{
\S4method{nonparboot}{unmarkedFit}(object, B = 0, keepOldSamples = TRUE, ...)
}

\arguments{
  \item{object}{A \code{unmarkedFit} object}
  \item{B}{Number of bootstrap samples to generate}
  \item{keepOldSamples}{Logical. Should previously generated bootstrap samples
        be kept?}
  \item{...}{Other arguments, currently ignored} 
}

\value{An \code{unmarkedFit} object with added bootstrap samples.}

\details{

  Calling \code{nonparboot} on an unmarkedFit returns the original
  unmarkedFit, with the bootstrap samples added on.  Then subsequent
  calls to \code{\link{vcov}} with the argument
  \code{method="nonparboot"} will use these bootstrap samples.
  Additionally, standard errors of derived estimates from either
  \code{\link{linearComb}} or \code{\link{backTransform}} can be
  instructed to use bootstrap samples by providing the argument
  \code{method = "nonparboot"}.

  For \code{\link{occu}} and \code{\link{occuRN}} both sites and
  occassions are re-sampled. For all other fitting functions, only sites
  are re-sampled.

}

\examples{
data(ovendata)
ovenFrame <- unmarkedFrameMPois(ovendata.list$data,
siteCovs=as.data.frame(scale(ovendata.list$covariates[,-1])), type = "removal")
(fm <- multinomPois(~ 1 ~ ufc + trba, ovenFrame))
fm <- nonparboot(fm, B = 20) # should use larger B in real life.
vcov(fm, method = "hessian")
vcov(fm, method = "nonparboot")
avg.abundance <- backTransform(linearComb(fm, type = "state", coefficients = c(1, 0, 0)))

## Bootstrap sample information propagates through to derived quantities.
vcov(avg.abundance, method = "hessian")
vcov(avg.abundance, method = "nonparboot")
SE(avg.abundance, method = "nonparboot")
}
