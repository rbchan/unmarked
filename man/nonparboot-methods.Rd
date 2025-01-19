\name{nonparboot-methods}
\docType{methods}
\alias{nonparboot}
\alias{nonparboot-methods}
\alias{nonparboot,unmarkedFit-method}
\alias{nonparboot,unmarkedFitColExt-method}
\alias{nonparboot,unmarkedFitDS-method}
\alias{nonparboot,unmarkedFitMPois-method}
\alias{nonparboot,unmarkedFitOccu-method}
\alias{nonparboot,unmarkedFitOccuPEN-method}
\alias{nonparboot,unmarkedFitOccuPEN_CV-method}
\alias{nonparboot,unmarkedFitOccuRN-method}
\alias{nonparboot,unmarkedFitPCount-method}
\alias{nonparboot,unmarkedFitGDS-method}
\alias{nonparboot,unmarkedFitGMM-method}
\alias{nonparboot,unmarkedFitOccuTTD-method}
\alias{nonparboot,unmarkedFitOccuMulti-method}
\alias{nonparboot,unmarkedFitNmixTTD-method}
\alias{nonparboot,unmarkedFitGDR-method}
\alias{nonparboot,unmarkedFitIDS-method}
\alias{nonparboot,unmarkedFitDailMadsen-method}
\alias{nonparboot,unmarkedFitOccuCOP-method}
\alias{nonparboot,unmarkedFitOccuRNMulti-method}

\title{ Nonparametric bootstrapping in unmarked }
\description{
Call \code{nonparboot} on an unmarkedFit to obtain non-parametric
bootstrap samples.  These can then be used by \code{vcov} in order to
get bootstrap estimates of standard errors.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "unmarkedFit")}}{ Obtain nonparametric
  bootstrap samples for a general unmarkedFit.  }

\item{\code{signature(object = "unmarkedFitColExt")}}{ Obtain nonparametric
  bootstrap samples for colext fits.  }

\item{\code{signature(object = "unmarkedFitDS")}}{ Obtain nonparametric
  bootstrap samples for a distsamp fits. }

\item{\code{signature(object = "unmarkedFitMPois")}}{ Obtain nonparametric
  bootstrap samples for a distsamp fits. }

\item{\code{signature(object = "unmarkedFitOccu")}}{  Obtain nonparametric
  bootstrap samples for a occu fits. }

\item{\code{signature(object = "unmarkedFitOccuPEN")}}{  Obtain nonparametric
  bootstrap samples for an occuPEN fit. }

\item{\code{signature(object = "unmarkedFitOccuPEN_CV")}}{  Obtain nonparametric
  bootstrap samples for occuPEN_CV fit. }

\item{\code{signature(object = "unmarkedFitOccuRN")}}{  Obtain nonparametric
  bootstrap samples for a occuRN fits. }

\item{\code{signature(object = "unmarkedFitPCount")}}{  Obtain nonparametric
  bootstrap samples for a pcount fits. }
}}
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
\keyword{methods}

