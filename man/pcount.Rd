\name{pcount}

\alias{pcount}

\title{Fit the N-mixture model of Royle (2004)}

\usage{pcount(formula, data, K, mixture=c("P", "NB"), starts, method="BFGS",
    control=list(), se=TRUE)}

\description{Fit the N-mixture model of Royle (2004)}

\arguments{
    \item{formula}{Double right-hand side formula describing covariates of 
        detection and abundance, in that order}
    \item{data}{an unmarkedFramePCount object supplying data to the model.}
    \item{K}{Integer upper index of integration for N-mixture. This should be 
        set high enough so that it does not affect the parameter estimates. Note
        that computation time will increase with K.}
    \item{mixture}{character specifying mixture: either "P" or "NB".}
    \item{starts}{vector of starting values}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{control}{Other arguments passed to \code{\link{optim}}.}
    \item{se}{logical specifying whether or not to compute standard errors.}}
    

\details{This function fits N-mixture model of Royle (2004) to spatially replicated count data.

See \code{\link{unmarkedFramePCount}} for a description of how to format data 
for \code{pcount}.

This function fits the latent N-mixture model for point count data
(Royle 2004, Kéry et al 2005).

The latent abundance distribution, \eqn{f(N | \mathbf{\theta})}{f(N |
theta)} can be set as either a Poisson or a negative binomial random
variable, depending on the setting of the \code{mixture} argument.
\code{mixture = "P"} or \code{mixture = "NB"} select the Poisson or
negative binomial distribution respectively.  The mean of \eqn{N_i} is
\eqn{\lambda_i}{lambda_i}.  If \eqn{N_i \sim NB}{N_i ~ NB}, then an
additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
\eqn{\alpha}{alpha} implies higher variance).

The detection process is modeled as binomial: \eqn{y_{ij} \sim
Binomial(N_i, p_{ij})}{y_ij ~ Binomial(N_i, p_ij)}.

Covariates of \eqn{\lambda_i}{lamdba_i} use the log link and
covariates of \eqn{p_{ij}}{p_ij} use the logit link.}

\value{unmarkedFit object describing the model fit.}

\author{Ian Fiske \email{ianfiske@gmail.com}}

\seealso{\code{\link{unmarkedFramePCount}}, \code{\link{pcountOpen}}}

\references{

Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from 
Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Kéry, M., Royle, J. A., and Schmid, H. (2005) Modeling Avaian Abundance from 
Replicated Counts Using Binomial Mixture Models. \emph{Ecological Applications} 
15(4), pp. 1450--1461.
}

\examples{

data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
obsCovs = mallard.obs)
(fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest, mallardUMF, K=30))
(fm.mallard.nb <- pcount(~ date + I(date^2) ~ length + elev, mixture = "NB", mallardUMF, K=30))

}

\keyword{models}