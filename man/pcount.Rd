\name{pcount}

\alias{pcount}

\title{Fit the N-mixture model of Royle (2004)}

\usage{pcount(formula, data, K, mixture=c("P", "NB", "ZIP"),
    starts, method="BFGS", se=TRUE, engine=c("C", "R"), ...)}

\description{Fit the N-mixture model of Royle (2004)}

\arguments{
    \item{formula}{Double right-hand side formula describing covariates of
        detection and abundance, in that order}
    \item{data}{an unmarkedFramePCount object supplying data to the model.}
    \item{K}{Integer upper index of integration for N-mixture. This should be
        set high enough so that it does not affect the parameter estimates. Note
        that computation time will increase with K.}
    \item{mixture}{character specifying mixture: "P", "NB", or "ZIP".}
    \item{starts}{vector of starting values}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{logical specifying whether or not to compute standard
      errors.}
    \item{engine}{Either "C" or "R" to use fast C++ code or native R
      code during the optimization.}
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }


\details{This function fits N-mixture model of Royle (2004) to spatially replicated count data.

See \code{\link{unmarkedFramePCount}} for a description of how to format data
for \code{pcount}.

This function fits the latent N-mixture model for point count data
(Royle 2004, Kery et al 2005).

The latent abundance distribution, \eqn{f(N | \mathbf{\theta})}{f(N |
theta)} can be set as a Poisson, negative binomial, or zero-inflated
Poisson random
variable, depending on the setting of the \code{mixture} argument,
\code{mixture = "P"}, \code{mixture = "NB"}, \code{mixture = "ZIP"}
respectively.  For the first two distributions, the mean of \eqn{N_i} is
\eqn{\lambda_i}{lambda_i}.  If \eqn{N_i \sim NB}{N_i ~ NB}, then an
additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
\eqn{\alpha}{alpha} implies higher variance). For the ZIP distribution,
the mean is \eqn{\lambda_i(1-\psi)}{lambda*(1-psi)}, where psi is the
zero-inflation parameter.

The detection process is modeled as binomial: \eqn{y_{ij} \sim
Binomial(N_i, p_{ij})}{y_ij ~ Binomial(N_i, p_ij)}.

Covariates of \eqn{\lambda_i}{lamdba_i} use the log link and
covariates of \eqn{p_{ij}}{p_ij} use the logit link.}

\value{unmarkedFit object describing the model fit.}

\author{Ian Fiske and Richard Chandler}

\seealso{\code{\link{unmarkedFramePCount}}, \code{\link{pcountOpen}},
\code{\link{ranef}}, \code{\link{parboot}}
}

\references{

Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from
Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Kery, M., Royle, J. A., and Schmid, H. (2005) Modeling Avaian Abundance from
Replicated Counts Using Binomial Mixture Models. \emph{Ecological Applications}
15(4), pp. 1450--1461.

Johnson, N.L, A.W. Kemp, and S. Kotz. (2005) Univariate Discrete
Distributions, 3rd ed. Wiley.
}

\examples{

# Simulate data
set.seed(35)
nSites <- 100
nVisits <- 3
x <- rnorm(nSites)               # a covariate
beta0 <- 0
beta1 <- 1
lambda <- exp(beta0 + beta1*x)   # expected counts at each site
N <- rpois(nSites, lambda)       # latent abundance
y <- matrix(NA, nSites, nVisits)
p <- c(0.3, 0.6, 0.8)            # detection prob for each visit
for(j in 1:nVisits) {
  y[,j] <- rbinom(nSites, N, p[j])
  }

# Organize data
visitMat <- matrix(as.character(1:nVisits), nSites, nVisits, byrow=TRUE)

umf <- unmarkedFramePCount(y=y, siteCovs=data.frame(x=x),
    obsCovs=list(visit=visitMat))
summary(umf)

# Fit a model
fm1 <- pcount(~visit-1 ~ x, umf, K=50)
fm1

plogis(coef(fm1, type="det")) # Should be close to p


# Empirical Bayes estimation of random effects
(fm1re <- ranef(fm1))
plot(fm1re, subset=site \%in\% 1:25, xlim=c(-1,40))
sum(bup(fm1re))         # Estimated population size
sum(N)                  # Actual population size



\dontrun{

# Real data
data(mallard)
mallardUMF <- unmarkedFramePCount(mallard.y, siteCovs = mallard.site,
obsCovs = mallard.obs)
(fm.mallard <- pcount(~ ivel+ date + I(date^2) ~ length + elev + forest, mallardUMF, K=30))
(fm.mallard.nb <- pcount(~ date + I(date^2) ~ length + elev, mixture = "NB", mallardUMF, K=30))

}

}

\keyword{models}