\name{ranef-methods}
\docType{methods}
\alias{ranef}
\alias{ranef-methods}
\alias{ranef,unmarkedFitOccu-method}
\alias{ranef,unmarkedFitOccuRN-method}
\alias{ranef,unmarkedFitPCount-method}
\alias{ranef,unmarkedFitMPois-method}
\alias{ranef,unmarkedFitDS-method}
\alias{ranef,unmarkedFitGMM-method}
\alias{ranef,unmarkedFitGDS-method}
\alias{ranef,unmarkedFitGMMorGDS-method}
\alias{ranef,unmarkedFitColExt-method}
\alias{ranef,unmarkedFitPCO-method}
\title{ Methods for Function \code{ranef} in Package \pkg{unmarked} }
\description{
Estimate posterior distributions of the random variables (latent abundance or
occurrence) using empirical Bayes methods. These methods return an
object storing the posterior distributions of the latent variables at
each site, and for each year (primary period) in the case of open
population models. See \link{unmarkedRanef-class} for methods used to
manipulate the returned object.
}
\section{Methods}{
\describe{

\item{\code{signature(object = "unmarkedFitOccu")}}{Computes the
  conditional distribution of occurrence given the data and the
  estimates of the
  fixed effects, \eqn{Pr(z_i=1 | y_{ij}, \hat{\psi}_i,
  \hat{p}_{ij})}{Pr(z(i)=1 | y(i,j), psi(i), p(i,j))}}
\item{\code{signature(object = "unmarkedFitOccuRN")}}{Computes the
  conditional abundance distribution given the data and the estimates of
  the fixed effects, \eqn{Pr(N_i=k | y_{ij}, \hat{\psi}_i, \hat{r}_{ij})
  k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,j), psi(i), r(i,j)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitPCount")}}{\eqn{Pr(N_i=k |
  y_{ij}, \hat{\lambda}_i, \hat{p}_{ij}) k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,j), lambda(i), p(i,j)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitMPois")}}{\eqn{Pr(N_i=k |
  y_{ij}, \hat{\lambda}_i, \hat{p}_{ij}) k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,j), lambda(i), p(i,j)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitDS")}}{\eqn{Pr(N_i=k |
  y_{i,1:J}, \hat{\lambda}_i, \hat{\sigma}_{i}) k =
  0,1,\dots,K}{Pr(N(i)=k |
  y(i,1:J), lambda(i), sigma(i)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitGMM")}}{\eqn{Pr(M_i=k |
  y_{i,1:J,t}, \hat{\lambda}_i, \hat{\phi}_{it}, \hat{p}_{ijt})  k =
  0,1,\dots,K}{Pr(N(i)=k |
  y(i,1:J,t), lambda(i), phi(i,t), p(i,j,t)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitGDS")}}{\eqn{Pr(M_i=k |
  y_{i,1:J,t}, \hat{\lambda}_i, \hat{\phi}_{it}, \hat{\sigma}_{it})
  k = 0,1,\dots,K}{Pr(N(i)=k |
  y(i,1:J,t), lambda(i), phi(i,t), sigma(i,t)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitColExt")}}{\eqn{Pr(z_{it}=1 |
  y_{ijt}, \hat{\psi}_i, \hat{\gamma}_{it}, \hat{\epsilon}_{it},
  \hat{p}_{ijt})}{Pr(z(i,t)=1 |
  y(i,j,t), psi(i), gamma(i,t), epsilon(i,t), p(i,j,t)) for k = 0,1,...,K}}
\item{\code{signature(object = "unmarkedFitPCO")}}{\eqn{Pr(N_{it}=k |
  y_{ijt}, \hat{\lambda}_i, \hat{\gamma}_{it}, \hat{\omega}_{it},
  \hat{p}_{ijt}) k = 0,1,...,K}{Pr(N(i,t)=k |
  y(i,j,t), lambda(i), gamma(i,t), omega(i,t), p(i,j,t)) for k =
  0,1,...,K}}

}
}

\note{
  From Carlin and Louis (1996): \dQuote{... the Bayesian approach to
  inference depends on a prior distribution for the model
  parameters. This prior can depend on unknown parameters which in turn
  may follow some second-stage prior. This sequence of parameters and
  priors consitutes a hierarchical model. The hierarchy must stop at
  some point, with all remaining prior parameters assumed known. Rather
  than make this assumption, the basic empirical Bayes approach uses the
  observed data to estimate these final stage parameters (or to estimate
  the Bayes rule), and proceeds as in a standard Bayesian analysis.}
  }


\section{Warning}{
  Empirical Bayes methods can underestimate the variance of the
  posterior distribution because they do not account for uncertainty in
  the hyperparameters (lambda or psi). However, we have not found this
  to be the case in a set of limited simulations studies that can be
  found in the inst/unitTests directory of unmarked. Simulations do,
  however, indicate that the posterior mode can exhibit slight (3-5
  percent) negative bias as a point
  estimator of site-specific abundance. It appears to be safer to use
  the posterior mean even though this will not be an integer in general.
}

\author{Richard Chandler \email{rchandler@usgs.gov}}


\references{
  Laird, N.M. and T.A. Louis. 1987. Empirical Bayes confidence intervals
  based on bootstrap samples. Journal of the American Statistical
  Association 82:739--750.

  Carlin, B.P and T.A Louis. 1996. Bayes and Empirical Bayes Methods for
  Data Analysis. Chapman and Hall/CRC.

  Royle, J.A and R.M. Dorazio. 2008. Hierarchical Modeling and Inference
  in Ecology. Academic Press.
  }

\seealso{
  \link{unmarkedRanef-class}
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
# Best Unbiased Predictors
bup(re, stat="mean")           # Posterior mean
bup(re, stat="mode")           # Posterior mode
confint(re, level=0.9) # 90\% CI

# Plots
plot(re, subset=site \%in\% c(1:10), layout=c(5, 2), xlim=c(-1,20))

# Compare estimates to truth
sum(N)
sum(bup(re))
colSums(confint(re))

# Extract all values in convenient formats
post.df <- as(re, "data.frame")
head(post.df)
post.arr <- as(re, "array")


}
\keyword{methods}

