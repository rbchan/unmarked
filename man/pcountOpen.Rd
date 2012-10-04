\name{pcountOpen}
\alias{pcountOpen}
\title{Fit the open N-mixture model of Dail and Madsen}
\description{Fit the model of Dail and Madsen, which is 
a generalized form of the Royle (2004) N-mixture model for open populations.}
\usage{
pcountOpen(lambdaformula, gammaformula, omegaformula, pformula, data, 
	mixture = c("P", "NB"), K, dynamics=c("constant", "autoreg", "notrend"),
    fix=c("none", "gamma", "omega"), starts, method = "BFGS", se = TRUE, ...)
}
\arguments{
  \item{lambdaformula}{Right-hand sided formula for initial abundance}
  \item{gammaformula}{Right-hand sided formula for recruitment rate}
  \item{omegaformula}{Right-hand sided formula for apparent survival probability}
  \item{pformula}{Right-hand sided formula for detection probability}
  \item{data}{An object of class \code{\link{unmarkedFramePCO}}. See details}
  \item{mixture}{character specifying mixture: either "P" or "NB" for the Poisson and negative binomial distributions.}
  \item{K}{Integer defining upper bound of discrete integration. This should be higher than the maximum observed count and high enough that it does not affect the parameter estimates. However, the higher the value the slower the compuatation.}
  \item{dynamics}{Character string describing the type of population dynamics. "constant" indicates that there is no relationship between omega and gamma. "autoreg" is an auto-regressive model in which recruitment is modeled as gamma*N[i,t-1]. "notrend" model gamma as lambda*(1-omega) such that there is no temporal trend.}
  \item{fix}{If "omega", omega is fixed at 1. If "gamma", gamma is fixed at 0.}
  \item{starts}{vector of starting values}
  \item{method}{Optimization method used by \code{\link{optim}}.}
  \item{se}{logical specifying whether or not to compute standard errors.}
  \item{\dots}{additional arguments to be passed to \code{\link{optim}}.}
}
\details{
This model generalizes the Royle (2004) N-mixture model by relaxing the 
closure assumption. The model includes two additional parameters: gamma, 
the recruitment rate (births and immigrations), and omega, 
the apparent survival rate (deaths and emigrations). Estimates of population 
size at each time period can be derived from these parameters, and thus so can 
trend estimates.

The latent initial abundance distribution, 
\eqn{f(N | \mathbf{\theta})}{f(N | theta)} can be set as either a 
Poisson or a negative binomial random variable, depending on the setting of the 
\code{mixture} argument. \code{mixture = "P"} or \code{mixture = "NB"} select 
the Poisson or negative binomial distribution respectively.  The mean of 
\eqn{N_i} is \eqn{\lambda_i}{lambda_i}.  If \eqn{N_i \sim NB}{N_i ~ NB}, then an
additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
\eqn{\alpha}{alpha} implies higher variance).

The latent abundance state following the initial sampling period arises from a
Markovian process in which survivors are modeled as \eqn{S_it \sim 
Binomial(N_it-1, omega_it)}{S_it ~ Binomial(N_it-1, omega_it)}, and recruits
follow \eqn{G_it \sim Poisson(gamma_it)}{G_it ~ Poisson(gamma_it)}.

The detection process is modeled as binomial: \eqn{y_{it} \sim
Binomial(N_it, p_it)}{y_it ~ Binomial(N_it, p_it)}.

Covariates of \eqn{\lambda_i}{lamdba_i} use the log link and
covariates of \eqn{p_it}{p_it} use the logit link.
 
}
\value{
An object of class unmarkedFitPCO.
}
\references{
Royle, J. A. (2004) N-Mixture Models for Estimating Population Size from Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Dail, D. and L. Madsen (In press) Models for Estimating Abundance from Repeated Counts of an Open Metapopulation. \emph{Biometrics}.

}

\author{Richard Chandler \email{rchandler@usgs.gov}}

\note{
When gamma or omega are modeled using observation-level covariates, the covariate data for the final survey occasion will be ignored; however, they must be supplied.

If time intervals vary among primary sampling occasions, an M by T matrix of integers should be supplied to \code{\link{unmarkedFramePCO}} using the \code{primaryPeriod} argument.

Secondary sampling periods are optional.
}

\section{Warning}{This function can be extremely slow, especially if there are covariates of gamma or omega. Consider testing the timing on a small subset of the data, perhaps with se=FALSE. Finding the lowest value of K that does not affect estimates will also help with speed. }

\seealso{
\code{\link{pcount}, \link{unmarkedFramePCO}}
}

\examples{

## Simulation 
## No covariates, constant time intervals between primary periods, and 
## no secondary sampling periods

set.seed(3)
M <- 50
T <- 5
lambda <- 1
gamma <- 0.5
omega <- 0.8
p <- 0.7
y <- N <- matrix(NA, M, T)
S <- G <- matrix(NA, M, T-1)
N[,1] <- rpois(M, lambda)
for(t in 1:(T-1)) {
	S[,t] <- rbinom(M, N[,t], omega)
	G[,t] <- rpois(M, gamma)
	N[,t+1] <- S[,t] + G[,t]
	}
y[] <- rbinom(M*T, N, p)

                            
# Prepare data                               
umf <- unmarkedFramePCO(y = y, numPrimary=T)

summary(umf)


# Fit model and backtransform
(m1 <- pcountOpen(~1, ~1, ~1, ~1, umf, K=10))

(lam <- coef(backTransform(m1, "lambda"))) # or
lam <- exp(coef(m1, type="lambda"))
gam <- exp(coef(m1, type="gamma"))
om <- plogis(coef(m1, type="omega"))
p <- plogis(coef(m1, type="det"))


# Calculate population size at each time period
N.hat <- matrix(NA, M, T)
N.hat[,1] <- lam
for(t in 2:T) {
    N.hat[,t] <- om*N.hat[,t-1] + gam
    }

rbind(N=colSums(N), N.hat=colSums(N.hat))

}

\keyword{models}
