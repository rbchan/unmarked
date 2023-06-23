\name{nmixTTD}
\alias{nmixTTD}
\title{Fit N-mixture Time-to-detection Models}

\usage{nmixTTD(stateformula= ~1, detformula = ~1, data, K=100,
    mixture = c("P","NB","ZIP"), ttdDist = c("exp", "weibull"), starts, method="BFGS", 
    se=TRUE, engine = c("C", "R"), threads = 1, ...)}

\arguments{
  \item{stateformula}{Right-hand sided formula for the abundance at each site.}
  \item{detformula}{Right-hand sided formula for mean time-to-detection.}
  \item{data}{\code{unmarkedFrameOccuTTD} object that supplies the data
    (see \code{\link{unmarkedFrameOccuTTD}}). Note that only single-season
    models are supported by \code{nmixTTD}.}
  \item{K}{The upper summation index used to numerically integrate out the
    latent abundance. This should be set high enough so that it does not
    affect the parameter estimates. Computation time will increase with K.}
  \item{mixture}{String specifying mixture distribution: "P" for Poisson, 
    "NB" for negative binomial, "ZIP" for zero-inflated Poisson.}
  \item{ttdDist}{Distribution to use for time-to-detection; either
    \code{"exp"} for the exponential, or \code{"weibull"} for the Weibull,
    which adds an additional shape parameter \eqn{k}.}
  \item{starts}{optionally, initial values for parameters in the optimization.}
  \item{method}{Optimization method used by \code{\link{optim}}.}
  \item{se}{logical specifying whether or not to compute standard errors.}
  \item{engine}{Either "C" or "R" to use fast C++ code or native R
    code during the optimization.}
  \item{threads}{Set the number of threads to use for optimization in C++, if
      OpenMP is available on your system. Increasing the number of threads
      may speed up optimization in some cases by running the likelihood 
      calculation in parallel. If \code{threads=1} (the default), OpenMP is disabled.} 
  \item{\dots}{Additional arguments to optim, such as lower and upper bounds}
}

\description{Fit N-mixture models with time-to-detection data.}

\value{unmarkedFitNmixTTD object describing model fit.}

\details{

This model extends time-to-detection (TTD) occupancy models to estimate site 
abundance using data from single or repeated visits. Latent abundance can be 
modeled as Poisson (\code{mixture="P"}) or negative binomial (\code{mixture="NB"}). 
Time-to-detection can be modeled as an exponential (\code{ttdDist="exp"}) or 
Weibull (\code{ttdDist="weibull"}) random variable with rate parameter \eqn{\lambda} 
and, for the Weibull, an additional shape parameter \eqn{k}. Note that 
\code{occuTTD} puts covariates on \eqn{\lambda} and not \eqn{1/\lambda}, i.e., 
the expected time between events. 

Assuming that there are \eqn{N} independent individuals at a site, and all
individuals have the same individual detection rate, the expected
detection rate across all individuals \eqn{\lambda} is equal to the the 
individual-level detection rate \eqn{r} multipled by the number of individuals
present \eqn{N}.

In the case where there are no detections before the maximum sample time at
a site (\code{surveyLength}) is reached, we are not sure if the site has 
\eqn{N=0} or if we just didn't wait long enough for a detection. We therefore 
must censor (\eqn{C} the exponential or Weibull distribution at the maximum survey 
length, \eqn{Tmax}. Thus, assuming true abundance at site \eqn{i} is 
\eqn{N_i}, and an exponential distribution for the TTD \eqn{y_i} (parameterized 
with the rate), then:

\deqn{y_i \sim Exponential(r_i * N_i) C(Tmax)}

Note that when \eqn{N_i = 0}, the exponential rate \eqn{lambda = 0} and the
scale is therefore \eqn{1 / 0 = Inf}, and thus the value will be censored
at \eqn{Tmax}.

Because in \code{unmarked} values of \code{NA} are typically used to indicate 
missing values that were a result of the sampling structure (e.g., lost data), 
we indicate a censored \eqn{y_i} in \code{nmixTTD} instead by setting 
\eqn{y_i = Tmax_i} in the \code{y} matrix provided to 
\code{\link{unmarkedFrameOccuTTD}}. You can provide either a single value of 
\eqn{Tmax} to the \code{surveyLength} argument of \code{unmarkedFrameOccuTTD}, 
or provide a matrix, potentially with a unique value of \eqn{Tmax} for each 
value of \code{y}. Note that in the latter case the value of \code{y} that will 
be interpreted by \code{nmixTTD} as a censored observation (i.e., \eqn{Tmax}) 
will differ between observations!

}

\references{
Strebel, N., Fiss, C., Kellner, K. F., Larkin, J. L., Kery, M., & Cohen, J (2021). 
  Estimating abundance based on time-to-detection data. Methods in Ecology and Evolution 12: 909-920.
}

\author{Ken Kellner \email{contact@kenkellner.com}}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuTTD}}}

\keyword{models}

\examples{

\dontrun{

# Simulate data
M = 1000 # Number of sites
nrep <- 3 # Number of visits per site
Tmax = 5 # Max duration of a visit
alpha1 = -1 # Covariate on rate
beta1 = 1 # Covariate on density
mu.lambda = 1 # Rate at alpha1 = 0
mu.dens = 1 # Density at beta1 = 0

covDet <- matrix(rnorm(M*nrep),nrow = M,ncol = nrep) #Detection covariate
covDens <- rnorm(M) #Abundance/density covariate
dens <- exp(log(mu.dens) + beta1 * covDens)
sum(N <- rpois(M, dens)) # Realized density per site
lambda <- exp(log(mu.lambda) + alpha1 * covDet) # per-individual detection rate
ttd <- NULL
for(i in 1:nrep) {
  ttd <- cbind(ttd,rexp(M, N*lambda[,i]))  # Simulate time to first detection per visit
}
ttd[N == 0,] <- 5 # Not observed where N = 0; ttd set to Tmax
ttd[ttd >= Tmax] <- 5 # Crop at Tmax

#Build unmarked frame
umf <- unmarkedFrameOccuTTD(y = ttd, surveyLength=5,
                            siteCovs = data.frame(covDens=covDens),
                            obsCovs = data.frame(covDet=as.vector(t(covDet))))

#Fit model
fit <- nmixTTD(~covDens, ~covDet, data=umf, K=max(N)+10)

#Compare to truth
cbind(coef(fit), c(log(mu.dens), beta1, log(mu.lambda), alpha1))

#Predict abundance/density values
head(predict(fit, type='state'))

}
}
