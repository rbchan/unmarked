\name{gdistsamp}
\alias{gdistsamp}
\title{
Fit the generalized distance sampling model of Chandler et al. (2011).
}
\description{
Extends the distance sampling model of Royle et al. (2004) to estimate
the probability of being available for detection. Also allows abundance
to be modeled using the negative binomial distribution.
}
\usage{
gdistsamp(lambdaformula, phiformula, pformula, data, keyfun =
c("halfnorm", "exp", "hazard", "uniform"), output = c("abund",
"density"), unitsOut = c("ha", "kmsq"), mixture = c("P", "NB"), K,
starts, method = "BFGS", se = TRUE, rel.tol=1e-4, ...)
}
\arguments{
  \item{lambdaformula}{
      A right-hand side formula describing the abundance covariates.
}
  \item{phiformula}{
      A right-hand side formula describing the availability covariates.
}
  \item{pformula}{
      A right-hand side formula describing the detection function covariates.
}
  \item{data}{
      An object of class \code{unmarkedFrameGDS}
}
  \item{keyfun}{
      One of the following detection functions: "halfnorm", "hazard", "exp",
          or "uniform." See details.
}
  \item{output}{
      Model either "density" or "abund"
}
  \item{unitsOut}{
      Units of density. Either "ha" or "kmsq" for hectares and square
          kilometers, respectively.
}
  \item{mixture}{
      Either "P" or "NB" for the Poisson and negative binomial models of
          abundance.
}
  \item{K}{
      An integer value specifying the upper bound used in the integration.
}
  \item{starts}{
      A numeric vector of starting values for the model parameters.
}
  \item{method}{
      Optimization method used by \code{\link{optim}}.
}
  \item{se}{
      logical specifying whether or not to compute standard errors.
}
\item{rel.tol}{relative accuracy for the integration of the detection function.
    See \link{integrate}. You might try adjusting this if you get an error
    message related to the integral. Alternatively, try providing
    different starting values.}
\item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}

}
\details{
  This model extends the model of Royle et al. (2004) by estimating the
  probability of being available for detection \eqn{\phi}{phi}. This
  effectively relaxes the assumption that \eqn{g(0)=1}. In other words,
  inividuals at a distance of 0 are not assumed to be detected with
  certainty. To estimate this additional parameter, replicate distance
  sampling data must be collected at each transect. Thus the data are
  collected at i = 1, 2, ..., R transects on t = 1, 2, ..., T
  occassions. As with the model of Royle et al. (2004), the detections
  must be binned into distance classes. These data must be formatted in
  a matrix with R rows, and JT columns where J is the number of distance
  classses. See \code{\link{unmarkedFrameGDS}} for more information.
}

\note{
  If you aren't interested in estimating phi, but you want to
  use the negative binomial distribution, simply set numPrimary=1 when
  formatting the data.
  }

\value{
    An object of class unmarkedFitGDS.
    }
\references{
    Royle, J. A., D. K. Dawson, and S. Bates. 2004. Modeling
    abundance effects in distance sampling. \emph{Ecology}
    85:1591-1597.

    Chandler, R. B, J. A. Royle, and D. I. King. 2011. Inference about
    density and temporary emigration in unmarked
    populations. \emph{Ecology}  92:1429--1435.
    }
\author{
    Richard Chandler \email{rbchan@uga.edu}
    }
\note{
    You cannot use obsCovs, but you can use yearlySiteCovs (a confusing name
    since this model isn't for multi-year data. It's just a hold-over
    from the colext methods of formatting data upon which it is based.)
    }
\seealso{
    \code{\link{distsamp}}
    }
\examples{


# Simulate some line-transect data

set.seed(36837)

R <- 50 # number of transects
T <- 5  # number of replicates
strip.width <- 50
transect.length <- 100
breaks <- seq(0, 50, by=10)

lambda <- 5 # Abundance
phi <- 0.6  # Availability
sigma <- 30 # Half-normal shape parameter

J <- length(breaks)-1
y <- array(0, c(R, J, T))
for(i in 1:R) {
    M <- rpois(1, lambda) # Individuals within the 1-ha strip
    for(t in 1:T) {
        # Distances from point
        d <- runif(M, 0, strip.width)
        # Detection process
        if(length(d)) {
            cp <- phi*exp(-d^2 / (2 * sigma^2)) # half-normal w/ g(0)<1
            d <- d[rbinom(length(d), 1, cp) == 1]
            y[i,,t] <- table(cut(d, breaks, include.lowest=TRUE))
            }
        }
    }
y <- matrix(y, nrow=R) # convert array to matrix

# Organize data
umf <- unmarkedFrameGDS(y = y, survey="line", unitsIn="m",
    dist.breaks=breaks, tlength=rep(transect.length, R), numPrimary=T)
summary(umf)


# Fit the model
m1 <- gdistsamp(~1, ~1, ~1, umf, output="density", K=50)

summary(m1)


backTransform(m1, type="lambda")
backTransform(m1, type="phi")
backTransform(m1, type="det")

\dontrun{
# Empirical Bayes estimates of abundance at each site
re <- ranef(m1)
plot(re, layout=c(10,5), xlim=c(-1, 20))
}

}
\keyword{ models }

