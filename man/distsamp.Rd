\name{distsamp}
\alias{distsamp}
\title{Fit the hierarchical distance sampling model of Royle et al. (2004)}
\usage{distsamp(formula, data, keyfun=c("halfnorm", "exp",
  "hazard", "uniform"), output=c("density", "abund"),
  unitsOut=c("ha", "kmsq"), starts, method="BFGS", se=TRUE,
  engine=c("C", "R"), rel.tol=0.001, ...)
}

\description{
  Fit the hierarchical distance sampling model of Royle et al. (2004)
  to line or point transect data recorded in discrete distance intervals.
}

\arguments{
  \item{formula}{Double right-hand formula describing detection
    covariates followed by abundance covariates. ~1 ~1 would be a null
    model.}
  \item{data}{object of class \code{unmarkedFrameDS}, containing response
    matrix, covariates, distance interval cut points, survey type ("line"
    or "point"), transect lengths (for survey = "line"), and units ("m"
    or "km") for cut points and transect lengths. See example for set up.}
  \item{keyfun}{One of the following detection functions:
    "halfnorm", "hazard", "exp", or "uniform." See details.}
  \item{output}{Model either "density" or "abund"}
  \item{unitsOut}{Units of density. Either "ha" or "kmsq" for hectares and
    square kilometers, respectively.}
  \item{starts}{Vector of starting values for parameters.}
  \item{method}{Optimization method used by \code{\link{optim}}.}
  \item{se}{logical specifying whether or not to compute standard errors.}
  \item{engine}{Use code written in C++ or R}
  \item{rel.tol}{Requested relative accuracy of the integral, see
    \code{\link{integrate}}}
  \item{\dots}{Additional arguments to optim, such as lower and upper
    bounds}
}

\details{
Unlike conventional distance sampling, which uses the 'conditional on
detection' likelihood formulation, this model is based upon the
unconditional likelihood and allows for modeling both abundance and
detection  function parameters.

The latent transect-level abundance distribution
\eqn{f(N | \mathbf{\theta})}{f(N | theta)} assumed to be
Poisson with mean \eqn{\lambda}{lambda} (but see \code{\link{gdistsamp}}
for alternatives).

The detection process is modeled as multinomial:
\eqn{y_{ij} \sim Multinomial(N_i, \pi_{ij})}{y_ij ~ Multinomial(N_i, pi_i1, pi_i2, ..., pi_iJ)},
where \eqn{\pi_{ij}}{pi_ij} is the multinomial cell probability for transect i in
distance class j. These are computed based upon a detection function
\eqn{g(x | \mathbf{\sigma})}{g(x | sigma)}, such as the half-normal,
negative exponential, or hazard rate.

Parameters \eqn{\lambda}{lambda} and \eqn{\sigma}{sigma} can be vectors
affected by transect-specific covariates using the log link.}

\value{unmarkedFitDS object (child class of \code{\link{unmarkedFit-class}})
describing the model fit.}

\note{
You cannot use obsCovs.
}

\author{Richard Chandler \email{rbchan@uga.edu}}

\references{

Royle, J. A., D. K. Dawson, and S. Bates (2004) Modeling
  abundance effects in distance sampling. \emph{Ecology} 85, pp. 1591-1597.

Sillett, S. and Chandler, R.B. and Royle, J.A. and Kery, M. and
  Morrison, S.A. In Press. Hierarchical distance sampling models to
  estimate population size and habitat-specific abundance of an island
  endemic. \emph{Ecological Applications}

}

\seealso{\code{\link{unmarkedFrameDS}},
    \code{\link{unmarkedFit-class}} \code{\link{fitList}},
    \code{\link{formatDistData}}, \code{\link{parboot}},
    \code{\link{sight2perpdist}}, \code{\link{detFuns}},
    \code{\link{gdistsamp}}, \code{\link{ranef}}.
    Also look at vignette("distsamp").
  }


\examples{
## Line transect examples

data(linetran)

ltUMF <- with(linetran, {
   unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4),
   siteCovs = data.frame(Length, area, habitat),
   dist.breaks = c(0, 5, 10, 15, 20),
   tlength = linetran$Length * 1000, survey = "line", unitsIn = "m")
   })

ltUMF
summary(ltUMF)
hist(ltUMF)

# Half-normal detection function. Density output (log scale). No covariates.
(fm1 <- distsamp(~ 1 ~ 1, ltUMF))

# Some methods to use on fitted model
summary(fm1)
backTransform(fm1, type="state")                # animals / ha
exp(coef(fm1, type="state", altNames=TRUE))     # same
backTransform(fm1, type="det")                  # half-normal SD
hist(fm1, xlab="Distance (m)")	# Only works when there are no det covars
# Empirical Bayes estimates of posterior distribution for N_i
plot(ranef(fm1, K=50))

# Effective strip half-width
(eshw <- integrate(gxhn, 0, 20, sigma=10.9)$value)

# Detection probability
eshw / 20 # 20 is strip-width


# Halfnormal. Covariates affecting both density and and detection.
(fm2 <- distsamp(~area + habitat ~ habitat, ltUMF))

# Hazard-rate detection function.
(fm3 <- distsamp(~ 1 ~ 1, ltUMF, keyfun="hazard"))

# Plot detection function.
fmhz.shape <- exp(coef(fm3, type="det"))
fmhz.scale <- exp(coef(fm3, type="scale"))
plot(function(x) gxhaz(x, shape=fmhz.shape, scale=fmhz.scale), 0, 25,
	xlab="Distance (m)", ylab="Detection probability")



## Point transect examples

# Analysis of the Island Scrub-jay data.
# See Sillett et al. (In press)

data(issj)
str(issj)

jayumf <- unmarkedFrameDS(y=as.matrix(issj[,1:3]),
 siteCovs=data.frame(scale(issj[,c("elevation","forest","chaparral")])),
 dist.breaks=c(0,100,200,300), unitsIn="m", survey="point")

(fm1jay <- distsamp(~chaparral ~chaparral, jayumf))




\dontrun{

data(pointtran)

ptUMF <- with(pointtran, {
	unmarkedFrameDS(y = cbind(dc1, dc2, dc3, dc4, dc5),
	siteCovs = data.frame(area, habitat),
	dist.breaks = seq(0, 25, by=5), survey = "point", unitsIn = "m")
	})

# Half-normal.
(fmp1 <- distsamp(~ 1 ~ 1, ptUMF))
hist(fmp1, ylim=c(0, 0.07), xlab="Distance (m)")

# effective radius
sig <- exp(coef(fmp1, type="det"))
ea <- 2*pi * integrate(grhn, 0, 25, sigma=sig)$value # effective area
sqrt(ea / pi) # effective radius

# detection probability
ea / (pi*25^2)

}
}

\keyword{models}
