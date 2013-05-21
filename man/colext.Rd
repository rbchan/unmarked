\name{colext}
\alias{colext}
\title{Fit the dynamic occupancy model of MacKenzie et. al (2003)}
\usage{colext(psiformula= ~1, gammaformula =  ~ 1, epsilonformula = ~ 1,
    pformula = ~ 1, data, starts, method="BFGS", se=TRUE, ...)}

\description{Estimate parameters of the colonization-extinction model,
including covariate-dependent rates and detection process.}

\details{This function fits the colonization-extinction model of
MacKenzie et al (2003).  The colonization and extinction rates can be
modeled with covariates that vary yearly at each site using a logit
link.  These covariates are supplied by special unmarkedMultFrame
\code{yearlySiteCovs} slot.  These parameters are specified using the
\code{gammaformula} and \code{epsilonformula} arguments.  The initial
probability of occupancy is modeled by covariates specified in the
\code{psiformula}.

The conditional detection rate can also be modeled as a function of
covariates that vary at the secondary sampling period (ie., repeat
visits).  These covariates are specified by the first part of the
\code{formula} argument and the data is supplied via the usual
\code{obsCovs} slot.

The projected and smoothed trajectories (Weir et al 2009) can be
obtained from the \code{smoothed.mean} and \code{projected.mean} slots
(see examples).
}
\seealso{\code{\link{nonparboot}}, \code{\link{unmarkedMultFrame}}, and \code{\link{formatMult}}}
\value{unmarkedFitColExt object describing model fit.}
\keyword{models}
\references{

MacKenzie, D.I. et al. (2002) Estimating Site Occupancy Rates When Detection Probabilities Are Less Than One. Ecology, 83(8), 2248-2255.

MacKenzie, D. I., J. D. Nichols, J. E. Hines, M. G. Knutson, and A. B. Franklin.
  2003. Estimating site occupancy, colonization, and local extinction when a
  species is detected imperfectly. \emph{Ecology} 84:2200--2207.

MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.Amsterdam: Academic Press.

Weir L. A., Fiske I. J., Royle J. (2009) Trends in Anuran Occupancy from
Northeastern States of the North American Amphibian Monitoring
Program. Herpetological Conservation and Biology. 4(3):389-402.

}
\arguments{\item{psiformula}{Right-hand sided formula for the initial
    probability of occupancy at each site.}
  \item{gammaformula}{Right-hand sided formula for colonization probability.}
  \item{epsilonformula}{Right-hand sided formula for extinction probability.}
  \item{pformula}{Right-hand sided formula for detection probability.}
\item{data}{unmarkedMultFrame object that supplies the data (see \code{\link{unmarkedMultFrame}}).}
\item{starts}{optionally, initial values for parameters in the optimization.}
\item{method}{Optimization method used by \code{\link{optim}}.}
\item{se}{logical specifying whether or not to compute standard errors.}
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
    }
\examples{

# Fake data
R <- 4 # number of sites
J <- 3 # number of secondary sampling occasions
T <- 2 # number of primary periods

y <- matrix(c(
   1,1,0,  0,0,0,
   0,0,0,  0,0,0,
   1,1,1,  1,1,0,
   1,0,1,  0,0,1), nrow=R, ncol=J*T, byrow=TRUE)
y

site.covs <- data.frame(x1=1:4, x2=factor(c('A','B','A','B')))
site.covs

yearly.site.covs <- list(
   year = matrix(c(
      'year1', 'year2',
      'year1', 'year2',
      'year1', 'year2',
      'year1', 'year2'), nrow=R, ncol=T, byrow=TRUE)
      )
yearly.site.covs

obs.covs <- list(
   x4 = matrix(c(
      -1,0,1,  -1,1,1,
      -2,0,0,  0,0,2,
      -3,1,0,  1,1,2,
      0,0,0,   0,1,-1), nrow=R, ncol=J*T, byrow=TRUE),
   x5 = matrix(c(
      'a','b','c',  'a','b','c',
      'd','b','a',  'd','b','a',
      'a','a','c',  'd','b','a',
      'a','b','a',  'd','b','a'), nrow=R, ncol=J*T, byrow=TRUE))
obs.covs

umf <- unmarkedMultFrame(y=y, siteCovs=site.covs,
    yearlySiteCovs=yearly.site.covs, obsCovs=obs.covs,
    numPrimary=2)                  # organize data
umf                                # look at data
summary(umf)                       # summarize
fm <- colext(~1, ~1, ~1, ~1, umf)  # fit a model
fm



\dontrun{
# Real data
data(frogs)
umf <- formatMult(masspcru)
obsCovs(umf) <- scale(obsCovs(umf))

## Use 1/4 of data just for run speed in example
umf <- umf[which((1:numSites(umf)) \%\% 4 == 0),]

## constant transition rates
(fm <- colext(psiformula = ~ 1,
gammaformula = ~ 1,
epsilonformula = ~ 1,
pformula = ~ JulianDate + I(JulianDate^2), umf, control = list(trace=1, maxit=1e4)))

## get the trajectory estimates
smoothed(fm)
projected(fm)

# Empirical Bayes estimates of number of sites occupied in each year
re <- ranef(fm)
modes <- colSums(bup(re, stat="mode"))
plot(1:7, modes, xlab="Year", ylab="Sites occupied", ylim=c(0, 70))

## Find bootstrap standard errors for smoothed trajectory
fm <- nonparboot(fm, B = 100)  # This takes a while!
fm@smoothed.mean.bsse

## try yearly transition rates
yearlySiteCovs(umf) <- data.frame(year = factor(rep(1:7, numSites(umf))))
(fm.yearly <- colext(psiformula = ~ 1,
gammaformula = ~ year,
epsilonformula = ~ year,
pformula = ~ JulianDate + I(JulianDate^2), umf,
	control = list(trace=1, maxit=1e4)))
}

}
