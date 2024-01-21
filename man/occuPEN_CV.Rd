\name{occuPEN_CV}

\alias{occuPEN_CV}

\title{Fit the MacKenzie et al. (2002) Occupancy Model with the penalized likelihood methods of Hutchinson et al. (2015) using cross-validation}

\usage{occuPEN_CV(formula, data, knownOcc=numeric(0), starts, method="BFGS",
    engine=c("C", "R"), lambdaVec=c(0,2^seq(-4,4)),
    pen.type = c("Bayes","Ridge"), k = 5, foldAssignments = NA,
    ...)}

\arguments{
    \item{formula}{Double right-hand side formula describing covariates of
        detection and occupancy in that order.}
    \item{data}{An \code{\link{unmarkedFrameOccu}} object}
    \item{knownOcc}{Vector of sites that are known to be occupied. These
    should be supplied as row numbers of the y matrix, eg, c(3,8) if
    sites 3 and 8 were known to be occupied a priori.}
    \item{starts}{Vector of parameter starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{engine}{Either "C" or "R" to use fast C++ code or native R
      code during the optimization.}
    \item{lambdaVec}{Vector of values to try for lambda.}
    \item{pen.type}{Which form of penalty to use.}
    \item{k}{Number of folds for k-fold cross-validation.}
    \item{foldAssignments}{Vector containing the number of the fold
      that each site falls into. Length of the vector should be equal
      to the number of sites, and the vector should contain k unique
      values. E.g. for 9 sites and 3 folds, c(1,2,3,1,2,3,1,2,3) or
      c(1,1,1,2,2,2,3,3,3).}
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{This function fits the occupancy model of MacKenzie et al
(2002) with the penalized methods of Hutchinson et al (2015) using
k-fold cross-validation to choose the penalty weight.}

\details{

See \code{\link{unmarkedFrameOccu}} for a description of how to supply data to 
the \code{data} argument.

This function wraps k-fold cross-validation around \code{occuPEN_CV}
for the "Bayes" and "Ridge" penalties of Hutchinson et al. (2015). The
user may specify the number of folds (\code{k}), the values to try
(\code{lambdaVec}), and the assignments of sites to folds
(\code{foldAssignments}). If \code{foldAssignments} is not provided,
the assignments are done pseudo-randomly, and the function attempts to
put some sites with and without positive detections in each fold. This
randomness introduces variability into the results of this function
across runs; to eliminate the randomness, supply foldAssignments. 
}

\value{unmarkedFitOccuPEN_CV object describing the model fit.}

\references{

Hutchinson, R. A., J. V. Valente, S. C. Emerson, M. G. Betts, and
T. G. Dietterich. 2015. Penalized Likelihood Methods Improve Parameter
Estimates in Occupancy Models. Methods in Ecology and Evolution. DOI:
10.1111/2041-210X.12368

MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege,
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy Rates
  When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.

}

\author{Rebecca A. Hutchinson}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccu}},
    \code{\link{occu}}, \code{\link{occuPEN}}, \code{\link{nonparboot}}}


\examples{

# Simulate occupancy data
set.seed(646)
nSites <- 60
nReps <- 2
covariates <- data.frame(veght=rnorm(nSites),
    habitat=factor(c(rep('A', 30), rep('B', 30))))

psipars <- c(-1, 1, -1)
ppars <- c(1, -1, 0)
X <- model.matrix(~veght+habitat, covariates) # design matrix
psi <- plogis(X \%*\% psipars)
p <- plogis(X \%*\% ppars)

y <- matrix(NA, nSites, nReps)
z <- rbinom(nSites, 1, psi)       # true occupancy state
for(i in 1:nSites) {
    y[i,] <- rbinom(nReps, 1, z[i]*p[i])
    }

# Organize data and look at it
umf <- unmarkedFrameOccu(y = y, siteCovs = covariates)
obsCovs(umf) <- covariates
head(umf)
summary(umf)

\dontrun{

# Fit some models
fmMLE <- occu(~veght+habitat ~veght+habitat, umf)
fmMLE@estimates

fm1penCV <- occuPEN_CV(~veght+habitat ~veght+habitat,
 umf,pen.type="Ridge", foldAssignments=rep(1:5,ceiling(nSites/5))[1:nSites])
fm1penCV@lambdaVec
fm1penCV@chosenLambda
fm1penCV@estimates

fm2penCV <- occuPEN_CV(~veght+habitat ~veght+habitat,
umf,pen.type="Bayes",foldAssignments=rep(1:5,ceiling(nSites/5))[1:nSites])
fm2penCV@lambdaVec
fm2penCV@chosenLambda
fm2penCV@estimates

# nonparametric bootstrap for uncertainty analysis:
# bootstrap is wrapped around the cross-validation
fm2penCV <- nonparboot(fm2penCV,B=10) # should use more samples
vcov(fm2penCV,method="nonparboot")

# Mean squared error of parameters:
mean((c(psipars,ppars)-c(fmMLE[1]@estimates,fmMLE[2]@estimates))^2)
mean((c(psipars,ppars)-c(fm1penCV[1]@estimates,fm1penCV[2]@estimates))^2)
mean((c(psipars,ppars)-c(fm2penCV[1]@estimates,fm2penCV[2]@estimates))^2)
}


}

\keyword{models}
