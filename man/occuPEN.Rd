\name{occuPEN}

\alias{occuPEN}

\title{Fit the MacKenzie et al. (2002) Occupancy Model with the penalized likelihood methods of Hutchinson et al. (2015)}

\usage{occuPEN(formula, data, knownOcc=numeric(0), starts, method="BFGS",
    engine=c("C", "R"), lambda=0, pen.type = c("Bayes","Ridge","MPLE"), se = FALSE, ...)}

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
    \item{lambda}{Penalty weight parameter.}
    \item{pen.type}{Which form of penalty to use.}
    \item{se}{Nonfunctional; always FALSE; for compatibility with
              other parts of \code{unmarked} }
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{This function fits the occupancy model of MacKenzie et al (2002) with the penalized methods of Hutchinson et al (2015).}

\details{

See \code{\link{unmarkedFrame}} and \code{\link{unmarkedFrameOccu}} for a
description of how to supply data to the \code{data} argument.

\code{occuPEN} fits the standard occupancy model based on
zero-inflated binomial models (MacKenzie et al. 2006, Royle and
Dorazio 2008) using the penalized likelihood methods described in
Hutchinson et al. (2015).  See \code{\link{occu}} for model
details. \code{occuPEN} returns parameter estimates that maximize a
penalized likelihood in which the penalty is specified by the
\code{pen.type} argument. The penalty function is weighted by
\code{lambda}. 

The MPLE method includes an equation for computing \code{lambda}
(Moreno & Lele, 2010). If the value supplied does not equal match the
one computed with this equation, the supplied value is used anyway
(with a warning).
}

\value{unmarkedFitOccuPEN object describing the model fit.}

\references{

Hutchinson, R. A., J. V. Valente, S. C. Emerson, M. G. Betts, and
T. G. Dietterich. 2015. Penalized Likelihood Methods Improve Parameter
Estimates in Occupancy Models. Methods in Ecology and Evolution. DOI:
10.1111/2041-210X.12368

MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege,
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy Rates
  When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.

MacKenzie, D. I. et al. 2006. \emph{Occupancy Estimation and Modeling}.
  Amsterdam: Academic Press.

Moreno, M. and S. R. Lele. 2010. Improved estimation of site occupancy
using penalized likelihood. Ecology 91: 341-346.

Royle, J. A. and R. Dorazio. 2008. \emph{Hierarchical Modeling and Inference
  in Ecology}. Academic Press.

}

\author{Rebecca A. Hutchinson}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccu}},
    \code{\link{occu}}, \code{\link{computeMPLElambda}}, 
    \code{\link{occuPEN_CV}}, \code{\link{nonparboot}}}


\examples{

# Simulate occupancy data
set.seed(344)
nSites <- 100
nReps <- 2
covariates <- data.frame(veght=rnorm(nSites),
    habitat=factor(c(rep('A', nSites/2), rep('B', nSites/2))))

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


# Fit some models
fmMLE <- occu(~veght+habitat ~veght+habitat, umf)
fm1pen <- occuPEN(~veght+habitat ~veght+habitat, umf,lambda=0.33,pen.type="Ridge")
fm2pen <- occuPEN(~veght+habitat ~veght+habitat, umf,lambda=1,pen.type="Bayes")

# MPLE:
fm3pen <- occuPEN(~veght+habitat ~veght+habitat, umf,lambda=0.5,pen.type="MPLE")
MPLElambda = computeMPLElambda(~veght+habitat ~veght+habitat, umf) 
fm4pen <- occuPEN(~veght+habitat ~veght+habitat, umf,lambda=MPLElambda,pen.type="MPLE")

# nonparametric bootstrap for uncertainty analysis:
fm1pen <- nonparboot(fm1pen,B=20) # should use more samples
vcov(fm1pen,method="nonparboot")




}

\keyword{models}
