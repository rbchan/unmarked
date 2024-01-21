\name{occu}

\alias{occu}

\title{Fit the MacKenzie et al. (2002) Occupancy Model}

\usage{occu(formula, data, knownOcc=numeric(0), linkPsi=c("logit", "cloglog"),
            starts, method="BFGS", se=TRUE, engine=c("C", "R", "TMB"),
            threads = 1, ...)}

\arguments{
    \item{formula}{Double right-hand side formula describing covariates of
        detection and occupancy in that order.}
    \item{data}{An \code{\link{unmarkedFrameOccu}} object}
    \item{knownOcc}{Vector of sites that are known to be occupied. These
    should be supplied as row numbers of the y matrix, eg, c(3,8) if
    sites 3 and 8 were known to be occupied a priori.}
    \item{linkPsi}{Link function for the occupancy model. Options are  
      \code{"logit"} for the standard occupancy model or \code{"cloglog"} 
      for the complimentary log-log link, which relates occupancy
      to site-level abundance. See details.}
    \item{starts}{Vector of parameter starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{Logical specifying whether or not to compute standard
      errors.}
    \item{engine}{Code to use for optimization. Either "C" for fast C++ code,
      "R" for native R code, or "TMB" for Template Model Builder. "TMB" is
      used automatically if your formula contains random effects.}
    \item{threads}{Set the number of threads to use for optimization in C++, if
      OpenMP is available on your system. Increasing the number of threads
      may speed up optimization in some cases by running the likelihood 
      calculation in parallel. If \code{threads=1} (the default), OpenMP is disabled.} 
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{This function fits the single season occupancy model of MacKenzie et al (2002).}

\details{

See \code{\link{unmarkedFrameOccu}} for a description of how to supply data to the 
\code{data} argument.

\code{occu} fits the standard occupancy model based on zero-inflated
binomial models (MacKenzie et al. 2006, Royle and Dorazio
2008).  The occupancy state process (\eqn{z_i}) of site \eqn{i} is
modeled as

\deqn{z_i \sim Bernoulli(\psi_i)}{z_i ~ Bernoulli(psi_i)}

The observation process is modeled as

\deqn{y_{ij}|z_i \sim Bernoulli(z_i p_{ij})}{y_ij | z_i ~ Bernoulli(z_i  * p_ij)}

By default, covariates of \eqn{\psi_i}{psi_i} and \eqn{p_{ij}}{p_ij} are modeled
using the logit link according to the \code{formula} argument.  The formula is a double right-hand sided formula
like \code{~ detform ~ occform} where \code{detform} is a formula for the detection process and \code{occform} is a
formula for the partially observed occupancy state.  See \link{formula} for details on constructing model formulae
in \R.

When \code{linkPsi = "cloglog"}, the complimentary log-log link 
function is used for \eqn{psi} instead of the logit link. The cloglog link
relates occupancy probability to the intensity parameter of an underlying
Poisson process (Kery and Royle 2016). Thus, if abundance at a site is 
can be modeled as \eqn{N_i ~ Poisson(\lambda_i)}, where 
\eqn{log(\lambda_i) = \alpha + \beta*x}, then presence/absence data at the 
site can be modeled as \eqn{Z_i ~ Binomial(\psi_i)} where 
\eqn{cloglog(\psi_i) = \alpha + \beta*x}. 

}

\value{unmarkedFitOccu object describing the model fit.}

\references{

Kery, Marc, and J. Andrew Royle. 2016. \emph{Applied Hierarchical Modeling in
  Ecology}, Volume 1. Academic Press. 

MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege,
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy Rates
  When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.

MacKenzie, D. I. et al. 2006. \emph{Occupancy Estimation and Modeling}.
  Amsterdam: Academic Press.

Royle, J. A. and R. Dorazio. 2008. \emph{Hierarchical Modeling and Inference
  in Ecology}. Academic Press.

}

\author{Ian Fiske}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccu}},
    \code{\link{modSel}}, \code{\link{parboot}}}


\examples{

data(frogs)
pferUMF <- unmarkedFrameOccu(pfer.bin)
plot(pferUMF, panels=4)
# add some fake covariates for illustration
siteCovs(pferUMF) <- data.frame(sitevar1 = rnorm(numSites(pferUMF)))

# observation covariates are in site-major, observation-minor order
obsCovs(pferUMF) <- data.frame(obsvar1 = rnorm(numSites(pferUMF) * obsNum(pferUMF)))

(fm <- occu(~ obsvar1 ~ 1, pferUMF))

confint(fm, type='det', method = 'normal')
confint(fm, type='det', method = 'profile')

# estimate detection effect at obsvars=0.5
(lc <- linearComb(fm['det'],c(1,0.5)))

# transform this to probability (0 to 1) scale and get confidence limits
(btlc <- backTransform(lc))
confint(btlc, level = 0.9)

# Empirical Bayes estimates of proportion of sites occupied
re <- ranef(fm)
sum(bup(re, stat="mode"))

}

\keyword{models}
