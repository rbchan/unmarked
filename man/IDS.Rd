\name{IDS}
\alias{IDS}
\alias{hist,unmarkedFitIDS-method}
\alias{names,unmarkedFitIDS-method}

\title{
Fit the integrated distance sampling model of Kery et al. (2022).
}

\description{Model abundance using a combination of distance sampling data (DS)
  and other similar data types, including simple point counts (PC) and 
  occupancy/detection-nondetection (OC/DND) data.}

\usage{
IDS(lambdaformula = ~1, detformulaDS = ~1, detformulaPC = NULL, detformulaOC = NULL,
    dataDS, dataPC = NULL, dataOC = NULL, availformula = NULL,
    durationDS = NULL, durationPC = NULL, durationOC = NULL, keyfun = "halfnorm",
    maxDistPC, maxDistOC, K = 100, unitsOut = "ha", 
    starts = NULL, method = "BFGS", ...)
}

\arguments{
  \item{lambdaformula}{Formula for abundance}
  \item{detformulaDS}{Formula for distance-based (DS) detection probability}
  \item{detformulaPC}{Formula for point count detection probability. If
    \code{NULL}, will share a model with DS detection probability}
  \item{detformulaOC}{Formula for occupancy/detection-nondetection detection 
    probability. If \code{NULL}, will share a model with DS detection probability}
  \item{dataDS}{An object of class \code{unmarkedFrameDS}. Required}
  \item{dataPC}{An object of class \code{unmarkedFramePCount}. If \code{NULL},
    no PC data will be used in the model}
  \item{dataOC}{An object of class \code{unmarkedFrameOccu}. If \code{NULL},
    no OC/DND data will be used in the model}
  \item{availformula}{Optional. If specified, formula for availability. Only possible to
    use if you have variable detection survey lengths (see below)}
  \item{durationDS}{Optional. Vector of survey durations at each distance sampling site} 
  \item{durationPC}{Optional. Vector of survey durations at each PC site} 
  \item{durationOC}{Optional. Vector of survey durations at each OC/DND site} 
  \item{keyfun}{Distance sampling key function; either "halfnorm" or "exp"}
  \item{maxDistPC}{Maximum observation distance for PC surveys; defaults to
    maximum of distance bins from the distance sampling data}
  \item{maxDistOC}{Maximum observation distance for OC/DND surveys; defaults to
    maximum of distance bins from the distance sampling data}
  \item{K}{Integer, upper bound for integrating out latent abundance. Only used if
    you have included OC/DND data}
  \item{unitsOut}{Units of density for output. Either "ha" or "kmsq" for 
    hectares and square kilometers, respectively}
  \item{starts}{A numeric vector of starting values for the model parameters}
  \item{method}{Optimization method used by \code{\link{optim}}}
  \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
}

\value{An object of class unmarkedFitIDS}

\details{
This function facilitates a combined analysis of distance sampling data (DS) with other similar
data types, including simple point counts (PC) and occupancy/detection-nondetection (OC/DND) data.
The combined approach capitalizes on the strengths and minimizes the weaknesses
of each type. The PC and OC/DND data are viewed as latent distance sampling surveys
with an underlying abundance model shared by all data types. All analyses
must include some distance sampling data, but can include either PC or DND data
or both. If surveys are of variable duration, it is also possible to estimate
availability.

Input data must be provided as a series of separate \code{unmarkedFrame}s: 
\code{unmarkedFrameDS} for the distance sampling data, \code{unmarkedFramePCount}
for the point count data, and \code{unmarkedFrameOccu} for OC/DND data.
See the help files for these objects for guidance on how to organize the data.
}

\note{
  Simulations indicated estimates of availability were very unreliable when
  including detection/non-detection data, so the function will not allow you
  to use DND data and estimate availability at the same time. 
  In general estimation of availability can be difficult; use simulations
  to see how well it works for your specific situation.
}

\references{
  Kery M, Royle JA, Hallman T, Robinson WD, Strebel N, Kellner KF. 2024. 
  Integrated distance sampling models for simple point counts. Ecology. 
}
\author{Ken Kellner \email{contact@kenkellner.com}}
\seealso{\code{\link{distsamp}}}

\examples{

\dontrun{

# Simulate data based on a real dataset

# Formulas for each model
formulas <- list(lam=~elev, ds=~1, phi=~1)

# Sample sizes
design <- list(Mds=2912, J=6, Mpc=506)

# Model parameters
coefs <- list(lam = c(intercept=3, elev=-0.5),
              ds = c(intercept=-2.5),
              phi = c(intercept=-1.3))

# Survey durations
durs <- list(ds = rep(5, design$Mds), pc=runif(design$Mpc, 3, 30))

set.seed(456)
sim_umf <- simulate("IDS", # name of model we are simulating for
                    nsim=1, # number of replicates
                    formulas=formulas, 
                    coefs=coefs,
                    design=design,
                    # arguments used by unmarkedFrameDS
                    dist.breaks = seq(0, 0.30, length.out=7),
                    unitsIn="km", 
                    # arguments used by IDS
                    # could also have e.g. keyfun here
                    durationDS=durs$ds, durationPC=durs$pc, durationOC=durs$oc,
                    maxDistPC=0.5, maxDistOC=0.5,
                    unitsOut="kmsq")

# Look at the results
lapply(sim_umf, head)

# Fit a model
(mod_sim <- IDS(lambdaformula = ~elev, detformulaDS = ~1,
                dataDS=sim_umf$ds, dataPC=sim_umf$pc,
                availformula = ~1, durationDS=durs$ds, durationPC=durs$pc,
                maxDistPC=0.5,
                unitsOut="kmsq"))

# Compare with known parameter values
# Note:  this is an unusually good estimate of availability
# It is hard to estimate in most cases
cbind(truth=unlist(coefs), est=coef(mod_sim))

# Predict density at each distance sampling site
head(predict(mod_sim, 'lam'))

}

}
