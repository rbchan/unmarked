\name{simulate-methods}
\docType{methods}
\alias{simulate-methods}
\alias{simulate,unmarkedFit-method}
\alias{simulate,unmarkedFrame-method}
\alias{simulate,unmarkedFitColExt-method}
\alias{simulate,unmarkedFitDS-method}
\alias{simulate,unmarkedFitMPois-method}
\alias{simulate,unmarkedFitOccu-method}
\alias{simulate,unmarkedFitOccuRN-method}
\alias{simulate,unmarkedFitOccuFP-method}
\alias{simulate,unmarkedFitOccuMulti-method}
\alias{simulate,unmarkedFitOccuMS-method}
\alias{simulate,unmarkedFitOccuTTD-method}
\alias{simulate,unmarkedFitNmixTTD-method}
\alias{simulate,unmarkedFitPCount-method}
\alias{simulate,unmarkedFitPCO-method}
\alias{simulate,unmarkedFitGMM-method}
\alias{simulate,unmarkedFitGDS-method}
\alias{simulate,unmarkedFitGPC-method}
\alias{simulate,unmarkedFitGDR-method}
\alias{simulate,unmarkedFitIDS-method}
\alias{simulate,unmarkedFitDailMadsen-method}
\alias{simulate,unmarkedFitGOccu-method}
\alias{simulate,unmarkedFitOccuCOP-method}
\alias{simulate,unmarkedFitOccuRNMulti-method}
\alias{simulate,character-method}

\title{Methods for Function simulate in Package `unmarked'}
\description{
Simulate data from a fitted model or an \code{unmarkedFrame}.
}
\usage{
\S4method{simulate}{unmarkedFit}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFrame}(object, nsim=1, seed=NULL, model = NULL,
                                   coefs = NULL, quiet = FALSE, ...)
}

\arguments{
\item{object}{Fitted model or \code{unmarkedFrame}.}
\item{nsim}{Number of simulations}
\item{seed}{Seed for random number generator. Not currently implemented}
\item{na.rm}{Logical, should missing values be removed?}
\item{model}{The model to use when \code{object} is an \code{unmarkedFrame}
  used for multiple model types. For example, if the \code{object} is an
  \code{unmarkedFrameOccu}, model should be set to \code{occu} or \code{occuRN}.}
\item{coefs}{List with one element per submodel. Each list element should be
  named with the corresponding submodel, and should be a numeric vector of 
  parameter values to use for that submodel when simulating the dataset.
  Note that parameter values should be on the inverse link scale.
  The number of parameter values in the vector depends on the model specified,
  covariates, etc. If you are not sure how to specify this list, set
  \code{coefs = NULL} and the function will return the correct structure.}
\item{quiet}{If TRUE, don't print informational messages.}
\item{...}{Used only for the \code{unmarkedFrame} method. Arguments to send to 
  the corresponding fitting function. Most importantly this will include
  formula arguments, but could also include distributions, key functions, etc. 
  For example, for simulating occupancy data, you must also supply the argument
  \code{formula = ~1~1} for a no-covariate model, \code{formula=~1~x} for a
  covariate effect of \code{x} on occupancy, etc. See examples below.}
}

\examples{

\dontrun{

# Simulation of an occupancy dataset from scratch

# First create an unmarkedFrame with the correct design

M <- 300 # number of sites
J <- 5   # number of occasions

# The values in the y-matrix don't matter as they will be simulated
# We can supply them as all NAs
y <- matrix(NA, M, J)

# Site covariate
x <- rnorm(M)

# Create unmarkedFrame
umf <- unmarkedFrameOccu(y = y, siteCovs = data.frame(x = x))

# Must specify model = occu since unmarkedFrameOccu is also used for occuRN
# the formula species the specific model structure we want to simulate
# If we don't specify coefs, unmarked will generate a template you can copy and use
simulate(umf, model = occu, formula = ~1~x)

# Now set coefs
# Here we imply a mean occupancy and mean detection of 0.5
# (corresponding to values of 0 on the inverse link scale) and a positive effect of x
s <- simulate(umf, model = occu, formula = ~1~x, 
              coefs = list(state = c(0,0.3), det = 0))

head(s[[1]])

occu(~1~x, s[[1]])

# For some models we can also include a random effect
# add a factor covariate
umf@siteCovs$x2 <- factor(sample(letters[1:10], M, replace=TRUE))

# The final value in coefs now represents the random effect SD for x2
s <- simulate(umf, model = occu, formula = ~1~x+(1|x2), 
              coefs = list(state = c(0,0.3, 1), det = 0))

head(s[[1]])

occu(~1~x+(1|x2), s[[1]])

# Here's a more complicated example simulating a gdistsamp dataset
# using a negative binomial distribution
M <- 100
J <- 3
T <- 2
y <- matrix(NA, M, J*T)
umf2 <- unmarkedFrameGDS(y=y, 
                         siteCovs=data.frame(x=rnorm(M)),
                         dist.breaks = c(0, 10, 20, 30), unitsIn='m',
                         numPrimary = T, survey="point")

cf <- list(lambda=c(1, 0.3), phi=0, det=c(log(20), 0), alpha=log(1))

# Note we now also supply another argument mixture="NB" to ... 
s2 <- simulate(umf2, coefs=cf, lambdaformula=~x, phiformula=~1, pformula=~x,
               mixture="NB")
head(s2[[1]])

gdistsamp(~x, ~1, ~x, s2[[1]], mixture="NB")

}
}
