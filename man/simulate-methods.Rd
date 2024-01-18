\name{simulate-methods}
\docType{methods}
\alias{simulate-methods}
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
\alias{simulate,character-method}

\title{Methods for Function simulate in Package `unmarked'}
\description{
Simulate data from a fitted model.
}
\usage{
\S4method{simulate}{unmarkedFitColExt}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitDS}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitMPois}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitOccu}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitOccuRN}(object, nsim, seed, na.rm)
\S4method{simulate}{unmarkedFitPCount}(object, nsim, seed, na.rm)
\S4method{simulate}{character}(object, nsim=1, seed=NULL, formulas, coefs=NULL,
  design, guide=NULL, ...)
}

\arguments{
\item{object}{Fitted model of appropriate S4 class}
\item{nsim}{Number of simulations}
\item{seed}{Seed for random number generator. Not currently implemented}
\item{na.rm}{Logical, should missing values be removed?}
\item{formulas}{
  A named list of formulas, one per submodel (e.g. a formula for occupancy
  \code{"state"} and a formula for detection \code{"det"}). To get the correct
  submodel names for a given model, fit an example for that model, and then
  call \code{names(fitted_model)}
}
\item{coefs}{
  A named list of vectors of coefficients associated with the regression
  intercepts and slopes for each submodel. List should be named as with
  \code{formulas} above. Each element of the list should be a named vector,
  where the names correspond to the names of the parameters in the model
  (intercept and covariates). If you are not sure how to structure this list,
  just run \code{simulate} with \code{coefs=NULL}; this will generate
  a template list you can copy and fill in.
}
\item{design}{
  A named list of components of the study design. Must include at least \code{M},
  the number of sites, and \code{J} the number of observations per site. If you
  are fitting a model with multiple primary periods you must also provide
  \code{T}, the number of primary periods.
}
\item{guide}{
  An optional list defining the format (continuous or categorical/factor) and distribution,
  if continuous, of covariates you want to simulate. By default all covariates
  are simulated from a standard normal. See example below for an example of
  how to specify entries in the \code{guide} list.
}
\item{...}{
  Additional arguments that are needed to fully specify the simulated dataset
  for a particular model. For example, \code{mixture} for \code{pcount} models
  or \code{keyfun} for \code{distsamp} models.
}
}

\section{Methods}{
\describe{
\item{object = "unmarkedFitColExt"}{A model fit by \code{\link{colext}}}
\item{object = "unmarkedFitDS"}{A model fit by \code{\link{distsamp}}}
\item{object = "unmarkedFitMPois"}{A model fit by \code{\link{multinomPois}}}
\item{object = "unmarkedFitOccu"}{A model fit by \code{\link{occu}}}
\item{object = "unmarkedFitOccuRN"}{A model fit by \code{\link{occuRN}}}
\item{object = "unmarkedFitPCount"}{A model fit by \code{\link{pcount}}}
\item{object = "character"}{An \code{unmarkedFrame} of the appropriate type}
}}
\keyword{methods}

\examples{

\dontrun{

# Simulation of an occupancy dataset from scratch

# Formulas for each submodel
# occupancy is a function of elevation, detection is intercept-only
forms <- list(state=~elev, det=~1)

# Specify list of coefficients - there must be a value for each
# covariate plus an intercept for each submodel
coefs <- list(state=c(intercept=0, elev=-0.4), det=c(intercept=0))

# Study design
design <- list(M=300, J=8) # 300 sites, 8 occasions per site

# If we don't specify coefs, unmarked will generate a template you can copy and use
simulate("occu", formulas=forms, design=design)

# Generate unmarkedFrameOccu
occu_umf <- simulate("occu", formulas=forms, coefs=coefs, design=design)
head(occu_umf) # note one covariate, elev

# What if we wanted to add a categorical/factor covariate or
# customize the distribution of elev?
# Use the guide argument

# Updated formulas with new covariate
forms2 <- list(state=~elev+landcover, det=~1)

# Guide
# landcover is factor, you must provide the levels
guide <- list(landcover=factor(levels=c("forest","grass")),  
              elev=list(dist=rnorm, mean=2, sd=0.5)) # custom distribution

# Updated coefficients list
coefs2 <- list(state=c(intercept=0, elev=-0.4, landcovergrass=0.2), det=c(intercept=0))

# Simulate new dataset
head(simulate("occu", formulas=forms2, coefs=coefs2, design=design, guide=guide))
# Note new categorical covariate

# For some models you may want to specify other arguments, such as 'mixture'
# for pcount or 'keyfun' for distsamp
# See the documentation for the associated fitting function and unmarkedFrame
# for what arguments are possible to include for a given model
head(simulate("pcount", formulas=forms, coefs=coefs, design=design, mixture="NB"))
}
}
