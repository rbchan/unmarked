\name{cooccu}
\alias{cooccu}
\title{
Co-occurrence model
}
\description{
A co-occurrence model allowing for (possibly) dependent occurrence,
false negative errors, and false positive errors.
}
\usage{
cooccu(psiformulaA = ~1, psiformulaB = ~1, gammaformula = NULL,
       pformulaA = ~1, pformulaB = ~1, fpformulaA = NULL,
       fpformulaB = NULL, data, interaction = c("dominance", "general"),
       fpmodel = c("full", "confusion"), starts, method = "BFGS",
       se = TRUE, ...)
}
\arguments{
  \item{psiformulaA}{
    Occupancy model for species A
}
  \item{psiformulaB}{
    Occupancy model for species B
}
  \item{gammaformula}{
    Optional interaction model describing the dependence between the
    occupancy of the two species
}
  \item{pformulaA}{
    Detection (false-negative probability) model for species A
}
  \item{pformulaB}{
    Detection (false-negative probability) model for species B
}
  \item{fpformulaA}{
    Optional false-positive error model for species B
}
  \item{fpformulaB}{
    Optional false-positive error model for species B
}
  \item{data}{
    An unmarkedFrameCo object
}
  \item{interaction}{
    Either "dominance" for the case where occurrence of species B
    depends upon species A but not vice-versa
    (Waddle et al. 2010, Richmond et al. 2010), or "general" for model in
    which species can attract or repel one another (MacKenzie et al. 2004).
}
  \item{fpmodel}{
    The "full" model allows false positives to occur anytime a species
    is not present. The "confusion" model only allows false positives to
    occurr when one species is confused with the other, i.e. they cannot
    happen when both species are absent.
}
  \item{starts}{
    Starting values for the parameters. The order of these is the same
    as the order of the formulas.
}
  \item{method}{
    Optimization method. See \code{\link{optim}}
}
  \item{se}{
    Should the standard errors be computed?
}
  \item{\dots}{
    Additional arguments passed to \code{\link{optim}}
}
}
\details{
  If gammaformula is NULL, then the two species are assumed to occur
  independently of one another. For the "dominance" model, an estimate
  of gamma=1 indicates no influence of the presence of species A on the
  occurrence of species B. Note that gamma is a probability in this
  model, so gamma=1 is the upper boundary. For the "general" model,
  gamma can take any positive real value, subject to the constraints
  described in MacKenzie et al. (2004). An estimate of gamma=1 indicates
  that the two species occur independently.

  Both the false-positive and the false-negative probabilities can be
  modeled as functions of the occurrence of the other species. For
  instance, if you think the false-negative rate of species B depends
  upon the presence of species A (perhaps because species A is its
  predator), then you can specify this using: pmodelB~SppA. The
  variables "SppA" and "SppB" do not need to be provided -- they are
  generated automatically.
}
\value{
  An object of class \code{unmarkedFitCo}.
}
\references{
  MacKenzie, D., L. Bailey, and J. Nichols, 2004. Investigating species
  co-occurrence patterns when species are detected imperfectly. Journal
  of Animal Ecology 73:546--555.

  Royle, J. and W. Link, 2006. Generalized site occupancy models
  allowing for false positive and false negative errors. Ecology
  87:835--841.

  Waddle, J., R. Dorazio, S. Walls, K. Rice, J. Beauchamp, M. Schuman,
  and F. Mazzotti, 2010. A new parameterization for estimating
  co-occurrence of interacting species. Ecological Applications
  20:1467--1475.

  Richmond, O., J. Hines, and S. Beissinger, 2010. Two-species occupancy
  models: a new parameterization applied to co-occurrence of secretive
  rails. Ecological Applications 20:2036--2046.

  Miller, D., J. Nichols, B. McClintock, E. Grant, L. Bailey, and L. Weir,
  2011. Improving occupancy estimation when two types of observational
  error occur: non-detection and species misidentification. Ecology
  92:1422--1428.
}
\author{
  Richard Chandler
}
\note{
%%  ~~further notes~~
}
\section{Warning}{
  The default is to assume no false positives. Caution is advised
  when trying to estimate both the false negative and false positive
  rates because the parameters are only weakly-identifiable (Royle and Link
  2006). This problem can be averted by assuming that the false-negative
  rate is higher than the false-positive rate (providing good starting
  values should help ensure this). Alternatively, if false-positives
  were known not to occur at some sites, then you can specify this when
  creating the \code{unmarkedFrameCo} object.
}
\seealso{
  \code{\link{occu}, \code{unmarkedFrameCo}}
}
\examples{
nSites <- 100
nReps <- 5
yA <- matrix(rbinom(nSites*nReps, 1, 0.1), nSites, nReps)
yB <- matrix(rbinom(nSites*nReps, 1, 0.1), nSites, nReps)
umf <- unmarkedFrameCo(yA=yA, yB=yB)

fm1 <- cooccu(psiformulaA=~1, psiformulaB=~1, #gammaformula=~1,
              data=umf)
fm1

fm2 <- cooccu(psiformulaA=~1, psiformulaB=~1, gammaformula=~1,
              data=umf)
fm2

fm3 <- cooccu(psiformulaA=~1, psiformulaB=~1, gammaformula=NULL,
              pformulaA=~SppB,
              data=umf)
fm3

fm4 <- cooccu(psiformulaA=~1, psiformulaB=~1, gammaformula=NULL,
              pformulaB=~SppA, #se=FALSE,
              data=umf)
fm4

fm5 <- cooccu(psiformulaA=~1, psiformulaB=~1, gammaformula=NULL,
              fpformulaA=~SppB, #se=FALSE,
              data=umf)
fm5

fm6 <- cooccu(psiformulaA=~1, psiformulaB=~1, gammaformula=NULL,
              fpformulaA=~SppB,
              fpformulaB=~SppA, se=FALSE,
              starts=c(4, 4, 0, 0, -3, 0, -3, 0),
              data=umf)
fm6

fm7 <- cooccu(psiformulaA=~1, psiformulaB=~1, gammaformula=NULL,
              fpformulaA=~1,
              fpformulaB=~1, #se=FALSE,
              data=umf)
fm7

debugonce(cooccu)
source("cooccu.R")
source("getDesign.R")
source("unmarkedEstimate.R")
source("unmarkedFrame.R")
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
