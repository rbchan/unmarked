\name{unmarked-package}
\alias{unmarked}
\alias{unmarked-package}
\title{Models for Data from Unmarked Animals}
\description{Estimate wildlife abundance or occurrence.}
\details{\tabular{ll}{
Package: \tab unmarked\cr
Type: \tab Package\cr
Version: \tab 0.8-7\cr
License: \tab GPL (>= 3)\cr
}

Unmarked estimates wildlife parameters for many popular sampling methods including occurrence and point count data.

\strong{Overview of Model-fitting Functions:}  Unmarked provides
several functions for fitting integrated likelihood models for wildlife
abundance and occurrence to replicated survey data.  

\code{\link{occu}} fits occurrence models with no linkage between abundance and detection (MacKenzie et al. 2006).  

\code{\link{occuRN}} fits abundance models to presence/absence data by exploiting the link between detection probability and abundance (Royle and Nichols 2003).

\code{\link{colext}} fits the mutli-season occupancy model of MacKenzie et al. 2003.

\code{\link{pcount}} fits N-mixture models for repeated count data
(Royle 2004a, Kéry et al 2005).  

\code{\link{distsamp}} fits the distsance sampling model of 
Royle et al. (2004) to distance data recorded in discrete intervals. 

\code{\link{multinomPois}} fits the multinomial-Poisson model of Royle 2004b.

\code{\link{gmultmix}} fits a generalized form of the multinomial-mixture model of Royle 2004a that allows for estimating availability and detection probability.  

All of these functions allow the user to specify covariates that affect the detection and state processes.

\strong{Data:} All data is passed to unmarked's estimation functions
as a formal S4 class called an unmarkedFrame, which has child classes for each
model type. This allows metadata (eg as distance interval cut points, 
measurement units, etc...) to be stored with the response and covariate data.
See \code{\link{unmarkedFrame}} for a detailed description of unmarkedFrames
and how to create them.

\strong{Model Specification:}  Most of \emph{unmarked}'s
model-fitting functions allow specification of covariates for both the
state process and the detection process.  Covariates for the detection
process (at the site or observation level) and the state process (at the site 
level) are specified with a double right-hand sided formula, in that 
order.  Such a formula looks like 
\eqn{~ x1 + x2 + ... + x3 ~ x_1 + x_2 + \ldots + x_n}{} where 
\eqn{x_1}{x1} through \eqn{x_n}{xn} are additive covariates of
the process of interest.  The meaning of these covariates or, what they
model, is full described in the help files for the individual functions
and is not the same for all functions.

\strong{Utility Functions:}  \emph{unmarked} contains several utility
functions for organizing data into the form required by its model-fitting
functions.  \code{\link{csvToUMF}} converts an appropriately
formated comma-separated values (.csv) file to a list containing the
components required by model-fitting functions.}
\references{

MacKenzie, D. I. et al. (2006) \emph{Occupancy Estimation and Modeling}.
Amsterdam: Academic Press.

Royle, J. A. and Nichols, J. D. (2003) Estimating Abundance from
Repeated Presence-Absence Data or Point Counts. \emph{Ecology}, 84(3)
pp. 777--790.

Royle, J. A. (2004a) N-Mixture Models for Estimating Population Size from
Spatially Replicated Counts. \emph{Biometrics} 60, pp. 108--105.

Royle, J. A., D. K. Dawson, and S. Bates (2004) Modeling 
abundance effects in distance sampling. \emph{Ecology} 85, pp. 1591--1597.

Kéry, M., Royle, J. A., and Schmid, H. (2005) Modeling Avaian Abundance from
Replicated Counts Using Binomial Mixture Models. \emph{Ecological
Applications} 15(4), pp. 1450--1461.

Royle, J. A. and Link W. A. (2005) A general class of multinomial
mixture models for anuran calling survey data. \emph{Ecology}, \bold{86}(9),
pp. 2505--2512.}
\docType{package}
\author{Ian Fiske, Richard Chandler, and Andy Royle}
\keyword{package}
