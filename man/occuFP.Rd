\name{occuFP}

\alias{occuFP}

\title{Fit occupancy models when false positive detections occur (e.g., Royle and Link [2006]
and Miller et al. [2011]) }

\usage{occuFP(detformula = ~ 1, FPformula = ~ 1, Bformula = ~ 1, 
stateformula = ~ 1, data, starts, method="BFGS", se = TRUE, engine = "R", ...)}

\arguments{
    \item{detformula}{formula describing covariates of
        detection.}
    \item{FPformula}{formula describing covariates of
        false positive detection probability.}
    \item{Bformula}{formula describing covariates of
        probability detections are certain.}
    \item{stateformula}{formula describing covariates of
        occupancy.}
    \item{data}{An \code{\link{unmarkedFrameOccuFP}} object}
    \item{starts}{Vector of parameter starting values.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{Logical specifying whether or not to compute standard
      errors.}
    \item{engine}{Currently only choice is R.}
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{This function fits the single season occupancy model while allowing for
false positive detections.}

\details{

See \code{\link{unmarkedFrame}} and \code{\link{unmarkedFrameOccuFP}} for a
description of how to supply data to the \code{data} argument.

\code{occuFP} fits an extension of the standard single-season occupancy model (MacKenzie et al. 2002), which allows
false positive detections. The occupancy status of a site is the same way as with the \code{\link{occu}} function, where stateformula
is used to specify factors that lead to differences in occupancy probabilities among sites.

The observation process differs in that both false negative and false positive errors are modeled for observations. The
function allows data to be of 3 types. These types are specified using in \code{unmarkedFrameOccuFP} as type. Occassions
are specified to belong to 1 of the 3 data types and all or a subset of the data types can be combined in the same model.

For type 1 data, the detection process is assumed to fit the assumptions of the standard MacKenzie model
where false negative probabilities are estimated but false positive detections are assumed not to occur. If all of your
data is of this type you should use code{occu} to analyze data. The detection parameter p, which is modeled using the
detformula is the only observation parameter for these data.

For type 2 data, both false negative and false positive detection probabilities are estimated. If all data is of this
type the likelihood follows Royle and Link (2006). Both p (the true positive detection probability) and fp (the false
positive detection probability described by fpformula) are estimated for occassions when this data type occurs

For type 3 data, observations are assumed to include both certain detections (false positives assumed not to occur)
and uncertain detections that may include false positive detections. When only this data type occurs, the estimator
is the same as the multiple detection state model described in Miller et al. (2011). Three observation parameters occur
for this data type: p - true positive detection probability, fp - false positive detection probability, and b - the
probability a true positive detection was designated as certain.

When both type 1 and type 2 data occur, the estimator is equivalent to the multiple detection method model described
in Miller et al. (2011). The frog data example in the same paper uses an analysis where type 1 (dipnet surveys) and
type 3 (call surveys) data were used.

Data in the y matrix of the unmarked frame should be all 0s and 1s for
type 1 and type 2 data. For type 3 data, uncertain detections are given
a value of 1 and certain detections a value of 2.  }


\value{unmarkedFitOccuFP object describing the model fit.}

\references{

MacKenzie, D. I., J. D. Nichols, G. B. Lachman, S. Droege,
  J. Andrew Royle, and C. A. Langtimm. 2002. Estimating Site Occupancy Rates
  When Detection Probabilities Are Less Than One. Ecology 83: 2248-2255.

Miller, D.A., J.D. Nichols, B.T. McClintock, E.H.C. Grant, L.L. Bailey,
and L.A. Weir. 2011. Improving occupancy estimation when two types of
observational error occur: non-detection and species
misidentification. Ecology 92:1422-1428.

Royle, J.A., and W.A. Link. 2006. Generalized site occupancy models
allowing for false positive and false negative errors. Ecology
87:835-841.


}

\author{David Miller}

\seealso{\code{\link{unmarked}}, \code{\link{unmarkedFrameOccuFP}},
    \code{\link{modSel}}, \code{\link{parboot}}}


\examples{

n = 100
o = 10
o1 = 5
y = matrix(0,n,o)
p = .7
r = .5
fp = 0.05
y[1:(n*.5),(o-o1+1):o] <- rbinom((n*o1*.5),1,p)
y[1:(n*.5),1:(o-o1)] <- rbinom((o-o1)*n*.5,1,r)
y[(n*.5+1):n,(o-o1+1):o] <- rbinom((n*o1*.5),1,fp)
type <- c((o-o1),o1,0)  ### vector with the number of each data type
site <- c(rep(1,n*.5*.8),rep(0,n*.5*.2),rep(1,n*.5*.2),rep(0,n*.8*.5))
occ <- matrix(c(rep(0,n*(o-o1)),rep(1,n*o1)),n,o)
site <- data.frame(habitat = site)
occ <- list(METH = occ)

umf1 <- unmarkedFrameOccuFP(y,site,occ, type = type)

m1 <- occuFP(detformula = ~ METH, FPformula = ~1,
             stateformula = ~ habitat, data = umf1)
predict(m1, type = 'fp')
coef(m1)
confint(m1, type = 'det')

}

\keyword{models}
