\name{goccu}
\alias{goccu}
\title{
Fit multi-scale occupancy models
}
\description{
Fit multi-scale occupancy models as described in Nichols et al. (2008) to
repeated presence-absence data collected using the robust design. This model
allows for inference about occupancy, availability, and detection probability.
}
\usage{
goccu(psiformula, phiformula, pformula, data, linkPsi = c("logit", "cloglog"),
      starts, method = "BFGS", se = TRUE, ...)
}
\arguments{
  \item{psiformula}{
    Right-hand sided formula describing occupancy covariates
}
  \item{phiformula}{
    Right-hand sided formula describing availability covariates
}
  \item{pformula}{
    Right-hand sided formula for detection probability covariates
}
  \item{data}{
    An object of class unmarkedFrameGOccu or unmarkedMultFrame
}
\item{linkPsi}{Link function for the occupancy model. Options are  
      \code{"logit"} for the standard occupancy model or \code{"cloglog"} 
      for the complimentary log-log link, which relates occupancy
      to site-level abundance.
}
\item{starts}{
    Starting values
}
  \item{method}{
    Optimization method used by \code{\link{optim}}
}
  \item{se}{
    Logical. Should standard errors be calculated?
}
\item{\dots}{
    Additional arguments to \code{\link{optim}}, such as lower and upper
    bounds
}
}
\details{
  Primary periods could represent spatial or temporal sampling replicates.
  For example, you could have several spatial sub-units within each site, where each
  sub-unit was then sampled repeatedly. This is a frequent design for eDNA studies.
  Or, you could have multiple primary periods of sampling at each site
  (conducted at different times within a season), each of which contains 
  several secondary sampling periods. In both cases the robust design structure
  can be used to estimate an availability probability in addition to
  detection probability. See Kery and Royle (2015) 10.10 for more details.
}
\value{
  An object of class unmarkedFitGOccu
}
\references{
  Kery, M., & Royle, J. A. (2015). Applied hierarchical modeling in ecology: 
  Volume 1: Prelude and static models. Elsevier Science.
  
  Nichols, J. D., Bailey, L. L., O'Connell Jr, A. F., Talancy, N. W., 
  Campbell Grant, E. H., Gilbert, A. T., Annand E. M., Husband, T. P., & Hines, J. E.
  (2008). Multi-scale occupancy estimation and modelling using multiple detection methods. 
  Journal of Applied Ecology, 45(5), 1321-1329.
}
\author{
  Ken Kellner \email{contact@kenkellner.com}
}

\seealso{
\code{\link{occu}}, \code{\link{colext}},
  \code{\link{unmarkedMultFrame}}, \code{\link{unmarkedFrameGOccu}}
}

\examples{

set.seed(123)
M <- 100
T <- 5
J <- 4

psi <- 0.5
phi <- 0.3
p <- 0.4

z <- rbinom(M, 1, psi)
zmat <- matrix(z, nrow=M, ncol=T)

zz <- rbinom(M*T, 1, zmat*phi)
zz <- matrix(zz, nrow=M, ncol=T)

zzmat <- zz[,rep(1:T, each=J)]
y <- rbinom(M*T*J, 1, zzmat*p)
y <- matrix(y, M, J*T)
umf <- unmarkedMultFrame(y=y, numPrimary=T)

\dontrun{
  mod <- goccu(psiformula = ~1, phiformula = ~1, pformula = ~1, umf)
  plogis(coef(mod))
}

}
