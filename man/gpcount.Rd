\name{gpcount}
\alias{gpcount}
\title{
Generalized binomial N-mixture model for repeated count data
}
\description{
Fit the model of Chandler et al. (2011) to repeated count data collected
using the robust design. This model allows for inference about
population size, availability, and detection probability.
}
\usage{
gpcount(lambdaformula, phiformula, pformula, data,
mixture = c("P", "NB"), K, starts, method = "BFGS", se = TRUE,
engine = c("C", "R"), ...)
}
\arguments{
  \item{lambdaformula}{
    Right-hand sided formula describing covariates of abundance.
}
  \item{phiformula}{
    Right-hand sided formula describing availability covariates
}
  \item{pformula}{
    Right-hand sided formula for detection probability covariates
}
  \item{data}{
    An object of class unmarkedFrameGPC
}
  \item{mixture}{
    Either "P" or "NB" for Poisson and negative binomial distributions
}
  \item{K}{
    The maximum possible value of M, the super-population size.
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
  \item{engine}{
    Either "C" or "R" for the C++ or R versions of the likelihood. The C++
    code is faster, but harder to debug.
}
  \item{\dots}{
    Additional arguments to \code{\link{optim}}, such as lower and upper
    bounds
}
}
\details{
  The latent transect-level super-population abundance distribution
  \eqn{f(M | \mathbf{\theta})}{f(M | theta)} can be set as either a
  Poisson or a negative binomial random variable, depending on the
  setting of the \code{mixture} argument. The expected value of
  \eqn{M_i} is \eqn{\lambda_i}{lambda_i}.  If \eqn{M_i \sim NB}{M_i ~ NB},
  then an additional parameter, \eqn{\alpha}{alpha}, describes
  dispersion (lower \eqn{\alpha}{alpha} implies higher variance).

  The number of individuals available for detection at time j
  is a modeled as binomial:
  \eqn{N_{ij} \sim Binomial(M_i, \mathbf{\phi_{ij}})}{N(i,j) ~
    Binomial(M(i), phi(i,j))}.

  The detection process is also modeled as binomial:
  \eqn{y_{ikj} \sim Binomial(N_{ij}, p_{ikj})}{y(i,k,j) ~
    Binomial(N(i,t), p(i,k,j))}.

  Parameters \eqn{\lambda}{lambda}, \eqn{\phi}{phi} and \eqn{p}{p} can be
  modeled as linear functions of covariates using the log, logit and logit
  links respectively.
}
\value{
  An object of class unmarkedFitGPC
}
\references{
  Royle, J. A. 2004. N-Mixture models for estimating population size from
  spatially replicated counts. \emph{Biometrics} 60:108--105.

  Chandler, R. B., J. A. Royle, and D. I. King. 2011. Inference about
  density and temporary emigration in unmarked populations. Ecology
  92:1429-1435.
}
\author{
  Richard Chandler \email{rbchan@uga.edu}
}
\note{
In the case where availability for detection is due to random temporary
emigration, population density at time j, D(i,j), can be estimated by
N(i,j)/plotArea.

This model is also applicable to sampling designs in which the local
population size is closed during the J repeated counts, and availability
is related to factors such as the probability of vocalizing. In this
case, density can be estimated by M(i)/plotArea.

If availability is a function of both temporary emigration and other
processess such as song rate, then density cannot be directly estimated,
but inference about the super-population size, M(i), is possible.

Three types of covariates can be supplied, site-level,
site-by-year-level, and observation-level. These must be formatted
correctly when organizing the data with \code{\link{unmarkedFrameGPC}}
}


\seealso{
\code{\link{gmultmix}}, \code{\link{gdistsamp}},
  \code{\link{unmarkedFrameGPC}}
}

\examples{
set.seed(54)

nSites <- 20
nVisits <- 4
nReps <- 3

lambda <- 5
phi <- 0.7
p <- 0.5

M <- rpois(nSites, lambda) # super-population size

N <- matrix(NA, nSites, nVisits)
y <- array(NA, c(nSites, nReps, nVisits))
for(i in 1:nVisits) {
    N[,i] <- rbinom(nSites, M, phi) # population available during vist j
}
colMeans(N)

for(i in 1:nSites) {
    for(j in 1:nVisits) {
        y[i,,j] <- rbinom(nReps, N[i,j], p)
    }
}

ym <- matrix(y, nSites)
ym[1,] <- NA
ym[2, 1:nReps] <- NA
ym[3, (nReps+1):(nReps+nReps)] <- NA
umf <- unmarkedFrameGPC(y=ym, numPrimary=nVisits)

\dontrun{
fmu <- gpcount(~1, ~1, ~1, umf, K=40, control=list(trace=TRUE, REPORT=1))

backTransform(fmu, type="lambda")
backTransform(fmu, type="phi")
backTransform(fmu, type="det")
}

}
