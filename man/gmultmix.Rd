\name{gmultmix}
\alias{gmultmix}
\title{Generalized multinomial N-mixture model}
\description{
A three level hierarchical model for designs involving
repeated counts that yield multinomial outcomes.  Possible data collection
methods include repeated removal sampling and double observer sampling.
The three model parameters are abundance, availability, and detection
probability.
}
\usage{
gmultmix(lambdaformula, phiformula, pformula, data,
mixture = c("P", "NB"), K, starts, method = "BFGS", se = TRUE, ...)
}
\arguments{
    \item{lambdaformula}{Righthand side (RHS) formula describing abundance
        covariates}
    \item{phiformula}{RHS formula describing availability covariates}
    \item{pformula}{RHS formula describing detection covariates}
    \item{data}{An object of class unmarkedFrameGMM}
    \item{mixture}{Either "P" or "NB" for Poisson and Negative Binomial mixing
        distributions.}
    \item{K}{The upper bound of integration}
    \item{starts}{Starting values}
    \item{method}{Optimization method used by \code{\link{optim}}}
    \item{se}{Logical. Should standard errors be calculated?}
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
    }
\details{
The latent transect-level super-population abundance distribution
\eqn{f(M | \mathbf{\theta})}{f(M | theta)} can be set as either a
Poisson or
a negative binomial random variable, depending on the setting of the
\code{mixture} argument. \code{mixture = "P"} or \code{mixture = "NB"}
select
the Poisson or negative binomial distribution respectively.  The mean of
\eqn{M_i} is \eqn{\lambda_i}{lambda_i}.  If \eqn{M_i \sim NB}{M_i ~ NB},
then an
additional parameter, \eqn{\alpha}{alpha}, describes dispersion (lower
\eqn{\alpha}{alpha} implies higher variance).

The number of individuals available for detection at time j
is a modeled as binomial:
\eqn{N_{ij} \sim Binomial(M_i, \mathbf{\phi_{ij}})}{N(i,j) ~
  Binomial(M(i), phi(i,j))}.

The detection process is modeled as multinomial:
\eqn{\mathbf{y_{it}} \sim
Multinomial(N_{it}, \pi_{it})}{y(i,1:J,t) ~
Multinomial(N(i,t), pi(i,1,t), pi(i,2,t), ..., pi(i,J,t))},
where \eqn{\pi_{ijt}}{pi(ijt)} is the multinomial cell probability for
plot i at time t on occasion j.

Cell probabilities are computed via a user-defined function related to the
sampling design. Alternatively, the default functions
\code{\link{removalPiFun}}
or \code{\link{doublePiFun}} can be used for equal-interval removal
sampling or
double observer sampling. Note that the function for computing cell
probabilites
is specified when setting up the data using \code{\link{unmarkedFrameGMM}}.

Parameters \eqn{\lambda}{lambda}, \eqn{\phi}{phi} and \eqn{p}{p} can be
modeled as linear functions of covariates using the log, logit and logit
links
respectively.

}
\value{
An object of class unmarkedFitGMM.
}
\references{
Royle, J. A. (2004) Generalized estimators of avian abundance from count
survey data. \emph{Animal Biodiversity and Conservation} 27,
pp. 375--386.

Chandler, R. B., J. A. Royle, and D. I. King. 2011. Inference about
density and temporary emigration in unmarked populations. Ecology
92:1429-1435.

}
\author{
Richard Chandler \email{rbchan@uga.edu} and Andy Royle
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
\code{\link{unmarkedFrameGMM}} for setting up the data and metadata.
\code{\link{multinomPois}} for surveys where no secondary sampling periods were
used. Example functions to calculate multinomial cell probabilities are
described \code{\link{piFuns}}
}
\examples{

# Simulate data using the multinomial-Poisson model with a
# repeated constant-interval removal design.

n <- 100  # number of sites
T <- 4    # number of primary periods
J <- 3    # number of secondary periods

lam <- 3
phi <- 0.5
p <- 0.3

#set.seed(26)
y <- array(NA, c(n, T, J))
M <- rpois(n, lam)          # Local population size
N <- matrix(NA, n, T)       # Individuals available for detection

for(i in 1:n) {
    N[i,] <- rbinom(T, M[i], phi)
    y[i,,1] <- rbinom(T, N[i,], p)    # Observe some
    Nleft1 <- N[i,] - y[i,,1]         # Remove them
    y[i,,2] <- rbinom(T, Nleft1, p)   # ...
    Nleft2 <- Nleft1 - y[i,,2]
    y[i,,3] <- rbinom(T, Nleft2, p)
    }

y.ijt <- cbind(y[,1,], y[,2,], y[,3,], y[,4,])


umf1 <- unmarkedFrameGMM(y=y.ijt, numPrimary=T, type="removal")

(m1 <- gmultmix(~1, ~1, ~1, data=umf1, K=30))

backTransform(m1, type="lambda")        # Individuals per plot
backTransform(m1, type="phi")           # Probability of being avilable
(p <- backTransform(m1, type="det"))    # Probability of detection
p <- coef(p)

# Multinomial cell probabilities under removal design
c(p, (1-p) * p, (1-p)^2 * p)

# Or more generally:
head(getP(m1))

# Empirical Bayes estimates of super-population size
re <- ranef(m1)
plot(re, layout=c(5,5), xlim=c(-1,20), subset=site\%in\%1:25)


}
\keyword{ model }
