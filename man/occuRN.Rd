\name{occuRN}

\alias{occuRN}

\title{Fit the occupancy model of Royle and Nichols (2003)}

\usage{occuRN(formula, data, K=25, starts, method="BFGS", se=TRUE, ...)}

\arguments{
    \item{formula}{double right-hand side formula describing covariates of
        detection and abundance, in that order.}
    \item{data}{Object of class \code{\link{unmarkedFrameOccu}} supplying data
        to the model.}
    \item{K}{the upper summation index used to numerically integrate out the
        latent abundance. This should be set high enough so that it does not
        affect the parameter estimates. Computation time will increase with K.}
    \item{starts}{initial values for the optimization.}
    \item{method}{Optimization method used by \code{\link{optim}}.}
    \item{se}{logical specifying whether or not to compute standard
      errors.}
    \item{\dots}{Additional arguments to optim, such as lower and upper
      bounds}
  }

\description{Fit the occupancy model of Royle and Nichols (2003)}

\details{

This function fits the latent abundance mixture model described in
Royle and Nichols (2003).

The latent abundance of site \eqn{i} is modelled as Poisson:

\deqn{N_i \sim Poisson(\lambda_i)}{N_i ~ Poisson(lambda_i)}

The detection of a single individual in site \eqn{i} during sample
\eqn{j} is modelled as Bernoulli:

\deqn{w_{ij} \sim Bernoulli(r_{ij})}{w_ij ~ Bernoulli(r_ij)}.

Thus, the detection probability for a single site is linked to the
detection probability for an individual by

\deqn{p_{ij} = 1 - (1 - r_{ij}) ^ {N_i}}{p_ij = 1 - (1 - r_ij) ^ N_i}

Covariates of \eqn{\lambda_i}{lambda_i} are modelled with the log link
and covariates of \eqn{r_{ij}}{r_ij} are modelled with the logit link.}

\value{unmarkedFit object describing the model fit.}

\author{Ian Fiske}

\references{

Royle, J. A. and Nichols, J. D. (2003) Estimating Abundance from
Repeated Presence-Absence Data or Point Counts. \emph{Ecology}, 84(3)
pp. 777--790.

}

\examples{

\dontrun{

data(birds)
woodthrushUMF <- unmarkedFrameOccu(woodthrush.bin)
# survey occasion-specific detection probabilities
(fm.wood.rn <- occuRN(~ obsNum ~ 1, woodthrushUMF))

# Empirical Bayes estimates of abundance at each site
re <- ranef(fm.wood.rn)
plot(re)


}

}

\keyword{models}
